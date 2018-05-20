# -*- coding: utf-8 -*-

import matplotlib as mpl

mpl.use(u'TkAgg')
mpl.interactive(False)

import sys
import numpy as np
from helixmc.pose import HelixPose
from helixmc.util import frames2params, frames2params_3dna
from helixmc.random_step import RandomStepSimple
from lmfit import Minimizer, Parameters, report_fit, minimize
import matplotlib.pyplot as plt
# ChromatinMC modules:
import FileIO as fileio
import NucleosomeMC as nMC

plt.interactive(False)
np.set_printoptions(formatter={'float': '{: 0.3f}, '.format})


def create_curved_linker(p=None):
    if p is None:
        p = Parameters()
        p.add('L_bp', value=400)
        p.add('n_nuc', value=2)
        p.add('nld_A', value=100)
        p.add('NRL', value=197)
        p.add('Unwrapped_bp', value=30)
        p.add('Pitch_bp', value=10.380450417509)
        p.add('G_kT', value=71.8913935387115)
        p.add('Phase_rad', value=3.75743188237752)
        p.add('Period', value=0.787822017675862)
        p.add('RelFreq', value=1.09559893586307)

    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')

    linkerlength_bp = p['NRL'] - 147 + p['Unwrapped_bp']

    params = np.tile(random_step.params_avg, (linkerlength_bp, 1))

    i = np.linspace(-linkerlength_bp // 2, linkerlength_bp // 2 - 1, linkerlength_bp)
    twist = 0 * i + 2 * np.pi / p['Pitch_bp']

    roll = p['Phase_rad'] + p['RelFreq'] * 2 * np.pi * i / p['Pitch_bp']
    roll = np.cos(roll)

    modulation = np.cos(p['Period'] * 2.0 * np.pi * i / linkerlength_bp)
    roll = roll * modulation

    energy = np.sqrt(np.sum(roll ** 2)) / 0.115  # 0.115 = approximate twist stiffness
    # Chou2014PloSCmpBiol: sigma_roll = 5' so k = kT/sigma^2 =4.1/(5*pi/180)^2 = 0.01857
    roll = roll * np.sqrt(p['G_kT']) / energy

    params[:, 5] = twist
    params[:, 4] = roll

    return params


def replace_linker(dna, dyads, linker_params):
    '''
    Replace basepair parameters for linker DNA

    Parameters
    ----------
    dna : basepair step parameters for DNA
    dyads : dyad index in DNA (bp)
    free_bps : array of bps; 0 = nucleosomal DNA, 1 = free DNA
    linker_params : basepair step parameters for linker DNA
    Returns
    -------
    dna_params : ndarray of (N, 6)
        The basepair step parameters of DNA
    free_bps : array of bps; 0 = nucleosomal DNA, 1 = free DNA
    '''
    params = dna.params
    for d1, d2 in zip(dyads, dyads[1:]):
        linker_mid = d1 + (d2 - d1) // 2
        linker_start = linker_mid - len(linker_params) // 2 - 1
        params[linker_start: linker_start + len(linker_params)] = linker_params
    dna = HelixPose(params)
    return dna


def insert_nucs(dna, dyads, pdb='1KX5'):
    '''
    Fix basepair parameters for nucleosomes in DNA at dyad positions

    Parameters
    ----------
    dna: HelixMC pose
    dyads : dyad positions in DNA (bp)
    Returns
    -------
    dna : HelixMC pose
    nucl : NucleosomeMC pose
    '''
    #    # Get nucleosome
    nucl = nMC.NucPose()
    nucl.from_file(fileio.change_extension(pdb, '3DNA'))

    params = dna.params
    length = len(nucl.dna.params)
    new_dyads = []
    for dyad in dyads:
        start = dyad - nucl.dyad
        end = start + length
        if (start >= 0 and (end < len(params))):
            params[start:end] = nucl.dna.params
            new_dyads.append(dyad)
    dna = HelixPose(params)
    return dna, nucl, np.asarray(new_dyads)


def create_unfolded_fiber(fiber_pars):
    '''
    Create DNA configuration, including nucleosomes

    '''

    n_bp = fiber_pars['L_bp']
    n_nucs = fiber_pars['n_nuc']
    NRL = fiber_pars['NRL']

    dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
    dyads = (dyads + n_bp // 2).astype(int)

    random_step = RandomStepSimple.load_gaussian_params('DNA_gau.npy')
    params = np.tile(random_step.params_avg, (n_bp - 1, 1))

    dna = HelixPose(params)
    dna, nucl, dyads = insert_nucs(dna, dyads)

    return dna, dyads, nucl


def get_casted_fiber_frames(par):
    '''
    Get the frame coordinates that define a stacked, folded fiber

    Returns
    -------
    n_ofs : ndarray of (N, 4, 3)
        The coordinates of the nucleosome frames
    d_ofs : ndarray of (N, 4, 3)
        The coordinates of the dyad frames
    '''
    fiber_start = par['fiber_start'].value
    r = par['diameter_A'].value / 2.0 - 65.0
    nld = par['nld_A'].value
    rise = par['rise_A'].value
    phi = np.sqrt(rise ** 2 - nld ** 2) / r
    if fiber_start is 2:
        phi = 0.5 * phi + np.pi
    if par['chirality'].value < 0:
        phi *= -1

    cm = []
    for i in np.arange(par['n_nuc'].value + fiber_start):
        cm.append(np.asarray([r * np.cos(i * phi), r * np.sin(i * phi), i * nld]))

    n_ofs = []
    for cm1, cm2 in zip(cm, cm[fiber_start:]):
        Nx = par['face'].value * cm1 * [-1, -1, 0]
        Nx /= np.linalg.norm(Nx)
        Nz = cm2 - cm1
        Nz /= np.linalg.norm(Nz)
        Ny = np.cross(Nx, Nz)
        Ny /= np.linalg.norm(Nz)
        Nx = np.cross(Nz, Ny)
        Nx /= np.linalg.norm(Nx)
        frame = np.array([Nx, Ny, -Nz])
        n_ofs.append(nMC.join_o_f(cm1, np.transpose(frame)))

    # n_coords = []
    # for n_of in n_ofs:
    #     n_coords.append(nMC.of2axis(n_of, length=20))
    #
    # filename = fileio.get_filename(root='array', ext='pov', incr=True)
    # fileio.create_pov(filename, n_coords, range_A=[500, 800], offset_A=[0, 0, 60], show=True, width_pix=500)

    return n_ofs


def init_fixed_stack_params2(par, nucl, plot=True):
    n_ofs = get_casted_fiber_frames(par)
    nucl_dyad_of = nMC.join_o_f(nucl.dna.coords[nucl.dyad], nucl.dna.frames[nucl.dyad])
    new_nucl_dyad_ofs = []
    for n_of in n_ofs:
        tf = nMC.get_transformation(nucl.of, target=n_of)
        new_nucl_dyad_ofs.append(nMC.apply_transformation(nucl_dyad_of, tf))
    stack_params = nMC.ofs2params(new_nucl_dyad_ofs[0], new_nucl_dyad_ofs[1], _3dna=False)

    if plot:
        dna, dyads, _ = create_unfolded_fiber(fiber_pars=par)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for n_of, dyad in zip(n_ofs, dyads):
            tf = nMC.get_transformation(nucl.of, target=n_of)

            dna_coords = (nMC.apply_transformation(nucl.dna.coords, tf))
            ax.scatter(dna_coords[:, 0], dna_coords[:, 1], dna_coords[:, 2])
        plt.show()

    return stack_params


def init_fixed_stack_params(par, nucl, plot=True):
    folded_fiber_dna, dyads, _ = create_casted_fiber(par, nucl)
    stack_params = get_stack_pars(folded_fiber_dna.coords, folded_fiber_dna.frames, dyads[0], dyads[1])
    return stack_params


def create_casted_fiber(par, nucl):
    dna, dyads, _ = create_unfolded_fiber(fiber_pars=par)
    params = dna.params * 0
    w = np.zeros(len(dna.coords))
    n_ofs = get_casted_fiber_frames(par)
    origin_of = np.concatenate(([np.zeros(3)], np.eye(3)))
    length = len(nucl.dna.params)
    unwrapped = int(par['Unwrapped_bp'])

    # from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    #
    # n_coords = []
    # n_coords.append(nMC.of2axis(nucl.of, length=30))
    # start_coords = []

    for n_of, dyad in zip(n_ofs, dyads):
        tf = nMC.get_transformation(nucl.of, target=n_of)

        # dna_coords = (nMC.apply_transformation(nucl.dna.coords, tf))
        # n_coords.append(dna_coords)
        # start_coords.append(dna_coords[0])
        # ax.scatter(dna_coords[:, 0], dna_coords[:, 1], dna_coords[:, 2])

        start_bp = dyad - nucl.dyad
        params[start_bp:start_bp + length] = nucl.dna.params

        dna = HelixPose(params)
        start_of = nMC.get_of(dna, start_bp)

        params[start_bp - 1] = nMC.ofs2params(origin_of, nMC.apply_transformation(start_of, tf), _3dna=True)

        dna = HelixPose(params)
        end_of = nMC.get_of(dna, start_bp + length)
        params[start_bp + length] = nMC.ofs2params(end_of, origin_of, _3dna=True)

        w[start_bp + unwrapped:start_bp + length - unwrapped] += 1
        # ax_coords = nMC.of2axis(nucl.of, length=30)
        # n_coords.append(nMC.apply_transformation(nucl.dna.coords, tf))
        # n_coords.append(ax_coords)
    w = np.transpose(np.asarray([w, w, w]))

    # n_coords.insert(0, start_coords)
    # filename = fileio.get_filename(root='cast_array', ext='pov', incr=True)
    # fileio.create_pov(filename, n_coords, range_A=[500, 800], offset_A=[0, 0, 60], show=True, width_pix=500)

    # fiber_dna_coords = HelixPose(params).coords*w
    # ax.scatter(fiber_dna_coords[:, 0], fiber_dna_coords[:, 1], fiber_dna_coords[:, 2], color='k', s=1)
    # r = 300
    # ax.set_xlim(-r / 2, r / 2)
    # ax.set_ylim(-r / 2, r / 2)
    # ax.set_zlim(0, r)
    # plt.show()

    return dna, dyads, w


def residual_stack_coords(fiber_par, dyad0_of, cast_coords, w):
    dna, dyads, _ = create_folded_fiber(fiber_par)

    tf = nMC.get_transformation(nMC.get_of(dna, dyads[0]), target=dyad0_of)
    coords = nMC.apply_transformation(dna.coords, tf)

    res = (coords - cast_coords) * w
    sys.stdout.write('.')
    return res


def residual_stack_params(fiber_par, fixed_stack_params, w):
    folded_fiber, dyads, nucl = create_folded_fiber(fiber_par)

    stack_params = get_stack_pars(folded_fiber.coords, folded_fiber.frames, dyads[0], dyads[1], nucl)

    res = (stack_params - fixed_stack_params) * w
    sys.stdout.write('.')
    return res


def get_stack_pars(coords, frames, dyad1, dyad2, nucl):
    n_of1 = nMC.get_nuc_of(coords, frames, dyad1, nucl)
    n_of2 = nMC.get_nuc_of(coords, frames, dyad2, nucl)
    stack_params = nMC.ofs2params(n_of2, n_of1, _3dna=True)
    return stack_params


def create_folded_fiber(fiber_pars):
    dna, dyads, nucl = create_unfolded_fiber(fiber_pars=fiber_pars)
    linker_params = create_curved_linker(p=fiber_pars)
    dna = replace_linker(dna, dyads, linker_params)

    return dna, dyads, nucl


def scan_fiber_param(fiber_pars, iter_param, iter_vals):
    setnr = 0
    image_files = []
    report_file = fileio.get_filename(incr=True)
    fileio.report_progress(len(iter_vals), title='scan_fiber_param', init=True)

    # from mpl_toolkits.mplot3d import Axes3D
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')

    for i in iter_vals:
        fileio.report_progress(setnr, title='scan_fiber_param\n')

        filename = fileio.get_filename(incr=True, sub=True)
        fiber_pars[iter_param].value = i

        dna, dyads, nucl = create_unfolded_fiber(fiber_pars=fiber_pars)
        fiber_dna, dyads, w = create_casted_fiber(fiber_pars, nucl)

        # ax.scatter(fiber_dna.coords[:, 0], fiber_dna.coords[:, 1], fiber_dna.coords[:, 2])
        # plt.show()

        fiber_pars['ref_frame'].value = dyads[0]
        ref_of = nMC.get_of(fiber_dna, fiber_pars['ref_frame'].value)

        try:
            out = minimize(residual_stack_coords, fiber_pars, args=(ref_of, fiber_dna.coords, w), method='nelder')
            report_fit(out)
            fiber_pars = out.params
            coords, dna, tf = create_folded_fiber(fiber_pars, ref_of)
            image_files.append(
                fileio.create_pov(filename, [dna.coords], range_A=[400, 400], offset_A=[0, 0, 150], show=True))
            fileio.write_xlsx_row(filename, str(setnr), fiber_pars, report_file=report_file)
            dna.write2disk(fileio.get_filename(sub=True, ext='npz'))
        except:
            print('Fit did not converge')
        setnr += 1
    # fileio.create_movie(image_files, filename=report_file)
    return


def main(pars):
    root = '{0}st{1}x{2}'.format(pars['fiber_start'].value, pars['n_nuc'].value, pars['NRL'].value)
    filename = fileio.get_filename(incr=True, root=root)

    # create fiber with straight linker DNA
    unfolded_fiber, dyads, nucl = create_unfolded_fiber(pars)
    # fileio.plot_dna(unfolded_fiber, range_nm=50, wait=1, origin_index=dyads[0])

    # create fiber with positioned nucleosomes, but missing linker DNA
    casted_fiber, _, w = create_casted_fiber(pars, nucl)
    # fileio.plot_dna(casted_fiber, range_nm=50, wait=5, origin_index=0)

    bases2 = []
    bases1 = []
    bases3 = []
    for dyad in dyads:
        bases1.append(casted_fiber.coords[dyad + nucl.fixed[0]])
        bases2.append(casted_fiber.coords[dyad + nucl.fixed[-1]])
        for i in nucl.fixed[1:-1]:
            bases3.append(casted_fiber.coords[dyad + i])

    print(fileio.create_pov(filename, [bases1, bases2, casted_fiber.coords, bases3], range_A=[500, 500],
                            offset_A=[0, 0, 150], show=True, width_pix=3000, colors='rbky', radius=[11, 11, 10, 11]))

    print('\n fixed stack_params:')
    for dyad0, dyad1 in zip(dyads, dyads[pars['fiber_start'].value:]):
        fixed_stack_params = get_stack_pars(casted_fiber.coords, casted_fiber.frames, dyad0, dyad1, nucl)
        print(fixed_stack_params)

    return

    # fit linker DNA parameters
    dyad0_of = nMC.get_of_2(casted_fiber.coords, casted_fiber.frames, dyads[0])
    out = minimize(residual_stack_coords, pars, args=(dyad0_of, casted_fiber.coords, w), method='nelder')
    # report_fit(out)
    print('\n fit stack_coords:')
    pars = out.params

    # create folded fiber with bent linker DNA
    folded_fiber, _, _ = create_folded_fiber(pars)
    fileio.plot_dna(folded_fiber, range_nm=50, wait=1, origin_index=dyads[0], save=True)
    print(get_stack_pars(folded_fiber.coords, folded_fiber.frames, dyads[0], dyads[1], nucl))

    # write to file
    fileio.write_xlsx_row(fileio.get_filename(sub=True, incr=True, ext='npz'), 0, pars, report_file=filename)
    folded_fiber.write2disk(fileio.get_filename(sub=True, ext='npz'))
    print(filename)
    return

    w = [1, 1, 1, 10, 10, 10]
    out = minimize(residual_stack_params, pars, args=(fixed_stack_params, w), method='nelder')
    print('\n fit stack_params:')
    pars = out.params
    folded_fiber, _, _ = create_folded_fiber(out.params)
    fileio.plot_dna(folded_fiber, range_nm=50, wait=10, origin_index=dyads[0])
    print(get_stack_pars(folded_fiber.coords, folded_fiber.frames, dyads[0], dyads[1], nucl))

    return

    iter_vals = np.linspace(20, 100, 11)
    scan_fiber_param(pars, 'nld_A', iter_vals)

    unfolded_fiber, dyads, nucl = create_unfolded_fiber(fiber_pars=pars)
    casted_fiber, dyads, w = create_casted_fiber(pars, nucl)
    print (fileio.create_pov(fileio.get_filename(incr=True), [casted_fiber.coords], range_A=[1000, 1500],
                             offset_A=[0, 0, 150], show=True))

    return


if __name__ == "__main__":
    # execute only if run as a script
    fiber_par = Parameters()
    # fiber parameters
    fiber_par.add('L_bp', value=4000, vary=False)
    fiber_par.add('n_nuc', value=4, vary=False)
    fiber_par.add('NRL', value=197, vary=False)
    fiber_par.add('Unwrapped_bp', value=20, min=0, max=60, vary=False)
    # linker DNA parameters
    fiber_par.add('Pitch_bp', value=10.780, min=9.4, max=11.4, vary=True)
    fiber_par.add('G_kT', value=27.891, min=0, max=100, vary=True)
    fiber_par.add('Phase_rad', value=3.757, min=0, max=2 * np.pi, vary=True)
    fiber_par.add('Period', value=0.4787, min=0, max=2, vary=True)
    fiber_par.add('RelFreq', value=1.01, min=0.95, max=1.05, vary=True)
    # stacking parameters
    fiber_par.add('diameter_A', value=330, vary=False)
    fiber_par.add('rise_A', value=90, vary=False)
    fiber_par.add('chirality', value=-1, vary=False)
    fiber_par.add('face', value=1, vary=False)

    fiber_par.add('fiber_start', value=2)
    if fiber_par['fiber_start'].value is 1:
        fiber_par.add('nld_A', value=55, vary=False)
    else:
        fiber_par.add('nld_A', value=35, vary=False)
    main(fiber_par)
