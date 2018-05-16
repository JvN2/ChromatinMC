import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose
import numpy as np
from lmfit import Parameters
import pandas as pd
import easygui as eg

# ChromatinMC modules:
import NucleosomeMC as nMC
import POVutils as pov
import FileIO as fileio
import RunMC as rMC
import FiberMC as fMC

default_step_file = 'C:\\Python27\\Lib\\site-packages\\helixmc\\data\\DNA_gau.npy'
default_folder = 'E:\\users\\'
kT = 41


def plot_gf(filename, calc=False):
    filename = fileio.change_extension(filename, 'xlsx')
    print('>>> Current file: {}'.format(filename))

    if calc:
        default_params = np.load(default_step_file)
        p0 = default_params[0]
        k = np.linalg.inv(default_params[1:])
        pars, _ = fileio.read_xlsx_row(filename, 0)
        dna, dyads, nucl = fMC.create_unfolded_fiber(fiber_pars=pars)

        # Get stack and wrap parameters
        fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coords, nucl.dna.frames, nucl.dyad, nucl.fixed)

        fiber_start = pars['fiber_start'].value
        fiber_dna, dyads, w = fMC.create_casted_fiber(pars, nucl)
        fixed_stack_params = fMC.get_stack_pars(fiber_dna.coords, fiber_dna.frames, dyads[0], dyads[fiber_start])

        e_stack_kT = pars['e_stack_kT'].value
        e_wrap_kT = pars['e_wrap_kT'].value
        fiber_start = pars['fiber_start'].value

        sets, files, _ = fileio.contents_xlsx(filename)
        fileio.report_progress(len(sets), title='plot_gf', init=True)
        i = 0
        g_nuc_kT_all = []
        for set, file in zip(sets, files):
            fileio.report_progress(i)
            pars, _ = fileio.read_xlsx_row(filename, set)
            force = pars['F_pN'].value
            dna = HelixPose.from_file(fileio.change_extension(file, 'npz'))
            g_nuc_kT, names = rMC.get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl,
                                                   e_wrap_kT, e_stack_kT, fiber_start, p0, k, force)
            g_nuc_kT_all.append(g_nuc_kT)
            for g, name in zip(g_nuc_kT, names):
                pars[name].value = g
            i += 1
    else:
        names = ['g_dna_kT', 'g_wrap_kT', 'g_stack_kT', 'g_work_kT']
        g_nuc_kT_all = []
        for name in names:
            g_nuc_kT_all.append(fileio.read_xlsx_collumn(filename, name))
        g_nuc_kT_all = np.transpose(g_nuc_kT_all)

    force = fileio.read_xlsx_collumn(filename, 'F_pN')
    plt.close()
    plt.figure(figsize=(4, 3))
    plt.axes([0.15, 0.15, .8, .75])

    plt.semilogx(force, np.asarray(g_nuc_kT_all))
    plt.xlim(0.08, 12)
    plt.xlabel('F (pN)')
    plt.ylim(-75, 75)
    plt.ylabel(r'$\Delta$G (kT)')
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    for i, name in enumerate(names):
        names[i] = name.split('_')[1]
    plt.legend(names)
    plt.draw()
    plt.pause(5)
    filename = fileio.change_extension(filename, '_gf.jpg')
    plt.savefig(filename, dpi=600, format='jpg')
    return


def create_dummy_dna():
    pars = Parameters()
    # Params that define the nucleosomal array
    pars.add('Unwrapped_bp', value=0)
    pars.add('g_dna_kT', value=0)
    pars.add('g_wrap_kT', value=0)
    pars.add('g_stack_kT', value=0)
    pars.add('g_work_kT', value=0)

    pars.add('L_bp', value=1000)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=4)
    pars.add('dyad0_bp', value=0)
    pars.add('Pitch_bp', value=10.380, min=9, max=11, vary=True)
    pars.add('G_kT', value=0, min=0, max=100, vary=True)
    pars.add('Phase_rad', value=3.757, min=0, max=2 * np.pi, vary=True)
    pars.add('Period', value=0.787, min=0, max=2, vary=True)
    pars.add('RelFreq', value=1.0, min=0.95, max=1.05, vary=True)

    # Params that define the folded fiber
    pars.add('diameter_A', value=330)
    pars.add('rise_A', value=100)
    pars.add('nld_A', value=25)
    pars.add('chirality', value=-1)
    pars.add('face', value=1)
    pars.add('e_wrap_kT', value=3)
    pars.add('e_stack_kT', value=0)

    # Params that are typically varied between simulations
    pars.add('NRL', value=197)
    pars.add('fiber_start', value=2)
    pars.add('ref_frame', value=0, vary=False)

    filename = fileio.get_filename(incr=True, root='dummy')

    for unwrap in range(0, 146):
        pars['Unwrapped_bp'].value = unwrap
        # dna, dyads, nucl = fMC.create_nuc_array(p=pars)
        _, dna, _, dyads = fMC.create_folded_fiber(pars, ref_of=None)
        pars['dyad0_bp'].value = dyads[0]
        fileio.plot_dna(dna, update=(unwrap > 0), title='Unwrapped = {:.1f} bp\n'.format(unwrap), save=True)
        pose_file = fileio.get_filename(sub=True, ext='npz', incr=True)
        fileio.write_xlsx_row(pose_file, unwrap, pars, report_file=filename)
        dna.write2disk(pose_file)
    return filename


def plot_energy(filename):
    params = np.load(default_step_file)
    p0 = params[0]
    cov = params[1:]
    k = np.linalg.inv(cov)

    sets, files, _ = fileio.contents_xlsx(filename)
    dna = HelixPose.from_file(fileio.change_extension(files[0], 'npz'))

    energy_kT = np.empty((len(files), len(dna.params)))
    i = 0
    fileio.report_progress(len(files), title='analyze_step_parameters', init=True)
    for f in files:
        dna = HelixPose.from_file(fileio.change_extension(f, 'npz'))
        fileio.report_progress(i + 1)
        j = 0
        for p in dna.params:
            energy_kT[i, j] = (np.sum(0.5 * (p - p0) * np.dot(k, p - p0)))
            j += 1
        i += 1

    F_data = fileio.read_xlsx_collumn(filename, 'F_pN')

    forces = np.linspace(10, 0, 6)
    energy_F = []
    legend = []
    for Fmin, Fmax in zip(forces[1:], forces):
        selected = (F_data > Fmin) & (F_data <= Fmax)
        try:
            energy_F.append(np.mean(energy_kT[selected], axis=0))
            legend.append('F = {:.1f} pN'.format(((Fmin + Fmax) / 2)))
        except:
            pass

    energy_F = np.asarray(energy_F)

    i = range(len(dna.params))

    plt.close()
    plt.figure(figsize=(12, 3))
    # plt.plot(i, energy_kT)
    plt.plot(i, energy_F.T)
    energy_thermal = np.ones(len(dna.params)) * 3
    plt.plot(i, energy_thermal, color='k', linestyle=':', linewidth=0.8)

    plt.xlim(0, len(dna.params))
    plt.ylim(-1, 6)
    plt.ylabel('G (kT)')
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    plt.xlabel('i (bp)')
    plt.legend(legend)
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    plt.pause(5)
    filename = fileio.change_extension(filename, '_Edna.jpg')
    plt.savefig(filename, dpi=600, format='jpg')

    return


def plot_step_params(filename, dataset, save=False, wait=0, plot_energy=True):
    filename = fileio.change_extension(filename, 'xlsx')
    sets, files, _ = fileio.contents_xlsx(filename)
    if dataset == -1:
        filename = sorted(files)[-1]
    else:
        filename = files[sets.index(dataset)]
    dna = HelixPose.from_file(fileio.change_extension(filename, 'npz'))

    p0 = np.load(default_step_file)[0]
    sigma2 = np.load(default_step_file)[1:]
    k = 1 / sigma2
    energy_kT = []
    for p in dna.params:
        energy_kT.append(0.5 * (p - p0) * np.dot(k, p - p0) / kT)
    energy_kT = np.asarray(energy_kT)
    energy_kT = np.sum(np.abs(energy_kT), axis=1)

    i = range(len(dna.params))

    plt.close()
    plt.figure(figsize=(12, 4))
    if plot_energy:
        plt.plot(i, energy_kT)
        plt.ylim(-1, 25)
        plt.ylabel('G (kT)')
    else:
        plt.plot(i, dna.params)
        plt.ylim(-2.5, 4.5)
        plt.ylabel('step parameters')
        plt.legend(['Shift', 'Slide', 'Rise', 'Tilt', 'Roll', 'Twist'], loc=1)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0], loc='left',
              fontdict={'fontsize': 10})
    plt.xlabel('i (bp)')

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    if wait != 0:
        plt.pause(wait)
    if save:
        filename = fileio.change_extension(filename, '_step.jpg')
        plt.savefig(filename, dpi=600, format='jpg')
    return


def plot_fz(filename):
    forces = fileio.read_xlsx_collumn(filename, 'F_pN')
    z = fileio.read_xlsx_collumn(filename, 'z_nm')

    selected = np.diff(np.append([-1], forces)) > 0

    pars, datafile = fileio.read_xlsx_row(filename, 0)

    forces = np.clip(forces, 1e-3, 1e2)
    wlc = 1 - 0.5 * np.sqrt(0.1 * kT / (forces * pars['P_nm'])) / selected
    grid = []
    for i in range(1, pars['n_nuc'] + 1):
        grid.append(wlc * (pars['L_bp'] - 80 * i) / 3)
        grid.append(wlc * (pars['L_bp'] - (pars['n_nuc'] * 80 + i * (147 - 80))) / 3)
    wlc *= pars['L_bp'] / 3

    filename = fileio.change_extension(filename, '_fz.jpg')

    fileio.save_plot((forces, z, wlc, selected), filename=filename,
                     ax_labels=['z (nm)', 'F (pN)'], grid=grid, yrange=[-0.5, 10.5],
                     transpose=True, xrange=[0, 1.1 * pars['L_bp'] / 3])
    return


def main():
    filename = (fileio.get_filename(ext='xlsx', current=True, wildcard='*197*'))
    # filename = 'E:\\users\\noort\\data\\20180516\\1st4x197_002.xlsx'
    # filename = eg.fileopenbox()
    # plot_energy(filename)
    plot_gf(filename)
    plot_fz(filename)
    fileio.create_pov_movie(filename, origin_frame=0, fps=5, reverse=False)
    return


if __name__ == "__main__":
    # execute only if run as a script
    main()
