import matplotlib as mpl

import warnings

warnings.filterwarnings("ignore")

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass
from matplotlib import pyplot as plt
from helixmc.pose import HelixPose
import numpy as np
from lmfit import Parameters
import easygui as eg

# ChromatinMC modules:
import NucleosomeMC as nMC
import FileIO as fileio
import RunMC as rMC
import FiberMC as fMC
import StatPhys as sp

default_step_file = 'C:\\Python27\\Lib\\site-packages\\helixmc\\data\\DNA_gau.npy'
default_folder = 'D:\\users\\'
kT = 41


def plot_gf(filename, calc=False, force_range=[0, 1e9]):
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
            for g_nuc_kT_all, name in zip(g_nuc_kT, names):
                pars[name].value = g_nuc_kT_all
            i += 1
    else:
        names = ['e_nuc_kT', 'e_wrap_kT', 'e_stack_kT']
        for name in names:
            print('{0:12s} = {1:7.1f}'.format(name, fileio.read_xlsx_collumn(filename, name)[0]))

        # e_stack_kT = fileio.read_xlsx_collumn(filename, 'e_stack_kT')

        names = ['g_work_kT', 'g_stack_kT', 'g_wrap_kT', 'g_dna_kT']
        names += ['g_shift_kT', 'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT', 'g_twist_kT']

        g_nuc_kT_all = []
        for name in names:
            g_nuc_kT_all.append(fileio.read_xlsx_collumn(filename, name))

        force = fileio.read_xlsx_collumn(filename, 'F_pN')
        pulling = (np.diff(np.append(0, force), axis=0) > 0)

        print('\naverages for F = {} pN:'.format(force_range))
        for name, g in zip(names, np.asarray(g_nuc_kT_all)):
            selection = g[(force > force_range[0]) * (force < force_range[1]) * (pulling > 0)]
            mean = np.mean(selection)
            sd = np.std(selection)
            print('{0:12s} = {1:7.1f} +/- {2:3.1f} '.format(name, mean, sd))

        plt.close()
        plt.figure(figsize=(4, 3))
        plt.axes([0.15, 0.15, .8, .75])
        colors = ['red', 'green', 'orange', 'steelblue']
        for g, color, name in zip(g_nuc_kT_all, colors, names):
            plt.semilogx(force[pulling], g[pulling], color=color, linewidth=1.5, label=name.split('_')[1])
        for g, color in zip(g_nuc_kT_all, colors):
            plt.semilogx(force[np.logical_not(pulling)], g[np.logical_not(pulling)], color=color, linestyle=':',
                         linewidth=1.5, markersize=15)

        plt.xlim(0.08, 12)
        plt.xlabel('F (pN)')
        plt.ylim(-75, 50)
        plt.ylabel(r'$\Delta$E (kT)')
        plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
        plt.title(filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0] + '\n', loc='left',
                  fontdict={'fontsize': 8})
        plt.legend(fontsize=8, loc='lower center', frameon=False)
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


def entropy_init(n_bps, n_bins=1000):
    params = np.load(default_step_file)
    p0 = params[0]
    sd = np.diag(params[1:])
    range_sd = 10

    bins = []
    for mean, sigma in zip(p0, sd):
        bins.append(np.linspace(mean - range_sd * sigma, mean + range_sd * sigma, n_bins))

    hist = np.zeros((n_bps, 6, len(bins[1]) - 1))
    return hist, np.asarray(bins)


def entropy_add(bp, param, hist, bins):
    for i, p in enumerate(param):
        hist[bp, i, :] += np.histogram([p], bins=bins[i])[0]
    return np.asarray(hist)


def entropy_calc(hist):
    s = np.zeros(len(hist))
    for bp, pars_hist in enumerate(hist):
        for par_hist in pars_hist:
            P = par_hist / par_hist.sum()
            s[bp] += -np.sum(P[P > 0] * np.log(P[P > 0]))
    return s


def test_entropy():
    p = np.random.rand(5, 6)
    p = [0, 1, 3.2, 0, 0, 0.6]
    hist, bins = entropy_init(5)
    hist = entropy_add(0, p, hist, bins)
    s = entropy_calc(hist)
    print(s)

    return


def plot_gi(filename, force_range=[0.15, 1.5], row_nr=0, report_file=None):
    dna_pars = Parameters()

    # print('>>> {}'.format(filename))
    params = np.load(default_step_file)
    p0 = params[0]
    cov = params[1:]
    k = np.linalg.inv(cov)

    sets, filenames, _ = fileio.contents_xlsx(filename, update_dir=True)
    dna = HelixPose.from_file(fileio.change_extension(filenames[0], 'npz'))

    n_nuc = fileio.read_xlsx_collumn(filename, 'n_nuc')[0]
    dyad0 = fileio.read_xlsx_collumn(filename, 'dyad0')[0]
    fiber_start = fileio.read_xlsx_collumn(filename, 'fiber_start')[0]
    e_wrap_kT = fileio.read_xlsx_collumn(filename, 'e_wrap_kT')[0] * -14.0
    e_stack_kT = -fileio.read_xlsx_collumn(filename, 'e_stack_kT')[0]
    NRL = fileio.read_xlsx_collumn(filename, 'NRL')[0]
    nuc_range = NRL * (n_nuc - 1)

    force = fileio.read_xlsx_collumn(filename, 'F_pN')
    pulling = (np.diff(np.append(0, force), axis=0) > 0)
    selection = (force > force_range[0]) * (force <= force_range[1] * (pulling))
    e_wrap_kT += np.mean(fileio.read_xlsx_collumn(filename, 'g_wrap_kT')[selection])
    e_stack_kT += np.mean(fileio.read_xlsx_collumn(filename, 'g_stack_kT')[selection])

    selected_files = []
    for par_file, selected in zip(filenames, selection):
        if selected:
            selected_files.append(par_file)
    if len(selected_files) == 0:
        print('>>> No files in force range')
        return

    hist, bins = entropy_init(len(dna.params))
    energy_kT_all = np.zeros((len(dna.params), 6))

    fileio.report_progress(len(selected_files), title='analyze_step_parameters', init=True)
    for i, par_file in enumerate(selected_files):
        fileio.report_progress(i + 1, title=fileio.change_extension(par_file, 'npz'))
        dna = HelixPose.from_file(fileio.change_extension(par_file, 'npz'))
        for j, p in enumerate(dna.params):
            energy_kT_all[j, :] += (0.5 * (p - p0) * np.dot(k, p - p0))
            hist = np.asarray(entropy_add(j, p, hist, bins))
    energy_kT_all /= (i + 1)
    energy_kT = np.sum(energy_kT_all, axis=1)

    entropy_kT = entropy_calc(hist) / kT

    TS_kT = np.sum(entropy_kT[dyad0:dyad0 + nuc_range]) / (n_nuc - 1)
    E_kT = np.sum(energy_kT[dyad0:dyad0 + nuc_range]) / (n_nuc - 1)
    G_kT = E_kT - TS_kT
    G_kT += e_wrap_kT + e_stack_kT

    energy_per_bp = np.sum(energy_kT_all[dyad0:dyad0 + nuc_range, :], axis=0) / (n_nuc - 1)

    dna_pars.add('n_nuc', value=n_nuc)
    dna_pars.add('NRL', value=NRL)
    dna_pars.add('fiber_start', value=fiber_start)

    names = ['Eshift', 'Eslide', 'Erise', 'Etilt', 'E_roll', 'Etwist']
    for name, value in zip(names, energy_per_bp):
        dna_pars.add(name, value=value)
    dna_pars.add('Edna_kT', value=E_kT)
    dna_pars.add('Ewrap_kT', value=e_wrap_kT)
    dna_pars.add('Estack_kT', value=e_stack_kT)
    dna_pars.add('dE_kT', value=E_kT - (NRL * 3) + e_stack_kT + e_wrap_kT)
    dna_pars.add('TS_kT', value=TS_kT)
    dna_pars.add('G_kT', value=G_kT)
    dna_pars.add('Fmin_pN', value=np.min(force[selection]))
    dna_pars.add('Fmax_pN', value=np.max(force[selection]))

    fileio.write_xlsx_row(filename, row_nr, dna_pars, report_file=report_file)

    result_title = filename.split('\\')[-2] + '\\' + filename.split('\\')[-1].split('.')[0]
    result_title += ', averages for F = {0} pN), per nuc; '.format(force_range)
    result_title += 'Edna = {0:.1f} kT; '.format(E_kT)
    result_title += 'Ewrap = {0:.1f} kT; '.format(e_wrap_kT)
    result_title += 'Estack = {0:.1f} kT; '.format(e_stack_kT)
    result_title += 'TS = {0:.1f} kT; '.format(TS_kT)
    result_title += 'G = {0:.1f} kT; '.format(G_kT)

    i = range(len(dna.params))

    plt.close()
    plt.figure(figsize=(12, 3))
    plt.plot(i, energy_kT, color='b')
    # plt.plot(i, entropy_kT, color='r')

    energy_thermal = np.ones(len(dna.params)) * 3
    plt.plot(i, energy_thermal, color='k', linestyle=':', linewidth=1.2)
    plt.plot(i, energy_thermal * 0, color='k', linestyle=':', linewidth=1.2)

    plt.xlim(0, len(dna.params) + 1)
    plt.ylim(-0.5, 5.5)
    plt.ylabel('G (kT)')
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.title(result_title + '\n', loc='left', fontdict={'fontsize': 8})
    plt.xlabel('i (bp)')
    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.draw()
    plt.pause(5)
    filename = fileio.change_extension(filename, '_Edna.jpg')
    plt.savefig(filename, dpi=600, format='jpg')
    return result_title


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
    if len(selected) == 0:
        return

    pars, datafile = fileio.read_xlsx_row(filename, 0)
    pars['L_bp'].value = 1800
    pars['n_nuc'].value = 8

    forces = np.clip(forces, 1e-3, 1e2)
    # wlc = 1 - 0.5 * np.sqrt(0.1 * kT / (forces * pars['P_nm'])) / selected
    wlc = sp.WLC(forces, 1)[0] / selected
    grid = []
    for i in range(1, pars['n_nuc'] + 1):
        grid.append(wlc * (pars['L_bp'] - 80 * i))
        grid.append(wlc * (pars['L_bp'] - (pars['n_nuc'] * 80 + i * (147 - 80))))
    grid.append(wlc * pars['L_bp'])
    wlc *= pars['L_bp']

    g1_kT = pars['e_stack_kT']
    g1_kT = 25
    e_wrap_kt = pars['e_wrap_kT']
    fiber_start = pars['fiber_start']
    # pars.pretty_print(columns=['value'])

    l2_bp = 86
    if g1_kT == 0:
        l1_bp = 147
    else:
        l1_bp = 120
    g2_kT = e_wrap_kt * (l1_bp - l2_bp - 10) / 10.0
    if e_wrap_kt == 0:
        n_nuc = 0
    else:
        n_nuc = pars['n_nuc'].value
    g2_kT = 5.5

    if pars['NRL'] == 167:
        fib = sp.tether(forces, pars['L_bp'].value, pars['NRL'].value, n_nuc, g1_kT=g1_kT,
                        g2_kT=g2_kT, l1_bp=l1_bp, l2_bp=85, fiber_start=2)
    else:
        fib = sp.tether(forces, pars['L_bp'].value, pars['NRL'].value, n_nuc, g1_kT=g1_kT,
                        g2_kT=g2_kT, l1_bp=l1_bp, l2_bp=85, fiber_start=1)
    # fib = sp.tether(forces, 1000, 167, 4, g1_kT=20,
    #                 g2_kT=5.8, l1_bp=98, l2_bp=85, fiber_start=2)
    # plt.plot(fib, forces)
    # plt.scatter(z, forces)
    # plt.show()
    # return

    filename = fileio.change_extension(filename, '_fz.jpg')
    print(filename)
    fileio.save_plot((forces, z, fib, selected), filename=filename,
                     ax_labels=['z (nm)', 'F (pN)'], grid=grid, yrange=[-0.5, 10.5],
                     transpose=True, xrange=[0, 1.1 * pars['L_bp'] / 3])
    return


def main(filenames):
    report_file = fileio.get_filename(incr=True, root='DNA_analysis', ext='xlsx')
    forces = [0.1, 1.5]
    for row, filename in enumerate(filenames):
        try:
            print('>>> {0}/{1}'.format(row, len(filenames)))
            # plot_gi(filename, force_range=forces, row_nr=row, report_file=report_file)
            # plot_gf(filename, force_range=forces)
            # plot_fz(filename)
            fileio.create_pov_movie(filename, fps=5, octamers=True, overwrite=False, frame=[60, 0, -100])
        except Exception as e:
            print('>>> Error in file', filename, e)
    if report_file is not None:
        print('Results stored in:')
        print(report_file)
    return


if __name__ == "__main__":
    filenames = fileio.get_filename(ext='xlsx', wildcard='8x*_001', date='20180613', list='all')
    # filenames += fileio.get_filename(ext='xlsx', wildcard='8x*', date='20180612', list='all')
    # print(filenames)
    # filenames.append(fileio.get_filename(ext='xlsx', wildcard='*', date='20180612', list='all'))
    # filenames = []
    # filenames.append(eg.fileopenbox(filetypes='*.xlsx'))
    # filenames.append('D:\\Users\\Noort\\data\\20180530\\2st4x167s0w25.xlsx')
    # filenames.append('D:\\Users\\Noort\\data\\20180530\\2st4x197s0w25.xlsx')
    main(filenames)
