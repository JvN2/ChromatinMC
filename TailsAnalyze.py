import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from helixmc.pose import HelixPose
from helixmc import util
import NucleosomeMC as nMC
import FileIO as fileio
import Tails as tMC
import POVenergy as POVe
import RunMC as rMC

def repulsion_exp():

    filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201110 repulsion d varies\20201112_cms_dists.xlsx"
    zero_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,B,C")
    one_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,D,E")
    two_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,F,G")

    zero_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,B,C")
    one_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,D,E")
    two_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,F,G")

    x_1 = range(1, 8)
    x_tick = one_start_167.index[:-1]

    y_0 = zero_start_167.iloc[:, 0]
    e_0 = zero_start_167.iloc[:, 1]

    y_1 = one_start_167.iloc[:, 0]
    e_1 = one_start_167.iloc[:, 1]

    y_2 = two_start_167.iloc[:, 0]
    e_2 = two_start_167.iloc[:, 1]

    y9_0 = zero_start_197.iloc[:, 0]
    e9_0 = zero_start_197.iloc[:, 1]

    y9_1 = one_start_197.iloc[:, 0]
    e9_1 = one_start_197.iloc[:, 1]

    y9_2 = two_start_197.iloc[:, 0]
    e9_2 = two_start_197.iloc[:, 1]

    fig, ax = plt.subplots()
    # left, bottom, width, height = [0.6, 0.55, 0.35, 0.35]
    # ax2 = fig.add_axes([left, bottom, width, height]) # inset
    #
    # ax.errorbar(x_tick, y_1.iloc[:-1], e_1.iloc[:-1], color=(0.75, 0, 0.25), marker='^', markersize=5, label='167 1-start', linewidth=0,
    #             ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    # ax.axhline(y=y_1['old'], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)
    ax.errorbar(x_tick, y_2.iloc[:-1], e_2.iloc[:-1], color=(0.75, 0, 0.25), marker='s', markersize=5, label='167 2-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y_2['old'], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)

    # ax.errorbar(x_tick, y_0.iloc[:-1], e_0.iloc[:-1], color=(0.75, 0, 0.25), marker='o', markersize=5, label='167 0-start', linewidth=0,
    #             ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    # ax.axhline(y=y_0['old'], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)

    # ax.errorbar(x_tick, y9_1.iloc[:-1], e9_1.iloc[:-1], color=(0, 0.75, 0.25), marker='^', markersize=5, label='197 1-start', linewidth=0,
    #             ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    # ax.axhline(y=y9_1['old'], color=(0.6, 0.6, 0.6), linestyle='--', lw=2)
    ax.errorbar(x_tick, y9_2.iloc[:-1], e9_2.iloc[:-1], color=(0, 0.75, 0.25), marker='s', markersize=5, label='197 2-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y9_2['old'], color=(0.6, 0.6, 0.6), linestyle='--', lw=2)

    # ax.errorbar(x_tick, y9_0.iloc[:-1], e9_0.iloc[:-1], color=(0, 0.75, 0.25), marker='o', markersize=5, label='197 0-start', linewidth=0,
    #             ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    # ax.axhline(y=y9_0['old'], color=(0.6, 0.6, 0.6), linestyle='--', lw=2)

    labels=['167 1-start fps', '167 2-start fps', '197 1-start fps', '197 2-start fps', '167 1-start', '167 2-start', '197 1-start', '197 2-start',]
    save_loc = fileio.change_extension(filename, '2-start.png')
    format_plot('decay length (nm)', 'distance (nm)', 'title', scale_page=(1.0/3),
                aspect=1, save=save_loc, yrange=[0,33], legend=None, ax=ax)

    return


def expo_decay ():

    kT = 41.0
    Amp = 2.5 * kT # amplitude pNA
    # Amp = [20, 61, 102, 143, 184, 205] # amplitude pNA
    decay_l = [5.0, 25.0, 45.0, 65.0, 85.0, 100.0]
    # decay_l = 50.0

    x = np.linspace(0,100, 500)
    y = {}
    for d in decay_l:
        y[d] = Amp * np.exp( - (1 / d) * x)
        y[d] /= kT
    # for A in Amp:
    #     y[A] = A * np.exp(- (1 / decay_l) * x)
    #     y[A] /= kT

    x /= 10 # A to nm
    #
    # plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()
    label = []
    for d in decay_l:
        ax.plot(x, y[d], color=(d / 100, 0, 1 - (d / 100)), linewidth=4, label='$\lambda$' + ' = ' + str(d/10))
        label.append('$\lambda$' + ' = ' + str(d/10))

    # for A in Amp:
    #     ax.plot(x, y[A], color=(A / 205.0, 0, 1 - (A / 205.0)), linewidth=5, label='A' + ' = ' + "{:3.1f}".format(A/kT))
    #
    sv_loc = r"D:\users\vlaar\data\repulsion_graphA2_5kTdvaries.png"
    format_plot('distance (nm)', 'energy (kT)', 'title', scale_page=(1.0/2.0),
                aspect=1, save=sv_loc, yrange=None, legend=label, ax=ax)

    return


def debye ():

    Lb = 0.71 #Bjerrum length

    x = np.linspace(0.00055,0.25, 500)
    y = 1/(np.sqrt(8*np.pi*Lb*x))

    fig, ax = plt.subplots()

    ax.plot(x, y, color=(0, 0, 0.75), linewidth=4)

    sv_loc = r"D:\users\vlaar\data\Debeylength_graph.png"
    format_plot('concentration (M)', 'Screening length (nm)', 'title', scale_page=(1.0/2.0),
                aspect=1, save=sv_loc, yrange=None, ax=ax)

    return


def coth(x):
    return np.cosh(x) / np.sinh(x)


def Langevin(x):
    return (coth(x) - 1.0 / x)


def tail_energy():

    # initial parameters
    L_nm = 6.8  # contour length H4 tial
    b_nm = 0.6  # kuhn length H4 tail
    S_pN = 630  # stiffness H4 tail
    kT = 4.10
    # possible force on H4 tials
    f_array = np.linspace(1e-4, 80, 1e6)
    # corresponding extension of H4 tials
    z_array = L_nm * (Langevin(b_nm * f_array / kT) + f_array / S_pN)
    # corresponding energy of H4 tials
    g_array = -(kT * L_nm / b_nm) * (np.log((b_nm * f_array) / (kT)) - np.log(
        np.sinh((b_nm * f_array) / (kT)))) + L_nm * f_array ** 2 / (2 * S_pN)

    g_array /= kT


    fig, ax = plt.subplots()

    ax.plot(z_array, g_array, color=(0.8, 0, 0.8), linewidth=5)
    sv_loc = r"C:\Users\Annelies\Documents\Natuurkunde\Master project 02\figures\Tail_extension_energy_graph.png"
    format_plot('distance (nm)', 'energy (kT)', 'titel', scale_page=(3.0/5.0),
                aspect=0.7, save=sv_loc, ax=ax)

    return


def get_nucl_dyads(L_bp, NRL, n_nucs):
    """

    Returns
    -------
    nucl: nucleosome pose
    dyads: indices of dyads
    """
    # create nucleosomepose
    nucl = nMC.NucPose()
    nucl.from_file('1KX5')
    # get list of dyad indices
    dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
    dyads = (dyads + L_bp // 2).astype(int)

    return nucl, dyads


def get_stack_params(filename):
    # get list of xlsx files in filename folder
    xlsx_f = glob.glob(fileio.change_extension(filename, '\*.xlsx'))
    xlsx = []
    nucl_stack_params = []
    nucl_stack_params_std = []
    for i, file in enumerate(xlsx_f):
        xlsx = pd.read_excel(file, usecols =['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'])
        nucl_stack_params.append(xlsx.iloc[1:-1].mean(axis=0))
        nucl_stack_params_std.append(xlsx.iloc[1:-1])

    nucl_stack_params = np.mean(nucl_stack_params, axis=0)
    nucl_stack_params_std = np.std(np.vstack(nucl_stack_params_std), axis=0)

    return nucl_stack_params, nucl_stack_params_std


def format_plot(xtitle='x (a.u.)', ytitle='y (a.u.)', title='', xrange=None, yrange=None,
                ylog=False, xlog=False, scale_page=1.0, aspect=0.5, save=None, boxed=True,
                GUI=False, ref='', legend=None, fig=None, ax=None):
    # adjust the format to nice looks
    from matplotlib.ticker import AutoMinorLocator
    import os, subprocess

    page_width = 7  # inches ( = A4 width - 1 inch margin)
    margins = (0.55, 0.45, 0.2, 0.2)  # inches
    fontsize = 14

    if fig == None:
        fig = plt.gcf()
    if ax == None:
        ax = plt.gca()

    # Set up figure
    fig_width = page_width * scale_page
    fig_height = (fig_width - (margins[0] + margins[2])) * aspect + margins[1] + margins[3]
    fig.set_size_inches(fig_width, fig_height)

    # Set up axes
    ax_rect = [margins[0] / fig_width]
    ax_rect.append(margins[1] / fig_height)
    ax_rect.append((fig_width - margins[0] - margins[2]) / fig_width)
    ax_rect.append((fig_height - margins[1] - margins[3]) / fig_height)

    # ax_rect = [margins[0] / fig_width,
    #            margins[1] / fig_height,
    #            (fig_width - margins[2] - margins[0]) / fig_width,
    #            (fig_height - margins[3] - margins[1]) / fig_height
    #            ]

    ax.set_position(ax_rect)

    # Add axis titles and frame label; use absolute locations, rather then leaving it up to Matplotlib
    if ref is not None:
        plt.text(ax_rect[1] * 0.15, ax_rect[-1] + ax_rect[1], ref, horizontalalignment='left',
                 verticalalignment='center',
                 fontsize=fontsize * 1.2, transform=fig.transFigure)
    plt.text(ax_rect[0] + 0.5 * ax_rect[2], 0, xtitle, horizontalalignment='center',
             verticalalignment='bottom', fontsize=fontsize, transform=fig.transFigure)
    plt.text(ax_rect[1] * 0.005, ax_rect[1] + 0.5 * ax_rect[3], ytitle, horizontalalignment='left',
             verticalalignment='center', fontsize=fontsize, transform=fig.transFigure, rotation=90)

    if legend is not None:
        plt.rcParams["legend.frameon"] = False
        plt.rcParams["legend.edgecolor"] = 'none'
        plt.rcParams["legend.labelspacing"] = 0.25
        plt.rcParams["legend.handlelength"] = 1
        plt.rcParams["legend.handletextpad"] = 0.25
        plt.legend(legend, prop={'size': fontsize * 0.8}, markerscale=100)

    # fig.canvas.mpl_connect("key_press_event", key_press_action)

    if not boxed:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', axis='both', bottom=True, top=boxed, left=True, right=boxed, direction='in')
    ax.tick_params(which='major', length=4, labelsize=fontsize * 0.8, width=1)
    ax.tick_params(which='minor', length=2, width=1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1)

    if xlog:
        ax.semilogx()
    else:
        ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, 3))
    if ylog:
        ax.semilogy()
    else:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 3))

    if xrange is not None:
        ax.set_xlim(xrange)
    if yrange is not None:
        ax.set_ylim(yrange)

    if not GUI and save is None: plt.show()

    if save is not None:
        if not os.path.exists(os.path.dirname(save)):
            os.makedirs(os.path.dirname(save))
        base, ext = os.path.splitext(save)
        if ext == '.emf':
            save = base + '.pdf'
        fig.savefig(save, dpi=1200, transparent=False)
        if ext == '.emf':
            try:
                subprocess.call(["C:\Program Files\Inkscape\inkscape.exe", "--file", save, "--export-emf",
                                 base + '.emf'])
                os.remove(save)
            except:
                print('Install Inkscape for conversion to emf.\nPlot saved as pdf in stead.')
        plt.close()
    return fig, ax

def stack_exp(param_name):

    filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201231 Nucleosome stacking.xlsx"

    fixed = pd.read_excel(filename, sheet_name='fixed', header=0, index_col=0, usecols="A,B,C,D,E,F,G")
    NRL167 = pd.read_excel(filename, sheet_name='167', header=0, index_col=0)
    NRL197 = pd.read_excel(filename, sheet_name='197', header=0, index_col=0)
    #
    fig, ax = plt.subplots()
    # ax.get_xaxis().set_visible(False)
    ax.set_xticks([5.5,15.5])
    ax.set_xticklabels(['167', '197'])
    # left, bottom, width, height = [0.6, 0.55, 0.35, 0.35]
    # ax2 = fig.add_axes([left, bottom, width, height]) # inset
    ax.axhline(y=fixed.loc['1-start', param_name], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)
    # ax.axhline(y=fixed.loc['2-start', param_name], color=(0.6, 0.6, 0.6), linestyle='--', lw=2)

    # FPS 167
    ax.errorbar(1, NRL167.loc['fps 0', param_name], NRL167.loc['fps 0', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='o', markersize=5, label='fps 0-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(2, NRL167.loc['fps 1', param_name], NRL167.loc['fps 1', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='^', markersize=5, label='fps 1-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(3, NRL167.loc['fps 2', param_name], NRL167.loc['fps 2', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='s', markersize=5, label='fps 2-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)

    # 1.31 nm 167
    ax.errorbar(5, NRL167.loc['1.31 0', param_name], NRL167.loc['1.31 0', param_name + ' std'], color=(0.75, 0, 0.25), marker='o', markersize=5, label='1.31 0-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(7, NRL167.loc['1.31 1', param_name], NRL167.loc['1.31 1', param_name + ' std'], color=(0.75, 0, 0.25), marker='^', markersize=5, label='1.31 1-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(9, NRL167.loc['1.31 2', param_name], NRL167.loc['1.31 2', param_name + ' std'], color=(0.75, 0, 0.25), marker='s', markersize=5, label='1.31 2-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)

    # 7.91 nm 167
    ax.errorbar(6, NRL167.loc['7.91 0', param_name], NRL167.loc['7.91 0', param_name + ' std'], mec=(0.75, 0, 0.25), marker='o', markersize=5, mfc=(1, 1, 1), label='7.91 0-start',
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(8, NRL167.loc['7.91 1', param_name], NRL167.loc['7.91 1', param_name + ' std'], mec=(0.75, 0, 0.25), marker='^', markersize=5, mfc=(1, 1, 1), label='7.91 1-start',
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(10, NRL167.loc['7.91 2', param_name], NRL167.loc['7.91 2', param_name + ' std'], mec=(0.75, 0, 0.25), marker='s', markersize=5, mfc=(1, 1, 1), label='7.91 2-start',
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)

    # # FPS 197
    ax.errorbar(12, NRL197.loc['fps 0', param_name], NRL197.loc['fps 0', param_name + ' std'], color=(0.2, 0.2, 0.2),
                marker='o', markersize=5, label='fps 0-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(13, NRL197.loc['fps 1', param_name], NRL197.loc['fps 1', param_name + ' std'], color=(0.2, 0.2, 0.2),
                marker='^', markersize=5, label='fps 1-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(14, NRL197.loc['fps 2', param_name], NRL197.loc['fps 2', param_name + ' std'], color=(0.2, 0.2, 0.2),
                marker='s', markersize=5, label='fps 2-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    #
    # # 1.31 nm 197
    ax.errorbar(16, NRL197.loc['1.31 0', param_name], NRL197.loc['1.31 0', param_name + ' std'], color=(0, 0.75, 0.25),
                marker='o', markersize=5, label='1.31 0-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(18, NRL197.loc['1.31 1', param_name], NRL197.loc['1.31 1', param_name + ' std'], color=(0, 0.75, 0.25),
                marker='^', markersize=5, label='1.31 1-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(20, NRL197.loc['1.31 2', param_name], NRL197.loc['1.31 2', param_name + ' std'], color=(0, 0.75, 0.25),
                marker='s', markersize=5, label='1.31 2-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    #
    # # 7.91 nm 197
    ax.errorbar(17, NRL197.loc['7.91 0', param_name], NRL197.loc['7.91 0', param_name + ' std'], mec=(0, 0.75, 0.25),
                marker='o', markersize=5, mfc=(1, 1, 1), label='7.91 0-start',
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(19, NRL197.loc['7.91 1', param_name], NRL197.loc['7.91 1', param_name + ' std'], mec=(0, 0.75, 0.25),
                marker='^', markersize=5, mfc=(1, 1, 1), label='7.91 1-start',
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(21, NRL197.loc['7.91 2', param_name], NRL197.loc['7.91 2', param_name + ' std'], mec=(0, 0.75, 0.25),
                marker='s', markersize=5, mfc=(1, 1, 1), label='7.91 2-start',
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    #
    # labels=['167 1-start fps', '167 2-start fps', '197 1-start fps', '197 2-start fps', '167 1-start', '167 2-start', '197 1-start', '197 2-start',]
    save_loc = fileio.change_extension(filename, (param_name + '.png'))
    format_plot(' ', param_name + ' ('u'\xb0'')', 'title', scale_page=(1.0/3),
    # format_plot(' ', param_name, 'title', scale_page=(1.0/3),
                aspect=1, save=save_loc, yrange=None, legend=None, ax=ax)
    return



def get_stack_params(filename):

    # get list of xlsx files in filename folder
    xlsx_f = glob.glob(fileio.change_extension(filename, '\*.xlsx'))

    params = []
    for f in xlsx_f:
        xlsx = pd.read_excel(f, header=0, index_col=0)
        params.append(xlsx.iloc[1:-1])

    df_stack = pd.concat(params, axis=0)
    df_stack['shift (A)'] = df_stack['shift (A)']/10
    df_stack['slide (A)'] = df_stack['slide (A)']/10
    df_stack['rise (A)'] = df_stack['rise (A)']/10
    df_stack['tilt'] = df_stack['tilt'] * 180 / np.pi
    df_stack['roll'] = df_stack['roll'] * 180 / np.pi
    df_stack['twist'] = df_stack['twist'] * 180 / np.pi

    df_stack.to_excel(fileio.change_extension(filename, '_stack.xlsx'))
    print(df_stack)

def read_npz(filename):

    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))

    # Get parameters from excel sheet
    pars = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0)

    # save parameters of bp of every dna pose in params
    params = []
    frames = []
    coords = []

    for f in npz_f:
        dna = HelixPose.from_file(f)
        params.append(dna.params)
        frames.append(dna.frames)
        coords.append(dna.coord)

    return params, frames, coords, pars

def g_dna_kT(filename):

    params, frames, coords, pars = read_npz(filename)

    e_wrap_kT = pars.iloc[0]['e_wrap_kT']
    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']

    # Initialize random steps
    p0 = np.load(util.locate_data_file('DNA_gau.npy'))[0]
    k = np.linalg.inv(np.load(util.locate_data_file('DNA_gau.npy'))[1:])

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)
    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coord, nucl.dna.frames, nucl.dyad, nucl.fixed)

    g_dna = np.zeros(6)
    g_dna_kT = []

    for i, p in enumerate(params):
        # create binary list of fixed and free bp's
        w = np.ones(len(p))
        for dyad in dyads:
            fixed_bps = rMC.score_wrapping(dyad + 1, coords[i], frames[i], dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                           half_nuc=True)[1]
            if len(fixed_bps) > 0:
                w[fixed_bps[0]:fixed_bps[1]] = 0

            # else:
            #     w = None
        g_dna_all = rMC.score_dna(p, p0, k, w=w)

        for dyad1, dyad2 in zip(dyads[:-1], dyads[1:]):
            g_dna += np.sum(g_dna_all[dyad1:dyad2], axis=0)


        g_dna /= (nucs - 1)

        g_dna_kT.append(np.sum(g_dna))

    print('DNA')
    print(np.mean(g_dna_kT))
    print(np.std(g_dna_kT))
    return g_dna_kT

def tabel(filename):

    g_dna = g_dna_kT(filename)
    xlsx = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0, index_col=0)

    print('stack (tails)')
    print(xlsx.iloc[:-1]['g_stack_kT'].mean())
    print(xlsx.iloc[:-1]['g_stack_kT'].std())
    print('tails')
    print(xlsx.iloc[:-1]['g_tails_kT'].mean())
    print(xlsx.iloc[:-1]['g_tails_kT'].std())
    print('rep')
    print(xlsx.iloc[:-1]['g_rep_kT'].mean())
    print(xlsx.iloc[:-1]['g_rep_kT'].std())
    print('wrap')
    print(xlsx.iloc[:-1]['g_wrap_kT'].mean())
    print(xlsx.iloc[:-1]['g_wrap_kT'].std())

    df = np.sum([g_dna, xlsx.iloc[:-1]['g_stack_kT'], xlsx.iloc[:-1]['g_tails_kT'], xlsx.iloc[:-1]['g_rep_kT'], xlsx.iloc[:-1]['g_wrap_kT']], axis=0)
    print('Total')
    if xlsx.iloc[1]['fiber_start'] == 0:
        print('0-start')
        print(np.mean(df) - 0.0)
    else:
        print(np.mean(df) - 25.0)
    print(np.std(df))

def wrap(filename):

    xlsx = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0, index_col=0)
    print(xlsx.iloc[1]['NRL'])
    print(xlsx.iloc[1]['fiber_start'])
    print(xlsx.iloc[1]['e_wrap_kT'])
    print('wrap')
    print(xlsx.iloc[:-1]['g_wrap_kT'].mean())
    print(xlsx.iloc[:-1]['g_wrap_kT'].std())
    return

def plot_wrap(NRL, fiberstart):

    filename = r"D:\Downloads\2021-01-02 E unwrap exp.xlsx"

    data = pd.read_excel(filename, sheet_name=str(NRL), header=0, index_col=0)
    x_tick = data.columns
    y_fps = data.loc['FPS ' + str(fiberstart) + '-start']
    e_fps = data.loc['FPS ' + str(fiberstart) + '-start' + ' std']

    y_131 = data.loc['131 ' + str(fiberstart) + '-start']
    e_131 = data.loc['131 ' + str(fiberstart) + '-start' + ' std']

    y_791 = data.loc['791 ' + str(fiberstart) + '-start']
    e_791 = data.loc['791 ' + str(fiberstart) + '-start' + ' std']


    if fiberstart == 1:
        marker = '^'
    elif fiberstart == 2:
        marker = 's'
    else:
        marker = 'o'

    if NRL == 167:
        color = (0.75, 0, 0.25)
    else:
        color = (0, 0.75, 0.25)


    fig, ax = plt.subplots()


    # FPS
    ax.errorbar(x_tick, y_fps, e_fps, color=(0.2, 0.2, 0.2),
                marker=marker, markersize=5, label='fps', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=1, capsize=3)

    ax.errorbar(x_tick, y_131, e_131, color=color,
                marker=marker, markersize=5, label='1.31', linewidth=0,
                ecolor=color, elinewidth=1, capsize=3)

    ax.errorbar(x_tick, y_791, e_791, mec=color,
                marker=marker, markersize=5, mfc=(1, 1, 1), label='7.91', linewidth=0,
                ecolor=color, elinewidth=1, capsize=3)

    save_loc = fileio.change_extension(filename, (str(NRL) + 's' + str(fiberstart) + '.png'))
    format_plot('E$_{unwrap, max}$ (K$_B$T)', 'E$_{unwrap}$ (K$_B$T)', 'title', scale_page=(1.0/3),
                aspect=1, save=save_loc, yrange=[-1.0,11], legend=None, ax=ax)
    return

def get_g_dna(filename):

    params, frames, coords, pars = read_npz(filename)
    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']
    e_wrap_kT = pars.iloc[0]['e_wrap_kT']

    # Initialize random steps
    p0 = np.load(util.locate_data_file('DNA_gau.npy'))[0]
    k = np.linalg.inv(np.load(util.locate_data_file('DNA_gau.npy'))[1:])

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)
    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coord, nucl.dna.frames, nucl.dyad, nucl.fixed)

    # get energies per bp in fiberpose
    g_dna = []

    for i, p in enumerate(params):
        # create binary list of fixed and free bp's
        w = np.ones(len(p))
        for dyad in dyads:
            fixed_bps = rMC.score_wrapping(dyad + 1, coords[i], frames[i], dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                           half_nuc=True)[1]
            if len(fixed_bps) > 0:
                w[fixed_bps[0]:fixed_bps[1]] = 0
            # else:
            #     w = None
        g_dna.append(rMC.score_dna(p, p0, k, w=w)) #npzsxL_bpX6

    g_dna_all = np.mean(g_dna, axis=0) #L_bpx6
    g_dna_std = np.std(np.sum(g_dna, axis=2), axis=0)
    g_dna_m = np.mean(np.sum(g_dna, axis=2), axis=0)

    df_std = pd.DataFrame(g_dna_std, columns=['g_total std'], index=range(1, L_bp - 1))
    df_m = pd.DataFrame(g_dna_m, columns=['g_total (kT)'], index=range(1, L_bp - 1))

    df_all = pd.DataFrame(g_dna_all, columns=['g_shift_kT',
                                                  'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT',
                                                  'g_twist_kT'], index=range(1, L_bp - 1))

    df_g_dna = pd.concat([df_all, df_m, df_std], axis=1)
    df_g_dna.to_excel(fileio.change_extension(filename, '_g_dna.xlsx'))

    begin = dyads[3] - 75
    end = dyads[5] + 75

    print(np.max(g_dna_m[begin:end]))

    return df_m


def get_mean_coords(filename):

    params = read_npz(filename)[0]
    pars = read_npz(filename)[3]
    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)

    # calculate mean parameters per bp
    params = np.mean(params, axis=0)

    # use 6 parameters to get coordinates of every basepair
    dr, frames = util.params2data(params)
    coords = util.dr2coords(dr)

    # get dyad_ofs to project histones in mean fiber pose
    of_d_fiber = []  # origin frame of dyad in fiber
    of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)  # origin frame dyad in nucl pose
    tf_d = []  # transformation matrix

    for i, d in enumerate(dyads):
        # get origin frame of dyad in fiber
        of_d_fiber.append(nMC.get_of_2(coords, frames, d))

        # get transformation matrix of nucleosome dyad onto fiber dyad
        tf_d.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

    # get center of masses of nucleosome
    nuc_cms = []
    for d, dyad in enumerate(dyads):
        nuc_cms.append(nMC.apply_transformation(nucl.of, tf_d[d]))

    # append histone positions to coordinates
    # tf_d ensures that histones are placed correct at nucleosome position
    coord_w_hist, radius, colors = tMC.get_histones(coords, dyads, nucl, tf=tf_d, tail=False)

    return dyads, nuc_cms, coord_w_hist, radius, colors


def dna_energy_display(filename, energy_kT='g_total (kT)'):

    file = get_g_dna(filename)
    # return

    E_max = 1.5 #kT
    file[energy_kT] = np.clip(file[energy_kT], 0 , E_max)

    colorwaaier = []
    rod = []
    for i, b in enumerate(np.arange(50)):
        p = b/50.0
        colorwaaier.append([0.6 + (0.4 * p) , 0.6 * (1 - p),  0.6 * (1 - p)])
        rod.append([0,0,i/10.0])


    # POVe.main(fileio.change_extension(filename, 'label.png'), rod, colorwaaier, radius=1, range_A=[25, 25], offset_A=[0, 0, 10], width_pix=500, showt=True)

    colors = []
    for e in file[energy_kT]:
        c = e/E_max
        colors.append([0.6 + (0.4 * c) , 0.6 * (1 - c),  0.6 * (1 - c)])


    dyads, nuc_cms, coords = get_mean_coords(filename)[0:3]

    # transform fiber to origin
    # origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]]) np.pi/6.
    # origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]]) 0.25*np.pi
    angle = 0.7*np.pi
    cos = np.cos(angle)
    sin = np.sin(angle)
    origin_of = np.asarray([[0, 0, 0], [-cos, -sin, 0], [-sin, cos, 0], [0, 0, -1]])
    tf_o = nMC.get_transformation(nuc_cms[3], target=origin_of)
    # Tranform coords where first nucleosome is placed in origin
    t_coord = nMC.apply_transformation(coords[0], tf_o)

    begin = dyads[3] - 75
    end = dyads[5] + 75

    POVe.main(fileio.change_extension(filename, 'Etot.png'), t_coord[begin:end], colors[begin:end], radius=10, range_A=[500, 500], offset_A=[0, 0, 125], width_pix=500, showt=True)

    return


def plot_npz(filename):

    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))

    # Get parameters from excel sheet
    pars = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0)

    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']
    # L_bp = 1364
    # nrl = 167
    # nucs = 8

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)

    for f in npz_f:
        dna = HelixPose.from_file(f)

        # use 6 parameters to get coordinates of every basepair
        coords = dna.coord
        frames = dna.frames

        # get dyad_ofs to project histones in mean fiber pose
        of_d_fiber = []  # origin frame of dyad in fiber
        of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)  # origin frame dyad in nucl pose
        tf_d = []  # transformation matrix

        for i, d in enumerate(dyads):
            # get origin frame of dyad in fiber
            of_d_fiber.append(nMC.get_of_2(coords, frames, d))

            # get transformation matrix of nucleosome dyad onto fiber dyad
            tf_d.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

        # get center of masses of nucleosome
        nuc_cms = []
        for d, dyad in enumerate(dyads):
            nuc_cms.append(nMC.apply_transformation(nucl.of, tf_d[d]))

        # append histone positions to coordinates
        # tf_d ensures that histones are placed correct at nucleosome position
        coord_w_hist, radius, colors = tMC.get_histones(coords, dyads, nucl, tf=tf_d, tail=False)
        # coord_w_hist, radius, colors = tMC.get_histones(coords[dyads[3] - 75:dyads[5] + 75], dyads, nucl, tf=tf_d[3:6], tail=False)

        # transform fiber to origin
        origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]])
        # origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]])
        tf_o = nMC.get_transformation(nuc_cms[4], target=origin_of)
        # origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]]) np.pi/6.
        # origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]]) 0.25*np.pi
        angle = 1.3 * np.pi
        cos = np.cos(angle)
        sin = np.sin(angle)
        origin_of = np.asarray([[0, 0, 0], [-cos, -sin, 0], [-sin, cos, 0], [0, 0, -1]])
        tf_o = nMC.get_transformation(nuc_cms[3], target=origin_of)
        t_coord = []  # transformed coords
        # Tranform coords where first nucleosome is placed in origin
        for c in coord_w_hist:
            t_coord.append(nMC.apply_transformation(c, tf_o))


        print(fileio.create_pov((fileio.change_extension(f, 'png')), t_coord, radius=radius, colors=colors, range_A=[1000, 1000],
                                    offset_A=[0, 0, 300], show=False, width_pix=1500))

    return


def de_grote_chromatine_show(filename, size):

    params = read_npz(filename)[0]
    pars = read_npz(filename)[3]
    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']
    fiber_start = pars.iloc[0]['fiber_start']

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)

    # calculate mean parameters per bp
    params = np.mean(params, axis=0)

    dyads_block_idx = np.arange(dyads[fiber_start], dyads[-2], fiber_start*nrl)
    params_block = []
    for idx in dyads_block_idx:
        params_block.append(params[idx:idx+fiber_start*nrl])

    params_block = np.mean(params_block, axis=0)
    params_full = []
    dyads_new = []
    for i in range(0,size):
        params_full.extend(params_block)
        if fiber_start == 2:
            dyads_new.append(2*i*nrl)
            dyads_new.append((2*i+1)*nrl)
        else:
            dyads_new.append(i*nrl)

    # delete dyad at position zero (the edge)
    dyads_new.pop(0)
    # use 6 parameters to get coordinates of every basepair
    dr, frames = util.params2data(np.asarray(params_full))
    coords = util.dr2coords(dr)

    # get dyad_ofs to project histones in mean fiber pose
    of_d_fiber = []  # origin frame of dyad in fiber
    of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)  # origin frame dyad in nucl pose
    tf_d = []  # transformation matrix

    for i, d in enumerate(dyads_new):
        # get origin frame of dyad in fiber
        of_d_fiber.append(nMC.get_of_2(coords, frames, d))

        # get transformation matrix of nucleosome dyad onto fiber dyad
        tf_d.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

    # append histone positions to coordinates
    # tf_d ensures that histones are placed correct at nucleosome position
    coord_w_hist, radius, colors = tMC.get_histones(coords[100:-100], dyads_new, nucl, tf=tf_d, tail=False)

    print(fileio.create_pov((fileio.change_extension(filename, '_2big.png')), coord_w_hist, radius=radius, colors=colors,
                            range_A=[1500, 1500],
                            offset_A=[0, 0, 500], show=False, width_pix=1500))

    return

def get_g_linker(filename):

    params, frames, coords, pars = read_npz(filename)
    L_bp = pars.iloc[0]['L_bp']
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']
    e_wrap_kT = pars.iloc[0]['e_wrap_kT']

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(L_bp, nrl, nucs)
    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coord, nucl.dna.frames, nucl.dyad, nucl.fixed)

    # Initialize random steps
    p0 = np.load(util.locate_data_file('DNA_gau.npy'))[0]
    k = np.linalg.inv(np.load(util.locate_data_file('DNA_gau.npy'))[1:])

    # get energies per bp in fiberpose
    g_dna = []

    for i, p in enumerate(params):
        # create binary list of fixed and free bp's
        w = np.ones(len(p))
        for dyad in dyads:
            fixed_bps = rMC.score_wrapping(dyad + 1, coords[i], frames[i], dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                           half_nuc=True)[1]
            if len(fixed_bps) > 0:
                w[fixed_bps[0]:fixed_bps[1]] = 0
            # else:
            #     w = None
        g_dna.append(rMC.score_dna(p, p0, k, w=w))  # 10xL_bpX6


    dyads_block_idx = np.arange(dyads[1], dyads[-2], nrl)
    g_dna_block = []
    for g, npz in enumerate(g_dna): # energy of each bp in each npz file
        for idx in dyads_block_idx:
            g_dna_block.append(g_dna[g][idx:idx + nrl])


    g_dna_std = np.std(g_dna_block, axis=0)
    g_dna_all = np.mean(g_dna_block, axis=0)
    g_dna_m = np.sum(g_dna_all, axis=1)

    df_std = pd.DataFrame(g_dna_std, columns=['g_shift_std',
                                              'g_slide_std', 'g_rise_std', 'g_tilt_std', 'g_roll_std',
                                              'g_twist_std'], index=range(nrl))

    df_m = pd.DataFrame(g_dna_m, columns=['g_total (kT)'], index=range(nrl))

    df_all = pd.DataFrame(g_dna_all, columns=['g_shift_kT',
                                              'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT',
                                              'g_twist_kT'], index=range(nrl))

    df_g_dna = pd.concat([df_all, df_m, df_std], axis=1)
    df_g_dna.to_excel(fileio.change_extension(filename, 'g_linker.xlsx'))

    return df_g_dna

def plot_g_linker(filename_1, filename_2):

    # param_name = 'g_shift_kT'
    # error_name = 'g_shift_std'
    # param_name = 'g_slide_kT'
    # error_name = 'g_slide_std'
    # param_name = 'g_rise_kT'
    # error_name = 'g_rise_std'
    # param_name = 'g_tilt_kT'
    # error_name = 'g_tilt_std'
    # param_name = 'g_roll_kT'
    # error_name = 'g_roll_std'
    param_name = 'g_twist_kT'
    error_name = 'g_twist_std'
    # df_g_dna_1 = get_g_linker(filename_1)
    # df_g_dna_2 = get_g_linker(filename_2)
    df_g_dna_1 = pd.read_excel(filename_1 + 'g_linker.xlsx', header=0, index_col=0)
    df_g_dna_2 = pd.read_excel(filename_2 + 'g_linker.xlsx', header=0, index_col=0)

    # 167
    # x_tick = df_g_dna_1.index[55:-55]
    # y_1 = df_g_dna_1.iloc[55:-55][param_name]
    # e_1 = df_g_dna_1.iloc[55:-55][error_name]
    # y_2 = df_g_dna_2.iloc[55:-55][param_name]
    # e_2 = df_g_dna_2.iloc[55:-55][error_name]
    # 197
    x_tick = df_g_dna_1.index[50:-50]
    y_1 = df_g_dna_1.iloc[50:-50][param_name]
    e_1 = df_g_dna_1.iloc[50:-50][error_name]
    y_2 = df_g_dna_2.iloc[50:-50][param_name]
    e_2 = df_g_dna_2.iloc[50:-50][error_name]

    fig, ax = plt.subplots()
    #
    ax.errorbar(x_tick, y_1, e_1, color=(1, 0, 0), marker='^', markersize=1, label='167 1-start', linewidth=0,
                ecolor=(1, 0, 0), elinewidth=0.0, capsize=0.0)
    # ax.axhline(y=y_1['old'], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)
    ax.errorbar(x_tick, y_2, e_2, color=(0, 0, 1), marker='s', markersize=1,
                label='167 2-start', linewidth=0,
                ecolor=(0, 0, 1), elinewidth=0.0, capsize=0.0)
    # ax.axhline(y=y_2['old'], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)

    save_loc = fileio.change_extension(filename_1, (param_name + '.png'))
    format_plot('basepair', '$\\Delta$G (k$_B$T)' , 'title', scale_page=(1.0/4.0),
                aspect=1, save=save_loc, yrange=[-1.0,1.0], legend=None, ax=ax)
    return

def plot_tail(filename, filename_2=None):

    tail = pd.read_excel(filename, header=0, index_col=0)
    x_tick = tail.index[:]
    t_up = tail.loc[:, 'tail up (nm)']
    t_down = tail.loc[:, 'tail down (nm)']

    fig, ax = plt.subplots()
    ax.plot(x_tick, t_up, color=(0.75, 0, 0.5), marker='^', label='up', markersize=0.05, linestyle='')
    ax.plot(x_tick, t_down, color=(1, 0, 1), marker='v', label='down', markersize=0.05, linestyle='')
    # print(np.mean(t_up))
    # print(np.std(t_up))
    # print(np.mean(t_down))
    # print(np.std(t_down))
    # return

    label = ['up', 'down']

    if filename_2 is not None:
        tail_2 = pd.read_excel(filename_2, header=0, index_col=0)
        t_up_2 = tail_2.loc[:, 'tail up (nm)']
        t_down_2 = tail_2.loc[:, 'tail down (nm)']
        ax.plot(x_tick, t_up_2, color=(0, 0, 0.5), marker='^', label='tail up', markersize=0.05, linestyle='')
        ax.plot(x_tick, t_down_2, color=(0, 0, 1), marker='v', label='tail down', markersize=0.05, linestyle='')
        label = []


    save_loc = fileio.change_extension(filename, '.png')
    format_plot('iteration', 'distance (nm)', 'title', scale_page=(1.0 / 2.0),
                aspect=0.5, save=save_loc, yrange=[0,24], #legend=label,
                ax=ax)

    return

def plot_tail2(filename, filename_2=None):

    tail = pd.read_excel(filename, header=0, index_col=0)
    x_tick = tail.index[:]
    t_1 = np.mean(tail, axis=1)

    tail_2 = pd.read_excel(filename_2, header=0, index_col=0)
    t_2 = np.mean(tail_2, axis=1)

    fig, ax = plt.subplots()
    ax.plot(x_tick, t_1, color=(0.75, 0, 0.25), marker='o', label='tail up', markersize=0.1, linestyle='')
    ax.plot(x_tick, t_2, color=(0, 0.75, 0.25), marker='o', label='tail up', markersize=0.1, linestyle='')

    label = ['167', '197']

    save_loc = fileio.change_extension(filename, '3.png')
    format_plot('iteration (#)', 'distance (nm)', 'title', scale_page=(1.0 / 2.8),
                aspect=0.8, save=save_loc, yrange=[0,20], legend=label, ax=ax)

    return

