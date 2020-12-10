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
    left, bottom, width, height = [0.6, 0.55, 0.35, 0.35]
    ax2 = fig.add_axes([left, bottom, width, height]) # inset
    #
    ax.errorbar(x_tick, y_1.iloc[:-1], e_1.iloc[:-1], color=(0.75, 0, 0.25), marker='^', markersize=5, label='167 1-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y_1['old'], color=(0.75, 0, 0.25), linestyle='-', lw=3)
    ax.errorbar(x_tick, y_2.iloc[:-1], e_2.iloc[:-1], color=(0.75, 0, 0.25), marker='s', markersize=5, label='167 2-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y_2['old'], color=(0.75, 0, 0.25), linestyle='--', lw=3)

    ax2.errorbar(x_tick, y_0.iloc[:-1], e_0.iloc[:-1], color=(0.75, 0, 0.25), marker='o', markersize=5, label='167 0-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax2.axhline(y=y_0['old'], color=(0.75, 0, 0.25), linestyle='-', lw=2)

    ax.errorbar(x_tick, y9_1.iloc[:-1], e9_1.iloc[:-1], color=(0, 0.75, 0.25), marker='^', markersize=5, label='197 1-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y9_1['old'], color=(0, 0.75, 0.25), linestyle='-', lw=3)
    ax.errorbar(x_tick, y9_2.iloc[:-1], e9_2.iloc[:-1], color=(0, 0.75, 0.25), marker='s', markersize=5, label='197 2-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.axhline(y=y9_2['old'], color=(0, 0.75, 0.25), linestyle='--', lw=3)

    ax2.errorbar(x_tick, y9_0.iloc[:-1], e9_0.iloc[:-1], color=(0, 0.75, 0.25), marker='o', markersize=5, label='197 0-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax2.axhline(y=y9_0['old'], color=(0, 0.75, 0.25), linestyle='-', lw=2)

    labels=['167 1-start fps', '167 2-start fps', '197 1-start fps', '197 2-start fps', '167 1-start', '167 2-start', '197 1-start', '197 2-start',]
    save_loc = fileio.change_extension(filename, '.png')
    format_plot('decay length (nm)', 'nucleosome distance (nm)', 'title', scale_page=1.0,
                aspect=0.5, save=save_loc, yrange=[5.5, 11.1], legend=None, ax=ax)
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
    sv_loc = r"C:\Users\Annelies\Documents\Natuurkunde\Master project 02\figures\repulsion_graphA2_5kTdvaries.png"
    format_plot('distance (nm)', 'energy (kT)', 'title', scale_page=(3.0/5.0),
                aspect=0.7, save=sv_loc, yrange=None, legend=label, ax=ax)

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

def get_nucl_dyads(coords, NRL, n_nucs):
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
    n_bp = len(coords)  # number of bp
    dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
    dyads = (dyads + n_bp // 2).astype(int)

    return nucl, dyads


def plot_npz(filename):

    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))

    pars = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0, usecols=['NRL',
                                                                           'n_nuc'])
    nrl = pars.iloc[0]['NRL']
    nucs = pars.iloc[0]['n_nuc']

    # save parameters of bp of every dna pose in params
    params = []
    for f in npz_f[:]:
        dna = HelixPose.from_file(f)
        params.append(dna.params)
    # # calculate mean parameters per bp
    params = np.mean(params, axis=0)
    # params = dna.params

    # use 6 parameters to get coordinates of every basepair
    dr, frames = util.params2data(params)
    coords = util.dr2coords(dr)

    # get nucleosome pose en list of dyads
    nucl, dyads = get_nucl_dyads(coords, nrl, nucs)

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
    # coord_w_hist, radius, colors = tMC.get_histones(coords, dyads, nucl, tf=tf_d)
    coord_w_hist, radius, colors = tMC.get_histones(coords[dyads[3] - 75:dyads[5] + 75], dyads, nucl, tf=tf_d[3:6])

    # transform fiber to origin
    origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]])
    # origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]])
    tf_o = nMC.get_transformation(nuc_cms[3], target=origin_of)
    t_coord = []  # transformed coords
    # Tranform coords where first nucleosome is placed in origin
    for c in coord_w_hist:
        t_coord.append(nMC.apply_transformation(c, tf_o))


        # print(fileio.create_pov((fileio.change_extension(f, 'png')), t_coord, radius=radius, colors=colors, range_A=[1000, 1000],
        #                         offset_A=[0, 0, 300], show=False, width_pix=1500))

    return t_coord[0]

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

def dna_energy_display(filename, energy_kT='g_total (kT)'):
    pars = pd.read_excel(fileio.change_extension(filename, '.xlsx'), sheet_name='params', header=0, usecols=['NRL',
                                                                                                 'L_bp', 'n_nuc'])
    NRL = pars.iloc[0]['NRL']
    n_nucs = pars.iloc[0]['n_nuc']
    nucl = nMC.NucPose()
    nucl.from_file('1KX5')
    # get list of dyad indices
    n_bp = pars.iloc[0]['L_bp']  # number of bp
    dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
    dyads = (dyads + n_bp // 2).astype(int)
    begin = dyads[3] - 75
    end = dyads[5] + 75

    file = pd.read_excel(filename, sheet_name='mean', header=0, usecols =['g_shift_kT',
               'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT', 'g_twist_kT', 'g_total (kT)']).iloc[begin:end]
    # coords = pd.read_excel(filename, sheet_name='mean', header=0, usecols =['x (nm)', 'y (nm)', 'z (nm)']).iloc[begin:end]
    coords = plot_npz(filename)
    coord_n =[]
    for c in coords:
        coord_n.append(c/10.0)

    # energy = np.sqrt(file[energy_kT]**2)/np.amax(file[energy_kT])
    # file[energy_kT] = (np.sqrt(file[energy_kT] ** 2))
    n_bins = 50
    bin_max = 1.5 #kT
    bin = np.linspace(0, bin_max, n_bins, dtype=float)
    bin_labels = np.arange(n_bins-1)
    file[energy_kT] = np.clip(file[energy_kT], 0 , bin_max)

    colorwaaier = []
    rod = []
    for i, b in enumerate(bin):
        p = b/bin_max
        colorwaaier.append([0.6 + (0.4 * p) , 0.6 * (1 - p)**2,  0.6 * (1 - p)**2])
        rod.append([0,0,i/10.0])

    POVe.main(fileio.change_extension(filename, 'label.png'), rod, colorwaaier, radius=1, range_A=[25, 25], offset_A=[0, 0, 10], width_pix=500, showt=True)

    file['bin'] = pd.cut(file[energy_kT], bins=bin, labels=bin_labels, include_lowest=True)
    colors = []

    for b in file['bin']:
        colors.append(colorwaaier[b])

    # POVe.main(fileio.change_extension(filename, 'Etot.png'), coords.to_numpy(), colors, radius=1, range_A=[50, 50], offset_A=[0, 0, 25], width_pix=500, showt=True)
    POVe.main(fileio.change_extension(filename, 'Etot.png'), coord_n, colors, radius=1, range_A=[50, 50], offset_A=[0, 0, 25], width_pix=500, showt=True)

    return

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
        plt.legend(legend, prop={'size': fontsize * 0.8}, )

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

    filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201210 Nucleosome stacking.xlsx"

    fixed = pd.read_excel(filename, sheet_name='fixed', header=0, index_col=0, usecols="A,B,C,D,E,F,G")
    NRL167 = pd.read_excel(filename, sheet_name='167', header=0, index_col=0)
    NRL197 = pd.read_excel(filename, sheet_name='197', header=0, index_col=0)
    # print(zero_start_167.loc['fps 0', param_name])
    # one_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,D,E")
    # two_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,F,G")
    #
    # zero_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,B,C")
    # one_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,D,E")
    # two_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,F,G")
    #
    # x_1 = range(1, 8)
    # x_tick = one_start_167.index[:-1]
    #
    # y_0 = zero_start_167.iloc[:, 0]
    # e_0 = zero_start_167.iloc[:, 1]
    #
    # y_1 = one_start_167.iloc[:, 0]
    # e_1 = one_start_167.iloc[:, 1]
    #
    # y_2 = two_start_167.iloc[:, 0]
    # e_2 = two_start_167.iloc[:, 1]
    #
    # y9_0 = zero_start_197.iloc[:, 0]
    # e9_0 = zero_start_197.iloc[:, 1]
    #
    # y9_1 = one_start_197.iloc[:, 0]
    # e9_1 = one_start_197.iloc[:, 1]
    #
    # y9_2 = two_start_197.iloc[:, 0]
    # e9_2 = two_start_197.iloc[:, 1]
    #
    fig, ax = plt.subplots()
    ax.get_xaxis().set_visible(False)
    # left, bottom, width, height = [0.6, 0.55, 0.35, 0.35]
    # ax2 = fig.add_axes([left, bottom, width, height]) # inset
    ax.axhline(y=fixed.loc['1-start', param_name], color=(0.6, 0.6, 0.6), linestyle='-', lw=2)
    ax.axhline(y=fixed.loc['2-start', param_name], color=(0.6, 0.6, 0.6), linestyle='--', lw=2)

    # FPS 167
    # ax.errorbar(1, NRL167.loc['fps 0', param_name], NRL167.loc['fps 0', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='o', markersize=5, label='fps 0-start', linewidth=0,
    #             ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(2, NRL167.loc['fps 1', param_name], NRL167.loc['fps 1', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='^', markersize=5, label='fps 1-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(3, NRL167.loc['fps 2', param_name], NRL167.loc['fps 2', param_name + ' std'], color=(0.2, 0.2, 0.2), marker='s', markersize=5, label='fps 2-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)

    # 1.31 nm 167
    # ax.errorbar(5, NRL167.loc['1.31 0', param_name], NRL167.loc['1.31 0', param_name + ' std'], color=(0.75, 0, 0.25), marker='o', markersize=5, label='1.31 0-start', linewidth=0,
    #             ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(6, NRL167.loc['1.31 1', param_name], NRL167.loc['1.31 1', param_name + ' std'], color=(0.75, 0, 0.25), marker='^', markersize=5, label='1.31 1-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(7, NRL167.loc['1.31 2', param_name], NRL167.loc['1.31 2', param_name + ' std'], color=(0.75, 0, 0.25), marker='s', markersize=5, label='1.31 2-start', linewidth=0,
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)

    # 7.91 nm 167
    # ax.errorbar(5, NRL167.loc['7.91 0', param_name], NRL167.loc['7.91 0', param_name + ' std'], mec=(0.75, 0, 0.25), marker='o', markersize=5, mfc=(1, 1, 1), label='7.91 0-start',
    #             ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(6, NRL167.loc['7.91 1', param_name], NRL167.loc['7.91 1', param_name + ' std'], mec=(0.75, 0, 0.25), marker='^', markersize=5, mfc=(1, 1, 1), label='7.91 1-start',
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(7, NRL167.loc['7.91 2', param_name], NRL167.loc['7.91 2', param_name + ' std'], mec=(0.75, 0, 0.25), marker='s', markersize=5, mfc=(1, 1, 1), label='7.91 2-start',
                ecolor=(0.75, 0, 0.25), elinewidth=2, capsize=3)

    # # FPS 197
    # ax.errorbar(9, NRL197.loc['fps 0', param_name], NRL197.loc['fps 0', param_name + ' std'], color=(0.2, 0.2, 0.2),
    #             marker='o', markersize=5, label='fps 0-start', linewidth=0,
    #             ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(10, NRL197.loc['fps 1', param_name], NRL197.loc['fps 1', param_name + ' std'], color=(0.2, 0.2, 0.2),
                marker='^', markersize=5, label='fps 1-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    ax.errorbar(11, NRL197.loc['fps 2', param_name], NRL197.loc['fps 2', param_name + ' std'], color=(0.2, 0.2, 0.2),
                marker='s', markersize=5, label='fps 2-start', linewidth=0,
                ecolor=(0.2, 0.2, 0.2), elinewidth=2, capsize=3)
    #
    # # 1.31 nm 197
    # ax.errorbar(13, NRL197.loc['1.31 0', param_name], NRL197.loc['1.31 0', param_name + ' std'], color=(0, 0.75, 0.25),
    #             marker='o', markersize=5, label='1.31 0-start', linewidth=0,
    #             ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(14, NRL197.loc['1.31 1', param_name], NRL197.loc['1.31 1', param_name + ' std'], color=(0, 0.75, 0.25),
                marker='^', markersize=5, label='1.31 1-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(15, NRL197.loc['1.31 2', param_name], NRL197.loc['1.31 2', param_name + ' std'], color=(0, 0.75, 0.25),
                marker='s', markersize=5, label='1.31 2-start', linewidth=0,
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    #
    # # 7.91 nm 197
    # ax.errorbar(13, NRL197.loc['7.91 0', param_name], NRL197.loc['7.91 0', param_name + ' std'], mec=(0, 0.75, 0.25),
    #             marker='o', markersize=5, mfc=(1, 1, 1), label='7.91 0-start',
    #             ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(14, NRL197.loc['7.91 1', param_name], NRL197.loc['7.91 1', param_name + ' std'], mec=(0, 0.75, 0.25),
                marker='^', markersize=5, mfc=(1, 1, 1), label='7.91 1-start',
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    ax.errorbar(15, NRL197.loc['7.91 2', param_name], NRL197.loc['7.91 2', param_name + ' std'], mec=(0, 0.75, 0.25),
                marker='s', markersize=5, mfc=(1, 1, 1), label='7.91 2-start',
                ecolor=(0, 0.75, 0.25), elinewidth=2, capsize=3)
    #
    # labels=['167 1-start fps', '167 2-start fps', '197 1-start fps', '197 2-start fps', '167 1-start', '167 2-start', '197 1-start', '197 2-start',]
    save_loc = fileio.change_extension(filename, '.png')
    format_plot(' ', param_name, 'title', scale_page=(1.0/3),
                aspect=1, save=save_loc, yrange=None, legend=None, ax=ax)
    return