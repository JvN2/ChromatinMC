import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
from helixmc.pose import HelixPose
from helixmc import util
import NucleosomeMC as nMC
import FileIO as fileio
import Tails as tMC

def plotten(x, y, xlabel, ylabel):
    plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()

    ax.plot(x, y, color=(0.75, 0, 0.25), linewidth=5)

    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    ax.set_xlim(left=0)
    # plt.legend()

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)


    plt.show()


def repulsion_exp():
    filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201110\20201112_cms_dists.xlsx"
    zero_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,B,C")
    one_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,D,E")
    two_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,F,G")

    zero_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,B,C")
    one_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,D,E")
    two_start_197 = pd.read_excel(filename, sheet_name='197', header=1, index_col=0, usecols="A,F,G")

    x_1 = range(1, 8)
    x_tick = one_start_167.index

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

    #
    #
    plt.rcParams.update({'font.size': 22})
    #
    fig, ax = plt.subplots()
    #
    ax.errorbar(x_1, y_1, e_1, color=(0.25, 0, 0.25), marker='o', markersize=10, label='NRL 167 1-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y_2, e_2, color=(0.5, 0, 0.25), marker='o', markersize=10, label='NRL 167 2-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y_0, e_0, color=(0.75, 0, 0.25), marker='o', markersize=10, label='NRL 167 0-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)

    ax.errorbar(x_1, y9_1, e9_1, color=(0, 0.25, 0.25), marker='o', markersize=10, label='NRL 197 1-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y9_2, e9_2, color=(0, 0.5, 0.25), marker='o', markersize=10, label='NRL 197 2-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y9_0, e9_0, color=(0, 0.75, 0.25), marker='o', markersize=10, label='NRL 197 0-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)

    #
    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    plt.xticks(x_1, x_tick)
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=0)
    ax.set_ylim(top=35)
    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.legend(loc='upper left')
    #
    # plt.title('Repulsion amplitude 2,5 kT', y=1.08)
    # plt.xlabel('decay length (nm)')
    plt.title('Decay length 5.0 nm', y=1.08)
    plt.xlabel('Amplitude (kT)')
    plt.ylabel('Distance between nucleosomes (nm)')
    plt.show()

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

    plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()

    for d in decay_l:
        ax.plot(x, y[d], color=(d / 100, 0, 1 - (d / 100)), linewidth=5, label='$\lambda$' + ' = ' + str(d/10))
    # for A in Amp:
    #     ax.plot(x, y[A], color=(A / 205.0, 0, 1 - (A / 205.0)), linewidth=5, label='A' + ' = ' + "{:3.1f}".format(A/kT))

    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    ax.set_xlim(left=0)
    plt.legend()
    # plt.title('Decay length 5.0 nm', y=1.08)
    plt.title('Amplitude 2.5 kT', y=1.08)
    plt.xlabel('distance (nm)')
    plt.ylabel('energy (kT)')

    plt.show()

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

    plotten(z_array, g_array, xlabel='distance (nm)', ylabel='energy (kT)')

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


def plot_npz(filename, nrl, nucs):

    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))
    # save parameters of bp of every dna pose in params
    params = []
    for f in npz_f[:]:
        dna = HelixPose.from_file(f)
        # params.append(dna.params)
        # # calculate mean parameters per bp
        # params = np.mean(params, axis=0)
        params = dna.params

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
        coord_w_hist, radius, colors = tMC.get_histones(coords[dyads[3] - 100:dyads[5] + 100], dyads, nucl, tf=tf_d[3:6])

        # transform fiber to origin
        # origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]])
        origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]])
        tf_o = nMC.get_transformation(nuc_cms[3], target=origin_of)
        t_coord = []  # transformed coords
        # Tranform coords where first nucleosome is placed in origin
        for c in coord_w_hist:
            t_coord.append(nMC.apply_transformation(c, tf_o))


        print(fileio.create_pov((fileio.change_extension(f, 'png')), t_coord, radius=radius, colors=colors, range_A=[1000, 1000],
                                offset_A=[0, 0, 300], show=False, width_pix=1500))

    return

def get_stack_params(filename):
    # get list of npz files in filename folder
    xlsx_f = glob.glob(fileio.change_extension(filename, '\*.xlsx'))
    xlsx = []
    nucl_stack_params = []
    for i, file in enumerate(xlsx_f):
        xlsx = pd.read_excel(file, usecols =['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'])
        nucl_stack_params.append(xlsx.iloc[1:-1].mean(axis=0))

    nucl_stack_params = np.mean(nucl_stack_params, axis=0)

    return nucl_stack_params
