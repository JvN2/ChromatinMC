import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from pynverse import inversefunc
from helixmc import util
import glob
import math

# ChromatinMC modules:
import NucleosomeMC as nMC
import FileIO as fileio

def coth(x):
    return np.cosh(x) / np.sinh(x)


def Langevin(x):
    return (coth(x) - 1.0 / x)


# initial parameters
L_nm = 6.8      # contour length H4 tial
b_nm = 0.6      # kuhn length H4 tail
S_pN = 630      # stiffness H4 tail
kT = 4.10
# possible force on H4 tials
f_array = np.linspace(1e-4, 4800, 1e6)
# corresponding extension of H4 tials
z_array = L_nm * (Langevin(b_nm * f_array / kT) + f_array / S_pN)
print('z max: ', np.max(z_array))
# corresponding energy of H4 tials
g_array = -(kT * L_nm / b_nm) * (np.log((b_nm * f_array) / (kT)) - np.log(
        np.sinh((b_nm * f_array) / (kT)))) + L_nm * f_array ** 2 / (2 * S_pN)

def tf_dyad(dyads, dna, nucl):
    '''
    get transformation matrix for every dyad

    Parameters
    ----------
    dyads:  indices of dyads in fiber
    dna:    HelixPose
    nucl:   nucleosome pose

    Returns
    -------
    tf: transformation matrix, for each dyad

    '''
    of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)  # origin frame dyad in nucl pose
    of_d_fiber = []  # origin frame of dyad in fiber
    tf = []  # transformation matrix
    for i, d in enumerate(dyads):
        # define origin frame of dyad in fiber
        of_d_fiber.append((nMC.get_of(dna, d)))
        # get transformation matrix of nucleosome dyad onto fiber dyad
        tf.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

    return tf


def get_histones(coord, dyads, nucl, dna=None, tf=None):

    '''
    projects histones coordinates onto every
    nucleosome in fiber pose, appends these
    coordinates to coord

    Parameters
    ----------
    coord:  coordinates of helixPose
    dyads:  indices of dyads in fiber
    dna:    HelixPose
    nucl:   nucleosome pose

    Returns
    -------
    coord:  [[DNA], n*[[histones],[linker-coord]]] n = number of nucleosomes
    radius: list
    color:  string

    '''
    # coordinates of histone proteins
    p_coord = []
    for chain in nucl.chains:
        if chain == 'DNA':
            pass
        # remove coordinates that are produce sphere in the middle of nowhere
        elif chain == 'H2B*':
            p_coord.append(nucl.chains[chain][2][:-3])
        elif chain == 'H2B':
            p_coord.append(nucl.chains[chain][2][:-3])
        else:
            p_coord.append(nucl.chains[chain][2])


    radius = [10]  # radius of DNA in POVray
    colors = 'v'  # color of DNA
    coord = [coord]

    if tf == None:
        tf = tf_dyad(dyads, dna, nucl)  # transformation matrix for every dyad

    for tfm in tf:
        # apply transformation on coordinates of histone proteins
        for c in p_coord:
            coord.append(nMC.apply_transformation(c, tfm))
        # apply transformation on linker coordinates
        for l in nucl.l_coord:
            coord.append(nMC.apply_transformation(l, tfm))
        # radius of histone proteins
        radius = np.append(radius, np.ones(8) * 4)
        # radius of linker-amino-acids
        radius = np.append(radius, np.ones(4) * 13)
        # colors of histone proteins and linker-amino-acids
        # [H3*, H3, H4*, H4, H2A*, H2B*, H2A, H2B]
        colors += 'bbggyryrpmpm'  # 8 histone proteins + H2A, H2A*, H4, H4*

    return coord, radius, colors


def tail_dist(dyad_1, dyad_2, dyads, dna, nucl, orientation=None):
    """

    Parameters
    ----------
    dyad_1:         index of first dyad
    dyad_2:         index of second dyad
    orientation:    orientation of two nucleosomes
    dyads:          indices of all dyads
    dna:            helixpose
    nucl:           nucleosomepose

    Returns
    -------
    d_up:           distance from dyad 1 to dyad 2 (nm)
    d_down:         distance from dyad 2 to dyad 1 (nm)

    """
    if orientation == None:
        orientation = '-*'
    l_coord = np.asarray(nucl.l_coord)
    n_l_coord = []  # new coordinates of linker coordinates after transformation
    tf_matrix = tf_dyad(dyads, dna, nucl)

    for tf in tf_matrix:
        # apply transformation on coordinates of linker coordinates
        n_l_coord.append(nMC.apply_transformation(l_coord, tf))

    # orientation of tails and patches
    # * (star) and - (no star) refers to hist_int
    # '*-': top dyad 1 connects to bottom dyad 2
    # '-*': bottom dyad 1 connects to top dyad 2
    # '**': top dyad 1 connects to top dyad 2
    # '--': bottom dyad 1 connects to bottom dyad 2

    hist_int = {'H2A': 0, 'H2A*': 1, 'H4': 2, 'H4*': 3}

    tail_star_1 = n_l_coord[dyad_1][hist_int['H4*']]
    tail_stripe_1 = n_l_coord[dyad_1][hist_int['H4']]
    tail_star_2 = n_l_coord[dyad_2][hist_int['H4*']]
    tail_stripe_2 = n_l_coord[dyad_2][hist_int['H4']]

    patch_star_1 = n_l_coord[dyad_1][hist_int['H2A*']]
    patch_stripe_1 = n_l_coord[dyad_1][hist_int['H2A']]
    patch_star_2 = n_l_coord[dyad_2][hist_int['H2A*']]
    patch_stripe_2 = n_l_coord[dyad_2][hist_int['H2A']]

    if orientation == '*-':
        d_up = np.sqrt(np.sum((tail_star_1 - patch_stripe_2) ** 2))
        d_down = np.sqrt(np.sum((tail_stripe_2 - patch_star_1) ** 2))
    elif orientation == '-*':
        d_up = np.sqrt(np.sum((tail_stripe_1 - patch_star_2) ** 2))
        d_down = np.sqrt(np.sum((tail_star_2 - patch_stripe_1) ** 2))
    elif orientation == '**':
        d_up = np.sqrt(np.sum((tail_star_1 - patch_star_2) ** 2))
        d_down = np.sqrt(np.sum((tail_star_2 - patch_star_1) ** 2))
    elif orientation == '--':
        d_up = np.sqrt(np.sum((tail_stripe_1 - patch_stripe_2) ** 2))
        d_down = np.sqrt(np.sum((tail_stripe_2 - patch_stripe_1) ** 2))

    # print('d up: ', d_up)
    # print('d down: ', d_down)
    d_up /= 10
    d_down /= 10

    return d_up, d_down


def tail_plot(filename, tails, save=False):
    """

    Parameters
    ----------
    tails: distance between tails of two nucleosomes

    Returns
    -------

    """
    dist_up_nm = [d[0] for d in tails]
    dist_down_nm = [d[1] for d in tails]

    plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()
    # # if orientation is '*-':
    # ax.plot(dist_up_nm, color=(1,0,1), marker='o', label='tail up', markersize=12, linestyle='')
    # ax.plot(dist_down_nm, color=(0.75,0,0.25), marker='o', label='tail down', markersize=12, linestyle='')
    # if orientation is '-*':
    ax.plot(dist_up_nm, color=(0.75, 0, 0.25), marker='o', label='tail up', markersize=2, linestyle='')
    ax.plot(dist_down_nm, color=(1, 0, 1), marker='o', label='tail down', markersize=2, linestyle='')
    # default plot parameters

    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # spines = ax.spines
    # [i.set_linewidth(2) for i in spines.values()]
    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    # ax.xaxis.set_tick_params(width=5, size=5)
    # ax.yaxis.set_tick_params(width=5, size=10)
    # ax.set_ylim(bottom=0, top=(max(dist_down_nm) + 5))

    plt.legend(frameon=False, loc=1, markerscale=6)
    plt.ylabel('Distance (nm)')
    plt.xlabel('Iteration (#)')

    if save:
        # save plot
        fig.set_size_inches(16, 9)
        fig.savefig(fileio.change_extension(filename, 'tails.png'), dpi=300)
        # # save tails in xlsx
        df = pd.DataFrame(np.array(tails), columns=['Tail up (nm)', 'Tail down (nm)'])
        df.to_excel(fileio.change_extension(filename, 'tails.xlsx'), index=False, header=True)

    plt.show()

    return


def dist_plot(filename, dist, save=False):
    """

    Parameters
    ----------
    dist: distance between center of mass of two nucleosomes

    Returns
    -------

    """

    plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()
    # # if orientation is '*-':
    # ax.plot(dist_up_nm, color=(1,0,1), marker='o', label='tail up', markersize=12, linestyle='')
    # ax.plot(dist_down_nm, color=(0.75,0,0.25), marker='o', label='tail down', markersize=12, linestyle='')
    # if orientation is '-*':
    ax.plot(dist, color=(1, 0, 1), marker='o', label='distance', markersize=5, linestyle='')
    # default plot parameters

    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # spines = ax.spines
    # [i.set_linewidth(2) for i in spines.values()]
    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    # ax.xaxis.set_tick_params(width=5, size=5)
    # ax.yaxis.set_tick_params(width=5, size=10)
    # ax.set_ylim(bottom=0, top=(max(dist) + 5))

    # plt.legend(frameon=False, loc=0, markerscale=6)
    plt.ylabel('Distance (nm)')
    plt.xlabel('Iteration (#)')

    if save:
        # save plot
        fig.set_size_inches(16, 9)
        fig.savefig(fileio.change_extension(filename, 'cms_dist.png'), dpi=300)
        # # save tails in xlsx
        df = pd.DataFrame(np.array(dist), columns=['cms distance (nm)'])
        df.to_excel(fileio.change_extension(filename, 'cms_dist.xlsx'), index=False, header=True)

    plt.show()

    return


def expected_value():

    f_pN = np.linspace(1e-9, 4000, 1000)
    z_m = []

    g_f = -(kT * L_nm / b_nm) * (np.log((b_nm * f_pN) / (kT)) - np.log(
        np.sinh((b_nm * f_pN) / (kT)))) + L_nm * f_pN ** 2 / (2 * S_pN)
    print('max g_f: ', np.max(g_f))

    z = (L_nm * (Langevin(b_nm * f_pN / kT) + f_pN / S_pN))
    print('expected function max z: ', np.max(z))

    for i, f in enumerate(f_pN):

        expo = (g_f - f * z)
        expo -= np.min(expo)
        # print("expo: ", expo)

        p_f_z = np.exp((-(expo)) / kT)
        p_f_z /= sum(p_f_z)
        # print('p_f_z: ', p_f_z)

        z_m.append(sum(z * p_f_z))


    plt.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots()

    ax.plot(z_m, f_pN, color=(0.75, 0, 0.75), linewidth=5)

    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    ax.set_xlim(left=0)
    # plt.legend()

    plt.ylabel('Force (pN)')
    plt.xlabel('<z>')
    # plt.xlabel('z (nm)')
    # plt.ylabel('z (nm)')
    # plt.ylabel('G(kt)')
    # plt.xlabel('Force (pN)')
    plt.show()


def fFJC(z_nm):

    z = lambda f: L_nm * (Langevin(b_nm * f / kT) + f / S_pN)
    f_pN = inversefunc(z, y_values=z_nm, domain=[(1e-7), 15000], image=[0, 15000])

    return f_pN


def find_nearest(array, value):
    """

    Parameters
    ----------
    array
    value

    Returns
    -------
    idx:            index in array of nearest value
    array[idx]:     nearest value in array

    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
        return idx - 1, array[idx - 1]
    else:
        return idx, array[idx]


def gFJC(z_nm):

    """

    Parameters
    ----------
    z_nm:       distance
    L_nm:       contour length
    b_nm:       Kuhnlength
    S_pN:       Stretch modulus (pN)

    Returns
    -------

    """
    indx, z = find_nearest(z_array, z_nm)
    g_pNnm = g_array[indx]

    #
    # g_pNnm = -(kT * L_nm / b_nm) * (np.log((b_nm * f_pN) / (kT)) - np.log(
    #     np.sinh((b_nm * f_pN) / (kT)))) + L_nm * f_pN ** 2 / (2 * S_pN)
    #
    # P_z = np.exp((z_array - z_nm)**4)
    # P_z /= np.sum(P_z)
    # g_pNnm_e = np.sum(P_z * g_array)


    return g_pNnm / kT


def score_tails(moving_bp, fiber_start, dyads, dna, nucl):
    left_dyad = np.argmax(dyads > moving_bp) - 1
    right_dyad = left_dyad + fiber_start
    g = 0

    if 0 <= left_dyad < len(dyads) - fiber_start:
        t_up, t_down = tail_dist(left_dyad, right_dyad, dyads, dna, nucl)
        g += gFJC(t_up)
        g += gFJC(t_down)

    if fiber_start is 2 and left_dyad >= 1:
        left_dyad -= 1
        right_dyad -= 1
        t_up, t_down = tail_dist(left_dyad, right_dyad, dyads, dna, nucl)
        g += gFJC(t_up)
        g += gFJC(t_down)
    # print('g: ', g)
    return g


def dist_cms(nuc_1, nuc_2, dna, dyads, nucl):
    """

    Parameters
    ----------
    dna:    helixpose
    dyads:  indices of all dyads
    nucl:   nucleosomepose

    Returns
    -------
    distance between center of masses of nucleosomes (nm)
    """
    nuc_cms = []
    for dyad in dyads:
        nuc_cms.append(nMC.get_nuc_of(dna.coord, dna.frames, dyad, nucl)[0])

    dist = np.sqrt(np.sum((nuc_cms[nuc_1] - nuc_cms[nuc_2]) ** 2))

    return dist/10


def origin(dna, dyads, nucl, coord, filename, axis=False):
    """
    first nucleosome is positioned in origin after transformation

    Parameters
    ----------
    dna:    Helixpose
    dyads:  indices of all dyads
    nucl:   Nucleosomepose
    coord:  coordinates of basepairs (and histones)

    Returns
    -------
    f_coord: transformed coordinates
    """

    nuc_cms = []
    for d, dyad in enumerate(dyads):
        nuc_cms.append(nMC.get_nuc_of(dna.coord, dna.frames, dyads[d], nucl))

    if axis == True:
        coord.append(nMC.of2axis(nuc_cms[0]))

    f_coord = []
    tf = nMC.get_transformation(nuc_cms[0], target=np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]]))
    for c in coord:
        f_coord.append(nMC.apply_transformation(c, tf))

    # save dna coords to s_coord
    s_coord = f_coord[0]
    df = pd.DataFrame(np.array(s_coord) / 10)
    df.to_excel(fileio.change_extension(filename, 'coord.xlsx'), index=False, header=True)

    nuc_cms_c = []
    for n, cms in enumerate(nuc_cms):
        nuc_cms_c.append(nMC.apply_transformation(nuc_cms[n][0], tf)[0])

    dfc = pd.DataFrame(np.array(nuc_cms_c) / 10)
    dfc.to_excel(fileio.change_extension(filename, 'coord_cms.xlsx'), index=False, header=True)

    return f_coord


def coord_mean(filename, dyads, nucl):
    """

    Parameters
    ----------
    filename
    dyads
    nucl

    Returns
    -------

    """


    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))

    params = []
    # params contains all 6 parameters for each basepair out of every npz file
    for f, file in enumerate(npz_f):
        data = np.load(file)
        params.append(data['params'])

    # get mean value of every parameter for each basepair
    params_m = np.mean(params,axis=0)

    # use 6 parameters to get coordinates of every basepair
    dr, frames = util.params2data(params_m)
    coords = util.dr2coords(dr)

    # get dyad_ofs to project histones in mean fiber pose
    of_d_fiber = []  # origin frame of dyad in fiber
    of_d_nucl = nMC.get_of(nucl.dna, nucl.dyad)  # origin frame dyad in nucl pose
    tf_d = []  # transformation matrix

    for i, d in enumerate(dyads):

        # get coord en frame of dyad
        dyad_coord = coords[dyads[i]]
        dyad_frame = frames[dyads[i]]

        # combine coords and frames to origin frame of dyad
        of_d_fiber.append(np.concatenate(([dyad_coord], dyad_frame.T+dyad_coord), axis=0))

        # get transformation matrix of nucleosome dyad onto fiber dyad
        tf_d.append(nMC.get_transformation(of_d_nucl, of_d_fiber[i]))

    # append histone positions to coordinates
    coord_w_hist, radius, colors = get_histones(coords, dyads, nucl, tf=tf_d)


    # transform fiber to origin
    nuc_cms = []
    for d, dyad in enumerate(dyads):
        nuc_cms.append(nMC.get_nuc_of(coords, frames, dyads[d], nucl))

    t_coord = [] # transformed coords
    # tf_o = nMC.get_transformation(nuc_cms[0], target=np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]]))
    tf_o = nMC.get_transformation(nuc_cms[0],
                                  target=np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]]))
    for c in coord_w_hist:
        t_coord.append(nMC.apply_transformation(c, tf_o))

    # save dna coords to s_coord
    s_coord = t_coord[0]
    df = pd.DataFrame(np.array(s_coord) / 10)
    df.to_excel(fileio.change_extension(filename, 'coord_m.xlsx'), index=False, header=True)

    nuc_cms_c = []
    for n, cms in enumerate(nuc_cms):
        nuc_cms_c.append(nMC.apply_transformation(nuc_cms[n][0], tf_o)[0])

    dfc = pd.DataFrame(np.array(nuc_cms_c) / 10)
    dfc.to_excel(fileio.change_extension(filename, 'coord_m_cms.xlsx'), index=False, header=True)
    print(fileio.create_pov((fileio.change_extension(filename, '_m.png')), t_coord, radius=radius, colors=colors, range_A=[750, 750],
                            offset_A=[0, 0, 150], show=True, width_pix=1500))

    return


def score_repulsion(moving_bp, fiber_start, dyads, dna, nucl):
    left_dyad = np.argmax(dyads > moving_bp) - 1
    right_dyad = left_dyad + fiber_start

    g = 0
    rep_dist = []

    if 0 <= left_dyad < len(dyads) - fiber_start:
        up_turn = dna.coord[left_dyad: left_dyad + 74]
        down_turn = dna.coord[right_dyad - 74: right_dyad]

        for n in up_turn:
            for m in down_turn:
                rep_dist.append(np.sqrt(np.sum((m - n) ** 2)))

    if fiber_start is 2 and left_dyad >= 1:
        left_dyad -= 1
        right_dyad -= 1
        up_turn = dna.coord[left_dyad: left_dyad + 74]
        down_turn = dna.coord[right_dyad - 74: right_dyad]

        for n in up_turn:
            for m in down_turn:
                rep_dist.append(np.sqrt(np.sum((m - n) ** 2)))