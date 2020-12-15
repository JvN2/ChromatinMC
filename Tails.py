import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from helixmc import util
from helixmc.random_step import RandomStepSimple, RandomStepAgg, symmetrize_WC
from helixmc.pose import HelixPose
import glob
import math

# ChromatinMC modules:
import NucleosomeMC as nMC
import FileIO as fileio
import RunMC as rMC
import FiberMC as fMC

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


def get_histones(coord, dyads, nucl, dna=None, tf=None, tail=True):

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
        elif chain == 'H2B*' or chain == 'H2B':
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
        # radius of histone proteins
        radius = np.append(radius, np.ones(8) * 4)

        if tail == True:
            # apply transformation on linker coordinates
            for l in nucl.l_coord:
                coord.append(nMC.apply_transformation(l, tfm))

            # radius of linker-amino-acids
            radius = np.append(radius, np.ones(4) * 13)

            # colors of histone proteins and linker-amino-acids
            # [H3*, H3, H4*, H4, H2A*, H2B*, H2A, H2B]
            colors += 'bbggyryrpmpm'  # 8 histone proteins + H2A, H2A*, H4, H4*

        else:
            # colors of histone proteins and linker-amino-acids
            # [H3*, H3, H4*, H4, H2A*, H2B*, H2A, H2B]
            colors += 'bbggyryr'  # 8 histone proteins


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
    d_up = 0
    d_down = 0

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


    return g_pNnm * 10 # g_pNA


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

    if fiber_start is 0:
        g = 0

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


def origin(dna, dyads, nucl, filename, axis=False):
    """
    first nucleosome is positioned in origin after transformation

    Parameters
    ----------
    filename
    axis
    dna:    Helixpose
    dyads:  indices of all dyads
    nucl:   Nucleosomepose
    coord:  coordinates of basepairs (and histones)

    Returns
    -------
    df_last_c: coordinates of last pose
    df_last_cms: coordinates cms of nucleosome in last pose
    """

    coord, radius, colors = get_histones(dna.coord, dyads, nucl, dna=dna, tail=False)

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
    s_coord = np.array(f_coord[0]) / 10
    df_last_c = pd.DataFrame(s_coord, columns=['x (nm)', 'y (nm)', 'z (nm)'])

    nuc_cms_c = []
    for n, cms in enumerate(nuc_cms):
        nuc_cms_c.append(nMC.apply_transformation(nuc_cms[n][0], tf)[0])

    df_last_cms = pd.DataFrame(np.array(nuc_cms_c) / 10)


    print(fileio.create_pov((fileio.change_extension(filename, '_org.png')), f_coord, radius=radius, colors=colors,
                            range_A=[750, 750], offset_A=[0, 0, 300], show=False, width_pix=1500))

    return df_last_c, df_last_cms


def coord_mean(filename, dyads, nucl, fiber_start, pars, fixed_wrap_params, p0, k):
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
    g_dna = []
    e_wrap_kT = pars['e_wrap_kT'].value
    # params contains all 6 parameters for each basepair out of every npz file
    for f, file in enumerate(npz_f):
        dna = HelixPose.from_file(file)
        params.append(dna.params)
        # frames.append(dna.frames)
        # coords.append(dna.coord)

        w = np.ones(len(dna.params))

        for dyad in dyads:
            fixed_bps = rMC.score_wrapping(dyad + 1, dna.coord, dna.frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                           half_nuc=True)[1]
            if len(fixed_bps) > 0:
                w[fixed_bps[0]:fixed_bps[1]] = 0
            else:
                w = None
        g_dna.append(rMC.score_dna(dna.params, p0, k, w=w))

    g_dna_all = np.mean(g_dna, axis=0) #L_bpx6
    g_dna_m = np.mean(np.sum(g_dna, axis=2), axis=0)


    df_m_g_all = pd.DataFrame(g_dna_all, columns=['g_shift_kT',
               'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT', 'g_twist_kT'], index=range(1,len(g_dna_all) + 1))

    df_m_g = pd.DataFrame(g_dna_m, columns=['g_total (kT)'], index=range(1,len(g_dna_all) + 1))

    # get mean value of every parameter for each basepair
    params_m = np.mean(params, axis=0)
    # save mean parameter values per bp in dataframe
    df_p_m = pd.DataFrame(params_m, columns=['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'],
                          index=range(1,len(params_m) + 1))

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
    coord_w_hist, radius, colors = get_histones(coords, dyads, nucl, tf=tf_d, tail=False)
    coord_3_nuc, radius3, colors3 = get_histones(coords[dyads[3]-100:dyads[5]+100], dyads, nucl, tf=tf_d[3:6], tail=False)

    # transform fiber to origin
    nuc_cms = []
    for d, dyad in enumerate(dyads):
        nuc_cms.append(nMC.apply_transformation(nucl.of, tf_d[d]))


    t_coord = [] # transformed coords
    coord3 = [] # 3 nucleosomes coords
    origin_of = np.asarray([[0, 0, 0], [0.707, 0.707, 0], [0.707, -0.707, 0], [0, 0, -1]])
    # origin_of = np.asarray([[0, 0, 0], [0.866, -0.5, 0], [-0.5, -0.866, 0], [0, 0, -1]])
    tf_o = nMC.get_transformation(nuc_cms[3], target=origin_of)
    # Tranform coords where 4th nucleosome is placed in origin
    for c in coord_w_hist:
        t_coord.append(nMC.apply_transformation(c, tf_o))

    for c in coord_3_nuc:
        coord3.append(nMC.apply_transformation(c, tf_o))


    # create separate list of transformed histone coords
    df_H2A, df_H2B, df_H3, df_H4, df_l_coord = histones_coords(nucl, tf_d, tf_o)

    # # save dna coords in dataframe
    df_m_c = pd.DataFrame(np.array(t_coord[0]) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])

    # transform coord of cms of nucleosomes the same as bp coords and save in dataframe
    nuc_cms_c = []
    nuc_params = []

    for n, cms in enumerate(nuc_cms):
        nuc_cms[n] = nMC.apply_transformation(nuc_cms[n], tf_o)
        nuc_cms_c.append(nuc_cms[n][0])
        if n >= fiber_start:
            if fiber_start != 0:
                nuc_params.append(nMC.ofs2params(nuc_cms[n], nuc_cms[n - fiber_start], _3dna=True, flipx=[0,0]))
            else:
                nuc_params.append(nMC.ofs2params(nuc_cms[n], nuc_cms[n - 1], _3dna=True, flipx=[0, 0]))


    df_cms_c = pd.DataFrame(np.array(nuc_cms_c) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])
    df_nucl_p = pd.DataFrame(nuc_params, columns=['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'],
                             index=range(fiber_start,len(dyads)))


    print(fileio.create_pov((fileio.change_extension(filename, '_m.png')), t_coord, radius=radius, colors=colors, range_A=[1000, 1000],
                            offset_A=[0, 0, 500], show=False, width_pix=1500))
    print(fileio.create_pov((fileio.change_extension(filename, '_3m.png')), coord3, radius=radius3, colors=colors3, range_A=[1000, 1000],
                            offset_A=[0, 0, 200], show=False, width_pix=1500))



    # combine dataframe of bp coords and parameters
    df_m = pd.concat([df_m_c, df_p_m, df_m_g_all, df_m_g], axis=1)
    df_nuc = pd.concat([df_cms_c, df_nucl_p], axis=1)


    return df_m, df_nuc, df_H2A, df_H2B, df_H3, df_H4, df_l_coord


def score_repulsion(moving_bp, fiber_start, dyads, dna, nucl, pars):

    g = 0
    # Repulsion parameters
    Amp = pars['Rep_Amp_pNA'].value  # amplitude (pNA)
    decay_l = pars['Rep_decay_A'].value  # decay length (A)
    left_dyad = np.argmax(dyads > moving_bp) - 1

    if fiber_start != 0:

        right_dyad = left_dyad + fiber_start

        if 0 <= left_dyad < len(dyads) - fiber_start:
            left_nucl_bps = dyads[left_dyad] + nucl.fixed
            right_nucl_bps = dyads[right_dyad] + nucl.fixed

            up_turn = dna.coord[left_nucl_bps]
            down_turn = dna.coord[right_nucl_bps]

            for i, n in enumerate(up_turn):
                for j, m in enumerate(down_turn[i:]):
                    dist = np.sqrt(np.sum((m - n) ** 2))
                    # calculate dist between surface dna
                    g += Amp * np.exp(- (1 / decay_l) * (dist))


        if fiber_start is 2 and left_dyad >= 1:
            left_dyad -= 1
            right_dyad -= 1
            left_nucl_bps = dyads[left_dyad] + nucl.fixed
            right_nucl_bps = dyads[right_dyad] + nucl.fixed

            up_turn = dna.coord[left_nucl_bps]
            down_turn = dna.coord[right_nucl_bps]

            for i, n in enumerate(up_turn):
                for j, m in enumerate(down_turn[i:]):
                    dist = np.sqrt(np.sum((m - n) ** 2))
                    # calculate dist between surface dna
                    g += Amp * np.exp(- (1 / decay_l) * (dist))

    elif fiber_start == 0 and left_dyad >= 0:
        dna_coord = dna.coord
        dna_frames = dna.frames

        nuc_cms = []
        nuc_dist = []
        for dyad in dyads:
            nuc_cms.append(nMC.get_nuc_of(dna_coord, dna_frames, dyad, nucl)[0])
        # check distance of nucleosomes
        for i, cm1 in enumerate(nuc_cms):
            nuc_dist.extend([np.sum((nuc_cms[left_dyad] - cm1) ** 2)])
        # nuc_sort_idx = np.argsort(nuc_dist)

        # determine indices of two dyads closest to left_dyad
        first_dyad = np.argsort(nuc_dist)[1]
        second_dyad = np.argsort(nuc_dist)[2]

        left_nucl_bps = dyads[left_dyad] + nucl.fixed
        first_nucl_bps = dyads[first_dyad] + nucl.fixed
        second_nucl_bps = dyads[second_dyad] + nucl.fixed

        up_turn = dna.coord[left_nucl_bps]
        first_turn = dna.coord[first_nucl_bps]
        second_turn = dna.coord[second_nucl_bps]

        for i, n in enumerate(up_turn):
            for j, m in enumerate(first_turn[i:]):
                dist = np.sqrt(np.sum((m - n) ** 2))
                # calculate dist between surface dna
                g += Amp * np.exp(- (1 / decay_l) * (dist))

        for i, n in enumerate(up_turn):
            for j, m in enumerate(second_turn[i:]):
                dist = np.sqrt(np.sum((m - n) ** 2))
                # calculate dist between surface dna
                g += Amp * np.exp(- (1 / decay_l) * (dist))

    return g


def sequence():
    dna_step_file = util.locate_data_file('DNA_default.npz')

    DNA_default = np.load(dna_step_file)
    DNA_default_p0_k = dict()
    seq_list = [
        'AA', 'AT', 'AG', 'AC', 'GA', 'GC', 'GT', 'GG', 'TT',
        'TA', 'TC', 'TG', 'CC', 'CG', 'CA', 'CT'
    ]

    for seq in seq_list:

        DNA_default_p0_k[seq] = dict()
        DNA_default_p0_k[seq]['params'] = DNA_default[seq]
        DNA_default_p0_k[seq]['avg'] = np.average(DNA_default[seq], axis=0)
        DNA_default_p0_k[seq]['cov'] = np.cov(DNA_default[seq], rowvar=False)



    print('DNA_default k: ', np.linalg.inv(DNA_default_p0_k['GT']['cov']))
    print('DNA_default p0: ', DNA_default_p0_k['TA']['avg'])

    Random_step = RandomStepAgg(data_file=dna_step_file)
    print('random step: ', Random_step.get_rand_step(identifier='TA').params_avg)


    return


def save_values(pars, filename, dyads, nucl, results, results_std, energy_all, fixed_wrap_params, p0, k):

    fiber_start = pars['fiber_start'].value

    mean_sheet, nucl_sheet, df_H2A, df_H2B, df_H3, df_H4, df_l_coord = coord_mean(filename, dyads, nucl, fiber_start, pars, fixed_wrap_params, p0, k)

    for key in energy_all:
        pars[key].value = np.mean(energy_all[key])
        energy_all[key] = np.std(energy_all[key])
    # save mean energy in results df
    results.loc['average'] = pars.valuesdict()
    results_std.loc['average'] = energy_all


    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(filename, engine='xlsxwriter')

    # Write each dataframe to a different worksheet.
    results.to_excel(writer, sheet_name='params')
    results_std.to_excel(writer, sheet_name='std')
    mean_sheet.to_excel(writer, sheet_name='mean')
    nucl_sheet.to_excel(writer, sheet_name='nucl')
    df_H2A.to_excel(writer, sheet_name='H2A')
    df_H2B.to_excel(writer, sheet_name='H2B')
    df_H3.to_excel(writer, sheet_name='H3')
    df_H4.to_excel(writer, sheet_name='H4')
    df_l_coord.to_excel(writer, sheet_name='l_coord')

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

    return


def energy_could_be_our_closest_friend(pars, energy, dyads, dna, nucl, fiber_start, fixed_wrap_params, fixed_stack_params, p0, k, force):

    dna_coord = dna.coord
    dna_params = dna.params
    dna_frames = dna.frames

    e_wrap_kT = pars['e_wrap_kT'].value
    e_stack_kT = pars['e_stack_kT'].value

    g_dna = np.zeros(6)
    g_wrap = 0
    g_stack = 0
    g_tails = 0
    g_rep = 0
    g_work = 0
    nucl_cms = 0
    tail = np.zeros(2)


    n_nucs = pars['n_nuc'].value
    tail_switch = pars['tail_switch'].value

    w = np.ones(len(dna_params))

    for dyad in dyads:
        fixed_bps = rMC.score_wrapping(dyad + 1, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                   half_nuc=True)[1]
        if len(fixed_bps) > 0:
            w[fixed_bps[0]:fixed_bps[1]] = 0
        # else:
        #     w = None
    g_dna_all = rMC.score_dna(dna_params, p0, k, w=w)

    for dyad1, dyad2 in zip(dyads[:-1], dyads[1:]):
        g_dna += np.sum(g_dna_all[dyad1:dyad2], axis=0)
        if tail_switch == True:
            g_tails += score_tails(dyad1 + 1, fiber_start, dyads, dna, nucl)
            g_rep   += score_repulsion(dyad1 + 1, fiber_start, dyads, dna, nucl, pars)
        else:
            g_stack += rMC.score_stacking(dyad1 + 1, dna_coord, dna_frames, dyads, fixed_stack_params, e_stack_kT, nucl,
                                  fiber_start)
        g_work += rMC.score_work(dna_coord, force, start_bp=dyad1, end_bp=dyad2, )
        g_wrap += rMC.score_wrapping(dyad1 + 2, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)[0]
        g_wrap += rMC.score_wrapping(dyad2 - 2, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)[0]

    # calculate distances between nucleosomes and distance of tails
    if fiber_start > 0:
        for n in range(fiber_start, n_nucs):
            nucl_cms += dist_cms(n, n - fiber_start, dna, dyads, nucl)
            tail += tail_dist(n - fiber_start, n, dyads, dna, nucl)
    else:
        for n in range(1, n_nucs):
            nucl_cms += dist_cms(n, n - 1, dna, dyads, nucl)
            tail += tail_dist(n - 1, n, dyads, dna, nucl)


    g_dna /= (n_nucs - 1)
    g_wrap /= (n_nucs - 1)
    g_work /= (n_nucs - 1)

    if fiber_start > 0:
        g_stack /= (n_nucs - fiber_start) * fiber_start
        g_tails /= (n_nucs - fiber_start) * fiber_start
        g_rep /= (n_nucs - fiber_start) * fiber_start
        nucl_cms /= (n_nucs - fiber_start)
        tail /= (n_nucs - fiber_start)
    else:
        g_rep /= (n_nucs - 1)
        nucl_cms /= (n_nucs - 1)
        tail /= (n_nucs - 1)

    energy['g_dna_kT'].extend([np.sum(g_dna)]) # Is already in kT
    energy['g_wrap_kT'].extend([g_wrap  / 41.0]) # kT
    energy['g_stack_kT'].extend([g_stack / 41.0]) # kT
    energy['g_tails_kT'].extend([g_tails / 41.0]) # kT
    energy['g_rep_kT'].extend([g_rep / 41.0])
    energy['g_work_kT'].extend([g_work / 41.0])
    energy['g_total'].extend([np.sum([np.sum(g_dna) * 41.0, g_wrap, g_stack, g_tails, g_rep, g_work]) / 41.0])

    energy['g_shift_kT'].extend([g_dna[0]])
    energy['g_slide_kT'].extend([g_dna[1]])
    energy['g_rise_kT'].extend([g_dna[2]])
    energy['g_tilt_kT'].extend([g_dna[3]])
    energy['g_roll_kT'].extend([g_dna[4]])
    energy['g_twist_kT'].extend([g_dna[5]])

    energy['nucl_cms_nm'].extend([nucl_cms])
    energy['tail_up_nm'].extend([tail[0]])
    energy['tail_down_nm'].extend([tail[1]])

    return


def which_energies(energy):

    g_names = ['g_tails_kT', 'g_rep_kT', 'g_total', 'g_dna_kT', 'g_wrap_kT', 'g_stack_kT', 'g_work_kT', 'g_shift_kT',
               'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT', 'g_twist_kT', 'nucl_cms_nm', 'tail_up_nm',
               'tail_down_nm']
    results_std = pd.DataFrame(columns=g_names)
    for names in g_names:
        energy[names] = []

    return results_std


def histones_coords(nucl, tf_d, tf_o):

    # create separate list of histone coords
    H2A_c = []
    H2B_c = []
    H3_c = []
    H4_c = []

    for chain in nucl.chains:
        if chain == 'DNA':
            pass
        # remove coordinates that are produce sphere in the middle of nowhere
        elif chain == 'H2A*' or chain == 'H2A':
            H2A_c.extend(nucl.chains[chain][2])
        elif chain == 'H2B*' or chain == 'H2B':
            H2B_c.extend(nucl.chains[chain][2][:-3])
        elif chain == 'H3*' or chain == 'H3':
            H3_c.extend(nucl.chains[chain][2])
        elif chain == 'H4*' or chain == 'H4':
            H4_c.extend(nucl.chains[chain][2])

    # transform histone coords to every nucl position
    H2A_ct = []
    H2B_ct = []
    H3_ct = []
    H4_ct = []
    l_coordt = []

    for tf in tf_d:
        H2A_ct.extend(nMC.apply_transformation(nMC.apply_transformation(np.array(H2A_c), tf), tf_o))
        H2B_ct.extend(nMC.apply_transformation(nMC.apply_transformation(np.array(H2B_c), tf), tf_o))
        H3_ct.extend(nMC.apply_transformation(nMC.apply_transformation(np.array(H3_c), tf), tf_o))
        H4_ct.extend(nMC.apply_transformation(nMC.apply_transformation(np.array(H4_c), tf), tf_o))
        l_coordt.extend(nMC.apply_transformation(nMC.apply_transformation(np.array(nucl.l_coord), tf), tf_o))


    df_H2A = pd.DataFrame(np.array(H2A_ct) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])
    df_H2B = pd.DataFrame(np.array(H2B_ct) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])
    df_H3 = pd.DataFrame(np.array(H3_ct) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])
    df_H4 = pd.DataFrame(np.array(H4_ct) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'])


    l_index = ['H2A', 'H2A*', 'H4', 'H4*'] * len(tf_d)

    df_l_coord = pd.DataFrame(np.array(l_coordt) / 10, columns=['x (nm)', 'y (nm)', 'z (nm)'], index=l_index)


    return df_H2A, df_H2B, df_H3, df_H4, df_l_coord


def nuc_pars(dna, dyads, nucl, fiber_start, datafile):

    nuc_pars=[]

    for i, d in enumerate(dyads):
        if fiber_start == 2:
            if i >= fiber_start:
                nuc_pars.append(fMC.get_stack_pars(dna.coord,dna.frames, dyads[i - fiber_start], dyads[i],
                                                   nucl, fiber_start))

        else:
            if i >= 1:
                nuc_pars.append(fMC.get_stack_pars(dna.coord,dna.frames, dyads[i - 1], dyads[i],
                                                       nucl, fiber_start))

    idx = 1
    if fiber_start == 2:
        idx = 2

    df_nuc_pars = pd.DataFrame(nuc_pars, columns=['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'],
                             index=range(idx,len(dyads)))

    df_nuc_pars.to_excel(fileio.change_extension(datafile, 'xlsx'))

    return


def fixed_pars2excel(fixed_stack_pars, datafile):

    df_nuc_pars = pd.DataFrame(np.reshape(fixed_stack_pars, (1,6)), columns=['shift (A)', 'slide (A)', 'rise (A)', 'tilt', 'roll', 'twist'],
                             index=['fixed'])
    df_nuc_pars.to_excel(fileio.change_extension(datafile, '_fixed.xlsx'))
    return