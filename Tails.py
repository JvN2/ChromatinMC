import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from pynverse import inversefunc
from helixmc import util
import glob

# ChromatinMC modules:
import NucleosomeMC as nMC
import FileIO as fileio

kT = 4.10


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


def get_histones(coord, dyads, dna, nucl):

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
        colors += 'bbggryrypmpm'  # 8 histone proteins + H2A, H2A*, H4, H4*

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

    return d_up/10, d_down/10


def tail_plot(filename, tails, save=False):
    """

    Parameters
    ----------
    tails: distance between tails of two nucleosomes

    Returns
    -------

    """
    dist_up_nm = tails[0]
    dist_down_nm = tails[1]

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


def coth(x):
    return np.cosh(x) / np.sinh(x)


def Langevin(x):
    return (coth(x) - 1.0 / x)

def expected_value():

    # z = np.linspace(0, 20, 100)
    f_pN = np.linspace(1e-9, 2000, 100)
    # f_pN = 50
    z_m = []
    L_nm = 6.8
    b_nm = np.linspace(1e-7, 0.6, 6)
    S_pN = 6300
    g_old = []
    # g = (L_nm * kT / b_nm) * (-np.log(f_pN) + np.log(np.tanh(b_nm * f_pN / kT))
    #                           + np.log(np.cosh(b_nm * f_pN / kT)) + (b_nm * f_pN**2) / (2 * kT * S_pN))
    # #
    for b in b_nm:
        g_old.append(-(kT * L_nm / b) * (np.log((b * f_pN) / (kT)) - np.log(
        np.sinh((b * f_pN) / (kT)))) + L_nm * f_pN ** 2 / (2 * S_pN))

    #
    # g_wiki = (kT * L_nm / b_nm)*(np.log(4*np.pi*np.sinh(f_pN * b_nm / kT))
    #                                - np.log(f_pN * b_nm / kT)) + L_nm * f_pN ** 2 / (2 * S_pN)
    #
    # g_meng = L_nm * (f_pN - np.sqrt(f_pN * kT/(2 * b_nm)) + (f_pN**2)/(2 * S_pN))

    z = []
    for b in b_nm:
        z.append(L_nm * (Langevin(b * f_pN / kT) + f_pN / S_pN))
    print(z)


    # coef = np.polyfit(z, f_pN, 6)
    # p = np.poly1d(coef)
    # print(coef)

    #
    # g = np.exp(f_pN/1000)
    # g -= np.min(g)
    # g /= kT
    g_old -= np.min(g_old)
    g_old /= kT
    # g_wiki -= np.min(g_wiki)
    # g_wiki /= kT
    # g_meng -= np.min(g_meng)
    # g_meng /= kT


    # for i, f in enumerate(f_pN):
    #
    #     # g_f = -(kT * L_nm / b_nm) * (np.log((b_nm * f) / (kT)) - np.log(
    #     #     np.sinh((b_nm * f) / (kT)))) + L_nm * f ** 2 / (2 * S_pN)
    #
    #
    #     #
    #     # print(np.log(4*np.pi*np.sinh(f * b_nm / kT))
    #     #                            - np.log(f * b_nm / kT))
    #     # g.append(g_f)
    #
    #     expo = (-(g - f * z))
    #
    #
    #     # fig, ax = plt.subplots()
    #     # ax.plot(expo, color=(0.75, 0, 0.25), linewidth=5)
    #     # plt.show()
    #     # expo = (-(g[i] - f * z))
    #
    #     # print(expo)
    #     # expo /= 10
    #     p_f_z = np.exp(expo / kT)
    #     p_f_z /= sum(p_f_z)
    #     # # p_f_z *= 1000
    #     #
    #
    #     #
    #     #
    #     z_m.append(sum(z * p_f_z))

    # return
    plt.rcParams.update({'font.size': 22})

    # f =0.1
    # g_f = [-(kT * L_nm / b_nm) * (np.log((b_nm * f) / (kT)) - np.log(
    #     np.sinh((b_nm * f) / (kT)))) + L_nm * f ** 2 / (2 * S_pN) for f in force]
    # g_f -= np.min(g_f)
    # g -= np.min(g)
    # g /= kT

    # p_f_z = np.exp(-(g_f - f * z) / kT)
    # p_f_z /= sum(p_f_z)
    # plt.plot(force, g_f)
    # plt.show()
    # return

    fig, ax = plt.subplots()
    # ax.plot(z, f_pN, color=(0.75, 0, 0.25), linewidth=5)
    # ax.plot(p(np.linspace(0,20)), color=(0.28, 0, 0.75), linewidth=5, label='poly')
    # ax.plot(z, g, color=(0.75, 0, 0.25), linewidth=5, label='intergral-z')
    for i, b in enumerate(b_nm):
        ax.plot(z[i], g_old [i], linewidth=5, label=str(b))
    # ax.plot(g_old[1], color=(0.28, 0, 0.10), linewidth=5, label='b = 0.1')
    # ax.plot(g_old[2], color=(0.28, 0, 0.20), linewidth=5, label='b = 0.2')
    # ax.plot(g_old[3], color=(0.28, 0, 0.30), linewidth=5, label='b = 0.3')
    # ax.plot(z, g_wiki, color=(0.28, 0, 0.28), linewidth=5, label='wiki-g')
    # ax.plot(f_pN, g_meng, color=(0.75, 0, 0.75), linewidth=5, label='meng-g')

    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    plt.legend()

    # plt.ylabel('Force (pN)')
    # # plt.xlabel('<z>')
    plt.xlabel('z (nm)')
    plt.ylabel('G(kt)')
    # plt.xlabel('Force (pN)')
    plt.show()

def fFJC(z_nm, L_nm=6.8, b_nm=0.22, S_pN=6300.0):

    z = lambda f: L_nm * (Langevin(b_nm * f / kT) + f / S_pN)
    f_pN = inversefunc(z, y_values=z_nm, domain=[(1e-7), 15000], image=[0, 15000])

    return f_pN


def gFJC(z_nm, L_nm=6.8, b_nm=0.01, S_pN=6300.0):
    """

    Parameters
    ----------
    z_nm:       distance
    k_pN__nm:   stiffness(pN/nm)
    L_nm:       contour length
    b_nm:       Kuhnlength
    S_pN:       Stretch modulus (pN)
    EFJC:       extended freely jointed chain

    Returns
    -------

    """
    z = lambda f: L_nm * (Langevin(b_nm * f / kT) + f / S_pN)
    f_pN = inversefunc(z, y_values=z_nm, domain=[(1e-7), 15000], image=[0, 15000])

    print('z_nm: ', z_nm)
    print('f_pn: ', f_pN)
    g_pNnm = -(kT * L_nm / b_nm) * (np.log((b_nm * f_pN) / (kT)) - np.log(
        np.sinh((b_nm * f_pN) / (kT)))) + L_nm * f_pN ** 2 / (2 * S_pN)

    # Remove offset at f = 0
    f_pN = 1e-9
    g_pNnm -= -(kT * L_nm / b_nm) * (np.log((b_nm * f_pN) / (kT)) - np.log(
        np.sinh((b_nm * f_pN) / (kT)))) + L_nm * f_pN ** 2 / (2 * S_pN)

    return g_pNnm / kT


def score_tails(moving_bp, fiber_start, dyads, dna, nucl):
    left_dyad = np.argmax(dyads > moving_bp) - 1
    right_dyad = left_dyad + fiber_start
    g = 0

    if 0 <= left_dyad < len(dyads) - fiber_start:
        t_up, t_down = tail_dist(left_dyad, right_dyad, dyads, dna, nucl)
        g += gFJC(t_up)
        print('gFJC t up: ', g)
        g += gFJC(t_down)

    if fiber_start is 2 and left_dyad >= 1:
        left_dyad -= 1
        right_dyad -= 1
        t_up, t_down = tail_dist(left_dyad, right_dyad, dyads, dna, nucl)
        g += gFJC(t_up)
        g += gFJC(t_down)
    print('g: ', g)
    return g


def dist_cms(nuc_1, nuc_2, dna,dyads,nucl):
    """

    Parameters
    ----------
    dna:    helixpose
    dyads:  indices of all dyads
    nucl:   nucleosomepose

    Returns
    -------
    distance between center of masses of nucleosomes
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


def params2coords(params):

    dr, frames = util.params2data(params)
    coords = [util.dr2coords(dr)]

    return coords


def coord_mean(filename, dyads, nucl):
    # get list of npz files in filename folder
    npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))
    params = []
    # params contains all 6 parameters for each basepair out of every npz file
    for f, file in enumerate(npz_f):
        data = np.load(file)
        params.append(data['params'])

    # get mean value of every parameter for each basepair
    params_m = np.mean(params,axis=0)

    # get cms of nucleosomes
    dr, frame = util.params2data(params_m)
    coords = params2coords(params_m)

    dyad_coord = coords[0][dyads[0]]
    dyad_frame = frame[dyads[0]]

    dyad_of = np.concatenate(([dyad_coord], dyad_frame.T+dyad_coord), axis=0)
    # target = np.concatenate(([dyad_coord], dyad_frame.T+dyad_coord), axis=0)
    # tf = nMC.get_transformation(dyad_of, target)
    #
    # ori = np.mean(coords[0][nucl.fixed[4:-4]], axis=0)
    # Nx = dyad_coord - ori
    # Nx = Nx / np.linalg.norm(Nx)
    # Nz = np.mean(coords[0][nucl.fixed[:7], :], axis=0) - np.mean(coords[0][nucl.fixed[7:], :], axis=0)
    # Ny = np.cross(Nx, Nz)
    # Ny = Ny / np.linalg.norm(Nz)
    # Nz = np.cross(Nx, Ny)
    # Nz = Nz / np.linalg.norm(Nz)
    # frame = np.array([Nx, Ny, Nz])
    # nucl_of = nMC.join_o_f(ori, np.transpose(frame))
    # nuc_cms = nMC.apply_transformation(nucl_of, tf)

    f_coord = []
    tf = nMC.get_transformation(dyad_of,
                                target=np.asarray([[0, 0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]]))
    for c in coords[0]:
        f_coord.append(nMC.apply_transformation(c, tf))

    # save dna coords to s_coord
    s_coord = f_coord[0]
    df = pd.DataFrame(np.array(s_coord) / 10)
    df.to_excel(fileio.change_extension(filename, 'coord.xlsx'), index=False, header=True)

    # nuc_cms_c = []
    # for n, cms in enumerate(nuc_cms):
    #     nuc_cms_c.append(nMC.apply_transformation(nuc_cms[n][0], tf)[0])
    #
    # dfc = pd.DataFrame(np.array(nuc_cms_c) / 10)
    # dfc.to_excel(fileio.change_extension(filename, 'coord_cms.xlsx'), index=False, header=True)

    return f_coord

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