# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 10:07:46 2018

@author: noort
"""

from __future__ import print_function
import warnings

warnings.filterwarnings("ignore")

import matplotlib as mpl
import math
import pandas as pd

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass

import numpy as np
import os as os
import re
from lmfit import Parameters
from helixmc import util
from helixmc.random_step import RandomStepSimple, RandomStepAgg, symmetrize_WC
from helixmc.score import ScoreTweezers
from helixmc.pose import HelixPose
from helixmc.util import locate_data_file
# ChromatinMC modules:
import FiberMC as fMC
import NucleosomeMC as nMC
import analyzeMC as aMC
import FileIO as fileio
import Tails as tMC

dna_step_file = locate_data_file('DNA_gau.npy')
kT = 41.0
np.set_printoptions(formatter={'float': '{: 0.3f}, '.format})


def score_dna(dna_params, p0, k, start_bp=0, end_bp=None, w=None):
    if end_bp is None:
        end_bp = len(dna_params) - 1
    if w is None:
        w = np.ones(len(dna_params))
    g = []
    for i in range(start_bp, end_bp, 1):
        g.append(w[i] * ((0.5 * (dna_params[i] - p0) * np.dot(k, dna_params[i] - p0)) - 0.5))
    return np.asarray(g)


def score_work(dna_coord, force, start_bp=0, end_bp=None):
    if end_bp is None:
        end_bp = len(dna_coord) - 1
    g_work = -(dna_coord[end_bp, 2] - dna_coord[start_bp, 2]) * force
    return g_work


def score_exclusion(dna_coord, dna_frames, dyads, nucl):
    nuc_cms = []
    for dyad in dyads:
        nuc_cms.append(nMC.get_nuc_of(dna_coord, dna_frames, dyad, nucl)[0])
    r_excl = 55
    g_excl = 0
    for i, cm1 in enumerate(nuc_cms[:-1]):
        for cm2 in nuc_cms[i + 1:]:
            if np.sum((cm2 - cm1) ** 2) < r_excl ** 2:
                g_excl += 1e7
    return g_excl


def get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT):
    sigma_trans = 1  # 1 Angstrom, similar to basepair steps
    sigma_rot = 5 * np.pi / 180  # 5 degree, similar to basepair steps
    k1 = kT / np.asarray([sigma_trans, sigma_trans, sigma_trans, sigma_rot, sigma_rot, sigma_rot]) ** 2
    k1 *= 10
    G_unwrap = 0.5 * k1 * (wrap_params - fixed_wrap_params) ** 2

    G_unwrap = np.clip(G_unwrap, 0, e_wrap_kT * kT / 6)

    # G_unwrap[0] = np.clip(G_unwrap[0], 0, 0.5 * e_wrap_kT * kT / 6)
    # G_unwrap[-1] = np.clip(G_unwrap[-1], 0, 0.5 * e_wrap_kT * kT / 6)
    #
    return np.sum(G_unwrap), G_unwrap


def score_wrapping(moving_bp, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT, half_nuc=False):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    start_bp = closest_dyad - nucl.dyad
    end_bp = start_bp + len(nucl.dna.params)

    if start_bp <= moving_bp < end_bp:
        wrap_params = nMC.get_wrap_param(dna_coord, dna_frames, closest_dyad, nucl.fixed)
        g, g_wrap_all = get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT)

        if e_wrap_kT is 0:
            fixed_bps = [closest_dyad, closest_dyad]
        else:
            fixed_bps = nucl.fixed[np.sum(g_wrap_all, axis=1) < 0.99 * e_wrap_kT * kT]
            if len(fixed_bps) > 0:
                fixed_bps = np.asarray([min(fixed_bps), max(fixed_bps)] + closest_dyad)

        if half_nuc is True:
            if moving_bp < closest_dyad:
                g = np.sum(g_wrap_all[nucl.fixed < 0])
            else:
                g = np.sum(g_wrap_all[nucl.fixed > 0])
    else:
        g = 0
        fixed_bps = []
    return g, fixed_bps


def score_stacking(moving_bp, coord, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start=1):
    left_dyad = np.argmax(dyads > moving_bp) - 1
    right_dyad = left_dyad + fiber_start
    g_min = 0

    sigma = np.asarray([1.0, 1.0, 1.0, 0.1, 0.1, 0.1])
    # sigma *= 2.0
    k = kT / sigma ** 2

    if 0 <= left_dyad < len(dyads) - fiber_start:
        stack_params = fMC.get_stack_pars(coord, frames, dyads[left_dyad], dyads[right_dyad],
                                          nucl, fiber_start)
        g = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g += 0.01 * g * g
        g_min += np.clip(g, 0, e_stack_kT * kT)

    if fiber_start is 2 and left_dyad >= 1:
        stack_params = fMC.get_stack_pars(coord, frames, dyads[left_dyad - 1], dyads[right_dyad - 1],
                                          nucl, fiber_start)
        g = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g += 0.01 * g * g
        g_min += np.clip(g, 0, e_stack_kT * kT)

    return g_min


def score_surface(dna_coord):
    surface = dna_coord[:, 2] < 0
    bead = dna_coord[:, 2] > dna_coord[-1, 2]
    if np.sum(surface + bead) > 0:
        return 1e7
    else:
        return 0


def get_new_step_params(moving_bp, prev_bp, dna, dyads, nucl, random_step):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    nucl_bp = moving_bp - closest_dyad + nucl.dyad
    nucl_prev_bp = prev_bp - closest_dyad + nucl.dyad
    if 1 <= nucl_bp < len(nucl.dna.params) - 1:
        if np.array_equal(dna.params[prev_bp], nucl.dna.params[nucl_prev_bp]):
            new_step_params = nucl.dna.params[nucl_bp]
        else:
            new_step_params = random_step()[0]
    else:
        new_step_params = random_step()[0]
    return new_step_params


def get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl, e_wrap_kT, e_stack_kT, e_nuc_kT,
                     fiber_start, p0, k, force):
    dna_coord = dna.coord
    dna_params = dna.params
    dna_frames = dna.frames
    g_wrap = 0
    g_dna = np.zeros(6)
    g_stack = 0
    g_work = 0
    n_nucs = len(dyads)

    w = np.ones(len(dna_params))
    for dyad in dyads:
        fixed_bps = score_wrapping(dyad + 1, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                   half_nuc=True)[1]
        if len(fixed_bps) > 0:
            w[fixed_bps[0]:fixed_bps[1]] = 0
        else:
            w = None
    g_dna_all = score_dna(dna_params, p0, k, w=w)

    for dyad1, dyad2 in zip(dyads[:-1], dyads[1:]):
        g_dna += np.sum(g_dna_all[dyad1:dyad2], axis=0)
        g_stack += score_stacking(dyad1 + 1, dna_coord, dna_frames, dyads, fixed_stack_params, e_stack_kT, nucl,
                                  fiber_start)
        g_work += score_work(dna_coord, force, start_bp=dyad1, end_bp=dyad2, )
        g_wrap += score_wrapping(dyad1 + 2, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)[0]
        g_wrap += score_wrapping(dyad2 - 2, dna_coord, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)[0]

    g_dna /= (n_nucs - 1)
    g_wrap /= (n_nucs - 1)
    g_work /= (n_nucs - 1)
    g_stack /= (n_nucs - fiber_start) / fiber_start

    g_nuc_kT = np.asarray([np.sum(g_dna) * kT, g_wrap, g_stack, g_work]) / kT
    names = ['g_dna_kT', 'g_wrap_kT', 'g_stack_kT', 'g_work_kT']

    g_nuc_kT = np.append(g_nuc_kT, g_dna)
    names += ['g_shift_kT', 'g_slide_kT', 'g_rise_kT', 'g_tilt_kT', 'g_roll_kT', 'g_twist_kT']
    return g_nuc_kT, names


def MC_move(dna, bp, previous_bp, force, fixed_wrap_params, fixed_stack_params, dyads, nucl,
            random_step, e_wrap_kT, e_stack_kT, fiber_start, Tail_switch = False):
    old_step_params = dna.params[bp]
    new_step_params = get_new_step_params(bp, previous_bp, dna, dyads, nucl, random_step)
    if np.array_equal(old_step_params, new_step_params):
        return False
    elif Tail_switch == False: # use old stacking
        coord = dna.coord
        frames = dna.frames
        old_score = [
            score_wrapping(bp, coord, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT)[0],
            score_stacking(bp, coord, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start),
            score_exclusion(coord, frames, dyads, nucl),
            score_work(coord, force),
            # score_surface(coord),
            0]
        dna.update(bp, new_step_params)
        coord = dna.coord
        frames = dna.frames
        new_score = [
            score_wrapping(bp, coord, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT)[0],
            score_stacking(bp, coord, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start),
            score_exclusion(coord, frames, dyads, nucl),
            score_work(coord, force),
            # score_surface(coord),
            0]
    else:
        coord = dna.coord
        frames = dna.frames
        old_score = [
            score_wrapping(bp, coord, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT)[0],
            tMC.score_tails(bp, fiber_start, dyads, dna, nucl),
            score_exclusion(coord, frames, dyads, nucl),
            score_work(coord, force),
            # score_surface(coord),
            0]
        dna.update(bp, new_step_params)
        coord = dna.coord
        frames = dna.frames
        new_score = [
            score_wrapping(bp, coord, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT)[0],
            tMC.score_tails(bp, fiber_start, dyads, dna, nucl),
            score_exclusion(coord, frames, dyads, nucl),
            score_work(coord, force),
            # score_surface(coord),
            0]
    if util.MC_acpt_rej(np.sum(old_score), np.sum(new_score)):
        return True
    else:
        dna.update(bp, old_step_params)
        return False


def main(n_steps, root):
    pars = Parameters()
    # Parameters that define the nucleosomal array
    pars.add('L_bp', value=428)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=4)
    pars.add('e_nuc_kT', value=34.7)

    # Parameters that define the folded fiber
    pars.add('rise_A', value=100)
    pars.add('nld_A', value=20)
    pars.add('chirality', value=1)
    pars.add('face', value=1)
    pars.add('diameter_A', value=330)

    pars.add('e_wrap_kT', value=2.1)
    pars.add('e_stack_kT', value=25)
    pars.add('NRL', value=187)
    pars.add('fiber_start', value=2)
    # Parameters for reporting results
    pars.add('F_pN', value=0)
    pars.add('z_nm', value=0)
    pars.add('g_dna_kT', value=0)
    pars.add('g_wrap_kT', value=0)
    pars.add('Unwrapped_bp', value=0)
    pars.add('g_stack_kT', value=0)
    pars.add('g_work_kT', value=0)
    pars.add('g_shift_kT', value=0)
    pars.add('g_slide_kT', value=0)
    pars.add('g_rise_kT', value=0)
    pars.add('g_tilt_kT', value=0)
    pars.add('g_roll_kT', value=0)
    pars.add('g_twist_kT', value=0)

    # Setup files and forces
    if root is None:
        root = '{1}x{2}x{0}s{3}w{4:0.1f}'.format(pars['fiber_start'].value, pars['n_nuc'].value, pars['NRL'].value,
                                                 pars['e_stack_kT'].value, pars['e_wrap_kT'].value).replace('.', '-')
    else:
        iterpar = []
        for par_txt in re.findall('\d+', root):
            iterpar.append(float(par_txt))
        iterpar[4] += iterpar[5] / 10.0

        pars['n_nuc'].value = iterpar[0]
        pars['NRL'].value = iterpar[1]
        pars['fiber_start'].value = int(iterpar[2])
        pars['e_stack_kT'].value = iterpar[3]
        pars['e_wrap_kT'].value = iterpar[4]

    # create optimal fiber length for each NRL, with 14 bp handles
    pars['L_bp'].value = int(pars['n_nuc'].value * pars['NRL'].value + 28)
    print('L_bp: ', pars['L_bp'].value)

    filename = fileio.get_filename(incr=True, root=root, ext='xlsx', )
    # print('\n>>> Current file: {}'.format(filename))
    n_samples = 250
    fmin_pN = 0
    fmax_pN = 0
    forces = np.linspace(fmin_pN, fmax_pN, n_steps / 2)

    sample_forces = np.logspace(np.log10(fmin_pN), np.log10(fmax_pN), n_samples / 2)
    sample_indices = np.searchsorted(forces, sample_forces)
    sample_indices = np.append(sample_indices, n_steps / 2 + (n_steps / 2 - sample_indices[::-1]) - 1)
    forces = np.append(forces, forces[::-1])

    dummy_steps = 10
    sample_indices += dummy_steps
    sample_indices[0] = 0
    forces = np.append(np.zeros(dummy_steps), forces)

    # Initialize random steps
    random_step = RandomStepSimple.load_gaussian_params(dna_step_file)
    p0 = np.load(dna_step_file)[0]
    k = np.linalg.inv(np.load(dna_step_file)[1:])

    # Initialize fiber pose
    dna, dyads, nucl = fMC.create_unfolded_fiber(fiber_pars=pars)
    pars.add('dyad0', value=dyads[0])
    e_nuc_kT = pars['e_nuc_kT'].value

    # Get from file
    if False:
        datafile = 'E:\\Users\\noort\\data\\20180513\\2x197_006\\2x197_006_0001.npz'
        dna = HelixPose.from_file(fileio.change_extension(datafile, 'npz'))

    # Get stack and wrap parameters
    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coord, nucl.dna.frames, nucl.dyad, nucl.fixed)
    e_wrap_kT = pars['e_wrap_kT'].value

    fiber_start = pars['fiber_start'].value

    n_ofs = fMC.get_casted_fiber_frames(pars)
    fixed_stack_params = nMC.ofs2params(n_ofs[fiber_start], n_ofs[0], _3dna=True)

    basepairs = np.asarray(range(pars['L_bp'] - 1))
    if pars['e_stack_kT'].value == 0:
        e_stack_kT = 0
    else:
        e_stack_kT = 1e3 * kT

    g_nuc_kT_all = []
    tails = []
    cms_dist = []
    # number of npz files that will be stored during simulation
    num_npz = 5
    idx = np.round(np.linspace(dummy_steps, len(forces) - 1, num_npz))
    # indices of nucleosomes of which distances will be calculated in tails and cms_dist
    nuc_1 = 0
    nuc_2 = 1
    Tail_switch = True # True: score tails, False: score_stacking


    pars['F_pN'].value = 0
    pars['z_nm'].value = dna.coord_terminal[2] / 10

    previous_bp = 0
    datafile = fileio.get_filename(sub=True, incr=True, ext='npz')
    paramsfile = fileio.change_extension(datafile, 'parms_0001.npz')

    fileio.report_progress(n_steps, title='RunMC', init=True)
    for i, force in enumerate(forces):
        if i == dummy_steps:
            e_stack_kT = pars['e_stack_kT'].value
            g_nuc_kT_all = []
            # Tail_switch = True

        g_nuc_kT, names = get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl, e_wrap_kT,
                                           e_stack_kT, e_nuc_kT, fiber_start, p0, k, force)
        g_nuc_kT_all.append(g_nuc_kT)

        tails.append(tMC.tail_dist(nuc_1, nuc_2, dyads, dna, nucl, orientation='-*'))
        cms_dist.append(tMC.dist_cms(nuc_1, nuc_2, dna, dyads, nucl))

        fileio.report_progress(i, title='Force = {0:.1f} pN {1}'.format(force, os.path.splitext(filename)[0]))

        if i in sample_indices:
            g_nuc_kT_all = np.mean(g_nuc_kT_all, axis=0)
            for g, name in zip(g_nuc_kT_all, names):
                pars[name].value = g
            pars['F_pN'].value = force
            pars['z_nm'].value = dna.coord_terminal[2] / 10
            fileio.write_xlsx_row(datafile, i - dummy_steps, pars, report_file=filename)
            # dna.write2disk(datafile)
            g_nuc_kT_all = []
            datafile = fileio.increment_file_nr(datafile)

        if i in idx:
            dna.write2disk(paramsfile)
            paramsfile = fileio.increment_file_nr(paramsfile)

        for bp in basepairs:
            MC_move(dna, bp, previous_bp, force, fixed_wrap_params, fixed_stack_params,
                    dyads, nucl, random_step, e_wrap_kT, e_stack_kT, fiber_start, Tail_switch)
            previous_bp = bp
        basepairs = basepairs[::-1]



    tMC.coord_mean(filename, dyads, nucl)

    coord, radius, colors = tMC.get_histones(dna.coord, dyads, nucl, dna=dna)

    print(fileio.create_pov(filename, coord, radius=radius, colors=colors, range_A=[750, 750], offset_A=[0, 0, 150],
                            show=True, width_pix=1500))
    #
    f_coord = tMC.origin(dna, dyads, nucl, coord, filename, axis=False)
    # # colors += 'z'
    # # radius = np.append(radius, 10)
    print(fileio.create_pov((fileio.change_extension(filename, '_org.png')), f_coord, radius=radius, colors=colors,
                            range_A=[750, 750], offset_A=[0, 0, 300], show=True, width_pix=1500))
    #
    tMC.dist_plot(filename, cms_dist[dummy_steps:], save=True)
    tMC.tail_plot(filename, tails[dummy_steps:], save=True)


    # aMC.plot_fz(filename)
    # aMC.plot_gi(filename, force_range=[0.1, 1.5])
    # fileio.create_pov_movie(filename, fps=5, octamers=True, overwrite=False, frame=[60, 0, -90])
    return


if __name__ == '__main__':
    # pars.pretty_print(columns=['value'])

    main(5e4, '8x197x1s21w2-1')
