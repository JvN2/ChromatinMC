# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 10:07:46 2018

@author: noort
"""

from __future__ import print_function
import matplotlib as mpl

try:
    mpl.use(u'TkAgg')
    mpl.interactive(False)
except:
    pass

import numpy as np
from lmfit import Parameters
from helixmc import util
from helixmc.random_step import RandomStepSimple, RandomStepAgg, symmetrize_WC
from helixmc.score import ScoreTweezers
from helixmc.pose import HelixPose
# ChromatinMC modules:
import FiberMC as fMC
import NucleosomeMC as nMC
import analyzeMC as aMC
import FileIO as fileio

dna_step_file = 'C:\\Python27\\Lib\\site-packages\\helixmc\\data\\DNA_gau.npy'
kT = 41.0
np.set_printoptions(formatter={'float': '{: 0.1f}, '.format})


def score_dna(start_bp, end_bp, dna, p0, k):
    g_dna = 0
    for i in range(start_bp, end_bp, 1):
        g_dna += np.sum(0.5 * (dna.params[i] - p0) * np.dot(k, dna.params[i] - p0)) - 3
    return g_dna * kT


def score_work(start_bp, end_bp, dna, force):
    g_work = -(dna.coords[end_bp, 2] - dna.coords[start_bp, 2]) * force
    return g_work


def get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT):
    sigma_trans = 1  # 1 Angstrom, similar to basepair steps
    sigma_rot = 5 * np.pi / 180  # 5 degree, similar to basepair steps
    k1 = kT / np.asarray([sigma_trans, sigma_trans, sigma_trans, sigma_rot, sigma_rot, sigma_rot]) ** 2
    k1 *= 10
    G_unwrap = 0.5 * k1 * (wrap_params - fixed_wrap_params) ** 2
    G_unwrap = np.clip(G_unwrap, 0, e_wrap_kT * kT / 6)
    return np.sum(G_unwrap), G_unwrap


def score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, e_wrap_kT, half_nuc=False):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    start_bp = closest_dyad - nucl.dyad
    end_bp = start_bp + len(nucl.dna.params)
    if start_bp <= moving_bp < end_bp:
        wrap_params = nMC.get_wrap_params(dna, closest_dyad, nucl.fixed)
        G, Gwrap_all = get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT)
        if half_nuc is True:
            if moving_bp < closest_dyad:
                G = np.sum(Gwrap_all[nucl.fixed < 0])
            else:
                G = np.sum(Gwrap_all[nucl.fixed > 0])
    else:
        G = 0
    return G


def score_stacking(moving_bp, dna, dyads, fixed_stack_params, e_stack_kT, fiber_start=1):
    start_dyad = np.argmax(dyads > moving_bp) - 1
    if 0 <= start_dyad < len(dyads) - fiber_start:
        stack_params = fMC.get_stack_pars(dna, [dyads[start_dyad], dyads[start_dyad + fiber_start]])
        sigma = np.asarray([1, 1, 1, 0.1, 0.1, 0.1])
        sigma /= 1
        k = kT / sigma ** 2
        g = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g = np.clip(g, 0, e_stack_kT * kT)
        # second potential
        sigma = np.asarray([10, 10, 10, 1e3, 1e3, 1e3])
        k = kT / sigma ** 2
        g2 = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g2 = np.clip(g2, 0, e_stack_kT * kT)
        g += g2
    else:
        g = 0
    return g


def score_surface(dna):
    surface = dna.coords[:, 2] < 0
    if np.sum(surface) > 0:
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


def get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl, e_wrap_kT, e_stack_kT, \
                     fiber_start, p0, k, force):
    e_nucl = score_dna(0, len(nucl.dna.params), nucl.dna, p0, k)
    g_wrap = 0
    g_dna = 0
    g_stack = 0
    g_work = 0
    n_nucs = len(dyads)
    for dyad1, dyad2 in zip(dyads[:-1], dyads[1:]):
        g_dna += score_dna(dyad1, dyad2, dna, p0, k)
        g_stack += score_stacking(dyad1 + 1, dna, dyads, fixed_stack_params, e_stack_kT, fiber_start)
        g_work += score_work(dyad1, dyad2, dna, force)
        g_wrap += score_wrapping(dyad1 + 2, dna, dyads, nucl, fixed_wrap_params, e_wrap_kT, half_nuc=True)
        g_wrap += score_wrapping(dyad2 - 2, dna, dyads, nucl, fixed_wrap_params, e_wrap_kT, half_nuc=True)

    g_wrap /= (n_nucs - 1)
    g_dna /= (n_nucs - 1)
    g_dna -= e_nucl
    g_work /= (n_nucs - 1)
    if fiber_start == 1:
        g_stack /= (n_nucs - 1)
    else:
        g_stack /= (n_nucs - 2)
    g_nuc_kT = np.asarray([g_dna, g_wrap, g_stack, g_work]) / kT
    names = ['g_dna_kT', 'g_wrap_kT', 'g_stack_kT', 'g_work_kT']
    return g_nuc_kT, names


from numba import jit


@jit
def MC_moves(dna, basepairs, scorefxn, fixed_wrap_params, fixed_stack_params, dyads, nucl,
             random_step, e_wrap_kT, e_stack_kT, fiber_start):
    previous_bp = 0
    old_dna_score = scorefxn(dna) + score_surface(dna)
    for moving_bp in basepairs:
        old_step_params = dna.params[moving_bp]
        new_step_params = get_new_step_params(moving_bp, previous_bp, dna, dyads, nucl, random_step)
        if not np.array_equal(old_step_params, new_step_params):
            old_score = old_dna_score \
                        + score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, e_wrap_kT) \
                        + score_stacking(moving_bp, dna, dyads, fixed_stack_params, e_stack_kT, fiber_start)
            dna.update(moving_bp, new_step_params)
            new_dna_score = scorefxn(dna) + score_surface(dna)
            new_score = new_dna_score \
                        + score_wrapping(moving_bp, dna, dyads, nucl, fixed_wrap_params, e_wrap_kT) \
                        + score_stacking(moving_bp, dna, dyads, fixed_stack_params, e_stack_kT, fiber_start)
            if util.MC_acpt_rej(old_score, new_score):
                old_dna_score = new_dna_score
            else:
                dna.update(moving_bp, old_step_params)
        previous_bp = moving_bp
    return


def main():
    # Params for reporting results
    pars = Parameters()
    pars.add('F_pN', value=0)
    pars.add('z_nm', value=0)

    pars.add('g_dna_kT', value=0)
    pars.add('g_wrap_kT', value=0)
    pars.add('g_stack_kT', value=0)
    pars.add('g_work_kT', value=0)

    # Params that define the nucleosomal array
    pars.add('L_bp', value=1000)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=4)
    pars.add('dyad0_bp', value=0)
    # Params that define the folded fiber
    pars.add('diameter_A', value=330)
    pars.add('rise_A', value=100)
    pars.add('nld_A', value=25)
    pars.add('chirality', value=-1)
    pars.add('face', value=1)

    # Params that are typically varied between simulations
    pars.add('NRL', value=197)
    pars.add('fiber_start', value=1)
    pars.add('Unwrapped_bp', value=30)
    pars.add('e_wrap_kT', value=3)
    pars.add('e_stack_kT', value=23)

    e_wrap_kT = pars['e_wrap_kT'].value

    # Setup files and forces
    filename = fileio.get_filename(incr=True, root='4x197')
    n_force_steps = 500
    n_samples = 250

    fmax_pN = 10
    fmin_pN = 0.1
    forces = np.linspace(fmin_pN, fmax_pN, n_force_steps / 2)

    sample_forces = np.logspace(np.log10(fmin_pN), np.log10(fmax_pN), n_samples / 2)
    sample_indices = (np.searchsorted(forces, sample_forces))
    sample_indices = (np.append(sample_indices, n_force_steps / 2 + (n_force_steps / 2 - sample_indices[::-1]) - 1))

    forces = np.append(forces, forces[::-1])

    get_from_file = False
    if get_from_file:
        # Get from file
        file_in = 'E:\\Users\\noort\\data\\20180321\\data_003.xlsx'
        dataset = 0
        pars, datafile = fileio.read_xlsx(file_in, dataset, pars=pars)
        dna, dyads, nucl = fMC.create_nuc_array(p=pars)
        dna = HelixPose.from_file(fileio.change_extension(datafile, 'npz'))
    else:
        # Initialize fiber pose
        dna, dyads, nucl = fMC.create_nuc_array(p=pars)

    pars['dyad0_bp'].value = dyads[0]
    # fileio.plot_dna(dna, title='Initial conformation\n', range_nm=100, save=True)
    dna.write2disk(fileio.get_filename(sub=True, ext='npz'))

    pars.pretty_print(columns=['value'])
    print('>>> Current file: {}'.format(filename))

    # Get stack and wrap parameters
    fixed_wrap_params = nMC.get_wrap_params(nucl.dna, nucl.dyad, nucl.fixed)
    fiber_dna, dyads, w = fMC.create_folded_fiber(pars, nucl)
    fixed_stack_params = fMC.get_stack_pars(fiber_dna, dyads)[0]
    fiber_start = pars['fiber_start'].value

    # Initialize random steps
    random_step = RandomStepSimple.load_gaussian_params(dna_step_file)
    p0 = np.load(dna_step_file)[0]
    k = np.linalg.inv(np.load(dna_step_file)[1:])

    basepairs = np.asarray(range(pars['L_bp'] - 1))
    accept = 0
    all_coord = np.empty((n_samples, 3))

    current_step = 0
    e_stack_kT = 1e6
    g_nuc_kT_all = []
    fileio.report_progress(n_force_steps, title='RunMC3', init=True)
    for i, force in enumerate(forces):
        fileio.report_progress(i, title='Force = {:.1f} pN'.format(force))
        scorefxn = ScoreTweezers(force)
        MC_moves(dna, basepairs, scorefxn, fixed_wrap_params, fixed_stack_params,
                 dyads, nucl, random_step, e_wrap_kT, e_stack_kT, fiber_start)
        basepairs = basepairs[::-1]

        g_nuc_kT, names = get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl, e_wrap_kT,
                                           e_stack_kT, fiber_start, p0, k, force)
        g_nuc_kT_all.append(g_nuc_kT)
        if i == 10:
            e_stack_kT = pars['e_stack_kT'].value

        if i in sample_indices:
            # fileio.plot_dna(dna, update=True, title='F = {:.1f} pN\n'.format(force), save=True)
            g_nuc_kT_all = np.sum(g_nuc_kT_all, axis=0)
            for g, name in zip(g_nuc_kT_all, names):
                pars[name].value = g
            g_nuc_kT_all = []
            pars['F_pN'].value = force
            pars['z_nm'].value = dna.coord_terminal[2] / 10
            fileio.write_xlsx_row(fileio.get_filename(sub=True, incr=True, ext='npz'), i, pars, report_file=filename)
            dna.write2disk(fileio.get_filename(sub=True, ext='npz'))
            all_coord[current_step] = dna.coord_terminal

    aMC.plot_fz(filename)
    try:
        fileio.create_pov_movie(fileio.get_filename(sub=True, folder=True), origin_frame=0, reverse=False)
    except Exception as e:
        print(Exception, e)


if __name__ == '__main__':
    main()
# Github works!!!!
