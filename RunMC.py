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
np.set_printoptions(formatter={'float': '{: 0.3f}, '.format})


def score_dna(dna_params, p0, k, start_bp=0, end_bp=None):
    if end_bp is None:
        end_bp = len(dna_params) - 1
    g_dna = 0
    for i in range(start_bp, end_bp, 1):
        g_dna += np.sum(0.5 * (dna_params[i] - p0) * np.dot(k, dna_params[i] - p0)) - 3
    return g_dna * kT


def score_work(dna_coords, force, start_bp=0, end_bp=None):
    if end_bp is None:
        end_bp = len(dna_coords) - 1
    g_work = -(dna_coords[end_bp, 2] - dna_coords[start_bp, 2]) * force
    return g_work


def get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT):
    sigma_trans = 1  # 1 Angstrom, similar to basepair steps
    sigma_rot = 5 * np.pi / 180  # 5 degree, similar to basepair steps
    k1 = kT / np.asarray([sigma_trans, sigma_trans, sigma_trans, sigma_rot, sigma_rot, sigma_rot]) ** 2
    k1 *= 10
    G_unwrap = 0.5 * k1 * (wrap_params - fixed_wrap_params) ** 2
    G_unwrap = np.clip(G_unwrap, 0, e_wrap_kT * kT / 6)
    return np.sum(G_unwrap), G_unwrap


def score_wrapping(moving_bp, dna_coords, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT, half_nuc=False):
    closest_dyad = dyads[np.abs(dyads - moving_bp).argmin()]
    start_bp = closest_dyad - nucl.dyad
    end_bp = start_bp + len(nucl.dna.params)

    if start_bp <= moving_bp < end_bp:
        wrap_params = nMC.get_wrap_param(dna_coords, dna_frames, closest_dyad, nucl.fixed)
        g, g_wrap_all = get_unwrap_energy(wrap_params, fixed_wrap_params, e_wrap_kT)
        if half_nuc is True:
            if moving_bp < closest_dyad:
                g = np.sum(g_wrap_all[nucl.fixed < 0])
            else:
                g = np.sum(g_wrap_all[nucl.fixed > 0])
    else:
        g = 0
    return g


def score_stacking(moving_bp, coords, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start=1):
    start_dyad = np.argmax(dyads > moving_bp) - 1
    if 0 <= start_dyad < len(dyads) - fiber_start:
        # print(dyads[start_dyad], dyads[start_dyad + fiber_start], -dyads[start_dyad]+ dyads[start_dyad + fiber_start])
        stack_params = fMC.get_stack_pars(coords, frames, dyads[start_dyad], dyads[start_dyad + fiber_start],
                                          nucl)
        sigma = np.asarray([10.0, 1.0, 1.0, 0.1, 0.1, 0.1])
        k = kT / sigma ** 2

        g = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g = np.clip(g, 0, e_stack_kT * kT)

        # second potential
        sigma = np.asarray([10, 10, 10, 1e3, 1e3, 1e3])
        sigma *= 1
        k = kT / sigma ** 2
        g2 = 0.5 * np.sum(k * (stack_params - fixed_stack_params) ** 2) / kT
        g += g2
    else:
        g = 0
    return g


def score_surface(dna_coords):
    surface = dna_coords[:, 2] < 0
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
    dna_coords = dna.coords
    dna_params = dna.params
    dna_frames = dna.frames
    e_nucl = score_dna(nucl.dna.params, p0, k)
    g_wrap = 0
    g_dna = 0
    g_stack = 0
    g_work = 0
    n_nucs = len(dyads)
    for dyad1, dyad2 in zip(dyads[:-1], dyads[1:]):
        g_dna += score_dna(dna_params, p0, k, start_bp=dyad1, end_bp=dyad2)
        g_stack += score_stacking(dyad1 + 1, dna_coords, dna_frames, dyads, fixed_stack_params, e_stack_kT, nucl,
                                  fiber_start)
        g_work += score_work(dna_coords, force, start_bp=dyad1, end_bp=dyad2, )
        g_wrap += score_wrapping(dyad1 + 2, dna_coords, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)
        g_wrap += score_wrapping(dyad2 - 2, dna_coords, dna_frames, dyads, nucl, fixed_wrap_params, e_wrap_kT,
                                 half_nuc=True)

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


def MC_move(dna, bp, previous_bp, force, fixed_wrap_params, fixed_stack_params, dyads, nucl,
            random_step, e_wrap_kT, e_stack_kT, fiber_start):
    old_step_params = dna.params[bp]
    new_step_params = get_new_step_params(bp, previous_bp, dna, dyads, nucl, random_step)
    if np.array_equal(old_step_params, new_step_params):
        return False
    else:
        coords = dna.coords
        frames = dna.frames
        old_score = [
            score_wrapping(bp, coords, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT),
            score_stacking(bp, coords, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start),
            score_work(coords, force),
            score_surface(coords),
            0]
        dna.update(bp, new_step_params)
        coords = dna.coords
        frames = dna.frames
        new_score = [
            score_wrapping(bp, coords, frames, dyads, nucl, fixed_wrap_params, e_wrap_kT),
            score_stacking(bp, coords, frames, dyads, fixed_stack_params, e_stack_kT, nucl, fiber_start),
            score_work(coords, force),
            score_surface(coords),
            0]
    if util.MC_acpt_rej(np.sum(old_score), np.sum(new_score)):
        return True
    else:
        dna.update(bp, old_step_params)
        return False


def main(pars):
    # Setup files and forces
    root = '{0}st{1}x{2}'.format(pars['fiber_start'].value, pars['n_nuc'].value, pars['NRL'].value)
    filename = fileio.get_filename(incr=True, root=root)
    n_force_steps = 25000
    n_samples = 250

    fmax_pN = 10
    fmin_pN = 0.1
    forces = np.linspace(fmin_pN, fmax_pN, n_force_steps / 2)

    sample_forces = np.logspace(np.log10(fmin_pN), np.log10(fmax_pN), n_samples / 2)
    sample_indices = np.searchsorted(forces, sample_forces)
    sample_indices = np.append(sample_indices, n_force_steps / 2 + (n_force_steps / 2 - sample_indices[::-1]) - 1)

    forces = np.append(forces, forces[::-1])
    duty_cycle = 1  # only Monte Carlo dna handles in 1 out of [dutycyle] forcesteps

    # Initialize fiber pose
    dna, dyads, nucl = fMC.create_unfolded_fiber(fiber_pars=pars)
    # Get from file
    if False:
        datafile = 'E:\\Users\\noort\\data\\20180513\\2x197_006\\2x197_006_0001.npz'
        dna = HelixPose.from_file(fileio.change_extension(datafile, 'npz'))

    # fileio.plot_dna(dna, title='Initial conformation\n', range_nm=100, save=False, wait=10)

    pars.pretty_print(columns=['value'])
    print('>>> Current file: {}'.format(filename))

    # Get stack and wrap parameters
    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coords, nucl.dna.frames, nucl.dyad, nucl.fixed)
    e_wrap_kT = pars['e_wrap_kT'].value

    casted_fiber, _, _ = fMC.create_casted_fiber(pars, nucl)
    # fileio.plot_dna(casted_fiber, range_nm=50, wait=1, origin_index=dyads[0])
    fixed_stack_params = fMC.get_stack_pars(casted_fiber.coords, casted_fiber.frames, dyads[0], dyads[1], nucl)
    fiber_start = pars['fiber_start'].value

    # Initialize random steps
    random_step = RandomStepSimple.load_gaussian_params(dna_step_file)
    p0 = np.load(dna_step_file)[0]
    k = np.linalg.inv(np.load(dna_step_file)[1:])

    basepairs = np.asarray(range(pars['L_bp'] - 1))
    fiber_basepairs = [dyads[0] + nucl.fixed[0], dyads[-1] + nucl.fixed[-1]]

    if pars['fiber_start'].value == 1:
        e_stack_kT = 1e6
    else:
        e_stack_kT = pars['e_stack_kT'].value
    g_nuc_kT_all = []

    pars['F_pN'].value = 0
    pars['z_nm'].value = dna.coord_terminal[2] / 10
    fileio.write_xlsx_row(fileio.get_filename(sub=True, incr=True, ext='npz'), -1, pars, report_file=filename)
    dna.write2disk(fileio.get_filename(sub=True, ext='npz'))

    previous_bp = 0
    fileio.report_progress(n_force_steps, title='RunMC', init=True)
    for i, force in enumerate(forces):
        if i % duty_cycle == 0:
            iter_range = [0, len(dna.params) - 1]
        else:
            iter_range = fiber_basepairs

        for bp in basepairs:
            MC_move(dna, bp, previous_bp, force, fixed_wrap_params, fixed_stack_params,
                    dyads, nucl, random_step, e_wrap_kT, e_stack_kT, fiber_start)
            previous_bp = bp
        basepairs = basepairs[::-1]

        g_nuc_kT, names = get_nuc_energies(dna, fixed_wrap_params, fixed_stack_params, dyads, nucl, e_wrap_kT,
                                           e_stack_kT, fiber_start, p0, k, force)
        g_nuc_kT_all.append(g_nuc_kT)

        fileio.report_progress(i, title='Force = {0:.1f} pN, g_stack = {1:.1f} kT'.format(force, g_nuc_kT[2]))

        if i == 100:
            e_stack_kT = pars['e_stack_kT'].value

        if i in sample_indices:
            # fileio.plot_dna(dna, update=True, title='F = {:.1f} pN\n'.format(force), save=True)
            g_nuc_kT_all = np.mean(g_nuc_kT_all, axis=0)
            for g, name in zip(g_nuc_kT_all, names):
                pars[name].value = g
            g_nuc_kT_all = []
            pars['F_pN'].value = force
            pars['z_nm'].value = dna.coord_terminal[2] / 10
            fileio.write_xlsx_row(fileio.get_filename(sub=True, incr=True, ext='npz'), i, pars, report_file=filename)
            dna.write2disk(fileio.get_filename(sub=True, ext='npz'))

    aMC.plot_fz(fileio.change_extension(filename, 'xlsx'))
    aMC.plot_gz(fileio.change_extension(filename, 'xlsx'))
    try:
        fileio.create_pov_movie(fileio.get_filename(sub=True, folder=True), origin_frame=0, reverse=False)
    except Exception as e:
        print(Exception, e)
    return


if __name__ == '__main__':
    # Parameters for reporting results
    pars = Parameters()
    pars.add('F_pN', value=0)
    pars.add('z_nm', value=0)
    pars.add('g_dna_kT', value=0)
    pars.add('g_wrap_kT', value=0)
    pars.add('g_stack_kT', value=0)
    pars.add('g_work_kT', value=0)

    # Parameters that define the nucleosomal array
    pars.add('L_bp', value=1000)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=4)

    # Parameters that define the folded fiber
    pars.add('rise_A', value=50)
    pars.add('nld_A', value=20)
    pars.add('chirality', value=-1)
    pars.add('face', value=1)

    # Parameters that are typically varied between simulations
    pars.add('diameter_A', value=330)
    pars.add('NRL', value=167)

    if pars['NRL'].value == 167:
        pars.add('fiber_start', value=2)
    else:
        pars.add('fiber_start', value=1)
    pars.add('Unwrapped_bp', value=0)
    pars.add('e_wrap_kT', value=3.0)
    pars.add('e_stack_kT', value=20)

    main(pars)
