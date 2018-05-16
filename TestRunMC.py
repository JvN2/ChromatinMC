# -*- coding: utf-8 -*-
"""
@author: John van Noort
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
from helixmc.random_step import RandomStepSimple
from helixmc.pose import HelixPose
# ChromatinMC modules:
import FiberMC as fMC
import NucleosomeMC as nMC
import FileIO as fileio
import RunMC as rMC

dna_step_file = 'C:\\Python27\\Lib\\site-packages\\helixmc\\data\\DNA_gau.npy'
dna_pose_file = 'E:\\Users\\noort\\data\\20180513\\2x197_006\\2x197_006_0001.npz'
kT = 41.0
np.set_printoptions(formatter={'float': '{: 0.4f}, '.format})


def main():
    pars = Parameters()
    pars.add('F_pN', value=0)
    pars.add('z_nm', value=0)
    pars.add('L_bp', value=400)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=2)
    pars.add('NRL', value=197)
    pars.add('fiber_start', value=1)
    pars.add('e_wrap_kT', value=3)
    pars.add('e_stack_kT', value=25)

    e_wrap_kT = pars['e_wrap_kT'].value
    e_stack_kT = pars['e_stack_kT'].value
    fiber_start = pars['fiber_start'].value

    n_steps = 1000

    dna, dyads, nucl = fMC.create_unfolded_fiber(fiber_pars=pars)
    dna = HelixPose.from_file(dna_pose_file)

    fixed_wrap_params = nMC.get_wrap_param(nucl.dna.coords, nucl.dna.frames, nucl.dyad, nucl.fixed)
    fixed_stack_params = fMC.get_stack_pars(dna.coords, dna.frames, dyads[0], dyads[fiber_start], nucl)

    print(fixed_stack_params)

    random_step = RandomStepSimple.load_gaussian_params(dna_step_file)
    basepairs = range(len(dna.params))
    previous_bp = 0
    force = 10.0

    fileio.plot_dna(dna, title='{1}: g_stack = {0:5.1f} kT'.format(0, 0), wait=1, range_nm=50, origin_index=199)
    fileio.report_progress(n_steps, title='TestRunMC', init=True)
    for i in range(n_steps):
        for bp in basepairs:
            rMC.MC_move(dna, bp, previous_bp, force, fixed_wrap_params, fixed_stack_params,
                        dyads, nucl, random_step, e_wrap_kT, e_stack_kT, fiber_start)
            previous_bp = bp
        basepairs = basepairs[::-1]

        if i % 10 == 0:
            g_stack_kt = rMC.score_stacking(dyads[0] + 1, dna.coords, dna.frames, dyads, fixed_stack_params, e_stack_kT,
                                            fiber_start)/kT
            fileio.report_progress(i, title='g_stack = {0:.1f} kT'.format(g_stack_kt))
            fileio.plot_dna(dna, update=True, title='{1}: g_stack = {0:5.1f} kT'.format(g_stack_kt, i), origin_index=199)

    return


if __name__ == '__main__':
    main()
