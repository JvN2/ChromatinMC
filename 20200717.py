import NucleosomeMC as nmc
import TailsAnalyze as TA

TA.plot_npz(r"D:\users\vlaar\data\20201117\8x197x1s25w2-1_002")
# #

# import numpy as np
# NRL = 197
# n_nucs = 2
# n_bp = 422
#
# dyads = np.asarray(NRL * (np.arange(0, n_nucs, 1) - (n_nucs - 1) / 2.0))
# dyads = (dyads + n_bp // 2).astype(int)
#
# print(dyads)
# nmc.main()
#
import pandas as pd
import matplotlib.pyplot as plt
import FileIO as fileio
import RunMC as runmc
# runmc.main()


# params = {}
# params['n_nuc'] = 5
# params['dummy_steps'] = 1
# params['num_npz'] = 50
# params['tail_switch'] = False
#
# results, filename = runmc.main(3, root=None, input=params)
# runmc.main(10, '2x197x1s500w50-0')
# data = pd.DataFrame(results, index=['test'])
# data.drop(index='test')
#
#
# NRL = [167, 197]
# start = [0,1,2]
# params['dummy_steps'] = 100
#
# for nrl in NRL:
#     for s in start:
#         params['NRL'] = nrl
#         params['fiber_start'] = s
#         results, filename = runmc.main(2e4, root=None, input=params)
#         data.loc[filename] = results
#
# data.to_excel(r'D:\users\vlaar\data\20201109\old_stack.xlsx', index=True, header=True)
