import NucleosomeMC as nmc
import TailsAnalyze as TA

TA.repulsion_exp()
#
# nmc.main()
#
import pandas as pd
import matplotlib.pyplot as plt
import FileIO as fileio
import RunMC as runmc


# params = {}
# params['n_nuc'] = 5
# params['dummy_steps'] = 1
# params['num_npz'] = 50
# params['tail_switch'] = False
#
# results, filename = runmc.main(3, root=None, input=params)
# runmc.main(20, '2x197x1s500w50-0')
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

# # # # # # #
#
import Tails as tMC

#
# Amp = 100 # amplitude pNA
# decay_l = 28.0 # decay length A
# #
# x = np.linspace(0,100)
# y = Amp * np.exp(- (1 / decay_l) * x)
#
# # for Amp in np.linspace(1200, 12000, 10):
# #     print('amp: ', Amp)
# #
# #     twin = Amp * np.exp(- (1 / decay_l) * 20)
# #     vijf = Amp * np.exp(- (1 / decay_l) * 35)
# #
# #     print('20: ', twin / 41.0)
# #     print('35: ', vijf / 41.0)
# #     print('verschil kT: ', (twin - vijf) / 41.0)
#
# tMC.plotten(x,y, xlabel='distance (A)', ylabel='Energy (pN$\AA$)')
#
# from helixmc import util
# import glob
# import FileIO as fileio
# import pandas as pd
# #
