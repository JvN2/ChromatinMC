# import NucleosomeMC as nmc
#
# nmc.main()
# #
import pandas as pd
import FileIO as fileio
import RunMC as runmc


# params = {}
# params['n_nuc'] = 5
# params['dummy_steps'] = 1
# params['num_npz'] = 50
# params['tail_switch'] = False
#
# results, filename = runmc.main(3, root=None, input=params)
runmc.main(3, '4x167x0s100w28-0')
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
# tMC.sequence()

# # tMC.expected_value()
# # z = np.linspace(0,20,1e4)
# z = 3.28
# # import time
# # start_time = time.time()
# # for z in z:
# tMC.gFJC(z)
#
# # print("--- %s seconds ---" % (time.time() - start_time))
#
# from helixmc import util
# import glob
# import FileIO as fileio
# import pandas as pd
#
# filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201016 repulsion exp\6x167x2s25w2-1_001"
#
# # get list of npz files in filename folder
# npz_f = glob.glob(fileio.change_extension(filename, '\*.npz'))
#
# params = []
# # params contains all 6 parameters for each basepair out of every npz file
# for f, file in enumerate(npz_f):
#     data = np.load(file)
#     params.append(data['params'])
#
# # get mean value of every parameter for each basepair
# params_m = np.mean(params, axis=0)
#
# # use 6 parameters to get coordinates of every basepair
# dr, frames = util.params2data(params_m)
# coords = util.dr2coords(dr)
#
# df = pd.DataFrame(np.array(coords) / 10)
# df.to_excel(fileio.change_extension(filename, 'coord_m.xlsx'), index=False, header=True)