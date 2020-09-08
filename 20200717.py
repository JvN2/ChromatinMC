# import NucleosomeMC as nmc
#
# nmc.main()
#
import RunMC as runmc

runmc.main(28, '3x197x1s25w2-1')

# import numpy as np
# import matplotlib.pyplot as plt
# from pynverse import inversefunc
# import matplotlib as mpl
#
# import Tails as Tmc
#
# z = np.linspace(0.1,20,200)
# f = Tmc.fFJC(z_nm=z)
# g = Tmc.gFJC(f)
#
# plt.rcParams.update({'font.size': 22})
#
# fig, ax = plt.subplots()
# # # if orientation is '*-':
# # ax.plot(dist_up_nm, color=(1,0,1), marker='o', label='tail up', markersize=12, linestyle='')
# # ax.plot(dist_down_nm, color=(0.75,0,0.25), marker='o', label='tail down', markersize=12, linestyle='')
# # if orientation is '-*':
# ax.plot(f,g, color=(0.75, 0, 0.25), markersize=6, linestyle='-', linewidth=6)
#
# # default plot parameters
#
# # ax.spines['top'].set_visible(False)
# # ax.spines['right'].set_visible(False)
# # spines = ax.spines
# # [i.set_linewidth(2) for i in spines.values()]
# plt.setp(ax.spines.values(), linewidth=2)
# ax.tick_params(which='both', width=2, length=5, top=True, right=True)
# # ax.xaxis.set_tick_params(width=5, size=5)
# # ax.yaxis.set_tick_params(width=5, size=10)
# # ax.set_ylim(bottom=0, top=(max(dist_down_nm) + 5))
#
# # plt.legend(frameon=False)
# plt.ylabel('Energy (pN*nm/kT)')
# plt.xlabel('Force (pN)')
# # plt.ylabel('Force (pN)')
# # plt.xlabel('Distance (nm)')
#
# plt.show()
#
