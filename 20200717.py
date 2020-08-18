# import NucleosomeMC as nmc
#
# nmc.main()
#
# import RunMC as runmc
#
# runmc.main(20000, '2x167x1s25w2-1')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file = r'D:\users\vlaar\data\20200818\2x197x1s25w2-1_001tails.xlsx'
tails = pd.read_excel(file)/10
dist_up_nm = tails['Tail up (nm)']
dist_down_nm = tails['Tail down (nm)']
#
# dist_up_nm = [d[0]/1 for d in tails]
# dist_down_nm = [d[1]/1 for d in tails]
#
plt.rcParams.update({'font.size': 22})

fig, ax = plt.subplots()
# # if orientation is '*-':
# ax.plot(dist_up_nm, color=(1,0,1), marker='o', label='tail up', markersize=12, linestyle='')
# ax.plot(dist_down_nm, color=(0.75,0,0.25), marker='o', label='tail down', markersize=12, linestyle='')
# if orientation is '-*':
ax.plot(dist_up_nm, color=(0.75,0,0.25), marker='o', label='tail up', markersize=2, linestyle='')
ax.plot(dist_down_nm, color=(1,0,1), marker='o', label='tail down', markersize=2, linestyle='')
# default plot parameters

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# spines = ax.spines
# [i.set_linewidth(2) for i in spines.values()]
plt.setp(ax.spines.values(), linewidth=2)
ax.tick_params(which='both', width=2, length=5, top=True, right=True)
# ax.xaxis.set_tick_params(width=5, size=5)
# ax.yaxis.set_tick_params(width=5, size=10)
ax.set_ylim(bottom=0, top=(max(dist_down_nm)+5))

plt.legend(frameon=False, loc=3, markerscale=6)
plt.ylabel('Distance (nm)')
plt.xlabel('Iteration (#)')


# # save plot
# fig.set_size_inches(16, 9)
# fig.savefig(fileio.change_extension(filename, 'tails.png'), dpi=300)
# # # save tails in xlsx
# df = pd.DataFrame(np.array(tails)/1, columns = ['Tail up (nm)', 'Tail down (nm)'])
# df.to_excel(fileio.change_extension(filename, 'tails.xlsx'), index = False, header=True)


plt.show()