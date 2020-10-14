# import NucleosomeMC as nmc
#
# nmc.main()
#
import RunMC as runmc

runmc.main(20, '2x167x1s25w2-1')
# # # # # #
#
import Tails as tMC
import numpy as np

# Amp = 41 # amplitude
# decay_l = 0.14 # decay length
#
# x = np.linspace(0,50)
# y = Amp * np.exp(- decay_l * x)
# print(y)

# tMC.plotten(x,y, xlabel='distance (A)', ylabel='Energy (pNnm)')



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