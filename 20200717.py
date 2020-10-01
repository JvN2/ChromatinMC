# import NucleosomeMC as nmc
#
# nmc.main()
#
import RunMC as runmc

runmc.main(20, '2x167x1s25w2-1')
#
# import Tails as tMC
# import numpy as np
# #
# tMC.expected_value()
# #
# z = 10
# f = tMC.fFJC(z)
# #
# print(f)
# L_nm = 6.8
# b_nm = 0.22
# S_pN = 6300
# kT = 4.1
# f_pN = f
#
# def coth(x):
#     return np.cosh(x) / np.sinh(x)
#
#
# def Langevin(x):
#     return (coth(x) - 1.0 / x)
# # #
# z = L_nm * (Langevin(b_nm * f_pN / kT) + f_pN / S_pN)
# print(z)
#
# # p = np.polyfit(x, y, 3)
# # print(p)