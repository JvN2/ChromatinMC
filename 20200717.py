# import NucleosomeMC as nmc
#
# nmc.main()
#
import RunMC as runmc

runmc.main(20, '2x177x1s25w2-1')
# # #
# import time
# start_time = time.time()
#
#
#
# import Tails as tMC
# import numpy as np
#
# z = np.linspace(0.01, 10, 1e4)
# for i in z:
#     tMC.gFJC(i)
# # #
# # tMC.expected_value()
# #
# z_nm = 4.6
# # f = tMC.fFJC(z)
# # # #
# # print(f)
# L_nm = 6.8
# b_nm = 0.6
# S_pN = 6300
# kT = 4.1
# # test = 49
# f_pN = np.linspace(0.01,4800, 1e6)
# # print('f_pN init: ', f_pN[test])
# # #
# def coth(x):
#     return np.cosh(x) / np.sinh(x)
#
# def Langevin(x):
#     return (coth(x) - 1.0 / x)
# # #
# z = L_nm * (Langevin(b_nm * f_pN / kT) + f_pN / S_pN)
# # print('x: ', b_nm * f_pN / kT)
# print('z(f): ', z[test])
# z_nm = z[test]
# y = z_nm/L_nm
# p = [2.14234, 4.22785, 3, 0.71716, 0.41103, 0.39165]
#
# # f_z_pN = (kT / b_nm) * ((p[0] * y**3 - p[1] * y**2 + p[2] * y) /
# #                         ((1 - y) * (p[3] * y**3 - p[4] * y**2 - p[5] * y + 1))) #+ S_pN * y**3
#
# f_z_pN = y * ( (3 * kT * S_pN) / (3 * kT + b_nm * S_pN))
#
# print('f(z): ', f_z_pN)

# from scipy import optimize
#
# def test_func(x, a, b):
#     return a * x**b
#
# params, params_covariance = optimize.curve_fit(test_func, f_pN, z)
#
# print(params)
# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return idx, array[idx]
#
# print(find_nearest(z, z_nm))
# print(f_pN[4291])
# f = f_pN[4291]
# print(L_nm * (Langevin(b_nm * f/ kT) + f / S_pN))

# print("--- %s seconds ---" % (time.time() - start_time))