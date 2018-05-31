import numpy as np
import math as math
from matplotlib import pyplot as plt

kT = 4.1


def WLC(f_pN, L_bp, P_nm=50, S_pN=1000):
    z = 1 - 0.5 * np.sqrt(kT / (f_pN * P_nm)) + f_pN / S_pN
    g = f_pN - np.sqrt(kT * f_pN / P_nm) + f_pN ** 2 / 2 * S_pN
    return z * L_bp * 0.34, -g * L_bp * 0.34


def fiber(f_pN, n_nuc, k_pN_nm=0.3, z0_nm=1.7):
    z = f_pN / k_pN_nm + z0_nm
    g = f_pN ** 2 / 2 * k_pN_nm + f_pN * z0_nm
    return z * n_nuc, -g * n_nuc


def tether(f_pN, L_bp, NRL_bp, n_nuc, g1_kT=17, g2_kT=5.5, g3_kT=40, l1_bp=100, l2_bp=80):
    z_nm = []
    d = []
    g_kT = []

    for n0 in range(n_nuc + 1):
        for n1 in range(n_nuc - n0 + 1):
            for n2 in range(n_nuc - n0 - n1 + 1):
                n3 = n_nuc - (n0 + n1 + n2)
                lt_bp = L_bp - n0 * NRL_bp - n1 * (NRL_bp - l1_bp) - n2 * (NRL_bp - l2_bp)
                zt_nm, gt_kT = WLC(f_pN, lt_bp)
                z_nm.append(zt_nm)

                gt_kT += n0 * (g1_kT + g2_kT + g3_kT)
                print(n0, n1, n2, n3, lt_bp, gt_kT)
                gt_kT += n1 * (g2_kT + g3_kT)
                gt_kT += n2 * (g3_kT)
                gt_kT -= f_pN * zt_nm / kT
                g_kT.append(gt_kT)

                dt = math.factorial(n_nuc)
                dt /= (math.factorial(n1) * math.factorial(n2) * math.factorial(n3))
                dt /= math.factorial(n_nuc - (n0 + n1 + n2))
                d.append(np.full((np.shape(f_pN)), dt))

    z_nm = np.asarray(z_nm)
    g_kT = np.asarray(g_kT)

    g_min = np.min(g_kT, axis=1)
    print(g_min)
    g_kT = np.subtract(g_kT.T, g_min)
    g_kT = g_kT.T

    d = np.asarray(d)

    z = np.sum(z_nm * d * np.exp(- g_kT / kT), axis=0)
    z /= np.sum(d * np.exp(- g_kT / kT), axis=0)
    print(np.shape(z))
    return g_kT.T


f = np.arange(0.05, 10, 0.5)
z = tether(f, 1000, 197, 1)
z, g = WLC(f, 1000)
plt.plot(f, g)
plt.show()
