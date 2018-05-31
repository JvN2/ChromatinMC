import numpy as np
import math as math
import scipy.special
from matplotlib import pyplot as plt

kT = 4.1


def WLC(f_pN, L_bp, P_nm=50, S_pN=1000):
    z = 1 - 0.5 * np.sqrt(kT / (f_pN * P_nm)) + f_pN / S_pN
    g = f_pN ** 2 / (2 * S_pN) + 0.5 * np.sqrt(f_pN * kT / P_nm)
    g *= L_bp * 0.34 / kT
    return z * L_bp * 0.34, g


def fiber(f_pN, n_nuc, k_pN_nm=0.3, z0_nm=1.7):
    z = f_pN / k_pN_nm + z0_nm
    g = (f_pN ** 2) / (2 * k_pN_nm * kT)
    return z * n_nuc, g * n_nuc


def tether(f_pN, L_bp, NRL_bp, n_nuc, g1_kT=23, g2_kT=5.5, g3_kT=80, l1_bp=100, l2_bp=80, degeneracy=1):
    z_nm = []
    d = []
    g = []

    for n3 in range(n_nuc + 1):
        for n2 in range(n_nuc - n3 + 1):
            for n1 in range(n_nuc - n3 - n2 + 1):
                n0 = n_nuc - (n3 + n2 + n1)
                lt_bp = L_bp - (n_nuc-n3) * NRL_bp + n1 * (NRL_bp - l1_bp) + n2 * (NRL_bp - l2_bp)
                zt_nm, gt_kT = WLC(f_pN, lt_bp)

                z_fib, g_fib = fiber(f_pN, n0)
                z_nm.append(zt_nm + z_fib)

                gt_kT += -n0 * (g1_kT + g2_kT + g3_kT) + g_fib
                gt_kT += -n1 * (g2_kT + g3_kT)
                gt_kT += -n2 * (g3_kT)
                gt_kT -= f_pN * zt_nm / kT
                g.append(gt_kT)

                dt = scipy.special.binom(n0 + n1, n0)
                dt = (dt - 1) * degeneracy + 1
                dt *= scipy.special.binom(n1 + n2, n2) * scipy.special.binom(n2 + n3, n3)
                d.append(np.full((np.shape(f_pN)), dt))

    z_nm = np.asarray(z_nm)
    g = np.asarray(g)
    d = np.asarray(d)

    g = np.subtract(g, np.min(g, axis=0))

    z = np.sum(z_nm * d * np.exp(- g), axis=0)
    z /= np.sum(d * np.exp(- g), axis=0)

    return z


L_bp = 1000
f = np.arange(0.01, 10, 0.05)

plt.close()
plt.figure(figsize=(4, 3))
plt.axes([0.15, 0.15, .8, .75])
plt.xlim([0, 1.1 * L_bp / 3])
plt.ylim([-0.5, 10.5])

z, g = WLC(f, L_bp)
plt.plot(z, f)

z = tether(f, L_bp, 197, 4)
plt.plot(z, f)

plt.show()
