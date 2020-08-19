# import NucleosomeMC as nmc
#
# nmc.main()
#
# import RunMC as runmc
#
# runmc.main(20000, '2x167x1s25w2-1')


import numpy as np
import matplotlib.pyplot as plt
import os as os
import pandas as pd
from pynverse import inversefunc

kT = 41.0
L=9.88
b=0.44

def coth(x):
    return np.cosh(x) / np.sinh(x)

def Langevin(x):
    return (coth(x) - 1.0 / x)

def fFJC(z_nm, L, b):
    zFJC = lambda f: L * Langevin(f * b / kT)
    f = inversefunc(zFJC, y_values=z_nm, domain=0, open_domain=[True,False], image=[0,L])
    return f
