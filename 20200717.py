import NucleosomeMC as nmc
import TailsAnalyze as TA
import numpy as np
import matplotlib.pyplot as plt
import glob
import FileIO as fileio

# get list of xlsx files in filename folder
xlsx_f = glob.glob(fileio.change_extension(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201110 repulsion d varies\20201110 New stacking 2200iter", '\*.xlsx'))

for f in xlsx_f:
    TA.plot_npz(f)

# TA.repulsion_exp()
# TA.expo_decay()
# TA.tail_energy()
# print(TA.get_stack_params(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 fixed parameters stacking\8x167x0s102w79-1_001.xlsx"))
# TA.dna_energy_display(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 new stacking two decay lengths\8x167x1s102w13-1_001.xlsx")

# TA.stack_exp('twist')
# print('('u'\xb0'')')
# TA.get_g_dna(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 fixed parameters stacking\8x197x1s102w79-1_001")
# TA.de_grote_chromatine_show(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 fixed parameters stacking\8x167x2s102w79-1_001.xlsx", 8)
# TA.plot_g_linker(r"D:\users\vlaar\data\20201216\8x167x1s2500000-0w0-5_001.xlsx",r"D:\users\vlaar\data\20201216\8x167x2s2500000-0w0-5_001.xlsx")

# TA.plot_tail2(r"D:\users\vlaar\data\20201216\8x167x1s2500000-0w0-5_001_tail.xlsx", r"D:\users\vlaar\data\20201216\8x167x2s2500000-0w0-5_001_tail.xlsx")