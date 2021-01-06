import NucleosomeMC as nmc
import TailsAnalyze as TA
import numpy as np
import matplotlib.pyplot as plt
import glob
import FileIO as fileio

# get list of xlsx files in filename folder
# xlsx_f = glob.glob(fileio.change_extension(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201110 repulsion d varies\20201110 New stacking 2200iter", '\*.xlsx'))
#
# for f in xlsx_f:
#     TA.plot_npz(f)

# TA.plot_npz(r"D:\Downloads\7_91nm\8x197x2s102w2-5_001")

# TA.repulsion_exp()
# TA.expo_decay()
# TA.debye()
# TA.tail_energy()
# TA.get_stack_params(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201218 Fixed Parameters Stacking\8x197x2s2500000w2-5_001")
# TA.get_stack_params(r"D:\Downloads\20201218 Fixed Parameters Stacking\8x197x0s2500000w2-5_001")
# print(TA.get_stack_params(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 fixed parameters stacking\8x167x0s102w79-1_001.xlsx"))
# TA.dna_energy_display(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201218 Fixed Parameters Stacking\8x197x2s2500000w2-5_001")
# TA.dna_energy_display(r"D:\users\vlaar\data\20201216\8x167x1s2500000-0w0-5_001.xlsx")

# TA.stack_exp('twist')
# print('('u'\xb0'')')
# TA.get_g_dna(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201208 fixed parameters stacking\8x197x1s102w79-1_001")
# TA.de_grote_chromatine_show(r"D:\Downloads\7_91nm\8x197x2s102w2-5_001", 16)
# TA.plot_g_linker(r"D:\Downloads\1_31nm\8x197x1s102w2-5_001",r"D:\Downloads\1_31nm\8x197x2s102w2-5_001")

# TA.plot_tail2(r"D:\users\vlaar\data\20201216\8x167x1s2500000-0w0-5_001_tail.xlsx", r"D:\users\vlaar\data\20201216\8x167x2s2500000-0w0-5_001_tail.xlsx")

# TA.g_dna_kT(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201225 Two salt concentrations exp\7_91nm\8x197x2s102w2-5_001")
# TA.wrap(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201216 Eunwrap\FPS\8x167x1s2500000w5-5_001")


# # get list of files
# files = glob.glob(r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201216 Eunwrap\7_91\8x197x2s102w*_001")
#
# for f in files:
#     # print(f)
#     TA.wrap(f)

TA.plot_wrap(197, 2)