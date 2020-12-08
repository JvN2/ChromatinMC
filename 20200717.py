import NucleosomeMC as nmc
import TailsAnalyze as TA
import numpy as np
import matplotlib.pyplot as plt


# TA.plot_npz(r"D:\users\vlaar\data\20201202\8x197x1s25w2-1_002")

# TA.repulsion_exp()
# TA.expo_decay()

# TA.get_stack_params(r"D:\users\vlaar\data\20201202\8x197x2s102w13-1_002")
# TA.dna_energy_display(r"D:\users\vlaar\data\20201202\8x197x2s102w13-1_001.xlsx")

plt.plot(np.arange(0, 10, 1), np.random.random(10))
TA.format_plot('x-as', 'y-ax', 'title', scale_page=1.0, aspect=0.5, save=r"D:\Downloads\test.png")