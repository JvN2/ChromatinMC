import pandas as pd
import matplotlib.pyplot as plt

def repulsion_exp():
    filename = r"C:\Users\Annelies\OneDrive\Documents\experimental data\20201110\20201112_cms_dists.xlsx"
    zero_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,B,C")
    one_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,D,E")
    two_start_167 = pd.read_excel(filename, sheet_name='167', header=1, index_col=0, usecols="A,F,G")

    x_1 = range(1, 8)
    x_tick = one_start_167.index

    y_0 = zero_start_167.iloc[:, 0]
    e_0 = zero_start_167.iloc[:, 1]

    y_1 = one_start_167.iloc[:, 0]
    e_1 = one_start_167.iloc[:, 1]

    y_2 = two_start_167.iloc[:, 0]
    e_2 = two_start_167.iloc[:, 1]

    #
    #
    plt.rcParams.update({'font.size': 22})
    #
    fig, ax = plt.subplots()
    #
    ax.errorbar(x_1, y_1, e_1, color=(0.25, 0, 0.75), marker='o', markersize=10, label='NRL 167 1-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y_2, e_2, color=(0.5, 0, 0.5), marker='o', markersize=10, label='NRL 167 2-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)
    ax.errorbar(x_1, y_0, e_0, color=(0.75, 0, 0.25), marker='o', markersize=10, label='NRL 167 0-start', linewidth=0,
                ecolor='r', elinewidth=5, capsize=5)

    #
    plt.setp(ax.spines.values(), linewidth=2)
    ax.tick_params(which='both', width=2, length=5, top=True, right=True)
    plt.xticks(x_1, x_tick)
    ax.set_xlim(left=0)
    # ax.set_ylim(bottom=0)
    plt.legend()
    #
    plt.title('Repulsion amplitude 2,5 kT', y=1.08)
    plt.xlabel('decay length (nm)')
    plt.ylabel('Distance between nucleosomes (nm)')
    plt.show()

    return