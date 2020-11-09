import multiprocessing
from lmfit import Parameters
import time
import RunMC
import warnings
import re

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    params = {}
    params['tail_switch'] = True
    params['dummy_steps'] = 3
    params['num_npz'] = 5

    iterpars = [
       [8, 197, 2, 100, 28.0],
       #  [8, 167, 2, 0, 2.5],
       #  [8, 167, 2, 0, 2.5],
       #  [8, 167, 2, 0, 2.5],
       #  [8, 168, 2, 25, 2.5],
       #  [8, 168, 2, 25, 2.5],
       #  [8, 168, 2, 25, 2.5],
       #  [8, 168, 1, 25, 2.5],
       #  [8, 168, 1, 25, 2.5],
       #  [8, 168, 1, 25, 2.5],
       #  [8, 169, 2, 25, 2.5],
       #  [8, 169, 2, 25, 2.5],
       #  [8, 169, 2, 25, 2.5],
       #  [8, 169, 1, 25, 2.5],
       #  [8, 169, 1, 25, 2.5],
       #  [8, 169, 1, 25, 2.5],
       #  [8, 170, 2, 25, 2.5],
       #  [8, 170, 2, 25, 2.5],
       #  [8, 170, 2, 25, 2.5],
       #  [8, 170, 1, 25, 2.5],
       #  [8, 170, 1, 25, 2.5],
       #  [8, 170, 1, 25, 2.5],
    ]

    if len(iterpars) > 22:
        print('>>> More than 22 processes, only the first 22 are executed')
        iterpars = iterpars[:22]

    n_steps = 28

    procs = len(iterpars)  # Number of processes to create

    # Create a list of jobs and then iterate through
    # the number of processes appending each process to
    # the job list
    jobs = []
    for i, iterpar in enumerate(iterpars):

        root = None

        params['n_nuc'] = iterpar[0]
        params['NRL'] = iterpar[1]
        params['fiber_start'] = iterpar[2]
        params['Rep_Amp_pNA'] = iterpar[3]
        params['Rep_decay_A'] = iterpar[4]

        out_list = list()
        process = multiprocessing.Process(target=RunMC.main, args=(n_steps, root, params))
        jobs.append(process)

    # Start the processes
    for j in jobs:
        print('\n>>> Start job {}'.format(j))
        j.start()
        time.sleep(5)

        # Ensure all of the processes have finished
    for j in jobs:
        j.join()

    print "List processing complete."
