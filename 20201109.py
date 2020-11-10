import multiprocessing
from lmfit import Parameters
import time
import RunMC
import warnings
import re

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    params = {}
    params['tail_switch'] = False
    params['dummy_steps'] = 100
    params['num_npz'] = 50

    iterpars = [
       # [8, 197, 2, 100, 28.0],
        [8, 167, 0, 100, 28.0],
        [8, 167, 1, 100, 28.0],
        # [8, 167, 2, 100, 28.0],
        # [8, 197, 0, 100, 28.0],
        # [8, 197, 1, 100, 28.0],
        # [8, 197, 2, 100, 28.0],

       #  [8, 167, 0, 102.5, 5.0],
       #  [8, 167, 0, 102.5, 25.0],
       #  [8, 167, 0, 102.5, 45.0],
       #  [8, 167, 0, 102.5, 65.0],
       #  [8, 167, 0, 102.5, 85.0],
       #  [8, 167, 0, 102.5, 100.0],
       #  [8, 167, 1, 102.5, 5.0],
       #  [8, 167, 1, 102.5, 25.0],
       #  [8, 167, 1, 102.5, 45.0],
       #  [8, 167, 1, 102.5, 65.0],
       #  [8, 167, 1, 102.5, 85.0],
       #  [8, 167, 1, 102.5, 100.0],
       #  [8, 167, 2, 102.5, 5.0],
       #  [8, 167, 2, 102.5, 25.0],
       #  [8, 167, 2, 102.5, 45.0],
       #  [8, 167, 2, 102.5, 65.0],
       #  [8, 167, 2, 102.5, 85.0],
       #  [8, 167, 2, 102.5, 100.0],
    ]

    if len(iterpars) > 22:
        print('>>> More than 22 processes, only the first 22 are executed')
        iterpars = iterpars[:22]

    n_steps = 500

    procs = len(iterpars)  # Number of processes to create

    # Create a list of jobs and then iterate through
    # the number of processes appending each process to
    # the job list
    jobs = []
    for i, iterpar in enumerate(iterpars):

        root = '{0}x{1}x{2}s{3}w{4:0.1f}'.format(iterpar[0], iterpar[1], iterpar[2], iterpar[3], iterpar[4]).replace('.', '-')
        #
        # params['n_nuc'] = iterpar[0]
        # params['NRL'] = iterpar[1]
        # params['fiber_start'] = iterpar[2]
        # params['Rep_Amp_pNA'] = iterpar[3]
        # params['Rep_decay_A'] = iterpar[4]

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
