import multiprocessing
from lmfit import Parameters
import time
import RunMC
import warnings
import re

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    pars = Parameters()
    # Parameters that define the nucleosomal array
    pars.add('L_bp', value=1800)
    pars.add('P_nm', value=50)
    pars.add('n_nuc', value=4)
    pars.add('e_nuc_kT', value=34.7)

    # Parameters that define the folded fiber
    pars.add('rise_A', value=100)
    pars.add('nld_A', value=17)
    pars.add('chirality', value=1)
    pars.add('face', value=1)
    pars.add('diameter_A', value=330)

    pars.add('e_wrap_kT', value=2.0)
    pars.add('e_stack_kT', value=0)
    pars.add('NRL', value=197)
    pars.add('fiber_start', value=1)

    iterpars = [
       # [4, 197, 2, 22, 2.5],
       #  [8, 167, 0, 102, 5.0],
       #  [8, 167, 0, 102, 25.0],
       #  [8, 167, 0, 102, 45.0],
       #  [8, 167, 0, 102, 65.0],
       #  [8, 167, 0, 102, 85.0],
       #  [8, 167, 0, 102, 100.0],
       #  [8, 197, 0, 102, 5.0],
        # [8, 197, 0, 102, 25.0],
        # [8, 197, 0, 102, 45.0],
        # [8, 197, 0, 102, 65.0],
        # [8, 197, 0, 102, 85.0],
        # [8, 197, 0, 102, 100.0],
        # [16, 167, 2, 102, 25.0],
        [16, 197, 1, 102, 25.0],
       #  [8, 167, 1, 102, 45.0],
       #  [8, 167, 1, 102, 65.0],
       #  [8, 167, 1, 102, 85.0],
       #  [8, 167, 1, 102, 100.0],
       #  [8, 167, 2, 102, 5.0],
       #  [8, 167, 2, 102, 25.0],
       #  [8, 167, 2, 102, 45.0],
       #  [8, 167, 2, 102, 65.0],
       #  [8, 167, 2, 102, 85.0],
       #  [8, 167, 2, 102, 100.0],

    ]

    if len(iterpars) > 22:
        print('>>> More than 22 processes, only the first 22 are executed')
        iterpars = iterpars[:22]

    n_steps = 2000

    procs = len(iterpars)  # Number of processes to create

    # Create a list of jobs and then iterate through
    # the number of processes appending each process to
    # the job list
    jobs = []
    for i, iterpar in enumerate(iterpars):
        root = '{0}x{1}x{2}s{3}w{4:0.1f}'.format(iterpar[0], iterpar[1], iterpar[2], iterpar[3], iterpar[4]).replace('.', '-')

        out_list = list()
        process = multiprocessing.Process(target=RunMC.main, args=(n_steps, root))
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
