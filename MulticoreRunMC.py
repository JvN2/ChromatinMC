import multiprocessing
from lmfit import Parameters
import time
import RunMC
import warnings

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    pars = Parameters()
    # Parameters that define the nucleosomal array
    pars.add('L_bp', value=1000)
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
        [4, 167, 1, 0, 2.5],
        [4, 167, 1, 25, 2.5],
        [4, 167, 1, 25, 2.5],
        [4, 167, 1, 25, 2.5],
        [4, 167, 2, 25, 2.5],
        [4, 167, 2, 25, 2.5],
        [4, 167, 2, 25, 2.5],
        [4, 197, 1, 0, 2.5],
        [4, 197, 1, 25, 2.5],
        [4, 197, 1, 25, 2.5],
        [4, 197, 1, 25, 2.5],
        [4, 197, 2, 25, 2.5],
        [4, 197, 2, 25, 2.5],
        [4, 197, 2, 25, 2.5],
    ]

    n_steps = 5e4

    procs = len(iterpars)  # Number of processes to create

    # Create a list of jobs and then iterate through
    # the number of processes appending each process to
    # the job list
    jobs = []
    for i, iterpar in enumerate(iterpars):
        pars['n_nuc'].value = iterpar[0]
        pars['NRL'].value = iterpar[1]
        pars['fiber_start'].value = iterpar[2]
        pars['e_stack_kT'].value = iterpar[3]
        pars['e_wrap_kT'].value = iterpar[4]

        root = '{1}x{2}x{0}s{3}w{4:0.1f}'.format(pars['fiber_start'].value, pars['n_nuc'].value, pars['NRL'].value,
                                                 pars['e_stack_kT'].value, pars['e_wrap_kT'].value).replace('.', '-')

        out_list = list()
        process = multiprocessing.Process(target=RunMC.main, args=(pars, n_steps, root))
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
