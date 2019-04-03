#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Various functions used throughout batch running of Python simulations. This script should be imported into other scripts
rather than being itself run.
"""

import time
import datetime
from shutil import which
from platform import system
from typing import Tuple, Union
from subprocess import Popen, PIPE
from multiprocessing import cpu_count


def process_line_output(line: str, spec_cycle: bool, print_crap: bool = True, verbose: bool = False) -> bool:
    """
    Process the output from a Python simulation and print something to screen. Very ugly! Sad!

    Parameters
    ----------
    line                str
                        The line to process
    spec_cycle          bool
                        If True then the line will be processed as a spectral cycle instead
    print_crap          bool, optional
                        If this is False, then all output to screen will be suppressed
    verbose             bool, optional
                        If this is True, then every line will be printed to screen

    Returns
    -------
    spec_cycle          bool
                        Indicates if the previously processes line was a spectral cycle line or not
    """

    if verbose:
        print(line)
    elif line.find("for defining wind") != -1 and print_crap:
        line = line.split()
        cycle = int(line[3])  # + 1
        ncycles = line[5]
        current_time = time.strftime("%H:%M")
        print("{} : Ionisation Cycle ....... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("to calculate a detailed spectrum") != -1 and print_crap:
        line = line.split()
        spec_cycle = True
        cycle = int(line[1])  # + 1
        ncycles = line[3]
        current_time = time.strftime("%H:%M")
        print("{} : Spectrum Cycle ......... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("per cent") != -1 and line.find("Photon") != -1 and print_crap:
        line = line.split()
        print("   - {}% of {} photons transported".format(line[-3], line[-5]))
    elif line.find("!!Check_converging:") != -1 and print_crap:
        line = line.split()
        nconverged = int(line[1])
        fconverged = line[2]
        print("   - {} cells converged {}".format(nconverged, fconverged))
    elif (line.find("Completed ionization cycle") != -1 or line.find("Completed spectrum cycle") != -1) and print_crap:
        line = line.split()
        elapsed_time_seconds = float(line[-1])
        elapsed_time = datetime.timedelta(seconds=elapsed_time_seconds // 1)
        print("   - Elapsed run time {} hrs:mins:secs".format(elapsed_time))
    elif line.find("photon transport completed in") != -1 and print_crap:
        line = line.split()
        transport_time_seconds = float(line[5])
        transport_time = datetime.timedelta(seconds=transport_time_seconds // 1)
        print("   - Photon transported in {} hrs:mins:secs".format(transport_time))
    elif line.find("Completed entire program.") != -1 and print_crap:
        line = line.split()
        tot_run_time_seconds = float(line[-1])
        tot_run_time = datetime.timedelta(seconds=tot_run_time_seconds // 1)
        print("\nSimulation completed in {} hrs:mins:secs".format(tot_run_time))

    return spec_cycle


def find_number_of_physical_cores_lscpu() -> float:
    """
    Find the physical number of CPU cores in a computer. This function looks at the number of available CPUs
    and the number of cores per CPU.

    Note that this will only work with systems where lscpu is installed, i.e. Linux systems.

    Parameters
    ----------
    None

    Returns
    -------
    ncores * nsockets       float
                            The number of physical CPU cores available
    """

    # Find the number of cores per CPU and number of CPUs using lscpu
    ncores_cmd = "lscpu | grep 'Core(s) per socket'"
    nsockets_cmd = "lscpu | grep 'Socket(s):'"
    ncores, stderr = Popen(ncores_cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    nsockets, stderr = Popen(nsockets_cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()

    if ncores == "":
        return 0
    if nsockets == "":
        return 0

    ncores = int(ncores.decode("utf-8").replace("\n", "").split()[-1])
    nsockets = int(nsockets.decode("utf-8").replace("\n", "").split()[-1])

    return ncores * nsockets


def get_num_procs(default_cores: int = 0) -> Tuple[bool, int]:
    """
    Determine the number of Python processes to run. For linux systems, and probably Windows :^), this will use
    lscpu to figure out the number of cores to use. On macOS, this will use multiprocessing.cpu_count() to figure
    out the number of cores to use - note that this will include virtual threads (hyperthreads) UGH.

    Parameters
    ----------
    default_cores       int, optional
                        If this is provided, the function will assume that this is the number of cores Python
                        will be run using and that the computer is capable of using this many cores

    Returns
    -------
    mpi                 bool
                        True will be returned if mpirun is found in the system
    n_cores             int
                        The number of cores to parallelise Python with
    """

    mpi = which("mpirun")
    if mpi:
        mpi = True
    else:
        mpi = False

    # If the user provided the number of cores to use, use that instead :-)
    if default_cores:
        return mpi, default_cores
    else:
        if system() == "Darwin":
            n_cores = cpu_count()
        else:
            n_cores = find_number_of_physical_cores_lscpu()
            if n_cores == 0:
                n_cores = cpu_count()

    return mpi, n_cores


def check_convergence(root_name: str, work_dir: str) -> Union[int, float]:
    """
    Check the convergence of a Python simulation.

    Parameters
    ----------
    root_name           str
                        The root name of the Python simulation
    work_dir            str
                        The directory containing the Python simulation

    Returns
    -------
    converge_fraction   int or float
                        If an int is return, there has been an error in finding the convergence fraction of
                        the simulation. Otherwise, the fraction of wind cells which converged is returned.
    """

    diag_path = "{}/diag_{}/{}_0.diag".format(work_dir, root_name, root_name)

    try:
        with open(diag_path, "r") as file:
            diag = file.readlines()
    except IOError:
        print(
            "py_util.read_convergence: Couldn't open read only copy of {}. Does the diag file exist?".format(diag_path))
        return -1

    converge_lines = []
    converge_fraction = None
    for line in diag:
        if line.find("!!Check_converging") != -1:
            c_string = line.split()[2].replace("(", "").replace(")", "")
            converge_fraction = float(c_string)
            converge_lines.append(line)

    if converge_fraction is None:
        print("py_util.read_convergence: unable to parse convergence fraction from diag file {}"
              .format(diag_path))
        return -1

    if 0 > converge_fraction > 1:
        print("py_util.read_convergence: the convergence in the simulation is negative or more than one")
        print("py_util.read_convergence: convergence_fraction = {}".format(converge_fraction))
        return -1

    return converge_fraction
