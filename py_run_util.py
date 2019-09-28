#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Various functions used throughout batch running of Python simulations. This
script should be imported into other scripts rather than being itself run.
"""


import time
import datetime
from platform import system
from typing import Tuple, Union, List
from shutil import which, copyfile
from subprocess import Popen, PIPE
from multiprocessing import cpu_count

LOGFILE = None


def init_logfile(logfile_name: str, use_global_log: bool = True):
    """
    Initialise a logfile global variable

    Parameters
    ----------
    logfile_name: str
        The name of the logfile to initialise
    use_global_log: bool, optional
        If this is false, a object for a logfile will be returned instead.
    """

    global LOGFILE

    if use_global_log:
        if LOGFILE:
            print("{}:{}: logfile already initialised as {}".format(__file__, init_logfile.__name__, LOGFILE.name))
            return
        LOGFILE = open(logfile_name, "a")
    else:
        logfile = open(logfile_name, "a")
        return logfile

    return


def close_logfile(logfile=None) -> None:
    """
    Close a log file for writing - this will either use the log file provided
    or will attempt to close the global log file.

    Parameters
    ----------
    logfile: io.TextIO, optional
        An externa log file object
    """

    global LOGFILE

    if logfile:
        logfile.close()
    elif LOGFILE:
        LOGFILE.close()
    else:
        print("{}:{}: No logfile to close? ahhh".format(__file__, close_logfile.__name__))

    return


def log(message: str, logfile=None) -> None:
    """
    Log a message to screen and to the log file provided or the global log file.

    Parameters
    ----------
    message: str
        The message to log to screen and file
    logfile: io.TextIO, optional
        An open file object which is the logfile to log to. If this is not
        provided, then the global logfile
    """

    print(message)

    if logfile:
        logfile.write("{}\n".format(message))
    elif LOGFILE:
        LOGFILE.write("{}\n".format(message))
    else:
        return

    return


def process_line_output(line: str, pcycle: bool, n_cores: int = 1, print_crap: bool = True, verbose: bool = False) \
        -> bool:
    """
    Process the output from a Python simulation and print something to screen.
    Very ugly! Very sad!

    Parameters
    ----------
    line: str
        The line to process
    pcycle: bool
        If True then the line will be processed as a spectral cycle instead
    n_cores: int, optional
        The number of cores the simulation is being run with. This is required
        to calculate the total photon number
    print_crap: bool, optional
        If this is False, then all output to screen will be suppressed
    verbose: bool, optional
        If this is True, then every line will be printed to screen

    Returns
    -------
    pcycle: bool
        Indicates if the previously processed line was a spectral cycle line
        or not
    """

    if verbose:
        print(line)
    elif line.find("for defining wind") != -1 and print_crap:
        line = line.split()
        cycle = int(line[3])  # + 1
        ncycles = line[5]
        current_time = time.strftime("%H:%M")
        log("{} : Ionisation Cycle ....... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("to calculate a detailed spectrum") != -1 and print_crap:
        line = line.split()
        pcycle = True
        cycle = int(line[1])  # + 1
        ncycles = line[3]
        current_time = time.strftime("%H:%M")
        log("{} : Spectrum Cycle ......... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("per cent") != -1 and line.find("Photon") != -1 and print_crap:
        line = line.split()
        log("      : {}% of {:1.2e} photons transported".format(line[-3], int(int(line[-5]) * n_cores)))
    elif line.find("!!Check_converging:") != -1 and print_crap:
        line = line.split()
        nconverged = int(line[1])
        fconverged = line[2]
        log("      : {} cells converged {}".format(nconverged, fconverged))
    elif (line.find("Completed ionization cycle") != -1 or line.find("Completed spectrum cycle") != -1) and print_crap:
        line = line.split()
        elapsed_time_seconds = float(line[-1])
        elapsed_time = datetime.timedelta(seconds=elapsed_time_seconds // 1)
        log("      : Elapsed run time {} hrs:mins:secs".format(elapsed_time))
    elif line.find("photon transport completed in") != -1 and print_crap:
        line = line.split()
        transport_time_seconds = float(line[5])
        transport_time = datetime.timedelta(seconds=transport_time_seconds // 1)
        log("      : Photons transported in {} hrs:mins:secs".format(transport_time))
    elif line.find("Completed entire program.") != -1 and print_crap:
        line = line.split()
        tot_run_time_seconds = float(line[-1])
        tot_run_time = datetime.timedelta(seconds=tot_run_time_seconds // 1)
        log("\nSimulation completed in {} hrs:mins:secs".format(tot_run_time))

    return pcycle


def ncores_lscpu() -> int:
    """
    Find the physical number of CPU cores in a computer. This function looks at
    the number of available CPUs and the number of cores per CPU.

    Note that this will only work with systems where lscpu is installed,
    i.e. Linux systems.

    Returns
    -------
    ncores * nsockets: int
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

    return int(ncores * nsockets)


def ncores(default_cores: int = 0) -> Tuple[bool, int]:
    """
    Determine the number of Python processes to run. For linux systems, and
    probably Windows :^), this will use lscpu to figure out the number of cores
    to use. On macOS, this will use multiprocessing.cpu_count() to figure
    out the number of cores to use - note that this will include virtual threads
    (hyperthreads) UGH.

    Parameters
    ----------
    default_cores: int, optional
        If this is provided, the function will assume that this is the number
        of cores Python will be run using and that the computer is capable of
        using this many cores

    Returns
    -------
    mpi: bool
        True will be returned if mpirun is found in the system
    n_cores: int
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
            n_cores = cpu_count() // 2  # divide by 2 because hyperthreading :^)
        else:
            n_cores = ncores_lscpu()
            if n_cores == 0:
                n_cores = cpu_count()

    return mpi, n_cores


def check_convergence(root: str, wd: str, verbose: bool = True) -> Union[int, float]:
    """
    Check the convergence of a Python simulation.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    wd: str
        The directory containing the Python simulation
    verbose: bool, optional
        Enable verbose logging if set to True

    Returns
    -------
    converge_fraction: int or float
        If an int is return, there has been an error in finding the convergence
        fraction of the simulation. Otherwise, the fraction of wind cells which
        converged is returned.
    """

    diag_path = "{}/diag_{}/{}_0.diag".format(wd, root, root)

    try:
        with open(diag_path, "r") as file:
            diag = file.readlines()
    except IOError:
        if verbose:
            print("{}:{}: could not open diag file {}".format(__file__, check_convergence.__name__, diag_path))
        return -1

    converge_lines = []
    converge_fraction = None
    for line in diag:
        if line.find("!!Check_convergence") != -1:
            if len(line.split()) != 11:
                continue
            c_string = line.split()[2].replace("(", "").replace(")", "")
            converge_fraction = float(c_string)
            converge_lines.append(line)

    if converge_fraction is None:
        if verbose:
            print("{}:{}: unable to parse convergence from diag {}"
                  .format(__file__, check_convergence.__name__, diag_path))
        return -1

    if 0 > converge_fraction > 1:
        if verbose:
            print("{}:{}: the convergence in the simulation is not sane with convergence {}"
                  .format(__file__, check_convergence.__name__, converge_fraction * 100.0))
        return -1

    return converge_fraction


def change_parameter(pf: str, parameter: str, value: str, verbose: bool = False):
    """
    Search a parameter file for a given parameter and replaces the current value
    with a new value. This script will change the parameter file, even if the
    old and new parameter values are the same :-).

    Parameters
    ----------
    pf: str
        The name of the parameter file to edit
    parameter: str
        The name of the parameter to be edited
    value: str
        The new value of the parameter
    verbose : bool, optional
        Enable verbose logging

    Returns
    -------
    Returns a non-zero integer on non-successful exit.
    """

    assert (type(pf) == str)
    assert (type(parameter) == str)
    assert (type(value) == str)

    if pf.find(".pf") == -1:
        pf += ".pf"

    old = ""
    new = ""

    # Create back up file, in case things go to shit
    copyfile(pf, pf + ".bak")

    # Open file, search for parameter and replace
    with open(pf, "r") as f:
        lines = f.readlines()
    nlines = len(lines)
    for l in range(nlines):
        if lines[l].find(parameter) != -1:
            old = lines[l]
            new = "{} {}\n".format(parameter, value)
            lines[l] = new
            break
    if old and new:
        if verbose:
            print("Changed parameter {} to value {}".format(parameter, value))
            print("OLD: {}".format(old.replace("\n", "")))
            print("NEW: {}".format(new.replace("\n", "")))
    else:
        print("Could not find parameter {} in {}".format(parameter, pf))
        exit(1)

    # Now write out modified lines to file
    with open(pf, "w") as f:
        f.writelines(lines)

    return


def print_error_summary(root: str, wd: str, logfile=None) -> List[str]:
    """
    Print the error summary from a completed or crashed Python run.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str
        The working directory to run the Python simulation in.
    logfile: io.TextIO, optional
        An open file object which is the logfile to log to. If this is None,
        the global LOGFILE in this file will be used instead

    Returns
    -------
    error_summary: List[str]
        A list of strings where each string is an error summary output from
        Python
    """

    diag_file = "{}/diag_{}/{}_0.diag".format(wd, root, root)
    with open(diag_file, "r") as f:
        diag = f.readlines()

    py_error_idx = 0
    error_start_idx = 0
    diaglen = len(diag)

    for i in range(diaglen):
        line = diag[-(i + 1)]
        if line == "Run py_error.py for full error report.\n":
            py_error_idx = diaglen - i
        if line == "Error summary: End of program, Thread 0 only\n":
            error_start_idx = diaglen - i - 1
            break

    error_summary = diag[error_start_idx:py_error_idx]
    n_unique_errors = len(error_summary) - 3
    if n_unique_errors < 0:
        n_unique_errors = 0
    n_total_errors = 0

    for i in range(n_unique_errors):
        idx = i + 2
        error_summary[idx] = error_summary[idx].replace("\n", "").lstrip()
        line = error_summary[idx]
        line_len = len(line)
        j = 0
        for j in range(line_len):
            if line[j].isdigit() is False:
                break
        if j != line_len - 1:
            n = int(line[0:j])
            n_total_errors += n

    log("Process 0:")
    log("Unique error messages ........ {}".format(n_unique_errors), logfile)
    log("Total error messages ......... {}".format(n_total_errors), logfile)
    log("\nError summary:\n")
    for i in range(n_unique_errors):
        idx = i + 2
        log("\t{}".format(error_summary[idx]), logfile)
    log("")

    return error_summary
