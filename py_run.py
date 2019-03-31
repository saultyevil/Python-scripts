#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run a batch of Python simulations. Currently, the script will search recursively for Python pfs (disregarding anything
which is py_wind.pf or .out.pf files) and execute a number of commands depending on what is requested by the user.

Currently, the script is able to:
    - Run simulations with -r
    - Resume simulations with -resume
    - Check the convergence of simulations with -c
    - Plot output of simulations with -p
    - Default action -d which is -r -c and -p

The script can also be run in a directory containing only one Python pf.

TODO: be able to run simulations from file input
TODO: flexible flag choices for py_plot.py - input option for the script
"""

import argparse
import datetime
import py_run_util
import py_plot_util
from sys import exit
from shutil import which
from typing import Tuple, Union, List
from subprocess import Popen, PIPE, CalledProcessError

CONVERGED = \
    """
                                                 _
      ___ ___  _ ____   _____ _ __ __ _  ___  __| |
     / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
    | (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
     \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                                  |___/
    """

NOT_CONVERGED = \
    """
                 _                                                 _
     _ __   ___ | |_    ___ ___  _ ____   _____ _ __ __ _  ___  __| |
    | '_ \ / _ \| __|  / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
    | | | | (_) | |_  | (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
    |_| |_|\___/ \__|  \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                                                    |___/
    """

ITS_A_MYSTERY = \
    """
     _ _   _                                   _
    (_) |_( )___    __ _   _ __ ___  _   _ ___| |_ ___ _ __ _   _
    | | __|// __|  / _` | | '_ ` _ \| | | / __| __/ _ \ '__| | | |
    | | |_  \__ \ | (_| | | | | | | | |_| \__ \ ||  __/ |  | |_| |
    |_|\__| |___/  \__,_| |_| |_| |_|\__, |___/\__\___|_|   \__, |
                                     |___/                  |___/
    """

DATE = datetime.datetime.now()
PY_VERSION = "py"
DRY_RUN = False
RUN_SIMS = False
RESUME_RUN = False
N_CORES = 0
CONVERGENCE = False
CLIM = 0.90
CREATE_PLOTS = False
TDE_PLOT = False
WMIN = None
WMAX = None
NOT_QUIET = True
VERBOSE = False


def get_run_mode() -> None:
    """
    Parse the different run modes for the script from the command line. The run options are returned via global
    variables, because why not? :^)

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    global VERBOSE
    global RUN_SIMS
    global RESUME_RUN
    global CONVERGENCE
    global CREATE_PLOTS
    global CLIM
    global TDE_PLOT
    global PY_VERSION
    global WMIN
    global WMAX
    global DRY_RUN
    global NOT_QUIET
    global N_CORES

    p = argparse.ArgumentParser(description="Enable or disable features")
    # Different run modes
    p.add_argument("-d", action="store_true", help="Run with -r -c -p")
    p.add_argument("-r", action="store_true", help="Run simulations")
    p.add_argument("-resume", action="store_true", help="Restart a previous run")
    p.add_argument("-c", action="store_true", help="Check the convergence of runs - calls convergence.py")
    p.add_argument("-p", action="store_true", help="Run plotting scripts")
    p.add_argument("-tde", action="store_true", help="Enable TDE plotting")
    # Script parameters
    p.add_argument("-py_ver", type=str, action="store", help="Name of the Python executable")
    p.add_argument("-clim", type=float, action="store", help="The convergence limit: c_value < 1")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose outputting")
    p.add_argument("-dry", action="store_true", help="Print the simulations found and exit")
    p.add_argument("-q", action="store_true", help="Enable quiet mode")
    p.add_argument("-n_cores", action="store", help="The number of processor cores to run Python with")
    # Plot parameters
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to display")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to display")
    args = p.parse_args()

    do_something = False

    # Set up different run modes
    if args.d:
        RUN_SIMS = True
        CONVERGENCE = True
        CREATE_PLOTS = True
        do_something = True
    if args.verbose:
        VERBOSE = True
    if args.r:
        RUN_SIMS = True
        do_something = True
    if args.resume:
        RESUME_RUN = True
    if args.c:
        CONVERGENCE = True
        do_something = True
    if args.p:
        CREATE_PLOTS = True
        do_something = True
    if args.tde:
        TDE_PLOT = True
    # Set up script parameters
    if args.py_ver:
        PY_VERSION = args.py_ver
    if args.clim:  # Bad variable names here :-) clim is from the arg list
        if 0 < args.c_value < 1:
            CLIM = args.clim
        else:
            print("Invalid value of c_value {}".format(args.c_value))
            exit(1)
    if args.dry:
        DRY_RUN = True
        do_something = True
    if args.q:
        NOT_QUIET = False
    if args.n_cores:
        N_CORES = int(args.n_cores)
    # Set up plotting parameters
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax

    print("--------------------------\n")
    print("Python version ............. {}".format(PY_VERSION))
    print("Show Verbose Output ........ {}".format(VERBOSE))
    print("Run Simulations ............ {}".format(RUN_SIMS))
    print("Resume Run ................. {}".format(RESUME_RUN))
    print("Show Convergence ........... {}".format(CONVERGENCE))
    if CONVERGENCE:
        print("Convergence Limit .......... {}".format(CLIM))
    print("Create Plots ............... {}".format(CREATE_PLOTS))
    if CREATE_PLOTS:
        print("wmin ....................... {}".format(WMIN))
        print("wmax ....................... {}".format(WMAX))
    if TDE_PLOT:
        print("Plot TDE ................... {}".format(TDE_PLOT))
    print("")

    if do_something is False:
        print("No run mode parameter provided, there is nothing to do!\n")
        p.print_help()
        print("\n--------------------------")
        exit(0)

    return


def py_run(root: str, work: str, use_mpi: bool, n_cores: int) -> Tuple[list, list]:
    """
    Function to control running and logging a Python simulation.
    """

    outf_name = "{}/{}_{}{}{}.txt".format(work, root, DATE.year, DATE.month, DATE.day)
    outf = open(outf_name, "w")

    # Construct shell command to run Python and use subprocess to run
    command = ""
    pf = root + ".pf"
    command += "cd {}; Setup_Py_Dir; ".format(work)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(PY_VERSION)
    if RESUME_RUN:
        command += "-r "
    command += "{} ".format(pf)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    print("{}\n".format(command))

    # Real time output of Python
    lines = []
    spec_cycle = False
    for stdout_line in iter(cmd.stdout.readline, ""):
        if not stdout_line:
            break
        line = stdout_line.decode("utf-8").replace("\n", "")
        outf.write("{}\n".format(line))
        lines.append(line)
        spec_cycle = py_run_util.process_line_output(line, spec_cycle, NOT_QUIET, VERBOSE)
    print("")
    outf.close()

    # Have to do this in case the buffer is too large and causes a deadlock
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        print("Captured from stderr:")
        print(err)
        errfname = "{}/err_{}{}{}.out" \
            .format(work, root, DATE.year, DATE.month, DATE.day)
        with open(errfname, "w") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        raise CalledProcessError(rc, command)

    # Append a file containing the Python version and commit hash
    version, hash = py_plot_util.get_python_version(PY_VERSION, VERBOSE)
    with open("version", "w") as f:
        f.write("{}\n{}".format(version, hash))

    return lines, err


def check_python_convergence(root: str, dir: str) -> Union[float, int]:
    """
    Check the convergence of a Python simulation by reading the diag file.
    """

    if 0 > CLIM > 1:
        print("py_run.check_python_convergence: convergence_limit ({}) should be between 0 and 1".
              format(CLIM))
        return -1

    convergence_fraction = py_run_util.check_convergence(root, dir)
    print("convergence limit ............ {}".format(CLIM))
    print("actual convergence ........... {}\n".format(convergence_fraction))

    if convergence_fraction < CLIM:
        print(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(dir, root, convergence_fraction))
    elif convergence_fraction >= CLIM :
        print(CONVERGED)
    else:
        print(ITS_A_MYSTERY)
    print("")

    return convergence_fraction


def plot_python_output(root_name: str, work_dir: str, extra_commands: str = None) -> None:
    """
    Call py_plot.py and plot the output from a Python simulation.
    """

    path = which("py_plot.py")
    if not path:
        print("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(work_dir, root_name)
    if WMIN:
        commands += " -wmin {}".format(WMIN)
    if WMAX:
        commands += " -wmax {}".format(WMAX)
    if extra_commands:
        commands += "{}".format(extra_commands)
    print(commands)

    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        print(output)
    if err:
        print("Captured from stderr:")
        print(err)
    print("")

    return


def plot_tde(root: str, dir: str, verbose: bool = False) -> None:
    """
    Do some tde plotting :-)
    """

    path = which("tde_spec_plot.py")
    if not path:
        print("tde_spec_plot.py not in $PATH and executable")
        return

    command = "cd {}; tde_spec_plot.py {} -tde iPTF15af -l".format(dir, root)
    print("{}\n".format(command))

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8")
    err = stdout.decode("utf-8")

    if verbose:
        print(out)
    if err:
        print("Captured from stderr:")
        print(err)

    return


def run_python_etc(pf_paths: List[str], n_sims: int, use_mpi: bool, n_cores: int) -> None:
    """
    Execute all of the different commands to run a Python simulation
    """

    if CONVERGENCE:
        open("not_converged.txt", "w").close()
    outfname = "py_run_{}{}{}.txt".format(DATE.year, DATE.month, DATE.day)
    f = open(outfname, "w")

    for i, path in enumerate(pf_paths):
        root_name, pf_relative_path = py_plot_util.get_root_name_and_path(path)

        print("--------------------------\n")
        print("Simulation {}/{}".format(i + 1, n_sims))
        print("\n--------------------------\n")
        print("Working directory ......... {}".format(pf_relative_path))
        print("Python root name .......... {}\n".format(root_name))

        f.write("--------------------------\n")
        f.write("Simulation {}/{}\n".format(i + 1, n_sims))
        f.write("--------------------------\n")
        f.write("Working directory ......... {}\n".format(pf_relative_path))
        f.write("Python root name .......... {}\n".format(root_name))

        if RUN_SIMS:
            print("Running the simulation: {}\n".format(root_name))
            python_output, python_err = py_run(root_name, pf_relative_path, use_mpi, n_cores)
            for line in python_output:
                f.write("{}\n".format(line))
            if python_err:
                f.write(python_err)
        if CONVERGENCE:
            print("Checking the convergence of the simulation:\n")
            convergence_fraction = check_python_convergence(root_name, pf_relative_path)
            if convergence_fraction < CLIM:
                f.write("{}\n".format(NOT_CONVERGED))
                f.write("cvalue {} < clim {}\n".format(convergence_fraction, CLIM))
            elif convergence_fraction >= CLIM:
                f.write("{}\n".format(CONVERGED))
                f.write("cvalue {} >= clim {}\n".format(convergence_fraction, CLIM))
            else:
                f.write("{}\n".format(ITS_A_MYSTERY))
        f.write("\n--------------------------\n")

        if CREATE_PLOTS:
            print("Creating plots for the simulation\n")
            plot_python_output(root_name, pf_relative_path, None)
            if TDE_PLOT:
                plot_tde(root_name, pf_relative_path, VERBOSE)

    f.close()

    return


def main() -> None:
    """
    Main control function of the script

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    # Determine which routines to run for each simulation
    get_run_mode()
    pf_paths = py_plot_util.find_pf(ignore_out_pf=True)
    n_sims = len(pf_paths)
    if not n_sims:
        print("No parameter files found, nothing to do!\n")
        print("--------------------------")
        exit(0)

    use_mpi, n_procs = py_run_util.get_num_procs(N_CORES)

    print("")

    if DRY_RUN:
        print("--------------------------")
        return

    # Now run Python, plotting and convergence procedures
    run_python_etc(pf_paths, n_sims, use_mpi, n_procs)

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
