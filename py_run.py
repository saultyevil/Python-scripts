#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run a batch of Python simulations. Currently, the script will search recursively
for Python pfs (disregarding anything which is py_wind.pf or .out.pf files) and
execute a number of commands depending on what is requested by the user.

The script can also be run in a directory containing only one Python pf.

usage: py_run.py [-h] [-d] [-s] [-sc] [-r] [-c] [-p] [-tde] [-path PATH]
                 [-py_ver PY_VER] [-f PY_FLAGS] [-clim CLIM] [-v] [--dry] [-q]
                 [-n_cores N_CORES] [-polar]

General script to run Python simulations

optional arguments:
  -h, --help            show this help message and exit
  -d                    Run with -s -c -p
  -s                    Run simulations
  -sc                   Run spectral cycles even if a simulation hasn't
                        converged
  -r                    Restart a previous run
  -c                    Check the convergence of runs - calls convergence.py
  -p                    Run plotting scripts
  -tde                  Enable TDE plotting
  -path PATH            Provide a list of directories of Python parameter
                        files to run
  -py_ver PY_VER        Name of the Python executable
  -f PY_FLAGS, --py_flags PY_FLAGS
                        Runtime flags to pass to Python
  -clim CLIM            The convergence limit: c_value < 1
  -v, --verbose         Verbose outputting
  --dry                 Print the simulations found and exit
  -q                    Enable quiet mode
  -n_cores N_CORES      The number of processor cores to run Python with
  -polar                Plot using a polar projection
"""

import argparse
from os import access, R_OK
import datetime
import py_rm_data as prd
import py_run_util as rutil
import py_plot_util as putil
from sys import exit
from shutil import which, copyfile
from typing import Union, List
from subprocess import Popen, PIPE
import py_change_parameter as pcp

CONVERGED = \
    r"""
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                                             _
  ___ ___  _ ____   _____ _ __ __ _  ___  __| |
 / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
| (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
 \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                              |___/
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                        \
                         \     ____________
                          \    |__________|
                              /           /\
                             /           /  \
                            /___________/___/|
                            |          |     |
                            |  ==   == |     |
                            |   O   O  | \ \ |
                            |     <    |  \ \|
                           /|          |   \ \
                          / |  \_____/ |   / /
                         / /|          |  / /|
                        /||\|          | /||\/
                            -------------|
                                | |    | |
                               <__/    \__>

"""

NOT_CONVERGED = \
    r"""
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
             _                                                 _
 _ __   ___ | |_    ___ ___  _ ____   _____ _ __ __ _  ___  __| |
| '_ \ / _ \| __|  / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \/ _` |
| | | | (_) | |_  | (_| (_) | | | \ V /  __/ | | (_| |  __/ (_| |
|_| |_|\___/ \__|  \___\___/|_| |_|\_/ \___|_|  \__, |\___|\__,_|
                                                |___/
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
                                            \                                 
                                             \     ____________
                                              \    |__________|
                                                  /           /\
                                                 /           /  \
                                                /___________/___/|
                                                |          |     |
                                                |  ==\ /== |     |
                                                |   O   O  | \ \ |
                                                |   () ()  |  \ \|
                                               /|          |   \ \
                                              / |   _____  |   / /
                                             / /|  /     \ |  / /|
                                            /||\|          | /||\/
                                                -------------|
                                                    | |    | |
                                                   <__/    \__>
"""

ITS_A_MYSTERY = \
    r"""
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _
                                                          _
  ___ ___  _ ____   _____ _ __ __ _  ___ _ __   ___ ___  (_)___
 / __/ _ \| '_ \ \ / / _ \ '__/ _` |/ _ \ '_ \ / __/ _ \ | / __|
| (_| (_) | | | \ V /  __/ | | (_| |  __/ | | | (_|  __/ | \__ \
 \___\___/|_| |_|\_/ \___|_|  \__, |\___|_| |_|\___\___| |_|___/
                              |___/
                                    _
         __ _   _ __ ___  _   _ ___| |_ ___ _ __ _   _
        / _` | | '_ ` _ \| | | / __| __/ _ \ '__| | | |
       | (_| | | | | | | | |_| \__ \ ||  __/ |  | |_| |
        \__,_| |_| |_| |_|\__, |___/\__\___|_|   \__, |
                          |___/                  |___/
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ __ _ _ _ _ _ _ _ _
"""

DATE = datetime.datetime.now()
PY_VERSION = "py"
PY_FLAGS = None
DRY_RUN = False
RUN_SIMS = False
RESUME_RUN = False
N_CORES = 0
CHECK_CONVERGENCE = False
CONV_LIMIT = 0.80
CREATE_PLOTS = False
TDE_PLOT = False
NOT_QUIET = True
VERBOSE = False
SPLIT_CYCLES = True
SIMS_FROM_FILE = False
POLAR = False


def plot_model(root: str, wd: str) -> None:
    """
    Use py_plot.py to create plots of the synthetic spectra and wind variables
    for a Python model.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str
        The working directory for the Python simulation.
    """

    path = which("py_plot.py")
    if path == "":
        rutil.log("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(wd, root)
    if POLAR:
        commands += " -p"

    rutil.log(commands)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        rutil.log("\n{}".format(output))
    if err:
        rutil.log("\nCaptured from stderr:")
        rutil.log(err)

    rutil.log("")

    return


def plot_spec_tde(root: str, wd: str) -> None:
    """
    Use tde_spec_plot.py to create a comparison plot to the TDE iPTF15af.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str
        The working directory for the Python simulation.
    """

    path = which("tde_spec_plot.py")
    if path == "":
        rutil.log("tde_spec_plot.py not in $PATH and executable")
        return

    command = "cd {}; tde_spec_plot.py {}".format(wd, root)

    rutil.log(command)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        rutil.log("\n{}".format(out))
    if err:
        rutil.log("Captured from stderr:")
        rutil.log(err)

    rutil.log("")

    return


def check_python_convergence(root: str, wd: str) -> Union[float, int]:
    """
    Check the convergence of a Python simulation by parsing the master diag file.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    wd: str
        The working directory containing the Python simulation

    Returns
    -------
    rc: bool
        If the simulation has converged, True is returned.
    """

    rc = False
    c_fraction = rutil.check_run_convergence(root, wd)

    rutil.log("convergence limit ............ {}".format(CONV_LIMIT))
    rutil.log("actual convergence ........... {}\n".format(c_fraction))

    if 0 > c_fraction > 1:
        rutil.log(ITS_A_MYSTERY)
    elif c_fraction < CONV_LIMIT:
        rutil.log(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root, c_fraction))
    elif c_fraction >= CONV_LIMIT:
        rutil.log(CONVERGED)
        rc = True

    rutil.log("")

    return rc


def python(root: str, wd: str, use_mpi: bool, n_cores: int, restart_run: bool = False, split_cycles: bool = False,
           restart_from_spec: bool = False) -> int:
    """
    The purpose of this function is to use the Subprocess library to call
    Python. Unfortunately, to cover a wide range of situations with how one
    may want to run Python, this function has become rather complicated and
    could benefit from being modularised further.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str
        The working directory to run the Python simulation in.
    use_mpi: bool
        If True, Python will be run using mpirun.
    n_cores: int
        If use_mpi is True, then Python will be run using the number of cores
        provided.
    restart_run: bool, optional
        If True, the -r flag will be passed to Python to restart a run from the
        previous cycle
    split_cycles: bool, optional
        If True, the -r flag will be passed to Python to restart a run from the
        first spectrum cycle with a reduced photon sample
    restart_from_spec: bool, optional
        If True, Python will probably run just the spectral cycles with a reduced
        photon number.

    Returns
    -------
    rc: int
        The return code from the Python simulation
    """

    logfile_name = "{}/{}_{}{:02d}{:02d}.txt".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
    logfile = open(logfile_name, "a")
    pf = root + ".pf"

    if split_cycles and restart_from_spec is False:
        pcp.change_python_parameter(pf, "Spectrum_cycles", "0", VERBOSE, bakup=True)
    if split_cycles and restart_from_spec:
        pcp.change_python_parameter(pf, "Spectrum_cycles", "5", VERBOSE, bakup=False)
        pcp.change_python_parameter(pf, "Photons_per_cycle", "1e6", VERBOSE, bakup=False)

    # Construct shell command to run Python and use subprocess to run
    command = "cd {}; Setup_Py_Dir; ".format(wd)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(PY_VERSION)

    if restart_run:
        command += "-r "

    if PY_FLAGS:
        if type(PY_FLAGS) != str:
            rutil.log("The provided additional flags for Python is not a string")
            exit(1)
        command += " {} ".format(PY_FLAGS)

    # Add the root name at the end of the call to Python
    command += "{}".format(pf)
    rutil.log("{}\n".format(command))

    # Use Popen to create a new Python process
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    # This next bit provides real time output of Python's output...
    if split_cycles:
        pcycle = True
    else:
        pcycle = False
    for stdout_line in iter(cmd.stdout.readline, ""):
        if not stdout_line:
            break
        line = stdout_line.decode("utf-8").replace("\n", "")
        logfile.write("{}\n".format(line))
        pcycle = rutil.process_line_output(line, pcycle, n_cores, NOT_QUIET, VERBOSE)

    rutil.log("")
    logfile.close()

    # Sometimes with Subprocess, if the output buffer is too large then subprocess
    # breaks and causes a deadlock. To get around this, one can use .communicate()
    # to flush the buffer or s/t
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        rutil.log("Captured from stderr:")
        rutil.log(err)
        errfname = "{}/err_{}{:02d}{:02d}.out".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
        with open(errfname, "a") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        print("Python exited with non-zero exit code: {}\n".format(rc))
        rutil.print_error_summary(root, wd)
        # raise CalledProcessError(rc, command)

    # Create a file containing the Python version and commit hash - helpful when
    # trying to figure out what version you need for py_wind or windsave2table
    version, hash = putil.get_python_version(PY_VERSION, VERBOSE)
    with open("version", "w") as f:
        f.write("{}\n{}".format(version, hash))

    return rc


def restore_bakup_pf(root: str, wd: str):
    """
    Copy a bakup parameter file back to the original parameter file destination.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str
        The working directory to run the Python simulation in.
    """

    opf = "{}/{}.pf".format(wd, root)
    bak = opf + ".bak"
    copyfile(bak, opf)

    return


def go(roots: List[str], use_mpi: bool, n_cores: int) -> None:
    """
    Run the parts of the scripts requested to by run by the user.

    Parameters
    ----------
    roots: List[str]
        A list containing the root names of the Python simulations to run.
    use_mpi: bool
        If True, MPI will be used to run Python.
    n_cores: int
        If use_mpi is True, this will be the number of cores to run Python with.
    """

    n_sims = len(roots)

    if CHECK_CONVERGENCE:
        open("not_converged.txt", "w").close()

    for i, path in enumerate(roots):
        rc = 0
        root, wd = putil.parse_root_name_and_path(path)
        rutil.log("------------------------\n")
        rutil.log("     Simulation {}/{}".format(i + 1, n_sims))
        rutil.log("\n------------------------\n")
        rutil.log("Working directory ......... {}".format(wd))
        rutil.log("Python root name .......... {}\n".format(root))

        if RUN_SIMS:
            rutil.log("Running the simulation: {}\n".format(root))
            rc = python(root, wd, use_mpi, n_cores, RESUME_RUN, SPLIT_CYCLES)

        if CHECK_CONVERGENCE:
            rutil.log("Checking the convergence of the simulation:\n")
            c = check_python_convergence(root, wd)
            if c and SPLIT_CYCLES:
                rc = python(root, wd, use_mpi, n_cores, True, True, True)
                restore_bakup_pf(root, wd)
            elif not c and SPLIT_CYCLES:
                print("Simulation has not converged, hence no spectral cycles will be run.")
                print("Use -sc to override this.\n")
                rutil.print_error_summary(root, wd)
                restore_bakup_pf(root, wd)
                continue

        if rc == 0:
            rutil.print_error_summary(root, wd)
        else:
            continue

        putil.run_windsave2table(wd, root, VERBOSE)

        if CREATE_PLOTS:
            rutil.log("Creating plots for the simulation\n")
            plot_model(root, wd)

        if TDE_PLOT:
            rutil.log("Creating TDE specific plot for simulation\n")
            plot_spec_tde(root, wd)

        rutil.log("Removing the data directory\n")
        prd.remove_data_dir(wd, VERBOSE)

    return


def get_pf_from_file() -> List[str]:
    """
    Read in a list of Python simulations to run from a file provided as a command
    line argument.

    Returns
    -------
    roots: List[str]
        A list containing all of the root names of the Python simulations
        to run.
    """

    assert(type(SIMS_FROM_FILE) == str)

    roots = []
    broken = []

    with open(SIMS_FROM_FILE, "r") as f:
        tmp = f.readlines()

    nsims = len(tmp)
    for i in range(nsims):
        file = tmp[i].replace("\r", "").replace("\n", "")
        if access(file, R_OK):
            roots.append(file)
        else:
            broken.append(file)

    if broken:
        rutil.log("\nSome provided parameter files could not be opened:")
        for i in range(len(broken)):
            print("\t- {}".format(broken[i]))
        rutil.log("\n------------------------")
        exit(1)

    return roots


def get_run_mode() -> None:
    """
    Determine the configuration to run the script in by parsing the command
    line for flags given by the user.
    """

    global VERBOSE
    global RUN_SIMS
    global RESUME_RUN
    global CHECK_CONVERGENCE
    global CREATE_PLOTS
    global CONV_LIMIT
    global TDE_PLOT
    global PY_VERSION
    global PY_FLAGS
    global DRY_RUN
    global NOT_QUIET
    global N_CORES
    global SPLIT_CYCLES
    global SIMS_FROM_FILE
    global POLAR

    p = argparse.ArgumentParser(description="General script to run Python simulations")
    # Different run modes
    p.add_argument("-d", action="store_true", help="Run with -s -c -p")
    p.add_argument("-s", action="store_true", help="Run simulations")
    p.add_argument("-sc", action="store_true", help="Run spectral cycles even if a simulation hasn't converged")
    p.add_argument("-r", action="store_true", help="Restart a previous run")
    p.add_argument("-c", action="store_true", help="Check the convergence of runs - calls convergence.py")
    p.add_argument("-p", action="store_true", help="Run plotting scripts")
    p.add_argument("-tde", action="store_true", help="Enable TDE plotting")
    p.add_argument("-path", action="store", help="Provide a list of directories of Python parameter files to run")
    # Script Parameters
    p.add_argument("-py_ver", type=str, action="store", help="Name of the Python executable")
    p.add_argument("-f", "--py_flags", type=str, action="store", help="Runtime flags to pass to Python")
    p.add_argument("-clim", type=float, action="store", help="The convergence limit: c_value < 1")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose outputting")
    p.add_argument("--dry", action="store_true", help="Print the simulations found and exit")
    p.add_argument("-q", action="store_true", help="Enable quiet mode")
    p.add_argument("-n_cores", action="store", help="The number of processor cores to run Python with")
    # Plot Parameters
    p.add_argument("-polar", action="store_true", help="Plot using a polar projection")

    args = p.parse_args()

    do_something = False

    # Set up different run modes
    if args.d:
        RUN_SIMS = True
        CHECK_CONVERGENCE = True
        CREATE_PLOTS = True
        do_something = True
    if args.verbose:
        VERBOSE = True
    if args.s:
        RUN_SIMS = True
        CHECK_CONVERGENCE = True
        do_something = True
    if args.path:
        SIMS_FROM_FILE = args.path
    if args.c:
        CHECK_CONVERGENCE = True
        do_something = True
    if args.sc:
        SPLIT_CYCLES = False
    if args.r:
        RESUME_RUN = True
        RUN_SIMS = True
        do_something = True
    if args.p:
        CREATE_PLOTS = True
        do_something = True
    if args.tde:
        TDE_PLOT = True
        do_something = True

    # Set up script parameters
    if args.py_ver:
        PY_VERSION = args.py_ver
    if args.py_flags:
        PY_FLAGS = args.py_flags
    if args.clim:
        if 0 < args.clim < 1:
            CONV_LIMIT = args.clim
        else:
            rutil.log("Invalid value for convergence limit {}".format(args.clim))
            exit(1)
    if args.dry:
        DRY_RUN = True
        do_something = True
    if args.q:
        NOT_QUIET = False
    if args.n_cores:
        N_CORES = int(args.n_cores)

    # Set up plotting parameters
    if args.polar:
        POLAR = True

    rutil.log("------------------------\n")
    rutil.log("Python version ................... {}".format(PY_VERSION))
    rutil.log("Run simulations .................. {}".format(RUN_SIMS))
    rutil.log("Split cycles ..................... {}".format(SPLIT_CYCLES))
    rutil.log("Resume run ....................... {}".format(RESUME_RUN))
    rutil.log("Convergence limit ................ {}".format(CONV_LIMIT))
    rutil.log("Number of cores .................. {}".format(N_CORES))
    rutil.log("Show convergence ................. {}".format(CHECK_CONVERGENCE))
    rutil.log("Create plots ..................... {}".format(CREATE_PLOTS))
    rutil.log("Polar projection ................. {}".format(POLAR))
    rutil.log("Plot TDE ......................... {}".format(TDE_PLOT))
    rutil.log("Don't suppress Python output ..... {}".format(NOT_QUIET))
    rutil.log("Show Verbose Output .............. {}".format(VERBOSE))

    if PY_FLAGS:
        rutil.log("\nUsing these extra python flags:\n\t{}".format(PY_FLAGS))

    if do_something is False:
        rutil.log("\nNo run mode parameter provided, there is nothing to do!\n")
        p.print_help()
        rutil.log("\n------------------------")
        exit(0)

    return


def main() -> None:
    """
    Main control function of the script.
    """

    outfname = "py_{}{:02d}{:02d}.txt".format(DATE.year, int(DATE.month), int(DATE.day))
    rutil.init_logfile(outfname)

    # Determine which routines to run for each simulation
    get_run_mode()
    if SIMS_FROM_FILE:
        pf_paths = get_pf_from_file()
    else:
        pf_paths = putil.find_pf_files()
    n_sims = len(pf_paths)
    if not n_sims:
        rutil.log("No parameter files found, nothing to do!\n")
        rutil.log("------------------------")
        exit(0)

    use_mpi, n_procs = rutil.get_num_procs(N_CORES)

    rutil.log("")

    rutil.log("The following parameter files were found:\n")
    for i in range(len(pf_paths)):
        rutil.log("{}".format(pf_paths[i]))
    rutil.log("")

    if DRY_RUN:
        rutil.log("------------------------")
        return

    # Now run Python, plotting and convergence procedures
    go(pf_paths, use_mpi, n_procs)

    rutil.log("------------------------")

    rutil.close_logfile()

    return


if __name__ == "__main__":
    main()
