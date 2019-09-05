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
import py_run_util as pru
import py_plot_util as ppu
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
        pru.log("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(wd, root)
    if POLAR:
        commands += " -p"

    pru.log(commands)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        pru.log("\n{}".format(output))
    if err:
        pru.log("\nCaptured from stderr:")
        pru.log(err)

    pru.log("")

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
        pru.log("tde_spec_plot.py not in $PATH and executable")
        return

    spec = ppu.read_spec_file(wd + root + ".spec", pandas_table=True)
    incs = ppu.spec_inclinations_pandas(spec)
    nincs = len(incs)
    print(incs)

    commands = ["cd {}; tde_spec_plot.py {}".format(wd, root)]
    for i in range(nincs):
        commands.append("cd {}; tde_spec_plot.py {} -i {}".format(wd, root, incs[i]))

    for command in commands:
        pru.log(command)
        cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = cmd.communicate()
        out = stdout.decode("utf-8")
        err = stderr.decode("utf-8")

        if VERBOSE:
            pru.log("\n{}".format(out))
        if err:
            pru.log("Captured from stderr:")
            pru.log(err)

    pru.log("")

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
    c_fraction = pru.check_run_convergence(root, wd)

    pru.log("convergence limit ............ {}".format(CONV_LIMIT))
    pru.log("actual convergence ........... {}\n".format(c_fraction))

    if 0 > c_fraction > 1:
        pru.log(ITS_A_MYSTERY)
    elif c_fraction < CONV_LIMIT:
        pru.log(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root, c_fraction))
    elif c_fraction >= CONV_LIMIT:
        pru.log(CONVERGED)
        rc = True

    pru.log("")

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
        pcp.change_python_parameter(pf, "Spectrum_cycles", "0", bakup=True, verbose=VERBOSE)
    if split_cycles and restart_from_spec:
        pcp.change_python_parameter(pf, "Spectrum_cycles", "5", bakup=False, verbose=VERBOSE)
        pcp.change_python_parameter(pf, "Photons_per_cycle", "1e6", bakup=False, verbose=VERBOSE)

    # Construct shell command to run Python and use subprocess to run
    command = "cd {}; Setup_Py_Dir; ".format(wd)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(PY_VERSION)

    if restart_run:
        command += "-r "

    if PY_FLAGS and restart_from_spec is False:
        if type(PY_FLAGS) != str:
            pru.log("The provided additional flags for Python is not a string")
            exit(1)
        command += " {} ".format(PY_FLAGS)

    # Add the root name at the end of the call to Python
    command += "{}".format(pf)
    pru.log("{}\n".format(command))

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
        pcycle = pru.process_line_output(line, pcycle, n_cores, NOT_QUIET, VERBOSE)

    if not NOT_QUIET:
        pru.log("")
    logfile.close()

    # Sometimes with Subprocess, if the output buffer is too large then subprocess
    # breaks and causes a deadlock. To get around this, one can use .communicate()
    # to flush the buffer or s/t
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        pru.log("Captured from stderr:")
        pru.log(err)
        errfname = "{}/err_{}{:02d}{:02d}.out".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
        with open(errfname, "a") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        print("Python exited with non-zero exit code: {}\n".format(rc))
        pru.print_error_summary(root, wd)

    version, hash = ppu.get_python_version(PY_VERSION, VERBOSE)
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
        rc = -1
        root, wd = ppu.get_root_wd(path)
        pru.log("------------------------\n")
        pru.log("     Simulation {}/{}".format(i + 1, n_sims))
        pru.log("\n------------------------\n")
        pru.log("Working directory ......... {}".format(wd))
        pru.log("Python root name .......... {}\n".format(root))

        if RUN_SIMS:
            pru.log("Running the simulation: {}\n".format(root))
            rc = python(root, wd, use_mpi, n_cores, RESUME_RUN, SPLIT_CYCLES)

        if CHECK_CONVERGENCE:
            pru.log("Checking the convergence of the simulation:\n")
            c = check_python_convergence(root, wd)
            if c and SPLIT_CYCLES and RUN_SIMS:
                rc = python(root, wd, use_mpi, n_cores, True, True, True)
                restore_bakup_pf(root, wd)
            elif not c and SPLIT_CYCLES and RUN_SIMS:
                print("Simulation has not converged, hence no spectral cycles will be run.")
                print("Use -sc to override this.\n")
                pru.print_error_summary(root, wd)
                restore_bakup_pf(root, wd)
                continue

        if rc == 0 or rc == 1:
            pru.print_error_summary(root, wd)
        elif rc > 0:
            continue

        if CREATE_PLOTS or RUN_SIMS:
            ppu.windsave2table(wd, root, VERBOSE)

        if CREATE_PLOTS:
            pru.log("Creating plots for the simulation\n")
            plot_model(root, wd)

        if TDE_PLOT:
            pru.log("Creating TDE specific plot for simulation\n")
            plot_spec_tde(root, wd)

        if CREATE_PLOTS or RUN_SIMS:
            pru.log("Removing the data directory\n")
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
        pru.log("\nSome provided parameter files could not be opened:")
        for i in range(len(broken)):
            print("\t- {}".format(broken[i]))
        pru.log("\n------------------------")
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
            pru.log("Invalid value for convergence limit {}".format(args.clim))
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

    pru.log("------------------------\n")
    pru.log("Python version ................... {}".format(PY_VERSION))
    pru.log("Run simulations .................. {}".format(RUN_SIMS))
    pru.log("Split cycles ..................... {}".format(SPLIT_CYCLES))
    pru.log("Resume run ....................... {}".format(RESUME_RUN))
    pru.log("Convergence limit ................ {}".format(CONV_LIMIT))
    pru.log("Number of cores .................. {}".format(N_CORES))
    pru.log("Show convergence ................. {}".format(CHECK_CONVERGENCE))
    pru.log("Create plots ..................... {}".format(CREATE_PLOTS))
    pru.log("Polar projection ................. {}".format(POLAR))
    pru.log("Plot TDE ......................... {}".format(TDE_PLOT))
    pru.log("Don't suppress Python output ..... {}".format(NOT_QUIET))
    pru.log("Show Verbose Output .............. {}".format(VERBOSE))

    if PY_FLAGS:
        pru.log("\nUsing these extra python flags:\n\t{}".format(PY_FLAGS))

    if do_something is False:
        pru.log("\nNo run mode parameter provided, there is nothing to do!\n")
        p.print_help()
        pru.log("\n------------------------")
        exit(0)

    return


def main() -> None:
    """
    Main control function of the script.
    """

    outfname = "py_{}{:02d}{:02d}.txt".format(DATE.year, int(DATE.month), int(DATE.day))
    pru.init_logfile(outfname)

    # Determine which routines to run for each simulation
    get_run_mode()
    if SIMS_FROM_FILE:
        pf_paths = get_pf_from_file()
    else:
        pf_paths = ppu.find_pf_files()
    n_sims = len(pf_paths)
    if not n_sims:
        pru.log("No parameter files found, nothing to do!\n")
        pru.log("------------------------")
        exit(0)

    use_mpi, n_procs = pru.get_num_procs(N_CORES)

    pru.log("\nThe following parameter files were found:\n")
    for i in range(len(pf_paths)):
        pru.log("{}".format(pf_paths[i]))
    pru.log("")

    if DRY_RUN:
        pru.log("------------------------")
        return

    # Now run Python, plotting and convergence procedures
    go(pf_paths, use_mpi, n_procs)

    pru.log("------------------------")

    pru.close_logfile()

    return


if __name__ == "__main__":
    main()
