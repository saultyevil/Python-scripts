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


import time
import argparse
from os import access, R_OK
import datetime
from sys import exit
from shutil import which, copyfile
from typing import Union, List
from subprocess import Popen, PIPE
from PyPython import Grid, Simulation, Utils, Log, SpectrumUtils, WindUtils


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
        Log.log("{} : Ionisation Cycle ....... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("to calculate a detailed spectrum") != -1 and print_crap:
        line = line.split()
        pcycle = True
        cycle = int(line[1])  # + 1
        ncycles = line[3]
        current_time = time.strftime("%H:%M")
        Log.log("{} : Spectrum Cycle ......... {}/{}".format(current_time, cycle, ncycles))
    elif line.find("per cent") != -1 and line.find("Photon") != -1 and print_crap:
        line = line.split()
        Log.log("      : {}% of {:1.2e} photons transported".format(line[-3], int(int(line[-5]) * n_cores)))
    elif line.find("!!Check_converging:") != -1 and print_crap:
        line = line.split()
        nconverged = int(line[1])
        fconverged = line[2]
        Log.log("      : {} cells converged {}".format(nconverged, fconverged))
    elif (line.find("Completed ionization cycle") != -1 or line.find("Completed spectrum cycle") != -1) and print_crap:
        line = line.split()
        elapsed_time_seconds = float(line[-1])
        elapsed_time = datetime.timedelta(seconds=elapsed_time_seconds // 1)
        Log.log("      : Elapsed run time {} hrs:mins:secs".format(elapsed_time))
    elif line.find("photon transport completed in") != -1 and print_crap:
        line = line.split()
        transport_time_seconds = float(line[5])
        transport_time = datetime.timedelta(seconds=transport_time_seconds // 1)
        Log.log("      : Photons transported in {} hrs:mins:secs".format(transport_time))
    elif line.find("Completed entire program.") != -1 and print_crap:
        line = line.split()
        tot_run_time_seconds = float(line[-1])
        tot_run_time = datetime.timedelta(seconds=tot_run_time_seconds // 1)
        Log.log("\nSimulation completed in {} hrs:mins:secs".format(tot_run_time))

    return pcycle


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
        Log.log("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(wd, root)
    if POLAR:
        commands += " -p"

    Log.log(commands)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        Log.log("\n{}".format(output))
    if err:
        Log.log("\nCaptured from stderr:")
        Log.log(err)

    Log.log("")

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
        Log.log("{}:{}: tde_spec_path.py not in $PATH and executable".format(__file__, plot_spec_tde.__name__))
        return

    incs = SpectrumUtils.spec_inclinations([wd + root + ".spec"])
    nincs = len(incs)

    commands = ["cd {}; tde_spec_plot.py {}".format(wd, root)]
    for i in range(nincs):
        commands.append("cd {}; tde_spec_plot.py {} -i {}".format(wd, root, incs[i]))

    for command in commands:
        Log.log(command)
        cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
        stdout, stderr = cmd.communicate()
        out = stdout.decode("utf-8")
        err = stderr.decode("utf-8")

        if VERBOSE:
            Log.log("\n{}".format(out))
        if err:
            Log.log("Captured from stderr:")
            Log.log(err)

    Log.log("")

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
    c_fraction = Simulation.check_convergence(root, wd)

    Log.log("convergence limit ............ {}".format(CONV_LIMIT))
    Log.log("actual convergence ........... {}\n".format(c_fraction))

    if 0 > c_fraction > 1:
        Log.log(ITS_A_MYSTERY)
    elif c_fraction < CONV_LIMIT:
        Log.log(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root, c_fraction))
    elif c_fraction >= CONV_LIMIT:
        Log.log(CONVERGED)
        rc = True

    Log.log("")

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
        Grid.change_parameter(pf, "Spectrum_cycles", "0", backup=True, verbose=VERBOSE)
    if split_cycles and restart_from_spec:
        Grid.change_parameter(pf, "Spectrum_cycles", "5", backup=False, verbose=VERBOSE)
        Grid.change_parameter(pf, "Photons_per_cycle", "1e6", backup=False, verbose=VERBOSE)

    # Construct shell command to run Python and use subprocess to run
    command = "cd {}; Setup_Py_Dir; ".format(wd)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(PY_VERSION)

    if restart_run:
        command += "-r "

    if PY_FLAGS and restart_from_spec is False:
        if type(PY_FLAGS) != str:
            Log.log("The provided additional flags for Python is not a string")
            exit(1)
        command += " {} ".format(PY_FLAGS)

    # Add the root name at the end of the call to Python
    command += "{}".format(pf)
    Log.log("{}\n".format(command))

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
        pcycle = process_line_output(line, pcycle, n_cores, NOT_QUIET, VERBOSE)

    if not NOT_QUIET:
        Log.log("")
    logfile.close()

    # Sometimes with Subprocess, if the output buffer is too large then subprocess
    # breaks and causes a deadlock. To get around this, one can use .communicate()
    # to flush the buffer or s/t
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        Log.log("Captured from stderr:")
        Log.log(err)
        errfname = "{}/err_{}{:02d}{:02d}.out".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
        with open(errfname, "a") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        print("Python exited with non-zero exit code: {}\n".format(rc))
        Simulation.error_summary(root, wd, VERBOSE)

    version, hash = Utils.get_python_version(PY_VERSION, VERBOSE)
    with open("version", "w") as f:
        f.write("{}\n{}".format(version, hash))

    return rc


def restore_bakup_pf(root: str, wd: str):
    """
    Copy a backup parameter file back to the original parameter file
    destination.

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
        root, wd = Utils.split_root_directory(path)
        Log.log("------------------------\n")
        Log.log("     Simulation {}/{}".format(i + 1, n_sims))
        Log.log("\n------------------------\n")
        Log.log("Working directory ......... {}".format(wd))
        Log.log("Python root name .......... {}\n".format(root))

        if RUN_SIMS:
            Log.log("Running the simulation: {}\n".format(root))
            rc = python(root, wd, use_mpi, n_cores, RESUME_RUN, SPLIT_CYCLES)

        if CHECK_CONVERGENCE:
            Log.log("Checking the convergence of the simulation:\n")
            c = check_python_convergence(root, wd)
            if c and SPLIT_CYCLES and RUN_SIMS:
                rc = python(root, wd, use_mpi, n_cores, True, True, True)
                restore_bakup_pf(root, wd)
            elif not c and SPLIT_CYCLES and RUN_SIMS:
                print("Simulation has not converged, hence no spectral cycles will be run.")
                print("Use -sc to override this.\n")
                Simulation.error_summary(root, wd, VERBOSE)
                restore_bakup_pf(root, wd)
                continue

        if rc == 0 or rc == 1:
            Simulation.error_summary(root, wd, VERBOSE)
        elif rc > 0:
            continue

        if CREATE_PLOTS or RUN_SIMS:
            WindUtils.windsave2table(root, wd, VERBOSE)

        if CREATE_PLOTS:
            Log.log("Creating plots for the simulation\n")
            plot_model(root, wd)

        if TDE_PLOT:
            Log.log("Creating TDE specific plot for simulation\n")
            plot_spec_tde(root, wd)

        if CREATE_PLOTS or RUN_SIMS:
            Log.log("Removing the data directory\n")
            Utils.remove_data_sym_links(wd, VERBOSE)

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
        Log.log("\nSome provided parameter files could not be opened:")
        for i in range(len(broken)):
            print("\t- {}".format(broken[i]))
        Log.log("\n------------------------")
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
            Log.log("Invalid value for convergence limit {}".format(args.clim))
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

    Log.log("------------------------\n")
    Log.log("Python version ................... {}".format(PY_VERSION))
    Log.log("Run simulations .................. {}".format(RUN_SIMS))
    Log.log("Split cycles ..................... {}".format(SPLIT_CYCLES))
    Log.log("Resume run ....................... {}".format(RESUME_RUN))
    Log.log("Convergence limit ................ {}".format(CONV_LIMIT))
    Log.log("Number of cores .................. {}".format(N_CORES))
    Log.log("Show convergence ................. {}".format(CHECK_CONVERGENCE))
    Log.log("Create plots ..................... {}".format(CREATE_PLOTS))
    Log.log("Polar projection ................. {}".format(POLAR))
    Log.log("Plot TDE ......................... {}".format(TDE_PLOT))
    Log.log("Don't suppress Python output ..... {}".format(NOT_QUIET))
    Log.log("Show Verbose Output .............. {}".format(VERBOSE))

    if PY_FLAGS:
        Log.log("\nUsing these extra python flags:\n\t{}".format(PY_FLAGS))

    if do_something is False:
        Log.log("\nNo run mode parameter provided, there is nothing to do!\n")
        p.print_help()
        Log.log("\n------------------------")
        exit(0)

    return


def main() -> None:
    """
    Main control function of the script.
    """

    outfname = "py_{}{:02d}{:02d}.txt".format(DATE.year, int(DATE.month), int(DATE.day))
    Log.init_logfile(outfname)

    # Determine which routines to run for each simulation
    get_run_mode()
    if SIMS_FROM_FILE:
        pf_paths = get_pf_from_file()
    else:
        pf_paths = Utils.find_parameter_files()
    n_sims = len(pf_paths)
    if not n_sims:
        Log.log("No parameter files found, nothing to do!\n")
        Log.log("------------------------")
        exit(0)

    use_mpi = False
    n_procs = Utils.get_cpu_count()
    if n_procs > 1:
        use_mpi = True

    Log.log("\nThe following parameter files were found:\n")
    for i in range(len(pf_paths)):
        Log.log("{}".format(pf_paths[i]))
    Log.log("")

    if DRY_RUN:
        Log.log("------------------------")
        return

    # Now run Python, plotting and convergence procedures
    go(pf_paths, use_mpi, n_procs)

    Log.log("------------------------")

    Log.close_logfile()

    return


if __name__ == "__main__":
    main()
