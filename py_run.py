#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Run a batch of Python simulations. Currently, the script will search recursively
for Python pfs (disregarding anything which is py_wind.pf or .out.pf files) and
execute a number of commands depending on what is requested by the user.

Currently, the script is able to:
    - Run simulations with -r
    - Resume simulations with -resume
    - Check the convergence of simulations with -c
    - Plot output of simulations with -p
    - Default action -d which is -r -c and -p

The script can also be run in a directory containing only one Python pf.

TODO: flexible flag choices for py_plot.py - input option for the script
"""

import argparse
import datetime
import py_rm_data
import py_run_util
import py_plot_util
from sys import exit
from shutil import which
from typing import Union, List
from subprocess import Popen, PIPE, CalledProcessError

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
    """

DATE = datetime.datetime.now()
PY_VERSION = "py"
PY_FLAGS = None
DRY_RUN = False
RUN_SIMS = False
RESUME_RUN = False
N_CORES = 0
CHECK_CONVERGENCE = False
CLIM = 0.85
CREATE_PLOTS = False
TDE_PLOT = False
WMIN = None
WMAX = None
NOT_QUIET = True
VERBOSE = False
SPEC_OVERRIDE = False
SIMS_FROM_FILE = False
POLAR = True


def get_run_mode() -> None:
    """
    Parse the different run modes for the script from the command line. The run
    options are returned via global variables, because why not? :^)
    """

    global VERBOSE
    global RUN_SIMS
    global RESUME_RUN
    global CHECK_CONVERGENCE
    global CREATE_PLOTS
    global CLIM
    global TDE_PLOT
    global PY_VERSION
    global PY_FLAGS
    global WMIN
    global WMAX
    global DRY_RUN
    global NOT_QUIET
    global N_CORES
    global SPEC_OVERRIDE
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
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to display")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to display")
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
        SPEC_OVERRIDE = True
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
    if args.clim:  # Bad variable names here :-) clim is from the arg list
        if 0 < args.clim < 1:
            CLIM = args.clim
        else:
            py_run_util.log("Invalid value of c_value {}".format(args.c_value))
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
    if args.polar:
        POLAR = True

    py_run_util.log("------------------------\n")
    py_run_util.log("Python version ................... {}".format(PY_VERSION))
    py_run_util.log("Run simulations .................. {}".format(RUN_SIMS))
    py_run_util.log("Spectral cycles .................. {}".format(SPEC_OVERRIDE))
    py_run_util.log("Resume run ....................... {}".format(RESUME_RUN))
    py_run_util.log("Convergence limit ................ {}".format(CLIM))
    py_run_util.log("Number of cores .................. {}".format(N_CORES))
    py_run_util.log("Show convergence ................. {}".format(CHECK_CONVERGENCE))
    py_run_util.log("Create plots ..................... {}".format(CREATE_PLOTS))
    py_run_util.log("Polar projection ................. {}".format(POLAR))
    py_run_util.log("wmin ............................. {}".format(WMIN))
    py_run_util.log("wmax ............................. {}".format(WMAX))
    py_run_util.log("Plot TDE ......................... {}".format(TDE_PLOT))
    py_run_util.log("Don't suppress Python output ..... {}".format(NOT_QUIET))
    py_run_util.log("Show Verbose Output .............. {}".format(VERBOSE))

    if PY_FLAGS:
        py_run_util.log("\nUsing these extra python flags:\n\t{}".format(PY_FLAGS))

    if do_something is False:
        py_run_util.log("\nNo run mode parameter provided, there is nothing to do!\n")
        p.print_help()
        py_run_util.log("\n------------------------")
        exit(0)

    return


def py_run(root: str, wd: str, use_mpi: bool, n_cores: int, spec_cycle: bool = False) -> None:
    """
    Function to control running and py_run_util.logging a Python simulation.

    If spec_cycle is False, then the script will assume that Python is being run
    for just the ionisation cycles therefore the number of spectrum cycles is
    set to 0. However, if it is True, then it is assumed that all ionisation
    cycles have been completed and the model is converged. Hence the model will
    be restarted from the spectrum cycles and the photon number will be normally
    be reduced to 1e6.

    TODO: holy fuck this is complicated now.. sort it out Ed
    """

    # if spec_cycle:
    #     py_run_util.log("Spectral cycles\n--------------\n")
    # else:
    #     py_run_util.log("Ionisation cycles\n----------------\n")

    pf = root + ".pf"

    # just_spec = False
    # if spec_cycle:
    #     py_run_util.change_parameter(wd + pf, "Spectrum_cycles", "5", VERBOSE)
    #     if RESUME_RUN:
    #         py_run_util.change_parameter(wd + pf, "Photons_per_cycle", "1e5", VERBOSE)
    #         just_spec = True
    #     else:
    #         py_run_util.change_parameter(wd + pf, "Photons_per_cycle", "1e2", VERBOSE)
    # else:
    #     if not spec_cycle:
    #         py_run_util.change_parameter(wd + pf, "Spectrum_cycles", "0", VERBOSE)


    outf_name = "{}/{}_{}{:02d}{:02d}.txt".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
    outf = open(outf_name, "a")

    # Construct shell command to run Python and use subprocess to run
    command = ""
    command += "cd {}; Setup_Py_Dir; ".format(wd)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(PY_VERSION)
    if RESUME_RUN:
        command += "-r "
    if PY_FLAGS:
        if type(PY_FLAGS) != str:
            py_run_util.log("PY_FLAGS is not a string")
            exit(1)
        command += " {} ".format(PY_FLAGS)
    command += "{} ".format(pf)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    py_run_util.log("{}\n".format(command))

    # Real time output of Python
    pcycle = False
    for stdout_line in iter(cmd.stdout.readline, ""):
        if not stdout_line:
            break
        line = stdout_line.decode("utf-8").replace("\n", "")
        outf.write("{}\n".format(line))
        pcycle = py_run_util.process_line_output(line, pcycle, n_cores, NOT_QUIET, VERBOSE)
    py_run_util.log("")

    outf.close()

    # Have to do this in case the buffer is too large and causes a deadlock
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        py_run_util.log("Captured from stderr:")
        py_run_util.log(err)
        errfname = "{}/err_{}{:02d}{:02d}.out".format(wd, root, DATE.year, int(DATE.month), int(DATE.day))
        with open(errfname, "a") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        print("Python exited with non-zero exit code: {}".format(rc))
        # raise CalledProcessError(rc, command)

    # Append a file containing the Python version and commit hash
    version, hash = py_plot_util.get_python_version(PY_VERSION, VERBOSE)
    with open("version", "w") as f:
        f.write("{}\n{}".format(version, hash))

    return


def check_python_convergence(root: str, wd: str) -> Union[float, int]:
    """
    Check the convergence of a Python simulation by reading the diag file.
    """

    if 0 > CLIM > 1:
        py_run_util.log("py_run.check_python_convergence: convergence_limit {} should be between 0 and 1".format(CLIM))
        return -1

    convergence_fraction = py_run_util.check_convergence(root, wd)
    py_run_util.log("convergence limit ............ {}".format(CLIM))
    py_run_util.log("actual convergence ........... {}\n".format(convergence_fraction))

    if convergence_fraction < 0:
        py_run_util.log(ITS_A_MYSTERY)
    elif convergence_fraction > 1:
        py_run_util.log(ITS_A_MYSTERY)
    elif convergence_fraction < CLIM:
        py_run_util.log(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root, convergence_fraction))
    elif convergence_fraction >= CLIM:
        py_run_util.log(CONVERGED)
    else:
        py_run_util.log(ITS_A_MYSTERY)
    py_run_util.log("")

    return convergence_fraction


def plot_python_output(root_name: str, wd: str, extra_commands: str = None) -> None:
    """
    Call py_plot.py and plot the output from a Python simulation.
    """

    path = which("py_plot.py")
    if not path:
        py_run_util.log("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(wd, root_name)
    if WMIN:
        commands += " -wmin {}".format(WMIN)
    if WMAX:
        commands += " -wmax {}".format(WMAX)
    if POLAR:
        commands += " -p"
    if extra_commands:
        commands += "{}".format(extra_commands)
    py_run_util.log(commands)

    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        py_run_util.log("\n{}".format(output))
    if err:
        py_run_util.log("\nCaptured from stderr:")
        py_run_util.log(err)
    py_run_util.log("")

    return


def plot_tde(root: str, wd: str) -> None:
    """
    Do some tde plotting :-)
    """

    path = which("tde_spec_plot.py")
    if not path:
        py_run_util.log("tde_spec_plot.py not in $PATH and executable")
        return

    command = "cd {}; tde_spec_plot.py {} -tde iPTF15af -l".format(wd, root)
    py_run_util.log(command)

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if VERBOSE:
        py_run_util.log("\n{}".format(out))
    if err:
        py_run_util.log("Captured from stderr:")
        py_run_util.log(err)
    py_run_util.log("")

    return


def run_python_etc(pf_paths: List[str], n_sims: int, use_mpi: bool, n_cores: int) -> None:
    """
    Execute all of the different commands to run a Python simulation
    """

    if CHECK_CONVERGENCE:
        open("not_converged.txt", "w").close()

    for i, path in enumerate(pf_paths):
        root, wd = py_plot_util.parse_root_name_and_path(path)

        py_run_util.log("------------------------\n")
        py_run_util.log("     Simulation {}/{}".format(i + 1, n_sims))
        py_run_util.log("\n------------------------\n")
        py_run_util.log("Working directory ......... {}".format(wd))
        py_run_util.log("Python root name .......... {}\n".format(root))

        # run_spec_cycles = False
        if RUN_SIMS:
            py_run_util.log("Running the simulation: {}\n".format(root))
            py_run(root, wd, use_mpi, n_cores, SPEC_OVERRIDE)

            # if SPEC_OVERRIDE:
            #     run_spec_cycles = True

        if CHECK_CONVERGENCE:
            py_run_util.log("Checking the convergence of the simulation:\n")
            convergence_fraction = check_python_convergence(root, wd)

            # if convergence_fraction >= CLIM and RUN_SIMS and not run_spec_cycles:
            #     py_run(root, wd, use_mpi, n_cores, True)
            # elif convergence_fraction < CLIM and RUN_SIMS:
            #     if SPEC_OVERRIDE and not run_spec_cycles:
            #         py_run_util.log("As the simulation hasn't converged, but SPEC_OVERRIDE is True, spec cycles will be"
            #                         " run as normal\n")
            #         py_run(root, wd, use_mpi, n_cores, True)
            #     elif RUN_SIMS:
            #         py_run_util.log("As the simulation hasn't converged, no spectrum cycles will be run\n")

        if CREATE_PLOTS:
            py_run_util.log("Creating plots for the simulation\n")
            plot_python_output(root, wd, None)

        if TDE_PLOT:
            py_run_util.log("Creating TDE specific plot for simulation\n")
            plot_tde(root, wd)

        py_rm_data.remove_data_dir(wd, VERBOSE)

    return


def get_pf_from_file() -> List[str]:
    """
    Read in a list of Python simulations to run from a file provided as a command
    line argument.
    """

    assert(type(SIMS_FROM_FILE) == str)

    pf_paths = []
    with open(SIMS_FROM_FILE, "r") as f:
        tmp = f.readlines()
    nsims = len(tmp)
    for i in range(nsims):
        pf_paths.append(tmp[i].replace("\r", "n").replace("\n", ""))

    return pf_paths


def main() -> None:
    """
    Main control function of the script
    """

    outfname = "py_{}{:02d}{:02d}.txt".format(DATE.year, int(DATE.month), int(DATE.day))
    py_run_util.init_logfile(outfname)

    # Determine which routines to run for each simulation
    get_run_mode()
    if SIMS_FROM_FILE:
        pf_paths = get_pf_from_file()
    else:
        pf_paths = py_plot_util.find_pf_files()
    n_sims = len(pf_paths)
    if not n_sims:
        py_run_util.log("No parameter files found, nothing to do!\n")
        py_run_util.log("------------------------")
        exit(0)

    use_mpi, n_procs = py_run_util.get_num_procs(N_CORES)

    py_run_util.log("")

    py_run_util.log("The following parameter files were found:\n")
    for i in range(len(pf_paths)):
        py_run_util.log("{}".format(pf_paths[i]))
    py_run_util.log("")

    if DRY_RUN:
        py_run_util.log("------------------------")
        return

    # Now run Python, plotting and convergence procedures
    run_python_etc(pf_paths, n_sims, use_mpi, n_procs)

    py_run_util.log("------------------------")

    py_run_util.close_logfile()

    return


if __name__ == "__main__":
    main()
