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

usage: py_run.py [-h] [-d] [-r] [-resume] [-c] [-p] [-tde] [-py_ver PY_VER]
                 [-clim CLIM] [-o] [-dry] [-q] [-n_cores N_CORES] [-wmin WMIN]
                 [-wmax WMAX]

Enable or disable features

optional arguments:
  -h, --help        show this help message and exit
  -d                Run with -r -c -p
  -r                Run simulations
  -resume           Restart a previous run
  -c                Check the convergence of runs - calls convergence.py
  -p                Run plotting scripts
  -tde              Enable TDE plotting
  -py_ver PY_VER    Name of the Python executable
  -clim CLIM        The convergence limit: c_value < 1
  -o                Verbose outputting
  -dry              Print the simulations found and exit
  -q                Enable quiet mode
  -n_cores N_CORES  The number of processor cores to run Python with
  -wmin WMIN        The smallest wavelength to display
  -wmax WMAX        The largest wavelength to display

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

DATE = datetime.datetime.now()     # The current data and time at launch of the script
VERSION = "py"                     # The name of the Python executable to use
N_CORES = 0                        # The number of processor cores to run Python over
RUN_SIMS = False                   # Flag for running a simulation
RESUME_RUN = False                 # Flag for resuming Python simulations
DRY_RUN = False                    # Flag for printing run information and then exiting
SHOW_OUTPUT = False                # If set to True, all stdout output will be shown
NOT_QUIET = True                   # If set to False, most output will not be shown
SHOW_CONVERGENCE = False           # Flag for checking the convergence of Python simulations
CLIM = 0.90                        # The limit at which a simulation is considered converged
CREATE_PLOTS = False               # Flag for running py_plot at the end of a simulation
TDE_PLOT = False                   # Flag for enabling overplotting of iPTF15af
WMIN = None                        # Smallest wavelength to plot
WMAX = None                        # Largest wavelength to plot


def get_run_mode()->None:
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

    global SHOW_OUTPUT
    global RUN_SIMS
    global RESUME_RUN
    global SHOW_CONVERGENCE
    global CREATE_PLOTS
    global CLIM
    global TDE_PLOT
    global VERSION
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
    p.add_argument("-o", action="store_true", help="Verbose outputting")
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
        SHOW_CONVERGENCE = True
        CREATE_PLOTS = True
        do_something = True
    if args.o:
        SHOW_OUTPUT = True
    if args.r:
        RUN_SIMS = True
        do_something = True
    if args.resume:
        RESUME_RUN = True
    if args.c:
        SHOW_CONVERGENCE = True
        do_something = True
    if args.p:
        CREATE_PLOTS = True
        do_something = True
    if args.tde:
        TDE_PLOT = True
    # Set up script parameters
    if args.py_ver:
        VERSION = args.py_ver
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
    print("Python version ............. {}".format(VERSION))
    print("Show Verbose Output ........ {}".format(SHOW_OUTPUT))
    print("Run Simulations ............ {}".format(RUN_SIMS))
    print("Resume Run ................. {}".format(RESUME_RUN))
    print("Show Convergence ........... {}".format(SHOW_CONVERGENCE))
    if SHOW_CONVERGENCE:
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


def py_run(root_name: str, work_dir: str, use_mpi: bool, n_cores: int, py: str = "py", resume: bool = False,
           not_quiet: bool = False, show_all: bool = False) ->Tuple[list, list]:
    """
    Function to control running and logging a Python simulation.

    Parameters
    ----------
    root_name: str
        The root name of the Python simulation
    work_dir: str
        The directory of the Python simulation.
    use_mpi: bool
        If set to True, MPI will be used to parallelise Python.
    n_cores: int
        The number of cores to parallelise Python over.
    py: str, optional
        The name of the Python executable to use.
    resume: bool, optional
        If set to True, Python will resume rather than start a run from scratch
    not_quiet: bool, optional
        If set to True, there will be less information printed to screen
    show_all: bool, optional
        If set to True, the Python output will be printed to screen

    Returns
    -------
    lines: list of strings
        The Python output which is sent to stdout.
    err: list of strings
        The Python output which is sent to stderr.
    """

    outf_name = "{}/{}_{}{}{}.txt".format(work_dir, root_name, DATE.year, DATE.month, DATE.day)
    outf = open(outf_name, "w")

    # Construct shell command to run Python and use subprocess to run
    command = ""
    pf = root_name + ".pf"
    command += "cd {}; Setup_Py_Dir; ".format(work_dir)
    if use_mpi:
        command += "mpirun -n {} ".format(n_cores)
    command += "{} ".format(py)
    if resume:
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
        spec_cycle = py_run_util.process_line_output(line, spec_cycle, not_quiet, show_all)
    print("")
    outf.close()

    # Have to do this in case the buffer is too large and causes a deadlock
    pystdout, pystderr = cmd.communicate()
    err = pystderr.decode("utf-8")
    if err:
        print("Captured from stderr:")
        print(err)
        errfname = "{}/err_{}{}{}.out"\
            .format(work_dir, root_name, DATE.year, DATE.month, DATE.day)
        with open(errfname, "w") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        raise CalledProcessError(rc, command)

    # Append a file containing the Python version and commit hash
    version, hash = py_plot_util.get_python_version(py, show_all)
    with open("version", "w") as f:
        f.write("{}\n{}".format(version, hash))

    return lines, err


def check_python_convergence(root_name:str, work_dir:str, convergence_limit:float)->Union[float, int]:
    """
    Check the convergence of a Python simulation by reading the diag file.

    Parameters
    ----------
    root_name: str
        The root name of the Python simulation.
    work_dir: str
        The directory containing the Python simulation
    convergence_limit: float
        The limit at which a Python simulation is considered considered converged.

    Returns
    -------
    convergence_fraction: float or int
        The fraction of converged cells in the Python simulation. If this is -1, this indicates an error.
    """

    if 0 > convergence_limit > 1:
        print("py_run.check_python_convergence: convergence_limit ({}) should be between 0 and 1".
              format(convergence_limit))
        return -1

    convergence_fraction = py_run_util.check_convergence(root_name, work_dir)
    print("convergence limit ............ {}".format(convergence_limit))
    print("actual convergence ........... {}\n".format(convergence_fraction))

    if convergence_fraction < convergence_limit:
        print(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(work_dir, root_name, convergence_fraction))
    elif convergence_fraction >= convergence_limit:
        print(CONVERGED)
    else:
        print(ITS_A_MYSTERY)
    print("")

    return convergence_fraction


def plot_python_output(root_name: str, work_dir: str, extra_commands: str = None, plot_tde: bool = False,
                       wmin: Union[float, int] = None, wmax: Union[float, int] = None, show_output: bool = False)->None:
    """
    Call py_plot.py and plot the output from a Python simulation.

    Parameters
    ----------
    root_name: str
        The root name of the Python simulation
    work_dir: str
        The directory containing the Python simulation
    show_output: bool, optional
        If set to True, output from the plotting script will be shown
    plot_tde: bool, optional
        If set to True, the iPTF15af UV spectrum will be overplotted
    extra_commands: str, optional
        Provide extra commands to the plotting script
    wmin: float or int, optional
        The smallest wavelength to plot
    wmax: float or int, optional
        The largest wavelength to plot

    Returns
    -------
    None
    """

    path = which("py_plot.py")
    if not path:
        print("py_plot.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot.py {}".format(work_dir, root_name)
    if plot_tde:
        commands += " -tde"
    if wmin:
        commands += " -wmin {}".format(WMIN)
    if wmax:
        commands += " -wmax {}".format(WMAX)
    if extra_commands:
        commands += "{}".format(extra_commands)
    print(commands)

    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)
    if err:
        print("Captured from stderr:")
        print(err)
    print("")

    return


def run_python_etc(pf_paths: List[str], n_sims: int, use_mpi: bool, n_cores: int, run_sims: bool = False,
                   resume_runs: bool = False, show_convergence: bool = False, clim: float = 0.90,
                   create_plots: bool = False, tde_plot: bool = False, wmin: float = None, wmax: float = None,
                   not_quiet: bool = True, show_output: bool = False) ->None:
    """
    Execute all of the different commands...

    Parameters
    ----------
    pf_paths: list of strings
        The directories of the Python simulations
    n_sims: int
        The number of simulations to be run
    use_mpi: bool
        If True, MPI will be used to parallelise Python
    n_cores: int
        The number of cores to run Python with
    clim: float
        The fraction of cells required to be converged for a Python simulation to be considered converged
    run_sims: bool, optional
        If True, Python will be called to run the simulation
    resume_runs: bool, optional
        If True, Python will resume runs rather than starting from scratch
    show_output: bool, optional
        If True, all of the output from the various commands will be printed to screen
    show_convergence: bool, optional
        If True, the convergence of the simulation will be checked
    create_plots: bool, optional
        If True, plots of the Python simulations will be created
    tde_plot: bool, optional
        If True, the iPTF15af UV spectrum will be overplotted on plotted spectra
    not_quiet: bool, optional
        If False, no command output will be printed to screen
    wmin: float, optional
        The smallest wavelength to plot
    wmax: float, optional
        The largest wavelength to plot

    Returns
    -------
    None
    """

    if show_convergence:
        open("not_converged.txt", "w").close()

    # Open up an output file
    outfname = "py_run_{}{}{}.txt".format(DATE.year, DATE.month, DATE.day)
    f = open(outfname, "w")

    for i, path in enumerate(pf_paths):
        root_name, pf_relative_path = py_plot_util.get_root_name_and_path(path)

        # Print some details to screen
        print("--------------------------\n")
        print("Simulation {}/{}".format(i+1, n_sims))
        print("\n--------------------------\n")
        print("Working directory ......... {}".format(pf_relative_path))
        print("Python root name .......... {}\n".format(root_name))

        # Print some details to the output file
        f.write("--------------------------\n")
        f.write("Simulation {}/{}\n".format(i+1, n_sims))
        f.write("--------------------------\n")
        f.write("Working directory ......... {}\n".format(pf_relative_path))
        f.write("Python root name .......... {}\n".format(root_name))

        if run_sims:
            print("Running the simulation: {}\n".format(root_name))
            python_output, python_err = py_run(root_name, pf_relative_path, use_mpi, n_cores, resume=resume_runs,
                                               not_quiet=not_quiet, show_all=show_output)
            # Write the Python output to file
            for line in python_output:
                f.write("{}\n".format(line))
            if python_err:
                f.write(python_err)

        if show_convergence:
            print("Checking the convergence of the simulation:\n")
            convergence_fraction = check_python_convergence(root_name, pf_relative_path, clim)

            # Write convergence information to file
            if convergence_fraction < clim:
                f.write("{}\n".format(NOT_CONVERGED))
                f.write("cvalue {} < clim {}\n".format(convergence_fraction, clim))
            elif convergence_fraction >= clim:
                f.write("{}\n".format(CONVERGED))
                f.write("cvalue {} >= clim {}\n".format(convergence_fraction, clim))
            else:
                f.write("{}\n".format(ITS_A_MYSTERY))

        # End of output to screen, don't write plot output to file
        f.write("\n--------------------------\n")

        if create_plots:
            print("Creating plots for the simulation\n")
            plot_python_output(root_name, pf_relative_path, tde_plot, wmin=wmin, wmax=wmax, show_output=show_output)

    f.close()

    return


def main()->None:
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
    run_python_etc(pf_paths, n_sims, use_mpi, n_procs, RUN_SIMS, RESUME_RUN, SHOW_CONVERGENCE, CLIM, CREATE_PLOTS,
                   TDE_PLOT, WMIN, WMAX, NOT_QUIET, SHOW_OUTPUT)

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
