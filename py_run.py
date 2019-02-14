#!/usr/bin/env python3

"""
Automatically run Python simulations, plot and check convergence by searching
for .pf files recursively from the calling directory

There are no required arguments to provide the script. However, if no "run time"
parameters are provided, then the script will exit as there is nothing to do.

The switches which can be used are as follows:

  -h, --help      show this help message and exit
  -d              Run with -r -c -p
  -r              Run simulations
  -c              Check the convergence of runs - calls convergence.py
  -p              Run plotting scripts
  -tde            Enable TDE plotting
  -py_ver PY_VER  Name of the Python executable
  -clim CLIM      The convergence limit: c_value < 1
  -o              Verbose outputting
  -dry            Print the simulations found and exit
  -wmin WMIN      The smallest wavelength to display
  -wmax WMAX      The largest wavelength to display


Note that -d, -p, -c, -dry or -r are required for the script to do something.

Example usage:
    python py_run.py -r -p -tde
    python py_run.py -d -tde -WMIN 200 -WMAX 5000 -o

TODO: be able to run simulations from file input
TODO: flexible flag choices for spec_plot.py
TODO: specify the number of MPI processes
"""


import time
import py_util
import argparse
import datetime
from sys import exit
from shutil import which
from platform import system
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, CalledProcessError


CLIM = 0.85
VERSION = "py"
SHOW_OUTPUT = False
RUN_SIMS = False
RESUME_RUN = False
SHOW_CONVERGENCE = False
CREATE_PLOTS = False
TDE_PLOT = False
WMIN = None
WMAX = None
DATE = datetime.datetime.now()
DRY_RUN = False
NOT_QUIET = True


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


def get_run_mode():
    """
    Read command line flags to determine the mode of operation
    """

    p = argparse.ArgumentParser(description="Enable or disable features")

    # Different run modes
    p.add_argument("-d", action="store_true", help="Run with -r -c -p")
    p.add_argument("-r", action="store_true", help="Run simulations")
    p.add_argument("-resume", action="store_true", help=
        "Restart a previous run")
    p.add_argument("-c", action="store_true", help=
        "Check the convergence of runs - calls convergence.py")
    p.add_argument("-p", action="store_true", help="Run plotting scripts")
    p.add_argument("-tde", action="store_true", help=
        "Enable TDE plotting")
    # Script parameters
    p.add_argument("-py_ver", type=str, action="store", help=
        "Name of the Python executable")
    p.add_argument("-clim", type=float, action="store", help=
        "The convergence limit: c_value < 1")
    p.add_argument("-o", action="store_true", help="Verbose outputting")
    p.add_argument("-dry", action="store_true", help=
        "Print the simulations found and exit")
    p.add_argument("-q", action="store_true", help="Enable quiet mode")
    # Plot parameters
    p.add_argument("-wmin", type=float, action="store",help=
        "The smallest wavelength to display")
    p.add_argument("-wmax", type=float, action="store", help=
        "The largest wavelength to display")

    args = p.parse_args()

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


def py_run(wd, root_name, mpi, ncores):
    """
    Do a standard Python run using either a single process or multiple
    """

    pf = root_name + ".pf"  # Don't actually need the .pf at the end
    command = "cd {};".format(wd)
    if mpi:
        command += " mpirun -n {}".format(ncores)
    command += " {}".format(VERSION)
    if RESUME_RUN:
        command += " -r"
    command += " {}".format(pf)
    print("{}\n".format(command))

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    outfname = "{}/{}_{}{}{}.txt"\
        .format(wd, root_name, DATE.year, DATE.month, DATE.day)
    outf = open(outfname, "w")

    lines = []
    spec_cycle = False
    for stdout_line in iter(cmd.stdout.readline, ""):
        if not stdout_line:
            break
        line = stdout_line.decode("utf-8").replace("\n", "")
        lines.append(line)
        outf.write("{}\n\n".format(line))

        #
        # This horrid bit of logic will determine from the output which
        # ionisation and wind cycle is currently being done and print out to
        # the screen. If SHOW_OUTPUT is True, then all lines will be printed
        # to the screen instead.
        #

        if SHOW_OUTPUT:
            print(line)
        elif line.find("for defining wind") != -1 and NOT_QUIET:
            line = line.split()
            cycle = int(line[3])  # + 1
            ncycles = line[5]
            current_time = time.strftime("%H:%M")
            print("{} : Ionisation Cycle ....... {}/{}".format(current_time,
                                                               cycle, ncycles))
        elif line.find("to calculate a detailed spectrum") != -1 and NOT_QUIET:
            line = line.split()
            spec_cycle = True
            cycle = int(line[1])  # + 1
            ncycles = line[3]
            current_time = time.strftime("%H:%M")
            print("{} : Spectrum Cycle ......... {}/{}".format(current_time,
                                                               cycle, ncycles))
        elif line.find("per cent") != -1 and line.find("Photon") != -1 and NOT_QUIET:
            line = line.split()
            print("   - {}% of {} photons transported".format(line[-3],
                                                              line[-5]))
        elif line.find("!!Check_converging:") != -1 and NOT_QUIET:
            line = line.split()
            nconverged = int(line[1])
            fconverged = line[2]
            print("           ----------           ")
            print("   - {} cells converged {}".format(nconverged, fconverged))
        elif (line.find("Completed ionization cycle") != -1 or
                line.find("Completed spectrum cycle") != -1) and NOT_QUIET:
            line = line.split()
            elapsed_time_seconds = float(line[-1])
            elapsed_time = datetime.timedelta(seconds=elapsed_time_seconds // 1)
            if spec_cycle:
                print("           ----------           ")
            print("   - Elapsed time: {} hours".format(elapsed_time))
            print("           ----------           ")
        elif line.find("PHOTON TRANSPORT COMPLETED") != -1 and NOT_QUIET:
            line = line.split()
            transport_time_seconds = float(line[4])
            transport_time = datetime.timedelta(seconds=transport_time_seconds
                                                        // 1)
            print("   - Photon transport: {} hours"
                  .format(transport_time))
        elif line.find("Completed entire program.") != -1 and NOT_QUIET:
            line = line.split()
            tot_run_time_seconds = float(line[-1])
            tot_run_time = datetime.timedelta(seconds=tot_run_time_seconds // 1)
            print("\nSimulation completed in {} hours".format(tot_run_time))

    print("")

    outf.close()

    #
    # When the stdout or stderr buffer becomes too large (>4KB), cmd.wait()
    # deadlocks. We can get around this by using cmd.communicate instead, I
    # think... it seems to be working...
    #

    pystdout, pystderr = cmd.communicate()

    err = pystderr.decode("utf-8")
    if err:
        print("Captured from stderr:")
        print(err)
        errfname = "{}/err_{}{}{}.out"\
            .format(wd, root_name, DATE.year, DATE.month, DATE.day)
        with open(errfname, "w") as f:
            f.writelines(err)
    rc = cmd.returncode
    if rc:
        raise CalledProcessError(rc, command)

    return lines, err


def get_convergence(wd, root_name):
    """
    Write to the screen if the simulation has or hasn't converged
    """

    convergence_fraction = py_util.check_convergence(wd, root_name)

    print("clim ............ {}".format(CLIM))
    print("convergence ..... {}\n".format(convergence_fraction))

    if convergence_fraction < CLIM:
        print(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root_name,
                                             convergence_fraction))
    elif convergence_fraction >= CLIM:
        print(CONVERGED)
    else:
        print(ITS_A_MYSTERY)

    print("")

    return convergence_fraction


def do_py_plot_output(wd, root_name):
    """
    Execute the standard py_plot_output routine located in py_progs
    """

    path = which("py_plot_output.py")
    if not path:
        print("py_plot_output.py not in $PATH and executable")
        return

    commands = "cd {}; py_plot_output.py {} all".format(wd, root_name)
    print(commands)

    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")

    if SHOW_OUTPUT:
        print(output)

    print("")

    return


def do_spec_plot(wd, root_name):
    """
    Execute my standard spec_plot script to plot spectra
    """

    path = which("spec_plot.py")
    if not path:
        print("spec_plot.py not in $PATH and executable")
        return

    commands = "cd {}; spec_plot.py {} all".format(wd, root_name)
    if TDE_PLOT:
        commands += " -tde"
    if WMIN:
        commands += " -wmin {}".format(WMIN)
    if WMAX:
        commands += " -wmax {}".format(WMAX)
    print(commands)

    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if SHOW_OUTPUT:
        print(output)
    if err:
        print("Captured from stderr:")
        print(err)
    print("")

    return


def check_for_mpi():
    """
    Check that mpirun exists on the system
    """

    path = which("mpirun")
    if not path:
        return False

    return True


def find_number_of_physical_cores_lscpu():
    """
    Created due to multiprocessing.cpu_count() returning physical and logical
    core count. Pathetic. (I don't want to use hyperthreads as this often can
    lead to slow down for a range of HPC problems)
    """

    ncores_cmd = "lscpu | grep 'Core(s) per socket'"
    nsockets_cmd = "lscpu | grep 'Socket(s):'"
    ncores, stderr = Popen(ncores_cmd, stdout=PIPE, stderr=PIPE,
                           shell=True).communicate()
    if ncores == "":
        return 0
    ncores = int(ncores.decode("utf-8").replace("\n", "").split()[-1])
    nsockets, stderr = Popen(nsockets_cmd, stdout=PIPE, stderr=PIPE,
                             shell=True).communicate()
    if nsockets == "":
        return 0
    nsockets = int(nsockets.decode("utf-8").replace("\n", "").split()[-1])

    return ncores * nsockets


def get_num_cores():
    """
    Determine the number of cores Python should be run using
    """

    n_cores = 1
    mpi = check_for_mpi()
    if RUN_SIMS and mpi:
        # We will suffer for macOS since lscpu doesn't exist
        if system() == "Darwin":
            n_cores = cpu_count()
        else:
            n_cores = find_number_of_physical_cores_lscpu()
            if n_cores == 0:
                n_cores = cpu_count()
        print("The following simulations will be run using {} cores:\n"
              .format(n_cores))
    elif RUN_SIMS and not mpi:
        print("The follow simulations will be run in single threaded mode:\n")
    else:
        print("Processing the following runs:\n")

    return mpi, n_cores


def run_choices(par_file_paths, n_sims, mpi, n_cores):
    """
    Controls the general flow of the simulation running, plotting etc. I did
    this to make the main function more readable.
    """

    # If the not_converged file exists, delete the contents :-)
    if SHOW_CONVERGENCE:
        open("not_converged.txt", "w").close()

    # Open up an output file
    outfname = "run_{}{}{}.txt".format(DATE.year, DATE.month, DATE.day,
                                       DATE.hour, DATE.minute)
    f = open(outfname, "w")

    # Iterate over the possible simulations
    for i, path in enumerate(par_file_paths):
        root_name, pf_relative_path = py_util.get_root_name_and_path(path)

        # Print some details to screen
        print("--------------------------\n")
        print("Simulation {}/{}".format(i+1, n_sims))
        print("Working dir ........ {}".format(pf_relative_path))
        print("Root name .......... {}\n".format(root_name))
        # Print some details to the output file
        f.write("--------------------------\n")
        f.write("Simulation {}/{}\n".format(i+1, n_sims))
        f.write("Working dir ........ {}\n".format(pf_relative_path))
        f.write("Root name .......... {}\n".format(root_name))

        if RUN_SIMS:
            if RESUME_RUN:
                print("Resuming the simulation: {}\n".format(root_name))
            else:
                print("Running the simulation: {}\n".format(root_name))

            python_output, python_err = py_run(pf_relative_path, root_name, mpi,
                                               n_cores)
            f.writelines(python_output)
            if python_err:
                f.write(python_err)

        if SHOW_CONVERGENCE:
            print("Checking the convergence of the simulation:\n")

            convergence_fraction = get_convergence(pf_relative_path, root_name)

            if convergence_fraction < CLIM:
                f.write("{}\n".format(NOT_CONVERGED))
                f.write("cvalue {} < clim {}\n".format(convergence_fraction,
                                                       CLIM))
            elif convergence_fraction >= CLIM:
                f.write("{}\n".format(CONVERGED))
                f.write("cvalue {} >= clim {}\n".format(convergence_fraction,
                                                        CLIM))
            else:
                f.write("{}\n".format(ITS_A_MYSTERY))

        # End of output to screen
        f.write("\n--------------------------\n")

        if CREATE_PLOTS:
            print("Creating plots for the simulation\n")
            do_py_plot_output(pf_relative_path, root_name)
            do_spec_plot(pf_relative_path, root_name)

    f.close()

    return


def main():
    """
    Main control function
    """

    # Determine which routines to run for each simulation
    get_run_mode()
    all_par_file_paths = py_util.find_pf(ignore_out_pf=True)
    mpi, n_cores = get_num_cores()

    # Remove any py_wind parameter files which might be lurking about
    par_file_paths = []
    for i, path in enumerate(all_par_file_paths):
        if path.find("py_wind.pf") == -1:
            par_file_paths.append(path)
            print("- {}".format(path))
    print("")

    if DRY_RUN:
        print("--------------------------")
        return

    # Now run Python, plotting and convergence procedures
    n_sims = len(par_file_paths)
    run_choices(par_file_paths, n_sims, mpi, n_cores)

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
