#!/usr/bin/env python3


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
SHOW_CONVERGENCE = False
CREATE_PLOTS = False
TDE_PLOT = False
WMIN = None
WMAX = None
DATE = datetime.datetime.now()
DRY_RUN = False

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


def get_run_mode():
    """
    Read command line flags to determine the mode of operation
    """

    p = argparse.ArgumentParser(description="Enable or disable features")

    # Different run modes
    p.add_argument("-d", action="store_true", help="Run with -r -c -p")
    p.add_argument("-r", action="store_true", help="Run simulations")
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
    # Plot parameters
    p.add_argument("-wmin", type=float, action="store",help=
        "The smallest wavelength to display")
    p.add_argument("-wmax", type=float, action="store", help=
        "The largest wavelength to display")

    args = p.parse_args()

    global SHOW_OUTPUT
    global RUN_SIMS
    global SHOW_CONVERGENCE
    global CREATE_PLOTS
    global CLIM
    global TDE_PLOT
    global VERSION
    global WMIN
    global WMAX
    global DRY_RUN

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
    if args.c:
        SHOW_CONVERGENCE = True
        do_something = True
    if args.p:
        CREATE_PLOTS = True
        do_something = True
    if args.tde:
        CREATE_PLOTS = True
        TDE_PLOT = True
        do_something = True
    # Set up script parameters
    if args.py_ver:
        VERSION = args.PY_VER
    if args.clim:  # Bad variable names here :-) clim is from the arg list
        if 0 < args.c_value < 1:
            CLIM = args.clim
        else:
            print("Invalid value of c_value {}".format(args.c_value))
            exit(1)
    if args.dry:
        DRY_RUN = True
        do_something = True
    # Set up plotting parameters
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax

    print("--------------------------\n")
    print("Python version ............. {}".format(VERSION))
    print("Show Verbose Output ........ {}".format(SHOW_OUTPUT))
    print("Run Simulations ............ {}".format(RUN_SIMS))
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
    Do a standard Python run using a single process
    """

    pf = root_name + ".pf"  # Don't actually need the .pf at the end
    if mpi:
        command = "cd {}; Setup_Py_Dir; mpirun -n {} {} {}; rm data"\
            .format(wd, ncores, VERSION, pf)
    else:
        command = "cd {}; Setup_Py_Dir; {} {}".format(wd, VERSION, pf)

    print("Working dir ........ {}".format(wd))
    print("Root name .......... {}".format(root_name))
    print("\n{}\n".format(command))

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    outfname = "{}/{}_{}{}{}.out"\
        .format(wd, root_name, DATE.year, DATE.month, DATE.day)
    outf = open(outfname, "w")

    lines = []
    for stdout_line in iter(cmd.stdout.readline, ""):
        if not stdout_line:
            break
        line = stdout_line.decode("utf-8").replace("\n", "")
        lines.append(line)
        outf.write("{}\n".format(line))

        #
        # This horrid bit of logic will determine from the output which
        # ionisation and wind cycle is currently being done and print out to
        # the screen. If SHOW_OUTPUT is True, then all lines will be printed
        # to the screen instead.
        #

        if line.find("for defining wind") != -1:
            line = line.split()
            cycle = int(line[3]) + 1
            ncycles = line[5]
            print("Ionisation Cycle ... {}/{}".format(cycle, ncycles))
        elif line.find("to calculate a detailed spectrum") != -1:
            line = line.split()
            cycle = int(line[1]) + 1
            ncycles = line[3]
            print("Spectrum Cycle ..... {}/{}".format(cycle, ncycles))
        elif SHOW_OUTPUT:
            print(line)
    print("")

    outf.close()

    #
    # When the stdout or stderr buffer becomes too large (>4KB), cmd.wait()
    # deadlocks. We can get around this by using cmd.communicate instead, I
    # think... it seems to be working...
    #

    pystdout, pystderr = cmd.communicate()
    rc = cmd.returncode
    if rc:
        raise CalledProcessError(rc, command)
    err = pystderr.decode("utf-8")
    if err:
        print("Captured from stderr:")
        print(err)
        errfname = "{}/err_{}{}{}.out"\
            .format(wd, root_name, DATE.year, DATE.month, DATE.day)
        with open(errfname, "w") as f:
            f.writelines(err)

    return lines, err


def check_convergence(wd, root_name, show_output):
    """
    Check the convergence of the simulation
    """

    commands = "cd {}; convergence.py {}".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)

    if err:
        print("Captured from stderr:")
        print(err)

    return output


def display_convergence(conv, wd, root_name):
    """
    Write to the screen if the simulation has or hasn't converged
    """

    final_cycle = 0
    for i, word in enumerate(conv):
        if word == "!!Check_converging:":
            final_cycle = i
    if final_cycle == 0:
        return 0

    final_cycle_conv = conv[final_cycle:]
    cvalue = float(final_cycle_conv[2].replace("(", "").replace(")", ""))

    print("dir ............. {}".format(wd))
    print("root name ....... {}".format(root_name))
    print("clim ............ {}".format(CLIM))
    print("convergence ..... {}\n".format(cvalue))
    if cvalue < CLIM:
        print(NOT_CONVERGED)
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root_name, cvalue))
    else:
        print(CONVERGED)
    print("")

    return cvalue


def do_py_plot_output(wd, root_name):
    """
    Execute the standard py_plot_output routine located in py_progs
    """

    path = which("py_plot_output.py")
    if not path:
        print("py_plot_output.py not in $PATH")
        return

    commands = "cd {}; Setup_Py_Dir; py_plot_output.py {} all; rm data"\
        .format(wd, root_name)
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
        print("spec_plot.py not in $PATH")
        return

    commands = "cd {}; spec_plot.py {} all".format(wd, root_name)
    if TDE_PLOT:
        commands += " -blag"
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


def find_number_of_physical_cores():
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


def main():
    """
    Main control function
    """

    # Determine which routines to run for each simulation
    get_run_mode()

    all_par_file_paths = py_util.find_pf(ignore_out_pf=True)
    mpi = check_for_mpi()

    n_cores = 1
    if RUN_SIMS and mpi:
        # We will suffer for macOS
        if system() == "Darwin":
            n_cores = cpu_count()
        else:
            n_cores = find_number_of_physical_cores()
            if n_cores == 0:
                n_cores = cpu_count()
        print("The following simulations will be run using {} cores:\n"
              .format(n_cores))
    elif RUN_SIMS and not mpi:
        print("The follow simulations will be run in single threaded mode:\n")
    else:
        print("Processing the following runs:\n")

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
    n_sims = len(par_file_paths)

    # If the not_converged file exists, delete the contents :-)
    if SHOW_CONVERGENCE:
        open("not_converged.txt", "w").close()

    # Write everything out to file
    outfname = "py_run_{}{}{}.txt"\
        .format(DATE.year, DATE.month, DATE.day, DATE.hour, DATE.minute)
    with open(outfname, "w") as f:
        for i, path in enumerate(par_file_paths):
            root_name, pf_relative_path = py_util.get_root_name_and_path(path)
            print("--------------------------\n")
            print("Simulation {}/{}\n".format(i+1, n_sims))

            # wark wark wark
            if RUN_SIMS:
                print("Running the simulation: {}\n".format(root_name))
                python_output, python_err = py_run(pf_relative_path, root_name,
                                                   mpi, n_cores)

                f.write("--------\nPython Output\n--------\n")
                f.write("Working dir ........ {}\n".format(pf_relative_path))
                f.write("Root name .......... {}\n".format(root_name))
                f.writelines(python_output)
                if python_err:
                    f.write(python_err)

            # Assumes plot_convergence.py is in $PATH
            if SHOW_CONVERGENCE:
                print("Checking the convergence of the simulation:\n")
                convergence_output = check_convergence(pf_relative_path,
                                                       root_name, SHOW_OUTPUT)
                cvalue = display_convergence(convergence_output.split(),
                                             pf_relative_path, root_name)
                f.write("--------\nConvergence\n--------\n")
                f.writelines(convergence_output)
                if cvalue < CLIM:
                    f.write("{}\n".format(NOT_CONVERGED))
                    f.write("cvalue {} < clim {}\n"
                            .format(cvalue, CLIM))
                else:
                    f.write("{}\n".format(CONVERGED))
                    f.write("cvalue {} >= clim {}\n"
                            .format(cvalue, CLIM))
                f.write("\n--------------------------")

            # Assumes plotting scripts are in $PATH
            if CREATE_PLOTS:
                print("Creating plots for the simulation\n")
                do_py_plot_output(pf_relative_path, root_name)
                do_spec_plot(pf_relative_path, root_name)

        print("--------------------------")

    return


if __name__ == "__main__":
    main()
