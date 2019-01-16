#!/usr/bin/env python3


import sys
import argparse
from shutil import which
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, CalledProcessError


CLIM = 0.85
VERSION = "py"
SHOW_OUTPUT = False
RUN_SIMS = False
SHOW_CONVERGENCE = False
non_converged = False
CREATE_PLOTS = False
TDE_PLOT = False


def get_run_mode():
    """
    Read command line flags to determine the mode of operation
    """

    p = argparse.ArgumentParser(description="Enable or disable features")
    p.add_argument("-D", "-d", action="store_true", help="Run with -R -C -P")
    p.add_argument("-R", "-r", action="store_true", help="Run simulations")
    p.add_argument("-C", "-c", action="store_true", help="Check the convergence of runs - calls convergence.py")
    p.add_argument("-C_LIM", type=float, action="store", help="The convergence limit: c_value < 1")
    p.add_argument("-P", "-p", action="store_true", help="Run plotting scripts")
    p.add_argument("-O", "-o", action="store_true", help="Verbose outputting")
    p.add_argument("-TP", "-tp", action="store_true", help="Enable TDE plotting")
    p.add_argument("-PY_VER", "-py_ver", type=str, action="store", help="Name of the Python executable")
    args = p.parse_args()

    global SHOW_OUTPUT
    global RUN_SIMS
    global SHOW_CONVERGENCE
    global CREATE_PLOTS
    global CLIM
    global TDE_PLOT
    global VERSION

    do_something = False

    if args.D:
        RUN_SIMS = True
        SHOW_CONVERGENCE = True
        CREATE_PLOTS = True
        do_something = True
    if args.O:
        SHOW_OUTPUT = True
        do_something = True
    if args.R:
        RUN_SIMS = True
        do_something = True
    if args.C:
        SHOW_CONVERGENCE = True
        do_something = True
    if args.P:
        CREATE_PLOTS = True
        do_something = True
    if args.TP:
        CREATE_PLOTS = True
        TDE_PLOT = True
        do_something = True
    if args.PY_VER:
        VERSION = args.PY_VER
    if args.C_LIM:  # Bad variable names here :-) C_LIM is from the arg list
        if 0 < args.c_value < 1:
            CLIM = args.C_LIM
        else:
            print("Invalid value of c_value {}".format(args.c_value))
            sys.exit(1)

    print("\nRun Simulations ............ {}".format(RUN_SIMS))
    print("Show Verbose Output ........ {}".format(SHOW_OUTPUT))
    print("Show Convergence ........... {}".format(SHOW_CONVERGENCE))
    print("Create Plots ............... {}".format(CREATE_PLOTS))
    print("Convergence Limit .......... {}".format(CLIM))
    print("Python version ............. {}".format(VERSION))
    if args.TP:
        print("Plot TDE ............... {}".format(TDE_PLOT))
    print("")

    if do_something is False:
        print("\nThere is nothing to do!\n")
        print("--------------------------")
        sys.exit(0)

    return


def py_run(wd, root_name):
    """
    Do a standard Python run using a single process
    """

    pf = root_name + ".pf"  # Don't actually need the .pf at the end
    command = "cd {}; Setup_Py_Dir; {} {}".format(wd, VERSION, pf)

    print("Working dir ........ {}".format(wd))
    print("Root name .......... {}".format(root_name))
    print("\n{}\n".format(command))

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    if SHOW_OUTPUT:
        for stdout_line in iter(cmd.stdout.readline, ""):
            if not stdout_line:
                break
            print(stdout_line.decode("utf-8").replace("\n", ""))
        print("")

    #
    # When the stdout or stderr buffer becomes too large (>4KB), cmd.wait() deadlocks.
    # We can get around this by using cmd.communicate instead, I think.
    #

    pystdout, pystderr = cmd.communicate()
    returncode = cmd.returncode
    if returncode:
        raise CalledProcessError(returncode, command)
    output = pystdout.decode("utf-8")
    err = pystderr.decode("utf-8")

    if err:
        print("Captured from stderr:")
        print(err)

    with open("{}/py_{}.out".format(wd, root_name), "w") as f:
        f.writelines(output)
    if err:
        with open("{}/py_err_{}.out".format(wd, root_name), "w") as f:
            f.writelines(err)

    return output, err


def mpi_py_run(wd, root_name, ncores):
    """
    Do a standard Python multi-processor run
    """

    pf = root_name + ".pf"  # Don't actually need the .pf at the end
    command = "cd {}; Setup_Py_Dir; mpirun -n {} {} {}".format(wd, ncores, VERSION, pf)

    print("Working dir ........ {}".format(wd))
    print("Root name .......... {}".format(root_name))
    print("\n{}\n".format(command))

    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)

    if SHOW_OUTPUT:
        for stdout_line in iter(cmd.stdout.readline, ""):
            if not stdout_line:
                break
            print(stdout_line.decode("utf-8").replace("\n", ""))
        print("")

    #
    # When the stdout or stderr buffer becomes too large (>4KB), cmd.wait() deadlocks.
    # We can get around this by using cmd.communicate instead, I think.
    #

    pystdout, pystderr = cmd.communicate()
    returncode = cmd.returncode
    if returncode:
        raise CalledProcessError(returncode, command)
    output = pystdout.decode("utf-8")
    err = pystderr.decode("utf-8")

    if err:
        print("Captured from stderr:")
        print(err)

    with open("{}/py_{}.out".format(wd, root_name), "w") as f:
        f.writelines(output)
    if err:
        with open("{}/py_err_{}.out".format(wd, root_name), "w") as f:
            f.writelines(err)

    return output, err


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

    conv = output.split()
    final_cycle = 0
    for i, word in enumerate(conv):
        if word == "!!Check_converging:":
            final_cycle = i
    if final_cycle == 0:
        return output, 0

    final_cycle_conv = conv[final_cycle:]
    cvalue = float(final_cycle_conv[2].replace("(", "").replace(")", ""))

    print("dir ............. {}".format(wd))
    print("root name ....... {}".format(root_name))
    print("clim ............ {}".format(CLIM))
    print("convergence ..... {}\n".format(cvalue))
    if cvalue < CLIM:
        print("NOT CONVERGED")
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root_name, cvalue))
    else:
        print("CONVERGED")
    print("")

    return output, cvalue


def do_py_plot_output(wd, root_name):
    """
    Execute the standard py_plot_output routine located in py_progs
    """

    commands = "cd {}; py_plot_output.py {} all".format(wd, root_name)
    print(commands)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")

    if SHOW_OUTPUT:
        print(output)

    return


def do_spec_plot(wd, root_name):
    """
    Execute my standard spec_plot script to plot spectra
    """

    commands = "cd {}; spec_plot.py {} all".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if SHOW_OUTPUT:
        print(output)
    if err:
        print("Captured from stderr:")
        print(err)

    return


def TDE_plot(wd, root_name):
    """
    Execute the standard Blag spec_plotTDE plotting script
    """

    commands = "cd {}; spec_plotTDE.py {} all -blag -wmin 1000 -wmax 2500".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if SHOW_OUTPUT:
        print(output)
    if err:
        print("Captured from stderr:")
        print(err)

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
    Created due to multiprocessing.cpu_count() returning physical and logical cores. Pathetic.
    """

    ncores_cmd = "lscpu | grep 'Core(s) per socket'"
    nsockets_cmd = "lscpu | grep 'Socket(s):'"
    ncores, stderr = Popen(ncores_cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    if ncores == "":
        return 0
    ncores = int(ncores.decode("utf-8").replace("\n", "").split()[-1])
    nsockets, stderr = Popen(nsockets_cmd, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    if nsockets == "":
        return 0
    nsockets = int(nsockets.decode("utf-8").replace("\n", "").split()[-1])

    return ncores * nsockets


def find_pf():
    """
    Find parameter files recursively
    """

    command = "find . -type f -name '*.pf'"
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    pfs = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")
    if err:
        print(err)

    for i, dir in enumerate(pfs):
        if dir.find(".out.pf") != -1:
            del(pfs[i])

    if len(pfs) == 0:
        print("No Python parameter files found")
        print("\n--------------------------")
        sys.exit(0)

    return pfs


def get_root_name_and_path(pf_path):
    """
    Split the path name into a directory and root name of the Python simulation
    """

    dot = 0
    slash = 0
    for l in range(len(pf_path)):
        if pf_path[l] == "/":
            slash = l + 1
        if pf_path[l] == ".":
            dot = l
    root_name = pf_path[slash:dot]
    sim = pf_path[:slash]

    return root_name, sim


def main():
    """
    Main control function
    """

    print("--------------------------\n")

    all_par_file_paths = find_pf()
    mpi = check_for_mpi()

    n_cores = 1
    if mpi:
        n_cores = find_number_of_physical_cores()  # I want the physical CPUs
        if n_cores == 0:
            n_cores = cpu_count()       # Not the physical CPU and SMT cores
        print("The following simulations will be run using {} cores:\n".format(n_cores))
    else:
        print("The follow simulations will be run in single threaded mode:\n")

    # Remove any py_wind parameter files which might be lurking about
    par_file_paths = []
    for i, path in enumerate(all_par_file_paths):
        if path.find("py_wind.pf") == -1:
            par_file_paths.append(path)
            print("- {}".format(path))

    # Determine which routines to run for each simulation
    get_run_mode()

    # If the not_converged file exists, delete the contents :-)
    if SHOW_CONVERGENCE:
        open("not_converged.txt", "w").close()

    # Write everything out to file
    with open("output.txt", "w") as f:
        for path in par_file_paths:
            root_name, pf_relative_path = get_root_name_and_path(path)
            print("--------------------------\n")

            # wark wark wark
            if RUN_SIMS:
                print("Running the simulation: {}\n".format(root_name))
                if mpi:
                    python_output, python_err = mpi_py_run(pf_relative_path, root_name, n_cores)
                else:
                    python_output, python_err = py_run(pf_relative_path, root_name)

                f.write("--------\nPython Output\n--------\n")
                f.write("Working dir ........ {}\n".format(pf_relative_path))
                f.write("Root name .......... {}\n".format(root_name))
                f.writelines(python_output)
                if python_err:
                    f.write(python_err)

            # Assumes plot_convergence.py is in $PATH
            if SHOW_CONVERGENCE:
                print("Checking the convergence of the simulation:\n")
                convergence_output, cvalue = check_convergence(pf_relative_path, root_name, SHOW_OUTPUT)
                f.write("--------\nConvergence\n--------\n")
                f.writelines(convergence_output)
                if cvalue < CLIM:
                    f.write("NOT CONVERGED, cvalue {} < clim {}\n".format(cvalue, CLIM))
                else:
                    f.write("CONVERGED, cvalue {} >= clim {}\n".format(cvalue, CLIM))
                f.write("\n--------------------------")

            # Assumes plotting scripts are in $PATH
            if CREATE_PLOTS:
                print("Creating plots for the simulation\n")
                do_py_plot_output(pf_relative_path, root_name)
                do_spec_plot(pf_relative_path, root_name)
                if TDE_PLOT:
                    TDE_plot(pf_relative_path, root_name)

        print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
