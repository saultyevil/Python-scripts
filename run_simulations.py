#!/usr/bin/env python3


import sys
import argparse
from multiprocessing import cpu_count
from subprocess import Popen, PIPE, CalledProcessError


clim = 0.85
show_output = False
run_sims = False
show_convergence = False
non_converged = False
create_plots = False


def mpi_py_run(wd, root_name, ncores, pyver):
    """
    Do a standard Python multi-processor run
    """

    pf = root_name + ".pf"
    command = "cd {}; Setup_Py_Dir; mpirun -n {} {} {}".format(wd, ncores, pyver, pf)
    print("Working dir ........ {}".format(wd))
    print("Root name .......... {}".format(root_name))
    print("\n{}\n".format(command))
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    if show_output:
        for stdout_line in iter(cmd.stdout.readline, ""):
            if not stdout_line:
                break
            print(stdout_line.decode("utf-8").replace("\n", ""))
    returncode = cmd.wait()
    if returncode:
        raise CalledProcessError(returncode, command)
    pystdout, pystderr = cmd.communicate()
    output = pystdout.decode("utf-8")
    err = pystderr.decode("utf-8")
    cmd.stdout.close()

    if err:
        print("Errors returned from stderr:")
        print(err)

    with open("{}/py_{}.out".format(wd, root_name), "w") as f:
        f.writelines(output)
    if err:
        with open("{}/py_err_{}.out".format(wd, root_name), "w") as f:
            f.writelines(err)

    return output


def convergence(wd, root_name, show_output):
    """
    Print to screen and plot the convergence of the simulation
    """

    commands = "cd {}; convergence.py {}".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)

    if err:
        print(err)

    return output


def list_non_converged(wd, root_name, clim):
    """
    Create a list of simulations which haven't reached a defined convergence limit
    """

    conv = convergence(wd, root_name, False).split()
    final_cycle = 0
    for i, word in enumerate(conv):
        if word == "!!Check_converging:":
            final_cycle = i
    if final_cycle == 0:
        return

    final_conv = conv[final_cycle:]
    cvalue = float(final_conv[2].replace("(", "").replace(")", ""))
    print("dir ............. {}".format(wd))
    print("root name ....... {}".format(root_name))
    print("clim ............ {}".format(clim))
    print("convergence ..... {}".format(cvalue))
    if cvalue < clim:
        print("NOT CONVERGED")
        with open("not_converged.txt", "a") as f:
            f.write("{}\t{}.pf\t{}\n".format(wd, root_name, cvalue))
    else:
        print("CONVERGED")

    return cvalue


def standard_py_progs(wd, root_name):
    """
    Execute the standard py_plot_output routine located in py_progs
    """

    commands = "cd {}; py_plot_output.py {} all".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)
    if err:
        print(err)

    return


def standard_spec_plot(wd):
    """
    Execute a slightly different TDE plot script
    """

    commands = "cd {}; spec_plot.py tde all -wmin 200 -wmax 2600".format(wd)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)
    if err:
        print(err)

    return


def TDE_plot(wd):
    """
    Execute the standard Blag spec TDE plotting script
    """

    # commands = "cd {}; spec_plotTDE.py tde all -dist 1.079987153448e+27 -wmin 1000 -wmax 2500".format(wd)
    commands = "cd {}; spec_plotTDE.py tde all -dist 1.079987153448e+27 -wmin 200 -wmax 2600".format(wd)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if show_output:
        print(output)
    if err:
        print(err)

    return


def clean_up_names(pf_path):
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


def find_physical_cpus():
    """
    Created due to multiprocessing.cpu_count() returning physical and logical threads. Pathetic.
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


def get_run_mode():
    """
    Read command line flags to determine the mode of operation
    """

    p = argparse.ArgumentParser(description="Enable or disable features")
    p.add_argument("-D", "-d", action="store_true", help="Run with -R -C -P")
    p.add_argument("-R", "-r", action="store_true", help="Run simulations")
    p.add_argument("-C", "-c", action="store_true", help="Check the convergence of runs - calls convergence.py")
    p.add_argument("-c_value", type=float, action="store", help="The convergence limit: c_value < 1")
    p.add_argument("-P", "-p", action="store_true", help="Run plotting scripts")
    p.add_argument("-V", "-v", action="store_true", help="Verbose outputting")
    args = p.parse_args()

    global show_output
    global run_sims
    global show_convergence
    global create_plots
    global clim

    do_something = False

    if args.D:
        run_sims = True
        show_convergence = True
        create_plots = True
        do_something = True
    if args.V:
        show_output = True
        do_something = True
    if args.R:
        run_sims = True
        do_something = True
    if args.C:
        show_convergence = True
        do_something = True
    if args.P:
        create_plots = True
        do_something = True
    if args.c_value:
        if 0 < args.c_value < 1:
            clim = args.c_value
        else:
            print("Invalid value of c_value {}".format(args.c_value))
            sys.exit(1)

    print("\nRun Simulations ............ {}".format(run_sims))
    print("Show Verbose Output ........ {}".format(show_output))
    print("Show Convergence ........... {}".format(show_convergence))
    print("Create Plots ............... {}".format(create_plots))
    print("Convergence Limit .......... {}\n".format(clim))

    if do_something is False:
        print("\nThere is nothing to do!")
        print("\n--------------------------")
        sys.exit(0)

    return


def main(py):
    """
    Main control function
    """

    print("--------------------------\n")

    sim_paths = find_pf()
    n_cores = find_physical_cpus()  # I want the physical CPUs
    if n_cores == 0:
        n_cores = cpu_count()       # Not the physical CPU and SMT cores

    print("The following simulations will be run using {} cores:\n".format(n_cores))

    # Remove any py_wind parameter files which might be lurking about
    py_sims = []
    for i, sim in enumerate(sim_paths):
        if sim.find("py_wind.pf") == -1:
            py_sims.append(sim)
            print("- {}".format(sim))

    # Determine which routines to run for each simulation
    get_run_mode()

    # If the not_converged file exists, delete the contents :-)
    if show_convergence:
        open("not_converged.txt", "w").close()

    # Write everything out to file
    with open("output.txt", "w") as f:
        for sim in py_sims:
            root_name, pf_path = clean_up_names(sim)
            print("--------------------------\n")

            # wark wark wark
            if run_sims:
                print("Running simulation")
                f.write("--------\nPython Output\n--------\n")
                f.write("Working dir ........ {}\n".format(pf_path))
                f.write("Root name .......... {}\n".format(root_name))
                python_output = mpi_py_run(pf_path, root_name, n_cores, py)
                f.writelines(python_output)

            # Assumes plot_convergence.py is in $PATH
            if show_convergence:
                f.write("--------\nConvergence\n--------\n")
                convergence_output = convergence(pf_path, root_name, show_output)
                f.writelines(convergence_output)
                cvalue = list_non_converged(pf_path, root_name, clim)
                if cvalue < clim:
                    f.write("NOT CONVERGED, cvalue {} < clim {}\n".format(cvalue, clim))
                else:
                    f.write("CONVERGED, cvalue {} >= clim {}\n".format(cvalue, clim))

            # Assumes plotting scripts are in $PATH
            if create_plots:
                print("Creating plots")
                standard_py_progs(pf_path, root_name)
                # standard_spec_plot(sim)
                TDE_plot(pf_path)

        print("\n--------------------------")

    return


if __name__ == "__main__":
    main("py")
