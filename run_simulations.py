#!/usr/bin/env python3

from subprocess import Popen, PIPE
from multiprocessing import cpu_count


show_output = False


def standard_python(wd, root_name, ncores, pyver):
    """
    Do a standard Python multi-processor run
    """

    pf = root_name + ".pf"
    command = "cd {}; Setup_Py_Dir; mpirun -n {} {} {}".format(wd, ncores, pyver, pf)
    print("Working dir ........ {}".format(wd))
    print("Root name .......... {}".format(root_name))
    print("\n{}\n".format(command))
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    pystdout, pystderr = cmd.communicate()
    output = pystdout.decode("utf-8")
    err = pystderr.decode("utf-8")

    if show_output:
        print("Output from Python:")
        print(output)

    if err:
        print("Errors returned from stderr:")
        print(err)

    with open("{}/py_{}.out".format(wd, root_name), "w") as f:
        f.writelines(output)
    if err:
        with open("{}/py_err_{}.out".format(wd, root_name), "w") as f:
            f.writelines(err)

    return output


def convergence(wd, root_name):
    """
    Print to screen and plot the convergence of the simulation
    """

    commands = "cd {}; convergence.py {}".format(wd, root_name)
    cmd = Popen(commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    print(output)
    if err:
        print(err)

    return output


def TDE_plot(wd):
    """
    Execute the standard Blag spec TDE plotting script
    """

    commands = "cd {}; spec_plotTDE.py tde all -dist 1.079987153448e+27 -wmin 1000 -wmax 2500".format(wd)
    cmd = Popen (commands, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    print(output)
    if err:
        print(err)

    return


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
        if dir.find(".out") != -1:
            del(pfs[i])

    return pfs


def main(py):
    """
    Main control function
    """

    show_convergence = True
    create_plots = True

    ncores = cpu_count()
    sims_to_run = find_pf()

    print("The following simulations will be run using {} cores:\n".format(ncores))
    for sim in sims_to_run:
        print("- {}".format(sim))

    with open("output.txt", "w") as f:
        for sim in sims_to_run:
            # Find root name and clean up sim path
            dot = 0
            slash = 0
            for l in range(len(sim)):
                if sim[l] == "/":
                    slash = l + 1
                if sim[l] == ".":
                    dot = l
            root_name = sim[slash:dot]
            sim = sim[:slash]

            # Do some wark
            print("\n--------------------------\n")
            f.write("--------\nPython Output\n--------\n")
            f.write("Working dir ........ {}\n".format(sim))
            f.write("Root name .......... {}\n".format(root_name))
            python_output = standard_python(sim, root_name, ncores, py)
            f.writelines(python_output)

            # Assumes plot_convergence.py is in $PATH
            if show_convergence:
                f.write("--------\nConvergence\n--------\n")
                convergence_output = convergence(sim, root_name)
                f.writelines(convergence_output)
            # Assumes plotting scripts are in $PATH
            if create_plots:
                TDE_plot(sim)

            f.write("\n--------------------------")
            print("\n--------------------------")

    return


if __name__ == "__main__":
    main("py")
