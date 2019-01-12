#!/usr/bin/env python3

import multiprocessing
from subprocess import Popen, PIPE


show_output = False


def standard_python(wd, root_name, ncores, pyver):
    """
    Do a standard Python multi-processor run
    """

    pf = root_name + ".pf"
    command = "cd {}; Setup_Py_Dir; mpirun -n {} {} {}".format(wd, ncores, pyver, pf)
    print(command)
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
    pfs = stdout.decode("utf-8")
    err = stderr.decode("utf-8")
    if err:
        print(err)

    return pfs


def main(py):
    """
    Main control function

    Parameters
    ----------
    py: str
        The command to call py from the command line
    """

    show_convergence = True
    create_plots = True

    sims = find_pf().split()
    ncores = multiprocessing.cpu_count()

    print("The following simulations will be run,")
    for sim in sims:
        sim = sim.strip("\r\n")
        print("- {}".format(sim))

    with open("output.txt", "w") as f:
        for sim in sims:
            # Find the root name of the simulation
            dot = 0
            slash = 0
            for i in range(len(sim)):
                if sim[i] == "/":
                    slash = i + 1
                if sim[i] == ".":
                    dot = i
            root_name = sim[slash:dot]
            if root_name.find(".out") != -1:
                root_name = root_name[:len(".out") - 1]

            # Wark
            sim = sim[:slash]
            print(sim)
            print("\n--------------------------\n\n{}".format(sim))
            python_output = standard_python(sim, root_name, ncores, py)
            f.write("{}\n".format(sim))
            f.write("--------\nPython Output\n--------\n")
            f.writelines(python_output)
            if show_convergence:
                convergence_output = convergence(sim, root_name)
                f.write("--------\nConvergence\n--------\n")
                f.writelines(convergence_output)
            if create_plots:
                TDE_plot(sim)
            f.write("\n--------------------------\n")
            print("\n--------------------------\n")

    return


if __name__ == "__main__":
    main("py")
