#!/usr/bin/env python3

from multiprocessing import cpu_count
from subprocess import Popen, PIPE, CalledProcessError


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

    run_sims = False
    show_convergence = False
    non_converged = False
    create_plots = True

    global show_output
    show_output = True

    clim = 0.85
    ncores = cpu_count()
    sims_to_run = find_pf()

    print("The following simulations will be run using {} cores:\n".format(ncores))
    sims = []
    for i, sim in enumerate(sims_to_run):
        if sim.find("py_wind.pf") == -1:
            sims.append(sim)
            print("- {}".format(sim))

    open("not_converged.txt", "w").close()

    with open("output.txt", "w") as f:
        for sim in sims:
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

            print("\n--------------------------\n")

            # Do some wark
            if run_sims:
                f.write("--------\nPython Output\n--------\n")
                f.write("Working dir ........ {}\n".format(sim))
                f.write("Root name .......... {}\n".format(root_name))
                python_output = standard_python(sim, root_name, ncores, py)
                f.writelines(python_output)

            # Assumes plot_convergence.py is in $PATH
            if show_convergence:
                f.write("--------\nConvergence\n--------\n")
                convergence_output = convergence(sim, root_name, show_output)
                f.writelines(convergence_output)
            if non_converged:
                cvalue = list_non_converged(sim, root_name, clim)
                if cvalue < clim:
                    f.write("NOT CONVERGED, cvalue {} < clim {}\n".format(cvalue, clim))
                else:
                    f.write("CONVERGED, cvalue {} >= clim {}\n".format(cvalue, clim))

            # Assumes plotting scripts are in $PATH
            if create_plots:
                standard_py_progs(sim, root_name)
                # standard_spec_plot(sim)
                TDE_plot(sim)

        print("\n--------------------------")

    return


if __name__ == "__main__":
    main("py")
