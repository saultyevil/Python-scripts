#!/usr/bin/env python3

from subprocess import Popen, PIPE

show_output = False
convergence = True
plots = True


def standard_python(wd, pf):
    """
    Initialise a standard Python multi-processor run
    """

    cd = "cd {}".format(wd)
    cmd = Popen(cd, stdout=PIPE, stderr=PIPE, shell=True)
    print("Currently in: {}".format(wd))
    cmd = Popen("Setup_Py_Dir", stdout=PIPE, stderr=PIPE, shell=True)
    cmd = Popen("mpirun py {}".format(pf), stdout=PIPE, stderr=PIPE, shell=True)
    pystdout, pystderr = cmd.communicate()
    output = pystdout.decode("utf-8")
    err = pystderr.decode("utf-8")
    if show_output:
        print("Output from Python:")
        print(output)
        if err:
            print("Errors returned from stderr:")
            print(err)
    with open("py_{}.out", "w") as f:
        f.writelines(output)
    if err:
        with open("py_err_{}.out", "w") as f:
            f.writelines(err)
    if convergence:
        cmd = Popen("convergence.py tde", stdout=PIPE, stderr=PIPE, shell=True)
    if plots:
        cmd = Popen("spec_plotTDE.py tde all -dist 1.079987153448e+27")

    return output


if __name__ == "__main__":
    with open("simulations", "r") as f:
        sims = f.readlines()
    f = open("script_output.txt", "w")
    for i in sims:
        output = standard_python(i, "tde.pf")
        f.write(i)
        f.writelines(output)
    f.close()