#!/usr/bin/env python3

import sys
import argparse
import numpy as np
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt

showplot = False  # If True, the plot will be drawn on screen
filetype = "png"  # The file type of the output


def plot_convergence(rootname, converged, converging):
    """
    Plot the fraction of converged and converging cells as a function of cycle.

    Parameters
    ----------
    converged: 2 x ncycles array of floats
        The fraction of converged cells for each Python ionisation cycle.
        Col 0: cycle
        Col 1: converged fraction
    converging: 2 x ncycles array of floats
        The fraction of converging cells for each Python ionisation cycle
        Col 0: cycle
        Col 1: converging fraction

    Returns
    -------
    None
    """

    ncycles = converged.shape[0]
    max_convergence = np.max(converged[:, 1])

    plt.figure(figsize=(8, 8))
    plt.plot(converged[:, 0], converged[:, 1], label="Fraction of converged cells")
    plt.plot(converging[:, 0], converging[:, 1], label="Fraction of converging cells")
    plt.ylim(0, 1)
    plt.xlim(0, ncycles - 1)
    plt.legend()
    plt.xlabel("Ionisation Cycle")
    plt.ylabel("Total Fraction")
    plt.title("Max convergence reached = {}".format(max_convergence))
    plt.savefig("{}_convergence.{}".format(rootname, filetype))
    plt.close()

    return


def grep_convergence(filename):
    """
    Grep the convergence fractions from the diag file given by filename.

    Parameters
    ----------
    filename: str
        The file path to the diag file to be grepped for the model convergence.

    Returns
    -------
    converged: 2 x ncycles array of floats
        The fraction of converged cells for each Python ionisation cycle.
        Col 0: cycle
        Col 1: converged fraction
    converging: 2 x ncycles array of floats
        The fraction of converging cells for each Python ionisation cycle
        Col 0: cycle
        Col 1: converging fraction
    """

    # Create a grep command and pass to subprocess to get stdout and stderr
    grep = "grep \!\!Check_converging {}".format(filename)
    stdout, stderr = Popen(grep, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    out = stdout.decode("utf-8").split()
    if len(out) == 0:
        print("No output from grep")
        sys.exit(0)
    # print(stdout.decode("utf-8"))

    # Loop over the grep output and separate converged and converging fractions into lists
    cycle = 0
    converged = []
    converging = []
    for i in range(len(out)):
        if out[i] == "converged":
            n_converged = out[i-1].replace("(", "").replace(")", "")
            converged.append([cycle, float(n_converged)])
        elif out[i] == "converging":
            n_converging = out[i-1].replace("(", "").replace(")", "")
            converging.append([cycle, float(n_converging)])
            cycle += 1

    ncycles = len(converged)
    for i in range(ncycles):
        print("Cycle {:2d}/{:2d}: {:3.0f}% cells converged and {:3.0f}% cells still converging".format(
                i + 1, ncycles, converged[i][1] * 100, converging[i][1] * 100))

    # Return the lists as arrays to make life easier later
    return np.array(converged), np.array(converging)


def parse_rootname():
    """
    Get the rootname from the command line argument.

    Parameters
    ----------
    None

    Returns
    -------
    rootname: str
from sys import argv
        The rootname of the Python simulation.
    """

    p = argparse.ArgumentParser(description="Please provide the rootname of the Python simulation")
    p.add_argument("rootname", type=str, help="The name of the Python pf")
    p.add_argument("-filetype", type=str, nargs="?", action="store",
                        help="The filetype for the output plot")
    p.add_argument("-show", action="store_true", help="Show the plot to screen")
    args = p.parse_args()
    if args.filetype:
        global filetype
        filetype = args.filetype
    if args.show:
        global showplot
        showplot = args.show

    return args.rootname


def main():
    """
    The main steering routine of the script.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    rootname = parse_rootname()
    print("Finding convergence for Python simulation {}.pf".format(rootname))
    filename = "diag_{}/{}_0.diag".format(rootname, rootname)
    converged, converging = grep_convergence(filename)
    if showplot:
        plot_convergence(rootname, converged, converging)

    return

if __name__ == "__main__":
    main()
