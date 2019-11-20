#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Check the progress and convergence of a Python simulation. This script parses
out text from the diag file of the Python simulation.

NOTE: this script will only work if grep is installed and callable from the
      command line - who knows if this works with Powershell or cmd.

Usage
    py_check_run.py [-h] [-filetype [FILETYPE]] [-show] root

    Please provide the rootname of the Python simulation

    positional arguments:
      rootname              The name of the Python pf

    optional arguments:
      -h, --help            show this help message and exit
      -filetype [FILETYPE]  The filetype for the output plot
      -show                 Show the plot to screen
"""

from os import access, R_OK
import sys
import argparse
import numpy as np
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from typing import List, Tuple

SHOW_PLOT = False  # If True, the plot will be drawn on screen
FILE_TYPE = "png"  # The file type of the output


def plot_convergence(root: str, converged: np.ndarray, converging: np.ndarray) -> None:
    """
    Plot the fraction of converged and converging cells as a function of cycle.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    converged: 2 x ncycles array of floats
        The fraction of converged cells for each Python ionisation cycle.
        Col 0: cycle
        Col 1: converged fraction
    converging: 2 x ncycles array of floats
        The fraction of converging cells for each Python ionisation cycle
        Col 0: cycle
        Col 1: converging fraction
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
    plt.savefig("{}_convergence.{}".format(root, FILE_TYPE))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def get_convergence(filename: str) -> Tuple[np.ndarray, np.ndarray]:
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
    grep = r"grep '\!\!Check_convergence:' {}".format(filename)
    stdout, stderr = Popen(grep, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    conv_out = stdout.decode("utf-8").split()
    if len(conv_out) == 0:
        print("No convergence output from grep")
        sys.exit(0)
    grep = r"grep '\!\!Python: Beginning cycle' {}".format(filename)
    stdout, stderr = Popen(grep, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    cycle_out = stdout.decode("utf-8").split()
    if len(cycle_out) == 0:
        print("No cycle output from grep")
        sys.exit(0)
    ncycles = int(cycle_out[-4])

    # Loop over the grep output and separate converged and converging fractions
    # into lists
    cycle = 0
    converged = []
    converging = []
    for i in range(len(conv_out)):
        if conv_out[i] == "converged":
            n_converged = conv_out[i - 1].replace("(", "").replace(")", "")
            converged.append([cycle, float(n_converged)])
        elif conv_out[i] == "converging":
            n_converging = conv_out[i - 1].replace("(", "").replace(")", "")
            converging.append([cycle, float(n_converging)])
            cycle += 1

    nlines = len(converged)
    for i in range(nlines):
        print("Cycle {:2d}/{:2d}: {:3.0f}% cells converged and {:3.0f}% cells still converging".format(
                i + 1, ncycles, converged[i][1] * 100, converging[i][1] * 100))

    # Return the lists as arrays to make life easier later
    return np.array(converged), np.array(converging)


def parse_rootname() -> str:
    """
    Prase the root name and various run time options from the command line.

    Returns
    -------
    root: str
        The root name of the Python simulation.
    """

    p = argparse.ArgumentParser(description="Please provide the root name of the Python simulation")
    p.add_argument("root", type=str, help="The root name of the Python simulation")
    p.add_argument("-show", action="store_true", help="Create a plot of convergence v cycle")
    p.add_argument("-filetype", type=str, nargs="?", action="store", help="The file type for the output plot")
    args = p.parse_args()

    if args.filetype:
        global FILE_TYPE
        FILE_TYPE = args.filetype
    if args.show:
        global SHOW_PLOT
        SHOW_PLOT = args.show

    return args.root


def main():
    """
    The main steering routine of the script.
    """

    root = parse_rootname()
    msg = "Convergence for Python simulation: {}".format(root)
    print("{}\n{}".format(msg, "-" * len(msg)))

    filename = "diag_{}/{}_0.diag".format(root, root)
    if not access(filename, R_OK):
        filename = "{}_0.diag".format(root)
        if not access(filename, R_OK):
            print("Could not find diag file for {} in calling directory\n".format(root))
            sys.exit(1)

    converged, converging = get_convergence(filename)
    plot_convergence(root, converged, converging)

    return


if __name__ == "__main__":
    main()
