#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Create and run a grid of Python simulations. Note that this script has to be
edited before being called as I didn't want to add one billion command line
arguments/switches to the script like with py_run.

Usage

    [python] py_run_grid.py [--dry-run]

    --dry-run: only create the grid of parameter files
"""


import os
import py_run
import shutil
import numpy as np
import py_run_util
from sys import argv
from typing import List


def create_grid(pf: str, parameter: str, grid: List[str]) -> List[str]:
    """
    The purpose of this function is to use a base parameter file and to create
    directories containing parameter files with a different value to the given
    parameter. These values are given in grid.

    This function will only work with one parameter and grids of multiple
    parameters cannot be created.

    A backup of the original base pf is made just in case :-).

    Parameters
    ----------
    pf              str
                    Path to the base Python pf
    parameter       str
                    The name of the parameter for which a grid will be created
    grid            List[str]
                    The values of the parameter to make a grid with

    Returns
    -------
    pfs             List[str]
                    A list containing the file paths to the newly created pfs.
    """

    if pf.find(".pf") == -1:
        pf += ".pf"

    sl = 0
    for i, letter in enumerate(pf):  # This iterates over the pf file path
        if letter == "/":
            sl = i  # This will find the index of the final /
    pl = pf.find(".pf")
    root = pf[sl:pl]  # Now we can extract just the root name from the file path
    shutil.copyfile(pf, pf + ".bak")

    with open(pf, "r") as f:
        pf_lines = f.readlines()

    base = []
    for letter in pf_lines:
        if letter[0] == "#" or letter[0] == "\n" or letter[0] == "\r":
            continue
        letter = letter.replace("\n", "").split()
        base.append([letter[0], letter[1]])

    # This bit creates the grid of pfs
    pfs = []
    npfs = len(grid)
    for i in range(npfs):
        new_pf = base.copy()
        for letter in range(len(new_pf)):
            if new_pf[letter][0] == parameter:
                new_pf[letter][1] = grid[i]
        try:
            os.mkdir(grid[i])
        except OSError:  # if the directory already exists an OS error is raised
            pass  # don't care for logging
        new = "{}/{}.pf".format(grid[i], root)
        pfs.append(new)
        with open(new, "w") as f:
            for par, val in new_pf:
                f.write("{}\t\t{}\n".format(par, val))

    return pfs


def run_grid() -> None:
    """
    Main controlling function of the script.

    Returns
    -------
    None
    """

    dry_run = False

    if len(argv) == 2:
        if argv[1] == "--dry-run":
            dry_run = True
        else:
            print("don't know command {}".format(argv[1]))
    elif len(argv) > 2:
        print("only takes 1 command ;-)")
        print(__doc__)
        return

    print("ENSURE THAT THE SCRIPT HAS BEEN EDITED APPROPRIATELY BEFORE RUNNING")
    input("Press a enter to continue...")

    # This is the parameter which will be changed
    root = "../tde_cv.pf"
    parameter = "Wind.mdot(msol/yr)"

    grid = []
    # tmp = np.linspace(0.01, 0.1, 10)
    tmp = [2, 3, 1/2, 1/3]
    for i in range(len(tmp)):
        grid.append("{:.4e}".format(2e-2 * tmp[i]))

    print("Running grid of {} simulations:".format(len(grid)))
    print("Parameter: {}".format(parameter))
    print("Grid values: {}".format(grid))

    pfs = create_grid(root, parameter, grid)

    # Set up global variables for py_run modules
    py_run.PY_VERSION = "py"
    py_run.PY_FLAGS = "-p 2"
    py_run.RUN_SIMS = True
    py_run.RESUME_RUN = False
    py_run.CHECK_CONVERGENCE = True
    py_run.CLIM = 0.85
    py_run.CREATE_PLOTS = True
    py_run.TDE_PLOT = True
    py_run.NOT_QUIET = True
    py_run.VERBOSE = False
    py_run.SPEC_OVERRIDE = True

    # Finally make a call to run_python to run all of the parameter files
    if dry_run:
        return
    else:
        nsims = len(pfs)
        mpi, ncores = py_run_util.get_num_procs()
        py_run.run_python_etc(pfs, nsims, mpi, ncores)

    return


if __name__ == "__main__":
    run_grid()
