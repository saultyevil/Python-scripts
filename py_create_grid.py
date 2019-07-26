#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Create and run a grid of Python simulations. Note that this script has to be
edited before being called as I didn't want to add one billion command line
arguments/switches to the script like with py_run.

Usage

    [python] py_create_grid.py [--run_grids]

    --full_run: as well as creating a grid of runs, also run them
"""


import os
import py_run
import shutil
import numpy as np
import py_run_util
from sys import argv
from typing import List
from consts import *


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
    for i, line in enumerate(pf):  # This iterates over the pf file path
        if line == "/":
            sl = i  # This will find the index of the final /
    pl = pf.find(".pf")
    root = pf[sl:pl]  # Now we can extract just the root name from the file path
    shutil.copyfile(pf, pf + ".bak")

    with open(pf, "r") as f:
        pf_lines = f.readlines()

    base = []
    for line in pf_lines:
        if line[0] == "#" or line[0] == "\n" or line[0] == "\r":
            continue
        line = line.replace("\n", "").split()
        base.append([line[0], line[1]])

    # This bit creates the grid of pfs
    pfs = []
    npfs = len(grid)
    for i in range(npfs):
        new_pf = base.copy()
        for line in range(len(new_pf)):
            if new_pf[line][0] == parameter:
                new_pf[line][1] = grid[i]
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

    dry_run = True
    full_run = False

    if len(argv) == 2:
        if argv[1] == "--run_grid":
            full_run = True
            dry_run = False
        else:
            print("don't know command {}".format(argv[1]))
    elif len(argv) > 2:
        print("only takes 1 command which is --run_grid ;-)")
        print(__doc__)
        return

    print("ENSURE THAT THE SCRIPT HAS BEEN EDITED APPROPRIATELY BEFORE RUNNING")
    input("Press a enter to continue...")

    # This is the parameter which will be changed
    root = "tde_spherical.pf"
    parameter = "Stellar_wind.v_infinity(cm)"

    grid = []
    # tmp = np.linspace(0.0, 2.0, 10)
    tmp = [0.07, 0.075, 0.08, 0.085, 0.090, 0.095]
    esc_vel = np.sqrt(2 * G * 3e7 * MSOL / 2.65e13)
    for i in range(len(tmp)):
        grid.append("{:.4e}".format(esc_vel * tmp[i]))

    print("Running grid of {} simulations:".format(len(grid)))
    print("Parameter: {}".format(parameter))
    print("Grid values: {}".format(grid))

    pfs = create_grid(root, parameter, grid)

    # Finally make a call to run_python to run all of the parameter files
    if dry_run:
        return
    if full_run:
        # Set up global variables for py_run modules
        py_run.PY_VERSION = "py"
        py_run.PY_FLAGS = "-p 2"
        py_run.RUN_SIMS = True
        py_run.RESUME_RUN = False
        py_run.CHECK_CONVERGENCE = True
        py_run.CONV_LIMIT = 0.85
        py_run.CREATE_PLOTS = True
        py_run.TDE_PLOT = True
        py_run.NOT_QUIET = True
        py_run.VERBOSE = False
        py_run.SPLIT_CYCLES = True
        # Run Python using py_run.py
        nsims = len(pfs)
        mpi, ncores = py_run_util.get_num_procs()
        py_run.go(pfs, mpi, ncores)

    return


if __name__ == "__main__":
    run_grid()
