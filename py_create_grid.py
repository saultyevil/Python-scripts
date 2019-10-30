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
import shutil
import numpy as np
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
        if line == "/":  # TODO: rfind exists u nonce
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
            elif new_pf[line][0] == "Wind.mdot(msol/yr)":
                new_pf[line][1] = str(0.1 * float(grid[i]))
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


def run_grid() -> List[str]:
    """
    Main controlling function of the script.

    Returns
    -------
    None
    """

    # This is the parameter which will be changed
    root = "../tde_cv.pf"
    parameter = "SV.acceleration_exponent"

    tmp = [0.7, 0.9, 1.1, 1.3]
    grid = []
    for i in range(len(tmp)):
        grid.append("{:.4e}".format(tmp[i]))
    print(parameter, grid)

    print("ENSURE THAT THE SCRIPT HAS BEEN EDITED APPROPRIATELY BEFORE RUNNING")
    input("Press a enter to continue...")

    print("Running grid of {} simulations:".format(len(grid)))
    print("Parameter: {}".format(parameter))
    print("Grid values: {}".format(grid))


    pfs = create_grid(root, parameter, grid)

    print(pfs)

    return pfs


if __name__ == "__main__":
    run_grid()
