#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The original purpose of this script was to plot various important quantities for
cells which are having trouble converging. However, there is nothing to stop one
from using this for any cell in a model. The input files in general are expected
to be in the form of rootXX.0.master.txt, where XX is a series of numbers.
"""


import glob
import argparse as ap
import pandas as pd
from matplotlib import pyplot as plt
from typing import Union
import numpy as np


def setup_script() -> Union[str, int, int]:
    """
    Parse the cell number from the command line and other setup options.

    Returns
    -------
    root: str
        The root name of the input files.
    icell: int
        The i index for the cell to query.
    jcell: int
        The j index for the cell to query.
    """

    p = ap.ArgumentParser(description=__doc__)

    p.add_argument("root",
                   help="The root name of master files to query")

    p.add_argument("icell",
                   type=int,
                   help="The i index number of the cell")

    p.add_argument("jcell",
                   type=int,
                   help="the j index number of the cell")

    args = p.parse_args()

    return args.root, args.icell, args.jcell


def main():
    """
    Main function of the script. Returns fuck all.
    """

    root, icell, jcell = setup_script()
    input_files = sorted(glob.glob("./{}*.0.master.txt".format(root)), key=str.lower)
    if len(input_files) == 0:
        print("Couldn't find any input you dumb fuck try again idiot")
        return

    converge = []
    converging = []
    te = []
    tr = []
    ip = []
    c4 = []
    ntot = []
    ne = []
    title = ["converge", "converging", "te", "tr", "ip", "c4 frac", "ntot", "ne"]
    plot = [converge, converging, te, tr, ip, c4, ntot, ne]

    for file in input_files:
        master = pd.read_csv(file, delim_whitespace=True)
        # Indexing and masking fuckery :)
        cell = master[master["i"] == icell]
        cell = cell[cell["j"] == jcell]
        cell = cell[cell["inwind"] == 0]
        # Pull out the relevant stuff from the dataframe
        converge.append(int(cell["converge"].values))
        te.append(float(cell["t_e"].values))
        tr.append(float(cell["t_r"].values))
        ip.append(float(cell["ip"].values))
        c4.append(float(cell["c4"].values))
        ntot.append(float(cell["ntot"].values))
        ne.append(float(cell["ne"].values))
        # This one may not exist in older master tables
        try:
            converging.append(int(cell["converging"].values))
        except KeyError:
            converging.append(0)

    nrows = 4
    ncols = 2
    fig, ax = plt.subplots(nrows, ncols, figsize=(13, 15), squeeze=False, sharex="col")
    cycles = np.arange(1, len(converge) + 1)

    ii = 0
    for i in range(nrows):
        for j in range(ncols):
            if ii == len(plot):
                break
            ax[i, j].plot(cycles, plot[ii], "k-")
            if title[ii] in ["te", "tr", "ip", "c4 frac", "ne"]:
                ax[i, j].set_yscale("log")
            ax[i, j].text(0.75, 0.1, title[ii], transform=ax[i, j].transAxes, fontsize=14)
            ax[i, j].set_xlim(1, len(converge) + 1)
            ii += 1

    fig.tight_layout(rect=[0.02, 0.05, 0.98, 0.95])
    fig.suptitle("icell = {} jcell = {}".format(icell, jcell), fontsize=14)
    plt.savefig("icell{}_jcell{}_convergence.png".format(icell, jcell))
    plt.close()

    return


if __name__ == "__main__":
    main()
