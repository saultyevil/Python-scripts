#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to print out detailed convergence information -
with a specific focus on cells which have failed the specified number of
convergence checks in Python. One should then be able to study which cells are
failing to converge, or study cells which are not converging. The cell information
contained in a regular master file should be printed out. Note that this will
not work as well with all master files, but only those which contain the
converging output as well.
"""


import argparse as ap
import pandas as pd
from PyPython.PythonUtils import windsave2table
from typing import Union


def setup_script() -> Union[str, int, int]:
    """
    Parse setup options from the command line.

    Returns
    -------
    root: str
        The root name of the wind save file.
    nfails: int
        The number of checks failed to look at.
    converging: int
        Either 0 or 1 depending on if you want a converging or non-converging cells.
    """

    pd.set_option('display.expand_frame_repr', False)  # Expand the Pandas table to streeeetch

    p = ap.ArgumentParser(description=__doc__)

    p.add_argument("root",
                   help="Root name of the wind_save file")

    p.add_argument("nfails",
                   type=int,
                   help="Extract cells which failed this number of convergence checks")

    p.add_argument("converging",
                   nargs="?",
                   type=int,
                   help="Extract cells with the corresponding converging flag")

    args = p.parse_args()

    return args.root, args.nfails, args.converging


def main():
    """
    Main function. Returns the screen output as a Pandas DataFrame and writes to file.
    """

    root, nfails, converging = setup_script()
    input_file = "{}.0.master.txt".format(root)

    try:
        master = pd.read_csv(input_file, delim_whitespace=True)
    except IOError:
        windsave2table(root, "./", no_ep_complete=True)
        master = pd.read_csv(input_file, delim_whitespace=True)

    master = master[master["inwind"] == 0]
    screen_output = master[master["converge"] == nfails]
    try:
        screen_output = screen_output[screen_output["converging"] == converging]
    except KeyError:
        print("Can't find converging flag: probably doesn't exist... so skipping.")
    screen_output.to_csv("{}_nfails{}_converging{}.txt".format(root, nfails, converging), sep=" ")
    print(screen_output)

    return screen_output


if __name__ == "__main__":
    main()
