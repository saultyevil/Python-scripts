#!/usr/bin/env python3

import py_util
import argparse
import numpy as np
from matplotlib import pyplot as plt
from subprocess import PIPE, Popen

verbose = False
showplot = False

def get_rootname():
    """


    Parameters
    ----------


    Returns
    -------

    """

    p = argparse.ArgumentParser(description="")
    p.add_argument("rootname", type=str, help="The rootname of the Python simulation.")
    p.add_argument("-v", "--verbose", help="Increase output to screen", action="store_true")
    p.add_argument("-s", "--show", help="Show plots to screen", action="store_true")
    args = p.parse_args()

    if args.verbose:
        global verbose
        verbose = True
    if args.show:
        global showplot
        showplot = True

    return args.rootname


def plot_temp_struc():
    """


    Parameters
    ----------


    Returns
    -------

    """

    return


if __name__ == "__main__":
    rootname = get_rootname()
    tefile = "{}.te.dat".format(rootname)
    tedata = py_util.read_file(tefile)
    headers = tedata[0, :]
    tedata = np.array(tedata[1:, :], dtype=float)

    plot_temp_struc()



