#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of the script is to create a root.ep.complete file which uses
windsave2table to first generate the master, heat and ion tables. The master
and heat table and then combined into a root.ep.complete file.

Usage:

    Two arguments are expected,

    [python] py_ep_complete.py root path

        - root      The root name of the simulation
        - path      The directory containing the simulation, this can
                    just be "./"
"""


from sys import argv, exit
import py_plot_util


def create(root, path):
    """
    A basic wrapper function to run windsave2table and create one of my
    root.ep.complete files :-)

    Parameters
    ----------
    root           str
                   The root name of the simulation
    path           str
                   The directory containing the simulation

    Returns
    -------
    None
    """

    verbose = False
    py_plot_util.windsave2table(root, path, verbose)

    return


if __name__ == "__main__":
    argc = len(argv)
    if argc == 3:
        root = argv[1]
        path = argv[2]
    else:
        print(__doc__)
        exit(1)
    create(root, path)
