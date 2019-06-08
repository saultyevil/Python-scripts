#!/usr/bin/env python3
# -*- coding: utf-8 -*-


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
    py_plot_util.run_windsave2table(path, root, verbose)

    return


if __name__ == "__main__":
    root = "fiducial_agn.pf"
    path = "./"
    create(root, path)
