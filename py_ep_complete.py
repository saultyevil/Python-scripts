#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Wrapper script for run_windsave2table
"""

import py_plot_util


def main():

    path = "./"
    root = "fiducial_agn"
    verbose = False
    py_plot_util.run_windsave2table(path, root, verbose)

    return


if __name__ == "__main__":
    main()
