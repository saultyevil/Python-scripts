#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to be able to parse out where luminosity is
being lost in a Python model by trawling through the diagnostic files.
"""

import argparse as ap
from glob import glob
from tqdm import tqdm


def setup_script():
    """
    Parse arguments from the command line.
    """

    p = ap.ArgumentParser(description=__doc__)

    p.add_argument("root",
                   help="The root name of the Python simulation.")

    p.add_argument("-wd",
                   "--working_directory",
                   default=".",
                   help="The directory containing the Python simulation.")

    args = p.parse_args()

    setup = (
        args.root,
        args.working_directory
    )

    return setup


def read_diag(sinks: tuple, fname: str) -> tuple:
    """
    Read the provided diag file for the sinks of luminosity.
    """

    try:
        with open(fname, "r") as f:
            lines = f.readlines()
    except IOError:
        return sinks

    for n in range(len(lines)):
        line = lines[n]

        if line.find("!!python: luminosity lost by adiabatic kpkt destruction") != -1:
            sinks["adiabatic"] += float(line.split()[-1])
            sinks["low_freq"] += float(lines[n].split()[-1])
            n += 1

        if line.find("!!python: luminosity lost by being completely absorbed") != -1:
            sinks["ncycles"] += 1
            sinks["absorbed"] += float(line.split()[-1])
            sinks["max_scatter"] += float(lines[n + 1].split()[-1])
            sinks["star"] += float(lines[n + 2].split()[-1])
            sinks["disk"] += float(lines[n + 3].split()[-1])
            sinks["errors"] += float(lines[n + 4].split()[-1])
            sinks["unknown"] += float(lines[n + 5].split()[-1])
            n += 5

    return sinks


def main():
    """
    Main function of the script.
    """

    sinks = {
        "ncycles": 0,
        "adiabatic": 0,
        "low_freq": 0,
        "absorbed": 0,
        "max_scatter": 0,
        "star": 0,
        "disk": 0,
        "errors": 0,
        "unknown": 0
    }

    root, wd = setup_script()
    glob_directory = "{}/diag_{}/{}_*.diag".format(wd, root, root)
    diag_files = glob(glob_directory)
    diag_files = sorted(diag_files)

    for diag in tqdm(diag_files, desc="Reading diag files for {}".format(root), unit="files", smoothing=0):
        sinks = read_diag(sinks, diag)

    for key in sinks:
        sinks[key] /= len(diag_files)

    print(sinks)

    return


if __name__ == "__main__":
    main()
