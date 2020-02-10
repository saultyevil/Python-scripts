#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to automatically generate a *.slurm file for a
Python simulation given sensible inputs. This script can also be used to
update an already existing .slurm file, for example if one wishes to restart
a Python simulation.
"""


import argparse
from typing import Tuple


def write_slurm_file(name: str, ncores: int, thours: int, tminutes: int, flags: str, wd: str = "./") -> None:
    """
    Create a slurm file in the directory wd with the name root.slurm. All
    of the script flags are passed using the flags variable.

    Parameters
    ----------
    name: str
        The name of the slurm file
    ncores: int
        The number of cores which to use
    thours: int
        The number of hours to allow
    tminutes: int
        The number of minutes to allow
    flags: str
        The run-time flags of which to execute Python with
    wd: str
        The directory to write the file to
    """

    slurm = \
        """#!/bin/bash
#SBATCH --mail-user=ejp1n17@soton.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --ntasks={}
#SBATCH --time={}:{}:00
#SBATCH --partition=batch
module load openmpi/3.0.0/gcc
module load conda/py3-latest
source activate PyPython
python /home/ejp1n17/PythonScripts/py_run.py -n {} {}
""".format(ncores, thours, tminutes, ncores, flags)

    if wd[-1] != "/":
        wd += "/"
    fname = wd + name + ".slurm"
    with open(fname, "w") as f:
        f.write("{}".format(slurm))

    return


def parse_arguments() -> Tuple[str, int, int, int, str]:
    """
    Parse arguments from the command line.

    Returns
    -------
    args.name: str
        The name of the slurm file
    args.root: str
        The root name of the Python simulation
    args.ncores: int
        The number of CPUs to use
    args.thours: int
        The maximum run time allowed + 1 hours
    args.flags: str
        Any flags to pass to Python
    args.vers: str
        The version of Python to use
    """

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("name", type=str, help="The name of the slurm file, i.e. name.slurm.")
    p.add_argument("ncores", type=int, help="The number of CPUs to use.")
    p.add_argument("thours", type=int, help="The number of hours of run time allowed.")
    p.add_argument("tminutes", type=int, help="The number of minutes of additional run time allowed.")
    p.add_argument("-f", "--flags", type=str, help="Any flags to pass to the py_run.py Python running script.")
    args = p.parse_args()

    return args.name, args.ncores, args.thours, args.tminutes, args.flags


def main() -> None:
    """
    Main function of the script. Parses the arguments from the command line and
    then executes the function to generate the slurm file.
    """

    name, ncores, thours, tminutes, flags = parse_arguments()

    if flags is None:
        flags = ""
    flags += " -t {} ".format(int(thours * 3600 + tminutes * 60))

    write_slurm_file(name, ncores, thours, tminutes, flags)

    return


if __name__ == "__main__":
    main()
