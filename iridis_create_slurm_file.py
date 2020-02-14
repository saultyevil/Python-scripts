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


def write_slurm_file(name: str, ncores: int, split_cycle: bool, thours: int, tminutes: int, root: str, flags: str, wd: str = "./") \
        -> None:
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
    split_cycle: bool
        If True, then py_run will use the split_cycle method
    flags: str
        The run-time flags of which to execute Python with
    root: str
        The root name of the model
    wd: str
        The directory to write the file to
    """

    split = ""
    if split_cycle:
        split = "-sc"

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
python /home/ejp1n17/PythonScripts/py_run.py -n {} {} -f="{}"
python /home/ejp1n17/PythonScripts/py_analyse_run.py {}
""".format(ncores, thours, tminutes, ncores, split, flags, root, root)

    if wd[-1] != "/":
        wd += "/"
    fname = wd + name + ".slurm"
    with open(fname, "w") as f:
        f.write("{}".format(slurm))

    return


def parse_arguments() -> Tuple[str, int, int, int, bool, str, str]:
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
    p.add_argument("root", type=str, help="The root name of the model.")
    p.add_argument("-f", "--flags", type=str, help="Any flags to pass to the py_run.py Python running script.")
    p.add_argument("-sc", "--split_cycle", action="store_true", help="Use the split cycle method for py_run.py")
    args = p.parse_args()

    split_cycle = False
    if args.split_cycle:
        split_cycle = True

    return args.name, args.ncores, args.thours, args.tminutes, split_cycle, args.root, args.flags


def main() -> None:
    """
    Main function of the script. Parses the arguments from the command line and
    then executes the function to generate the slurm file.
    """

    name, ncores, thours, tminutes, split_cycle, root, flags = parse_arguments()

    if flags is None:
        flags = ""
    flags += " -t {} ".format(int(thours * 3600 + tminutes * 60))

    write_slurm_file(name, ncores, split_cycle, thours, tminutes, root, flags)

    return


if __name__ == "__main__":
    main()
