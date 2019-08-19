#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to automatically generate a *.slurm file for a
Python simulation given sensible inputs. This script can also be used to
update an already existing .slurm file, for example if one wishes to restart
a Python simulation.

usage: iridis_create_slurm_file.py [-h] [-flags FLAGS] [-vers VERS]
                                   name root ncores thours

positional arguments:
  name          The name of the slurm file, i.e. name.slurm
  root          The root name of the Python simulation
  ncores        The number of CPUs to use
  thours        The maximum run time to run

optional arguments:
  -h, --help    show this help message and exit
  -flags FLAGS  Any flags to pass to Python
  -vers VERS    The version of Python to use
"""


import argparse
from typing import Tuple


def write_slurm_file(name: str, root: str, ncores: int, thours: int, flags: str, vers: str = "py") -> None:
    """
    Create a slurm file in the directory wd with the name root.slurm.

    Parameters
    ----------
    name: str
        The name of the slurm file
    root: str
        The root name of the Python simulation
    ncores: int
        The number of cores which to use
    thours: int
        The number of hours to execute for
    flags: str
        The flags of which to execute Python with
    vers: str, optional
        The version of Python to use, i.e. py83c
    """

    slurm = \
        """#!/bin/bash
#SBATCH --mail-user=ejp1n17@soton.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --ntasks={}
#SBATCH --time={}:00:00
#SBATCH --partition=batch
module load openmpi/3.0.0/gcc
N_TASKS="{}"
DIR=$(pwd)
PY_VER="{}"
PY_FLAGS="-t {} {}"
ROOT="{}"
export PYTHON="$HOME/python"
export PYTHON_BIN="$HOME/python/bin"
cd $DIR
$PYTHON_BIN/Setup_Py_Dir
mpirun -np $N_TASKS $PYTHON_BIN/$PY_VER $PY_FLAGS $ROOT >> $ROOT_out.txt
$HOME/Scripts/py_check_run.py $ROOT""".format(ncores, thours, ncores, vers, thours * 3600 - 60, flags, root)

    fname = name + ".slurm"
    with open(fname, "w") as f:
        f.write("{}\n".format(slurm))

    return


def parse_arguments() -> Tuple[str, str, int, int, str, str]:
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

    p = argparse.ArgumentParser(description="Create a slurm file to submit to the Iridis queue")
    p.add_argument("name", type=str, help="The name of the slurm file, i.e. name.slurm")
    p.add_argument("root", type=str, help="The root name of the Python simulation")
    p.add_argument("ncores", type=int, help="The number of CPUs to use")
    p.add_argument("thours", type=int, help="The maximum run time allowed")
    p.add_argument("-flags", type=str, help="Any flags to pass to Python")
    p.add_argument("-vers", type=str, help="The version of Python to use")
    args = p.parse_args()

    return args.name, args.root, args.ncores, args.thours, args.flags, args.vers


def main() -> None:
    """
    Main function of the script. Parses the arguments from the command line and
    then executes the function to generate the slurm file.
    """

    name, root, ncores, thours, flags, vers = parse_arguments()
    if vers is None:
        vers = "py"
    if flags is None:
        flags = ""
    write_slurm_file(name, root, ncores, thours, flags, vers)

    return


if __name__ == "__main__":
    main()
