#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Add a parameter to an already existing parameter file.

The script will search recursively from the calling directory for parameter
files. If a root name is provided, however, then the script will only operate
 on pf files which have the same root name.

The script expects 2 arguments and 1 optional argument, as documented below:

Usage
   $ python py_change_multiple_pf.py parameter value [root] [-h] [-checkpf]

        - parameter       The name of the parameter to update
        - value           The updated parameter value
        - root            Optional: the name of the parameter files to edit
        - h               Print this error message and exit
        - checkpf         Prints to the screen the pfs which will be updated
"""


from sys import argv, exit
from typing import List
from PyPython import PythonUtils as Utils
from PyPython import Grid


def add_parameter(wdpf: List[str], parameter: str, value: str):
    """
    Iterate over a list of pfs, and add the parameter given to the end of
    the parameter file. This function will also print out verbose, because it
    seems the most sensible to be loud about this.

    Parameters
    ----------
    wdpf: List[str]
        A list containing the directories of multiple pf files.
    parameter: str
        The parameter name of the parameter which is being updated.
    value: str
        The updated parameter value.
    """

    for i in range(len(wdpf)):
        Grid.add_parameter(wdpf[i], parameter, value, verbose=True)

    return


def get_pfs(root: str = None) -> List[str]:
    """
    Search recursively from the calling directory for Python pfs. If root is
    specified, then only pfs with the same root name as root will be returned.

    Parameters
    -------
    root: str, optional
        If this is set, then any pf which is not named with this root will be
        removed from the return pfs

    Returns
    -------
    pfs: List[str]
        A list containing the relative paths of the pfs to be updated.
    """

    pfs = []
    ppfs = Utils.find_parameter_files("./")

    for i in range(len(ppfs)):
        pf, wd = Utils.split_root_directory(ppfs[i])
        if root:
            if root == pf:
                pfs.append(ppfs[i])
        else:
            pfs.append(ppfs[i])

    return pfs


def main(argc: int, argv: List[str]):
    """
    Main function.

    Parameters
    ----------
    argc: int
        The number of command line arguments provided.
    argv: List[str]
        The command line arguments provided.
    """

    root = None
    if argc == 2 and argv[1] == "-h":
        print(__doc__)
        exit(0)
    elif argc == 2 and argv[1] == "-checkpf":
        pfs = get_pfs()
        print("Will operate on the following pfs:\n", pfs)
        exit(0)
    elif argc == 3:
        parameter = argv[1]
        value = argv[2]
    elif argc == 4:
        parameter = argv[1]
        value = argv[2]
        root = argv[3]
    else:
        print("Unknown arguments provided: ", argv[1:])
        print(__doc__)
        exit(1)
    add_parameter(get_pfs(root), parameter, value)

    return


if __name__ == "__main__":
    main(len(argv), argv)
