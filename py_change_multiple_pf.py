#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is basically an extension of the script py_change_parameter.py. The key
difference here is that the script will search recursively from the calling
directory for parameter files. If a root name is provided, however, then the
script will only operate on pf files which have the same root name.

The script expects 2 arguments and 1 optional argument, as documented below:

Usage
   $ python py_change_multiple_pf.py parameter value [root]

        - parameter       The name of the parameter to update
        - value           The updated parameter value
        - root            Optional: the name of the parameter files to edit
"""

import py_change_parameter as pcp
from sys import argv, exit
import py_plot_util as ppu
from typing import List


def change_pfs(wdpf: List[str], parameter: str, value: str) -> None:
    """
    Iterate over a list of pfs, and update the parameter given by the variable
    parameter with the new value given by value. This function will also
    print out verbose, because it seems most sensible to be loud about this.

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
        pcp.change_python_parameter(wdpf[i], parameter, value, verbose=True)

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
    ppfs = ppu.find_pf_files("./")

    print(root)

    for i in range(len(ppfs)):
        pf, wd = ppu.get_root_wd(ppfs[i])
        if root:
            if root == pf:
                pfs.append(ppfs[i])
        else:
            pfs.append(ppfs[i])

    return pfs


if __name__ == "__main__":
    root = None
    argc = len(argv)
    if argc == 3:
        parameter = argv[1]
        value = argv[2]
    elif argc == 4:
        parameter = argv[1]
        value = argv[2]
        root = argv[3]
    else:
        print(__doc__)
        exit(1)
    change_pfs(get_pfs(root), parameter, value)