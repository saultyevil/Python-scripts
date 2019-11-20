#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The purpose of this is to mass delete a bunch of wind_save files using a
standard name, i.e. python_01.wind_save
"""

from sys import argv
from subprocess import Popen, PIPE
from typing import List
from PyPython import Utils


def delete_windsaves(wdpf: List[str], root) -> None:
    """
    """

    for i in range(len(wdpf)):
        pf, wd = Utils.split_root_directory(wdpf[i])
        cmd = "cd {}; rm {}*.wind_save".format(wd, root)
        sh = Popen(cmd, stdout=True, stderr=True, shell=True)
        stdout, stderr = sh.communicate()
        if stderr:
            print(stderr.decode("utf-8"))
    
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

if __name__ == "__main__":
    delete_windsaves(get_pfs(), argv[1])
