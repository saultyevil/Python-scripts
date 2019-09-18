#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from typing import List
import py_run_util as pru
import py_plot_util as ppu
from os import access, R_OK


def check_multiple(wdpf: List[str]):
    """

    Parameters
    ----------
    wdpf: List[str]
        A list containing the directories of multiple pf files.
    """

    convergence = []

    for i in range(len(wdpf)):
        root, path = ppu.get_root_name(wdpf[i])
        c = pru.check_convergence(root, path, verbose=False)
        convergence.append(c)

    print("The following simulations have completed at least 1 cycle")
    for i in range(len(wdpf)):
        c = convergence[i]
        if c > 0:
            print(wdpf[i], c)

    print("The following simulations have no cycle information")
    for i in range(len(wdpf)):
        c = convergence[i]
        if c <= 0:
            print(wdpf[i], c)

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
    ppfs = ppu.find_pf("./")

    for i in range(len(ppfs)):
        pf, wd = ppu.get_root_name(ppfs[i])
        if root:
            if root == pf:
                pfs.append(ppfs[i])
        else:
            pfs.append(ppfs[i])

    return pfs


if __name__ == "__main__":
    check_multiple(get_pfs())
