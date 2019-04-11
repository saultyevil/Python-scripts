#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import py_run
import shutil
import py_run_util


def create_grid(pf, parameter, new_values):
    """

    Parameters
    ----------
    pf
        path to the original Python pf

    Returns
    -------

    """

    pfs = []

    if pf.find(".pf") == -1:
        pf += ".pf"

    sl = 0
    for i, l in enumerate(pf):
        if l == "/":
            sl = i
    pl = pf.find(".pf")
    root = pf[sl:pl]
    shutil.copyfile(pf, pf + ".bak")

    print(root)

    with open(pf, "r") as f:
        pf_lines = f.readlines()

    base = []
    for l in pf_lines:
        if l[0] == "#" or l[0] == "\n" or l[0] == "\r":
            continue
        line = l.replace("\n", "").split()
        base.append([line[0], line[1]])

    npfs = len(new_values)
    for i in range(npfs):
        new_pf = base.copy()
        for l in range(len(new_pf)):
            if new_pf[l][0] == parameter:
                new_pf[l][1] = new_values[i]
        try:
            os.mkdir(new_values[i])
        except OSError:
            pass
        new = "{}/{}.pf".format(new_values[i], root)
        print(new)
        pfs.append(new)
        with open(new, "w") as f:
            for par, val in new_pf:
                f.write("{}\t\t{}\n".format(par, val))

    return pfs


def run_grid(pfs):
    """

    """

    nsims = len(pfs)
    mpi, ncores = py_run_util.get_num_procs()
    py_run.run_python_etc(pfs, nsims, mpi, ncores)

    return


import numpy as np  ## fuck uuuuuuuuuuuuuuuuuuuuuuuuuuu


def main():
    """

    """

    parameter = "Wind.mdot(msol/yr)"
    values = []
    baa = np.linspace(0.01, 0.1, 10)
    for i in range(len(baa)):
        values.append("{:.4e}".format(2e-1 * baa[i]))

    print(values)

    pfs = create_grid("tde.pf", parameter, values)

    py_run.RUN_SIMS = True
    py_run.CHECK_CONVERGENCE = True
    py_run.CREATE_PLOTS = True
    py_run.TDE_PLOT = True
    py_run.PY_FLAGS = "-p 2"

    run_grid(pfs)

    return


if __name__ == "__main__":
    main()
