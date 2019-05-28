#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
from matplotlib import pyplot as plt
import py_plot_util
from consts import *


def plot_spec(root, dir: str = "./", xlims=None, plot_edges: bool = False, semilogy: bool = True) -> plt.Axes:
    """

    Parameters
    ----------
    root
    dir
    xlims
    plot_edges
    semilogy

    Returns
    -------

    """

    filename = "{}/{}.tau_spec.diag".format(dir, root, root)
    spec = np.loadtxt(filename, skiprows=1)
    with open(filename, "r") as f:
        cols = f.readline().split()
    if cols[0] == "#":
        cols = cols[1:]
    nangles = len(cols) - 1

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    xxlims = [spec[0, 0], spec[-1, 0]]
    if xlims:
        xxlims[0] = xlims[0]
        xxlims[1] = xlims[1]

    ymin = 1e+99
    ymax = 1e-99

    if plot_edges:
        scale = 1e6
    else:
        scale = 10

    for i in range(nangles):

        if semilogy:
            ax.semilogy(spec[:, 0], spec[:, i + 1], label=cols[i + 1])
        else:
            ax.plot(spec[:, 0], spec[:, i + 1], label=cols[i + 1])

        if plot_edges or semilogy is False:
            tymax, tymin = py_plot_util.get_ylimits(spec[:, 0], spec[:, i + 1], xxlims[0], xxlims[1], scale=scale)
            if tymax > ymax:
                ymax = tymax
            if tymin < ymin and tymin != 0:
                ymin = tymin

    ax.set_xlabel(r"Wavelength, $\AA$")
    ax.set_ylabel(r"Optical Depth, $\tau$")
    ax.set_xlim(xxlims[0], xxlims[1])

    if plot_edges or semilogy is False:
        ax.set_ylim(ymin, ymax)

    ax.legend()

    edges = {
        "He II Edge": 229,
        "Lyman Edge": 912,
        "Ly-delta": 950,
        "Ly-gamma": 973,
        "Ly-beta": 1026,
        "Ly-alpha": 1216,
        "Balmer Edge": 3646,
        "Ba-delta": 4101,
        "Ba-gamma": 4340,
        "Ba-beta": 4861,
        "Ba-alpha": 6564,
    }

    if plot_edges:
        py_plot_util.plot_line_ids(ax, edges, "vertical")

    plt.savefig("{}_tau_spec.png".format(root))
    plt.show()

    return ax


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Pls provide only the root name")
        sys.exit(1)
    else:
        root = sys.argv[1]

    plot_spec(root)
