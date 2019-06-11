#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
from matplotlib import pyplot as plt
import py_plot_util
from typing import Tuple
from consts import *


def plot_optical_depth_spec(root: str, dir: str = "./", plot_freq: bool = False, plot_edges: bool = True,
                            xlims: Tuple[int, int] = None, semilogy: bool = True, loglog: bool = False,
                            show_plot: bool = False) -> plt.Axes:
    """
    Plot an optical depth spectrum of either optical depth vs wavelength (Angstroms) or
    optical depth vs frequency (Hz).

    Parameters
    ----------
    root         str
                 The root name of the optical depth spectrum to plot
    dir          str
                 The directory containing the simulation
    plot_freq    bool, optional
                 If True, the optical depth will be printed in frequency space
    plot_edges   bool, optional
                 Label various absorption edges
    xlims        (float, flot), optional
                 The x-limits of the spectrum
    semilogy     bool, optional
                 If True, then the y axis will be log scaled
    loglog       bool, optional
                 If True, then the x and y axis will be log scaled
    show_plot    bool, optional
                 If True, then the plot will be showed before being saved to disk

    Returns
    -------
    ax           plt.Axes
                 A matplotlib Axes object of the optical depth spectrum
    """

    filename = "{}/{}.tau_spec.diag".format(dir, root, root, root)

    try:
        spec = np.loadtxt(filename, skiprows=1)
    except IOError:
        try:
            filename = "{}/diag_{}/".format(dir, root) + "{}.tau_spec.diag".format(root)
            spec = np.loadtxt(filename, skiprows=1)
        except IOError:
            print("Unable to find the optical depth spectrum: {}".format(root))
            sys.exit(1)

    with open(filename, "r") as f:
        cols = f.readline().split()
    if cols[0] == "#":
        cols = cols[1:]
    nangles = len(cols) - 1

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    if plot_freq:
        spec[:, 0] = C / (spec[:, 0] * ANGSTROM)

    if spec[0, 0] < spec[-1, 0]:
        xxlims = [spec[0, 0], spec[-1, 0]]
    else:
        xxlims = [spec[-1, 0], spec[0, 0]]

    if xlims:
        xxlims[0] = xlims[0]
        xxlims[1] = xlims[1]

    ymin = 1e+99
    ymax = 1e-99

    if plot_edges:
        scale = 50000
    else:
        scale = 10

    for i in range(nangles):
        the_label = r"$i$ = " + str(float(cols[i + 1][1:3])) + r"$^{\circ}$"
        if loglog:
            ax.semilogy(np.log10(spec[:, 0]), spec[:, i + 1], label=the_label)
        elif semilogy:
            ax.semilogy(spec[:, 0], spec[:, i + 1], label=the_label)
        else:
            ax.plot(spec[:, 0], spec[:, i + 1], label=the_label)
        ax.tick_params(axis="both", which="major", labelsize=13)

        if plot_edges or semilogy is False:
            tymax, tymin = py_plot_util.get_ylimits(spec[:, 0], spec[:, i + 1], xxlims[0], xxlims[1], scale=scale)
            if tymax > ymax:
                ymax = tymax
            if tymin < ymin and tymin != 0:
                ymin = tymin

    if plot_freq:
        if loglog:
            ax.set_xlabel(r"Log(Frequency), Hz")
        else:
            ax.set_xlabel(r"Frequency, Hz")
    else:
        ax.set_xlabel(r"Wavelength, $\AA$", fontsize=15)
    ax.set_ylabel(r"Optical Depth, $\tau$", fontsize=15)

    if plot_edges or semilogy is False:
        ax.set_ylim(ymin, ymax)

    ax.legend()

    if plot_edges:
        py_plot_util.plot_line_ids(ax, py_plot_util.get_common_absorption_edges(plot_freq, loglog), "horizontal")

    plt.savefig("{}_tau_spec.png".format(root))

    if show_plot:
        plt.show()
    else:
        plt.close()

    return ax


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Pls provide only the root name")
        sys.exit(1)
    else:
        root = sys.argv[1]

    plot_optical_depth_spec(root)
