#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from PyPython import SpectrumUtils
from matplotlib import pyplot as plt
import numpy as np
from consts import *


def plot_line_id(ax: plt.Axes, LINES) -> plt.Axes:
    """
    Plot labels and vertical lines to indicate important atomic transitions.

    Parameters
    ----------
    ax: plt.Axes
        The Axes object to add line ID labels to.


    Returns
    -------
    ax: plt.Axes
        The input Axes object with additional line IDs.
    """

    xlims = ax.get_xlim()
    nlines = len(LINES)

    for i in range(nlines):
        lab = LINES[i][0]
        x = LINES[i][1]
        if x < xlims[0]:
            continue
        elif x > xlims[1]:
            continue
        ax.axvline(x, ymax=1, linestyle="--", linewidth=0.45, color="k")
        xnorm = x - 0.05 * x
        ax.text(xnorm, 2e8, lab, ha="center", va="center", rotation="vertical", fontsize=13) #, transform=ax.transAxes,)

    return ax


def plot_optical_depth(fname):
    """Plot the optical depth as a function of wavelength"""

    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    spec = np.loadtxt(fname)
    wl = spec[:, 0]

    freq = C / (wl * ANGSTROM)

    with open(fname, "r") as f:
        angles = f.readline().split()
    if angles[0] == "#":
        angles = angles[2:]
    nangles = len(angles)
    for i in range(nangles):
        angles[i] = angles[i][1:3]
    angles = sorted(angles)
    ninclinations = len(angles)

    for i in range(ninclinations - 1):
        ax.loglog(freq, spec[:, i + 2], label="i = {}".format(int(angles[i + 1])) + r"$^{\circ}$")
    ax.legend(loc="upper right")
    ax.set_ylabel(r"Optical Depth $\tau_{\nu}}$", fontsize=15)
    ax.set_xlabel(r"Frequency $\nu$ [Hz]", fontsize=15)
    ax.tick_params(axis="x", labelsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_xlim(np.min(freq), np.max(freq))
    ax.set_ylim(9e-6, 2e10)

    edges = SpectrumUtils.absorption_edges(True)
    for i in range(len(edges)):
        print("{} {:3.2e}".format(edges[i][0], edges[i][1]))

    ax = plot_line_id(ax, edges)

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("tde_cv_optical_depth.png")
    plt.show()

    return


def main():

    fname = "/home/saultyevil/Dropbox/DiskWinds/PySims/tde/paper_models/clump/1e-1/cv/solar/diag_tde_cv/" \
            "tde_cv.tau_spec.diag"
    plot_optical_depth(fname)

    return


if __name__ == "__main__":
    main()
