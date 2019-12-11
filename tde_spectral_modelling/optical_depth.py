#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from PyPython import SpectrumUtils, Utils
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
        ax.text(xnorm, 2e8, lab, ha="center", va="center", rotation="vertical", fontsize=13)

    return ax


def plot_optical_depth(tau_fname: str, disk_fname: str):
    """
    Plot the optical depth as a function of wavelength.

    Parameters
    ----------
    tau_fname: str
        The file name and path for the optical depth spectrum.
    disk_fname: str
        The file name and path for the disc spectrum.

    """

    root, wd = Utils.split_root_directory(tau_fname)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    disc_spec = SpectrumUtils.read_spec(disk_fname)
    diskhz = disc_spec["Freq."].values.astype(float)
    diskflux = disc_spec["60"].values.astype(float)
    # https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf
    #          Converts Flambda -> Fnu                                               so now it's nu * F_nu
    diskflux = (disc_spec["Lambda"].values.astype(float) ** 2 * diskflux * 3.34e-19) * diskhz # * 4 * PI * (100 * PARSEC) ** 2

    spec = np.loadtxt(tau_fname)
    wl = spec[:, 0]
    freq = C / (wl * ANGSTROM)  # Because the optical depth thing doesn't output frequency

    with open(tau_fname, "r") as f:
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

    ax2 = ax.twinx()
    ax2.loglog(diskhz, SpectrumUtils.smooth_spectrum(diskflux, 75), "--", zorder=0)
    ax2.set_ylabel(r"$\nu$ F$_{\nu}$ [ergs s$^{-1}$ cm$^{-2}$]", fontsize=15)

    ax.legend(loc="upper right")
    ax.set_ylabel(r"Optical Depth $\tau_{\nu}}$", fontsize=15)
    ax.set_xlabel(r"Frequency $\nu$ [Hz]", fontsize=15)
    ax.tick_params(axis="x", labelsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_xlim(np.min(freq), np.max(freq))
    ax.set_ylim(9e-6, 2e10)
    plot_line_id(ax, SpectrumUtils.absorption_edges(True))

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("{}_optical_depth.png".format(root))
    plt.show()

    return fig, ax


def main():
    """
    Main function of the script.
    """

    root = "/home/saultyevil/Dropbox/DiskWinds/PySims/tde"

    fname = root + "/paper_models/clump/1e-1/cv/solar/diag_tde_cv/tde_cv.tau_spec.diag"
    disk_name = root + "/matrix_bb/paper_models/disk_spectra/cv/tde_cv.spec"
    plot_optical_depth(fname, disk_name)

    return


if __name__ == "__main__":
    main()
