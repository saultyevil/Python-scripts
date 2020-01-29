#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from PyPython import SpectrumUtils, Utils
from matplotlib import pyplot as plt
import numpy as np
from consts import *

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15


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
        # ax.axvline(x, ymin=0.75, ymax=0.98, linestyle="-", linewidth=2, color="k")
        xnorm = x - 0.07 * x
        ax.text(xnorm, 4e5, lab, ha="center", va="bottom", rotation="vertical", fontsize=13)

    return ax


def plot_optical_depth(tau_fname: str, disk_fname: str, extra: str, lims=None, los="60"):
    """
    Plot the optical depth as a function of wavelength.

    Parameters
    ----------
    tau_fname: str
        The file name and path for the optical depth spectrum.
    disk_fname: str
        The file name and path for the disc spectrum.
    extra: str
        An extra bit to add to the file name
    """

    root, wd = Utils.split_root_directory(tau_fname)
    fig, ax = plt.subplots(1, 1, figsize=(10, 7))

    disc_spec = SpectrumUtils.read_spec(disk_fname)
    disklamb = disc_spec["Lambda"].values.astype(float)
    diskhz = disc_spec["Freq."].values.astype(float)
    diskflux = disc_spec["Emitted"].values.astype(float) * 4 * PI * (100 * PARSEC) ** 2 * disklamb

    spectrum = SpectrumUtils.read_spec("solar/tde_cv.spec")
    spectrumlamb = spectrum["Lambda"].values.astype(float)
    spectrumhz = spectrum["Freq."].values.astype(float)
    spectrumflux = spectrum["Emitted"].values.astype(float) * 4 * PI * (100 * PARSEC) ** 2 * spectrumlamb

    # https://hea-www.harvard.edu/~pgreen/figs/Conversions.pdf
    #          Converts Flambda -> Fnu                                               so now it's nu * F_nu   and now it's nu * L nu
    # diskflux = (disc_spec["Lambda"].values.astype(float) ** 2 * diskflux * 3.34e-19) * diskhz * 4 * PI * (100 * PARSEC) ** 2

    spec = np.loadtxt(tau_fname)
    freq = spec[:, 1]
    # freq = C / (wl * ANGSTROM)  # Because the optical depth thing doesn't output frequency

    with open(tau_fname, "r") as f:
        angles = f.readline().split()
    if angles[0] == "#":
        angles = angles[2:]
    nangles = len(angles)
    for i in range(nangles):
        angles[i] = angles[i][1:3]
    angles = sorted(angles)
    if angles[-1] == "re":
        angles = angles[:-1]
    ninclinations = len(angles)

    ax2 = ax.twinx()
    ax2.loglog(diskhz, SpectrumUtils.smooth_spectrum(diskflux, 75), "k--", zorder=0, alpha=0.5)
    ax2.loglog(spectrumhz, SpectrumUtils.smooth_spectrum(spectrumflux, 75), "k-", zorder=1, alpha=0.5)
    ax2.set_ylabel(r"$\nu$ L$_{\nu}$ [ergs s$^{-1}$]", fontsize=15)

    columns = [0, 4.78e+25, 9.34e+25, 1.64e+26, 3.17e+26, 5.40e+26]
    for i in range(ninclinations - 1):
        ax.loglog(freq, spec[:, i + 3], label="i = {}"
                  .format(int(angles[i + 1])) + r"$^{\circ}$: $N_{H} \sim$" + "{:3.2e}".format(columns[i + 1])
                                              + " cm$^{-2}$")

    ax.legend(loc="center left")
    ax.set_ylabel(r"Optical Depth $\tau_{\nu}}$", fontsize=15)
    ax.set_xlabel(r"Frequency $\nu$ [Hz]", fontsize=15)
    ax.tick_params(axis="x")
    ax.tick_params(axis="y")
    ax.set_xlim(np.min(freq), np.max(freq))
    ax.text(5e14, 4e1, "Electron Scattering", fontsize=13)

    ax.set_zorder(ax2.get_zorder() + 1)
    ax.patch.set_visible(False)

    if lims:
        ax.set_ylim(lims[0], lims[1])
    else:
        ax.set_ylim(9e-6, 2e10)

    plot_line_id(ax, SpectrumUtils.absorption_edges(True))

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("{}_{}_optical_depth.png".format(root, extra))
    plt.show()

    return fig, ax


def main():
    """
    Main function of the script.
    """

    root = "/home/saultyevil/Dropbox/DiskWinds/PySims/tde/"

    fname = "solar/diag_tde_cv/tde_cv.tau_spec.diag"
    disk_name = root + "models/disc/cv/tde_cv.spec"
    plot_optical_depth(fname, disk_name, "clumpy", (1, 2e7))

    return


if __name__ == "__main__":
    main()
