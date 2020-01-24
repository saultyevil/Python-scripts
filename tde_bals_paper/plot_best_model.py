#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

"""
Plot a comparison between the best tde model against the data for a BAL and
BEL TDE.
"""

import sys
from sys import argv
from platform import system
from matplotlib import pyplot as plt
from typing import List
from PyPython import SpectrumUtils

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

import tde_spectra as tu
from consts import *

SMOOTH = 10

LINES = [
    ["O VI", 1032],
    ["P V", 0],
    ["N III]", 0],
    [r"Ly$\alpha$/N V", 1216],
    ["", 1242],
    ["O I", 0],
    ["O V/Si IV", 1371],
    ["", 1400],
    ["N IV]", 1489],
    ["C IV",  1549],
    ["He II", 1640],
    ["O III]", 0],
    ["N III]", 1750],
    ["C III]", 1908],
    ["Fe II", 0],
    ["Fe II / CII]", 0],
    ["Fe II", 0],
    ["Fe II", 0],
    ["Mg II", 2798],
    ["He II", 4686],
    [r"H$_{\beta}$", 4861],
    [r"H$_{\alpha}$", 6564],
]


def plot_line_id(ax: plt.Axes) -> plt.Axes:
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

        ax.axvline(x, ymin=0.8, ymax=0.98, linestyle="-", linewidth=1.5, color="k")
        x = x - 25
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        ax.text(xnorm, 0.90, lab, ha="center", va="center", rotation="vertical", transform=ax.transAxes, fontsize=13)

    return ax


def plot_against_data(dfname: str, disk: str, inclinations: List[str], name: str):
    """

    Parameters
    ----------
    dfname: str
        The file name of the models to plot
    disk: str
        The file names of the disk conts. to plot
    inclinations: List[str]
        The inclinations to plot - this should be the inclination for the BAL
        followed by the BEL model
    name: str
        A identification name for the plot to make. Allowed values are smooth
        and clump.
    """

    alpha = 0.8
    wmin = 950
    wmax = 3050

    iptf = tu.iptf15af_spec(SMOOTH)
    assa = tu.asassn14li_spec(SMOOTH)
    fig, ax = plt.subplots(2, 1, figsize=(9.5, 12), sharex="col")

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"

    dfname = pdir + dfname
    spec1 = SpectrumUtils.read_spec(dfname)
    disc_spec = SpectrumUtils.read_spec(pdir + disk)

    inclination = inclinations[0]
    spec1_wl = spec1["Lambda"].values.astype(float)
    spec1_flux = spec1[inclination].values.astype(float)
    spec1_flux *= (100 * PARSEC) ** 2 / (350 * 1e6 * PARSEC) ** 2
    disk_wl = disc_spec["Lambda"].values.astype(float)
    disc_flux = disc_spec[inclination].values.astype(float)
    disc_flux *= (100 * PARSEC) ** 2 / (350 * 1e6 * PARSEC) ** 2

    ax[0].semilogy(spec1_wl, SpectrumUtils.smooth_spectrum(spec1_flux, SMOOTH), color="C0", linewidth=3,
                   label=r"Model: i = {}".format(inclination) + r"$^{\circ}$", alpha=alpha, zorder=4)
    ax[0].semilogy(disk_wl, SpectrumUtils.smooth_spectrum(disc_flux, 15), "--", color="C0", linewidth=1, alpha=alpha)
    ax[0].semilogy(iptf[:, 0] / (0.07897 + 1), SpectrumUtils.smooth_spectrum(iptf[:, 1], SMOOTH),
                   "k", label=r"iPTF15af $\Delta t = $52 d", zorder=2)
    ax[0].set_xlim(wmin, wmax)
    ax[0].set_ylim(2e-17, 4e-14)
    ax[0] = plot_line_id(ax[0])
    ax[0].legend(loc="lower right", fontsize=13)
    ax[0].text(0.05, 0.05, "In-wind", ha="left", va="center", rotation="horizontal", fontsize=15,
                  transform=ax[0].transAxes)

    colors = ["C0", "C1"]

    for i in range(len(inclinations[1])):
        inclination = inclinations[1][i]
        spec1_wl = spec1["Lambda"].values.astype(float)
        spec1_flux = spec1[inclination].values.astype(float)
        spec1_flux *= (100 * PARSEC) ** 2 / (90 * 1e6 * PARSEC) ** 2
        disk_wl = disc_spec["Lambda"].values.astype(float)
        disc_flux = disc_spec[inclination].values.astype(float)
        disc_flux *= (100 * PARSEC) ** 2 / (90 * 1e6 * PARSEC) ** 2

        ax[1].semilogy(spec1_wl, SpectrumUtils.smooth_spectrum(spec1_flux, SMOOTH), color=colors[i], linewidth=3,
                       label=r"Model: i = {}".format(inclination) + r"$^{\circ}$", alpha=alpha, zorder=4)
        ax[1].semilogy(disk_wl, SpectrumUtils.smooth_spectrum(disc_flux, 15), "--", color=colors[i], linewidth=1,
                       alpha=0.8)

    ax[1].semilogy(assa[:, 0] / (0.02058 + 1), SpectrumUtils.smooth_spectrum(assa[:, 1], SMOOTH),
                   "k", label=r"ASASSN14li $\Delta t = $60 d", zorder=2)
    ax[1].set_xlim(wmin, wmax)
    ax[1].set_ylim(2e-16, 3e-12)
    ax[1] = plot_line_id(ax[1])
    ax[1].legend(loc="lower right", fontsize=13)
    ax[1].text(0.05, 0.05, "Outside wind", ha="left", va="center", rotation="horizontal", fontsize=15,
                  transform=ax[1].transAxes)

    fig.text(0.5, 0.02, r"Rest Wavelength $\lambda$ [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=15)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)
    plt.savefig("polar_clumped_wind.png")
    plt.show()

    return


def main(argc: int, argv: List[str]):
    """
    Main function of the script.

    Parameters
    ----------
    argc: int
        The number of command line arguments provided
    argv: List[str]
        The command line arguments provided
    """

    plot_against_data("models/clump/1e-1/cv/solar/tde_cv.spec",
                      "models/disc/cv/tde_cv.spec",
                      ["60", ["10", "75"]], "clumped_wind")

    return


if __name__ == "__main__":
    main(len(argv), argv)
