#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

"""
Plot a comparison between the best tde model against the data for a BAL and
BEL TDE.
"""

import sys
from sys import exit, argv
from platform import system
from matplotlib import pyplot as plt
from typing import List

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

import tde_util as tu
import py_plot_util as ppu
from consts import *

SMOOTH = 10

LINES = [
    ["O VI", 0],
    ["P V", 1118],
    ["N III]", 0],
    [r"Ly$\alpha$/N V", 1216],
    ["", 1240],
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
    ["Mg II", 2798]
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

        ax.axvline(x, ymax=1, linestyle="--", linewidth=0.45, color="k")
        x = x - 25
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        ax.text(xnorm, 0.92, lab, ha="center", va="center", rotation="vertical", transform=ax.transAxes, fontsize=13)

    return ax


def plot_against_data(dfname: List[str], inclinations: List[str]):
    """

    Parameters
    ----------
    dfname: List[str]
        The file names of the models to plot - typically the solar and cno
        abundance models are given.
    inclinations: List[str]
        The inclinations to plot - this should be the inclination for the BAL
        followed by the BEL model
    """

    wmin = 1000
    wmax = 3000

    iptf = tu.iptf15af_spec(SMOOTH)
    assa = tu.asassn14li_spec(SMOOTH)

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"
    for i in range(len(dfname)):
        dfname[i] = pdir + dfname[i]
    mspec = ppu.read_spec(dfname[0], numpy=True)
    cnospec = ppu.read_spec(dfname[1], numpy=True)

    fig, ax = plt.subplots(2, 1, figsize=(9.5, 11), sharex="col")

    inclination = inclinations[0]
    mspec_wl = mspec["Lambda"].values.astype(float)
    mspec_fl = mspec[inclination].values.astype(float)
    mspec_fl *= (100 * PARSEC) ** 2 / (350 * 1e6 * PARSEC) ** 2
    cnospec_wl = cnospec["Lambda"].values.astype(float)
    cnospec_fl = cnospec[inclination].values.astype(float)
    cnospec_fl *= (100 * PARSEC) ** 2 / (350 * 1e6 * PARSEC) ** 2

    ax[0].semilogy(iptf[:, 0] / (0.07897 + 1), ppu.smooth(iptf[:, 1], SMOOTH),
                   label=r"iPTF15af $\Delta t = $52 d")
    ax[0].semilogy(mspec_wl, ppu.smooth(mspec_fl, SMOOTH), label="Solar Abundance")
    ax[0].semilogy(cnospec_wl, ppu.smooth(cnospec_fl, SMOOTH), label="CNO Abundance")
    ax[0].set_xlim(wmin, wmax)
    ax[0].set_ylim(3e-17, 9e-15)
    ax[0] = plot_line_id(ax[0])
    ax[0].text(0.68, 0.06, r"$i = $" + inclination + "$^{\circ}$", transform=ax[0].transAxes, fontsize=15)
    ax[0].legend(loc="lower right", fontsize=9)

    inclination = inclinations[1]
    mspec_wl = mspec["Lambda"].values.astype(float)
    mspec_fl = mspec[inclination].values.astype(float)
    mspec_fl *= (100 * PARSEC) ** 2 / (55 * 1e6 * PARSEC) ** 2
    cnospec_wl = cnospec["Lambda"].values.astype(float)
    cnospec_fl = cnospec[inclination].values.astype(float)
    cnospec_fl *= (100 * PARSEC) ** 2 / (55 * 1e6 * PARSEC) ** 2

    ax[1].semilogy(assa[:, 0] / (0.02058 + 1), ppu.smooth(assa[:, 1], SMOOTH),
                   label=r"ASASSN14li $\Delta t = $60 d")
    ax[1].semilogy(mspec_wl, ppu.smooth(mspec_fl, SMOOTH), label="Solar Abundance")
    ax[1].semilogy(cnospec_wl, ppu.smooth(cnospec_fl, SMOOTH), label="CNO Abundance")
    ax[1].set_xlim(wmin, wmax)
    ax[1].set_ylim(9e-16, 8e-14)
    ax[1] = plot_line_id(ax[1])
    ax[1].text(0.67, 0.06, r"$i = $" + inclination + "$^{\circ}$", transform=ax[1].transAxes, fontsize=15)
    ax[1].legend(loc="lower right", fontsize=9)

    fig.text(0.5, 0.02, r"Rest Wavelength [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=13)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=13)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    plt.savefig("best_model_comparison.png")
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

    plot_against_data(["nv_fitting/zzzgood/tde_agn_5rwind_mdot_ratio_5e-2/tde_agn.spec",
                       "running/tde_agn_5rwind_mdot_ratio_5e-2_cno/tde_agn.spec"],
                      ["80", "85"])

    return


if __name__ == "__main__":
    main(len(argv), argv)
