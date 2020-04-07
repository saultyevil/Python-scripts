#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from sys import exit
from platform import system
from sys import argv
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
import numpy as np
from PyPython import SpectrumUtils, WindUtils
from PyPython import PythonUtils as Utils

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

SMOOTH = 5
VERBOSITY = False
SKIP = -1

LINES = [
    ["O VI", 1032],
    ["P V", 0],
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
    ["Mg II", 2798],
    ["He II", 4686],
    [r"H$_{\beta}$", 4861],
    [r"H$_{\alpha}$", 6564],
]


def plot_line_id(ax: plt.Axes) -> plt.Axes:
    """
    Plot labels and vertical lines to indicate important atomic transitions.
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


def model_comparison(directories: List[str], line_colours: List[str], wmin: float, wmax: float, ylims: List[float],
                     labels: List[str], continuum_spectrum_directory: str = None, extra_file_name: str = "")\
        -> Union[None, Tuple[plt.Figure, plt.Axes]]:
    """
    Plot a 1 x 3 grid of model comparisons.
    """

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde_optical/grid/round1/"
    else:
        pdir = "/home/saultyevil/PySims/tde_optical/grid/round1/"

    modelspecs = []
    for i in range(len(directories)):
        directories[i] = pdir + directories[i] + "/tde_uv.spec"
        try:
            modelspecs.append(SpectrumUtils.read_spec(directories[i]))
        except:
            modelspecs.append(SKIP)

    ncols = 3
    nrows = 1
    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 5), sharey="row")

    lstyle = [
        "-", "--", "-."
    ]

    incl = [
        "10", "60", "75"
    ]

    alpha = 0.75

    for i in range(ncols):
        for j in range(len(modelspecs)):
            if type(modelspecs[j]) == type(SKIP):
                continue
            wl = modelspecs[j]["Lambda"].values.astype(float)
            try:
                fl = modelspecs[j][incl[i]].values.astype(float)
            except KeyError:
                continue
            if i == 0:
                ax[i].semilogy(wl, SpectrumUtils.smooth(fl, SMOOTH), lstyle[j], label=labels[j], linewidth=3,
                               color=line_colours[j], alpha=alpha)
            else:
                ax[i].semilogy(wl, SpectrumUtils.smooth(fl, SMOOTH), lstyle[j], linewidth=3, color=line_colours[j],
                               alpha=alpha)
        ax[i].set_xlim(wmin, wmax)
        ax[i].set_ylim(ylims)
        ax[i] = plot_line_id(ax[i])

    ax[0].legend(loc="lower left", fontsize=13)
    fig.text(0.5, 0.02, r"Wavelength $\lambda$ [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=15)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0.1)

    comment = [":\n outside wind", ":\n in-wind", ":\n outside wind"]
    for i in range(ncols):
        tstr = r"$i = $" + incl[i] + r"$^{\circ}$" + comment[i]
        ax[i].text(0.83, 0.90, tstr, ha="center", va="center", rotation="horizontal", fontsize=13,
                   transform=ax[i].transAxes)

    fname = "comparison"
    if extra_file_name:
        fname += extra_file_name
    fname += ".png"

    plt.savefig(fname)
    plt.close()

    return fig, ax


def main(argc: int, argv: List[str]) -> None:
    """
    Main function of the script.

    Parameters
    ----------
    argc: int
        The number of command line arguments provided
    argv: List[str]
        The command line arguments provided
    """

    optical_wmin = 3900
    optical_wmax = 8100

    uv_wmin = 900
    uv_wmax = 2400

    mbh_grid = [
        "Mbh/1.0000e+06_temp_spec",
        "Mbh/1.0000e+07",
        "Mbh/1.0000e+08"
    ]

    mbh_labels = [
        r"M$_{BH}$ = 10$^6$ M$_{\odot}$",
        r"M$_{BH}$ = 10$^7$ M$_{\odot}$",
        r"M$_{BH}$ = 10$^8$ M$_{\odot}$",
    ]

    rmin_grid = [
        "Rmin/5.0000e+00",
        "Rmin/1.0000e+01",
        "Rmin/1.5000e+01"
    ]

    rmin_labels = [
        r"R$_{min}$ = 5 R$_{ISCO}$",
        r"R$_{min}$ = 10 R$_{ISCO}$",
        r"R$_{min}$ = 15 R$_{ISCO}$",
    ]

    vinf_grid = [
        "Vinf/1.0000e-01",
        "Vinf/5.0000e-01",
        "Vinf/8.0000e-01"
    ]

    vinf_labels = [
        r"V$_{\infty}$ = 0.1 V$_{esc}$",
        r"V$_{\infty}$ = 0.5 V$_{esc}$",
        r"V$_{\infty}$ = 0.8 V$_{esc}$"
    ]

    colors = ["C0", "C1", "C2"]

    # Optical figures

    fig, ax = model_comparison(mbh_grid.copy(), colors, optical_wmin, optical_wmax, [7e-6, 3e-2], mbh_labels,
                               extra_file_name="Mbh_optical")
    fig, ax = model_comparison(rmin_grid.copy(), colors, optical_wmin, optical_wmax, [5e-5, 3e-2], rmin_labels,
                               extra_file_name="Rmin_optical")
    fig, ax = model_comparison(vinf_grid.copy(), colors, optical_wmin, optical_wmax, [5e-5, 3e-2], vinf_labels,
                               extra_file_name="Vinf_optical")

    # UV figures

    fig, ax = model_comparison(mbh_grid.copy(), colors, uv_wmin, uv_wmax, [1e-6, 3], mbh_labels,
                               extra_file_name="Mbh_uv")
    fig, ax = model_comparison(rmin_grid.copy(), colors, uv_wmin, uv_wmax, [1e-3, 3], rmin_labels,
                               extra_file_name="Rmin_uv")
    fig, ax = model_comparison(vinf_grid.copy(), colors, uv_wmin, uv_wmax, [1e-3, 3], vinf_labels,
                               extra_file_name="Vinf_uv")

    return


if __name__ == "__main__":
    main(len(argv), argv)
