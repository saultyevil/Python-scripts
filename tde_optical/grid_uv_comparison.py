#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

"""
Create a comparison figure for three TDE geometries.

Usage: [python] plot_model_comparison.py [wmin] [wmax]

    wmin: optional, the smallest wavelength to show
    wmax: optional, the largest wavelength to show

Note that wmin and wmax have to be both provided at the same time.
"""

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


def model_comparison(direcs: List[str], colour, disk_specs: List[str] = None, extrafname: str = "", wmin: float = 4000,
                     wmax: float = 8000, label: str = None, return_figure: bool = False,
                     figure: Tuple[plt.Figure, plt.Axes] = None) -> Union[None, Tuple[plt.Figure, plt.Axes]]:
    """
    Plot a 3 x 3 grid of model comparisons for the three TDE geometries, AGN,
    CV and spherical. Each model will have a low, a medium and high inclination
    component to the plot - these will be hardcoded as I don't expect to use
    this function ever again.
    """

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde_optical/grid/round1/"
    else:
        pdir = "/home/saultyevil/PySims/tde_optical/grid/round1/"

    modelspecs = []
    for i in range(len(direcs)):
        direcs[i] = pdir + direcs[i] + "/tde_uv.spec"
        try:
            modelspecs.append(SpectrumUtils.read_spec(direcs[i]))
        except:
            modelspecs.append(SKIP)

    ncols = 3
    nrows = 3

    if figure:
        fig, ax = figure[0], figure[1]
        if return_figure:
            lstyle = "--"
        else:
            lstyle = "-."
    else:
        fig, ax = plt.subplots(nrows, ncols, figsize=(13.5, 15.5), sharex="col", sharey="row")
        lstyle = "-"

    incl = [
        "10", "60", "75"
    ]

    iidx = 0
    ylims = [(1e-3, 1), (1e-3, 1), (1e-3, 1)]
    alpha = 1.0

    for i in range(nrows):
        for j in range(ncols):
            if type(modelspecs[j]) == type(SKIP):
                continue
            wl = modelspecs[j]["Lambda"].values.astype(float)
            try:
                fl = modelspecs[j][incl[i]].values.astype(float)
                iidx += 1
            except KeyError:
                iidx += 1
                continue
            ax[i, j].semilogy(wl, SpectrumUtils.smooth(fl, SMOOTH), lstyle, label=label, linewidth=3, color=colour, alpha=alpha)
            ax[i, j].set_xlim(wmin, wmax)
            ax[i, j].set_ylim(ylims[i])
            ax[i, j] = plot_line_id(ax[i, j])

    # if disk_specs and not return_figure:
    #     dspecs = []
    #     for i in range(len(disk_specs)):
    #         dspecs.append(SpectrumUtils.read_spec(pdir + disk_specs[i]))
    #     iidx = 0
    #     for i in range(nrows):
    #         for j in range(ncols):
    #             dwl = dspecs[j]["Lambda"].values.astype(float)
    #             try:
    #                 dfl = dspecs[j][incl[iidx]].values.astype(float)
    #             except KeyError:
    #                 print("Unable to find inclination {} in disk spectrum {}".format(incl[iidx], disk_specs[j]))
    #                 iidx += 1
    #                 continue
    #             ax[i, j].semilogy(dwl, SpectrumUtils.smooth_spectrum(dfl, 15), "k-.", alpha=0.8, zorder=2,
    #                               label="Disc Model")
    #             plot_line_id(ax[i, j])
    #             iidx += 1

    if figure:
        ax[-1, -1].legend(loc="lower left", fontsize=13)

    fig.text(0.5, 0.02, r"Wavelength $\lambda$ [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=17)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=17)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    mnames = ["Mbh", "Rmin", "Vinf"]
    for i in range(ncols):
        ax[0, i].text(0.5, 1.1, mnames[i], va="center", ha="center", rotation="horizontal", fontsize=15,
                      transform=ax[0, i].transAxes)

    incls = incl
    comment = [" : outside wind", " : in-wind", " : outside wind"]
    for i in range(nrows):
        tstr = r"$i = $" + incls[i] + r"$^{\circ}$" + comment[i]
        ax[i, 0].text(0.5, 0.10, tstr, ha="center", va="center", rotation="horizontal", fontsize=15,
                    transform=ax[i, 0].transAxes)

    fname = "comparison"
    if extrafname:
        fname += extrafname
    fname += ".png"

    if not return_figure:
        plt.savefig(fname)
        plt.show()
    else:
        return fig, ax

    return


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

    wmin = 1000
    wmax = 3000

    if argc == 3:
        try:
            wmin = float(argv[1])
        except ValueError:
            print("Could not convert wmin = {} into float".format(argv[1]))
            exit(1)
        try:
            wmax = float(argv[2])
        except ValueError:
            print("Could not convert wmax = {} into float".format(argv[2]))
            exit(1)
    elif argc != 1:
        print(__doc__)
        exit(1)


    plot1 = [
        "",
        "Rmin/1.0000e+01",
        "Vinf/1.0000e-01"
    ]

    plot2 = [
        "Mbh/1.0000e+07",
        "Rmin/1.5000e+01",
        "Vinf/5.0000e-01"
    ]

    plot3 = [
        "Mbh/1.0000e+08",
        "Rmin/5.0000e+00",
        "Vinf/8.0000e-01"
    ]

    disk = None

    print("1")
    fig, ax = model_comparison(plot1.copy(), "C2", disk, "_uv_aaa", wmin, wmax, label="Mbh=N/A,Rmin=10,Vinf=0.1", return_figure=True)
    print("2")
    fig, ax = model_comparison(plot2.copy(), "C1", disk, "_uv_aaa", wmin, wmax, return_figure=True, figure=(fig, ax), label="Mbh=1e7,Rmin=15,Vinf=0.5")
    print("3")
    model_comparison(plot3.copy(), "C0", disk, "_uv_aaa", wmin, wmax, return_figure=False, figure=(fig, ax), label="Mbh=1e8,Rmin=5,Vinf=0.8")

    return


if __name__ == "__main__":
    main(len(argv), argv)
