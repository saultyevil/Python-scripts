#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

"""
Create a comparison figure for three TDE geometries.

Usage: [python] plot_model_comparison.py [wmin] [wmax]

    wmin: optional, the smallest wavelength to show
    wmax: optional, the largest wavelength to show

Note that wmin and wmax have to be both provided at the same time.
"""

import sys
from sys import exit
from platform import system
from sys import argv
from matplotlib import pyplot as plt
from typing import List, Tuple

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

import py_plot_util as ppu

SMOOTH = 5
VERBOSE = False

LINES = [
    ["O VI", 0],
    ["P V", 1118],
    ["N III]", 0],
    [r"Ly$\alpha$", 1216],
    ["N V", 1240],
    ["O I", 0],
    ["Si IV", 1400],
    ["N IV]", 1489],
    ["C IV",  1549],
    ["He II", 0],
    ["O III]", 0],
    ["N III]", 1750],
    ["C III]", 1908],
    ["Fe II", 0],
    ["Fe II / CII]", 0],
    ["Fe II", 0],
    ["Fe II", 0],
    ["Mg II", 2798]
]


def plot_line_id(ax: plt.Axes, xlims: Tuple[float, float]) -> plt.Axes:
    """
    Plot labels and vertical lines to indicate important atomic transitions.

    Parameters
    ----------
    ax: plt.Axes
        The Axes object to add line ID labels to.
    xlims: Tuple[float, float]
        The x coordinate limits for the figure. Line labels outside of this
        range will not be added to the plot.

    Returns
    -------
    ax: plt.Axes
        The input Axes object with additional line IDs.
    """

    xlims = ax.get_xlim()
    nlines = len(LINES)

    iloc = 0
    for i in range(nlines):
        lab = LINES[i][0]
        x = LINES[i][1]
        if x < xlims[0]:
            continue
        elif x > xlims[1]:
            continue
        ax.axvline(x, ymin=0.15, ymax=0.85, linestyle="--", linewidth=0.45, color="k")
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        if iloc % 2 == 0:
            yloc = 0.95
        else:
            yloc = 0.05
        iloc += 1
        ax.text(xnorm, yloc, lab, ha="center", va="center", rotation="vertical", transform=ax.transAxes,
                fontsize=9)

    return ax


def model_comparison(direcs: List[str], extrafname: str = "", wmin: float = 800, wmax: float = 3600) -> None:
    """
    Plot a 3 x 3 grid of model comparisons for the three TDE geometries, AGN,
    CV and spherical. Each model will have a low, a medium and high inclination
    component to the plot - these will be hardcoded as I don't expect to use
    this function ever again.

    Parameters
    ----------
    direcs: List[str]
        The directories containing the spectra to plot
    extrafname: str, optional
        Append an extra string to the end of the output filename
    wmin: float, optional
        The smallest wavelength to show on the figure
    wmax: float, optional
        The largest wavelength to show on the figure
    legend: bool, optional
        If True, include a legend for each panel in the figure
    """

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"

    modelspecs = []
    for i in range(len(direcs)):
        direcs[i] = pdir + direcs[i]
        modelspecs.append(ppu.read_spec_file(direcs[i], pandas_table=True))
        print("Loaded spec: ", direcs[i])

    ncols = 3
    nrows = 3

    fig, ax = plt.subplots(nrows, ncols, figsize=(9.5, 11), sharex="col", sharey="row")
    #       CV    AGN   Spherical
    incl = ["20", "70", "30",
            "62", "75", "60",
            "75", "85", "80"]
    iidx = 0

    ylims = [(5e-3, 0.13), (1e-3, 0.07), (1e-3, 0.06)]
    xlims = (wmin, wmax)

    for i in range(nrows):
        for j in range(ncols):
            wl = modelspecs[j]["Lambda"].values.astype(float)
            try:
                fl = modelspecs[j][incl[iidx]].values.astype(float)
                iidx += 1
            except KeyError:
                print("Inclination {} w/ iidx {} not found for model with x,y indices {},{}: {}".format(incl[iidx], iidx, i, j, direcs[j]))
                iidx += 1
                continue
            wminidx = 0
            wmaxidx = 0
            for k in range(len(wl)):
                if wl[k] > wmin:
                    wminidx = k - 1
                if wl[k] > wmax:
                    wmaxidx = k
            wl = wl[wmaxidx:wminidx]
            fl = fl[wmaxidx:wminidx]

            # TODO: line ids look a bit fucked up
            # ax[i, j].set_xlim(wmin - 40, wmax + 40)
            # ax[i, j] = plot_line_id(ax[i, j], xlims)

            ax[i, j].semilogy(wl, ppu.smooth_1d_array(fl, SMOOTH, VERBOSE))
            ax[i, j].set_ylim(ylims[i])
            tstr = r"$i = $" + incl[iidx - 1] + r"$^{\circ}$"
            ax[i, j].text(0.85, 0.93, tstr, ha="center", va="center", rotation="horizontal", fontsize=12,
                          transform=ax[i, j].transAxes)

    fig.text(0.5, 0.02, r"Rest Wavelength [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=15)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    fname = "tde_model_comparison"
    if extrafname:
        fname += extrafname
    fname += ".png"
    plt.savefig(fname)
    plt.show()

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

    wmin = 900
    wmax = 2100

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

    solar = ["cv_macro/zorig/tde_cv.spec", "agn_macro/zorig/tde_agn.spec", "spherical_macro/zorig/tde_spherical.spec"]
    cno = ["cno_processed/cv_macro_cno/zorig/tde_cv.spec", "cno_processed/agn_macro_cno/zorig/tde_agn.spec",
           "cno_processed/spherical_macro_cno/zorig/tde_spherical_cno.spec"]

    model_comparison(solar.copy(), "_solar", wmin, wmax)
    # model_comparison(cno.copy(), "_cno", wmin, wmax, legend=False)
    # model_comparison(solar.copy() + cno.copy(), "_solar_cno", wmin, wmax, legend=True)

    return


if __name__ == "__main__":
    main(len(argv), argv)
