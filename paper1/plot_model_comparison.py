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
from typing import List, Tuple, Union
import numpy as np
from PyPython import SpectrumUtils, Utils, WindUtils

SMOOTH = 10
VERBOSITY = False

LINES = [
    # ["O VI", 0],
    ["P V", 1118],
    # ["N III]", 0],
    [r"Ly$\alpha$/N V", 1216],
    ["", 1240],
    ["O V/Si IV", 1371],
    ["", 1400],
    # ["N IV]", 1489],
    ["C IV",  1549],
    ["He II", 1640],
    # ["O III]", 0],
    # ["N III]", 1750],
    # ["C III]", 1908],
    # ["Fe II", 0],
    # ["Fe II / CII]", 0],
    # ["Fe II", 0],
    # ["Fe II", 0],
    # ["Mg II", 2798]
]


def plot_line_id(ax: plt.Axes, labels: bool) -> plt.Axes:
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

    for i in range(nlines):
        lab = LINES[i][0]
        x = LINES[i][1]
        if x < xlims[0]:
            continue
        elif x > xlims[1]:
            continue

        ax.axvline(x, ymax=1, linestyle="--", linewidth=0.45, color="k")

        if labels is False:
            continue

        x = x - 25
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        ax.text(xnorm, 0.92, lab, ha="center", va="center", rotation="vertical", transform=ax.transAxes, fontsize=10)

    return ax


def sightline_coords(x: np.ndarray, theta: float):
    """
    Return the vertical coordinates for a sightline given the x coordinates
    and the inclination of the sightline.

    Parameters
    ----------
    x: np.ndarray[float]
        The x-coordinates of the sightline
    theta: float
        The opening angle of the sightline
    
    Returns
    -------
    z: np.ndarray[float]
        The z-coordinates of the sightline
    """

    return x * np.tan(np.pi / 2 - theta)


def model_comparison(direcs: List[str], extrafname: str = "", wmin: float = 800, wmax: float = 3600, label: str = None,
                     return_figure: bool = False, figure: Tuple[plt.Figure, plt.Axes] = None) -> Union[None, Tuple[plt.Figure, plt.Axes]]:
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
    label: str, optional
        A label to give to the data
    return_figure: bool, optional
        If True, the plt.Figure and plt.Axes objects will be returned
    figure: Tuple[plt.Figure, plt.Axes], optional
        If a plt.Figure and plt.Axes object are provided, then these will be
        used instead of generating a new one.

    Returns
    -------
    fig: plt.Figure
        If return_figure is True, then the plt.Figure and plt.Axes object which
        are the figure objects are returned instead of the figure being created.
        This particular return is the plt.Figure object.
    ax: plt.Axes
        If return_figure is True, then the plt.Figure and plt.Axes object which
        are the figure objects are returned instead of the figure being created.
        This particular return is the plt.Axes object.
    """

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"

    modelspecs = []
    for i in range(len(direcs)):
        direcs[i] = pdir + direcs[i]
        modelspecs.append(SpectrumUtils.read_spec(direcs[i]))

    ncols = 3
    nrows = 3

    if figure:
        fig, ax = figure[0], figure[1]
        lstyle = "--"
    else:
        fig, ax = plt.subplots(nrows, ncols, figsize=(9.5, 11), sharex="col", sharey="row")
        lstyle = "-"

    #       CV    AGN   Spherical
    incl = ["30", "30", "30",
            "60", "60", "60",
            "75", "75", "75"]
    iidx = 0

    ylims = [(2e-3, 0.5), (0.8e-3, 0.2), (6e-4, 0.1)]

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
            ax[i, j].semilogy(wl, SpectrumUtils.smooth_spectrum(fl, SMOOTH), lstyle, label=label)
            ax[i, j].set_xlim(wmin, wmax)
            ax[i, j].set_ylim(ylims[i])
            # if i != 0:
            #     labels = False
            # else:
            #     labels = True
            # ax[i, j] = plot_line_id(ax[i, j], labels)

    mnames = ["Conical Wind Model", "Equatorial Wind Model", "Spherical Wind Model"]
    for i in range(ncols):
        ax[0, i].text(0.5, 1.1, mnames[i], va="center", ha="center", rotation="horizontal", fontsize=13,
                      transform=ax[0, i].transAxes)
        ax[-1, i].tick_params(axis="x", rotation=30, labelsize=12)

    incls = ["30", "60", "75"]
    for i in range(nrows):
        ax[i, 0].tick_params(axis="y", labelsize=12)
        tstr = r"$i = $" + incls[i] + r"$^{\circ}$"
        ax[i, -1].text(0.85, 0.93, tstr, ha="center", va="center", rotation="horizontal", fontsize=12,
                    transform=ax[i, j].transAxes)

    if figure:
        ax[-1, 0].legend(loc="lower center")

    fig.text(0.5, 0.02, r"Rest Wavelength [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=15)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    fname = "tde_model_comparison"
    if extrafname:
        fname += extrafname
    fname += ".png"

    if not return_figure:
        plt.savefig(fname)
        plt.show()
    else:
        return fig, ax

    return


def wind_geos(direcs) -> Union[plt.Figure, list]:
    """
    Create a plot of the three different wind geometries.

    Parameters
    ----------
    direcs: List[str]
        The directories containing the root.ep.complete files.

    Returns
    -------
    fig: plt.Figure
        The figure object for the figure
    ax: list
        A list containing plt.Axes for each panel of the figure
    """

    nplots = 6
    incls = ["30", "60", "75"]
    lstyle = ["k--", "k-.", "k:"]

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"

    modelwind = []
    coords = ["rectilinear", "rectilinear", "polar"]
    for i in range(len(direcs)):
        direcs[i] = pdir + direcs[i]
        root, wd = Utils.split_root_directory(direcs[i])
        wx, wz, ww = WindUtils.extract_wind_var(root, "rho", "wind", wd, coords[i])
        modelwind.append([wx, wz, np.log10(ww)])

    fig = plt.figure(figsize=(15, 8))
    ax1 = plt.subplot(2, 3, 1)
    ax2 = plt.subplot(2, 3, 2)
    ax3 = plt.subplot(2, 3, 3, projection="polar")
    ax4 = plt.subplot(2, 3, 4)
    ax5 = plt.subplot(2, 3, 5)
    # ax6 = plt.subplot(2, 3, 6, projection="polar")
    ax = [ax1, ax2, ax3, ax4, ax5]

    for i in range(nplots):
        # Log-Log plots
        if i < 3:
            wx, wz, ww = modelwind[i]
            if coords[i] != "polar":
                im = ax[i].pcolor(np.log10(wx), np.log10(wz), ww)
                ax[i].set_xlim(np.log10(wx[1, 1]), np.log10(wx[-1, -1]))
                ax[i].set_ylim(np.log10(wz[1, 1]), np.log10(wz[-1, -1]))
                ax[i].set_xlabel("Log[x]")
                ax[i].set_ylabel("Log[z]")
            else:
                ax[i].set_theta_zero_location("N")
                ax[i].set_theta_direction(-1)
                im = ax[i].pcolor(np.deg2rad(wz), np.log10(wx), ww)
                ax[i].set_thetamin(0)
                ax[i].set_thetamax(90)
                ax[i].set_rlim(np.log10(wx[1][0]), np.log10(wx[-2][0]))
                ax[i].set_xticklabels([])
                ax[i].set_ylabel("Log[R]")
                ax[i].set_rlabel_position(90)

            plt.colorbar(im, ax=ax[i])

            for j in range(len(incls)):
                if coords[i] != "polar":
                    xsight = np.linspace(0, np.max(wx), 1e5)
                    zsight = sightline_coords(xsight, np.deg2rad(float(incls[j])))
                    ax[i].plot(np.log10(xsight), np.log10(zsight), lstyle[j],
                               label=incls[j] + r"$^{\circ}$ line of sight")
                else:
                    xsight = np.linspace(0, 1e17, 1e5)
                    zsight = sightline_coords(xsight, np.deg2rad(90 - float(incls[j])))
                    rsight = np.sqrt(xsight ** 2 + zsight ** 2)
                    thetasight = np.arctan(zsight / xsight)
                    ax[i].plot(thetasight, np.log10(rsight), lstyle[j], label=incls[j] + r"$^{\circ}$ line of sight")
        # Lin-Lin plots
        else:
            ii = i - 3
            wx, wz, ww = modelwind[ii]
            if coords[ii] != "polar":
                im = ax[i].pcolor(wx, wz, ww)
                ax[i].set_xlim(wx[1, 1], 5.5e17)
                ax[i].set_ylim(wz[1, 1], 5.5e17)
                ax[i].set_xlabel("x")
                ax[i].set_ylabel("z")
            else:
                continue
                # ax[i].set_theta_zero_location("N")
                # ax[i].set_theta_direction(-1)
                # im = ax[i].pcolor(np.deg2rad(wz), wx, ww)
                # ax[i].set_thetamin(0)
                # ax[i].set_thetamax(90)
                # ax[i].set_xticklabels([])
                # ax[i].set_ylabel("R")
                # ax[i].set_rlabel_position(90)
                # ax[i].set_rlim(1e13, 1e17)

            plt.colorbar(im, ax=ax[i])

            for j in range(len(incls)):
                if coords[ii] != "polar":
                    xsight = np.linspace(0, np.max(wx), 1e5)
                    zsight = sightline_coords(xsight, np.deg2rad(float(incls[j])))
                    ax[i].plot(xsight, zsight, lstyle[j], label=incls[j] + r"$^{\circ}$ line of sight")
                # else:
                #     xsight = np.linspace(0, 1e17, 1e5)
                #     zsight = sightline_coords(xsight, np.deg2rad(90 - float(incls[j])))
                #     rsight = np.sqrt(xsight ** 2 + zsight ** 2)
                #     thetasight = np.arctan(zsight / xsight)
                #     ax[i].plot(thetasight, rsight, lstyle[j], label=incls[j] + r"$^{\circ}$ line of sight")

    ax[0].legend()

    mnames = ["Conical Wind Model", "Equatorial Wind Model", "Spherical Wind Model"]
    for i in range(nplots - 3):
        ax[i].text(0.5, 1.1, mnames[i], va="center", ha="center", rotation="horizontal", fontsize=13,
                   transform=ax[i].transAxes)

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    
    plt.savefig("tde_model_comparison_wind_geo.png")
    plt.show()

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
            print("Could not convert wmax = {} intoig, ax,  float".format(argv[2]))
            exit(1)
    elif argc != 1:
        print(__doc__)
        exit(1)

    solar_smooth = ["paper_models/smooth/cv/solar/tde_cv.spec",
                    "paper_models/smooth/agn/solar/tde_agn.spec",
                    "paper_models/smooth/spherical/solar/tde_spherical.spec"]

    cno_smooth = ["paper_models/smooth/cv/cno/tde_cv.spec",
                  "paper_models/smooth/agn/cno/tde_agn.spec",
                  "paper_models/smooth/spherical/cno/tde_spherical.spec"]

    fig, ax = model_comparison(solar_smooth.copy(), "_solar_cno_smooth_abundances", wmin, wmax, label="Solar Abundance", return_figure=True)
    model_comparison(cno_smooth.copy(), "_solar_cno_smooth_abundances", wmin, wmax, return_figure=False, figure=(fig, ax),
                     label="CNO Processed Abundance\nHe = 2 x Solar\nC = 0.5 x Solar\nN = 7 x Solar")

    solar_clump = ["paper_models/clump/1e-1/cv/solar/tde_cv.spec",
                    "paper_models/clump/1e-1/agn/solar/tde_agn.spec",
                    "paper_models/clump/1e-1/spherical/solar/tde_spherical.spec"]

    fig, ax = model_comparison(solar_smooth.copy(), "_solar_clump", wmin, wmax, label="f = 1", return_figure=True)
    model_comparison(solar_clump.copy(), "_solar_clump", wmin, wmax, label="f = 0.1", return_figure=False,
                     figure=(fig, ax))

    wind_geos(solar_smooth.copy())

    return


if __name__ == "__main__":
    main(len(argv), argv)
