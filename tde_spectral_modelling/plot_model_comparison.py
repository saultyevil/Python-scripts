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
from PyPython import SpectrumUtils, Utils, WindUtils

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


def model_comparison(direcs: List[str], colour, disk_specs: List[str] = None, extrafname: str = "", wmin: float = 800,
                     wmax: float = 3600, label: str = None, return_figure: bool = False,
                     figure: Tuple[plt.Figure, plt.Axes] = None) -> Union[None, Tuple[plt.Figure, plt.Axes]]:
    """
    Plot a 3 x 3 grid of model comparisons for the three TDE geometries, AGN,
    CV and spherical. Each model will have a low, a medium and high inclination
    component to the plot - these will be hardcoded as I don't expect to use
    this function ever again.

    Parameters
    ----------
    direcs: List[str]
        The directories containing the spectra to plot
    disk_specs: List[str]
        The directories containing disk spectrum to plot
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
        try:
            modelspecs.append(SpectrumUtils.read_spec(direcs[i]))
        except FileNotFoundError:
            modelspecs.append(SKIP)

    ncols = 2
    nrows = 3

    if figure:
        fig, ax = figure[0], figure[1]
        lstyle = "--"
    else:
        fig, ax = plt.subplots(nrows, ncols, figsize=(13.5, 15.5), sharex="col", sharey="row")
        lstyle = "-"

    incl = ["10", "10",
            "60", "60",
            "75", "75"]

    iidx = 0
    ylims = [(1e-2, 0.4), (2e-4, 0.35), (2e-4, 0.35)]
    alpha = 1.0

    for i in range(nrows):
        for j in range(ncols):
            if type(modelspecs[j]) == type(SKIP):
                continue
            wl = modelspecs[j]["Lambda"].values.astype(float)
            try:
                fl = modelspecs[j][incl[iidx]].values.astype(float)
                iidx += 1
            except KeyError:
                print("Inclination {} w/ iidx {} not found for model with x,y indices {},{}: {}".format(incl[iidx], iidx, i, j, direcs[j]))
                print("Current headings: ", modelspecs[j].columns.values)
                iidx += 1
                continue
            ax[i, j].semilogy(wl, SpectrumUtils.smooth_spectrum(fl, SMOOTH), lstyle, label=label, linewidth=3, zorder=3,
                              color=colour, alpha=alpha)
            ax[i, j].set_xlim(wmin, wmax)
            ax[i, j].set_ylim(ylims[i])

    if disk_specs and not return_figure:
        dspecs = []
        for i in range(len(disk_specs)):
            dspecs.append(SpectrumUtils.read_spec(pdir + disk_specs[i]))
        iidx = 0
        for i in range(nrows):
            for j in range(ncols):
                dwl = dspecs[j]["Lambda"].values.astype(float)
                try:
                    dfl = dspecs[j][incl[iidx]].values.astype(float)
                except KeyError:
                    print("Unable to find inclination {} in disk spectrum {}".format(incl[iidx], disk_specs[j]))
                    iidx += 1
                    continue
                ax[i, j].semilogy(dwl, SpectrumUtils.smooth_spectrum(dfl, 15), "k-.", alpha=0.8, zorder=2,
                                  label="Disc Model")
                plot_line_id(ax[i, j])
                iidx += 1

    mnames = ["Polar Wind Model", "Equatorial Wind Model", "Spherical Wind Model"]
    for i in range(ncols):
        ax[0, i].text(0.5, 1.1, mnames[i], va="center", ha="center", rotation="horizontal", fontsize=15,
                      transform=ax[0, i].transAxes)

    incls = ["10", "60", "75"]
    comment = [" : outside wind", " : in-wind", " : outside wind"]
    for i in range(nrows):
        tstr = r"$i = $" + incls[i] + r"$^{\circ}$" + comment[i]
        ax[i, 0].text(0.77, 0.10, tstr, ha="center", va="center", rotation="horizontal", fontsize=15,
                    transform=ax[i, 0].transAxes)
    comment = [" : outside wind", " : outside wind", " : in-wind"]
    for i in range(nrows):
        tstr = r"$i = $" + incls[i] + r"$^{\circ}$" + comment[i]
        ax[i, 1].text(0.77, 0.10, tstr, ha="center", va="center", rotation="horizontal", fontsize=15,
                    transform=ax[i, 1].transAxes)

    if figure:
        ax[-1, 0].legend(loc="lower left", fontsize=13)

    fig.text(0.5, 0.02, r"Wavelength $\lambda$ [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=17)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=17)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    fname = "comparison"
    if extrafname:
        fname += extrafname
    fname += ".png"

    if not return_figure:
        plt.savefig(fname)
        plt.close()
    else:
        return fig, ax

    return


import numpy.ma as ma


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

    nplots = 4
    incls = ["10", "60", "75"]
    lstyle = ["k-", "k--", "k-.", "k:"]

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde/"
    else:
        pdir = "/home/saultyevil/PySims/tde/"

    colormap = "t_e"

    modelwind = []
    coords = ["rectilinear", "rectilinear"]
    for i in range(len(coords)):
        direcs[i] = pdir + direcs[i]
        root, wd = Utils.split_root_directory(direcs[i])
        wx, wz, ww = WindUtils.extract_wind_var(root, colormap, "wind", wd, coords[i])
        for j in range(ww.shape[0]):
            for k in range(ww.shape[1]):
                if ww[j, k] < 100 and type(ww[j, k]) == np.float64:
                    ww[j, k] = 100
        modelwind.append([wx, wz, np.log10(ww)])

    fig = plt.figure(figsize=(14, 10))
    ax1 = plt.subplot(2, 2, 1)
    ax2 = plt.subplot(2, 2, 2)
    ax3 = plt.subplot(2, 2, 3)
    ax4 = plt.subplot(2, 2, 4)
    ax = [ax1, ax2, ax3, ax4]

    for i in range(nplots):
        # Log-Log plots
        if i < len(coords):
            wx, wz, ww = modelwind[i]
            im = ax[i].pcolor(np.log10(wx), np.log10(wz), ww)
            ax[i].set_xlim(np.log10(wx[1, 1]), np.log10(wx[-1, -1]))
            ax[i].set_ylim(np.log10(wz[1, 1]), np.log10(wz[-1, -1]))
            ax[i].set_xlabel("Log[x] [cm]", fontsize=15)
            ax[i].set_ylabel("Log[z] [cm]", fontsize=15)
            plt.colorbar(im, ax=ax[i])
            for j in range(len(incls)):
                xsight = np.linspace(0, np.max(wx), 1e5)
                zsight = sightline_coords(xsight, np.deg2rad(float(incls[j])))
                ax[i].plot(np.log10(xsight), np.log10(zsight), lstyle[j],
                           label=incls[j] + r"$^{\circ}$ sightline")
        # Lin-Lin plots
        else:
            ii = i - 3
            wx, wz, ww = modelwind[ii]
            im = ax[i].pcolor(wx / 1e17, wz / 1e17, ww)
            ax[i].set_xlim(wx[1, 1] / 1e17, 5.5e17 / 1e17)
            ax[i].set_ylim(wz[1, 1] / 1e17, 5.5e17 / 1e17)
            ax[i].set_xlabel("x ($10^{17}$) [cm]", fontsize=15)
            ax[i].set_ylabel("z ($10^{17}$) [cm]", fontsize=15)
            plt.colorbar(im, ax=ax[i])
            for j in range(len(incls)):
                xsight = np.linspace(0, np.max(wx), 1e5)
                zsight = sightline_coords(xsight, np.deg2rad(float(incls[j])))
                ax[i].plot(xsight, zsight, lstyle[j], label=incls[j] + r"$^{\circ}$ sightline")

    ax[0].legend(loc="lower right")

    mnames = ["Polar Wind Model", "Equatorial Wind Model"]
    for i in range(len(mnames)):
        ax[i].text(0.5, 1.1, mnames[i], va="center", ha="center", rotation="horizontal", fontsize=13,
                   transform=ax[i].transAxes)

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    plt.savefig("comparison_wind_geo_{}.png".format(colormap))
    plt.show(block=True)

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

    wmin = 950
    wmax = 2050

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

    solar_smooth = ["models/smooth/cv/solar/tde_cv.spec",
                    "models/smooth/agn/solar/tde_agn.spec"]
                    # "ztodo_iridis/models/smooth/spherical/solar/tde_spherical.spec"]

    cno_smooth = ["models/smooth/cv/cno/tde_cv.spec",
                  "models/smooth/agn/cno/tde_agn.spec"]
                  # "ztodo_iridis/models/smooth/spherical/cno/tde_spherical.spec"]

    solar_clump = ["models/clump/1e-1/cv/solar/tde_cv.spec",
                   "models/clump/1e-1/agn/solar/tde_agn.spec"]
                   # "ztodo_iridis/models/clump/1e-1/spherical/solar/tde_spherical.spec"]

    cno_clump = ["models/clump/1e-1/cv/cno/tde_cv.spec",
                 "models/clump/1e-1/agn/cno/tde_agn.spec"]
                 # "ztodo_iridis/models/clump/1e-1/spherical/cno/tde_spherical.spec"]

    disk = ["models/disc/cv/tde_cv.spec",
            "models/disc/agn/tde_agn.spec"]
            # "ztodo_iridis/models/disc/spherical/tde_spherical.spec"]

    # fig, ax = model_comparison(solar_smooth.copy(), "C2", disk, "_solar_cno_smooth_abundances", wmin, wmax,
    #                            label="Solar Abundances", return_figure=True)
    # model_comparison(cno_smooth.copy(), "C1", disk, "_solar_cno_smooth_abundances", wmin, wmax, return_figure=False,
    #                  figure=(fig, ax), label="CNO Abundances\nHe = 2 x Solar\nC = 0.5 x Solar\nN = 7 x Solar")
    #
    # fig, ax = model_comparison(solar_smooth.copy(), "C2", disk, "_solar_clump", wmin, wmax, label=r"$f_{v}$ = 1", return_figure=True)
    # model_comparison(solar_clump.copy(), "C0", disk, "_solar_clump", wmin, wmax, label=r"$f_{v}$ = 0.1", return_figure=False,
    #                  figure=(fig, ax))
    #
    # fig, ax = model_comparison(solar_clump.copy(), "C0", disk, "_solar_cno_clump", wmin, wmax,
    #                            label=r"$f_{v}$ = 0.1" + "\nSolar Abundances", return_figure=True)
    # model_comparison(cno_clump.copy(), "C1", disk, "_solar_cno_clump", wmin, wmax,
    #                  label=r"$f_{v}$ = 0.1" + "\nCNO Abundances\nHe = 2 x Solar\nC = 0.5 x Solar\nN = 7 x Solar", return_figure=False,
    #                  figure=(fig, ax))

    solar_smooth = ["models/smooth/cv/solar/tde_cv.spec",
                    "models/smooth/agn/solar/tde_agn.spec"]
    wind_geos(solar_smooth.copy())

    plt.close()

    return


if __name__ == "__main__":
    main(len(argv), argv)
