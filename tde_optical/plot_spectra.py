#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from sys import argv
from platform import system
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
from PyPython import SpectrumUtils
from path import *

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

DISPLAY = False
if len(argv) > 1:
    DISPLAY = True
SMOOTH = 5
VERBOSITY = False

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
    ["C IV", 1549],
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
    [r"H$\beta$", 4861],
    [r"H$\alpha$", 6564],
]

colors = ["C0", "C1", "C2"]

alpha = 0.75

lstyle = [
    "-", "--", "-."
]

incl = [
    "10", "60", "75"
]

comment = [
    ":\n outside wind",
    ":\n in-wind",
    ":\n outside wind"
]


def plot_line_id(ax: plt.Axes, offset: Union[float, int] = 55) -> plt.Axes:
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
        x = x - offset
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        ax.text(xnorm, 0.90, lab, ha="center", va="center", rotation="vertical", transform=ax.transAxes, fontsize=13)

    return ax


def add_axes_labels(fig, ax):
    """Add the axes labels"""

    fig.text(0.5, 0.03, r"Wavelength $\lambda$ [$\AA$]", ha="center", va="center", rotation="horizontal", fontsize=15)
    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)
    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])

    return fig, ax


def panel_all_grids(directories: List[str], line_colours: List[str], wmin: float, wmax: float, ylims: List[float],
                    labels: List[str], extra_file_name: str):
    """Creates a 3 x 3 figure for showing the entire grid in one figure."""

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde_optical/grid/round1/"
    else:
        pdir = "/home/saultyevil/PySims/tde_optical/grid/round1/"
    t = get_fiducial_uv_model()

    incl = [
        "10", "60", "75"
    ]

    ncols = 3
    nrows = 3
    fig, ax = plt.subplots(nrows, ncols, figsize=(13.3, 14), sharex="col", sharey="row")
    ii = 0

    if extra_file_name == "the_uv_grid":
        offset = 25
    else:
        offset = 55

    for i in range(nrows):
        for j in range(ncols):

            if j == 0:
                fid_label = r"M$_{BH}$ = 5$\times$10$^6$ M$_{\odot}$"
            elif j == 1:
                fid_label = r"R$_{min}$ = 1 R$_{ISCO}$"
            elif j == 2:
                fid_label = r"V$_{\infty}$ = 1 V$_{esc}$"
            else:
                fid_label = None

            for k in range(len(directories[j])):
                mt = SpectrumUtils.read_spec(pdir + directories[j][k] + "/tde_uv.spec")
                wl = mt["Lambda"].values.astype(float)
                uv_fid_wl = t["Lambda"].values.astype(float)
                try:
                    fl = mt[incl[i]].values.astype(float)
                    uv_fid_fl = t[incl[i]].values.astype(float)
                except KeyError:
                    continue

                ax[i, j].semilogy(wl, SpectrumUtils.smooth(fl, SMOOTH), label=labels[j][k], linewidth=3,
                                  color=line_colours[k], linestyle=lstyle[k], alpha=alpha)
                if k == 0:
                    ax[i, j].semilogy(uv_fid_wl, SpectrumUtils.smooth(uv_fid_fl, SMOOTH), "k-", linewidth=3,
                                      alpha=alpha, zorder=0, label=fid_label)
                ax[i, j].set_xlim(wmin, wmax)
                if i == 0:  # Lazy programming :-)
                    if extra_file_name == "the_uv_grid":
                        ax[i, j].set_ylim(ylims)
                    else:
                        ax[i, j].set_ylim(9e-5, 4e-2)
                else:
                    ax[i, j].set_ylim(ylims)

                ax[i, j] = plot_line_id(ax[i, j], offset=offset)
            if i == nrows - 1:
                ax[i, j].legend(loc="lower center", fontsize=11, ncol=2)
            ii += 1

    fig, ax = add_axes_labels(fig, ax)
    fig.subplots_adjust(hspace=0, wspace=0)

    comment = [": outside wind", ": inside wind", ": outside wind"]
    for i in range(nrows - 1):
        ax[i, -1].text(0.5, 0.1, r"$i = $" + incl[i] + r"$^{\circ}$" + comment[i], ha="center", va="center",
                       rotation="horizontal", fontsize=15, transform=ax[i, -1].transAxes)
    ax[2, -1].text(0.5, 0.2, r"$i = $" + incl[2] + r"$^{\circ}$" + comment[2], ha="center", va="center",
                   rotation="horizontal", fontsize=15, transform=ax[2, -1].transAxes)

    fname = "spectra/"
    if extra_file_name:
        fname += extra_file_name

    fig.savefig(fname + ".pdf", dpi=300)
    fig.savefig(fname + ".png", dpi=300)

    if DISPLAY:
        plt.show()
    else:
        plt.close()

    return fig, ax


def panel_1_by_3(directories: List[str], line_colours: List[str], optical_wmin: float, optical_wmax: float,
                 optical_ylims: List[float], labels: List[str], fid_label: str, extra_file_name: str) \
        -> Union[None, Tuple[plt.Figure, plt.Axes]]:
    """
    Plot a 1 x 3 grid of model comparisons.
    """

    modelspecs = get_the_models(directories, "tde_uv.spec")
    t = get_fiducial_uv_model()

    ncols = 3
    nrows = 1
    fig, ax = plt.subplots(nrows, ncols, figsize=(16, 5.5), sharey="row")

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
            uv_fid_wl = t["Lambda"].values.astype(float)
            try:
                fl = modelspecs[j][incl[i]].values.astype(float)
                uv_fid_fl = t[incl[i]].values.astype(float)
            except KeyError:
                continue
            if j == 0:
                ax[i].semilogy(uv_fid_wl, SpectrumUtils.smooth(uv_fid_fl, SMOOTH), "k-", linewidth=2, alpha=alpha,
                               zorder=0, label=fid_label)
            ax[i].semilogy(wl, SpectrumUtils.smooth(fl, SMOOTH), lstyle[j], label=labels[j], linewidth=3,
                           color=line_colours[j], alpha=alpha)
        ax[i].set_xlim(optical_wmin, optical_wmax)
        ax[i].set_ylim(optical_ylims)
        ax[i] = plot_line_id(ax[i])

    ax[0].legend(loc="lower center", fontsize=11, ncol=2)

    add_axes_labels(fig, ax)
    fig.subplots_adjust(hspace=0, wspace=0)

    comment = [": outside wind", ": inside wind", ": outside wind"]
    for i in range(ncols):
        ax[i].set_title(r"$i = $" + incl[i] + r"$^{\circ}$" + comment[i], fontsize=15)

    fname = "spectra/"
    if extra_file_name:
        fname += extra_file_name

    fig.savefig(fname + ".pdf", dpi=300)
    fig.savefig(fname + ".png", dpi=300)

    if DISPLAY:
        plt.show()
    else:
        plt.close()

    return fig, ax


def single_panel(directories: List[str], line_colours: List[str], labels: List[str], extra_file_name: str = "") \
        -> Union[None, Tuple[plt.Figure, plt.Axes]]:
    """Plot a single panel"""

    modelspecs = get_the_models(directories, "tde_uv.spec")
    t = get_fiducial_uv_model()

    nrows = ncols = 1
    fig, ax = None, None

    sm = 50

    for j in range(len(incl)):
        fig, ax = plt.subplots(nrows, ncols, figsize=(15, 7))
        for i in range(len(modelspecs)):
            if type(modelspecs[i]) == type(SKIP):
                continue
            wl = modelspecs[i]["Lambda"].values.astype(float)
            try:
                fl = modelspecs[i][incl[j]].values.astype(float)
            except KeyError:
                continue
            ax.loglog(wl, SpectrumUtils.smooth(fl, sm), "-", label=labels[i], linewidth=3,
                      color=line_colours[i], alpha=alpha)
        ax.text(0.83, 0.9, r"$i = $" + incl[j] + r"$^{\circ}$" + comment[j], ha="center", va="center",
                rotation="horizontal", fontsize=15, transform=ax.transAxes)

        # Plot the fiducial model
        # fid_wl = t["Lambda"].values.astype(float)
        # fid_fl = t[incl[j]].values.astype(float)
        # ax.loglog(fid_wl, SpectrumUtils.smooth(fid_fl, sm), "k-", linewidth=2, alpha=0.3, zorder=0)

        ax.set_xlabel(r"Wavelength $\lambda$ [$\AA$]")
        ax.set_ylabel(r"Flux $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", )
        ax.legend(loc="upper left", fontsize=13)
        ax.set_ylim(1e-5, 1e0)
        # ax = SpectrumUtils.plot_line_ids(ax, LINES, logx=True, offset=55)

        fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])

        fname = "spectra/"
        if extra_file_name:
            fname += extra_file_name
        fname += "_i{}".format(incl[j])

        fig.savefig(fname + ".pdf", dpi=300)
        fig.savefig(fname + ".png", dpi=300)

    if DISPLAY:
        plt.show()
    else:
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
    optical_wmax = 7350
    uv_wmin = 900
    uv_wmax = 2400

    # Plot the "best" lines

    panel_1_by_3(best_lines_grid, colors, optical_wmin, optical_wmax, [1e-4, 2e-2], best_lines_labels,
                 "Original", "optical_fiducial")

    # Plot the grids on a single panel, full wavelength range, for each inclination angle

    single_panel(mbh_grid.copy(), colors, mbh_labels, extra_file_name="sed_Mbh")
    single_panel(rmin_grid.copy(), colors, rmin_labels, extra_file_name="sed_Rmin")
    single_panel(vinf_grid.copy(), colors, vinf_labels, extra_file_name="sed_Vinf")

    # Individual grids for just optical

    panel_1_by_3(mbh_grid.copy(), colors, optical_wmin, optical_wmax, [1e-5, 4e-2], mbh_labels,
                 r"M$_{BH}$ = 5$\times$10$^6$ M$_{\odot}$", "optical_Mbh")
    panel_1_by_3(rmin_grid.copy(), colors, optical_wmin, optical_wmax, [2e-4, 2e-2], rmin_labels,
                 r"R$_{min}$ = 1 R$_{ISCO}$", "optical_Rmin")
    panel_1_by_3(vinf_grid.copy(), colors, optical_wmin, optical_wmax, [2e-4, 1e-2], vinf_labels,
                 r"V$_{\infty}$ = 1 V$_{esc}$", "optical_Vinf")

    # Optical + UV together

    the_grid = [mbh_grid.copy(), rmin_grid.copy(), vinf_grid.copy()]
    the_labels = [mbh_labels, rmin_labels, vinf_labels]

    panel_all_grids(the_grid.copy(), colors, optical_wmin, optical_wmax, [9e-6, 2e-2], the_labels, "the_optical_grid")
    panel_all_grids(the_grid.copy(), colors, uv_wmin, uv_wmax, [1e-6, 3], the_labels, "the_uv_grid")

    return


if __name__ == "__main__":
    main(len(argv), argv)
