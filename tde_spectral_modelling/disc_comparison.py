#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from PyPython import SpectrumUtils
from matplotlib import pyplot as plt

C = 299792458
ANGSTROM = 1e-10


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
        # xnorm = x - 0.05 * x
        xnorm = x
        ax.text(xnorm, 2e8, lab, ha="center", va="center", rotation="vertical", fontsize=13)

    return ax


def plot_for_model(root, model_fname, disc_fname, los, smooth=10, xlims=([100, 1000], [3000, 6000]),
                   ylims=([1e-7, 1e3], [1e0, 1e2])):
    """
    Create a figure showing comparing the disc spectrum and model spectrum for
    short and long wavelengths.

    Parameters
    ----------
    root: str
        Root name of the simulation.
    model_fname: str
        The file name of the model spectrum.
    disc_fname: str
        The file name of the disc spectrum.
    los: str
        The inclination angle to consider
    smooth: int
        The amount of smoothing to use
    xlims: tuple(list, list)
        The xlimits to use for the two panels
    ylims: tuple(list, list)
        The ylimits to use for the two panels
    """

    nrows = 2
    fig, ax = plt.subplots(nrows, 1, figsize=(12, 8))

    model = SpectrumUtils.read_spec(model_fname)
    disc = SpectrumUtils.read_spec(disc_fname)

    model_wl = model["Lambda"].values.astype(float)
    model_fl = model[los].values.astype(float)
    disc_wl = disc["Lambda"].values.astype(float)
    disc_fl = disc[los].values.astype(float)

    model_x = model["Freq."].values.astype(float)
    model_y = model_fl * model_wl
    disc_x = disc["Freq."].values.astype(float)
    disc_y = disc_fl * disc_wl

    for i in range(nrows):
        wmin = C / xlims[i][0] / ANGSTROM
        wmax = C / xlims[i][1] / ANGSTROM
        ax[i].loglog(model_x, SpectrumUtils.smooth_spectrum(model_y, smooth), label="Model")
        ax[i].loglog(disc_x, SpectrumUtils.smooth_spectrum(disc_y, smooth), label="Disc")
        ax[i].set_xlabel(r"Frequency $\nu$", fontsize=15)
        ax[i].set_ylabel(r"$\nu$ F$_{\nu}$ [ergs s$^{-1}$ cm$^{-2}]$", fontsize=15)
        ax[i].set_ylim(ylims[i][0], ylims[i][1])
        ax[i].legend()
        plot_line_id(ax[i], SpectrumUtils.common_lines(True))
        ax[i].set_xlim(wmin, wmax)

    # fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("{}_disc_plot_i{}.png".format(root, los))
    plt.show()

    return


def main():
    """
    Main function for the script
    """

    root_directory = "/home/saultyevil/Dropbox/DiskWinds/PySims/tde"

    model_spec = root_directory + "/ztodo_iridis/models/clump/1e-1/cv/solar/tde_cv.spec"
    disc_spec = root_directory + "/ztodo_iridis/models/disc/cv/tde_cv.spec"
    for los in ["10", "30", "45", "60", "75", "85"]:
        if los == "85":
            plot_for_model("tde_cv", model_spec, disc_spec, los, ylims=([1e-7, 1e3], [1e-2, 1e2]))
        else:
            plot_for_model("tde_cv", model_spec, disc_spec, los)

    return


if __name__ == "__main__":
    main()
