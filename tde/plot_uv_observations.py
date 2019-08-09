#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import sys
sys.path.append("..")

import numpy as np
import py_plot_util as ppu
import tde_util as tu
from matplotlib import pyplot as plt


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


def plot_line_id(ax: plt.Axes, yloc: float) -> plt.Axes:
    """
    Plot labels and vertical lines to indicate important atomic transitions.

    Parameters
    ----------
    ax: plt.Axes
        The Axes object to add line ID labels to.

    yloc: float
        The y coordinate where to place the labels. Note that this is value is
        then padded by 0.2 :^).

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
        ax.axvline(x, ymax=0.93, linestyle="--", linewidth=0.45, color="k")
        ax.text(x, yloc + 0.2, lab, ha="center", va="center", rotation="vertical")

    return ax


def normalise_flux(input_flux: np.ndarray):
    """
    Normalise the flux between 0 and 1.

    Parameters
    ----------
    input_flux: np.ndarray
        The input flux to rescale.

    Returns
    -------
    output_flux: np.ndarray
        The rescaled input flux.
    """

    assert(type(input_flux) == np.ndarray), "input flux must be a Numpy array"

    fmax = np.max(input_flux)
    output_flux = input_flux / fmax

    return output_flux


def plot_uv_observations():
    """
    Plot four UV observations of TDE around ~55d. Also plot the SDSS composite
    QSO as a base for comparison.
    """

    # Load the spectra into memory
    iptf15af = tu.iPTF15af_spec(SMOOTH, VERBOSE)
    asassn14li = tu.ASSASN14li_spec(SMOOTH, VERBOSE)
    iptf16fnl = tu.iPTF16fnl_spec(SMOOTH, VERBOSE)
    at2018zr = tu.AT2018zr_spec(SMOOTH, VERBOSE)
    composite_qso = np.loadtxt("/home/saultyevil/PySims/tde/observed_spec/sdss_composite_qso.dat")

    # Book keeping lists for plotting the spectra in a loop
    nspec = 5
    spec_list = [composite_qso, asassn14li, iptf15af, iptf16fnl, at2018zr]
    spec_names = ["SDSS Composite QSO", "ASASSN14li: 60d", "iPTF15af: 52d",  "iPTF16fnl: 51d", "AT2018zr: 59d"]
    spec_z = [0,  0.02058, 0.07897, 0.0163, 0.071]

    fig, ax = plt.subplots(1, 1, figsize=(9.5, 11))
    ax.set_xlim(1000, 3000)
    ax = plot_line_id(ax, nspec + 0.17)
    for i in range(nspec):
        offset = 1 * i
        wlength = spec_list[i][:, 0] / (1 + spec_z[i])
        flux = normalise_flux(ppu.smooth_1d_array(spec_list[i][:, 1], SMOOTH, VERBOSE))
        ax.plot(wlength, flux + offset, label=spec_names[i])
    ax.legend(loc="lower right")
    ax.set_ylim(0, nspec + 0.55)
    ax.set_xlabel(r"Rest Wavelength [$\AA$]")
    ax.set_ylabel(r"Normalised Flux $F_{\lambda}$ + Offset")

    fig.tight_layout()
    plt.savefig("tde_uv_observations.png")
    plt.show()

    return


if __name__ == "__main__":
    plot_uv_observations()
