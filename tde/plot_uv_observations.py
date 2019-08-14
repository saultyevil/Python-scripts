#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import sys
from platform import system
import numpy as np
from matplotlib import pyplot as plt

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

from consts import *
import py_plot_util as ppu
import tde_util as tu

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
    ["N III]", 1759],
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
    fmin = np.min(input_flux)
    output_flux = (input_flux - fmin) / (fmax - fmin)

    return output_flux


def blackbody_flux(T: float, lamda: np.ndarray) -> np.ndarray:
    """
    Return the blackbody flux as a function of wavelength in Angstroms.

    Parameters
    ----------
    T: float
        The temperature of the blackbody
    lamda: np.ndarray[float]
        The wavelength range to calculate the blackbody flux over.

    Returns
    -------
    The monochromatic intensity for a black body at a wavelength lamda and
    temperature t, in units of ergs s^-1 cm^-2 A^-1.
    """

    # convert lambda into cm
    lcm = lamda * ANGSTROM

    x = H * C / (lcm * BOLTZMANN * T)
    y = 2 * H * C ** 2 / lcm ** 5
    b_lambda = y / (np.exp(x) - 1)

    return b_lambda * ANGSTROM


def plot_uv_observations() -> None:
    """
    Plot four UV observations of TDE around ~55d. Also plot the SDSS composite
    QSO as a base for comparison.
    """

    iptf15af = tu.iptf15af_spec(SMOOTH, VERBOSE)
    asassn14li = tu.asassn14li_spec(SMOOTH, VERBOSE)
    iptf16fnl = tu.iptf16fnl_spec(SMOOTH, VERBOSE)
    at2018zr = tu.at2018zr_spec(SMOOTH, VERBOSE)
    composite_qso = tu.sdss_qso_spec(VERBOSE)

    nspec = 5
    spec_list = [composite_qso, asassn14li, iptf15af, iptf16fnl, at2018zr]
    spec_names = ["SDSS Composite QSO", r"ASASSN14li $\Delta t = $60 d", r"iPTF15af $\Delta t = $52 d",
                  r"iPTF16fnl $\Delta t = $51 d", r"AT2018zr $\Delta t = $59 d"]
    name_x = [0.28, 1.25, 2.23, 3.27, 4.55]
    spec_z = [0,  0.02058, 0.07897, 0.0163, 0.071]
    bbT = [0, 35000, 43300, 19000, 22000]

    wmin = 1000
    wmax = 3000

    fig, ax = plt.subplots(1, 1, figsize=(9.5, 11))
    ax.set_xlim(wmin, wmax)
    ax = plot_line_id(ax, nspec + 0.17)
    for i in range(nspec):
        offset = 1 * i
        wlength = spec_list[i][:, 0] / (1 + spec_z[i])
        flux = normalise_flux(ppu.smooth_1d_array(spec_list[i][:, 1], SMOOTH, VERBOSE))
        ax.plot(wlength, flux + offset, label=spec_names[i])
        if bbT[i] != 0:
            twl = wlength
            bbfl = blackbody_flux(bbT[i], twl)
            bbfl = normalise_flux(bbfl)
            ax.plot(twl, bbfl + offset, linestyle="--", alpha=0.5, color="k")
        ax.text(2200, name_x[i], spec_names[i], fontsize=13)
    ax.set_ylim(0, nspec + 0.55)
    ax.set_xlabel(r"Rest Wavelength [$\AA$]")
    ax.set_ylabel(r"Normalised Flux $F_{\lambda}$ + Offset")

    fig.tight_layout()
    plt.savefig("tde_uv_observations.png")
    plt.show()

    return


if __name__ == "__main__":
    plot_uv_observations()
