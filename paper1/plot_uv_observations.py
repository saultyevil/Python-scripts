#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import sys
from platform import system
import numpy as np
from matplotlib import pyplot as plt
from PyPython import SpectrumUtils

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

from consts import *
import tde_util as tu

SMOOTH = 15
VERBOSITY = False

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


NORM = 6.652348e-15


def normalise_flux(input_flux: np.ndarray):
    """
    Rescale the input flux between 0 and 1.

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

    output_flux = input_flux / NORM

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

    iptf15af = tu.iptf15af_spec(SMOOTH, VERBOSITY)
    asassn14li = tu.asassn14li_spec(SMOOTH, VERBOSITY)
    iptf16fnl = tu.iptf16fnl_spec(SMOOTH, VERBOSITY)
    at2018zr = tu.at2018zr_spec(SMOOTH, VERBOSITY)
    # composite_qso = tu.sdss_qso_spec(VERBOSE)
    lobal = tu.lobal_qso_spec(VERBOSITY)
    composite_qso = lobal

    nspec = 5
    spec_list = [lobal, asassn14li, iptf16fnl, iptf15af, at2018zr]
    spec_names = ["Composite LoBALQSO",
                  r"ASASSN14li $\Delta t = $60 d" + "\n" + r"$T_{bb} = 35,000$K",
                  r"iPTF16fnl $\Delta t = $51 d" + "\n" + r"$T_{bb} = 19,000$K",
                  r"iPTF15af $\Delta t = $52 d" + "\n" + r"$T_{bb} = 43,300$K",
                  r"AT2018zr $\Delta t = $59 d" + "\n" + r"$T_{bb} = 22,000$K"]
    name_x = [0.8, 6.1, 100, 1200, 60000]
    spec_z = [0, 0.02058, 0.0163, 0.07897, 0.071]
    bb_temp = [0, 35000, 19000, 43300, 22000]
    bb_radius = [0, 1.35e14, 1.1e14, 1.35e14, 4e14]
    dl = np.array([0, 90, 67, 358, 337]) * 1e6 * PARSEC

    wmin = 1000
    wmax = 3000

    offsets = [1, 10, 700, 2e4, 5e5]

    fig, ax = plt.subplots(1, 1, figsize=(9.5, 11))
    ax.set_xlim(wmin, wmax)
    for i in range(nspec):
        offset = offsets[i]
        wlength = spec_list[i][:, 0] / (1 + spec_z[i])
        flux = SpectrumUtils.smooth_spectrum(spec_list[i][:, 1], SMOOTH)
        flux = normalise_flux(flux)
        ax.semilogy(wlength, flux * offset, label=spec_names[i])

        if bb_temp[i]:
            twl = wlength
            bbfl = blackbody_flux(bb_temp[i], twl)
            bbfl *= np.pi * bb_radius[i] ** 2 / dl[i] ** 2
            bbfl = normalise_flux(bbfl)
            ax.semilogy(twl, bbfl * offset, linestyle="--", alpha=0.5, color="k")
        ax.text(2200, name_x[i], spec_names[i], fontsize=13)

    ax.set_ylim(0.05, 4e6)
    ax = plot_line_id(ax)

    ax.set_xlabel(r"Rest Wavelength [$\AA$]", fontsize=15)
    ax.set_ylabel(r"Normalised Flux $F_{\lambda}$ + Offset", fontsize=15)
    ax.tick_params(axis="x", labelsize=13)
    ax.tick_params(axis="y", labelsize=13)

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("tde_uv_observations.png")
    plt.show()

    return


if __name__ == "__main__":
    plot_uv_observations()
