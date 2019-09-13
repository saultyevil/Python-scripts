#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from platform import system
import numpy as np
from matplotlib import pyplot as plt
from consts import *
import py_plot_util as ppu
import tde_util as tu
from scipy.optimize import curve_fit
from astropy.modeling import models
from astropy import units as u
from astropy.modeling.blackbody import FLAM

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
        ax.axvline(x, ymax=0.93, linestyle="--", linewidth=0.45, color="k", zorder=1)
        ax.text(x, yloc + 0.2, lab, ha="center", va="center", rotation="vertical")

    return ax


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


def plot_uv_observations() -> None:
    """
    Plot four UV observations of TDE around ~55d. Also plot the SDSS composite
    QSO as a base for comparison.
    """

    iptf15af = tu.iptf15af_spec(SMOOTH, VERBOSE)
    iptf16fnl = tu.iptf16fnl_spec(SMOOTH, VERBOSE)

    nspec = 2
    spec_list = [iptf15af, iptf16fnl]
    spec_names = [r"iPTF15af $\Delta t = $52 d", r"iPTF16fnl $\Delta t = $51 d"]
    # name_x = [0.28, 1.25, 2.23, 3.27, 4.55]
    spec_z = [0.07897, 0.0163]
    bbT = [43300, 19000]
    dl = [350 * 1e6 * PARSEC, 66.6 * 1e6 * PARSEC]

    wmin = 1000
    wmax = 3000

    fig, ax = plt.subplots(2, 1, figsize=(12, 9))
    # ax.set_xlim(wmin, wmax)
    # ax = plot_line_id(ax, nspec + 0.17)
    for i in range(nspec):
        offset = 0  # * i
        wlength = spec_list[i][:, 0] / (1 + spec_z[i])
        flux = normalise_flux(ppu.smooth(spec_list[i][:, 1], SMOOTH, VERBOSE))
        ax[i].plot(wlength, flux + offset, label=spec_names[i])

        wl = wlength.copy()
        bb = models.BlackBody1D(temperature=bbT[i] * u.K, bolometric_flux=1e50)
        wav = np.linspace(wl[0], wl[-1], 500) * u.AA
        flux = normalise_flux(np.array(bb(wav).to(FLAM, u.spectral_density(wav))))
        ax[i].plot(wav, flux, label="bb @ {} K".format(bbT[i]))

        # exit(1)

        ax[i].set_xlim(wlength[0], wlength[-1])
        ax[i].set_ylim(0, 1)
        ax[i].set_xlabel(r"Rest Wavelength [$\AA$]", fontsize=13)
        ax[i].set_ylabel(r"Flux $F_{\lambda}$", fontsize=13)
        ax[i].legend(loc="upper right")

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("tde_uv_observations.png")
    plt.show()

    return


def fit_bb():
    """
    """

    def func(wl, t):
        lcm = wl * ANGSTROM
        x = H * C / (lcm * BOLTZMANN * t)
        y = 2 * H * C ** 2 / lcm ** 5
        b_lambda = y / (np.exp(x) - 1)
        return b_lambda

    iptf15af = tu.iptf15af_spec(SMOOTH, VERBOSE)
    wl = iptf15af[:, 0]
    fl = iptf15af[:, 1] * np.pi * 4 * (350 * 1e6 * PARSEC) ** 2
    idx = np.where(wl > 2000)[0]

    assert(len(wl[idx]) == len(fl[idx]))
    print(wl[idx], fl[idx])

    popt, pcov = curve_fit(func, wl[idx], fl[idx], p0=(43300))
    print(popt)

    fix, ax = plt.subplots(1, 1, figsize=(12, 5))
    ax.semilogy(wl, fl, label="data")
    ax.semilogy(wl, func(wl.copy(), popt), label="bb model")
    ax.legend()
    plt.show()

    return


if __name__ == "__main__":
    # fit_bb()
    plot_uv_observations()
