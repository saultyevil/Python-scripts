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
import tde_spectra as tu

SMOOTH = 5
VERBOSITY = False


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

    spec_list = [asassn14li, iptf15af, iptf16fnl, at2018zr]
    spec_names = [r"ASASSN14li $\Delta t = $60 d", r"iPTF15af $\Delta t = $52 d", r"iPTF16fnl $\Delta t = $51 d",
                  r"AT2018zr $\Delta t = $59 d"]
    spec_z = [0.02058, 0.07897, 0.0163, 0.071]
    bb_temp = [35000, 43300, 19000, 22000]
    bb_radius = [1.35e14, 1.35e14, 1.1e14, 4e14]
    dl = np.array([90, 358, 67, 337]) * 1e6 * PARSEC

    wmin = 1000
    wmax = 3000

    fig, ax = plt.subplots(2, 2, figsize=(15, 11))

    ii = 0
    for i in range(2):
        for j in range(2):
            if ii > 3:
                break
            wlength = spec_list[ii][:, 0] / (1 + spec_z[ii])
            flux = SpectrumUtils.smooth_spectrum(spec_list[ii][:, 1], SMOOTH)
            ax[i, j].semilogy(wlength, flux, label=spec_names[ii])
            if bb_temp[ii]:
                twl = wlength
                bbfl = blackbody_flux(bb_temp[ii], twl)
                bbfl *= np.pi * bb_radius[ii] ** 2 / dl[ii] ** 2
                ax[i, j].semilogy(twl, bbfl, linestyle="--", alpha=0.5, color="k",
                                  label="T = {:1.3e} K R = {:1.3e} cm".format(bb_temp[ii], bb_radius[ii]))
            ax[i, j].legend()
            ax[i, j].set_xlim(wmin, wmax)
            ii += 1

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("tde_uv_observations.png")
    plt.show()

    return


if __name__ == "__main__":
    plot_uv_observations()
