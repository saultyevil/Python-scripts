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

from consts import H, C, BOLTZMANN, PARSEC
import tde_spectra as tu

C_SI = 299792458
C_CM = 1e2 * C_SI
ANGSTROM_TO_M = 1e-10
MS_TO_KMS = 1e-3
SMOOTH = 5
VERBOSITY = False

bb_temp = [43300,
           35000,
           19000,
           22000]
bb_radius = [1.35e14,
             1.35e14,
             1.1e14,
             4e14]
dl = np.array([358,
               90,
               67,
               337]) * 1e6 * PARSEC


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
    lcm = lamda * 1e-8

    x = H * C / (lcm * BOLTZMANN * T)
    y = 2 * H * C ** 2 / lcm ** 5
    b_lambda = y / (np.exp(x) - 1)

    return b_lambda * 1e-8


def create_line_comparison():
    """Compare line profiles in velocity space"""

    iptf15af = tu.iptf15af_spec(SMOOTH, VERBOSITY)
    asassn14li = tu.asassn14li_spec(SMOOTH, VERBOSITY)
    at2018zr = tu.at2018zr_spec(SMOOTH, VERBOSITY)
    iptf16fnl = tu.iptf16fnl_spec(SMOOTH, VERBOSITY)
    iptf16fnl = iptf16fnl[5:, :]

    spec = [iptf15af, asassn14li, at2018zr, iptf16fnl]
    spec = [iptf15af]

    ncols = 2
    nrows = len(spec)
    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 8), squeeze=False)

    carbon4_wavelength = 1549 * ANGSTROM_TO_M
    carbon_wl_range = [[1460, 1590]]
    silicon4_wavelength = 1400 * ANGSTROM_TO_M
    siiv_range = [[1290, 1400]]
    nitrogen5_wavelength = 1240 * ANGSTROM_TO_M
    nv_range = [[1140, 1340]]

    ii = 0
    for i in range(nrows):

        # Extract C IV
        wl_carbon = spec[i][:, 0]
        id1 = SpectrumUtils.get_wavelength_index(wl_carbon, carbon_wl_range[i][0])
        id2 = SpectrumUtils.get_wavelength_index(wl_carbon, carbon_wl_range[i][1])
        wl_carbon = wl_carbon[id1:id2]
        fl_carbon = spec[i][id1:id2, 1]
        # Extract Si IV
        wl_silicon = spec[i][:, 0]
        id1 = SpectrumUtils.get_wavelength_index(wl_silicon, siiv_range[i][0])
        id2 = SpectrumUtils.get_wavelength_index(wl_silicon, siiv_range[i][1])
        wl_silicon = wl_silicon[id1:id2]
        fl_silicon = spec[i][id1:id2, 1]
        # Extract N V
        wl_nitrogen = spec[i][:, 0]
        id1 = SpectrumUtils.get_wavelength_index(wl_nitrogen, nv_range[i][0])
        id2 = SpectrumUtils.get_wavelength_index(wl_nitrogen, nv_range[i][1])
        wl_nitrogen = wl_nitrogen[id1:id2]
        fl_nitrogen = spec[i][id1:id2, 1]

        # Convert the line profiles into velocity space
        carbon4_vel = C_SI * (wl_carbon * ANGSTROM_TO_M / carbon4_wavelength - 1) * MS_TO_KMS
        silicon4_vel = C_SI * (wl_silicon * ANGSTROM_TO_M / silicon4_wavelength - 1) * MS_TO_KMS
        nitrogen5_vel = C_SI * (wl_nitrogen * ANGSTROM_TO_M / nitrogen5_wavelength - 1) * MS_TO_KMS

        ax[i, 0].plot(carbon4_vel, fl_carbon, label=r"C IV")
        ax[i, 0].plot(silicon4_vel, fl_silicon, label=r"Si IV")
        ax[i, 0].legend()
        ax[i, 0].set_xlabel("Velocity km/s")
        ax[i, 0].set_ylabel("Normalised Flux")

        ax[i, 1].plot(carbon4_vel, fl_carbon, label=r"C IV")
        ax[i, 1].plot(nitrogen5_vel, fl_nitrogen, label=r"N V")
        ax[i, 1].legend()
        ax[i, 1].set_xlabel("Velocity km/s")
        ax[i, 1].set_ylabel("Normalised flux")

    plt.show()

    return


if __name__ == "__main__":
    create_line_comparison()
