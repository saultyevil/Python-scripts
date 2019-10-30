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
import tde_util as tu

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

    nrows = 1
    ncols = 2
    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 8), squeeze=False)
    temperature = 43300

    civ = 1549 * ANGSTROM_TO_M
    civ_range = [[1460, 1590]]
    civ_Hz = C_SI / (civ * ANGSTROM_TO_M)
    siiv = 1400
    siiv_range = [[1290, 1400]]
    siiv_Hz = C_SI / (siiv * ANGSTROM_TO_M)
    nv = 1240
    nv_range = [[1140, 1340]]
    nv_Hz = C_SI / (nv * ANGSTROM_TO_M)

    ii = 0
    for i in range(nrows):
        for j in range(ncols):
            if ii > len(spec) - 1:
                break
            
            # Extract relevant wavlength range
            wl = spec[ii][:, 0]
            id1 = SpectrumUtils.get_wavelength_index(wl, civ_range[ii][0])
            id2 = SpectrumUtils.get_wavelength_index(wl, civ_range[ii][1])
            wl = wl[id1:id2] * ANGSTROM_TO_M

            wll = spec[ii][:, 0]
            id11 = SpectrumUtils.get_wavelength_index(wll, siiv_range[ii][0])
            id21 = SpectrumUtils.get_wavelength_index(wll, siiv_range[ii][1])
            wll = wll[id11:id21] * ANGSTROM_TO_M

            # Extract the flux and normalise the line profile
            # norm = blackbody_flux(bb_temp[ii], civ) 
            # norm *= np.pi * bb_radius[ii] ** 2 / dl[ii] ** 2
            # norm = spec[ii][id1, 1]
            fl = spec[ii][id1:id2, 1]
            fll = spec[ii][id11:id21, 1]

            # Convert to velocity space
            vel = C_SI * (wl / civ - 1)
            vel *= MS_TO_KMS   
            vell = C_SI * (wll / siiv - 1)
            vell *= MS_TO_KMS         

            # Finally plot :^)
            ax[i, j].plot(wl / ANGSTROM_TO_M, fl, label="civ")
            # ax[i, j].plot(wll / ANGSTROM_TO_M, fll, label="siiv")
            ax[i, j + 1].plot(vel, fl / 1e-15)
            # ax[i, j + 1].plot(vell, fll / 1e-15)
            ax[i, j].legend()
            ii += 1

    plt.show()

    return


if __name__ == "__main__":
    create_line_comparison()
