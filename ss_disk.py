#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script contains multiple routines which are used to calculate various
properties of a Shakura and Sunyaev alpha accretion disk. Running the script
will plot a spectrum of a CV-esque accretion disk
"""

import numpy as np
from consts import *
from matplotlib import pyplot as plt


def planck_nu(t, nu):
    """
    The Planck function for a given frequency nu and temperature t.

    Parameters
    ----------
    t               float
                    Temperature, in Kelvin
    mu              float
                    Frequency in Hz

    Returns
    --------
    The monochromatic intensity for a black body at a frequency nu and temperature t,
    in units of ergs s^-1 cm^-2 Hz^-1
    """

    x = np.exp(H * nu / (BOLTZMANN * t))
    y = 2 * H * nu ** 3 / C ** 2

    return y * (x - 1) ** -1


def planck_lambda(t, lamda):
    """
    The Planck function for a given wavelength lamda and temperature t.

    Parameters
    ----------
    t               float
                    Temperature, in Kelvin
    lamda           float
                    Wavelength, in Angstroms

    Returns
    --------
    The monochromatic intensity for a black body at a wavelength lamda and temperature t,
    in units of ergs s^-1 cm^-2 A^-1
    """

    x = np.exp(H * C / (lamda * BOLTZMANN * t))
    y = 2 * H * C ** 2 / lamda ** 5

    return y * (x - 1) ** -1


def t_eff(r: float, m: float, mdot: float, rin: float) -> float:
    """
    Return the effective temperature at a position R in an accretion disk. See Frank, King, Raine 1995 Chapter 5.

    Parameters
    ----------
    r               float
                    The position on the disk given in cm
    m               float
                    The mass of the central object given in g
    rin             float
                    The radius of the central object given in cm
    mdot            float
                    The mass accretion rate of the accretion disk given in g s^-1

    Returns
    --------
    teff            float
                    The effective temperature of the accretion disk at point R on the accretion disk given in K
    """

    if r < rin:
        print("ss_disk.t_eff: r < rin")
        exit(1)

    teff = (3 * G * m * mdot) / (8 * np.pi * r ** 3 * STEFAN_BOLTZMANN) * (1 - (rin / r) ** 0.5)
    teff = teff ** 0.25

    return teff


def disk_spectrum(min: float, max: float, m: float, mdot: float, rin: float, rout: float, nu_or_lambda: str = "lambda",
                  npoints: int = 500, nrings: int = 500, verbose: bool = False) -> np.array:
    """
    Calculate the spectrum of an accretion disk of inner radius rin and outer radius rout with accretion rate mdot.
    Note that the spectrum of an accretion disk is an ensemble of black body distributions.

    Parameters
    ----------
    min, max        float
                    Wavelength or frequency limits for the spectrum
    m               float
                    Mass of the central object in units of msol
    mdot            float
                    Mass accretion rate of the disk in units of msol/yr
    rin, rout       float
                    Inner and outer radius of the accretion disk (assume rin = r_co) given in cm
    nu_or_lambda    str, optional
                    Sets whether the flux will be per unit wavelength or per unit frequency
    npoints         int, optional
                    The number of points on the wavelength of frequency grid
    nrings          int, optional
                    The number of disk annuli which makes up the spectrum
    verbose         bool, optional
                    Enable verbose output

    Returns
    -------
    disk_spec       np.array[float]
                    The disk spectrum
                        - disk_spec[:, 0] the wavelength or frequency range
                        - disk_spec[:, 1] the monochromatic flux per unit wavelength or frequency
    """

    # Convert into CGS units
    m *= MSOL
    mdot *= MSOL_PER_YEAR

    max_temperature = 0
    r_grid = np.linspace(rin, rout, nrings)
    unit_grid = np.linspace(min, max, npoints)
    disk_spec = np.zeros((npoints, 2))
    disk_spec[:, 0] = unit_grid

    for i in range(nrings - 1):
        # Use midpoint of annulus as point on r grid
        r = (r_grid[i + 1] + r_grid[i]) / 2.0
        annulus_area = np.pi * (r_grid[i + 1] ** 2 - r_grid[i] ** 2)
        temperature = t_eff(r, m, mdot, rin)
        if temperature > max_temperature:
            max_temperature = temperature
        # Determine which function to get the bb flux from
        bb = 0
        if nu_or_lambda == "lambda":
            bb = planck_lambda(temperature, unit_grid)
        elif nu_or_lambda == "nu":
            bb = planck_nu(temperature, unit_grid)
        else:
            print("value of {} for nu_or_lambda invalid".format(nu_or_lambda))
            exit(1)
        # Finally calculate the monochromatic flux and increment the spectrum
        disk_spec[:, 1] += bb * annulus_area

    if verbose:
        print("ss_disk.get_disk_spectrum: maximum effective disk temperature {} K".format(max_temperature))

    return disk_spec


def main():
    """
    Main function

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    object = "cv"
    if object == "cv":
        min = 1e11
        max = 1e18
        m = 1.0
        mdot = 1e-8
        rin = 7e8
        rout = rin * 1e5
        nu_or_lambda = "nu"
    elif object == "tde":
        min = 1000 * 1e-10
        max = 3000 * 1e-10
        rin = 2.65e13
        # rout = 5.00e14
        rout = 1e4 * rin
        m = 3e7
        mdot = 4e-2
        nu_or_lambda = "lambda"
    else:
        print("Unknown object")
        exit(1)

    # Get spectrum from above functions
    i = 0
    dist = 1
    disk_spec = disk_spectrum(min, max, m, mdot, rin, rout, nu_or_lambda, npoints=5000, nrings=5000)
    flux = disk_spec[:, 1] * (4 * PI * np.cos(i)) / dist ** 2
    lamda = disk_spec[:, 0]

    # Plot spectrum
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.plot(np.log10(lamda), np.log10(lamda*flux))
    ax.set_ylim(10, 40)
    ax.set_xlim(11, 18)
    ax.set_xlabel(r"log[$\nu$]", fontsize=15)
    ax.set_ylabel(r"log[$\nu$L$_{\nu}$]", fontsize=15)
    plt.show()

    return


if __name__ == "__main__":
    main()
