#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from consts import *
import numpy as np
from matplotlib import pyplot as plt


def planck_lambda(t, lamda):
    """
    The Planck function for a given wavelength lamda and temperature t.

    Parameters
    ----------
    t       float
            Temperature, in Kelvin
    lamda   float
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
    r       float
            The position on the disk given in cm
    m       float
            The mass of the central object given in g
    rin     float
            The radius of the central object given in cm
    mdot    float
            The mass accretion rate of the accretion disk given in g s^-1

    Returns
    --------
    teff    float
            The effective temperature of the accretion disk at point R on the accretion disk given in K
    """

    if r < rin:
        print("ss_disk.t_eff: r < rin")
        exit(1)

    teff = (3 * G * m * mdot) / (8 * np.pi * r ** 3 * STEFAN_BOLTZMANN)
    teff *= (1 - (rin / r) ** 0.5)
    teff = teff ** 0.25

    return teff


def disk_spectrum(wmin: float, wmax: float, m: float, mdot: float, rin: float, rout: float, nlamda: int = 500,
                  nrings: int = 500, verbose: bool = False) -> np.array:
    """
    Calculate the spectrum of an accretion disk of inner radius rin and outer radius rout with accretion rate mdot.
    Note that the spectrum of an accretion disk is an ensemble of black body distributions.

    Parameters
    ----------
    wmin, wmax      float
                    Wavelength limits for spectrum in Angstroms
    m               float
                    Mass of the central object in units of msol
    mdot            float
                    Mass accretion rate of the disk in units of msol/yr
    rin, rout       float
                    Inner and outer radius of the accretion disk (assume rin = r_co) given in cm
    nlamda          int, optional
                    The number of points on the wavelength grid
    nrings          int, optional
                    The number of disk annuli
    verbose         bool, optional
                    Enable verbose output

    Returns
    -------
    disk_spec       2 x nlamda array
                    The disk spectrum
                        - disk_spec[:, 0] the wavelength range in Angstrom
                        - disk_spec[:, 1] the monochromatic flux ergs s^-1 cm^-2 A^-1
    """

    # Convert into CGS units
    m *= MSOL
    mdot *= MSOL_PER_YEAR

    # Create the radius and wavelength grid
    rgrid = np.linspace(rin, rout, nrings)
    lamda = np.linspace(wmin, wmax, nlamda)
    disk_spec = np.zeros((nlamda, 2))
    disk_spec[:, 0] = lamda

    # Loop over each annulus
    teff_max = 0
    for i in range(nrings - 1):
        # Calculate midpoint for annulus, area of annulus and the effective temperature at the midpoint
        r = (rgrid[i + 1] + rgrid[i]) / 2.0
        teff = t_eff(r, m, mdot, rin)
        if teff > teff_max:
            teff_max = teff
        area = np.pi * (rgrid[i + 1] ** 2 - rgrid[i] ** 2)
        # Calculate monochromatic flux for each wavelength
        disk_spec[:, 1] += planck_lambda(teff, lamda) * area * 2 * np.pi

    if verbose:
        print("ss_disk.get_disk_spectrum: maximum effective disk temperature {} K".format(teff_max))

    return disk_spec


def main():
    """
    Nothing happens :-)

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    object = "cv"

    if object == "cv":
        wmin = 0.1
        wmax = 5e10
        rin = 7e8
        rout = 2.4e10
        m = 0.8
        mdot = 1e-8
    elif object == "tde":
        wmin = 5e2
        wmax = 1e4
        rin = 2.65e13
        rout = 5.00e14
        m = 3e7
        mdot = 4e-2
    else:
        print("Unknown object :-)")
        exit(1)


    # lamda = np.linspace(1e-5, 3e-4, int(5e3))
    # t = 5e3
    # f = planck_lambda(t, lamda)
    # plt.plot(lamda, f)
    # plt.show()
    # exit(1)

    disk_spec = disk_spectrum(wmin, wmax, m, mdot, rin, rout, verbose=True)

    flux = disk_spec[:, 1]
    lamda = disk_spec[:, 0]

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.plot(np.log10(lamda), np.log10(lamda*flux))
    ax.set_xlabel("Wavelength, $\AA$")
    ax.set_ylabel("Flux")
    plt.show()

    return


if __name__ == "__main__":
    main()
