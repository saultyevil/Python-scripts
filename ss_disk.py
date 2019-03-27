#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from consts import *
import numpy as np
from matplotlib import pyplot as plt


def planck_lambda(lamda, t):
    """
    The Planck function for a given wavelength lamda and temperature t.

    Parameters
        lamda           wavelength in Angstroms
        t               the temperature in Kelvin

    Returns
        The monochromatic intensity for a blackbody at a wavelength lamda and temperature t, in units of ergs s^-1 A^-1
    """

    x = np.exp(H * C / (lamda * BOLTZMANN * t))
    y = 2 * H * C ** 2 / lamda ** 5

    return y * 1 / (x - 1)


def t_eff(r: float, m: float, mdot: float, rin: float) -> float:
    """
    Return the effective temperature at a position R in an accretion disk. See Frank, King, Raine 1995 Chapter 5.

    Parameters
        r               the position on the disk in cm
        m               the mass of the central object in g
        rin             the radius of the central object in cm
        mdot            the mass accretion rate of the accretion disk in g s^-1

    Returns
        teff            the effective temperature of the accretion disk at a point R in K
    """

    if r < rin:
        print("ss_disk_teff.t_eff: r < rin")
        return 0

    teff = (3 * G * m * mdot) / (8 * np.pi * r ** 3 * STEFAN_BOLTZMANN)
    teff *= (1 - (rin / r) ** 0.5)
    teff = teff ** 0.25

    return teff


def bb_disk_spectrum(wmin: float, wmax: float, m: float, mdot: float, rin: float, rout: float, nlamda: int = 500,
                     nrings: int = 500, verbose: bool = False) -> np.array:
    """
    Calculate the spectrum of an accretion disk of inner radius rin and outer radius rout and accretion rate mdot.

    Parameters
        wmin, wmax      wavelength limits for spectrum in Angstroms
        m               mass of the central object in msol
        mdot            mass accretion rate of the disk in msol/yr
        rin, rout       inner and outer radius of the accretion disk (assume rin = r_co) in cm
        nlamda          [optional] the number of points on the wavelength grid
        nrings          [optional] the number of disk annuli
        verbose         [optional] enable verbose output

    Returns
        disk_spec       2 x nlamda array - the disk spectrum
                            - disk_spec[:, 0] the wavelength range in Angstrom
                            - disk_spec[:, 1] the monochromatic flux ergs s^-1 cm^-2 A^-1
    """

    # Convert into CGS units
    m *= MSOL
    mdot *= MSOL_PER_YEAR

    # Crate the radius and wavelength grid
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
        disk_spec[:, 1] += planck_lambda(lamda, teff) * area

    if verbose:
        print("tde_plot_spec.bb_disk_spectrum: max effective disk temperature: {} K".format(teff_max))

    return disk_spec


def main():
    """
    Nothing happens :-)

    Parameters
        None

    Returns
        None
    """

    wmin = 0.1
    wmax = 1e4
    rin = 7e8
    rout = 2.4e10
    m = 0.8
    mdot = 1e-8

    disk_spec = bb_disk_spectrum(wmin, wmax, m, mdot, rin, rout, verbose=True)

    flux = disk_spec[:, 1]
    lamda = disk_spec[:, 0]

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.semilogy(lamda, flux)
    ax.set_xlabel("Wavelength, $\AA$")
    ax.set_ylabel("Flux")
    plt.show()

    return


if __name__ == "__main__":
    main()
