#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
from consts import *
from matplotlib import pyplot as plt


def t_eff(r: np.ndarray, m: float, mdot: float, rco: float) -> float:
    """
    Return the effective temperature at a position R for an alpha disc.

    Parameters
    ----------
    r: np.ndarry
        The radial positions to calculate t_eff at, in cm.
    m: float
        The black hole mass in g.
    mdot: float
        The accretion rate in g/s.
    rco: float
        The radius of the central object in cm.

    Returns
    -------

    """

    return ((3 * G * m * mdot) / (8 * np.pi * r ** 3 * STEFAN_BOLTZMANN) * (1 - (rco / r) ** 0.5)) ** 0.25


if __name__ == "__main__":
    r = np.logspace(np.log10(8.9e11), np.log10(1e15), 500)
    mbh = 1e6 * MSOL
    mdot = 1e-2 * MSOL_PER_YEAR
    rco = 8.86e11
    teff = t_eff(r, mbh, mdot, rco)

    print("Max temperature in disc {:e}".format(np.max(teff)))

    plt.loglog(r, teff)
    plt.xlabel("Radius (cm)")
    plt.ylabel("Effective Temperature (k)")
    plt.show()
