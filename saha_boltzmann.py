#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to create a qualitative plot of the Saha-Boltzmann
ionisation equation. 

This is modelling the ion density ratio of He I and He II by default.
"""


import numpy as np
from matplotlib import pyplot as plt
from consts import *


def saha_equation(ne: float, g_upper: float, g_lower: float, energy_upper: float, energy_lower: float, 
                  temperature: float):
    """
    Calculate the ratio of n_i+1 / n_i, using the Saha-Boltzman equation.

    Parameters
    ----------
    ne: float
        The electron density of the plasma, in cm^-2.
    g_upper: float
        The statistical weight of the upper level.
    g_lower: float
        The statistical weight of the lower level.
    energy_upper: float
        The ionisation potential of the upper level, in ergs.
    energy_lower: upper
        The ionisation potential of the lwoer level, in ergs.
    temperature: float
        The temperature of the plasma in K.

    Returns
    -------
    N_i+1 / N_i: float
        The ratio of the population of the upper ionisation and ground state
        of the atom.
    """

    gratio = 2 * g_upper / g_lower
    saha = ((2 * PI * MELEC * BOLTZMANN * temperature) / H ** 2) ** (3 / 2)

    return saha * gratio * np.exp(-(energy_upper - energy_lower) / (BOLTZMANN * temperature)) / ne


def helium_example():
    """
    A Helium example for the Saha-Boltzman ionisation equation.
    """

    gi_upper = 2
    gi_lower = 1
    energy_upper = 54.418 * EV2ERGS
    energy_lower = 24.588 * EV2ERGS
    temperature = np.linspace(11000, 15000, 5000)
    densities = np.linspace(1e18, 1e21, 5000)

    temp_ratio = saha_equation(1e10, gi_upper, gi_lower, energy_upper, energy_lower, temperature)
    dens_ratio = saha_equation(densities, gi_upper, gi_lower, energy_upper, energy_lower, 40000)

    nplots = 2
    fig, ax = plt.subplots(1, nplots, figsize=(12, 5))

    ax[0].semilogy(temperature, temp_ratio)
    ax[0].set_xlabel("Temperature [K]")
    ax[0].set_xlim(11000, 15000)
    ax[0].text(0.2, 0.05, "Ionisation grows with temperature", fontsize=13, transform=ax[0].transAxes)
    ax[0].set_title("Ionisation ratio as a function of temperature")

    ax[1].loglog(densities, dens_ratio)
    ax[1].set_xlabel(r"Electron Density [cm$^{\circ}$]")
    ax[1].set_xlim(1e18, 1e21)
    ax[1].text(0.05, 0.05, "Ionistation drops with electron density\n(recombination)", fontsize=13, 
            transform=ax[1].transAxes)
    ax[1].set_title("Ionisation ratio as a function of electron density")

    for i in range(nplots):
        ax[i].set_ylabel(r"$N_{i+1} / N_{i}$")

    plt.show()

    return

if __name__ == "__main__":
    helium_example()
