#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from constants import *


def convert_to_wavelength(energy: float):
    """
    Convert a photon energy given in eV into Angstroms.

    E = hc / lambda

    Parameters
    ----------
    energy: float
         The energy of the photon in eV

    Returns
    -------
    The wavelength of the photon in Angstroms
    """

    wl = (H * C / (energy * EV2ERGS)) / ANGSTROM
    print("Wavelength for a photon of energy {:.2e} eV is {:.2f} Angstroms".format(energy, wl))

    return wl


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Please provide a photon energy in eV")
    else:
        wavelength = convert_to_wavelength(float(sys.argv[1]))
