#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from consts import *


def convert_to_wavelength(phot_energy: float):
    """
    Convert a photon energy given in eV into Angstroms.

    E = hc / lambda

    Parameters
    ----------
    phot_energy: float
         The energy of the photon in eV

    Returns
    -------
    The wavelength of the photon in Angstroms
    """

    wavelength = (H * C / (phot_energy * EV2ERGS)) / ANGSTROM

    print("Wavelength for a photon of energy {:.2e} eV is {:.2f} Angstroms".format(
            phot_energy, wavelength))

    return wavelength


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Please provide a photon energy in eV")
        sys.exit(1)
    else:
        phot_energy = float(sys.argv[1])
        wavelength = convert_to_wavelength(phot_energy)
