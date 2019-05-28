#!/usr/bin/env python 

import sys
from consts import *


def convert_to_frequency(wavelength):
	"""

	Parameters
	----------
	wavelength          float
						The wavelength in Angstroms

	Returns
	-------
	The frequency of a photon of the given wavelength.
	"""

	frequency = C / (wavelength * ANGSTROM)
	print("Photon of wavelength {:.0f} A has frequency {} Hz".format(wavelength, frequency))

	return frequency


if __name__ == "__main__":

	if len(sys.argv) != 2:
		print("Provide a wavelength pls")
		sys.exit(1)
	else:
		convert_to_frequency(sys.argv[1])
