#!/usr/bin/env python 
# -*- coding: utf-8 -*-

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

	if type(wavelength) != float:
		try:
			wavelength = float(wavelength)
		except ValueError:
			print("Could not convert the provided wavelength to a float!!! :-(")
			sys.exit(1)

	frequency = C / (wavelength * ANGSTROM)

	print("Photon of wavelength {:.0f} A has frequency {:e} Hz".format(wavelength, frequency))

	return frequency


if __name__ == "__main__":

	if len(sys.argv) != 2:
		print("Provide a wavelength pls")
		sys.exit(1)
	else:
		convert_to_frequency(sys.argv[1])
