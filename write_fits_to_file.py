#!/usr/bin/env python3

import argparse
from astropy.io import fits


def parse():
    """
    Parse run time options from the command line

    Parameters
    ----------
    None

    Returns
    -------
    input_filename: str
        The input fits file to be processes.
    output_filename: str
        The filename of the output dat file.
    """

    input_filename = "spec.txt"
    output_filename = "aa"
    return input_filename, output_filename


def process_fits(input_filename, output_filename):
    """
    Open and process the fits file

    Parameters
    ----------

    Returns
    -------
    None

    """

    with fits.open(input_filename) as fitsfile:
        fitsfile.info()
        fitsdata = fitsfile[1].data
        wav = fitsdata['WAVELENGTH']
        flux = fitsdata['FLUX']
    with open("{}.dat".format(output_filename), "w") as f:
        f.write("Wavelength\tFlux\n")
        for i in range(len(flux)):
            f.write("{}\t{}\n".format(wav[i], flux[i]))

    return


def main():
    """
    The main steering function

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    input_filename, output_filename = parse()
    process_fits(input_filename, output_filename)

    return


if __name__ == "__main__":
    main()
