#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import argparse
import sys
from astropy.io import fits
from astropy.table import Table
from typing import Tuple, List

from matplotlib import pyplot as plt
import py_plot_util as ppu


def parse_command_line() -> Tuple[int, List[str]]:
    """
    Parse the command line for various file names.

    Returns
    -------
    fnames: List[str]
        The fnames provided on the command line.
    """

    p = argparse.ArgumentParser(description="Combined multiple spectrum files into one spectrum")
    p.add_argument("fnames", nargs="+",  help="The file names of the spectra to combine")
    args = p.parse_args()

    return len(args.fnames), args.fnames


def read_fits_file(fname: str) -> np.ndarray:
    """
    Read in a 1D spectrum from a FITS file. This function is fucking shite, but
    I got fed up of the fucking god awful documentation which Astropy provides
    for their shitty and vile package. f.data is most definitely not a fucking
    numpy array as the cretin who wrote the documentation claims. It's actually
    some god awful internal class which should be converted into a cunting Astropy
    table for '''convenience''' :^^^^^^^^^^^^^^^^^^^^^^^^^^^^). Fuck you.

    Maybe I should just learn how Astropy works rather than getting angry at it.

    Parameters
    ----------
    fname: str
        The file name of the FITS file to open.

    Returns
    -------
    spec: np.ndarray
        A Numpy array containing the wavelength (0), flux (1) and error (2)
        for the spectrum which was read in.
    """

    f = fits.open(fname)
    f.info()

    f = f["SCI"]
    d = Table(f.data)
    wl = d["WAVELENGTH"].tolist()[0]
    fl = d["FLUX"].tolist()[0]
    er = d["ERROR"].tolist()[0]

    spec = np.zeros((len(wl), 3))
    spec[:, 0] = wl
    spec[:, 1] = fl
    spec[:, 2] = er

    return spec


def combine(fnames: List[str], nfiles: int):
    """
    Combine multiple spectra into one spectrum. It's assumed that the file
    extensions are .fits.

    Parameters
    ----------
    fnames: List[str]
        A list containing the file names of the spectra to combine.
    nfiles: int
        The number of spectra to combine.

    Returns
    -------

    """

    specs = []

    print("Reading in the various spec files:")
    for i in range(nfiles):
        specs.append(read_fits_file(fnames[i]))

    for i in range(nfiles):
        name = "{}.txt".format(fnames[i])
        np.savetxt(name, specs[i])

    # Sort the spectra by their smallest wavelengths
    specs.sort(key=lambda x: x[0, 0])

    print("Combining the spectra:")

    tmp = []

    for i in range(nfiles - 1):
        w1 = specs[i][0, 0]
        w2 = specs[i + 1][0, 0]
        print(w2)
        nele = len(specs[i][:, 0])
        for j in range(nele):
            if w2 < specs[i][j, 0]:
                print(j, specs[i][j - 1, 0], specs[i][j, 0], w2)
                break

    tlen = j + len(specs[1][:, 0])
    spec = np.zeros((tlen, 3))

    spec[:j, 0] = specs[0][:j, 0]
    spec[:j, 1] = specs[0][:j, 1]
    spec[:j, 2] = specs[0][:j, 2]
    spec[j:, 0] = specs[1][:, 0]
    spec[j:, 1] = specs[1][:, 1]
    spec[j:, 2] = specs[1][:, 2]

    return spec


def main():
    """
    Main function.
    """

    print("Combining the current spectra:")
    nfiles, fnames = parse_command_line()
    if len(fnames) < 2:
        print("Not enough filenames were provided, nothing to combine.")
        sys.exit(1)
    fnames = list(dict.fromkeys(fnames))
    if len(fnames) < 2:
        print("Not enough filenames provided after removing any duplicates")
        sys.exit(1)
    for i in range(nfiles):
        print("  - {}".format(fnames[i]))
    print()

    spec = combine(fnames, nfiles)

    print(spec[:, 1])

    plt.semilogy(spec[:, 0], ppu.smooth_1d_array(spec[:, 1], 10))
    # plt.semilogy(spec[:, 0], spec[:, 1])
    plt.show()

    return


if __name__ == "__main__":
    main()
