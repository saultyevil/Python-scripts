#!/usr/bin/env python
# -*- coding: utf-8 -*-


from sys import exit
import numpy as np
import argparse
import astropy.units as u
from dust_extinction.parameter_averages import F99, CCM89
from matplotlib import pyplot as plt
import py_plot_util as ppu


"""
iPTF16fnl:
    Rv = 3.1
    E(B-V) = 0.0777
AT2018zr:
    Rv = 3.06
    E(B-V) = 0.0404
"""


SMOOTH = 5


def parse_arguments():
    """
    Parse run time parameters from the command line.

    Returns
    -------
    fname: str
        The filename of the spectrum to deredden.
    Rv: float
        The selective extinction parameter
    Ebv: float
        The colour excess parameter
    skip: int
        The number of lines to skip before reading in the spectrum
    ext_curve: str
        The extinction curve to use
    """

    p = argparse.ArgumentParser(description="Deredden spectra given Av, Rv and E(B-V)")
    p.add_argument("fname", type=str, help="The filename of the spectrum to reredden")
    p.add_argument("Rv", type=float, help="The selective extinction parameter")
    p.add_argument("Ebv", type=float, help="The colour excess parameter")
    p.add_argument("-ext_curve", type=str, help="The extinction curve to use")
    p.add_argument("-skiplines", type=int, help="Skip a certain number of lines before reading in the spectrum")
    args = p.parse_args()

    skip = 0
    if args.skiplines:
        skip = args.skiplines
    ext_curve = "CCM89"
    if args.ext_curve:
        ext_curve = args.ext_curve

    return args.fname, args.Rv, args.Ebv, skip, ext_curve


def deredden_spectrum(spec: np.ndarray, Rv: float, Ebv: float, ext_curve: str = "CCM89"):
    """
    Deredden an input spectrum given the selective extinction and colour exccess.
    By default, this uses a CCM89 (Cardelli et al. 1989) extinction curve.

    Parameters
    ----------
    spec: np.ndarray
        The spectrum to be dereddened.
    Rv: float
        The selective extinction parameter.
    Ebv: float
        The colour excess parameter.
    ext_curve: str, optional
        The type of extinction curve to fit with.

    Returns
    -------
    dered: np.ndarray
        The dereddened input spectrum.
    """

    dered = np.copy(spec)

    wl = dered[:, 0] * u.angstrom

    if ext_curve == "CCM89":
        extinc = CCM89(Rv=Rv)
    elif ext_curve == "F99":
        extinc = F99(Rv=Rv)
    else:
        print("Extinction curve {} unknown".format(ext_curve))
        exit(1)

    dered[:, 1] /= extinc.extinguish(wl, Ebv=Ebv)

    return dered


def main():
    """
    Main control function for the script.
    """

    fname, Rv, Ebv, skip, ext_curve = parse_arguments()
    spec = np.loadtxt(fname, skiprows=skip)
    dered = deredden_spectrum(spec, Rv, Ebv, ext_curve)
    np.savetxt("{}_dered.txt".format(fname), dered)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.semilogy(spec[:, 0], ppu.smooth_1d_array(spec[:, 1], SMOOTH), label="Original spectrum")
    ax.semilogy(dered[:, 0], ppu.smooth_1d_array(dered[:, 1], SMOOTH), label="The dereddened spectrum")
    ax.set_xlabel(r"Observed Wavelength [$\AA$]")
    ax.set_ylabel(r"Flux")
    ax.legend()
    plt.show()

    return


if __name__ == "__main__":
    main()
