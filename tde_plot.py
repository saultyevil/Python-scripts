#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to plot a more TDE inclined spectrum. To do this,
the script includes options to over plot various TDE spectrum onto the plot, as well
as additionally plotting emission and absorption line IDs for common lines we
would expect in TDE spectrum with high velocity wind outflows.
"""

from sys import exit
import argparse
import tde_spectra
import numpy as np
from consts import *
from typing import Tuple
from matplotlib import pyplot as plt
from PyPython import SpectrumUtils
from PyPython import PythonUtils as Utils
from PyPython import Quotes

SMOOTH = 5
WMIN = 1000
WMAX = 3000
TDE_OBJ = "iPTF15af"
DEFAULT_DIST = 100 * PARSEC
FTYPE = "png"
PLOT_LINE_IDS = True
VERBOSITY = False


def get_tde_spectrum() -> Tuple[np.array, float, str]:
    """
    Return an array containing the TDE spectrum as well as returning the distance
    of the object and a string containing the reference for the observation.

    Returns
    -------
    tde_spec: np.array[float]
        The spectrum for the TDE
            - Column 0 is wavelength,
            - Column 1 is flux,
            - Column 2 is error (sometimes)
    observse_dist: float
        The distance of the TDE given in cm
    reference: str
        The journal reference for the observation
    """

    if TDE_OBJ == "iPTF15af" or TDE_OBJ == "iptf15af":
        z = 0.07897
        reference = "Blagorodnova et al. (2019)"
        observe_dist = 350 * 1e6 * PARSEC
        tde_spec = tde_spectra.iptf15af_spec(SMOOTH, VERBOSITY)
    elif TDE_OBJ == "ASASSN14li" or TDE_OBJ == "asassn14li":
        z = 0.02058
        reference = "Cenko et al. (2016)"
        observe_dist = 90 * 1e6 * PARSEC
        tde_spec = tde_spectra.asassn14li_spec(SMOOTH, VERBOSITY)
    else:
        print("tde_spec_plot.get_tde_spectrum: can't find spectrum for object {}".format(TDE_OBJ))
        exit(1)

    tde_spec[:, 0] /= (z + 1)

    return tde_spec, observe_dist, reference


def plot_single_inclination(root: str, inc: str) -> None:
    """
    Plot the spectrum for a single inclination for a Python TDE simulation.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    inc: str
        The inclination angle to plot the spectrum for
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()
        ax.semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ + " " + reference)

    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)

    try:
        spectrum = SpectrumUtils.read_spec(spec_file)
    except IOError:
        print("tde_spec_spec.spec_plot_one: could not open file {}".format(spec_file))
        return

    wavelength = spectrum["Lambda"].values.astype(float)
    try:
        raw_flux = spectrum[inc].values.astype(float)
    except KeyError:
        print("tde_spec_plot.spec_plot_one: could not find inclination angle of {}".format(inc))
        return
    flux = SpectrumUtils.smooth(raw_flux, SMOOTH)
    flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
    ax.semilogy(wavelength, flux, label="Model at i = {}".format(inc) + r"$^{\circ}$")

    ymin, ymax = SpectrumUtils.ylims(wavelength, flux, WMIN, WMAX)
    ax.set_xlim(WMIN, WMAX)
    ax.set_ylim(ymin, ymax)
    if PLOT_LINE_IDS:
        SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.legend(loc="best")

    fig.suptitle(root)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_{}_spectra.{}".format(root, inc, FTYPE))
    plt.close()

    return


def plot_all_inclinations(root: str) -> None:
    """
    Plot the spectrum for each inclination angle for a Python TDE simulation on
    a single plot.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    """

    tde = None
    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()

    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)

    try:
        spectrum = SpectrumUtils.read_spec(spec_file)
    except IOError:
        print("tde_spec_plot.spec_plot_multiple: could not open file {}".format(spec_file))
        return

    inclinations = spectrum.columns.values[9:]
    n_spec = len(inclinations)
    wavelength = spectrum["Lambda"].values.astype(float)
    nrows, ncols = Utils.subplot_dims(n_spec)
    fig, ax = plt.subplots(nrows, ncols, figsize=(17, 12), squeeze=False, sharex="col")

    index = 0
    for i in range(nrows):
        for j in range(ncols):
            if index > n_spec - 1:
                break
            if i == nrows - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")

            if TDE_OBJ:
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)

            inc = inclinations[index]
            raw_flux = spectrum[inc].values.astype(float)
            flux = SpectrumUtils.smooth(raw_flux, SMOOTH)
            flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
            ax[i, j].semilogy(wavelength, flux, label=r"$i$ = {}".format(inc) + r"$^{\circ}$")

            ymin, ymax = SpectrumUtils.ylims(wavelength, flux, WMIN, WMAX)
            ax[i, j].set_ylim(ymin, ymax)
            ax[i, j].set_xlim(WMIN, WMAX)
            if PLOT_LINE_IDS:
                ax[i, j] = SpectrumUtils.plot_line_ids(ax[i, j], SpectrumUtils.common_lines())
            ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
            ax[i, j].legend(loc="best")
            index += 1

    fig.suptitle(root)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_spectra.{}".format(root, FTYPE))
    plt.close()

    return


def plot_model_comparisons(name: str, inc: str = None):
    """
    Create a comparison plot similar to spec_plot for multiple spectra which are
    search for recursively in the working directory.

    Parameters
    ----------
    name: str
        The base name of the output plot
    inc: str, optional
        The inclination angle to be plotted
    """

    spec_files = SpectrumUtils.find_specs()
    if len(spec_files) == 0:
        print("No spec files found")
        print("\n--------------------------")
        return

    print("Creating spectra with the following simulations:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    tde = None
    observe_dist = 1
    if TDE_OBJ and TDE_OBJ != "None".lower():
        tde, observe_dist, ref = get_tde_spectrum()

    if inc:
        if inc.isdigit() is False:
            print("Provided inclination angle is not a number: {}".format(inc))
            exit(1)
        inclination = inc
        size = (12, 5)
        n_specs = 1
    else:
        inclination = []
        for f in spec_files:
            spec = SpectrumUtils.read_spec(f)
            inclination += list(spec.columns.values[9:])
        inclination = sorted(list(dict.fromkeys(inclination)))
        size = (20, 20)
        n_specs = len(inclination)

    nrows, ncols = Utils.subplot_dims(n_specs)
    fig, ax = plt.subplots(nrows, ncols, figsize=size, squeeze=False, sharex="col")

    index = 0
    for i in range(nrows):
        for j in range(ncols):
            if index > n_specs - 1:
                break
            if i == nrows - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            if j == 0:
                ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
            if TDE_OBJ and TDE_OBJ != "None".lower():
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)

            ymin = +1e99
            ymax = -1e99

            # If a specific inclination angle has been provided, then use this
            if inc:
                ii = inc
            else:
                ii = inclination[index]

            # Loop over each spectrum from each simulation
            for file in spec_files:
                root, dir = Utils.split_root_directory(file)

                spectrum = SpectrumUtils.read_spec(file)
                wavelength = spectrum["Lambda"].values.astype(float)

                try:
                    raw_flux = spectrum[ii].values.astype(float)
                except KeyError:
                    continue
                flux = SpectrumUtils.smooth(raw_flux, SMOOTH)
                flux *= DEFAULT_DIST ** 2 / observe_dist ** 2

                ax[i, j].semilogy(wavelength, flux, label=r"{}/{}: $i$: {}".format(dir, root, ii) + r"$^{\circ}$")

                tymin, tymax = SpectrumUtils.ylims(wavelength, flux, WMIN, WMAX)
                if tymin is not None and tymin < ymin:
                    ymin = tymin
                if tymax is not None and tymax > ymax:
                    ymax = tymax

            if ymin == +1e99:
                ymin = None
            if ymax == -1e99:
                ymax = None

            ax[i, j].set_xlim(WMIN, WMAX)
            ax[i, j].set_ylim(ymin, ymax)
            if PLOT_LINE_IDS:
                SpectrumUtils.plot_line_ids(ax[i, j], SpectrumUtils.common_lines())
            ax[i, j].legend(loc="best")

            index += 1

    fig.suptitle("Model Comparison")
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    fname = "{}_comparison".format(name)
    if inc:
        fname += "_i{}".format(inc)
    plt.savefig("{}.{}".format(fname, FTYPE))
    plt.close()

    return


def main() -> None:
    """
    Main controlling function of the script. Also contains all of the code which
    parses the command line for run time parameters.
    """

    global SMOOTH
    global FTYPE
    global VERBOSITY
    global WMIN
    global WMAX
    global TDE_OBJ
    global PLOT_LINE_IDS

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=str, help="The root name of the Python simulation to plot.")
    p.add_argument("-i", "--inclination", type=str, action="store", help="Plot a single inclination angle.")
    p.add_argument("-s", "--smooth", type=int, action="store", help="Window size for the boxcar smooth - default is 15.")
    p.add_argument("-l", "--lineid", action="store_true", help="Plot line IDs for common absorption and emission lines.")
    p.add_argument("-c", "--comparison", action="store_true", help="Plot a comparison of spectra from multiple runs.")
    p.add_argument("-wmin", type=float, action="store", help="Shortest wavelength to plot.")
    p.add_argument("-wmax", type=float, action="store", help="Longest wavelength to plot.")
    p.add_argument("-tde", type=str,  action="store", help="The name of the TDE object to overplot.")
    p.add_argument("-ftype", type=str,action="store", help="The file type of the output image.")
    p.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")
    args = p.parse_args()

    root = args.root

    if args.smooth:
        SMOOTH = args.smooth
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax
    if args.tde:
        TDE_OBJ = args.tde
    if args.ftype:
        FTYPE = args.ftype
    if args.verbose:
        VERBOSE = True
    if args.lineid:
        PLOT_LINE_IDS = False

    print("--------------------------\n")

    Quotes.random_quote()

    if args.inclination and not args.comparison:
        print("Plotting {} spectrum for inclination {}".format(root, args.inclination))
        plot_single_inclination(root, args.inclination)
    elif args.comparison:
        print("Plotting comparison spectrum")
        if args.inclination:
            plot_model_comparisons(root, args.inclination)
        else:
            plot_model_comparisons(root)
    else:
        print("Plotting {} spectra for all inclination angles".format(root))
        plot_all_inclinations(root)

    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
