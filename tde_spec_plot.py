#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to plot a more TDE inclined spectrum. To do this,
the script includes options to over plot various TDE spectrum onto the plot, as well
as additionally plotting emission and absorption line IDs for common lines we
would expect in TDE spectrum with high velocity wind outflows.
"""

import argparse
import tde_util
import numpy as np
import py_plot_util
from consts import *
from typing import Tuple
from matplotlib import pyplot as plt

SMOOTH = 15
WMIN = 800
WMAX = 3500
TDE_OBJ = "iPTF15af"
DEFAULT_DIST = 100 * PARSEC
FTYPE = "png"
PLOT_LINE_IDS = True
VERBOSE = False


def get_tde_spectrum() -> Tuple[np.array, float, str]:
    """
    Return an array containing the TDE spectrum as well as returning the distance
    of the object and a string containing the reference for the observation.

    Returns
    -------
    tde_spec            np.array
                        The spectrum for the TDE
                            - Column 0 is wavelength,
                            - Column 1 is flux,
                            - Column 2 is error (sometimes)
    observse_dist       float
                        The distance of the TDE given in cm
    reference           str
                        The journal reference for the observation
    """

    if TDE_OBJ == "iPTF15af" or TDE_OBJ == "iptf15af":
        z = 0.07897
        reference = "Blagorodnova et al. (2019)"
        observe_dist = 350 * 1e6 * PARSEC
        tde_spec = tde_util.iPTF15af_spec(SMOOTH, VERBOSE)
    elif TDE_OBJ == "ASSASN14li" or TDE_OBJ == "assasn14li":
        z = 0.02058
        reference = "Cenko et al. (2016)"
        observe_dist = 90 * 1e6 * PARSEC
        tde_spec = tde_util.ASSASN14li_spec(SMOOTH, VERBOSE)
    else:
        print("tde_spec_plot.get_tde_spectrum: can't find spectrum for object {}".format(TDE_OBJ))
        exit(1)

    tde_spec[:, 0] /= (z + 1)

    return tde_spec, observe_dist, reference


def spec_plot_inclination(root: str, inc: str) -> None:
    """
    Plot the spectrum for a single inclination for a Python TDE simulation.

    Parameters
    ----------
    root            str
                    The root name of the Python simulation
    inc             str
                    The inclination angle to plot the spectrum for

    Returns
    -------
    None
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()
        ax.semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ + " " + reference)

    # Get the spectrum from file for the TDE model
    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)
    try:
        spectrum = py_plot_util.read_spec_file(spec_file, pandas_table=True)
    except IOError:
        print("tde_spec_spec.spec_plot_one: could not open file {}".format(spec_file))
        return
    wavelength = spectrum["Lambda"].values.astype(float)
    try:
        raw_flux = spectrum[inc].values.astype(float)
    except KeyError:
        print("tde_spec_plot.spec_plot_one: could not find inclination angle of {}".format(inc))
        return
    flux = py_plot_util.smooth_1d_array(raw_flux, SMOOTH, VERBOSE)
    flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
    ymax, ymin = py_plot_util.get_ylimits(wavelength, flux, WMIN, WMAX)

    # Plot the TDE spectrum for the specific inclination angle
    ax.semilogy(wavelength, flux, label="i = {}".format(inc) + r"$^{\circ}$")
    ax.set_xlim(WMIN, WMAX)
    ax.set_ylim(ymin, ymax)
    if PLOT_LINE_IDS:
        py_plot_util.plot_line_ids(ax, py_plot_util.get_common_line_ids())
    ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.legend(loc="lower right")

    # Finishing touches to the plot, including title and tight layout
    fig.suptitle(root)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_{}_spectra.{}".format(root, inc, FTYPE))
    plt.close()

    return


def spec_plot_multiple(root: str) -> None:
    """
    Plot the spectrum for each inclination angle for a Python TDE simulation on
    a single plot.

    Parameters
    ----------
    root            str
                    The root name of the Python simulation

    Returns
    -------
    None
    """

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()

    # Read the spectrum for the TDE into memory
    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)
    try:
        spectrum = py_plot_util.read_spec_file(spec_file, pandas_table=True)
    except IOError:
        print("tde_spec_plot.spec_plot_multiple: could not open file {}".format(spec_file))
        return
    wavelength = spectrum["Lambda"].values.astype(float)
    inclinations = py_plot_util.spec_inclinations_pandas(spectrum)

    # Figure out the shape of the subplot grid and create the plotting object
    n_spec = len(inclinations)
    nrows, ncols = py_plot_util.subplot_dims(n_spec)
    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 8), squeeze=False, sharex="col")

    index = 0
    for i in range(nrows):
        for j in range(ncols):
            if index > n_spec - 1:
                break
            # This ensures that the wavelength label will only be put on the bottom plots
            if i == nrows - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            if TDE_OBJ:
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)

            # Extract the flux from the spectrum and scale for distance for a specific inclination
            inc = inclinations[index]
            raw_flux = spectrum[inc].values.astype(float)
            flux = py_plot_util.smooth_1d_array(raw_flux, SMOOTH, VERBOSE)
            flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
            ymax, ymin = py_plot_util.get_ylimits(wavelength, flux, WMIN, WMAX)

            # Plot the TDE spectrum
            ax[i, j].semilogy(wavelength, flux, label=r"$i$ = {}".format(inc) + r"$^{\circ}$")
            ax[i, j].set_ylim(ymin, ymax)
            ax[i, j].set_xlim(WMIN, WMAX)
            if PLOT_LINE_IDS:
                ax[i, j] = py_plot_util.plot_line_ids(ax[i, j], py_plot_util.get_common_line_ids())
            ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
            ax[i, j].legend(loc="upper right")

            # Increment index for inclination array
            index += 1

    # Finishing touches to the plot, including title and tight layout
    fig.suptitle(root)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_spectra.{}".format(root, FTYPE))
    plt.close()

    return


def spec_plot_multiple_comparison(name: str, inc: str = None):
    """
    Create a comparison plot similar to spec_plot for multiple spectra which are
    search for recursively in the working directory.

    Parameters
    ----------
    name                str
                        The base name of the output plot
    inc                 str
                        The inclination angle

    Returns
    -------
    None
    """

    print("USE WITH CAUTION. THIS SHITTY THING DOESN'T WORK SOMETIMES.")

    spec_files = py_plot_util.find_spec_files()
    if len(spec_files) == 0:
        print("No spec files found")
        print("\n--------------------------")
        return

    print("Creating spectra with the following simulations:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, ref = get_tde_spectrum()

    # Figure out the inclinations in the spec files for the simulations
    if inc:
        inclination = [inc]
        size = (12, 8)
    else:
        inclination = []
        for f in spec_files:
            spec = py_plot_util.read_spec_file(f, pandas_table=True)
            inclination += py_plot_util.spec_inclinations_pandas(spec)
        inclination = sorted(list(dict.fromkeys(inclination)))
        size = (20, 12)

    n_specs = len(inclination)
    nrows, ncols = py_plot_util.subplot_dims(n_specs)
    fig, ax = plt.subplots(nrows, ncols, figsize=size, squeeze=False, sharex="col")

    index = 0
    for i in range(nrows):
        for j in range(ncols):
            if i == nrows - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            if index > n_specs - 1:
                break

            if TDE_OBJ:
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)

            # If a specific inclination angle has been provided, then use this
            if inc:
                ii = inc
            else:
                ii = inclination[index]

            ymin = +1e99
            ymax = -1e99
            # Loop over each spectrum from each simulation
            for file in spec_files:
                # Get the spectrum, flux and scale the flux for distance
                root, dir = py_plot_util.parse_root_name_and_path(file)

                spectrum = py_plot_util.read_spec_file(file, pandas_table=True)
                wavelength = spectrum["Lambda"].values.astype(float)
                try:
                    raw_flux = spectrum[ii].values.astype(float)
                except KeyError:
                    continue
                flux = py_plot_util.smooth_1d_array(raw_flux, SMOOTH, VERBOSE)
                flux *= DEFAULT_DIST ** 2 / observe_dist ** 2

                # Plot the spectrum for a model
                ax[i, j].semilogy(wavelength, flux, label=r"$i$: {}".format(ii) + r"$^{\circ}$: " + dir)

                tmax, tmin = py_plot_util.get_ylimits(wavelength, flux, WMIN, WMAX)
                if tmax > ymax:
                    ymax = tmax
                if tmin < ymin:
                    ymin = tmin

            # Finishing touches to plot
            ax[i, j].set_xlim(WMIN, WMAX)
            ax[i, j].set_ylim(ymin, ymax)
            if PLOT_LINE_IDS:
                py_plot_util.plot_line_ids(ax[i, j], py_plot_util.get_common_line_ids())
            ax[i, j].legend(loc="upper right")

            # Increment inclination index
            index += 1

    fig.suptitle("Model Comparison")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_comparison.{}".format(name, FTYPE))
    plt.close()

    return


def main() -> None:
    """
    Main controlling function of the script.

    Returns
    -------
    None
    """

    global SMOOTH
    global FTYPE
    global VERBOSE
    global WMIN
    global WMAX
    global TDE_OBJ
    global PLOT_LINE_IDS

    p = argparse.ArgumentParser(description="Plot a more TDE inclined spectrum or spectra.")
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
        PLOT_LINE_IDS = True

    print("--------------------------\n")

    if args.inclination and not args.comparison:
        print("Plotting {} spectrum for inclination {}".format(root, args.inclination))
        spec_plot_inclination(root, args.inclination)
    elif args.comparison:
        print("Plotting comparison spectrum")
        if args.inclination:
            spec_plot_multiple_comparison(root, args.inclination)
        else:
            spec_plot_multiple_comparison(root)
    else:
        print("Plotting {} spectra for all inclination angles".format(root))
        spec_plot_multiple(root)

    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
