#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to plot a more TDE inclined spectrum. To do this,
the script includes options to over plot various TDE spectrum onto the plot, as well
as additionally plotting emission and absorption line IDs for common lines we
would expect in TDE spectrum with high velocity wind outflows.
"""

import argparse
import numpy as np
import pandas as pd
import py_plot_util
from consts import *
from typing import Tuple
from matplotlib import pyplot as plt

SMOOTH = 15
WMIN = 800
WMAX = 3000
TDE_OBJ = None
DEFAULT_DIST = 100 * PARSEC
FTYPE = "png"
PLOT_LINE_IDS = False
VERBOSE = False


def spectrum_inclinations(spectrum: pd.DataFrame) -> list:
    """
    Return a list of the viewing angles within a spectrum

    Parameters
    ----------
    spectrum            Pandas DataFrame
                        The spectrum to extract the viewing angles from

    Returns
    -------
    viewing_angles      list
                        The inclinations within the spectrum file
    """

    viewing_angles = []

    if type(spectrum) != pd.DataFrame:
        print("viewing_angles_in_spec: spectrum is not a pandas dataframe")
        exit(1)

    column_headers = spectrum.columns.values
    for header in column_headers:
        if header.isdigit():
            viewing_angles.append(header)

    return viewing_angles


def subplot_dims(n_plots: int) -> tuple:
    """
    Determine the dimensions for a plot with multiple subplot panels. Two columns
     of subplots will always be used.

    TODO: if n_plots is large, use 3 columns instead

    Parameters
    ----------
    n_plots         int
                    The number of subplots which will be plotted

    Returns
    -------
    dims            tuple
                    The dimensions of the subplots returned as (nrows, ncols)
    """

    if n_plots > 2:
        ncols = 2
        nrows = (1 + n_plots) // ncols
    else:
        ncols = 1
        nrows = n_plots

    dims = (nrows, ncols)

    return dims


def plot_line_ids(ax: plt.axes, lines: dict, rotation: str = "horizontal") -> plt.axes:
    """
    Plot major absorption and emission line IDs onto a spectrum.

    Note, this should be used after setting the x limits on a figure.

    Parameters
    ----------
    ax              plt.axes
                    The plot object to add line IDs to
    lines           dict
                    A dictionary containing the line name and wavelength in Angstroms (ordered by wavelength)
    rotation        str, optional
                    Vertical or horizontal rotation for text ids

    Returns
    -------
    ax              plt.axes
                    The plot object now with lines IDs :-)
    """

    xlims = ax.get_xlim()

    for i, l in enumerate(lines):
        # check to make sure that a line isn't off the actual axis of the plot
        if float(lines[l]) > xlims[1]:
            continue
        if float(lines[l]) < xlims[0]:
            continue
        # plot the text in an alternating up-down pattern
        if i % 2 == 0:
            line = ax.axvline(lines[l], 0.8, 0.9, color="k")
            x, y = line.get_data()
            coords = ax.transLimits.transform((x[0], y[1]))
            ax.text(coords[0], 0.95, l, transform=ax.transAxes, ha="center", va="center", rotation=rotation)
        else:
            line = ax.axvline(lines[l], 0.1, 0.2, color="k")
            x, y = line.get_data()
            coords = ax.transLimits.transform((x[0], y[0]))
            ax.text(coords[0], 0.05, l, transform=ax.transAxes, ha="center", va="center", rotation=rotation)

    return ax


def get_tde_spectrum() -> Tuple[np.array, float, str]:
    """
    Return an array containing the TDE spectrum as well as returning the distance
    of the object and a string containing the reference for the observation.

    Parameters
    ----------
    None

    Returns
    -------
    tde_spec            np.array
                        The spectrum for the TDE
                            Column 0 is wavelength,
                            Column 1 is flux,
                            Column 2 is error (sometimes)
    observse_dist       float
                        The distance of the TDE given in cm
    reference           str
                        The journal reference for the observation
    """

    if TDE_OBJ is None:
        return
    elif TDE_OBJ == "iPTF15af" or TDE_OBJ == "iptf15af":
        z = 0.07897
        reference = "Blagorodnova et al. (2019)"
        observe_dist = 350 * 1e6 * PARSEC
        tde_spec = py_plot_util.get_iPTF15af_spec(SMOOTH, VERBOSE)
    elif TDE_OBJ == "ASSASN14li" or TDE_OBJ == "assasn14li":
        z = 0.02058
        reference = "Cenko et al. (2016)"
        observe_dist = 90 * 1e6 * PARSEC
        tde_spec = py_plot_util.get_ASSASN14li_spec(SMOOTH, VERBOSE)
    else:
        print("tde_spec_plot.get_tde_spectrum: can't find spectrum for object {}".format(TDE_OBJ))
        exit(1)

    tde_spec[:, 0] /= (z + 1)

    return tde_spec, observe_dist, reference


def spec_plot_one(root: str, inc: str) -> None:
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

    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)

    try:
        spectrum = py_plot_util.read_spec_file(spec_file, pandas_table=True)
    except IOError:
        print("tde_spec_spec.spec_plot_one: could not open file {}".format(spec_file))
        exit(1)

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()
    if PLOT_LINE_IDS:
        lines = py_plot_util.get_common_line_ids()

    wavelength = spectrum["Lambda"].values.astype(float)
    try:
        raw_flux = spectrum[inc].values.astype(float)
    except KeyError:
        print("tde_spec_plot.spec_plot_one: could not find inclination {}".format(inc))
        exit(1)
    flux = py_plot_util.smooth_flux(raw_flux, SMOOTH, VERBOSE)
    flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
    fmax, fmin = py_plot_util.get_ylims(wavelength, flux, WMIN, WMAX)

    # plot the spectrum
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    if TDE_OBJ:
        ax.semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ + " " + reference)
    ax.semilogy(wavelength, flux, label="i = {}".format(inc) + r"$^{\circ}$")
    ax.set_xlim(WMIN, WMAX)
    ax.set_ylim(fmin, fmax)
    ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    ax.set_xlabel(r"Wavelength ($\AA$)")
    ax.legend()
    if PLOT_LINE_IDS:
        plot_line_ids(ax, lines)

    # finish up the plot
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

    idx = root.find(".pf")
    if idx != -1:
        root = root[:idx]
    spec_file = "{}.spec".format(root)

    try:
        spectrum = py_plot_util.read_spec_file(spec_file, pandas_table=True)
    except IOError:
        print("tde_spec_plot.spec_plot_multiple: could not open file {}".format(spec_file))
        exit(1)

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, reference = get_tde_spectrum()
    if PLOT_LINE_IDS:
        lines = py_plot_util.get_common_line_ids()

    wavelength = spectrum["Lambda"].values.astype(float)
    inclinations = spectrum_inclinations(spectrum)
    n_spec = len(inclinations)
    shape = subplot_dims(n_spec)
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(12, 8), squeeze=False, sharex="col")

    index = 0
    for i in range(shape[0]):
        for j in range(shape[1]):
            # ensure that the x label is shown only on the bottom plots
            if i == shape[0] - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")

            if index > n_spec - 1:
                break

            if TDE_OBJ:
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)  #  + " " + reference)

            # extract the flux from the spectrum array and scale for distance
            inc = inclinations[index]
            raw_flux = spectrum[inc].values.astype(float)
            flux = py_plot_util.smooth_flux(raw_flux, SMOOTH, VERBOSE)
            flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
            fmax, fmin = py_plot_util.get_ylims(wavelength, flux, WMIN, WMAX)

            # now plot the spectrum :-)
            ax[i, j].semilogy(wavelength, flux, label=r"$i$ = {}".format(inc) + r"$^{\circ}$")
            ax[i, j].set_xlim(WMIN, WMAX)
            ax[i, j].set_ylim(fmin, fmax)
            ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
            ax[i, j].legend(loc="upper right")

            if PLOT_LINE_IDS:
                ax[i, j] = plot_line_ids(ax[i, j], lines)

            index += 1


    fig.suptitle(root)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_spectra.{}".format(root, FTYPE))
    plt.close()

    return


def spec_plot_comparison(name):
    """
    Create a comparison plot similar to spec_plot for multiple spectra which are
    search for recursively in the working directory.

    Parameters
    ----------
    name                str
                        The base name of the output plot

    Returns
    -------
    None
    """

    spec_files = py_plot_util.find_spec_files()
    if len(spec_files) == 0:
        print("No spec files found")
        print("\n--------------------------")
        exit(1)

    print("Creating spectra with the following simulations:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    observe_dist = 1
    if TDE_OBJ:
        tde, observe_dist, ref = get_tde_spectrum()
    if PLOT_LINE_IDS:
        lines = py_plot_util.get_common_line_ids()

    # Use the first spectrum as a template for what all the other spectra are like
    test_spec = py_plot_util.read_spec_file(spec_files[0], pandas_table=True)
    wavelength = test_spec["Lambda"].values.astype(float)
    inclination = spectrum_inclinations(test_spec)
    n_specs = len(inclination)

    shape = subplot_dims(n_specs)
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(18, 12), squeeze=False, sharex="col")

    index = 0
    for i in range(shape[0]):
        for j in range(shape[1]):
            if index > n_specs - 1:
                break

            if TDE_OBJ:
                ax[i, j].semilogy(tde[:, 0], tde[:, 1], label=TDE_OBJ)

            inc = inclination[index]

            # Loop over each spectrum from each simulation
            for file in spec_files:
                root, dir = py_plot_util.get_root_name_and_path(file)

                # Get the spectrum
                spectrum = py_plot_util.read_spec_file(file, pandas_table=True)
                raw_flux = spectrum[inc].values.astype(float)
                flux = py_plot_util.smooth_flux(raw_flux, SMOOTH, VERBOSE)
                flux *= DEFAULT_DIST ** 2 / observe_dist ** 2
                fmax, fmin = py_plot_util.get_ylims(wavelength, flux, WMIN, WMAX)

                # Plot the spectrum
                ax[i, j].semilogy(wavelength, flux, label=r"$i$ = {}".format(inc) + r"$^{\circ}$ "+dir)
                ax[i, j].set_xlim(WMIN, WMAX)
                ax[i, j].set_ylim(fmin, fmax)
                ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
                ax[i, j].legend(loc="upper right")
                if PLOT_LINE_IDS:
                    plot_line_ids(ax[i, j], lines)

            if i == shape[0] - 1:
                ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            index += 1

    fig.suptitle("comparison plot")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}_comparison.{}".format(name, FTYPE))
    plt.close()

    return


def main() -> None:
    """
    Main controlling function of the script.

    Parameters
    ----------
    None

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

    if args.inclination:
        print("Plotting {} spectrum for inclination {}".format(root, args.inclination))
        spec_plot_one(root, args.inclination)
    elif args.comparison:
        print("Plotting comparison spectrum")
        spec_plot_comparison(root)
    else:
        print("Plotting {} spectra for all inclination angles".format(root))
        spec_plot_multiple(root)

    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
