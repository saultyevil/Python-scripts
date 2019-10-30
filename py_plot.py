#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to plot the output from a Python simulation.

The script can be called within a directory containing a single Python simulation,
or in a directory containing multiple Python simulations. In the case of being
called in a directory with multiple Python simulations, only spectra will be
plotted. However, if called in a directory containing a single simulation, then
the spectra, integration flux components, various wind quantities and wind ions
will be plotted.

You do not have to include the Python root name as an argument for this script,
as it uses the Unix find command to search for .spec files recursively to find
Python simulations.

usage: py_plot.py [-h] [-wmin WMIN] [-wmax WMAX] [-filetype FILETYPE]
                  [-smooth SMOOTH] [-dist DIST] [-p] [--dim_1d] [-v] [-s]
                  [-log] [-r]
                  output_name [plots]

positional arguments:
  output_name         The base name for the output
  plots               The type of plot to create

optional arguments:
  -h, --help          show this help message and exit
  -wmin WMIN          The smallest wavelength to show
  -wmax WMAX          The largest wavelength to show
  -filetype FILETYPE  The file format of the output
  -smooth SMOOTH      The amount of smoothing of the spectra
  -dist DIST          Distance of the observer
  -p, --polar         Project the wind on polar axes
  --dim_1d            Plot 1D data
  -v, --verbose       Increase output to screen
  -s, --show          Show the plots on screen
  -log                Enable log log axes
  -r, --root          Enables providing a root name instead
"""

import os
import pandas as pd
import argparse
import numpy as np
from consts import *
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
from PyPython import SpectrumUtils, Utils, WindUtils

PLOTS = "all"
POLAR_PROJECTION = False
VERBOSITY = False
SHOW_PLOT = False
PLOT_LOG = False
DIMS = "2d"
WMIN = None
WMAX = None
FILETYPE = "png"
SMOOTH = 5
DEFAULT_DIST = 100 * PARSEC
OBSERVE_DIST = 100 * PARSEC
MIN_FLUX = 1e-20
PLOT_OBS_LOS = False
ROOT = False


def get_script_arguments() -> str:
    """
    Parse the various global parameters from the command line.

    Returns
    -------
    args.output_name: str
       The output base name for plots
    """

    global PLOTS
    global POLAR_PROJECTION
    global VERBOSITY
    global SHOW_PLOT
    global PLOT_LOG
    global WMIN
    global WMAX
    global FILETYPE
    global SMOOTH
    global OBSERVE_DIST
    global DIMS
    global ROOT
    global PLOT_OBS_LOS

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("plots", nargs="?", type=str, help="The type of plot to create")
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help="The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help="The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help="Distance of the observer")
    p.add_argument("-los", action="store_true", help="Plot the observer line of sights in the simulation")
    p.add_argument("-p", "--polar", action="store_true", help="Project the wind on polar axes")
    p.add_argument("--dim_1d", action="store_true", help="Plot 1D data")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output to screen")
    p.add_argument("-s", "--show", action="store_true", help="Show the plots on screen")
    p.add_argument("-log", action="store_true", help="Enable log log axes")
    p.add_argument("-r", "--root", action="store_true", help="Enables providing a root name instead")
    args = p.parse_args()

    # Assign the optional arguments to their global vars
    if args.plots:
        PLOTS = args.plots
    if args.verbose:
        VERBOSE = True
    if args.show:
        SHOW_PLOT = True
    if args.log:
        PLOT_LOG = True
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax
    if args.filetype:
        FILETYPE = args.filetype
    if args.smooth:
        SMOOTH = int(args.smooth)
    if args.dist:
        OBSERVE_DIST = args.dist
    if args.polar:
        POLAR_PROJECTION = True
    if args.dim_1d:
        DIMS = "1d"
    if args.root:
        ROOT = True
    if args.los:
        PLOT_OBS_LOS = True

    return args.output_name


def all_inclinations(spec_names: List[str], delim: str = None) -> np.array:
    """
    Get all of the unique inclination angles for a set of Python .spec files.
    Parameters
    ----------
    spec_name: List[str]
        The directory path to Python .spec files
    delim: str, optional
        The delimiter in the .spec files, assumed to be spaces by default
    Returns
    -------
    inclinations: List[int]
        All of the unique inclination angles found in the Python .spec files
    """

    inclinations = []

    for i in range(len(spec_names)):
        if spec_names[i].find(".pf") != -1:
            spec_names[i] = spec_names[i].replace(".pf", ".spec")

    # Find the viewing angles in each .spec file
    for i in range(len(spec_names)):
        spec = SpectrumUtils.read_spec(spec_names[i], delim)

        if type(spec) == pd.core.frame.DataFrame:
            col_names = spec.columns.values
        elif type(spec) == np.ndarray:
            col_names = spec[0, :]
        else:
            raise TypeError("{}:{}: unknown data type {} for function"
                            .format(__file__, all_inclinations.__name__, type(spec)))

        # Go over the columns and look for viewing angles
        for j in range(len(col_names)):
            if col_names[j].isdigit() is True:
                angle = int(col_names[j])
                duplicate_flag = False
                for va in inclinations:  # Check for duplicate angle
                    if angle == va:
                        duplicate_flag = True
                if duplicate_flag is False:
                    inclinations.append(angle)

    return inclinations


def sightline_coords(x: np.ndarray, theta: float):
    """
    Return the vertical coordinates for a sightline given the x coordinates
    and the inclination of the sightline.

    Parameters
    ----------
    x: np.ndarray[float]
        The x-coordinates of the sightline
    theta: float
        The opening angle of the sightline

    Returns
    -------
    z: np.ndarray[float]
        The z-coordinates of the sightline
    """

    return x * np.tan(np.pi / 2 - theta)


def rectilinear_wind_plot(fig: plt.Figure, ax: plt.Axes, x: np.ndarray, z: np.ndarray, w: np.ndarray, i: int, j: int,
                          wvar: str, wvar_t: str, loglog_scale: bool = True, inclination_angles: List[float] = None) \
                          -> Tuple[plt.Figure, plt.axes]:
    """
    Create a wind plot in rectilinear coordinates.

    Parameters
    ----------
    fig: plt.Figure
        A matplotlib.pyplot figure object which contains the axes objects and
        etc to create the desired figure
    ax: plt.Axes
        A matplotlib.pyplot plt.Axes object which contains the subplot
        panels of the desired figure
    x: np.array[float]
        The x coordinate points for the wind variable
    z: np.array[float]
        The z coordinate points for the wind vairbale
    w: np.ndarray[float]
        The wind variable to be plotted - this should be a 2d masked array
        The quantity of interest to be plotted
    i: int
        The i-th (row) index for subplot panels
    j: int
        The j-th (column) index for the subplot panels
    wvar: List[str]
        The name of the wind variable being plotted in the panel
    wvar_t: List[str]
        The type of the wind variable being plotted in the panel
    loglog_scale: bool, optional
        If this is True then the figure will be created on a loglog scale
    inclination_angles: List[float], optional
        If this is provided, line of sights for each inclination angle
        will be overplotted the wind

    Returns
    -------
    fig: plt.Figure
        The plt.Figure object creating this figure.
    ax: plt.Axes
        The plt.Axes object which contains the panel of the figure.
    """

    with np.errstate(divide="ignore"):
        if wvar == "converge" or wvar == "converging":
            im = ax[i, j].pcolor(x, z, w)
        elif wvar_t == "ion":
            im = ax[i, j].pcolor(x, z, np.log10(w), vmin=-5, vmax=0)
        elif wvar_t == "wind":
            im = ax[i, j].pcolor(x, z, np.log10(w))
        else:
            print("{}:{}: type {} not recognised".format(__file__, rectilinear_wind_plot.__name__, wvar_t))
            return fig, ax

    if inclination_angles:
        xsight = np.linspace(0, np.max(x), int(1e5))
        for inc in inclination_angles:
            zsight = sightline_coords(xsight, np.deg2rad(inc))
            ax[i, j].plot(xsight, zsight, label="i = {}".format(inc) + r"$^{\circ}$ sightline")
        ax[i, j].legend()

    fig.colorbar(im, ax=ax[i, j])
    ax[i, j].set_xlabel("x")
    ax[i, j].set_ylabel("z")
    if wvar == "converge" or wvar == "converging":
        ax[i, j].set_title(wvar)
    else:
        ax[i, j].set_title(r"$\log_{10}$(" + wvar + ")")

    if loglog_scale:
        ax[i, j].set_xscale("log")
        ax[i, j].set_yscale("log")

    ax[i, j].set_xlim(x[1, 1], x[-1, -1])
    ax[i, j].set_ylim(z[1, 1], z[-1, -1])

    return fig, ax


def polar_wind_plot(r: np.ndarray, theta: np.ndarray, w: np.ndarray, index: int, wvar: str, wvar_t: str,
                    subplot_dims: Tuple[int, int], inclination_angles: List[float] = None) -> None:
    """
    Create a subplot panel for a polar plot. By itself, this can actually make
    a single panel. However, the intended use of this function is to be
    used from plot_wind.

    Parameters
    ----------
    r: np.array[float]
        The r coordinate points for the wind variable
    theta: np.array[float]
        The theta coordinate points for the wind variable
    w: np.ndarray[float]
        The wind variable to be plotted - this should be a 2d masked array
    index: int
        The index of the subplot panel in the overall figure,
        (nrows, ncols, index) - note that this uses 0th indexing
    wvar: List[str]
        The name of the wind variable being plotted in the panel
    wvar_t: List[str]
        The type of the wind variable being plotted in the panel
    subplot_dims: Tuple[int, int]
        The number of rows and columns of subplot panels, (nrows, ncols)
    inclination_angles: List[float], optional
        If this is provided, line of sights for each inclination angle
        will be overplotted the wind
    """

    ax = plt.subplot(subplot_dims[0], subplot_dims[1], index + 1, projection="polar")

    with np.errstate(divide="ignore"):
        if wvar == "converge":
            im = ax.pcolor(theta, np.log10(r), w)
        elif wvar_t == "wind":
            im = ax.pcolor(theta, np.log10(r), np.log10(w))
        elif wvar_t == "ion":
            im = ax.pcolor(theta, np.log10(r), np.log10(w), vmin=-5, vmax=0)

    if inclination_angles:
        xsight = np.linspace(0, np.max(r), int(1e5))
        for inc in inclination_angles:
            zsight = sightline_coords(xsight, np.deg2rad(90 - inc))
            rsight = np.sqrt(xsight ** 2 + zsight ** 2)
            thetasight = np.arctan2(zsight, xsight)
            ax.plot(thetasight, np.log10(rsight), label="i = {}".format(inc) + r"$^{\circ}$ sightline")
        # ax.legend()

    plt.colorbar(im, ax=ax)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_rlabel_position(90)
    ax.set_ylabel("Log[R]")
    rmin = r[1][0]
    rmax = r[-2][0]
    ax.set_rlim(np.log10(rmin), np.log10(rmax))

    if wvar == "converge":
        ax.set_title(r"convergence")
    else:
        ax.set_title(r"$\log_{10}$(" + wvar + ")")

    return


def plot_wind(root: str, output_name: str, wvars: List[str], wvar_types: List[str], path: str = "./",
              ndims: str = "2d", projection: str = "rectilinear", subplot_dims: Tuple[int, int] = None,
              fig_size: Tuple[int, int] = (10, 15), plot_title: str = None, loglog_scale: bool = True,
              filetype: str = "png", show_plot: bool = False, verbose: bool = False,
              input_file: str = None, inclination_angles: bool = None) -> None:
    """
    Creates a 2D figure with subplot panels for each provided wind variable
    given in the wvars list; for each wvar, there should be a provided variable
    type in the list wvar_types. The subplot panel dimensions in the variable
    subplot_dims should multiply together to be larger than the number of
    wind variables passed in wvars.

    TODO: polar plot doesn't take into account using log-log scale or not

    Notes
    -----
    Currently, this will only work with 2D data and only with rectilinear (x, y)
    and polar (r, theta) coordinate systems.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    output_name: str
        The base output file name for the figure
    wvars: List[str]
        A list containing the name of the wind variables to show in the figure.
        This list must be the same length as the wvar_types list
    wvar_types: List[str]
        A list containing the types of the wind variables, this is either going
        to be "wind" or "ion" depending on if it is a wind/plasma variable or
        an ion. This list must be the same length as the wvars list
    path: str, optional
        An absolute or relative path to the directory containing the simulation
    ndims: str, optional
        The dimensionality of the simulation, currently only allowed is "2d"
    projection: str, optional
        The coordinate system of the simulation grid, allowed values are
        "rectilinear" and "polar"
    subplot_dims: Tuple[int, int], optional
        The number of rows and columns of the subplot panels
    fig_size: Tuple[int, int], optional
        The size of the figure in inches given by the tuple (width, height)
    plot_title: str, optional
        If this is provided, the figure will include this title at the top
    loglog_scale: bool, optional
        When True, the subplot panels will use a log-log scale. This is the
        default behaviour.
    filetype: str, optional
        The file type of the output figure, set to png by default
    show_plot: bool, optional
        Show the figure after saving the figure to disk
    verbose: bool, optional
        Enable verbose logging
    input_file: str
        If a specific input file is required to use, then use this file path
        instead.
    """

    if ndims.lower() != "2d":
        raise NotImplementedError("{}:{}: only understand ndims 2d at the moment".format(__file__, plot_wind.__name__))

    if len(wvars) != len(wvar_types):
        print("{}:{}: vars and types should be of the same length".format(__file__, plot_wind.__name__,))
        return

    # If subplot dimensions are not provided, then assume a default shape of (4, 2)
    if subplot_dims is None:
        subplot_dims = (4, 2)
    if subplot_dims[0] * subplot_dims[1] < len(wvars):
        print("{}:{}: not enough subplot panels to plot all the provided vars".format(__file__, plot_wind.__name__,))
        return

    allowed_projections = ["rectilinear", "polar"]
    if projection not in allowed_projections:
        print("{}:{}: projection {} not allowed, allowed values: rectilinear or polar"
              .format(__file__, plot_wind.__name__, projection))
        return

    if projection == "rectilinear":
        fig, ax = plt.subplots(subplot_dims[0], subplot_dims[1], figsize=fig_size, squeeze=False)
    elif projection == "polar":
        fig = plt.figure(figsize=fig_size)

    index = 0
    for i in range(subplot_dims[0]):
        for j in range(subplot_dims[1]):
            if index > len(wvars) - 1:
                break
            wvar = wvars[index]
            wvar_t = wvar_types[index]
            if verbose:
                print("\tPlotting {} of type {}".format(wvar, wvar_t))
            try:
                x, z, w = WindUtils.extract_wind_var(root, wvar, wvar_t, path, projection, input_file=input_file)
            except Exception:
                print("Exception occured: Unable to plot {} for some god forsaken fucking cunting bastard reason".format(wvar))
                index += 1
                continue

            if projection == "rectilinear":
                fig, ax = rectilinear_wind_plot(fig, ax, x, z, w, i, j, wvar, wvar_t, loglog_scale, inclination_angles)
            elif projection == "polar":
                polar_wind_plot(x, z, w, index, wvar, wvar_t, subplot_dims, inclination_angles=inclination_angles)

            index += 1

    # Finishing touches of the plot
    if plot_title:
        fig.suptitle(plot_title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}/{}.{}".format(path, output_name, filetype))

    if show_plot:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_tau_spec(root: str, dir: str = "./", plot_freq: bool = False, plot_edges: bool = True,
                  wmin: Union[float, bool] = None, wmax: Union[float, bool] = None, semilogy: bool = True,
                  loglog: bool = False, show_plot: bool = False) -> None:
    """
    Create an optical depth spectrum for a Python simulation. Requires a
    root.tau_spec.diag file or something. The figure can be created as a
    function of both wavelength and frequency.

    Parameters
    ----------
    root: str
        The root name of the optical depth spectrum to plot
    dir: str
        The directory containing the simulation
    plot_freq: bool, optional
        If True, the optical depth will be printed in frequency space
    plot_edges: bool, optional
        Label various absorption edges
    wmin: float, optional
        The minimum wavelength to plot
    wmax: float, optional
        The maximum wavelength to plot
    semilogy: bool, optional
        If True, then the y axis will be log scaled
    loglog: bool, optional
        If True, then the x and y axis will be log scaled
    """

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    filename = "{}/{}.tau_spec.diag".format(dir, root, root, root)

    try:
        spec = np.loadtxt(filename, skiprows=1)
    except IOError:
        try:
            filename = "{}/diag_{}/".format(dir, root) + "{}.tau_spec.diag".format(root)
            spec = np.loadtxt(filename, skiprows=1)
        except IOError:
            print("{}:{}: unable to find the optical depth spectrum: {}".format(__file__, plot_tau_spec.__name__, root))
            return

    with open(filename, "r") as f:
        cols = f.readline().split()
    if cols[0] == "#":
        cols = cols[1:]
    nangles = len(cols) - 1

    if plot_freq:
        spec[:, 0] = C / (spec[:, 0] * ANGSTROM)

    if spec[0, 0] < spec[-1, 0]:
        xxlims = [spec[0, 0], spec[-1, 0]]
    else:
        xxlims = [spec[-1, 0], spec[0, 0]]

    if wmin:
        xxlims[0] = wmin
    if wmax:
        xxlims[1] = wmax

    for i in range(nangles):
        the_label = r"$i$ = " + str(float(cols[i + 1][1:3])) + r"$^{\circ}$"

        # skip if there is no optical depth for this angle
        n_non_zeros = np.count_nonzero(spec[:, i + 1])
        if n_non_zeros == 0:
            continue

        if loglog:
            ax.semilogy(np.log10(spec[:, 0]), spec[:, i + 1], label=the_label)
        elif semilogy:
            ax.semilogy(spec[:, 0], spec[:, i + 1], label=the_label)
        else:
            ax.plot(spec[:, 0], spec[:, i + 1], label=the_label)
        ax.tick_params(axis="both", which="major", labelsize=13)

    if plot_freq:
        if loglog:
            ax.set_xlabel(r"Log[Frequency], Hz")
        else:
            ax.set_xlabel(r"Frequency, Hz")
    else:
        ax.set_xlabel(r"Wavelength, $\AA$", fontsize=15)
    ax.set_ylabel(r"Optical Depth, $\tau$", fontsize=15)

    ylims = ax.get_ylim()
    ax.set_ylim(ylims[0] / 10, ylims[1] * 10)

    if plot_edges:
        SpectrumUtils.plot_line_ids(ax, SpectrumUtils.absorption_edges(freq=plot_freq, log=loglog))

    ax.legend()

    plt.savefig("{}_tau_spec.png".format(root))

    if show_plot:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_spec_comps(input_file: str, output_name: str, semilogy_scale: bool = False, smooth: int = 5,
                    wmin: Union[float, bool] = None, wmax: Union[float, bool] = None, filetype: str = "png",
                    show_plot: bool = False, verbose: bool = False) -> None:
    """
    Create a figure of the different spectrum components which contribute to
    the emitted spectrum in Python. Note that all the spectrum components added
    together DO NOT equal the emitted spectrum.

    From Knox:
    0 - The spectra with their original weights (before transmission through the
        wind)
    1 - The emitted spectra with their weights revised by transmission through
        the wind
    2 - The emitted spectra from photons that originated from the central object
    3 - The emitted spectra from photons that originated0.2 from the disk
    4 - The emitted spectra from photons that originated in the wind  (Note that
        with modification that we made recently to account for adiabatic
        heating, there are actually photons, generally quite few, that can be in
        this spectrum even in the macro-atom case).
    5 - The spectra of photons which hit (and were absorbed by the star) at the
        time they were absorbed.  This spectrum should be all zero’s for in the
        case where we scatter at a surface as this is not the flux that hit the
        surface, but the flux that is lost.
    6 - The flux of photons that scattered at least once in the wind. Note that
        this appears to be all photons that have scattered, regardless of
        whether they escape. The weights used are the final weights.

    Parameters
    ----------
    input_file: str
        A relative or absolute path to the .spec file
    output_name: str
        The base output file name for the figure
    semilogy_scale: bool, optional
        Use a log-linear scale for the plot
    smooth: int, optional
        The size of the window for the boxcar smoother
    wmin: float, optional
        The smallest wavelength to plot
    wmax: float, optional
        The largest wavelength to plot
    filetype: str, optional
        The file type of the figure saved to disk; default is png
    show_plot: bool, optional
        Show the figure after saving the figure to disk
    verbose: bool, optional
        Enable verbose logging
    """

    if type(input_file) is not str:
        print("{}:{}: can only plot spectrum components for one spectrum at a time"
              .format(__file__, plot_spec_comps.__name__))
        return

    if input_file.find(".spec") == -1:
        if input_file.find(".pf") != -1:
            input_file = input_file.replace(".pf", ".spec")
        else:
            print("{}:{}: can't fix file path {}".format(__file__, plot_spec_comps.__name__, input_file))

    fig, ax = plt.subplots(2, 1, figsize=(12, 10))

    # These are the headers which should be in Python spec files. They are the different components
    # of the Python spectra
    headers_top = ["Created", "Emitted"]
    headers_bot = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    try:
        spec = SpectrumUtils.read_spec(input_file)
    except Exception:
        print("{}:{}: could not open input file {}".format(__file__, plot_spec_comps.__name__, input_file))
        return

    wavelength = spec["Lambda"].values.astype(float)

    # First panel, plot the created and emitted emission
    for i in range(len(headers_top)):
        print("\tPlotting {}".format(headers_top[i]))
        flux = SpectrumUtils.smooth_spectrum(spec[headers_top[i]].values.astype(float), smooth)
        if len(flux[flux < MIN_FLUX]) > 0.7 * len(flux):
            print("\t\t!!Skipping {}".format(headers_top[i]))
            continue
        if semilogy_scale:
            ax[0].semilogy(wavelength, flux, label=headers_top[i])
        else:
            ax[0].plot(wavelength, flux, label=headers_top[i])

    xlims = [wavelength.min(), wavelength.max()]
    if wmin:
        xlims[0] = wmin
    if wmax:
        xlims[1] = wmax
    ax[0].set_xlim(xlims[0], xlims[1])
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[0].legend()

    # Second panel, plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bot)):
        print("\tPlotting {}".format(headers_bot[i]))
        flux = SpectrumUtils.smooth_spectrum(spec[headers_bot[i]].values.astype(float), smooth)
        if len(flux[flux < MIN_FLUX]) > 0.7 * len(flux):
            print("\t\t!!Skipping {}".format(headers_bot[i]))
            continue
        if semilogy_scale:
            ax[1].semilogy(wavelength, flux, label=headers_bot[i])
        else:
            ax[1].plot(wavelength, flux, label=headers_bot[i])

    ax[1].set_xlim(xlims[0], xlims[1])
    ax[1].set_ylim(MIN_FLUX)
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.{}".format(output_name, filetype))

    if show_plot:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_spectra(input_files: List[str], figure_inclinations: Union[List, np.array], output_name: str,
                 wmin: float = None, wmax: float = None, smooth: int = 5, filetype: str = "png",
                 show_plot: bool = False, verbose: bool = False) -> None:
    """
    Plot the spectrum for each provided inclination angle for a Python .spec
    file. The directory to a .pf or a .spec file can be passed, as the function
    should be able to cope with either...

    Parameters
    ----------
    input_files: str
        The relative or absolute path to the .spec or .pf file
    figure_inclinations: List[str]
        A list of the inclination angles to be plotted
    output_name: str
        The base output file name for the figure
    wmin: float, optional
        The smallest wavelength to plot
    wmax: float, optional
        The largest wavelength to plot
    smooth: int, optional
        The size of the window for the boxcar smoother
    filetype: str, optional
        The file type of the figure saved to disk; default is png
    show_plot: bool, optional
        Show the figure after saving the figure to disk
    verbose: bool, optional
        Enable verbose logging
    """

    nfiles = len(input_files)

    for i in range(nfiles):
        path = input_files[i]
        if path.find(".spec") != -1:
            continue
        elif path.find(".pf") != -1:
            path = path.replace(".pf", ".spec")
        input_files[i] = path

    # Loop over each possible viewing angle as provided in inclinations
    for inclination in figure_inclinations:
        ymin = +1e99
        ymax = -1e99
        inclination = str(inclination)
        fig, ax = plt.subplots(1, 1, figsize=(12, 5))

        for file in input_files:
            if file.find(".spec") == -1 and file.find(".pf") == -1:
                file += ".spec"
            root, filepath = Utils.split_root_directory(file)
            legend = filepath + root
            print("\tPlotting {} for i = {}°".format(legend, inclination))

            # Read in the spectrum and check that it can be plotted for the
            # current inclination
            spec = SpectrumUtils.read_spec(file)
            allowed_inclination = SpectrumUtils.check_inclination(inclination, spec)
            if not allowed_inclination:
                continue

            # Extract the wavelength range and flux and scale the flux appropriately
            wavelength = spec["Lambda"].values.astype(float)
            flux = SpectrumUtils.smooth_spectrum(spec[inclination].values.astype(float), smooth)
            flux *= (DEFAULT_DIST ** 2 / OBSERVE_DIST ** 2)  # TODO: should this be a function input?
            ax.semilogy(wavelength, flux, label=legend)

            # This is basically only required when the wavelength range is restricted
            tymax, tymin = SpectrumUtils.ylims(wavelength, flux, wmin, wmax)
            if tymin is not None and tymin < ymin:
                ymin = tymin
            if tymax is not None and tymax > ymax:
                ymax = tymax

        # If these values haven't been updated to something sensible, set to none
        # otherwise matplotlib complains
        if ymin == +1e99:
            ymin = None
        if ymax == -1e99:
            ymax = None

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(wmin, wmax)
        ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=15)
        ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)
        ax.legend(loc="upper right")
        ax.set_title(r"{} $i$ = {}".format(root, inclination) + r"$^{\circ}$", fontsize=20)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig("{}_i{}.{}".format(output_name, inclination, filetype))

        if show_plot:
            plt.show(block=False)
        else:
            plt.close()

    return


def main() -> None:
    """
    The main controlling function of the script. Has a complex system of routes
    which it can follow depending on the type of plot requested.
    """

    allowed_plots = ["spec", "spec_comps", "wind", "ions", "tau_spec", "all", "help"]

    # Parse the running options from the command line
    output_name = get_script_arguments()

    if POLAR_PROJECTION:
        projection = "polar"
    else:
        projection = "rectilinear"

    print("--------------------------\n")

    global PLOTS
    PLOTS = PLOTS.lower()

    if PLOTS not in allowed_plots or PLOTS == "help":
        if PLOTS != "help":
            print("Don't know how to plot {}".format(PLOTS))
        print("Allowed plots are: spec, spec_comp, wind, ion, tau_spec or all")
        print("\n--------------------------")
        return

    if ROOT:
        input_files = [output_name]
    else:
        if PLOTS == "spec" or PLOTS == "spec_comps":
            input_files = SpectrumUtils.find_specs()
        else:
            input_files = Utils.find_parameter_files()
        if len(input_files) == 0:
            print("No input files found")
            print("\n--------------------------")
            return

    print("Creating {} plots for the following simulations:\n".format(PLOTS))
    for i in range(len(input_files)):
        print("\t- {}".format(input_files[i]))

    print("\n--------------------------")

    # CREATE SPECTRA FOR EACH INCLINATION ANGLE
    if PLOTS == "spec" or PLOTS == "all":
        inclinations = all_inclinations(input_files)
        print("\nPlotting spectra".format(input_files))
        plot_spectra(input_files, inclinations, output_name, wmin=WMIN, wmax=WMAX, smooth=SMOOTH, filetype=FILETYPE,
                     show_plot=SHOW_PLOT, verbose=VERBOSITY)

    # IF THIS IS CALLED IN A DIRECTORY WITH A SINGLE SIMULATION
    if len(input_files) == 1:
        root, path = Utils.split_root_directory(input_files[0])
        if ROOT:
            root = output_name

        if PLOT_OBS_LOS:
            spec_path = "{}/{}.spec".format(path, root)
            inclination_angles = SpectrumUtils.spec_inclinations(spec_path)
        else:
            inclination_angles = None

        # CREATE SPECTRUM COMPONENTS FILE
        if PLOTS == "spec_comps" or PLOTS == "all":
            print("\nPlotting spectrum components")
            input_files[0] = input_files[0].replace(".pf", ".spec")
            plot_spec_comps(input_files[0], output_name, semilogy_scale=PLOT_LOG, smooth=SMOOTH, filetype=FILETYPE,
                            show_plot=SHOW_PLOT, verbose=VERBOSITY, wmin=WMIN, wmax=WMAX)

        # CREATE OPTICAL DEPTH SPECTRUM
        if PLOTS == "tau_spec" or PLOTS == "all":
            print("\nPlotting optical depth spectrum")
            plot_tau_spec(root, path, wmin=WMIN, wmax=WMAX, show_plot=SHOW_PLOT,)

        # RUN WINDSAVE2TABLE IF ROOT.EP.COMPLETE FILE IS MISSING
        if PLOTS == "wind" or PLOTS == "ions" or PLOTS == "all":
            if not os.path.isfile("{}.ep.complete".format(root)):
                try:
                    Utils.windsave2table(root, path, VERBOSITY)
                except:
                    return

        # CREATE WIND VARIABLE PLOT
        if PLOTS == "wind" or PLOTS == "all":
            print("\nPlotting wind quantities")
            vars = ["t_e", "t_r", "ne", "rho", "c4", "converge", "ip", "ntot"]
            var_types = ["wind"] * len(vars)
            plot_wind(root, output_name + "_wind", vars, var_types, path, projection=projection, filetype=FILETYPE,
                      ndims=DIMS, show_plot=SHOW_PLOT, verbose=True, inclination_angles=inclination_angles)

        # CREATE ION PLOTS
        if PLOTS == "ions" or PLOTS == "all":
            print("\nPlotting wind ions")
            inames = ["", "H", "He", "C", "N", "O", "Si"]
            dims = [(4, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 2), (5, 3)]
            size = [(15, 20), (15, 5), (15, 10), (15, 15), (15, 20), (15, 20), (22.5, 25)]
            # inames = ["", "H", "He", "C", "N", "O", "Al", "Si", "S"]
            # dims = [(4, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 2), (5, 3), (5, 3), (6, 3)]
            # size = [(15, 20), (15, 5), (15, 10), (15, 15), (15, 20), (15, 20), (22.5, 25), (22.5, 25), (22.5, 30)]
            ions = [
                ["O_i05", "Si_i04", "Si_i05", "N_i04", "N_i05", "N_i06", "C_i04", "C_i05"],
                ["H_i01", "H_i02"],
                ["He_i01", "He_i02", "He_i03"],
                ["C_i01", "C_i02", "C_i03", "C_i04", "C_i05", "C_i06"],
                ["N_i01", "N_i02", "N_i03", "N_i04", "N_i05", "N_i06", "N_i07", "N_i08"],
                ["O_i01", "O_i02", "O_i03", "O_i04", "O_i05", "O_i06", "O_i07", "O_i08"],
                # ["Al_i01", "Al_i02", "Al_i03", "Al_i04", "Al_i05", "Al_i06", "Al_i07", "Al_i08", "Al_i09", "Al_i10",
                #  "Al_i11", "Al_i12", "Al_i13", "Al_i14"],
                ["Si_i01", "Si_i02", "Si_i03", "Si_i04", "Si_i05", "Si_i06", "Si_i07", "Si_i08", "Si_i09", "Si_i10",
                 "Si_i11", "Si_i12", "Si_i13", "Si_i14", "Si_i15"],
                # ["S_i01", "S_i02", "S_i03", "S_i04", "S_i05", "S_i06", "S_i07", "S_i08", "S_i09", "S_i10", "S_i11",
                #  "S_i12", "S_i13", "S_i14", "S_i15", "S_i16", "S_i17"]
            ]
            for i in range(len(ions)):
                print("\tCreating ion plot {}".format(inames[i]))
                vars = ions[i]
                name = "_" + inames[i] + "_ions"
                var_types = ["ion"] * len(vars)
                plot_wind(root, output_name + name, vars, var_types, path, projection=projection, filetype=FILETYPE,
                          ndims=DIMS, verbose=VERBOSITY, show_plot=SHOW_PLOT, subplot_dims=dims[i], fig_size=size[i],
                          inclination_angles=inclination_angles)

        # REMOVE DATA SYMBOLIC LINK TO KEEP THINGS CLEAN FOR DROPBOX
        Utils.remove_data_sym_links(path, VERBOSITY)

    elif len(input_files) > 1 and PLOTS != "spec":
        print("Can only plot {} when one root in folder :^)".format(PLOTS))

    print("")  # spacer :-)

    if SHOW_PLOT:
        input("Press enter to exit")

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
