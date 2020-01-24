#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to act as a entry way into plotting a bunch of
the output from Python, as well as being able to make a bunch of default plots
which are helpful in the analysis of a Python simulation.

This script mostly acts as a wrapped around PyPython.SpectrumPlot and
PyPython.WindPlot, and as such all of the individual plotting routines can be
found in this module.
"""

import os
import argparse
import numpy as np
from consts import *
from matplotlib import pyplot as plt
from typing import List, Tuple, Union

from PyPython import SpectrumUtils
from PyPython import SpectrumPlot
from PyPython import Utils
from PyPython import WindUtils
from PyPython import WindPlot
from PyPython import Quotes

PLOTS = "all"
SHOW_PLOT = False
USE_POLAR_PROJECTION = False
VERBOSITY = False
PLOT_SCALE = False
DIMS = "2d"
WMIN = None
WMAX = None
OUTPUT_FILE_TYPE = "png"
SMOOTH = 5
DEFAULT_DIST = 100 * PARSEC
OBSERVE_DIST = 100 * PARSEC
MIN_FLUX = 1e-20
PLOT_LINE_OF_SIGHTS = False
USE_ROOT = False


def get_script_arguments() -> str:
    """
    Parse the various global parameters from the command line.

    Returns
    -------
    args.output_name: str
       The output base name for plots
    """

    global PLOTS
    global USE_POLAR_PROJECTION
    global VERBOSITY
    global SHOW_PLOT
    global PLOT_SCALE
    global WMIN
    global WMAX
    global OUTPUT_FILE_TYPE
    global SMOOTH
    global OBSERVE_DIST
    global DIMS
    global USE_ROOT
    global PLOT_LINE_OF_SIGHTS

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("plots", nargs="?", type=str, help="The type of plot to create")
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-otype", type=str, action="store", help="The file format of the output")
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
    if args.show:
        SHOW_PLOT = True
    if args.log:
        PLOT_SCALE = True
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax
    if args.otype:
        OUTPUT_FILE_TYPE = args.filetype
    if args.smooth:
        SMOOTH = int(args.smooth)
    if args.dist:
        OBSERVE_DIST = args.dist
    if args.polar:
        USE_POLAR_PROJECTION = True
    if args.dim_1d:
        DIMS = "1d"
    if args.root:
        USE_ROOT = True
    if args.los:
        PLOT_LINE_OF_SIGHTS = True

    return args.output_name


def all_inclinations(spec_names: List[str], delim: str = None) -> np.array:
    """
    Get all of the unique inclination angles for a set of Python .spec files.

    Parameters
    ----------
    spec_name: List[str]
        The directory path to Python .spec files
    delim: str [optional]
        The delimiter in the .spec files, assumed to be spaces by default

    Returns
    -------
    inclinations: List[int]
        All of the unique inclination angles found in the Python .spec files
    """

    n = all_inclinations.__name__

    # This is a bit of a unnecessary safety precuation :-)
    for i in range(len(spec_names)):
        if spec_names[i].find(".pf") != -1:
            spec_names[i] = spec_names[i].replace(".pf", ".spec")

    inclinations = []
    for i in range(len(spec_names)):
        spec = SpectrumUtils.read_spec(spec_names[i], delim)
        columns = spec.columns.values
        for j in range(len(columns)):
            if columns[j].isdigit() and int(columns[j]) not in inclinations:
                inclinations.append(int(columns[j]))

    inclinations.sort()

    return inclinations


def wind_plot(root: str, output_name: str, wind_names: List[str], wind_types: List[str], wd: str = "./",
              input_file: str = None, projection: str = "rectilinear", scale: str = "loglog",
              line_of_sights: bool = None, plot_indices: bool = False, ndims: str = "2d",
              subplot_dims: Tuple[int, int] = None, fig_size: Tuple[int, int] = (10, 15), plot_title: str = None,
              filetype: str = "png", verbose: bool = False) -> None:
    """
    Creates a figure containing 2D plots of each provided wind variable. For
    each variable in wind_names, there should be the corresponding variable
    type which can either be "wind" or "ion". The subplot dimensions should
    multiply together to be larger than the number of variables provided to
    the function.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    output_name: str
        The base output file name for the figure
    wind_names: List[str]
        A list containing the name of the wind variables to show in the figure.
        This list must be the same length as the wvar_types list
    wind_types: List[str]
        A list containing the types of the wind variables, this is either going
        to be "wind" or "ion" depending on if it is a wind/plasma variable or
        an ion. This list must be the same length as the wvars list
    wd: str [optional]
        An absolute or relative path to the directory containing the simulation
    ndims: str [optional]
        The dimensionality of the simulation, currently only allowed is "2d"
    projection: str [optional]
        The coordinate system of the simulation grid, allowed values are
        "rectilinear" and "polar"
    subplot_dims: Tuple[int, int] [optional]
        The number of rows and columns of the subplot panels
    fig_size: Tuple[int, int] [optional]
        The size of the figure in inches given by the tuple (width, height)
    plot_title: str [optional]
        If this is provided, the figure will include this title at the top
    loglog_scale: bool [optional]
        When True, the subplot panels will use a log-log scale. This is the
        default behaviour.
    filetype: str [optional]
        The file type of the output figure, set to png by default
    show_plot: bool [optional]
        Show the figure after saving the figure to disk
    verbose: bool [optional]
        Enable verbose logging
    input_file: str
        If a specific input file is required to use, then use this file path
        instead.
    plot_indices: bool [optional]
        Instead of using x, z coorindates plot the spatial components in terms
        of the cell's i j indices
    """

    n = wind_plot.__name__
    allowed_projections = ["rectilinear", "polar"]

    if ndims.lower() != "2d":
        raise NotImplementedError("{}: only understand ndims 2d at the moment".format(n))

    # Check to make sure each wind variable has a type corresponding type
    if len(wind_names) != len(wind_types):
        print("{}: wind_names and wind_types should be the same length".format(n))
        return

    # Check to ensure the subplot dimensions are provided and valid
    if subplot_dims is None:
        subplot_dims = (4, 2)
    if subplot_dims[0] * subplot_dims[1] < len(wind_names):
        print("{}: not enough subplot panels to plot all the provided wind variables".format(n))
        return

    # Check to see if the projection is understood and set up the figure and axes objects
    if projection not in allowed_projections:
        print("{}: projection {} not allowed, allowed values: rectilinear or polar".format(n, projection))
        return
    if projection == "rectilinear":
        fig, ax = plt.subplots(subplot_dims[0], subplot_dims[1], figsize=fig_size, squeeze=False)
    elif projection == "polar":
        fig = plt.figure(figsize=fig_size)
    else:
        return

    if plot_indices:
        scale = "linlin"

    index = 0
    for i in range(subplot_dims[0]):
        for j in range(subplot_dims[1]):
            if index > len(wind_names) - 1:
                break

            wind_name = wind_names[index]
            wind_type = wind_types[index]
            if verbose:
                print("\tPlotting {} of type {}".format(wind_name, wind_type))

            try:
                x, z, w = WindUtils.extract_wind_var(root, wind_name, wind_type, wd, projection,
                                                     input_file=input_file, return_indices=plot_indices)
            except Exception:  # Can't figure out what to catch here...?
                print("Exception occured >:(. Can't plot {} for some reason".format(wind_name))
                index += 1
                continue

            if projection == "rectilinear":
                fig, ax = WindPlot.create_rectilinear_wind_plot(x, z, w, wind_type, wind_name, fig, ax, i, j,
                                                                scale=scale, obs_los=line_of_sights)
            elif projection == "polar":
                ax = plt.subplot(subplot_dims[0], subplot_dims[1], index + 1, projection="polar")
                WindPlot.create_polar_wind_plot(x, z, w, wind_type, wind_name, ax, index + 1, line_of_sights, scale,
                                                fig_size)
            index += 1

    if plot_title:
        fig.suptitle(plot_title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}/{}.{}".format(wd, output_name, filetype))

    if SHOW_PLOT:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_tau_spec(root: str, wd: str = "./", plot_freq: bool = True, plot_edges: bool = True,
                  wmin: Union[float, bool] = None, wmax: Union[float, bool] = None, semilogy: bool = False,
                  loglog: bool = True) -> None:
    """
    Create an optical depth spectrum for a Python simulation. Requires a
    root.tau_spec.diag file or something. The figure can be created as a
    function of wavelength or frequency.

    Parameters
    ----------
    root: str
        The root name of the optical depth spectrum to plot
    wd: str
        The directory containing the simulation
    plot_freq: bool [optional]
        If True, the optical depth will be printed in frequency space
    plot_edges: bool  [optional]
        Label various absorption edges
    wmin: float [optional]
        The minimum wavelength to plot
    wmax: float [optional]
        The maximum wavelength to plot
    semilogy: bool [optional]
        If True, then the y axis will be log scaled
    loglog: bool [optional]
        If True, then the x and y axis will be log scaled
    """

    fig, ax, = SpectrumPlot.spec_plot_tau_spec(root, wd, wmin=wmin, wmax=wmax, logy=semilogy, loglog=loglog,
                                               show_absorption_edge_labels=plot_edges, frequency_space=plot_freq)
    plt.savefig("{}_tau_spec.png".format(root))

    if SHOW_PLOT:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_spec_comps(root: str, wd:str = "./", logy: bool = True, frequency_space: bool = False) -> None:
    """
    Create a figure of the different spectrum components which contribute to
    the emitted spectrum in Python. Note that all the spectrum components added
    together DO NOT equal the emitted spectrum.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    wd: str
        The absolute or relative path containing the Python simulation
    wmin: float [optional]
        The lower wavelength boundary for the figure
    wmax: float [optional]
        The upper wavelength boundary for the figure
    smooth: int [optional]
        The size of the boxcar filter to smooth the spectrum components
    logy: bool [optional]
        Use a log scale for the y axis of the figure
    frequency_space: bool [optional]
        Create the figure in frequency space instead of wavelength space
    verbose: bool [optional]
        Enable verbose output to screen
    """

    fig, ax = SpectrumPlot.spec_plot_components(root, wd, WMIN, WMAX, SMOOTH, logy, frequency_space, verbose=VERBOSITY)
    plt.savefig("{}_spec_comps.{}".format(root, OUTPUT_FILE_TYPE))

    if SHOW_PLOT:
        plt.show(block=False)
    else:
        plt.close()

    return


def plot_spectra(input_files: List[str], figure_inclinations: Union[List, np.array], output_name: str) -> None:
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
    wmin: float [optional]
        The smallest wavelength to plot
    wmax: float [optional]
        The largest wavelength to plot
    smooth: int [optional]
        The size of the window for the boxcar smoother
    filetype: str [optional]
        The file type of the figure saved to disk; default is png
    show_plot: bool [optional]
        Show the figure after saving the figure to disk
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
            flux = SpectrumUtils.smooth_spectrum(spec[inclination].values.astype(float), SMOOTH)
            flux *= (DEFAULT_DIST ** 2 / OBSERVE_DIST ** 2)  # TODO: should this be a function input?
            ax.semilogy(wavelength, flux, label=legend)

            # This is basically only required when the wavelength range is restricted
            tymax, tymin = SpectrumUtils.ylims(wavelength, flux, WMIN, WMAX)
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
        ax.set_xlim(WMIN, WMAX)
        ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=15)
        ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)
        ax.legend(loc="upper right")
        ax.set_title(r"{} $i$ = {}".format(root, inclination) + r"$^{\circ}$", fontsize=20)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig("{}_i{}.{}".format(output_name, inclination, OUTPUT_FILE_TYPE))

        if SHOW_PLOT:
            plt.show(block=False)
        else:
            plt.close()

    return


def main() -> None:
    """
    The main controlling function of the script. Has a complex system of routes
    which it can follow depending on the type of plot requested.
    """

    output_name = get_script_arguments()

    global PLOTS
    PLOTS = PLOTS.lower()
    allowed_plots = ["spec", "spec_comps", "wind", "ions", "tau_spec", "all", "help"]
    if PLOTS not in allowed_plots or PLOTS == "help":
        if PLOTS != "help":
            print("Don't know how to plot {}".format(PLOTS))
        print("Allowed plots are: spec, spec_comp, wind, ion, tau_spec or all")
        print("\n--------------------------")
        return

    if USE_POLAR_PROJECTION:
        projection = "polar"
    else:
        projection = "rectilinear"

    print("--------------------------\n")
    Quotes.random_quote()

    # If a specific root name has been provided, we will only plot that one
    # instead or something - I can't quite remember the functionality of this
    # but I am scared to remove it :-(
    if USE_ROOT:
        input_simulations = [output_name]
    else:
        if PLOTS == "spec" or PLOTS == "spec_comps":
            input_simulations = SpectrumUtils.find_specs()
        else:  # For when plotting the wind as well as the spectrum
            input_simulations = Utils.find_parameter_files()
        if len(input_simulations) == 0:
            print("No input files found")
            print("\n--------------------------")
            return
    print("Creating {} plots for the following simulations:\n".format(PLOTS))
    for i in range(len(input_simulations)):
        print("\t- {}".format(input_simulations[i]))
    print("\n--------------------------")

    # CREATE SPECTRA FOR EACH INCLINATION ANGLE
    if PLOTS == "spec" or PLOTS == "all":
        print("\nPlotting spectra".format(input_simulations))
        inclinations_in_spectra = all_inclinations(input_simulations)
        plot_spectra(input_simulations, inclinations_in_spectra, output_name)

    # The following plots can only be done when the script is called from a
    # directory containing a single simulation
    #    - spectrum components
    #    - continuum optical depth spectrum
    #    - wind variables, i.e. electron temperature
    #    - ion fractions or densities
    if len(input_simulations) == 1:
        root, path = Utils.split_root_directory(input_simulations[0])

        # Again, I can't remember the functionality of this... :-(
        if USE_ROOT:
            root = output_name

        # SPECTRUM COMPONENTS
        if PLOTS == "spec_comps" or PLOTS == "all":
            print("\nPlotting spectrum components")
            plot_spec_comps(root, path)

        # OPTICAL DEPTH SPECTRUM
        if PLOTS == "tau_spec" or PLOTS == "all":
            print("\nPlotting optical depth spectrum")
            plot_tau_spec(root, path)

        # RUN WINDSAVE2TABLE IF ROOT.EP.COMPLETE FILE IS MISSING
        if PLOTS == "wind" or PLOTS == "ions" or PLOTS == "all":
            if not os.path.isfile("{}.ep.complete".format(root)):
                Utils.windsave2table(root, path, VERBOSITY)

        # For plotting the line of sights on the wind plots
        if PLOT_LINE_OF_SIGHTS:
            spec_path = "{}/{}.spec".format(path, root)
            inclination_angles = SpectrumUtils.spec_inclinations(spec_path)
        else:
            inclination_angles = None

        # WIND VARIABLE
        if PLOTS == "wind" or PLOTS == "all":
            print("\nPlotting wind quantities")
            wind_names = ["t_e", "t_r", "ne", "rho", "c4", "converge", "ip", "ntot"]
            wind_types = ["wind"] * len(wind_names)
            wind_plot(root, output_name + "_wind", wind_names, wind_types, path, projection=projection,
                      line_of_sights=inclination_angles, ndims=DIMS, verbose=True)

        # WIND IONS
        if PLOTS == "ions" or PLOTS == "all":
            print("\nPlotting wind ions")
            dims = [(4, 2), (1, 2), (2, 2), (3, 2), (4, 2), (4, 2), (5, 3)]
            size = [(15, 20), (15, 5), (15, 10), (15, 15), (15, 20), (15, 20), (22.5, 25)]
            elements = ["ImportantBoys", "H", "He", "C", "N", "O", "Si"]
            ions = [
                ["O_i05", "Si_i04", "Si_i05", "N_i04", "N_i05", "N_i06", "C_i04", "C_i05"],
                ["H_i01", "H_i02"],
                ["He_i01", "He_i02", "He_i03"],
                ["C_i01", "C_i02", "C_i03", "C_i04", "C_i05", "C_i06"],
                ["N_i01", "N_i02", "N_i03", "N_i04", "N_i05", "N_i06", "N_i07", "N_i08"],
                ["O_i01", "O_i02", "O_i03", "O_i04", "O_i05", "O_i06", "O_i07", "O_i08"],
                ["Si_i01", "Si_i02", "Si_i03", "Si_i04", "Si_i05", "Si_i06", "Si_i07", "Si_i08", "Si_i09", "Si_i10",
                 "Si_i11", "Si_i12", "Si_i13", "Si_i14", "Si_i15"]
            ]
            for i in range(len(elements)):
                print("\tCreating ion plot for {}".format(elements[i]))
                ion_names = ions[i]
                ion_types = ["ion"] * len(ion_names)
                extra_name = "_" + elements[i] + "_ions"
                wind_plot(root, output_name + extra_name, ion_names, ion_types, path, projection=projection,
                          line_of_sights=inclination_angles, ndims=DIMS, subplot_dims=dims[i], fig_size=size[i])

        # REMOVE DATA SYMBOLIC LINK TO KEEP THINGS CLEAN FOR DROPBOX
        # Utils.remove_data_sym_links(path, VERBOSITY)

    elif len(input_simulations) > 1 and PLOTS != "spec":
        print("Can only plot {} when one root in folder :^)".format(PLOTS))

    print("")

    if SHOW_PLOT:
        input("Press enter to exit")

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
