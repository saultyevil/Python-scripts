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

TODO: switch to using pfs instead of spec files
TODO: add check to see if wind files already exist to avoid issues with windsave2table
"""

import os
import argparse
import py_rm_data
import numpy as np
import py_plot_util
from consts import *
from matplotlib import pyplot as plt
from typing import List, Tuple, Union

PLOTS = "all"
POLAR_PROJECTION = False
VERBOSE = False
SHOW_PLOT = False
SPEC_LOGLOG = False
WMIN = None
WMAX = None
FILETYPE = "png"
SMOOTH = 15
DEFAULT_DIST = 100 * PARSEC
OBSERVE_DIST = 100 * PARSEC


def get_script_arguments() -> Tuple[str, str]:
    """
    Parse the various global parameters from the command line.

    Returns
    -------
    args.output_name       str
                           The output base name for plots
    """

    global PLOTS
    global POLAR_PROJECTION
    global VERBOSE
    global SHOW_PLOT
    global SPEC_LOGLOG
    global WMIN
    global WMAX
    global FILETYPE
    global SMOOTH
    global CHANGE_DIST
    global OBSERVE_DIST
    global TDE_PLOT

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("plots", nargs="?", type=str, help="The type of plot to create")
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help="The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help="The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help="Distance of the observer")
    p.add_argument("-p", "--polar", action="store_true", help="Project the wind on polar axes")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output to screen")
    p.add_argument("-s", "--show", action="store_true", help="Show the plots on screen")
    p.add_argument("-loglog", action="store_true", help="Enable log log axes")
    args = p.parse_args()

    # Assign the optional arguments to their global vars
    if args.plots:
        PLOTS = args.plots
    if args.verbose:
        VERBOSE = True
    if args.show:
        SHOW_PLOT = True
    if args.loglog:
        SPEC_LOGLOG = True
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

    return args.output_name


def plot_rectilinear_wind(fig, ax, x, z, qoi, i, j, var, var_t, loglog_scale) -> Tuple[plt.axes, plt.axes]:
    """

    Parameters
    ----------
    fig
    ax
    x
    z
    qoi
    i
    j
    var
    var_t
    loglog_scale

    Returns
    -------
    None
    """

    with np.errstate(divide="ignore"):
        if var == "converge":
            im = ax[i, j].pcolor(x, z, qoi)
        elif var_t == "ion":
            im = ax[i, j].pcolor(x, z, np.log10(qoi), vmin=-5, vmax=0)
        elif var_t == "wind":
            im = ax[i, j].pcolor(x, z, np.log10(qoi))
        else:
            print("py_plot.plot_python_wind: type {} not recognised".format(type))
            return fig, ax

    fig.colorbar(im, ax=ax[i, j])
    ax[i, j].set_xlabel("x")
    ax[i, j].set_ylabel("z")
    if var == "converge":
        ax[i, j].set_title(r"convergence")
    else:
        ax[i, j].set_title(r"$\log_{10}$(" + var + ")")

    if loglog_scale:
        ax[i, j].set_xscale("log")
        ax[i, j].set_yscale("log")
        ax[i, j].set_xlim(x[1, 1], x[-1, -1])
        ax[i, j].set_ylim(z[1, 1], z[-1, -1])

    return fig, ax


def plot_polar_wind(r, theta, qoi, index, var, var_t, subplot_dims, verbose):
    """

    Parameters
    ----------
    r
    theta
    qoi
    index
    var
    var_t
    subplot_dims
    verbose

    Returns
    -------
    None
    """

    theta = np.deg2rad(theta)
    ax = plt.subplot(subplot_dims[0], subplot_dims[1], index + 1, projection="polar")

    with np.errstate(divide="ignore"):
        if var == "converge":
            im = ax.pcolor(theta, np.log10(r), qoi)
        elif var_t == "wind":
            im = ax.pcolor(theta, np.log10(r), np.log10(qoi))
        elif var_t == "ion":
            im = ax.pcolor(theta, np.log10(r), np.log10(qoi), vmin=-5, vmax=0)

    plt.colorbar(im, ax=ax)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_rlim(np.log10(2.65e13), 17)  # TODO: remove hard coded limit

    if var == "converge":
        ax.set_title(r"convergence")
    else:
        ax.set_title(r"$\log_{10}$(" + var + ")")

    return


def plot_wind(root_name: str, output_name: str, path: str = "./", vars: List[str] = None,
              var_type: List[str] = None, projection: str = "rectilinear", subplot_dims: Tuple[int, int] = None,
              plot_title: str = None, loglog_scale: bool = True, filetype: str = "png",
              data_ndims: str = "2d", verbose: bool = False) -> None:
    """
    Create a 2D wind plot of the wind variables given in the list vars, of var
    type given in the list var_type. This function will only work with 2d Python
     simulations - who even uses 1d in Python anyway? Other than for sn runs.
    There is some error checking going on for inputs - this should be enough to
    make sure that something is plotted but still may cause the script to fall
    over, and will force exit in some situations.

    Parameters
    ----------
    root_name           str
                        The root name of the Python simulation
    output_name         str
                        The base name of the output plot
    path                str
                        The directory containing the Python simulation
    vars                list[str], optional
                        The Python wind variables to plot
    var_type            list[str], optional
                        The type of the variables to be plotted, allowable
                        values are wind and ion
    coords              str
                        The coordinate system in use for the wind data
    subplot_dims        tuple[int, int], optional
                        The number of rows and columns of subplot panels
    plot_title          str, optional
                        The title of the plot
    loglog_scale        bool, optional
                        Plot using a log log scale, set to True by default
    filetype            str, optional
                        The file type of the output plot saved to disk, set to
                        png by default
    data_ndims          str, optional
                        The dimensionality of the Python simulation - note only
                        2d is supported for now
    verbose             bool, optional
                        Enable verbose logging

    Returns
    -------
    None
    """

    if data_ndims != "2d":
        print("py_plot.plot_wind: only understand ndims = 2d")
        return

    # If no vars are provided, then use some default ones
    if vars is None:
        vars = ["t_e", "t_r", "ne", "v_x", "v_y", "v_z", "ip", "c4"]
    if var_type is None:
        var_type = ["wind"] * len(vars)
    if len(vars) != len(var_type):
        print("py_plot.plot_wind: vars and types should be of the same length")
        return

    # If subplot dimensions are not provided, then assume a default shape of (4, 2)
    if subplot_dims is None:
        subplot_dims = (4, 2)
    if subplot_dims[0] * subplot_dims[1] < len(vars):
        print("py_plot.plot_python_wind: not enough panels to plot all the provided vars!")
        return

    allowed_projections = ["rectilinear", "polar"]
    if projection not in allowed_projections:
        print("py_plot.plot_wind: projection {} not understood, allowed values: rectilinear or polar".format(projection))
        return

    if projection == "rectilinear":
        fig, ax = plt.subplots(subplot_dims[0], subplot_dims[1], figsize=(10, 15), squeeze=False)
    elif projection == "polar":
        fig = plt.figure(figsize=(10, 15))

    index = 0
    for i in range(subplot_dims[0]):
        for j in range(subplot_dims[1]):
            if index > len(vars) - 1:
                break
            var = vars[index]
            var_t = var_type[index]
            if verbose:
                print("\tPlotting {} of type {}".format(var, var_t))
            x, z, qoi = py_plot_util.get_wind_data(root_name, var, var_t, path, projection)

            # If the data return is only of length 1, then something has gone wrong in the
            # get_wind_data function, so increment the index counter and skip plotting this
            if len(qoi) == 1 and len(x) == 1 and len(z) == 1:
                index += 1
                continue

            if projection == "rectilinear":
                fig, ax = plot_rectilinear_wind(fig, ax, x, z, qoi, i, j, var, var_t, loglog_scale)
            elif projection == "polar":
                plot_polar_wind(x, z, qoi, index, var, var_t, subplot_dims, verbose)

            index += 1

    # Finishing touches of the plot
    if plot_title:
        fig.suptitle(plot_title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}/{}.{}".format(path, output_name, filetype))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def plot_spec_comps(spec_patch: str, output_name: str, loglog_scale: bool = False, smooth: int = 15,
                    filetype: str = "png", verbose: bool = False) -> None:
    """
    Plot the different integrated flux components within a Python spec file.

    Parameters
    ----------
    spec_path       str
                    A directory path to Python .spec file
    output_name     str
                    The base name for the plot which is saved to disk
    loglog_scale    bool, optional
                    Use a log-log scale for the plot
    smooth          int, optional
                    The size of the window for the boxcar smoother. Larger
                    numbers result in more smoothing
    filetype        str, optional
                    The file type of the plot saved to disk, by default this is
                    png
    verbose         bool, optional
                    Enable verbose logging

    Returns
    -------
    None
    """

    if type(spec_patch) is not str:
        print("py_plot.plot_spec_comps: can only plot spectrum components for one spectrum at a time")
        return

    fig, ax = plt.subplots(2, 1, figsize=(12, 10))

    # These are the headers which should be in Python spec files. They are the different components
    # of the Python spectra
    headers_top = ["Created", "Emitted"]
    headers_bot = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    # Get the spectrum for the model
    spec = py_plot_util.read_spec_file(spec_patch, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    # First, plot the created and emitted emission
    for i in range(len(headers_top)):
        print("\tPlotting {}".format(headers_top[i]))
        flux = py_plot_util.smooth_1d_array(np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float), smooth,
                                            verbose)
        # I'm using plot by default as semilogy can make some very ugly plots
        if loglog_scale:
            ax[0].loglog(wavelength, flux, label=headers_top[i])
        else:
            ax[0].plot(wavelength, flux, label=headers_top[i])
    ax[0].set_xlim(wavelength.min(), wavelength.max())
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[0].legend()

    # Second, plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bot)):
        print("\tPlotting {}".format(headers_bot[i]))
        flux = py_plot_util.smooth_1d_array(np.array(spec[1:, spec[0, :] == headers_bot[i]], dtype=float), smooth,
                                            verbose)
        if loglog_scale:
            ax[1].loglog(wavelength, flux, label=headers_bot[i])
        else:
            ax[1].plot(wavelength, flux, label=headers_bot[i])
    ax[1].set_xlim(wavelength.min(), wavelength.max())
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.{}".format(output_name, filetype))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def plot_spectra(spec_path: List[str], inclinations: Union[List, np.array], output_name: str, wmin: float = None,
                 wmax: float = None, smooth: int = 15, filetype: str = "png", verbose: bool = False) -> None:
    """
    Plot the spectrum for each provided inclination angle for a Python .spec file.

    Parameters
    ----------
    spec_path       str
                    A directory path to Python .spec file
    inclinations    list[str]
                    A list of inlcination angles to plot
    output_name     str
                    The base name for the plot which is saved to disk
    wmin            float, optional
                    The smallest wavelength to show on the plot
    wmax            float, optional
                    The largest wavelength to show on the plot
    smooth          int, optional
                    The size of the window for the boxcar smoother. Larger
                    numbers result in more smoothing
    filetype        str, optional
                    The file type of the plot saved to disk, by default this is
                    png
    show_plot       bool, optional
                    Show the plot before saving to disk
    verbose         bool, optional
                    Enable verbose logging

    Returns
    -------
    None
    """

    root = "spec"

    # ymin = 0.1
    # ymax = 1

    # Loop over each possible viewing angle as provided in inclinations
    for angle in inclinations:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        # Loop over each possible .spec file provided in spec_path
        ymin = +1e99
        ymax = -1e99
        for file in spec_path:
            root, filepath = py_plot_util.parse_root_name_and_path(file)
            legend = filepath + root
            if verbose:
                print("\tPlotting {} {}°".format(legend, angle))

            # Read in the spectrum and check that it can be plotted for the current viewing angle
            spec = py_plot_util.read_spec_file(file)
            allowed = py_plot_util.check_inclination_present(angle, spec)
            if not allowed:
                continue

            # Weird hacky and stupid code to find the index for the inclination
            # This is only required from read_spec_file returns a numpy array of strings
            idx = 0
            for i in range(spec.shape[1]):
                if spec[0, i].isdigit() and float(spec[0, i]) == float(angle):
                    spec[0, i] = float(spec[0, i])
                    idx = i
                    break

            # Extract the wavelength range and flux and scale the flux appropriately
            wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
            flux = py_plot_util.smooth_1d_array(np.array(spec[1:, idx], dtype=float), smooth, verbose)
            flux *= (DEFAULT_DIST ** 2 / OBSERVE_DIST ** 2)

            ax.semilogy(wavelength, flux, label=legend)
            ax.set_xlim(wmin, wmax)

            tymin, tymax = ax.get_ylim()

            if tymin < ymin:
                ymin = tymin / 2
            if tymax > ymax:
                ymax = tymax * 2

        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=15)
        ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)
        ax.legend(loc="lower right")
        ax.set_title("{} {}".format(root, angle) + r"$^{\circ}$", fontsize=20)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig("{}_{}.{}".format(output_name, angle, filetype))

        if SHOW_PLOT:
            plt.show()
        else:
            plt.close()

    return


def main() -> None:
    """
    The main controlling function of the script.

    Returns
    -------
    None
    """

    allowed_plots = ["spec", "spec_comp", "wind", "ions", "all"]

    # Parse the running options from the command line
    outname = get_script_arguments()

    if POLAR_PROJECTION:
        projection = "polar"
    else:
        projection = "rectilinear"

    print("--------------------------\n")

    if PLOTS not in allowed_plots:
        print("Don't know how to plot {}".format(PLOTS))
        print("Allowed plots are: spec, spec_comp, wind, ion or all")
        print("\n--------------------------")
        return

    # Get the output name, the viewing angles and the file paths to the .spec files
    spec_files = py_plot_util.find_spec_files()
    if len(spec_files) == 0:
        print("No spec files found")
        print("\n--------------------------")
        exit(1)

    print("Creating {} plots for the following simulations:\n".format(PLOTS))
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    print("\n--------------------------")

    # Plot spectra
    spec_angles = py_plot_util.spec_inclinations_numpy(spec_files)
    if PLOTS == "spec" or PLOTS == "all":
        print("\nPlotting spectra".format(spec_files))
        plot_spectra(spec_files, spec_angles, outname, wmin=WMIN, wmax=WMAX, smooth=SMOOTH, filetype=FILETYPE,
                     verbose=VERBOSE)

    # If this is being run in an individual folder, then we can plot the spectrum components and wind parameters
    if len(spec_files) == 1:
        root, path = py_plot_util.parse_root_name_and_path(spec_files[0])
        if PLOTS == "spec_comp" or PLOTS == "all":
            print("\nPlotting spectrum components")
            plot_spec_comps(spec_files[0], outname, loglog_scale=SPEC_LOGLOG, smooth=SMOOTH, filetype=FILETYPE,
                            verbose=VERBOSE)

        # Run windsave2table to extract data from the wind_save file
        if os.path.isfile("{}.ep.complete".format(root)) is False:
            rc = py_plot_util.run_windsave2table(path, root, VERBOSE)
            if rc:
                print("py_plot.main: windsave2table failed to run")
                exit(1)

        # Plot some wind quantities first
        if PLOTS == "wind" or PLOTS == "all":
            print("\nPlotting wind quantities")
            vars = ["t_e", "t_r", "ne", "converge", "w", "ntot", "ip", "c4"]
            var_types = ["wind"] * len(vars)
            plot_wind(root, outname + "_wind", path, vars, var_types, projection=projection, filetype=FILETYPE,
                      verbose=VERBOSE)

        # Now plot some ions
        if PLOTS == "ions" or PLOTS == "all":
            print("\nPlotting wind ions")
            vars = ["H_i01", "H_i02", "C_i03", "C_i04", "C_i05", "Si_i04", "N_i05", "O_i06"]
            var_types = ["ion"] * len(vars)
            plot_wind(root, outname + "_ions", path, vars, var_types, projection=projection, filetype=FILETYPE,
                      verbose=VERBOSE)

        print("")  # spacer :-)
        py_rm_data.remove_data_dir(path)

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
