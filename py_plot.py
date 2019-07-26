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
"""

from sys import exit
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
PLOT_LOG = False
DIMS = "2d"
WMIN = None
WMAX = None
FILETYPE = "png"
SMOOTH = 15
DEFAULT_DIST = 100 * PARSEC
OBSERVE_DIST = 100 * PARSEC
MIN_FLUX = 1e-20
ROOT = False


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
    global PLOT_LOG
    global WMIN
    global WMAX
    global FILETYPE
    global SMOOTH
    global OBSERVE_DIST
    global DIMS
    global ROOT

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("plots", nargs="?", type=str, help="The type of plot to create")
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help="The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help="The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help="Distance of the observer")
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

    return args.output_name


def plot_rectilinear_wind(fig, ax, x, z, qoi, i, j, var, var_t, loglog_scale) -> Tuple[plt.axes, plt.axes]:
    """
    Create a wind plot in rectilinear coordinates.

    Parameters
    ----------
    fig                 plt.Figure
                        A Pyplot figure object containing the axes
    ax                  plt.Axes
                        The plt.Axes object containing the subplots
    x                   np.array[float]
                        The x coordinates to plot
    z                   np.array[float]
                        The z coordinates to plot
    qoi                 np.array[float] (masked array)
                        The quantity of interest to be plotted
    i                   int
                        The i (row) index for the plt.Axes object
    j                   int
                        The j (column) index for the plt.Axes object
    vars                list[str], optional
                        The Python wind variables to plot
    var_type            list[str], optional
                        The type of the variables to be plotted, allowable
                        values are wind and ion
    loglog_scale        bool, optional
                        If True, plot on a loglog scale

    Returns
    -------
    fig

    ax
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


def plot_polar_wind(r, theta, qoi, index, var, var_t, subplot_dims):
    """
    Create a wind plot in a polar coordinates.

    Parameters
    ----------
    r                   np.array[float]
                        The x coordinates to plot
    theta               np.array[float]
                        The z coordinates to plot
    qoi                 np.array[float] (masked array)
                        The quantity of interest to be plotted
    index               int
                        The subplot index, i.e. (nrows, ncols, index)
    vars                list[str], optional
                        The Python wind variables to plot
    var_type            list[str], optional
                        The type of the variables to be plotted, allowable
                        values are wind and ion
    subplot_dims        tuple[int, int], optional
                        The number of rows and columns of subplot panels
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
    ax. set_thetamax(90)
    rmin = r[1][0]
    rmax = r[-2][0]
    ax.set_rlim(np.log10(rmin), np.log10(rmax))

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
        print(
                "py_plot.plot_wind: projection {} not understood, allowed values: rectilinear or polar".format(
                        projection))
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
                plot_polar_wind(x, z, qoi, index, var, var_t, subplot_dims)

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


def plot_tau_spec(root: str, dir: str = "./", plot_freq: bool = False, plot_edges: bool = True,
                  wmin: Union[float, bool] = None, wmax: Union[float, bool] = None,
                  semilogy: bool = True, loglog: bool = False) -> plt.Axes:
    """
    Plot an optical depth spectrum of either optical depth vs wavelength (Angstroms) or
    optical depth vs frequency (Hz).

    Parameters
    ----------
    root         str
                 The root name of the optical depth spectrum to plot
    dir          str
                 The directory containing the simulation
    plot_freq    bool, optional
                 If True, the optical depth will be printed in frequency space
    plot_edges   bool, optional
                 Label various absorption edges
    wmin         float, optional
                 The minimum wavelength to plot
    wmax         float, optional
                 The maximum wavelength to plot
    semilogy     bool, optional
                 If True, then the y axis will be log scaled
    loglog       bool, optional
                 If True, then the x and y axis will be log scaled

    Returns
    -------
    ax           plt.Axes
                 A matplotlib Axes object of the optical depth spectrum
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
            print("Unable to find the optical depth spectrum: {}".format(root))
            return fig, ax

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

    ymin = 1e+99
    ymax = 1e-99

    if plot_edges:
        scale = 50000
    else:
        scale = 10

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

        if plot_edges or semilogy is False:
            tymax, tymin = py_plot_util.get_ylimits(spec[:, 0], spec[:, i + 1], xxlims[0], xxlims[1], scale=scale)
            if tymax > ymax:
                ymax = tymax
            if tymin < ymin and tymin != 0:
                ymin = tymin

    if plot_freq:
        if loglog:
            ax.set_xlabel(r"Log(Frequency), Hz")
        else:
            ax.set_xlabel(r"Frequency, Hz")
    else:
        ax.set_xlabel(r"Wavelength, $\AA$", fontsize=15)
    ax.set_ylabel(r"Optical Depth, $\tau$", fontsize=15)

    if plot_edges:
        py_plot_util.plot_line_ids(ax, py_plot_util.get_common_absorption_edges(plot_freq, loglog), "horizontal")

    ax.legend()

    if plot_edges:
        py_plot_util.plot_line_ids(ax, py_plot_util.get_common_absorption_edges(plot_freq), "horizontal")

    plt.savefig("{}_tau_spec.png".format(root))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return ax


def plot_spec_comps(spec_path: str, output_name: str, semilogy_scale: bool = False, smooth: int = 15,
                    wmin: Union[float, bool] = None, wmax: Union[float, bool] = None, filetype: str = "png",
                    verbose: bool = False) -> None:
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
    wmin            float, optional
                    The minimum wavelength to plot
    wmax            float, optional
                    The maximum wavelength to plot
                    png
    verbose         bool, optional
                    Enable verbose logging

    Returns
    -------
    None
    """

    if type(spec_path) is not str:
        print("py_plot.plot_spec_comps: can only plot spectrum components for one spectrum at a time")
        return

    if spec_path.find(".spec") == -1:
        if spec_path.find(".pf") != -1:
            spec_path = spec_path.replace(".pf", ".spec")
        else:
            print("py_plot.plot_spec_comps: can't fix file path {}".format(spec_path))

    fig, ax = plt.subplots(2, 1, figsize=(12, 10))

    # These are the headers which should be in Python spec files. They are the different components
    # of the Python spectra
    headers_top = ["Created", "Emitted"]
    headers_bot = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    if spec_path.find(".spec") == -1:
        spec_path += ".spec"

    # Get the spectrum for the model
    spec = py_plot_util.read_spec_file(spec_path, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    # First, plot the created and emitted emission
    for i in range(len(headers_top)):
        print("\tPlotting {}".format(headers_top[i]))
        flux = py_plot_util.smooth_1d_array(np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float), smooth,
                                            verbose)
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

    # Second, plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bot)):
        print("\tPlotting {}".format(headers_bot[i]))
        flux = py_plot_util.smooth_1d_array(np.array(spec[1:, spec[0, :] == headers_bot[i]], dtype=float), smooth,
                                            verbose)
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
    nfiles = len(spec_path)

    for i in range(nfiles):
        path = spec_path[i]
        if path.find(".spec") != -1:
            continue
        elif path.find(".pf") != -1:
            path = path.replace(".pf", ".spec")
        spec_path[i] = path

    # Loop over each possible viewing angle as provided in inclinations
    for angle in inclinations:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        # Loop over each possible .spec file provided in spec_path
        ymin = +1e99
        ymax = -1e99
        for file in spec_path:
            if file.find(".spec") == -1:
                file += ".spec"
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
            tymax, tymin = py_plot_util.get_ylimits(wavelength, flux, wmin, wmax, scale=5)

            if tymin == 0 or tymax == 0:
                ymin = None
                ymax = None
            else:
                if tymin < ymin:
                    ymin = tymin
                if tymax > ymax:
                    ymax = tymax

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(wmin, wmax)
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

    allowed_plots = ["spec", "spec_comps", "wind", "ions", "tau_spec", "all", "help"]

    # Parse the running options from the command line
    outname = get_script_arguments()

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
        files = [outname]
    else:
        if PLOTS == "spec" or PLOTS == "spec_comps":
            files = py_plot_util.find_spec_files()
            if len(files) == 0:
                print("No spec files found")
                print("\n--------------------------")
                exit(1)
        else:
            files = py_plot_util.find_pf_files()
            if len(files) == 0:
                print("No pf files found")
                print("\n--------------------------")
                exit(1)

    print("Creating {} plots for the following simulations:\n".format(PLOTS))
    for i in range(len(files)):
        print("\t- {}".format(files[i]))

    print("\n--------------------------")

    # Plot spectra for inclination angles
    if PLOTS == "spec" or PLOTS == "all":
        spec_angles = py_plot_util.spec_inclinations_numpy(files)
        print("\nPlotting spectra".format(files))
        plot_spectra(files, spec_angles, outname, wmin=WMIN, wmax=WMAX, smooth=SMOOTH, filetype=FILETYPE,
                     verbose=VERBOSE)

    # If this is being run in an individual folder, then we can plot the spectrum components and wind parameters
    if len(files) == 1:
        root, path = py_plot_util.parse_root_name_and_path(files[0])

        # Plot the spectrum components
        if PLOTS == "spec_comps" or PLOTS == "all":
            print("\nPlotting spectrum components")
            plot_spec_comps(files[0], outname, semilogy_scale=PLOT_LOG, smooth=SMOOTH, filetype=FILETYPE,
                            verbose=VERBOSE)

        # Plot the optical depth spectrum
        if PLOTS == "tau_spec" or PLOTS == "all":
            print("\nPlotting optical depth spectrum")
            plot_tau_spec(root, path, wmin=WMIN, wmax=WMAX)

        # Run windsave2table to extract data from the wind_save file
        if PLOTS == "wind" or PLOTS == "ions" or PLOTS == "all":
            if os.path.isfile("{}.ep.complete".format(root)) is False:
                rc = py_plot_util.run_windsave2table(path, root, VERBOSE)
                if rc:
                    print("py_plot.main: windsave2table failed to run")
                    exit(1)

        # Plot some wind quantities first
        if PLOTS == "wind" or PLOTS == "all":
            print("\nPlotting wind quantities")
            vars = ["t_e", "t_r", "ne", "rho", "w", "converge", "ip", "ntot"]
            var_types = ["wind"] * len(vars)
            plot_wind(root, outname + "_wind", path, vars, var_types, projection=projection, filetype=FILETYPE,
                      data_ndims=DIMS, verbose=VERBOSE)

        # Now plot some ions
        if PLOTS == "ions" or PLOTS == "all":
            print("\nPlotting wind ions")
            # vars = ["He_i01", "He_i02", "He_i03"]
            # vars = ["C_i01", "C_i02", "C_i03", "C_i04", "C_i05", "C_i06"]
            # vars = ["N_i01", "N_i02", "N_i03", "N_i04", "N_i05", "N_i06", "N_i07", "N_i08"]
            vars = ["H_i01", "H_i02", "He_i01", "He_i02", "He_i03", "Si_i04", "N_i05", "C_i04"]
            var_types = ["ion"] * len(vars)
            plot_wind(root, outname + "_ions", path, vars, var_types, projection=projection, filetype=FILETYPE,
                      data_ndims=DIMS, verbose=VERBOSE)

        py_rm_data.remove_data_dir(path)
    elif len(files) > 1 and PLOTS != "spec":
        print("Can only plot {} when one root in folder :^)".format(PLOTS))
    print("")  # spacer :-)

    print("--------------------------")

    return


if __name__ == "__main__":
    main()
