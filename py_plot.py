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
import tde_spec_plot

VERBOSE = False
SHOW_PLOT = False
SPEC_LOGLOG = False
WMIN = None
WMAX = None
FILETYPE = "png"
SMOOTH = 15
OBSERVE_DIST = 100 * PARSEC


def get_script_arguments() -> str:
    """
    Parse the various global parameters from the command line.

    Parameters
    ----------
    None

    Returns
    -------
    args.output_name       str
                           The output base name for plots
    """

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
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help="The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help="The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help="Distance of the observer")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output to screen")
    p.add_argument("-s", "--show", action="store_true", help="Show the plots on screen")
    p.add_argument("-loglog", action="store_true", help="Enable log log axes")
    args = p.parse_args()

    # Assign the optional arguments to their global vars
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

    return args.output_name


def plot_python_wind(root_name: str, output_name: str, path: str = "./", vars: List[str] = None,
                     var_type: List[str] = None, subplot_dims: Tuple[int, int] = None, plot_title: str = None,
                     loglog_scale: bool = True, filetype: str = "png", show_plot: bool = False, data_ndims: str = "2d",
                     verbose: bool = False) -> None:
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
    subplot_dims        tuple[int, int], optional
                        The number of rows and columns of subplot panels
    plot_title          str, optional
                        The title of the plot
    loglog_scale        bool, optional
                        Plot using a log log scale, set to True by default
    filetype            str, optional
                        The file type of the output plot saved to disk, set to
                        png by default
    show_plot           bool, optional
                        Show the plot before saving if set to True
    data_ndims          str, optional
                        The dimensionality of the Python simulation - note only
                        2d is supported for now
    verbose             bool, optional
                        Enable verbose logging

    Returns
    -------
    None
    """

    if data_ndims == "1d":
        print("py_plot.plot_python_wind: 1d Python runs are not supported yet")
        exit(1)
    elif data_ndims != "2d":
        print("py_plot.plot_python_wind: only understand ndims = 1d or ndims = 2d")
        exit(1)

    if vars is None:
        vars = ["t_e", "t_r", "ne", "v_x", "v_y", "v_z", "ip", "c4"]
    if var_type is None:
        var_type = ["wind"] * len(vars)
    if subplot_dims is None:
        subplot_dims = (4, 2)

    if subplot_dims[0] * subplot_dims[1] < len(vars):
        print("py_plot.plot_python_wind: not enough panels to plot all the provided vars!")
        return

    if len(vars) != len(var_type):
        print("py_plot.plot_python_wind: vars and types should be of the same length")
        return

    idx = 0
    fig, ax = plt.subplots(subplot_dims[0], subplot_dims[1], figsize=(10, 15), squeeze=False)
    for i in range(subplot_dims[0]):
        for j in range(subplot_dims[1]):
            var = vars[idx]
            var_t = var_type[idx]
            if verbose:
                print("\tPlotting {} of type {}".format(var, var_t))
            x, z, qoi = py_plot_util.get_wind_data(root_name, var, var_t)
            if len(qoi) == 1 and len(x) == 1 and len(z) == 1:
                idx += 1
                continue
            with np.errstate(divide="ignore"):
                if var == "converge":
                    im = ax[i, j].pcolor(x, z, qoi)
                elif var_t == "ion":
                    im = ax[i, j].pcolor(x, z, np.log10(qoi), vmin=-5, vmax=0)
                elif var_t == "wind":
                    im = ax[i, j].pcolor(x, z, np.log10(qoi))
                else:
                    print("py_plot.plot_python_wind: type {} not recognised".format(type))
                    continue
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
            idx += 1

    if plot_title:
        fig.suptitle(plot_title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}/{}.{}".format(path, output_name, filetype))

    if show_plot:
        plt.show()
    else:
        plt.close()

    return


def plot_spec_comps(spec_patch: str, output_name: str, loglog_scale: bool = False, smooth: int = 15,
                    filetype: str = "png", show_plot: bool = False, verbose: bool = False) -> None:
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
    show_plot       bool, optional
                    Show the plot before saving to disk
    verbose         bool, optional
                    Enable verbose logging

    Returns
    -------
    None
    """

    if type(spec_patch) is not str:
        print("py_plot.plot_spec_comps: can only plot spectrum components for one spectrum at a time")
        return

    headers_top = ["Created", "Emitted"]
    headers_bottom = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    spec = py_plot_util.read_spec_file(spec_patch, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    # Plot created and emitted emission
    for i in range(len(headers_top)):
        print("\tPlotting {}".format(headers_top[i]))
        flux = py_plot_util.smooth_flux(np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float), smooth, verbose)
        if loglog_scale:
            ax[0].loglog(wavelength, flux, label=headers_top[i])
        else:
            ax[0].plot(wavelength, flux, label=headers_top[i])
    ax[0].set_xlim(wavelength.min(), wavelength.max())
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[0].legend()

    # Plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bottom)):
        print("\tPlotting {}".format(headers_bottom[i]))
        flux = py_plot_util.smooth_flux(np.array(spec[1:, spec[0, :] == headers_bottom[i]], dtype=float), smooth,
                                        verbose)
        if loglog_scale:
            ax[1].loglog(wavelength, flux, label=headers_bottom[i])
        else:
            ax[1].plot(wavelength, flux, label=headers_bottom[i])
    ax[1].set_xlim(wavelength.min(), wavelength.max())
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.{}".format(output_name, filetype))

    if show_plot:
        plt.show()
    else:
        plt.close()

    return


def plot_spectra(spec_path: List[str], inclinations: Union[List, np.array], output_name: str, wmin: float = None,
                 wmax: float = None, smooth: int = 15, filetype: str = "png", show_plot: bool = False,
                 verbose: bool = False) -> None:
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

    # Loop over all the possible viewing angles
    for angle in inclinations:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        # Loop over each possible .spec file
        for file in spec_path:
            root, filepath = py_plot_util.get_root_name_and_path(file)
            legend = filepath + root

            if VERBOSE:
                print("\tPlotting {} {}°".format(legend, angle))

            # Read in .spec file and check that it can be plotted for the current viewing angle
            spec = py_plot_util.read_spec_file(file)
            allowed = py_plot_util.check_inclination_angle(angle, spec)
            if not allowed:
                continue

            # Weird hacky code to find the correct index because I refuse to use Pandas tables?
            idx = 0
            for i in range(spec.shape[1]):
                if spec[0, i].isdigit() and float(spec[0, i]) == float(angle):
                    spec[0, i] = float(spec[0, i])
                    idx = i
                    break

            wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
            flux = py_plot_util.smooth_flux(np.array(spec[1:, idx], dtype=float), smooth, verbose)

            # Scale the flux
            default_dist = 100 * 3.08567758128e18  # 100 pc
            flux *= (default_dist ** 2 / OBSERVE_DIST ** 2)
            yupper, ylower = py_plot_util.get_ylims(wavelength, flux, WMIN, WMAX, VERBOSE)

            ax.semilogy(wavelength, flux, label=legend)
            ax.set_xlim(wmin, wmax)
            ax.set_ylim(ylower, yupper)
            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=15)
            ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)

        ax.legend(loc="lower right")
        ax.set_title("{} {}".format(root, angle) + r"$^{\circ}$", fontsize=20)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig("{}_{}.{}".format(output_name, angle, filetype))

        if show_plot:
            plt.show()
        else:
            plt.close()

    return


def main() -> None:
    """
    The main controlling function of the script.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    # Parse the running options from the command line
    outname = get_script_arguments()

    print("--------------------------\n")

    # Get the output name, the viewing angles and the file paths to the .spec files
    spec_files = py_plot_util.find_spec_files()
    if len(spec_files) == 0:
        print("No spec files found")
        print("\n--------------------------")
        exit(1)

    print("Creating plots for the following simulations:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    print("\n--------------------------")

    # Plot spectra
    print("\nPlotting spectra".format(spec_files))
    spec_angles = py_plot_util.get_spec_viewing_angles(spec_files)
    plot_spectra(spec_files, spec_angles, outname, wmin=WMIN, wmax=WMAX, smooth=SMOOTH, filetype=FILETYPE,
                 show_plot=SHOW_PLOT, verbose=VERBOSE)

    # If this is being run in an individual folder, then we can plot the spectrum components and wind parameters
    if len(spec_files) == 1:
        root, path = py_plot_util.get_root_name_and_path(spec_files[0])
        print("\nPlotting spectrum components")
        plot_spec_comps(spec_files[0], outname, loglog_scale=SPEC_LOGLOG, smooth=SMOOTH, filetype=FILETYPE,
                        show_plot=SHOW_PLOT, verbose=VERBOSE)
        # Run windsave2table to extract data from the wind_save file
        if os.path.isfile("{}.ep.complete".format(root)) is False:
            rc = py_plot_util.run_windsave2table(path, root, VERBOSE)
            if rc:
                print("py_plot.main: windsave2table failed to run")
                exit(1)
        # Plot some wind quantities first
        print("\nPlotting wind quantities")
        vars = ["t_e", "t_r", "ne", "converge", "w", "ntot", "ip", "c4"]
        var_types = ["wind"] * len(vars)
        plot_python_wind(root, outname + "_wind", path, vars, var_types, filetype=FILETYPE, show_plot=SHOW_PLOT,
                         verbose=VERBOSE)
        # Now plot some ions
        print("\nPlotting wind ions")
        vars = ["H_i01", "H_i02", "C_i03", "C_i04", "C_i05", "Si_i04", "N_i05", "O_i06"]
        var_types = ["ion"] * len(vars)
        plot_python_wind(root, outname + "_ions", path, vars, var_types, filetype=FILETYPE, show_plot=SHOW_PLOT,
                         verbose=VERBOSE)
        print("")  # spacer :-)
        py_rm_data.remove_data_dir(path)

    print("\nDone!")
    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
