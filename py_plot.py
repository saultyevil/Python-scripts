#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plot the output from a Python simulation.

The script can be called within a directory containing a single Python simulation, or in a directory containing multiple
Python simulations. In the case of being called in a directory with multiple Python simulations, only spectra will be
plotted. However, if called in a directory containing a single simulation, then the spectra, integration flux components,
various wind quantities and wind ions will be plotted.

You do not have to include the Python root name as an argument for this script, as it uses the Unix find command to search
for .spec files recursively to find Python simulations. TODO: switch to using pfs instead of spec files

usage: py_plot.py [-h] [-wmin WMIN] [-wmax WMAX] [-filetype FILETYPE]
                  [-smooth SMOOTH] [-dist DIST] [-v] [-s] [-z Z] [-tde]
                  [-loglog]
                  output_name

positional arguments:
  output_name         The base name for the output

optional arguments:
  -h, --help          show this help message and exit
  -wmin WMIN          The smallest wavelength to show
  -wmax WMAX          The largest wavelength to show
  -filetype FILETYPE  The file format of the output
  -smooth SMOOTH      The amount of smoothing of the spectra
  -dist DIST          Distance of the observer
  -v, --verbose       Increase output to screen
  -s, --show          Show the plots on screen
  -z Z                The redshift of the object
  -tde                Overplot iPTF15af UV spectrum and scale flux appropriately
  -loglog             Enable log log axes

TODO: add check to see if wind files already exist to avoid issues with windsave2table
"""


import argparse
import numpy as np
import py_plot_util
from matplotlib import pyplot as plt
from typing import List, Tuple, Union


VERBOSE = False                # More info will be printed to screen if True
SHOW_PLOT = False              # If True, the plot will be drawn on screen
TDE_PLOT = False               # Enable default TDE plotting
SPECLOGLOG = False             # Enable loglog axes on spectra
WMIN = None                    # The smallest wavelength to show on the spectra
WMAX = None                    # The largest wavelength to show on the spectra
FILETYPE = "png"               # The file type of the output spectra
SMOOTH = 15                    # The amount of smoothing for the spectra
OBSERVE_DIST = 100 * 3.086e18  # 100 pc in cm - the default Python distance
Z = 0                          # Redshift


def get_script_arguments() -> str:
    """
    Parse the various global parameters from the command line.

    Parameters
    ----------
    None

    Returns
    -------
    str:
        The output base name of the created plots.
    """

    global VERBOSE
    global SHOW_PLOT
    global SPECLOGLOG
    global WMIN
    global WMAX
    global FILETYPE
    global SMOOTH
    global CHANGE_DIST
    global OBSERVE_DIST
    global TDE_PLOT
    global Z

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("-wmin", type=float, action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help="The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help="The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help="Distance of the observer")
    p.add_argument("-v", "--verbose", action="store_true", help="Increase output to screen")
    p.add_argument("-s", "--show", action="store_true", help="Show the plots on screen")
    p.add_argument("-z", type=float, action="store", help="The redshift of the object")
    p.add_argument("-tde", action="store_true", help="Overplot iPTF15af UV spectrum and scale flux appropriately")
    p.add_argument("-loglog", action="store_true", help="Enable log log axes")
    args = p.parse_args()

    # Assign the optional arguments to their global vars
    if args.verbose:
        VERBOSE = True
    if args.show:
        SHOW_PLOT = True
    if args.loglog:
        SPECLOGLOG = True
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
    if args.z:
        Z = args.z
    if args.tde:
        TDE_PLOT = True
        OBSERVE_DIST = 1.079987153448e+27  # 350 Mpc
        Z = 0.07897
        if not WMIN:
            WMIN = 800
        if not WMAX:
            WMAX = 3000

    return args.output_name


def plot_python_wind(root_name:str, output_name:str, path:str="./", vars:List[str]=None, types:List[str]=None,
                     shape:Tuple[int, int]=None, title:str=None, loglog:bool=True, filetype:str="png",
                     plot_show:bool=False, verbose:bool=False, ndims:str="2d") -> None:
    """
    Create a 2D wind plot of the wind variables given in vars.

    Parameters
    ----------
    root_name: str
        The root name of the Python simulation
    output_name: str
        The base name for the output plot
    path: str, optional
        The directory containing the Python simulation
    vars: list of str, optional
        The variables to be plotted. If none are provided, default ones are used instead
    types: list of str, optional
        The type of the variables to be plotted. Allowed values are wind or ion. One is required for each variable
    shape: tuple of two ints, optional
        The number of subplots panels to create, given as (n_rows, n_cols)
    title: str, optional
        The title of the plot, placed at the very top
    loglog: bool, optional
        If True, log log scaling will be used for the x and y axes
    filtype: str, optional
        The file type of the image to create, by default this is png
    plot_show: bool, optional
        If True, the plot will be show before saved
    verbose: bool, optional
        If True, more verbose output will be shown
    ndims: str, optional
        The number of dimensions of the Python simulation. Currently only 2d is supported

    Returns
    -------
    None
    """

    if ndims == "1d":
        print("py_plot.plot_python_wind: 1d Python runs are not supported yet")
        exit(0)
    elif ndims != "2d":
        print("py_plot.plot_python_wind: only understand ndims = 1d or ndims = 2d")

    if vars is None:
        vars = ["t_e", "t_r", "ne", "v_x", "v_y", "v_z", "ip", "c4"]
    if types is None:
        types = ["wind", "wind", "wind", "wind", "wind", "wind", "wind", "wind"]
    if shape is None:
        shape = (4, 2)

    if shape[0] * shape[1] < len(vars):
        print("py_plot.plot_python_wind: not enough panels to plot all the provided vars!")
        return

    if len(vars) != len(types):
        print("py_plot.plot_python_wind: vars and types should be of the same length")
        return

    idx = 0
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(8, 15), squeeze=False)
    for i in range(shape[0]):
        for j in range(shape[1]):
            var = vars[idx]
            var_type = types[idx]
            x, z, qoi = py_plot_util.get_wind_data(root_name, var, var_type)
            if qoi.all() == 0 and x.all() == 0 and z.all() == 0:
                idx += 1
                continue
            with np.errstate(divide="ignore"):
                if var_type.lower() == "ion":
                    im = ax[i, j].pcolor(x, z, np.log10(qoi), vmin=-5, vmax=0)
                elif var_type.lower() == "wind":
                    im = ax[i, j].pcolor(x, z, np.log10(qoi))
                else:
                    print("py_plot.plot_python_wind: type {} not recognised".format(type))
            fig.colorbar(im, ax=ax[i, j])
            ax[i, j].set_xlabel("x")
            ax[i, j].set_ylabel("z")
            ax[i, j].set_title(r"$\log_{10}$("+var+")")
            if loglog:
                ax[i, j].set_xscale("log")
                ax[i, j].set_yscale("log")
                ax[i, j].set_xlim(x[1, 1], x[-1, -1])
                ax[i, j].set_ylim(z[1, 1], z[-1, -1])
            idx += 1

    if title:
        fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("{}/{}.{}".format(path, output_name, filetype))

    if plot_show:
        plt.show()
    else:
        plt.close()

    return


def plot_spec_comps(spec_file:str, outname:str, loglog:bool=False, smooth:int=15, verbose:bool=False,
                    filetype:str="png", plot_show:bool=False) -> None:
    """
    Plot the components which makes the integrated flux of a Python Simulation

    Parameters
    ----------
    spec_file: str
        The spec file of the Python simulation to plot. Include the full path.
    outname: str
        The base name of the output plot
    loglog: bool, optional
        If True, log log scales will be used instead of semilogy
    smooth: int, optional
        The number of sample points to use in the boxcar smoother
    verbose: bool, optional
        If True, more verbose output will be shown
    filetype: str, optional
        The file type of the output image. By default, this is png
    plot_show: bool, optional
        If True, show the created plot before saving it

    Returns
    -------
    None
    """

    if type(spec_file) is not str:
        print("py_plot.plot_spec_comps: can only plot spectrum components for one spectrum at a time")
        return

    headers_top = ["Created", "Emitted"]
    headers_bottom = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    spec = py_plot_util.read_spec_file(spec_file, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    # Plot created and emitted emission
    for i in range(len(headers_top)):
        flux = py_plot_util.smooth_flux(np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float), smooth, verbose)
        if loglog:
            ax[0].loglog(wavelength, flux, label=headers_top[i])
        else:
            ax[0].semilogy(wavelength, flux, label=headers_top[i])
    ax[0].set_xlim(wavelength.min(), wavelength.max())
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[0].legend()

    # Plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bottom)):
        flux = py_plot_util.smooth_flux(np.array(spec[1:, spec[0, :] == headers_bottom[i]], dtype=float), smooth,
                                        verbose)
        if loglog:
            ax[1].loglog(wavelength, flux, label=headers_bottom[i])
        else:
            ax[1].semilogy(wavelength, flux, label=headers_bottom[i])
    ax[1].set_xlim(wavelength.min(), wavelength.max())
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.{}".format(outname, filetype))

    if plot_show:
        plt.show()
    else:
        plt.close()

    return


def plot_spectra(spec_files:List[str], spec_angles:Union[List, np.array], outname:str, wmin:float=None, wmax:float=None,
                 plot_iPTF15af:bool=False, smooth:int=15, plot_show:bool=False, filetype:str="png",
                 verbose:bool=False) -> None:
    """
    Plot the spectrum of each viewing angle of a Python simulation.

    Parameters
    ----------
    spec_file: str
        The spec file of the Python simulation to plot. Include the full path.
    spec_angles: list of ints
        The spectrum viewing angles which will be plotted
    outname: str
        The base name of the output plot
    wmin: float, optional
        The smallest wavelength to show
    wmax: float, optional
        The longest wavelength to show
    plot_iPTF15af: bool, optional
        If True, then the iPTF15af UV spectrum will be overplotted
    smooth: int, optional
        The number of sample points to use in the boxcar smoother
    plot_show: bool, optional
        If True, show the created plot before saving it
    verbose: bool, optional
        If True, more verbose output will be shown
    filetype: str, optional
        The file type of the output image. By default, this is png

    Returns
    -------
    None
    """

    if plot_iPTF15af:
        iPTF15af_spec = py_plot_util.get_iPTF15af_spec(smooth, verbose)

    # Loop over all the possible viewing angles
    for angle in spec_angles:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        if plot_iPTF15af:
            ax.semilogy(iPTF15af_spec[:, 0] / (Z + 1), iPTF15af_spec[:, 1], label="iPTF15af: Blagorodnova et al. (2019)")

        # Loop over each possible .spec file
        for file in spec_files:
            root, filepath = py_plot_util.get_root_name_and_path(file)
            legend = filepath + root

            if VERBOSE:
                print("\tPlotting {} {}°".format(legend, angle))

            # Read in .spec file and check that it can be plotted for the current viewing angle
            spec = py_plot_util.read_spec_file(file)
            allowed = py_plot_util.check_viewing_angle(angle, spec)
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
            flux *= (default_dist**2 / OBSERVE_DIST**2)
            yupper, ylower = py_plot_util.get_ylims(wavelength, flux)

            ax.semilogy(wavelength, flux, label=legend)
            ax.set_xlim(wmin, wmax)
            ax.set_ylim(ylower, yupper)
            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
            ax.tick_params(labelsize=17)

        ax.legend(loc="best")
        ax.set_title("{} {}".format(root, angle) + "$^{\circ}$", fontsize=20)

        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig("{}_{}.{}".format(outname, angle, filetype))

        if plot_show:
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
    plot_spectra(spec_files, spec_angles, outname, wmin=WMIN, wmax=WMAX, plot_iPTF15af=TDE_PLOT, smooth=SMOOTH,
                 verbose=VERBOSE, plot_show=SHOW_PLOT, filetype=FILETYPE)

    # If this is being run in an individual folder, then we can plot the spectrum components and wind parameters
    if len(spec_files) == 1:
        root, path = py_plot_util.get_root_name_and_path(spec_files[0])
        print("\nPlotting spectrum components")
        plot_spec_comps(spec_files[0], outname, loglog=SPECLOGLOG, smooth=SMOOTH, verbose=VERBOSE, filetype=FILETYPE,
                        plot_show=SHOW_PLOT)
        # Run windsave2table to extract data from the wind_save file
        rc = py_plot_util.run_windsave2table(path, root, VERBOSE)
        if rc:
            print("py_plot.main: windsave2table failed, skipping wind and ion plots")
        else:
            # Plot some wind quantities first
            print("\nPlotting wind quantities")
            vars = ["t_e", "t_r", "ne", "v_x", "v_y", "v_z", "ip", "c4"]
            var_types = ["wind"] * len(vars)
            plot_python_wind(root, outname+"_wind", path, vars, var_types, filetype=FILETYPE, plot_show=SHOW_PLOT,
                             verbose=VERBOSE)
            # Now plot some ions
            print("\nPlotting wind ions")
            vars = ["H_i01", "H_i02", "C_i03", "C_i04", "C_i05", "Si_i04", "N_i05", "O_i06"]
            var_types = ["ion"] * len(vars)
            plot_python_wind(root, outname+"_ions", path, vars, var_types, filetype=FILETYPE, plot_show=SHOW_PLOT,
                             verbose=VERBOSE)

    print("\nDone!")
    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()
