#!/usr/bin/env python3

"""
Create spectra from the .spec files from a MCRT Python simulation.

The script requires you to provide 2 arguments:
    - The base output name of the spectra which are going to be produced
    - The viewing angles to create spectra for: possible choices for this
      include all, a single number or a list of numbers separated by a comma,
      i.e. 20,30,40

There are also some optional arguments:
    - dist: the distance of the object from the observer, by default Python
            uses 100 pc
    - wmin: the smallest wavelength to plot
    - wmax: the largest wavelength to plot
    - filetype: the file type of the output spectra plots, by default this is png
    - tde: default plotting parameter choices for modelling iPTF15af (TDE)
    - loglog: enable log log axes, otherwise semilogy

There are also two optional switches:
    - -v or --verbose which will enable more verbose output
    - -s or --show which will show the output spectra as well as saving to disk

This script works by using the unix find command to find .spec files recursively
in the current directory and deeper. Thus, it's possible to plot spectra for a
single simulation or by plotting multiple simulations at once on the same plot
for comparison purposes.

Example usage:
    python spec_plot.py qDNe_test all -wmin 1250 -wmax 3000 -v
"""


import py_util
import argparse
import numpy as np
from matplotlib import pyplot as plt


VERBOSE = False                # More info will be printed to screen if True
SHOW_PLOT = False              # If True, the plot will be drawn on screen
TDE_PLOT = False               # Enable default TDE plotting
LOGLOG = False                 # Enable loglog axes
WMIN = None                    # The smallest wavelength to show on the spectra
WMAX = None                    # The largest wavelength to show on the spectra
FILETYPE = "png"               # The file type of the output spectra
SMOOTH = 11                    # The amount of smoothing for the spectra
OBSERVE_DIST = 100 * 3.086e18  # 100 pc in cm - the default Python distance
Z = 0                          # Redshift


def get_script_arguments(specfiles):
    """
    Parse the various global parameters from the command line.
    """

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output")
    p.add_argument("angles", type=str, help=
        "The viewing angles to plot: all, a single angle or a comma "
        "separated list")
    p.add_argument("-wmin", type=float, action="store", help=
        "The smallest wavelength to show")
    p.add_argument("-wmax", type=float, action="store", help=
        "The largest wavelength to show")
    p.add_argument("-filetype", type=str, action="store", help=
        "The file format of the output")
    p.add_argument("-smooth", type=float, action="store", help=
        "The amount of smoothing of the spectra")
    p.add_argument("-dist", type=float, action="store", help=
        "Distance of the observer")
    p.add_argument("-v", "--verbose", action="store_true", help=
        "Increase output to screen")
    p.add_argument("-s", "--show", action="store_true", help=
        "Show the plots on screen")
    p.add_argument("-z", type=float, action="store", help=
        "The redshift of the object")
    p.add_argument("-tde", action="store_true", help=
        "Plot for the Blagodovnova (?) UV TDE spec")
    p.add_argument("-loglog", action="store_true", help="Enable log log axes")
    args = p.parse_args()

    global VERBOSE
    global SHOW_PLOT
    global LOGLOG
    global WMIN
    global WMAX
    global FILETYPE
    global SMOOTH
    global CHANGE_DIST
    global OBSERVE_DIST
    global TDE_PLOT
    global Z

    # Assign the optional arguments to their global vars
    if args.verbose:
        VERBOSE = True
    if args.show:
        SHOW_PLOT = True
    if args.loglog:
        LOGLOG = True
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
        OBSERVE_DIST = 1.079987153448e+27  # 350 MPc
        Z = 0.07897
        if not WMIN:
            WMIN = 1100
        if not WMAX:
            WMAX = 2500

    # Get the viewing angles to plot
    angles = []
    if args.angles == "all":
        angles = py_util.get_spec_viewing_angles(specfiles)
    elif args.angles.isdigit() is True:
        angles.append(int(args.angles))
    else:
        ang = args.angles.replace(",", " ").split()
        for i in range(len(ang)):
            angles.append(int(ang[i]))

    return args.output_name, angles


def plot_ions(spec_file, outname, user_ions=None, user_shape=None, loglog=True):
    """
    Create 2D plots of the different ions in a Python simulation
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot ion components for one simulation at a time")
        return
    spec_file = spec_file[0]

    root, path = py_util.get_root_name_and_path(spec_file)

    ions = default_ions = ["H_i01", "H_i02", "C_i03", "C_i04", "C_i05",
                           "Si_i04", "N_i05", "O_i06"]
    shape = default_shape = (4, 2)  # 4 rows, 2 cols

    fig, ax = plt.subplots(shape[0], shape[1], figsize=(12, 16), squeeze=False)

    var_idx = 0
    for i in range(shape[0]):
        for j in range(shape[1]):
            var = ions[var_idx]
            x, z, ion = py_util.get_ion_data(path, root, var, VERBOSE)
            im = ax[i, j].pcolor(x, z, np.log10(ion), vmin=-5, vmax=0)
            # Set the scales to log if that is the wish of the user
            if loglog:
                ax[i, j].set_xscale("log")
                ax[i, j].set_yscale("log")
                # Hacky fix for weird x and y limiting when using loglog
                ax[i, j].set_xlim(x[1, 1], x[-1, -1])
                ax[i, j].set_ylim(z[1, 1], z[-1, -1])
            # Add helpful labels and colourbar
            ax[i, j].set_xlabel("x")
            ax[i, j].set_ylabel("z")
            ax[i, j].set_title(r"$\log_{10}$("+var+")")
            fig.colorbar(im, ax=ax[i, j])
            var_idx += 1

    fig.tight_layout()
    plt.savefig("{}_ions_plot.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def plot_wind(spec_file, outname, user_vars=None, user_shape=None, loglog=True):
    """
    Create 2D plots of some wind parameters in a Python simulation. I
    essentially wrote this as I grew tired of the default py_progs plotting
    routine and this one uses windsave2table instead of py_wind so should
    actually be simpler.
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot wind components for one simulation at a time")
        return
    spec_file = spec_file[0]

    root, path = py_util.get_root_name_and_path(spec_file)
    wind = py_util.get_master_data(path, root, VERBOSE)

    # If no user vars or shape is provided, use the default
    vars = default_vars = ["ne", "t_e", "t_r", "ntot"]
    shape = default_shape = (2, 2)  # 2 rows, 2 columns
    if user_vars and user_shape:
        vars = user_vars
        shape = user_shape
    # In the case where only user_vars is provided, just use the default and
    # inform the user of their mistake
    elif user_vars and user_shape is None:
        print("plot_wind: as only user_vars has been provided with no shape"
              " the default vars and shape are being used instead as"
              " dynamic subplotting hasn't been implemented yet")
        vars = default_vars
        shape = default_shape

    # If we use squeeze=False, then ax is always 2d :-)
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(12, 12), squeeze=False)

    var_idx = 0
    for i in range(shape[0]):
        for j in range(shape[1]):
            var = vars[var_idx]
            x, z, qoi = py_util.get_wind_quantity(wind, var, VERBOSE)
            im = ax[i, j].pcolor(x, z, np.log10(qoi))
            # Set the scales to log if that is the wish of the user
            if loglog:
                ax[i, j].set_xscale("log")
                ax[i, j].set_yscale("log")
                # Hacky fix for weird x and y limiting when using loglog
                ax[i, j].set_xlim(x[1, 1], x[-1, -1])
                ax[i, j].set_ylim(z[1, 1], z[-1, -1])
            # Add helpful labels and colourbar
            ax[i, j].set_xlabel("x")
            ax[i, j].set_ylabel("z")
            ax[i, j].set_title(r"$\log_{10}$("+var+")")
            fig.colorbar(im, ax=ax[i, j])

            var_idx += 1

    fig.tight_layout()
    plt.savefig("{}_wind_plot.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def plot_spec_comps(spec_file, outname):
    """
    Plot the different components of the spectra
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot spectrum components for one spectrum at a time")
        return

    spec_file = spec_file[0]
    headers_top = ["Created", "Emitted"]
    headers_bottom = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    fig, ax = plt.subplots(1, 2, figsize=(12, 8))
    spec = py_util.read_spec_file(spec_file, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    for i in range(len(headers_top)):
        flux = np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float)
        smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
        ax[0].plot(wavelength, smoothflux, label=headers_top[i])
    ax[0].set_xlim(wavelength.min(), wavelength.max())
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)",
                          fontsize=17)
    ax[0].legend()

    for i in range(len(headers_bottom)):
        flux = np.array(spec[1:, spec[0, :] == headers_bottom[i]], dtype=float)
        smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
        ax[1].plot(wavelength, smoothflux, label=headers_bottom[i])
    ax[1].set_xlim(wavelength.min(), wavelength.max())
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)",
                          fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    return


def print_info (spec_files, angles):
    """
    Print extra information to the screen
    """

    print("Creating spectra for the following .spec files:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    print("\nSpectra will be created for the following inclinations:\n")
    for i in range(len(angles)):
        print("\t- {}".format(angles[i]))

    print("")

    return


def get_ylims(wlength, flux):
    """
    Figure out the ylims given a wavelength range
    """

    wmin_flux = flux.min()
    wmax_flux = flux.max()

    # As WMIN/WMAX and wlength are floats, we have to be smarter to find where
    # WMIN and WMAX are in an array of floats
    if WMIN:
        wmin_idx = np.abs(wlength - float(WMIN)).argmin()
        if VERBOSE:
            print("wmin_idx: {}".format(wmin_idx))
        wmin_flux = flux[wmin_idx]
    if WMAX:
        wmax_idx = np.abs(wlength - float(WMAX)).argmin()
        if VERBOSE:
            print("wmax_idx: {}".format(wmax_idx))
        wmax_flux = flux[wmax_idx]

    if wmin_flux > wmax_flux:
        yupper = wmin_flux
        ylower = wmax_flux
    else:
        yupper = wmax_flux
        ylower = wmin_flux

    yupper *= 10
    ylower /= 10

    # Revolting hack to ensure the Blag TDE spectrum isn't cut off by the lims
    if TDE_PLOT:
        if yupper < 1e-14:
            yupper = 1e-14

    return yupper, ylower


def plot_spectra(spec_files, spec_angles, outname):
    """
    Plot the spectra for the give .spec files and for the given viewing angles.
    """

    print_info(spec_files, spec_angles)

    if len(spec_angles) > 1:
        print("Plotting spectra now...\n")
    else:
        print("Plotting spectrum now...")

    if TDE_PLOT:
        blag_spec = py_util.get_blagordnova_spec(SMOOTH, VERBOSE)

    # Loop over all the possible viewing angles
    for angle in spec_angles:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        if TDE_PLOT:
            ax.semilogy(blag_spec[:, 0], blag_spec[:, 1],
                        label="Blagorodnova et al. 2018 (in prep.)")

        # Loop over each possible .spec file
        for file in spec_files:
            root, filepath = py_util.get_root_name_and_path(file)
            legend = filepath + root

            if VERBOSE:
                print("\t- {}: {}°".format(legend, angle))

            # Read in .spec file and check that it can be plotted for the
            # current viewing angle
            spec = py_util.read_spec_file(file)
            allowed = py_util.check_viewing_angle(angle, spec)
            if not allowed:
                if VERBOSE:
                    print("Error: {}° not found in .spec file {}"
                          .format(angle, file))
                continue

            # Read the wavelength and flux and smooth the flux
            # There is something weird happening to determine which index data
            # should be extracted from. I did this because all of the data in
            # the read in .spec file are strings
            # NOTE: this won't work super well with different phase angles

            idx = 0
            for i in range(spec.shape[1]):
                if spec[0, i].isdigit() and float(spec[0, i]) == float(angle):
                    # Found the index for the current angle to be plotted :-)
                    spec[0, i] = float(spec[0, i])
                    idx = i
                    break

            flux = np.array(spec[1:, idx], dtype=float)
            smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
            wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

            # Plot and scale flux for observer distance
            default_dist = 100 * 3.08567758128e18  # 100 pc
            flux_dist = smoothflux * (default_dist**2 / OBSERVE_DIST**2)
            z_wav = wavelength * (Z + 1)

            if LOGLOG:
                ax.loglog(z_wav, flux_dist, label=legend)
            else:
                ax.semilogy(z_wav, flux_dist, label=legend)

            ax.set_xlim(WMIN, WMAX)
            yupper, ylower = get_ylims(wavelength, flux_dist)
            ax.set_ylim(ylower, yupper)

            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)",
                          fontsize=17)
            ax.tick_params(labelsize=17)

        ax.legend(loc="best")
        ax.set_title("Viewing angle = {}".format(angle) + "$^{\circ}$",
                     fontsize=20)

        plt.savefig("{}_{}.{}".format(outname, angle, FILETYPE))

        if SHOW_PLOT:
            plt.show()
        else:
            plt.close()

    if len(spec_angles) > 1:
        print("Spectra plotted...")
    else:
        print("Spectrum plotted...")

    return


def main():
    """
    Main controlling function
    """

    print("--------------------------\n")

    # Get the output name, the viewing angles and the file paths to the .spec files
    spec_files = py_util.find_spec_files()
    if len(spec_files) == 0:
        print("No spec files found")
        exit(1)

    outname, spec_angles = get_script_arguments(spec_files)
    if len(spec_angles) == 0:
        print("No angles provided")
        exit(2)

    plot_spectra(spec_files, spec_angles, outname)
    plot_spec_comps(spec_files, outname)
    plot_wind(spec_files, outname)
    plot_ions(spec_files, outname)

    if VERBOSE:
        print("\nAll done :-)\n")
    print("--------------------------")

    return


if __name__ == "__main__":
    main()
