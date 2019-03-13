#!/usr/bin/env python3

"""
Create spectra from the .spec files from a MCRT Python simulation.

The script requires you to provide 1 argument:
    - The base output name of the spectra which are going to be produced

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
    python py_plot.py qDNe_test -wmin 1250 -wmax 3000 -v

TODO: the wind and ion plotting are basically the same - make a better plotting function for these
"""


import py_util
import argparse
import numpy as np
from matplotlib import pyplot as plt


VERBOSE = False                # More info will be printed to screen if True
SHOW_PLOT = False              # If True, the plot will be drawn on screen
TDE_PLOT = False               # Enable default TDE plotting
SPECLOGLOG = False             # Enable loglog axes on spectra
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
    p.add_argument("-tde", action="store_true", help="Plot for the Blagodovnova (?) UV TDE spec")
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
            WMIN = 1100
        if not WMAX:
            WMAX = 2500

    # Get the viewing angles to plot
    angles = py_util.get_spec_viewing_angles(specfiles)

    return args.output_name, angles


def plot_python_ions(spec_file, outname, user_ions=None, user_shape=None, loglog=True):
    """
    Create 2D plots of wind ions for a Python simulation. Note that your
    own ions and shape can be provided instead of the default ones. Note
    that both must be provided otherwise the default ones will  be used.

    NOTE: this uses windsave2table which must have been compiled from the same
    git commit otherwise odd things can happen
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot win ions for one simulation at a time")
        return
    spec_file = spec_file[0]

    root, path = py_util.get_root_name_and_path(spec_file)

    # Default ions to plot and the shape for the resulting figure
    ions = default_ions = ["H_i01", "H_i02", "C_i03", "C_i04", "C_i05", "Si_i04", "N_i05", "O_i06"]
    shape = default_shape = (4, 2)  # 4 rows, 2 cols

    # If user quantities and figure shape is provided, then use those instead.
    # However, if only quantaties are provided, then just use the default values
    # until I implement smarter subplotting
    # TODO: improve subplotting
    if user_ions and user_shape:
        ions = user_ions
        shape = user_shape
    elif user_ions and user_shape is None:
        print("plot_wind: as only user_ions has been provided with no shape"
              " the default vars and shape are being used instead as"
              " dynamic subplotting hasn't been implemented yet")
        ions = default_ions
        shape = default_shape

    print("\nPlotting the following wind ions...")
    print("\n\t- {}\n".format(ions))

    var_idx = 0
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(12, 16), squeeze=False)
    for i in range(shape[0]):
        for j in range(shape[1]):
            var = ions[var_idx]
            x, z, ion = py_util.get_ion_data(path, root, var, VERBOSE)
            if VERBOSE:
                print("np.mean({}) = {}".format(var, np.mean(ion)))
            with np.errstate(divide="ignore"):
                im = ax[i, j].pcolor(x, z, np.log10(ion), vmin=-5, vmax=0)
            fig.colorbar(im, ax=ax[i, j])
            ax[i, j].set_xlabel("x")
            ax[i, j].set_ylabel("z")
            ax[i, j].set_title(r"$\log_{10}$("+var+")")
            if loglog:
                ax[i, j].set_xscale("log")
                ax[i, j].set_yscale("log")
                ax[i, j].set_xlim(x[1, 1], x[-1, -1])
                ax[i, j].set_ylim(z[1, 1], z[-1, -1])
            var_idx += 1

    fig.tight_layout()
    plt.savefig("{}_ions_plot.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    print("Wind ions plotted...")

    return


def plot_python_complete(spec_file, outname, user_vars=None, user_shape=None, loglog=True):
    """
    Create 2D plots of wind parameters for a Python simulation. Note that your
    own quantities and shape can be provided instead of the default ones. Note
    that both must be provided otherwise the default ones will  be used.

    NOTE: this uses windsave2table which must have been compiled from the same
    git commit otherwise odd things can happen
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot wind quantities for one simulation at a time")
        return
    spec_file = spec_file[0]

    root, path = py_util.get_root_name_and_path(spec_file)
    wind = py_util.get_master_data(path, root, VERBOSE)

    # Default wind quantities to plot and the shape for the resulting figure
    vars = default_vars = ["t_e", "t_r", "ne", "v_x", "v_y", "v_z", "ip", "c4"]
    shape = default_shape = (4, 2)  # 2 rows, 2 columns

    # If user quantities and figure shape is provided, then use those instead.
    # However, if only quantaties are provided, then just use the default values
    # until I implement smarter subplotting
    # TODO: improve subplotting
    if user_vars and user_shape:
        vars = user_vars
        shape = user_shape
    elif user_vars and user_shape is None:
        print("plot_wind: as only user_vars has been provided with no shape"
              " the default vars and shape are being used instead as"
              " dynamic subplotting hasn't been implemented yet")
        vars = default_vars
        shape = default_shape

    print("\nPlotting the following wind quantities...")
    print("\n\t- {}\n".format(vars))

    var_idx = 0
    fig, ax = plt.subplots(shape[0], shape[1], figsize=(12, 16), squeeze=False)
    for i in range(shape[0]):
        for j in range(shape[1]):
            var = vars[var_idx]
            x, z, qoi = py_util.get_wind_quantity(wind, var, VERBOSE)
            if VERBOSE:
                print("np.mean({}) = {}".format(var, np.mean(qoi)))
            with np.errstate(divide="ignore"):
                im = ax[i, j].pcolor(x, z, np.log10(qoi))
            fig.colorbar(im, ax=ax[i, j])
            ax[i, j].set_xlabel("x")
            ax[i, j].set_ylabel("z")
            ax[i, j].set_title(r"$\log_{10}$("+var+")")
            if loglog:
                ax[i, j].set_xscale("log")
                ax[i, j].set_yscale("log")
                ax[i, j].set_xlim(x[1, 1], x[-1, -1])
                ax[i, j].set_ylim(z[1, 1], z[-1, -1])
            var_idx += 1

    fig.tight_layout()
    plt.savefig("{}_wind_plot.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    print("Wind quantities plotted...")

    return


def plot_spec_comps(spec_file, outname):
    """
    Plot the different components of the spectra
    """

    if len(spec_file) > 1:
        if VERBOSE:
            print("Can only plot spectrum components for one spectrum at a time")
        return

    print("\nPlotting spectrum components...")

    spec_file = spec_file[0]
    headers_top = ["Created", "Emitted"]
    headers_bottom = ["CenSrc", "Disk", "Wind", "HitSurf", "Scattered"]

    fig, ax = plt.subplots(1, 2, figsize=(20, 8))
    spec = py_util.read_spec_file(spec_file, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    # Plot created and emitted emission
    for i in range(len(headers_top)):
        flux = np.array(spec[1:, spec[0, :] == headers_top[i]], dtype=float)
        smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
        if SPECLOGLOG:
            ax[0].loglog(wavelength, smoothflux, label=headers_top[i])
        else:
            ax[0].plot(wavelength, smoothflux, label=headers_top[i])
    ax[0].set_xlim(wavelength.min(), wavelength.max())
    ax[0].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[0].legend()

    # Plot CenSrc, Disk, Wind, HitSurf and Scattered emission
    for i in range(len(headers_bottom)):
        flux = np.array(spec[1:, spec[0, :] == headers_bottom[i]], dtype=float)
        smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
        if SPECLOGLOG:
            ax[1].loglog(wavelength, smoothflux, label=headers_bottom[i])
        else:
            ax[1].plot(wavelength, smoothflux, label=headers_bottom[i])
    ax[1].set_xlim(wavelength.min(), wavelength.max())
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
    ax[1].legend()

    plt.savefig("{}_spec_comps.png".format(outname))

    if SHOW_PLOT:
        plt.show()
    else:
        plt.close()

    print("\nSpectrum components plotted...")

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

    print("\nSpectra will be created for the following inclinations:\n")
    print("\t- {}".format(spec_angles))
    print("")

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
            ax.semilogy(blag_spec[:, 0] / (Z + 1), blag_spec[:, 1], label="iPTF15af: Blagorodnova et al. (2019)")

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
                    print("Error: {}° not found in .spec file {}".format(angle, file))
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

            if SPECLOGLOG:
                ax.loglog(wavelength, flux_dist, label=legend)
            else:
                ax.semilogy(wavelength, flux_dist, label=legend)

            ax.set_xlim(WMIN, WMAX)
            yupper, ylower = get_ylims(wavelength, flux_dist)
            ax.set_ylim(ylower, yupper)

            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
            ax.tick_params(labelsize=17)

        ax.legend(loc="best")
        ax.set_title("Viewing angle = {}".format(angle) + "$^{\circ}$", fontsize=20)

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
        print("\n--------------------------")
        exit(1)

    print("Creating plots for the following simulations:\n")
    for i in range(len(spec_files)):
        print("\t- {}".format(spec_files[i]))

    # Angles aren't actually supplied via the command line anymore - but I couldn't
    # be bothered to refactor this bit of the code so I have left it like this
    outname, spec_angles = get_script_arguments(spec_files)
    if len(spec_angles) == 0:
        print("No angles provided")
        exit(2)

    plot_spectra(spec_files, spec_angles, outname)
    plot_spec_comps(spec_files, outname)
    plot_python_complete(spec_files, outname)
    plot_python_ions(spec_files, outname)

    print("\n--------------------------")

    return


if __name__ == "__main__":
    main()