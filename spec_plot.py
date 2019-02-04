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
    - blag: default plotting parameter choices for modelling iPTF15af (TDE)

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
from socket import gethostname
from matplotlib import pyplot as plt
from scipy.signal import convolve, boxcar


VERBOSE = False                # More info will be printed to screen if True
SHOW_PLOT = False              # If True, the plot will be drawn on screen
TDE_PLOT = False               # Enable default TDE plotting
WMIN = None                    # The smallest wavelength to show on the spectra
WMAX = None                    # The largest wavelength to show on the spectra
FILETYPE = "png"               # The file type of the output spectra
SMOOTH = 11                    # The amount of smoothing for the spectra
OBSERVE_DIST = 100 * 3.086e18  # 100 pc in cm
Z = 0                          # Redshift


def get_outname_and_angles(specfiles):
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
    p.add_argument("-blag", action="store_true", help=
        "Plot for the Blagodovnova (?) UV TDE spec")
    args = p.parse_args()

    global VERBOSE
    global SHOW_PLOT
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
    if args.wmin:
        WMIN = args.wmin
    if args.wmax:
        WMAX = args.wmax
    if args.filetype:
        FILETYPE = args.filetype
    if args.smooth:
        SMOOTH = args.smooth
    if args.dist:
        OBSERVE_DIST = args.dist
    if args.z:
        Z = args.z
    if args.blag:
        TDE_PLOT = True
        OBSERVE_DIST = 1.079987153448e+27
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


def load_blag_spec():
    """
    Load the Blagorodnova iPTF15af UV spectrum into
    """

    blag_dir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        blag_dir = "/home/saultyevil/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "excession":
        blag_dir = "/home/ejp1n17/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "REXBOOK-AIR.local":
        blag_dir = "/Users/saultyevil/Dropbox/DiskWinds/PySims/TDE/" \
                   "Blagorodnova_iPTF15af.dat"
    elif hostname == "REXBUNTU":
        blag_dir = "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/" \
                   "Blagorodnova_iPTF15af.dat"
    else:
        print("Unknown hostname, update script with directory for the "
              "Blagorodnova spectrum")
        exit(1)

    if VERBOSE:
        print("Hostname: {}".format(hostname))
        print("Blagordnova spectra being read in from {}".format(blag_dir))

    blagorodnovaspec = np.loadtxt(blag_dir, skiprows=36)
    sm_blagorodnovaspec = np.copy(blagorodnovaspec)
    sm_blagorodnovaspec[:, 1] = convolve(
        sm_blagorodnovaspec[:, 1], boxcar(SMOOTH) / float(SMOOTH), mode="same")

    return sm_blagorodnovaspec


def print_info (specfiles, angles):
    """
    Print extra information to the screen
    """

    print("Creating spectra for the following .spec files:\n")
    for i in range(len(specfiles)):
        print("\t- {}".format(specfiles[i]))

    print("\nSpectra will be created for the following viewing angles:\n")
    for i in range(len(angles)):
        print("\t- {}".format(angles[i]))

    print("")

    return


def plot_spectra():
    """
    Plot the spectra for the give .spec files and for the given viewing angles.
    """

    print("--------------------------\n")

    # Get the output name, the viewing angles and the file paths to the .spec files
    specfiles = py_util.find_spec_files()
    if len(specfiles) == 0:
        print("No spec files found")
        exit(1)

    outname, angles = get_outname_and_angles(specfiles)
    if len(angles) == 0:
        print("No angles provided")
        exit(1)

    print_info(specfiles, angles)

    print("Beginning plotting...\n")

    # Loop over all the possible viewing angles
    for angle in angles:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))

        if TDE_PLOT:
            blag_spec = load_blag_spec()
            ax.semilogy(blag_spec[:, 0], blag_spec[:, 1],
                        label="Blagorodnova et al. 2018")

        # Loop over each possible .spec file
        for file in specfiles:
            rootname, filepath = py_util.get_root_name_and_path(file)
            legend = filepath + rootname
            print("\t- {}: {}°".format(legend, angle))

            # Read in .spec file and check that it can be plotted for the
            # current viewing angle
            spec = py_util.read_file(file)
            allowed = py_util.check_viewing_angle(angle, spec)
            if not allowed:
                if VERBOSE:
                    print("Error: {}° not found in .spec file {}"
                          .format(angle, file))
                continue

            #
            # Read the wavelength and flux and smooth the flux
            # There is something weird happening to determine which index data
            # should be extracted from. I did this because all of the data in
            # the read in .spec file are strings
            #

            idx = 0
            for i in range(spec.shape[1]):
                if spec[0, i].isdigit() and float(spec[0, i]) == float(angle):
                    # Found the index for the current angle to be plotted :-)
                    spec[0, i] = float(spec[0, i])
                    idx = i
                    break

            flux = np.array(spec[1:, idx], dtype=float)
            flux = np.reshape(flux, (len(flux),))
            smoothflux = convolve(flux, boxcar(SMOOTH) / float(SMOOTH),
                                  mode="same")
            wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

            # Plot and scale flux for observer distance
            default_dist = 100 * 3.08567758128e18  #  default Python distance
            ax.semilogy(wavelength * (Z + 1), smoothflux *
                        (default_dist**2 / OBSERVE_DIST**2), label=legend)
            ax.set_xlim(WMIN, WMAX)
            # TODO: implement smart way of determining correct ylim value
            # ax.set_ylim(1e-18, 1e-13)
            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel("$F_{\lambda}$ (ergs/s/cm$^{2}$/$\AA$)", fontsize=17)
            ax.tick_params(labelsize=17)

        ax.legend(loc="best")
        ax.set_title("Viewing angle = {}".format(angle) + "$^{\circ}$",
                     fontsize=20)
        plt.savefig("{}_{}.{}".format(outname, angle, FILETYPE))

        if SHOW_PLOT:
            plt.show()
        else:
            plt.close()

    print("\nAll done :-)")
    print("\n--------------------------")

    return


if __name__ == "__main__":
    plot_spectra()
