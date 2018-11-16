#!/usr/bin/env python3

"""
Create spectra from the .spec files from a MCRT Python simulation.

The script requires you to provide 2 arguments:
    - The base output name of the spectra which are going to be produced
    - The viewing angles to create spectra for: possible choices for this include all, a single number or a list of
      numbers separated by a comma, i.e. 20,30,40
There are also some optional arguments:
    - wmin: the smallest wavelength to plot
    - wmax: the largest wavelength to plot
    - filetype: the file type of the output spectra plots, by default this is png
There are two optional switches also:
    - -v or --verbose which will enable more verbose output
    - -s or --show which will show the output spectra as well as saving to disk

This script works by using the unix find command to find .spec files recursively in the current directory and deeper.
Thus, it's possible to plot spectra for a single simulation or by plotting multiple simulations at once on the same
plot for comparison purposes.

Example usage:
    python spec_plot.py qDNe all -v -show
"""

import py_util
import argparse
import numpy as np
from matplotlib import pyplot as plt


verbose = False      # More info will be printed to screen if True
show_plot = False    # If True, the plot will be drawn on screen
wmin = 850           # The smallest wavelength to show on the spectra
wmax = 2000          # The largest wavelength to show on the spectra
filetype = "png"    # The file type of the output spectra


def get_outname_angles(specfiles):
    """
    Parse the various global parameters from the command line.

    Parameters
    ----------
    None

    Returns
    -------
    outputname: str
        The base filename for the output spectra.
    angles: nangles list of ints.
        The viewing angles of the spectra to plot.
    """

    p = argparse.ArgumentParser(description="")
    p.add_argument("output_name", type=str, help="The base name for the output spectra")
    p.add_argument("angles", type=str,
                   help="The viewing angles to plot: all, a single angle or a comma separated list")
    p.add_argument("-wmin", type=float, nargs="?", action="store", help="The smallest wavelength to show")
    p.add_argument("-wmax", type=float, nargs="?", action="store", help="The largest wavelength to show")
    p.add_argument("-filetype", type=str, nargs="?", action="store",
                   help="The file format of the output mspectra")
    p.add_argument("-v", "--verbose", help="Increase output to screen", action="store_true")
    p.add_argument("-s", "--show", help="Show the plots on screen", action="store_true")
    args = p.parse_args()

    # Assign the optional arguments to their global vars
    if args.verbose:
        global verbose
        verbose = True
    if args.show:
        global show_plot
        show_plot = True
    if args.wmin:
        global wmin
        wmin = args.wmin
    if args.wmax:
        global wmax
        wmax = args.wmax
    if args.filetype:
        global filetype
        filetype = args.filetype

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


def plot_spectra():
    """
    Plot the spectra for the give .spec files and for the given viewing angles.

    Parameters
    ----------
    files: nfiles list of strings.
        The file pathes to the .spec files.
    angles: nangles list of ints.
        The viewing angles to be plotted.
    filename: str
        The base filename for the output spectra.

    Returns
    -------
    None
    """

    # Get the output name, the viewing angles and the file paths to the .spec files
    specfiles = py_util.find_spec_files()
    if len(specfiles) == 0:
        print("No spec files found, exiting.")
        return
    outname, angles = get_outname_angles(specfiles)
    if len(angles) == 0:
        print("No angles were provided")
        return

    print("Creating spectra for the following .spec files:")
    for i in range(len(specfiles)):
        print("\t+ {}".format(specfiles[i]))
    print("Spectra will be created for the following viewing angles:")
    for i in range(len(angles)):
        print("\t+ {}".format(angles[i]))

    # Loop over the viewing angles
    for angle in angles:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        # Loop over each .spec file
        for file in specfiles:
            # Read in the data, this could probably be hardcoded instead...
            # I don't think the .spec standard is changing anytime soon.
            spec = py_util.read_file(file)
            allowed = py_util.check_viewing_angle(angle, spec)
            if allowed is False:
                print("Error: viewing angle {} not found in .spec file {}".format(angle, file))
                continue
            wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
            flux = np.array(spec[1:, spec[0, :] == "{}".format(angle)], dtype=float)

            # Find the final slash and final dot and assume between this slash and
            # dot is the rootname of the Python pf
            slashIdx = 0
            dotIdx = len(file) - 1
            for i in range(len(file)):
                if file[i] == "/":
                    slashIdx = i
                elif file[i] == ".":
                    dotIdx = i
            rootname = file[slashIdx+1:dotIdx]
            legend = rootname.replace("_", " ")

            # Plot the spectrum
            ax.plot(wavelength, flux, label=legend)
            ax.set_xlim(wmin, wmax)
            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel("Flux", fontsize=17)
            ax.tick_params(labelsize=17)
            ax.legend(loc="best")

        title = "Viewing angle = {}".format(angle) + "$^{\circ}$"
        ax.set_title(title, fontsize=20)
        plt.savefig("{}_{}.{}".format(outname, angle, filetype))

        if show_plot:
            plt.show()
        else:
            plt.close()

    return


if __name__ == "__main__":
    plot_spectra()
