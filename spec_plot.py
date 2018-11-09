#!/usr/bin/env python3

import argparse
import numpy as np
from matplotlib import pyplot as plt
from subprocess import Popen, PIPE

verbose = False      # More info will be printed to screen if True
show_plot = False    # If True, the plot will be drawn on screen
wmin = 850           # The smallest wavelength to show on the spectra
wmax = 2000          # The largest wavelength to show on the spectra
filetype = "png"    # The file type of the output spectra


def read_file(filename, delim=" "):
    """
    Read in data from an external file, line by line whilst ignoring comments.
        - Comments begin with #
        - The default delimiter is assumed to be a space

    Parameters
    ----------
    filename: str
        The path to the file to be read
    delim: str
        The delimiter between values in the file. By default, space is assumed

    Returns
    -------
    lines: ncols x nlines array of strings
        The file as a numpy array of strings for each column and line
    """

    # Try to open the file, otherwise return an error
    try:
        f = open(filename, "r")
        flines = f.readlines()
        f.close()
    except IOError:
        print("Can't open file {}".format(filename))
        return -1

    # Now read in the line one by one and append to the list, lines
    lines = []
    for i in range(len(flines)):
        line = flines[i].strip()
        if delim == " ":
            line = line.split()
        else:
            line = line.split(delim)
        if len(line) > 0:
            if line[0] == "Freq.":  # Clean up the inclination angle names; makes life easier later
                for j in range(len(line)):
                    if line[j][0] == "A":
                        line[j] = line[j].replace("P0.50", "").replace("A", "")            
            if line[0][0] != "#":  # Don't add lines which are comments
                lines.append(line)

    return np.array(lines)


def plot_spectra(files, angles, filename):
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

    assert(len(files) >= 1)
    assert(len(angles) >= 1)

    # Loop over the viewing angles
    for angle in angles:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        # Loop over each .spec file
        for file in files:
            # Read in the data, this could probably be hardcoded instead...
            # I don't think the .spec standard is changing anytime soon.
            data = read_file(file)
            if data == -1:
                return -1
            wavelength = np.array(data[1:, data[0, :] == "Lambda"], dtype=float)
            flux = np.array(data[1:, data[0, :] == "{}".format(angle)], dtype=float)

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

            # Plot the spectra
            ax.plot(wavelength, flux, label=legend)
            ax.set_xlim(wmin, wmax)
            ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
            ax.set_ylabel("Flux", fontsize=17)
            ax.tick_params(labelsize=17)
            ax.legend(loc="best")

        title = "Viewing angle = {}".format(angle) + "$^{\circ}$"
        ax.set_title(title, fontsize=20)
        plt.savefig("{}_{}.{}".format(filename, angle, filetype))

        if show_plot:
            plt.show()
        else:
            plt.close()

    return


def parse_args():
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
        angles = [20, 40, 60, 70, 75, 80, 85, 89]
    elif args.angles.isdigit() is True:
        angles.append(int(args.angles))
    else:
        ang = args.angles.replace(",", " ").split()
        for i in range(len(ang)):
            angles.append(int(ang[i]))

    # Check that the provided viewing angles are legal
    allowed_angles = [20, 40, 60, 70, 75, 80, 85, 89]
    for i in range(len(angles)):
        allowed = False
        for j in range(len(allowed_angles)):
            if angles[i] == allowed_angles[j]:
                allowed = True
                break
        if allowed is False:
            print("Incorrect viewing angle: {}".format(angles[i]))
            return -1, -1

    return args.output_name, angles


def find_spec_files():
    """
    Use the unix find command to find spec files in the current working directory and in directories below

    Parameters
    ----------
    None

    Returns
    -------
    spec_files: nfiles list of str
        The file paths for the .spec files.
    """

    find = "find . -name '*.spec'"
    stdout, stderr = Popen(find, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    spec_files = sorted(stdout.decode("utf-8").split(), key=str.lower)
    if len(spec_files) == 0:
        return -1

    return spec_files


def main():
    """
    The main steering function of the script.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    outname, angles = parse_args()
    if outname == -1 or angles == -1:
        return

    spec_files = find_spec_files()
    if spec_files == -1:
        print("No spec files found, exiting.")
        return

    print("Creating spectra for the following .spec files:")
    for i in range(len(spec_files)):
        print("\t-{}".format(spec_files[i]))
    print("Spectra will be created for the viewing angles:")
    for i in range(len(angles)):
        print("\t-{}°".format(angles[i]))
    err = plot_spectra(spec_files, angles, outname)
    if err == -1:
        print("Error when loading .spec file.")

    return


if __name__ == "__main__":
    main()
