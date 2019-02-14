#!/usr/bin/env python3

"""
Various common routines which are used in scripts concerned with MCRT Python.
"""


import sys
import numpy as np
from sys import exit
from subprocess import Popen, PIPE
from scipy.signal import convolve, boxcar


def tests():
    """
    Warns the user that this script is a utility script and not meant to be run.
    TODO: add unit tests
    """

    print("This script is not designed to be run. Instead, import it using "
          "import py_util.")
    return


def read_spec_file(filename, delim=" "):
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

    try:
        with open(filename, "r") as f:
            flines = f.readlines()
    except IOError:
        print("Can't open file {}".format(filename))
        exit(1)

    # Now read in the line one by one and append to the list, lines
    lines = []
    for i in range(len(flines)):
        line = flines[i].strip()
        if delim == " ":
            line = line.split()
        else:
            line = line.split(delim)
        if len(line) > 0:
            if line[0] == "Freq.":  # Clean up the inclination angle names
                for j in range(len(line)):
                    if line[j][0] == "A":
                        # TODO: this will fail for other phase angles
                        line[j] = line[j].replace("P0.50", "").replace("A", "")
            if line[0][0] != "#":  # Don't add lines which are comments
                lines.append(line)

    return np.array(lines)


def find_spec_files():
    """
    Use the unix find command to find spec files in the current working
    directory and in directories below

    Parameters
    ----------
    None

    Returns
    -------
    spec_files: nfiles list of str
        The file paths for the .spec files.
    """

    find = "find . -name '*.spec'"
    stdout, stderr = Popen(find, stdout=PIPE, stderr=PIPE, shell=True)\
        .communicate()
    spec_files = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")
    spec_files = sorted(spec_files, key=str.lower)

    if err:
        print("ERROR: py_util.find_spec_files: find returning something to "
              "stderr...")
        print("Captured from stderr:")
        print(err)
        return []

    if len(spec_files) == 0:
        print("py_util.find_spec_files: No .spec files found")
        exit(2)

    return spec_files


def find_pf(ignore_out_pf):
    """
    Find parameter files recursively
    """

    command = "find . -type f -name '*.pf'"
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    pfs = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")

    if err:
        print("ERROR: py_util.find_pf: find returned something to stderr...")
        print("Captured from stderr:")
        print(err)
        sys.exit(1)

    if ignore_out_pf:
        for i, dir in enumerate(pfs):
            if dir.find(".out.pf") != -1:
                del(pfs[i])

    if len(pfs) == 0:
        print("py_util.find_pf: no Python parameter files found\n")
        print("--------------------------")
        sys.exit(0)

    return pfs


def get_spec_viewing_angles(specfiles, delim=" "):
    """
    Get all of the unique viewing angles for a set of .spec files.

    TODO: in the future, it may be beneficial to deal with phase angle

    Parameters
    ----------
    specfiles: list of str
        The file path to various .spec files
    delim: str
        The delimiter in the .spec files

    Returns
    -------
    vangles: list of ints
        All of the unique viewing angles found in the provided .spec files
    """

    vangles = []
    # Find the viewing angles in each .spec file
    for i in range(len(specfiles)):
        spec_data = read_spec_file(specfiles[i], delim)
        col_names = spec_data[0, :]
        # Go over the columns and look for viewing angles
        for i in range(len(col_names)):
            if col_names[i].isdigit() is True:
                angle = int(col_names[i])
                duplicate_flag = False
                for va in vangles:  # Check for duplicate angle
                    if angle == va:
                        duplicate_flag = True
                if duplicate_flag is False:
                    vangles.append(angle)

    return np.array(vangles)


def check_viewing_angle(angle, spec):
    """
    Check that a viewing angle is legal

    Parameters
    ----------
    angle: int
        The viewing angle to check
    spec: ncols by nrow array of str
        The raw spectrum read in from file, hence why the data type if still str

    Returns
    -------
    allowed: bool
        If True, angle is a legal angle. Otherwise will return False to indicate
        illegal angle.
    """

    allowed = False
    headers = spec[0, :]
    for i in range(len(headers)):
        if headers[i].isdigit() is True:
            if float(angle) == float(headers[i]):
                allowed = True

    return allowed


def get_root_name_and_path(pf_path):
    """
    Split the path name into a directory and root name of the Python simulation

    Parameters
    ----------
    pf_path: str
        The file path to a .pf file.

    Returns
    -------
    root_name: str
        The root name of the simulation.
    sim: str
        The absolute directory containing the provided .pf file.
    """

    dot = 0
    slash = 0
    for l in range(len(pf_path)):
        if pf_path[l] == "/":
            slash = l + 1
        if pf_path[l] == ".":
            dot = l
    root_name = pf_path[slash:dot]
    sim = pf_path[:slash]

    return root_name, sim


def check_convergence(wd, root):
    """
    Determine the convergence of a Python simulation

    Parameters
    ----------
    wd: str
        The directory of the Python simulation
    root: str
        The root name of the Python simulation

    Returns
    -------
    converge_fraction: float
        The fraction of cells which are converged, between 0 and 1. A value of 0
        means that no cells have converged, whereas a value of 1 means all cells
        have converged.
    """

    diag_path = "{}/diag_{}/{}_0.diag".format(wd, root, root)

    try:
        with open(diag_path, "r") as file:
            diag = file.readlines()
    except IOError:
        print("ERROR: py_util.read_convergence: Couldn't open ready only copy "
              "of {}. Does the diag file exist?".format(diag_path))
        return -1

    converge_fraction = None
    for line in diag:
        if line.find("!!Check_converging") != -1:
            c_string = line.split()[2].replace("(", "").replace(")","")
            converge_fraction = float(c_string)

    if converge_fraction is None:
        print("ERROR: py_util.read_convergence: unable to parse convergence "
              "fraction from diag file {}".format(diag_path))
        return -1

    if 0 > converge_fraction > 1:
        print("ERROR: py_util.read_convergence: the convergence in the "
              "simulation is negative or more than one")
        print("ERROR: py_util.read_convergence: convergence_fraction = {}"
              .format(converge_fraction))
        return -1

    return converge_fraction


def smooth_spectra(flux, smooth, verbose=False):
    """
    Smooth the data flux using a Boxcar averaging algorithm.

    Parameters
    ----------
    flux: list or array of floats/ints
        The data which is to be smoothed
    smooth: int
        Number of points to smooth by
    verbose: bool
        Enable verbose printing

    Returns
    -------
    smooth_flux: numpy array of floats
        The smoothed data.
    """

    if verbose:
        print("pu_util.smooth_spectra: type(flux) = {}".format(type(flux)))

    if type(flux) is not list and type(flux) is not np.ndarray:
        print("ERROR: py_util.smooth_spectra: data to be smoothed is not a "
              "Python list or numpy ndarray")
        print("ERROR: py_util.smooth_spectra: type(flux) = {}"
              .format(type(flux)))
        return flux

    if type(flux) is list:
        if verbose:
            print("py_util.smooth_spectra: Converting Python list to numpy "
                  "ndarray")
        flux = np.array(flux, dtype=float)

    if len(flux.shape) > 2:
        print("ERROR: py_util.smooth_spectra: The data to be smoothed should"
              " be one dimensional and not of dimension {}".format(flux.shape))
        return flux

    if type(smooth) is not int:
        try:
            if verbose:
                print("py_util.smooth_spectra: Converting smooth to an int")
            smooth = int(smooth)
        except ValueError:
            print("ERROR: py_util.smooth_spectra: could not convert smooth {} "
                  "into an integer".format(smooth))
            return flux

    flux = np.reshape(flux, (len(flux), ))
    smooth_flux = convolve(flux, boxcar(smooth) / float(smooth), mode="same")

    return smooth_flux


if __name__ == "__main__":
    tests()
