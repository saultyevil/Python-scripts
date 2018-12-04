#!/usr/bin/env python3

"""
Various common routines which are used in scripts concerned with MCRT Python.
"""

import numpy as np
from sys import exit
from subprocess import Popen, PIPE


def tests():
    """
    Warns the user that this script is a utility script and not meant to be run.
    TODO: add unit tests

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    print("This script is not designed to be run. Instead, import it using import py_util.")
    exit(0)


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
            if line[0] == "Freq.":  # Clean up the inclination angle names; makes life easier later
                for j in range(len(line)):
                    if line[j][0] == "A":
                        line[j] = line[j].replace("P0.50", "").replace("A", "")
            if line[0][0] != "#":  # Don't add lines which are comments
                lines.append(line)

    return np.array(lines)


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
    specfiles = stdout.decode("utf-8").split()
    specfiles = sorted(specfiles, key=str.lower)
    if len(specfiles) == 0:
        print("No .spec files found")
        exit(2)

    return specfiles


def get_spec_viewing_angles(specfiles, delim=" "):
    """
    Get all of the unique viewing angles for a set of .spec files.

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
        specdata = read_file(specfiles[i], delim)
        colnames = specdata[0, :]
        # Go over the columns and look for viewing angles
        for i in range(len(colnames)):
            if colnames[i].isdigit() is True:
                ang = int(colnames[i])
                dup = False
                for va in vangles:  # Check for duplicate angle
                    if ang == va:
                        dup = True
                if dup is False:
                    vangles.append(ang)

    return vangles


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
        If True, angle is a legal angle. Otherwise will return False to indicate illegal angle.
    """

    headers = spec[0, :]
    allowed = False
    for i in range(len(headers)):
        if headers[i].isdigit() is True:
            if float(angle) == float(headers[i]):
                allowed = True

    return allowed


if __name__ == "__main__":
    tests()
