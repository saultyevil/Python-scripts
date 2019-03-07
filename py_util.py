#!/usr/bin/env python3

"""
Various common routines which are used in scripts concerned with MCRT Python.
Most usually used when plotting output from Python.
"""


import os
import numpy as np
import pandas as pd
from sys import exit
from shutil import which
from socket import gethostname
from subprocess import Popen, PIPE
from scipy.signal import convolve, boxcar


def tests():
    """
    Warns the user that this script is a utility script and not meant to be run.
    TODO: add unit tests
    """

    print("This script is not designed to be run. Instead, import it using "
          "import py_util.")
    get_python_version(verbose=True)

    return


def get_python_version(py="py", verbose=False):
    """
    Get the Python version and commit hash for the provided Python binary. This
    should also work with windsave2table.

    Parameters
    ----------
    py: str
        The name of the Python executable in $PATH whose version will be
        queried
    verbose: bool
        If True, enable more verbose output

    Returns:
    --------
    version: str
        The version of the Python executable
    commit_hash: str
        The commit hash of the Python executable
    """

    version = ""
    commit_hash = ""

    path = which(py)
    if not path:
        print("ERROR: py_util.get_python_version: {} is not in $PATH".format(py))
        exit(1)

    command = "{} --version".format(py)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")

    if err:
        print("py_util.get_python_version: captured from stderr")
        print(stderr)

    for i in range(len(out)):
        if out[i] == "Version":
            version = out[i+1]
        if out[i] == "hash":
            commit_hash = out[i+1]

    if version == "":
        print("ERROR: py_util.get_python_version: couldn't find Python version "
              "for {}".format(py))
        exit(1)
    if commit_hash == "":
        print("ERROR: py_util.get_python_version: couldn't find commit hash for"
              " {}".format(py))
        exit(1)

    if verbose:
        print("Python Version      {}".format(version))
        print("Git commit hash     {}".format(commit_hash))
        print("Short commit hash   {}".format(commit_hash[:7]))

    return version, commit_hash


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

    # TODO: add some input error checking

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
    stdout, stderr = Popen(find, stdout=PIPE, stderr=PIPE, shell=True).communicate()
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
        exit(1)

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
        exit(1)

    if ignore_out_pf:
        for i, dir in enumerate(pfs):
            if dir.find(".out.pf") != -1:
                del(pfs[i])

    if len(pfs) == 0:
        print("py_util.find_pf: no Python parameter files found\n")
        print("--------------------------")
        exit(0)

    pfs = sorted(pfs, )

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

    # TODO: add some input error checking

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

    # TODO: add some input error checking

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

    TODO: the loop could probably be replaced by a str.find() instead

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

    if type(pf_path) != str:
        print("ERROR: py_util.get_root_name_and_path: Provided a {} when "
              "expecting a string".format(type(pf_path)))
        exit(1)

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
        print("ERROR: py_util.read_convergence: Couldn't open read only copy "
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


def get_blagordnova_spec(smooth, verbose=False):
    """
    Load the Blagorodnova iPTF15af UV spectrum into

    Parameters
    ----------
    smooth: int
        The amount of smoothing to be used by the boxcar smoother
    verbose: bool
        Enable verbose printing

    Returns
    -------
    blagorodnova_spec: numpy array of floats
        The UV spectrum of iPTF15af from N. Blagorodnova et al. 2018.
        ArXiv e-prints, arXiv:1809.07446v1 [astro-ph.HE].
            Column 0: Wavelength in Angstroms
            Column 1: Flux per unit wavelength in erg s^-1 cm^-2 A^-1
    """

    try:
        smooth = int(smooth)
    except:
        print("ERROR: py_util.get_blagordnvoa_spec: Unable to convert smooth "
              "into an integer")
        exit(1)

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
    elif hostname == "REX":
        blag_dir = "/home/saultyevil/PySims/TDE/Blagorodnova_iPTF15af.dat"
    else:
        print("Unknown hostname, update py_util with directory for the "
              "Blagordnova spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("Blagordnova spectra being read in from {}".format(blag_dir))

    try:
        blagorodnova_spec = np.loadtxt(blag_dir)
    except IOError:
        print("ERROR: py_util.get_blagordnova_spec: Unable to open the "
              "Blagordnova spectrum from the following path {}".format(blag_dir))
        print("ERROR: py_util.get_blagordnova_spec: check the directories "
              "provided in the script")
        exit(1)

    blagorodnova_spec[:, 1] = smooth_spectra(blagorodnova_spec[:, 1], smooth)

    return blagorodnova_spec


def run_windsave2table(path, root, verbose=False):
    """
    Use windsave2table to get the details of the wind

    Parameters
    ----------
    root: str
        The rootname of the Python simulation.

    verbose: bool
        If True, enable verbose output

    Returns
    -------

    """

    # Look for a version file in the directory to see if windsave2table is the
    # on the correct git hash to read the wind_save file
    try:
        with open("version", "r") as f:
            lines = f.readlines()
        run_version = lines[0]
        run_hash = lines[1]
        wind_version, wind_hash = get_python_version("windsave2table", verbose)
        if verbose:
            print("wind_save: version {} hash {}".format(run_version, run_hash))
            print("windsave2table: version {} hash {}".format(wind_version, wind_hash))
        if run_version != wind_version and run_hash != wind_hash:
            print("py_util.run_windsave2table: windsave2table git hash "
                  "and the git hash of the wind_save file are different")
            return True
    except IOError:
        if verbose:
            print("py_util.run_windsave2table: no version file assuming "
                  "everything will be a-ok")

    check = which("windsave2table")
    if not check:
        print("ERROR: py_util.get_wind_details: windsave2table not in $PATH")
        return False

    command = "cd {}; windsave2table {}".format(path, root)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if verbose:
        print(output)

    if err:
        print("py_util.get_wind_details: the following was sent to stderr:")
        print(err)

    return True


def get_ion_data(path, root, ion, verbose=False):
    """
    Get data for a certain ion from a Python simulation using windsave2table
    """

    if type(root) is not str:
        print("ERROR: py_util.get_ion_data: the root name provided is not a"
              " string :-(")
        return

    ele_idx = ion.find("_")
    element = ion[:ele_idx]
    ion_level = ion[ele_idx+1:]
    ion_file = "{}.0.{}.txt".format(root, element)

    # Add check to avoid situations where we re-create data with the wrong commit
    file_exists = os.path.isfile(ion_file)

    if not file_exists:
        if verbose:
            print("py_util.get_ion_data: running windsave2table")
        worked = run_windsave2table(path, root, verbose)
        if not worked:
            print("ERROR: py_util.get_ion_data: could not use windsave2table")
            return
    elif verbose:
        print("py_util.get_ion_data: required files already exist hence"
              " windsave2table will not be run")

    try:
        ions = pd.read_table(ion_file, delim_whitespace=True)
    except IOError:
        print("ERROR: py_util.get_ion_data: Could not find file for element {}"
              .format(element))
        return

    # Figure out the number of cells in the x and z directions of the simulation
    xi = ions["i"]
    zj = ions["j"]
    nx_cells = int(np.max(xi) + 1)
    nz_cells = int(np.max(zj) + 1)

    try:
        ion_wanted = ions[ion_level].values.reshape(nx_cells, nz_cells)
    except KeyError:
        print("ERROR: py_util.get_ion_data: Could not find ion {} for element "
              "{}".format(ion_level, element))
        return

    x = ions["x"].values.reshape(nx_cells, nz_cells)
    z = ions["z"].values.reshape(nx_cells, nz_cells)
    inwind = ions["inwind"].values.reshape(nx_cells, nz_cells)

    # Create a masked array where only the quantities in the wind are returned
    mask = (inwind < 0)
    ion_mask = np.ma.masked_where(mask, ion_wanted)

    return x, z, ion_mask


def get_master_data(path, root, verbose=False):
    """
    Get the "important" information for windsave2table

    Reads in root.0.master.txt and root.0.heat.txt and puts them into one
    numpy array
    """

    if type(root) is not str:
        print("ERROR: py_util.get_master_data: the root name provided is not a"
              " string :-(")
        return

    heat_file = "{}.0.heat.txt".format(root)
    master_file = "{}.0.master.txt".format(root)

    # Add check to avoid situations where we re-create data with the wrong commit
    heat_exists = os.path.isfile(heat_file)
    master_exists = os.path.isfile(master_file)

    if not heat_exists and not master_exists:
        if verbose:
            print("py_util.get_master_data: running windsave2table")
        worked = run_windsave2table(path, root, verbose)
        if not worked:
            print("ERROR: py_util.get_master_data: could not use windsave2table")
            return
    elif verbose:
        print("py_util.get_master_data: required files already exist hence"
              " windsave2table will not be run")


    try:
        heat = pd.read_table(heat_file, delim_whitespace=True)
        master = pd.read_table(master_file, delim_whitespace=True)
        # master = master.drop(columns=["xcen", "zcen"])
    except IOError:
        print("ERROR: py_uilt.get_master_data: Could not find master or heat "
              "file from windsave2table")
        return

    # Completely hideous way of appending columns, but join, merge and concat
    # were doing bad things and this worked so whatever
    append = heat.columns.values[14:]
    for i, col in enumerate(append):
        master[col] = pd.Series(heat[col])

    if verbose:
        print("Headers for data read in:")
        print(master.columns.values)
        print("")

    master.to_csv("{}.ep.complete".format(root))

    return master


def get_wind_quantity(wind, quantity, verbose):
    """
    Get and return x, z coordinates and a quantity of interest which is a wind
    quantity from a Python simulation. This assumes that wind is a Pandas
    data frame gotten from something like py_util.get_master_data or
    py_util.get_ion_data

    Parameters
    ----------
    wind: pandas dataframe
        An pandas dataframe of the wind quantities from windsave2table
    quantity: str
        The quantity wanting to be masked in question
    verbose: bool
        If True, enable verbose logging

    Returns
    -------


    """

    if type(quantity) is not str:
        print("ERROR: py_util.get_wind_quantity: the quantity desired needs to"
              " be given as a string")

    # Figure out the number of cells in the x and z directions of the simulation
    xi = wind["i"]
    zj = wind["j"]
    nx_cells = int(np.max(xi) + 1)
    nz_cells = int(np.max(zj) + 1)

    # Read in the quantity of interest, qoi, and x and z values
    try:
        qoi = wind[quantity].values.reshape(nx_cells, nz_cells)
    except KeyError:
        print("Quantity {} not found in table, try again!".format(quantity))
        return

    x = wind["x"].values.reshape(nx_cells, nz_cells)
    z = wind["z"].values.reshape(nx_cells, nz_cells)
    inwind = wind["inwind"].values.reshape(nx_cells, nz_cells)

    # Create a masked array where only the quantities in the wind are returned
    mask = (inwind < 0)
    qoi_mask = np.ma.masked_where(mask, qoi)

    return x, z, qoi_mask


if __name__ == "__main__":
    tests()
