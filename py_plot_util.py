#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Various functions used throughout plotting the output from Python. This script
should be imported into other scripts rather than being itself run.
"""

import os
import numpy as np
import pandas as pd
from sys import exit
from shutil import which
from socket import gethostname
from subprocess import Popen, PIPE
from typing import Tuple, List, Union
from scipy.signal import convolve, boxcar


def tests():
    """
    Warns the user that this script is a utility script and not meant to be run. It will also print the current
    Python version to screen.

    TODO: add unit tests

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    print("This script is not designed to be run. Instead, import it using import py_util.")
    get_python_version(verbose=True)

    return


def get_python_version(py: str = "py", verbose: bool = False) -> Tuple[str, str]:
    """
    Get the Python version and commit hash for the provided Python binary. This should also work with windsave2table.

    Parameters
    ----------
    py              str, optional
                    The name of the Python executable in $PATH whose version will be queried
    verbose         bool, optional
                    Enable verbose logging

    Returns:
    --------
    version         str
                    The version number of Python
    commit_hash     str
                    The commit hash of Python
    """

    version = ""
    commit_hash = ""

    path = which(py)
    if not path:
        print("py_util.get_python_version: {} is not in $PATH".format(py))
        return version, commit_hash

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
            version = out[i + 1]
        if out[i] == "hash":
            commit_hash = out[i + 1]

    if version == "" and verbose:
        print("py_util.get_python_version: couldn't find version for {}".format(py))
    if commit_hash == "" and verbose:
        print("py_util.get_python_version: couldn't find commit hash for {}".format(py))

    if verbose and version and commit_hash:
        print("{} version {}".format(py, version))
        print("Git commit hash    {}".format(commit_hash))
        print("Short commit hash  {}".format(commit_hash[:7]))

    return version, commit_hash


def read_spec_file(file_name: str, delim: str = " ", pandas_table: bool = False) -> Union[np.array, pd.DataFrame]:
    """
    Read in data from an external file, line by line whilst ignoring comments.
        - Comments begin with #
        - The default delimiter is assumed to be a space

    Parameters
    ----------
    file_name       str
                    The directory path to the spec file to be read in
    delim           str, optional
                    The delimiter between values in the file, by default a space is assumed
    pandas_table    bool, optional
                    Return the spectrum as a Pandas DataFrame instead of a Numpy array

    Returns
    -------
    lines           Numpy array of strings or a Pandas DataFrame
                    The .spec file as a Numpy array or a Pandas DataFrame
    """

    try:
        with open(file_name, "r") as f:
            flines = f.readlines()
    except IOError:
        print("Can't open file {}".format(file_name))
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
                        # TODO: this will fail for other phase angles, maybe?
                        line[j] = line[j].replace("P0.50", "").replace("A", "")
            if line[0][0] != "#":  # Don't add lines which are comments
                lines.append(line)

    if pandas_table:
        return pd.DataFrame(lines[1:], columns=lines[0])
    else:
        return np.array(lines)


def find_spec_files() -> List[str]:
    """
    Use the unix find command to find spec files in the current working directory and in directories below

    Parameters
    ----------
    None

    Returns
    -------
    spec_files:     list[str]
                    The file paths of various .spec files
    """

    find = "find . -name '*.spec'"
    stdout, stderr = Popen(find, stdout=PIPE, stderr=PIPE, shell=True).communicate()
    spec_files = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")
    spec_files = sorted(spec_files, key=str.lower)

    if err:
        print("py_util.find_spec_files: find returning something to stderr...")
        print("Captured from stderr:")
        print(err)
        return []

    return spec_files


def find_pf(ignore_out_pf: bool = True):
    """
    Find parameter files recursively from the directory this function is called in.

    Parameters
    ----------
    ignore_out_pf           bool, optional
                            Ignore Python .out.pf files

    Returns
    -------
    pfs                     list
                            The file path for any Python pf files founds
    """

    command = "find . -type f -name '*.pf'"
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")

    if err:
        print("py_util.find_pf: find returned something to stderr...")
        print("Captured from stderr:")
        print(err)
        exit(1)

    pfs = []
    if ignore_out_pf is True:
        for file in out:
            if file.find(".out.pf") == -1 and file is not "py_wind.pf":
                pfs.append(file)

    pfs = sorted(pfs, key=str.lower)

    return pfs


def get_spec_viewing_angles(spec_names: List[str], delim: str = " ") -> np.array:
    """
    Get all of the unique inclination angles for a set of Python .spec files.

    TODO: add better support for binary system phase angles

    Parameters
    ----------
    spec_name       list[str]
                    The directory path to Python .spec files
    delim           str, optional
                    The delimiter in the .spec files, assumed to be spaces by default

    Returns
    -------
    iangles         list[int]
                    All of the unique inclination angles found in the Python .spec files
    """

    # TODO: add some input error checking
    # TODO: why did I write it like this? It's hideous

    iangles = []
    # Find the viewing angles in each .spec file
    for i in range(len(spec_names)):
        spec_data = read_spec_file(spec_names[i], delim)
        col_names = spec_data[0, :]
        # Go over the columns and look for viewing angles
        for i in range(len(col_names)):
            if col_names[i].isdigit() is True:
                angle = int(col_names[i])
                duplicate_flag = False
                for va in iangles:  # Check for duplicate angle
                    if angle == va:
                        duplicate_flag = True
                if duplicate_flag is False:
                    iangles.append(angle)

    return iangles


def check_inclination_angle(inclination: int, spec: np.array) -> bool:
    """
    Check that an inclination angle is in the spectrum array.

    This is a terrible function.

    Parameters
    ----------
    inclination         int
                        The inclination angle to check
    spec                np.arrayp[str]
                        The spectrum array to read -- assume that it is a np.array of strings
                        Note that tde_spec_plot has a similar routine for pd.DataFrame's, whoops!
    Returns
    -------
    allowed             bool
                        If True, angle is a legal angle, otherise fault
    """

    # TODO: add some input error checking

    allowed = False
    headers = spec[0, :]
    for i in range(len(headers)):
        if headers[i].isdigit() is True:
            if float(inclination) == float(headers[i]):
                allowed = True

    return allowed


def get_root_name_and_path(pf_path: str) -> Tuple[str, str]:
    """
    Split a path name into a directory path and root name for a Python simulation.

    TODO: the loop could probably be replaced by a str.find() instead

    Parameters
    ----------
    pf_path         str
                    The directory path to a Python .pf file

    Returns
    -------
    root_name       str
                    The root name of the Python simulation
    path            str
                    The directory path containing the provided Python .pf file
    """

    if type(pf_path) != str:
        print("py_util.get_root_name_and_path: Provided a {} when expecting a string".format(type(pf_path)))
        exit(1)

    dot = 0
    slash = 0
    for l in range(len(pf_path)):
        if pf_path[l] == "/":
            slash = l + 1
        if pf_path[l] == ".":
            dot = l
    root_name = pf_path[slash:dot]
    path = pf_path[:slash]

    return root_name, path


def smooth_flux(flux: np.array, smooth: Union[int, float], verbose: bool = False) -> np.array:
    """
    Smooth 1d data using a boxcar smoother.

    Parameters
    ----------
    flux                np.array[float]
                        The data to smooth using the boxcar filter
    smooth              int
                        The size of the window for the boxcar filter
    verbose             bool
                        Enable verbose logging

    Returns
    -------
    smoothed            np.array[float]
                        The smoothed data
    """

    if type(flux) is not list and type(flux) is not np.ndarray:
        print("py_util.smooth_spectra: data to be smoothed is not a Python list or numpy ndarray")
        print("py_util.smooth_spectra: type(flux) = {}".format(type(flux)))
        return flux

    if type(flux) is list:
        if verbose:
            print("py_util.smooth_spectra: Converting Python list to numpy ndarray")
        flux = np.array(flux, dtype=float)

    if len(flux.shape) > 2:
        print("py_util.smooth_spectra: The data to be smoothed should be one dimensional and not of dimension "
              "{}".format(flux.shape))
        return flux

    if type(smooth) is not int:
        try:
            if verbose:
                print("py_util.smooth_spectra: Converting smooth to an int")
            smooth = int(smooth)
        except ValueError:
            print("py_util.smooth_spectra: could not convert smooth {} into an integer".format(smooth))
            return flux

    flux = np.reshape(flux, (len(flux),))
    smoothed = convolve(flux, boxcar(smooth) / float(smooth), mode="same")

    return smoothed


def get_ylims(wavelength: np.array, flux: np.array, wmin: float, wmax: float, verbose: bool = False) -> Tuple[float, float]:
    """
    Find more appropriate y limits to use when the wavelength range has been limited.

    Parameters
    ----------
    wavelength          np.array[float]
                        An array containing all wavelengths in a spectrum
    flux                np.array[float]
                        An array containing the flux at each wavelength point
    wmin                float
                        The shortest wavelength which is being plotted
    wmax                float
                        The longest wavelength which is being plotted
    verbose             bool, optional
                        Enable verbose logging

    Returns
    -------
    yupper          float
                    The upper y limit to use with the wavelength range
    ylower          float
                    The lower y limit to use with the wavelength range
    """

    if wavelength.shape != flux.shape:
        print("py_plot_util.get_y_lims: wavelength and flux are of different dimensions!")
        exit(1)

    wmin_flux = flux.min()
    wmax_flux = flux.max()

    if wmin:
        wmin_idx = np.abs(wavelength - float(wmin)).argmin()
        if verbose:
            print("wmin_idx: {}".format(wmin_idx))
        wmin_flux = flux[wmin_idx]
    if wmax:
        wmax_idx = np.abs(wavelength - float(wmax)).argmin()
        if verbose:
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

    return yupper, ylower


def get_iPTF15af_spec(smooth: int, verbose: bool = False) -> np.array:
    """
    Return an array containing the UV spectrum for iPTF15af as in Blagordonova et al. (2019)

    Parameters
    ----------
    smooth              int
                        The size of the window for the boxcar filter
    verbose             bool, optional
                        Enable verbose logging

    Returns
    --------
    iPTF15af_spec       np.array[float]
                        the UV spectrum of iPTF15af from N. Blagorodnova et al. (2019)
                            column 0: wavelength in Angstroms
                            column 1: flux per unit wavelength in erg s^-1 cm^-2 A^-1
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.get_iPTF15af_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        spec_dir = "/home/saultyevil/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "excession":
        spec_dir = "/home/ejp1n17/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "REXBOOK-AIR.local":
        spec_dir = "/Users/saultyevil/Dropbox/DiskWinds/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "REXBUNTU":
        spec_dir = "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/Blagorodnova_iPTF15af.dat"
    elif hostname == "REX":
        spec_dir = "/home/saultyevil/PySims/TDE/Blagorodnova_iPTF15af.dat"
    else:
        print("Unknown hostname, update py_util with directory for the Blagordnova spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("iPTF15af spectra being read in from {}".format(spec_dir))

    try:
        iPTF15af_spec = np.loadtxt(spec_dir)
    except IOError:
        print("py_util.get_iPTF15af_spec: Unable to open the iPTF15af UV spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    iPTF15af_spec[:, 1] = smooth_flux(iPTF15af_spec[:, 1], smooth)

    return iPTF15af_spec


def get_ASSASN14li_spec(smooth: int, verbose: bool = False) -> np.array:
    """
    Return an array containing the UV spectrum for ASSASN14li as in Cenko et al. (2016)

    Parameters
    ----------
    smooth              int
                        The size of the window for the boxcar filter
    verbose             bool, optional
                        Enable verbose logging

    Returns
    -------
    ASSASN_14li_spec    np.array[float]
                        the UV spectrum of ASASSN-14li from Cenko et al. 2016.
                            column 0: wavelength in Angstroms
                            column 1: flux per unit wavelength in erg s^-1 cm^-2 A^-1
                            column 2: error in flux
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.get_ASSASN_14li_spec: Unable to convert smooth into an integer")
        exit(1)

    cenk_dir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        cenk_dir = "/home/saultyevil/PySims/TDE/ASASSN-14li_spec_Cenko.dat"
    elif hostname == "excession":
        cenk_dir = "/home/ejp1n17/PySims/TDE/ASASSN-14li_spec_Cenko.dat"
    elif hostname == "REXBOOK-AIR.local":
        cenk_dir = "/Users/saultyevil/Dropbox/DiskWinds/PySims/TDE/ASASSN-14li_spec_Cenko.dat"
    elif hostname == "REXBUNTU":
        cenk_dir = "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/ASASSN-14li_spec_Cenko.dat"
    elif hostname == "REX":
        cenk_dir = "/home/saultyevil/PySims/TDE/ASASSN-14li_spec_Cenko.dat"
    else:
        print("Unknown hostname, update py_util with directory for the Cenko spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("Cenko spectra being read in from {}".format(cenk_dir))

    try:
        ASSASN_14li_spec = np.loadtxt(cenk_dir)
    except IOError:
        print("py_util.get_ASSASN_14li_spec: Unable to open the ASSASSN_14li spectrum from the following path {}. "
              "Update the directories in the script".format(cenk_dir))
        exit(1)

    ASSASN_14li_spec[:, 1] = smooth_flux(ASSASN_14li_spec[:, 1], smooth)

    return ASSASN_14li_spec


def get_common_line_ids() -> dict:
    """
    Return a dictionary containing the major absorption and emission lines which I'm interested in. The wavelengths
    of the lines are in Angstrom.

    Parameters
    ----------
    None

    Returns
    -------
    common_lines        dict
                        A dictionary where the keys are the line names and the values are the wavelength of the lines
                        in Angstroms
    """

    common_lines = {
        "P V": 1118,
        r"L$_{\alpha}$": 1216,
        "N V": 1240,
        "Si IV": 1400,
        "N IV": 1489,
        "C IV": 1549,
        "N III]": 1750,
        "Al III": 1854,
        "C III]": 1908,
        "Mg II": 2798,
        "He II": 4686,
        # "FeII": ,  # can't find a fucking value for this cunt
    }

    return common_lines


def run_windsave2table(path: str, root: str, verbose: bool = False) -> Union[int, None]:
    """
    Run windsave2table in the directory given by path. This function will also create a *.ep.complete file which
    combines both the heat and master data tables together.

    Parameters
    ----------
    path                str
                        The directory of the Python simulation where windsave2table will be run
    root                str
                        The root name of the Python simulation
    verbose             bool, optional
                        Enable verbose logging

    Returns
    -------
    Returns an int on error.
    """

    # Look for a version file in the directory to see if windsave2table is the
    # on the correct git hash to read the wind_save file
    try:
        wind_version, wind_hash = get_python_version("windsave2table", verbose)
        with open("version", "r") as f:
            lines = f.readlines()
        run_version = lines[0]
        run_hash = lines[1]
        if verbose:
            print("wind_save: version {} hash {}".format(run_version, run_hash))
            print("windsave2table: version {} hash {}".format(wind_version, wind_hash))
        if run_version != wind_version and run_hash != wind_hash:
            print("py_util.run_windsave2table: windsave2table version and the wind_save version file are different. "
                  "Results may be garbage!")
    except IOError:
        if verbose:
            print("py_util.run_windsave2table: no version file assuming everything will be a-ok")

    # Check that windsave2table can be used
    in_path = which("windsave2table")
    if not in_path:
        print("py_util.run_windsave2table: windsave2table not in $PATH and executable")
        return 1

    # Run windsave2table
    command = "cd {}; Setup_Py_Dir; windsave2table {}".format(path, root)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if verbose:
        print(output)

    if err:
        print("py_util.get_wind_details: the following was sent to stderr:")
        print(err)

    # Now create a "complete" file which is the master and heat put together into one csv
    heat_file = "{}.0.heat.txt".format(root)
    master_file = "{}.0.master.txt".format(root)

    try:
        heat = pd.read_table(heat_file, delim_whitespace=True)
        master = pd.read_table(master_file, delim_whitespace=True)
    except IOError:
        print("py_util.run_windsave2table: Could not find master or heat file for Python simulation")
        return 1

    # Terrible hack to append the columns I want to the end of the table :-)
    append = heat.columns.values[14:]
    for i, col in enumerate(append):
        master[col] = pd.Series(heat[col])

    master.to_csv("{}/{}.ep.complete".format(path, root), sep=" ")

    return


def get_wind_data(root_name: str, var: str, var_type: str, path: str = "./") -> Tuple[np.array, np.array, np.array]:
    """
    Read in variables contained within a windsave2table file. Requires the user to have already run windsave2table
    already so the data is in the directory.

    Parameters
    ----------
    root_name           str
                        The root name of the Python simulation
    var                 str
                        The name of the quantity from the Python simulation
    var_type            str
                        The type of quantity, this can be wind or ion
    path                str, optional
                        The directory containing the Python simulation

    Returns
    -------
    x                   np.array[float]
                        The x coordinates of the wind in the Python simulation
    z                   np.array[float]
                        The z coordinates of the wind in the Python simulation
    qoi_mask            np.array[float]
                        A numpy masked array for the quantity defined in var
    """

    if type(root_name) is not str:
        print("py_util.get_wind_data: the root name provided is not a string")
        exit(1)
    if type(var_type) is not str:
        print("py_util.get_wind_data: the type of data is not a string")
        exit(1)
    if type(var) is not str:
        print("py_util.get_wind_data: the type of var is not a string")
        exit(1)

    # Construct the name of the data file
    if var_type.lower() == "ion":
        ele_idx = var.find("_")
        element = var[:ele_idx]
        ion_level = var[ele_idx + 1:]
        file = "{}/{}.0.{}.txt".format(path, root_name, element)
    elif var_type.lower() == "wind":
        file = "{}/{}.ep.complete".format(path, root_name)
    else:
        print("py_util.get_wind_data: type {} not recognised for var {}".format(var_type, var))
        exit(1)

    # Check if the file already exists -- to avoid using windsave2table on wind_saves from different versions
    file_exists = os.path.isfile(file)
    if not file_exists:
        print("py_util.get_wind_data: file {} doesn't exist for var {}".format(file, var))
        rc = np.zeros((10, 10))
        return rc, rc, rc

    # Try to open the data file
    try:
        data = pd.read_table(file, delim_whitespace=True)
    except IOError:
        print("py_util.get_wind_data: could not open file {} for var {}".format(file, var))
        exit(1)

    # Read the wind data out of the file
    try:
        xi = data["i"]
        zj = data["j"]
        nx_cells = int(np.max(xi) + 1)
        nz_cells = int(np.max(zj) + 1)
        x = data["x"].values.reshape(nx_cells, nz_cells)
        z = data["z"].values.reshape(nx_cells, nz_cells)
        if var_type.lower() == "ion":
            qoi = data[ion_level].values.reshape(nx_cells, nz_cells)
        elif var_type.lower() == "wind":
            qoi = data[var].values.reshape(nx_cells, nz_cells)
        else:  # for safety, I guess?
            print("py_util.get_wind_data: type {} not recognised".format(var_type))
            exit(1)
    except KeyError:
        print("py_util.get_wind_data: could not find var {} or another key".format(var))
        exit(1)

    # construct mask for data wanted
    inwind = data["inwind"].values.reshape(nx_cells, nz_cells)
    mask = (inwind < 0)
    qoi_mask = np.ma.masked_where(mask, qoi)

    return x, z, qoi_mask


if __name__ == "__main__":
    tests()
