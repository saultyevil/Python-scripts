#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Various functions used throughout plotting the output from Python. This script
should be imported into other scripts rather than being itself run.
"""


import os
import numpy as np
import pandas as pd
from shutil import which
from subprocess import Popen, PIPE
from matplotlib import pyplot as plt
from typing import Tuple, List, Union
from scipy.signal import convolve, boxcar
from pathlib import Path
from consts import C, ANGSTROM


def tests() -> None:
    """
    Warns the user that this script is a utility script and not meant to be run.
    It will also print the current Python version to screen.
    """

    print("This script is not designed to be run. Instead, import it using import py_plot_util.")
    get_python_version(verbose=True)

    return


def get_python_version(py: str = "py", verbose: bool = False) -> Tuple[str, str]:
    """
    Get the Python version and commit hash for the provided Python binary. This should also work with windsave2table.

    Parameters
    ----------
    py: str, optional
        The name of the Python executable in $PATH whose version will be queried
    verbose: bool, optional
        Enable verbose logging

    Returns
    --------
    version: str
        The version number of Python
    commit_hash: str
        The commit hash of Python
    """

    version = ""
    commit_hash = ""

    path = which(py)
    if not path:
        raise Exception("{}:{}: {} is not in $PATH".format(__file__, get_python_version.__name__, py))

    command = "{} --version".format(py)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    out = stdout.decode("utf-8").split()
    err = stderr.decode("utf-8")

    if err:
        print("{}:{}: captured from stderr".format(__file__, get_python_version.__name__))
        print(stderr)

    for i in range(len(out)):
        if out[i] == "Version":
            version = out[i + 1]
        if out[i] == "hash":
            commit_hash = out[i + 1]

    if version == "" and verbose:
        print("{}:{}: couldn't find version for {}".format(__file__, get_python_version.__name__, py))
    if commit_hash == "" and verbose:
        print("{}:{}: couldn't find commit hash for {}".format(__file__, get_python_version.__name__, py))

    if verbose and version and commit_hash:
        print("{} version {}".format(py, version))
        print("Git commit hash    {}".format(commit_hash))
        print("Short commit hash  {}".format(commit_hash[:7]))

    return version, commit_hash


def get_root_name(pf_path: str) -> Tuple[str, str]:
    """
    Split a path name into a directory path and root name for a Python simulation.

    Parameters
    ----------
    pf_path: str
        The directory path to a Python .pf file

    Returns
    -------
    root: str
        The root name of the Python simulation
    path: str
        The directory path containing the provided Python .pf file
    """

    if type(pf_path) != str:
        raise TypeError("{}:{}: expected string as input".format(__file__, get_root_name.__name__))

    dot = 0
    slash = 0
    for i in range(len(pf_path)):
        l = pf_path[i]
        if l == ".":
            dot = i
        elif l == "/":
            slash = i + 1

    root = pf_path[slash:dot]
    path = pf_path[:slash]
    if path == "":
        path = "./"

    return root, path


def find_specs(path: str = "./") -> List[str]:
    """
    Find root.spec files recursively in provided directory

    Parameters
    ----------
    path: str
        The path to recursively search from

    Returns
    -------
    spec_files: List[str]
        The file paths of various .spec files
    """

    spec_files = []
    for filename in Path(path).glob("**/*.spec"):
        spec_files.append(str(filename))

    return spec_files


def find_pf(path: str = "./") -> List[str]:
    """
    Find parameter files recursively from the directory this function is called
    in.

    Parameters
    ----------
    path: str, optional
        The directory to search for pf files

    Returns
    -------
    pfs: List[str]
        The file path for any Python pf files founds
    """

    pfs = []
    for filename in Path(path).glob("**/*.pf"):
        fname = str(filename)
        if fname.find("out.pf") != -1:
            continue
        if fname == "py_wind.pf":
            continue
        if fname[0] == "/":
            fname = fname[1:]
        pfs.append(fname)
    pfs = sorted(pfs, key=str.lower)

    return pfs


def read_spec(file_name: str, delim: str = None, numpy: bool = False) -> Union[np.ndarray, pd.DataFrame]:
    """
    Read in data from an external file, line by line whilst ignoring comments.
        - Comments begin with #
        - The default delimiter is assumed to be a space

    Parameters
    ----------
    file_name: str
        The directory path to the spec file to be read in
    delim: str, optional
        The delimiter between values in the file, by default a space is assumed
    numpy:bool, optional
        If True, a Numpy array of strings will be used instead :-(

    Returns
    -------
    lines: np.ndarray or pd.DataFrame
        The .spec file as a Numpy array or a Pandas DataFrame
    """

    try:
        with open(file_name, "r") as f:
            flines = f.readlines()
    except IOError:
        raise Exception("{}:{}: cannot open spec file {}".format(__file__, read_spec.__name__, file_name))

    lines = []
    for i in range(len(flines)):
        line = flines[i].strip()
        if delim:
            line = line.split(delim)
        else:
            line = line.split()
        if len(line) > 0:
            if line[0] == "#":
                continue
            if line[0] == "Freq.":
                for j in range(len(line)):
                    if line[j][0] == "A":
                        index = line[j].find("P")
                        line[j] = line[j][1:index]
            lines.append(line)

    if numpy:
        return np.array(lines)

    return pd.DataFrame(lines[1:], columns=lines[0])


def spec_inclinations(spec: Union[pd.DataFrame, np.ndarray]) -> List[str]:
    """
    Return a list of the inclination angles in the spectrum as strings. This is
    done differently depending on if the provided spectrum is a Numpy array or
    a Pandas DataFrame.

    Parameters
    ----------
    spec: pd.DataFrame or np.ndarray
        The spectrum to extract the inclination angles from

    Returns
    -------
    inclinations: List[str]
        The possible inclination angles in the spectrum
    """

    inclinations = []
    arr_type = type(spec)

    if arr_type == pd.core.frame.DataFrame:
        column_headers = spec.columns.values
        for header in column_headers:
            if header.isdigit():
                inclinations.append(header)
    elif arr_type == np.ndarray:
        col_names = spec[0, :]
        for i in range(len(col_names)):
            if col_names[i].isdigit() is True:
                angle = col_names[i]
                duplicate_flag = False
                for va in inclinations:
                    if angle == va:
                        duplicate_flag = True
                if not duplicate_flag:
                    inclinations.append(angle)
    else:
        raise TypeError("{}:{}: unknown data type {}, require pd.DataFrame or np.ndarray"
                        .format(__file__, spec_inclinations.__name__, arr_type))

    return inclinations


def all_inclinations(spec_names: List[str], delim: str = None) -> np.array:
    """
    Get all of the unique inclination angles for a set of Python .spec files.
    Parameters
    ----------
    spec_name: List[str]
        The directory path to Python .spec files
    delim: str, optional
        The delimiter in the .spec files, assumed to be spaces by default
    Returns
    -------
    inclinations: List[int]
        All of the unique inclination angles found in the Python .spec files
    """

    inclinations = []

    for i in range(len(spec_names)):
        if spec_names[i].find(".pf") != -1:
            spec_names[i] = spec_names[i].replace(".pf", ".spec")

    # Find the viewing angles in each .spec file
    for i in range(len(spec_names)):
        spec = read_spec(spec_names[i], delim)

        if type(spec) == pd.core.frame.DataFrame:
            col_names = spec.columns.values
        elif type(spec) == np.ndarray:
            col_names = spec[0, :]
        else:
            raise TypeError("{}:{}: unknown data type {} for function"
                            .format(__file__, all_inclinations.__name__, type(spec)))

        # Go over the columns and look for viewing angles
        for j in range(len(col_names)):
            if col_names[j].isdigit() is True:
                angle = int(col_names[j])
                duplicate_flag = False
                for va in inclinations:  # Check for duplicate angle
                    if angle == va:
                        duplicate_flag = True
                if duplicate_flag is False:
                    inclinations.append(angle)

    return inclinations


def check_inclination(inclination: str, spec: Union[pd.DataFrame, np.ndarray]) -> bool:

    """
    Check that an inclination angle is in the spectrum array.

    Parameters
    ----------
    inclination: str
        The inclination angle to check
    spec: np.ndarray
        The spectrum array to read -- assume that it is a np.array of strings
        Note that tde_spec_plot has a similar routine for pd.DataFrame's, whoops!

    Returns
    -------
    allowed: bool
        If True, angle is a legal angle, otherise fault
    """

    allowed = False

    if type(spec) == pd.core.frame.DataFrame:
        headers = spec.columns.values
    elif type(spec) == np.ndarray:
        headers = spec[0, :]
    else:
        raise TypeError("{}:{}: unknown data type {} for function"
                        .format(__file__, check_inclination.__name__, type(spec)))

    if type(inclination) != str:
        try:
            inclination = str(inclination)
        except ValueError:
            raise TypeError("{}:{}: could not convert {} into string"
                            .format(__file__, check_inclination.__name__, inclination))

    if inclination in headers:
        allowed = True

    return allowed


def smooth(flux: np.ndarray, smooth: Union[int, float], verbose: bool = False) -> np.ndarray:
    """
    Smooth a 1D array of data using a boxcar filter of width smooth pixels.

    Parameters
    ----------
    flux: np.array[float]
        The data to smooth using the boxcar filter
    smooth: int
        The size of the window for the boxcar filter
    verbose: bool
        Enable verbose logging

    Returns
    -------
    smoothed: np.ndarray
        The smoothed data
    """

    if type(flux) != list and type(flux) != np.ndarray:
        raise TypeError("{}:{}: expecting list or np.ndarray".format(__file__, smooth.__name__))

    if type(flux) == list:
        flux = np.array(flux, dtype=float)

    if len(flux.shape) > 2:
        raise Exception("{}:{}: data is not 1 dimensional but has shape {}"
                        .format(__file__, smooth.__name__, flux.shape))

    if type(smooth) != int:
        try:
            smooth = int(smooth)
        except ValueError:
            raise Exception("{}:{}: could not convert smooth = {} into an integer"
                            .format(__file__, smooth.__name__, smooth))

    flux = np.reshape(flux, (len(flux),))
    smoothed = convolve(flux, boxcar(smooth) / float(smooth), mode="same")

    return smoothed


def windsave2table(root: str, path: str, verbose: bool = False) -> None:
    """
    Run windsave2table in the directory given by path. This function will also
    create a *.ep.complete file which combines both the heat and master data
    tables together.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    path: str
        The directory of the Python simulation where windsave2table will be run
    verbose: bool, optional
        Enable verbose logging
    """

    version, hash = get_python_version("windsave2table", verbose)
    try:
        with open("version", "r") as f:
            lines = f.readlines()
        run_version = lines[0]
        run_hash = lines[1]
        if run_version != version and run_hash != hash:
            print("{}:{}: windsave2table and wind_save versions are likey different. Output may be wank."
                  .format(__file__, windsave2table.__name__))
    except IOError:
        print("{}:{}: unable to determine windsave2table version from py_run.py version file"
              .format(__file__, windsave2table.__name__))

    in_path = which("windsave2table")
    if not in_path:
        raise Exception("{}:{}: windsave2table not in $PATH and executable", __file__, windsave2table.__name__)

    command = "cd {}; Setup_Py_Dir; windsave2table {}".format(path, root)
    cmd = Popen(command, stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = cmd.communicate()
    output = stdout.decode("utf-8")
    err = stderr.decode("utf-8")

    if verbose:
        print(output)
    if err:
        print("{}:{}: the following was sent to stderr:".format(__file__, windsave2table.__name__))
        print(err)

    # Now create a "complete" file which is the master and heat put together into one csv
    heat_file = "{}/{}.0.heat.txt".format(path, root)
    master_file = "{}/{}.0.master.txt".format(path, root)

    try:
        heat = pd.read_csv(heat_file, delim_whitespace=True)
        master = pd.read_csv(master_file, delim_whitespace=True)
    except IOError:
        raise IOError("{}: could not open master or heat file for root {}".format(__file__, windsave2table.__name__, root))

    # This merges the heat and master table together :-)
    append = heat.columns.values[14:]
    for i, col in enumerate(append):
        master[col] = pd.Series(heat[col])
    master.to_csv("{}/{}.ep.complete".format(path, root), sep=" ")

    return


def extract_wind_var(root: str, var_name: str, var_type: str, path: str = "./", coord: str = "rectilinear") \
        -> Tuple[np.array, np.array, np.array]:
    """
    Read in variables contained within a windsave2table file. Requires the user
    to have already run windsave2table so the data is in the directory. This
    also assumes that a .ep.complete file exists which contains both the
    heat and master data. This will also only work for 2d models :^).

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    var_name: str
        The name of the quantity from the Python simulation
    var_type: str
        The type of quantity to extract, this can either be wind or ion
    path: str, optional
        The directory containing the Python simulation
    coord: str
        The coordinate system in use. Currently this only works for polar
        and rectilinear coordinate systems

    Returns
    -------
    x: np.array[float]
        The x coordinates of the wind in the Python simulation
    z: np.array[float]
        The z coordinates of the wind in the Python simulation
    var_mask: np.array[float]
        A numpy masked array for the quantity defined in var
    """

    assert(type(root) == str), "{}:{}: root must be a string".format(__file__, extract_wind_var.__name__)
    assert(type(var_name) == str), "{}:{}: var_name must be a string".format(__file__, extract_wind_var.__name__)
    assert(type(var_type) == str), "{}:{}: var_type must be a string".format(__file__, extract_wind_var.__name__)

    # Create string containing the name of the data file required to be loaded
    nx_cells = 0
    nz_cells = 0
    key = var_name
    if var_type.lower() == "ion":
        ele_idx = var_name.find("_")
        element = var_name[:ele_idx]
        key = var_name[ele_idx + 1:]
        file = "{}/{}.0.{}.txt".format(path, root, element)
    elif var_type.lower() == "wind":
        file = "{}/{}.ep.complete".format(path, root)
    else:
        raise Exception("{}:{}: var type {} not recognised for var {}"
                        .format(__file__, extract_wind_var.__name__, var_type, var_name))

    file_exists = os.path.isfile(file)
    if not file_exists:
        raise IOError("{}:{}: file {} doesn't exist. var {} var type {}".
                      format(__file__, extract_wind_var.__name__, file, var_name, var_type))

    # Open the file and remove any garbage cells for theta grids
    try:
        data = pd.read_csv(file, delim_whitespace=True)
        if coord == "polar" and var_type != "ion":
            data = data[~(data["theta"] > 90)]
    except IOError:
        raise IOError("{}:{}: could not open file {} for some reason"
                      .format(__file__, extract_wind_var.__name__, file, var_name))

    # Now we can try and read the data out of the file...
    try:
        xi = data["i"]
        zj = data["j"]
        nx_cells = int(np.max(xi) + 1)
        nz_cells = int(np.max(zj) + 1)
        if coord == "rectilinear":
            x = data["x"].values.reshape(nx_cells, nz_cells)
            z = data["z"].values.reshape(nx_cells, nz_cells)
        elif coord == "polar":
            # Transform from cartesian to polar coordinates manually for ion tables
            if var_type.lower() == "ion":
                x = data["x"].values.reshape(nx_cells, nz_cells)
                z = data["z"].values.reshape(nx_cells, nz_cells)
                r = np.sqrt(x ** 2 + z ** 2)
                theta = np.rad2deg(np.arctan(z / x))
                x = r
                z = theta
            else:
                try:
                    x = data["r"].values.reshape(nx_cells, nz_cells)
                    z = data["theta"].values.reshape(nx_cells, nz_cells)
                except KeyError:
                    print("{}:{}: trying to read r theta  points in non-polar model"
                          .format(__file__, extract_wind_var.__name__))
        else:
            raise Exception("{}:{}: unknown projection {}: use rectilinear or polar"
                            .format(__file__, extract_wind_var.__name__, coord))
    except KeyError:
        print("{}:{}: could not find var {} or another key".format(__file__, extract_wind_var.__name__, var_name))

    try:
        var = data[key].values.reshape(nx_cells, nz_cells)
    except KeyError:
        raise Exception("{}:{}: unable to use {} as key for table".format(__file__, extract_wind_var.__name__, key))

    # Construct mask for variable
    inwind = data["inwind"].values.reshape(nx_cells, nz_cells)
    mask = (inwind < 0)
    var_mask = np.ma.masked_where(mask, var)

    return x, z, var_mask


def subplot_dims(n_plots: int) -> Tuple[int, int]:
    """
    Determine the dimensions for a plot with multiple subplot panels. Two columns
    of subplots will always be used.

    Parameters
    ----------
    n_plots: int
        The number of subplots which will be plotted

    Returns
    -------
    dims: Tuple[int, int]
        The dimensions of the subplots returned as (nrows, ncols)
    """

    if n_plots > 2:
        ncols = 2
        nrows = (1 + n_plots) // ncols
    elif n_plots > 9:
        ncols = 3
        nrows = (1 + n_plots) // ncols
    else:
        ncols = 1
        nrows = n_plots

    dims = (nrows, ncols)

    return dims


def ylims(wavelength: np.array, flux: np.array, wmin: float, wmax: float, scale: float = 10,
          verbose: bool = False) -> Union[Tuple[float, float], Tuple[bool, bool]]:
    """
    Return appropriate y limits to use when the wavelength range has been limited.

    Parameters
    ----------
    wavelength: np.array[float]
        An array containing all wavelengths in a spectrum
    flux: np.array[float]
        An array containing the flux at each wavelength point
    wmin: float
        The shortest wavelength which is being plotted
    wmax: float
        The longest wavelength which is being plotted
    scale: float, optional
        The scaling factor for white space around the data
    verbose: bool, optional
        Enable verbose logging

    Returns
    -------
    yupper: float
        The upper y limit to use with the wavelength range
    ylower: float
        The lower y limit to use with the wavelength range
    """

    if wavelength.shape[0] != flux.shape[0]:
        raise Exception("{}:{}: wavelength and flux are of different dimensions {} {}"
                        .format(__file__, ylims.__name__, wavelength.shape, flux.shape))

    if type(wavelength) != np.ndarray or type(flux) != np.ndarray:
        raise TypeError("{}:{}: wavelength or flux array not a numpy arrays {} {}"
                        .format(__file__, ylims.__name__, type(wavelength), type(flux)))

    yupper = ylower = None

    if wmin and wmax:
        idx_wmin = wavelength < wmin
        idx_wmax = wavelength > wmax
        flux_lim_wav = np.where(idx_wmin == idx_wmax)[0]
        yupper = np.max(flux[flux_lim_wav]) * scale
        ylower = np.min(flux[flux_lim_wav]) / scale

    return yupper, ylower


def common_lines(freq: bool = False, log: bool = False) -> list:
    """
    Return a dictionary containing the major absorption and emission lines which
    I'm interested in. The wavelengths of the lines are in Angstrom in the
    rest frame.

    Parameters
    ----------
    freq: bool, optional
       If True, return the dict in frequency space
    log: bool, optional
       If this and use_freq are true, the lines will be returned to be plotted
       on a log10 scale

    Returns
    -------
    line: List[List[str, float]]
        A list of lists where each element of the list is the name of the
        transition/edge and the rest wavelength of that transition in Angstroms.
    """

    lines = [
        ["HeII Edge", 229],
        ["Lyman Edge", 912],
        ["P V", 1118],
        [r"Ly$\alpha$/N V", 1216],
        ["", 1240],
        ["O V/Si IV", 1371],
        ["", 1400],
        ["N IV", 1489],
        ["C IV", 1549],
        ["He II", 1640],
        ["N III]", 1750],
        ["Al III", 1854],
        ["C III]", 1908],
        ["Mg II", 2798],
        ["Balmer Edge", 3646],
        ["He II", 4686],
        [r"H$_{\beta}$", 4861],
        [r"H$_{\alpha}$", 6564],
        ["Paschen Edge", 8204]
    ]

    if freq:
        for i in range(len(lines)):
            if log:
                lines[i][1] = np.log10(C / (lines[i][1] * ANGSTROM))
            else:
                lines[i][1] = C / (lines[i][1] * ANGSTROM)

    return lines


def absorption_edges(freq: bool = False, log: bool = False) -> list:
    """
    Return a dictionary containing major absorption edges which I am interested
    in. The wavelengths of the lines are in Angstroms or in frequency in Hz
    if use_freq is True.

    Parameters
    ----------
    freq: bool, optional
       If True, return the dict in frequency space
    log: bool, optional
       If this and use_freq are true, the lines will be returned to be plotted
       on a log10 scale

    Returns
    -------
    edges: dict
        A list of lists where each element of the list is the name of the
        transition/edge and the rest wavelength of that transition in Angstroms.
    """

    edges = [
        ["HeII Edge", 229],
        ["Lyman Edge", 912],
        ["Balmer Edge", 3646],
        ["Paschen Edge", 8204],
    ]

    if freq:
        for i in range(len(edges)):
            if log:
                edges[i][1] = np.log10(C / (edges[i][1] * ANGSTROM))
            else:
                edges[i][1] = C / (edges[i][1] * ANGSTROM)

    return edges


def plot_line_ids(ax: plt.Axes, lines: list, rotation: str = "vertical", fontsize: int = 10) -> plt.Axes:
    """
    Plot line IDs onto a figure. This should probably be used after the x-limits
    have been set on the figure which these labels are being plotted onto. 

    Parameters
    ----------
    ax: plt.Axes
        The plot object to add line IDs to
    lines: list
        A list containing the line name and wavelength in Angstroms 
        (ordered by wavelength)
    rotation: str, optional
        Vertical or horizontal rotation for text ids
    fontsize: int, optional
        The fontsize of the labels

    Returns
    -------
    ax: plt.Axes
        The plot object now with lines IDs :-)
    """

    nlines = len(lines)
    xlims = ax.get_xlim()

    for i in range(nlines):
        x = lines[i][1]
        if x < xlims[0]:
            continue
        if x > xlims[1]:
            continue
        label = lines[i][0]
        ax.axvline(x, linestyle="--", linewidth=0.5, color="k", zorder=1)
        x = x - 25
        xnorm = (x - xlims[0]) / (xlims[1] - xlims[0])
        ax.text(xnorm, 0.92, label, ha="center", va="center", rotation=rotation, fontsize=fontsize,
                transform=ax.transAxes)

    return ax


if __name__ == "__main__":
    tests()
