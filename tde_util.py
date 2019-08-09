#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from sys import exit
import numpy as np
import py_plot_util
from socket import gethostname


def iptf15af_spec(smooth: int, verbose: bool = False) -> np.array:
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
                            - Column 0: wavelength in Angstroms
                            - Column 1: flux per unit wavelength in erg s^-1 cm^-2 A^-1
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.get_iPTF15af_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        spec_dir = "/home/saultyevil/PySims/tde/observed_spec/Blagorodnova_iPTF15af.dat"
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

    iPTF15af_spec[:, 1] = py_plot_util.smooth_1d_array(iPTF15af_spec[:, 1], smooth)

    return iPTF15af_spec


def asassn14li_spec(smooth: int, verbose: bool = False) -> np.array:
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
                            - Column 0: wavelength in Angstroms
                            - Column 1: flux per unit wavelength in erg s^-1 cm^-2 A^-1
                            - Column 2: error in flux
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.get_ASSASN_14li_spec: Unable to convert smooth into an integer")
        exit(1)

    cenk_dir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        cenk_dir = "/home/saultyevil/PySims/tde/observed_spec/ASASSN-14li_spec_Cenko.dat"
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

    ASSASN_14li_spec[:, 1] = py_plot_util.smooth_1d_array(ASSASN_14li_spec[:, 1], smooth)

    return ASSASN_14li_spec


def iptf16fnl_spec(smooth: int, verbose: bool = False) -> np.array:
    """
    Parameters
    ----------
    smooth              int
                        The size of the window for the boxcar filter
    verbose             bool, optional
                        Enable verbose logging
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.iPTF16fnl_spec: Unable to convert smooth into an integer")
        exit(1)

    fdir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        fdir = "/home/saultyevil/PySims/tde/observed_spec/iPTF16fnl_52d.dat"
    else:
        print("Unknown hostname, update py_util with directory for the iPTF16fnl spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("iPTF16fnl spectra being read in from {}".format(fdir))

    try:
        iPTF16fnl_spec = np.loadtxt(fdir)
    except IOError:
        print("py_util.iPTF16fnl_spec: Unable to open the iPTF16fnl spectrum from the following path {}. "
              "Update the directories in the script".format(fdir))
        exit(1)

    iPTF16fnl_spec[:, 1] = py_plot_util.smooth_1d_array(iPTF16fnl_spec[:, 1], smooth)

    return iPTF16fnl_spec


def at2018zr_spec(smooth: int, verbose: bool = False) -> np.array:
    """
    Parameters
    ----------
    smooth              int
                        The size of the window for the boxcar filter
    verbose             bool, optional
                        Enable verbose logging
    """

    try:
        smooth = int(smooth)
    except ValueError:
        print("py_util.AT2018zr_spec: Unable to convert smooth into an integer")
        exit(1)

    fdir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        fdir = "/home/saultyevil/PySims/tde/observed_spec/at2018zr_59d.dat"
    else:
        print("Unknown hostname, update py_util with directory for the AT2018zr spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("AT2018zr spectra being read in from {}".format(fdir))

    try:
        AT2018zr_spec = np.loadtxt(fdir)
    except IOError:
        print("py_util.AT2018zr_spec: Unable to open the AT2018zr spectrum from the following path {}. "
              "Update the directories in the script".format(fdir))
        exit(1)

    AT2018zr_spec[:, 1] = py_plot_util.smooth_1d_array(AT2018zr_spec[:, 1], smooth)

    return AT2018zr_spec


def sdss_qso_spec(verbose: bool = False) -> np.array:
    """
    Parameters
    ----------
    verbose             bool, optional
                        Enable verbose logging
    """

    fdir = ""
    hostname = gethostname()
    if hostname == "ASTRO-REX":
        fdir = "/home/saultyevil/PySims/tde/observed_spec/sdss_composite_qso.dat"
    else:
        print("Unknown hostname, update py_util with directory for the SDSS QSO spectrum")
        exit(1)

    if verbose:
        print("Hostname: {}".format(hostname))
        print("SDSS QSO spectra being read in from {}".format(fdir))

    try:
        return np.loadtxt(fdir)
    except IOError:
        print("py_util.SDSS_QSO: Unable to open the SDSS QSO spectrum from the following path {}. "
              "Update the directories in the script".format(fdir))
        exit(1)

    return


def plot_iPTF15af_photometry(ax):
    """
    Plot photometric data for the UV spectrum of iPTF15af from Blagordnova et al. 2019.

    Note that the Swift optical filters (U, B, V) are commented out due to contamination from the host galaxy.

    Returns
    -------

    """

    def mag2flux(mag, zero):
        return zero * 10 ** (-mag / 2.5)

    mags = np.array([
        19.48,  # UVW1
        19.42,  # UVM2
        19.16,  # UVW2
        # 19.27,   # U
        # 18.81,   # B
        # 17.91,   # V
        20.09,  # g
        20.39,  # r
        20.31,  # i
    ])

    # Taken from http://svo2.cab.inta-csic.es/svo/theory/fps/
    zero_points = np.array([
        1.624e-8,
        2.193e-8,
        2.640e-8,
        # 8.880e-9,
        # 5.810e-9,
        # 3.730e-9,
        5.055e-9,
        2.904e-9,
        7.738e-9
    ])

    # Taken from http://svo2.cab.inta-csic.es/svo/theory/fps/
    central_waves = np.array([
        2604.57,
        2246.00,
        1941.22,
        # 3463.14,
        # 4371.22,
        # 5441.20,
        4700.33,
        6174.48,
        3684.29
    ])

    labels = [
        "Swift (+52d): UVW1",
        "Swift (+52d): UWM2",
        "Swift (+52d): UVW2",
        # "Swift (+52d): U",
        # "Swift (+52d): B",
        # "Swift (+52d): V",
        "P60 (+53.9d): g",
        "P60 (+53.9d): r",
        "LasCumbres (+49.5d): I"
    ]

    fluxs = mag2flux(mags, zero_points)
    for i in range(len(fluxs)):
        leg_label = labels[i] + r" ({:4.0f} $\AA$)".format(central_waves[i])
        ax.plot(central_waves[i], fluxs[i], "d", label=leg_label)

    return ax
