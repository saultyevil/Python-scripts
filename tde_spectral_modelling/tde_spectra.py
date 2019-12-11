#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This contains various functions for getting TDE spectra from somewhere on my
hard disk.
"""


from sys import exit
import numpy as np
from PyPython import SpectrumUtils
from socket import gethostname
from platform import system


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
        print("py_util.iptf15af_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.iptf15af_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/Blagorodnova_iPTF15af.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("System: {}".format(sys))
        print("iPTF15af spectra being read in from {}".format(spec_dir))

    try:
        spec = np.loadtxt(spec_dir)
    except IOError:
        print("py_util.get_iPTF15af_spec: Unable to open the iPTF15af UV spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    spec[:, 1] = SpectrumUtils.smooth_spectrum(spec[:, 1], smooth)

    return spec


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
        print("py_util.asassn14li_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.iptf15af_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/ASASSN-14li_spec_Cenko.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("Cenko spectra being read in from {}".format(spec_dir))

    try:
        spec = np.loadtxt(spec_dir)
    except IOError:
        print("py_util.get_ASSASN_14li_spec: Unable to open the ASSASSN_14li spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    spec[:, 1] = SpectrumUtils.smooth_spectrum(spec[:, 1], smooth)

    return spec


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
        print("py_util.iptf16fnl_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.iptf16fnl_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/iPTF16fnl_52d.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("iPTF16fnl spectra being read in from {}".format(spec_dir))

    try:
        spec = np.loadtxt(spec_dir)
    except IOError:
        print("py_util.iptf16fnl_spec: Unable to open the iPTF16fnl spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    spec[:, 1] = SpectrumUtils.smooth_spectrum(spec[:, 1], smooth)

    return spec


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
        print("py_util.at2018zr_spec: Unable to convert smooth into an integer")
        exit(1)

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.iptf16fnl_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/at2018zr_59d.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("AT2018zr spectra being read in from {}".format(spec_dir))

    try:
        spec = np.loadtxt(spec_dir)
    except IOError:
        print("py_util.at2018zr_spec: Unable to open the AT2018zr spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    spec[:, 1] = SpectrumUtils.smooth_spectrum(spec[:, 1], smooth)

    return spec


def sdss_qso_spec(verbose: bool = False) -> np.array:
    """
    Parameters
    ----------
    verbose             bool, optional
                        Enable verbose logging
    """

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.sdss_qso_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/sdss_composite_qso.dat"


    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("SDSS QSO spectra being read in from {}".format(spec_dir))

    try:
        return np.loadtxt(spec_dir)
    except IOError:
        print("py_util.sdss_qso_spec: Unable to open the SDSS QSO spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    return


def lobal_qso_spec(verbose: bool = False) -> np.array:
    """
    Parameters
    ----------
    verbose             bool, optional
                        Enable verbose logging
    """

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.lobal_qso_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/LoBALQSO.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("LoBALQSO spectra being read in from {}".format(spec_dir))

    try:
        return np.loadtxt(spec_dir, skiprows=1)
    except IOError:
        print("py_util.lobal_qso_spec: Unable to open the LoBALQSO spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    return


def n_rich_qso_spec(verbose: bool = False) -> np.array:
    """
    UV spectrum for the Nitrogen Rich Quasar SDSS J164148.20+223225.22.

    Parameters
    ----------
    verbose             bool, optional
                        Enable verbose logging
    """

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.n_rich_qso_spec: unknown system type {}".format(sys))
        exit(1)
    spec_dir += "PySims/tde/observed_spec/sdss_J164148.20+223225.22_n_rich_QSO.dat"
    # spec_dir += "PySims/tde/observed_spec/nrichqso.dat"

    if verbose:
        print("Hostname: {}".format(gethostname()))
        print("N Rich QSO spectra being read in from {}".format(spec_dir))

    try:
        return np.loadtxt(spec_dir)
    except IOError:
        print("py_util.n_rich_qso_spec: Unable to open the N Rich QSO spectrum from the following path {}. "
              "Update the directories in the script".format(spec_dir))
        exit(1)

    return


def get_tde_spec(name: str, plot: bool = False, verbose: bool = False):
    """

    Parameters
    ----------
    name
    plot
    verbose

    Returns
    -------

    """

    spec_dir = ""
    sys = system()
    if sys == "Linux":
        spec_dir = "/home/saultyevil/"
    elif sys == "Darwin":
        spec_dir = "/Users/saultyevil/"
    else:
        print("py_util.lobal_qso_spec: unknown system type {}".format(sys))
        exit(1)
    

    return


def plot_iPTF15af_photometry(ax):
    """
    Plot photometric data for the UV spectrum of iPTF15af from Blagordnova
    et al. 2019.

    Note that the Swift optical filters (U, B, V) are commented out due to
    contamination from the host galaxy.

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
