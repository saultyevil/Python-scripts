#!/usr/bin/env python

from platform import system
from typing import List
from PyPython import SpectrumUtils

SKIP = -1


def get_fiducial_uv_model():
    """Load the spectrum for the fiducial uv model"""

    if system() == "Darwin":
        original_model = "/Users/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.spec"
    else:
        original_model = "/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.spec"

    t = SpectrumUtils.read_spec(original_model)

    return t


def get_the_models(directories: List[str], generic_file_name: str):
    """Get the spectra of interest from file"""

    if system() == "Darwin":
        pdir = "/Users/saultyevil/PySims/tde_optical/grid/round1/"
    else:
        pdir = "/home/saultyevil/PySims/tde_optical/grid/round1/"

    modelspecs = []
    for i in range(len(directories)):
        directories[i] = pdir + directories[i] + "/" + generic_file_name
        if generic_file_name.find(".spec") != -1:
            try:
                modelspecs.append(SpectrumUtils.read_spec(directories[i]))
            except:
                modelspecs.append(SKIP)
        else:
            try:
                modelspecs.append(directories[i])
            except:
                modelspecs.append(SKIP)

    return modelspecs


# BLACK HOLE MASS

mbh_grid = [
    "Mbh/1.0000e+06",
    "Mbh/1.0000e+07",
    "Mbh/1.0000e+08"
]

mbh_labels = [
    r"M$_{BH}$ = 10$^6$ M$_{\odot}$",
    r"M$_{BH}$ = 10$^7$ M$_{\odot}$",
    r"M$_{BH}$ = 10$^8$ M$_{\odot}$",
]

# INNER RADIUS

rmin_grid = [
    "Rmin/5.0000e+00",
    "Rmin/1.0000e+01",
    "Rmin/1.5000e+01"
]

rmin_labels = [
    r"R$_{min}$ = 5 R$_{ISCO}$",
    r"R$_{min}$ = 10 R$_{ISCO}$",
    r"R$_{min}$ = 15 R$_{ISCO}$",
]

# TERMINAL VELOCITY

vinf_grid = [
    "Vinf/1.0000e-01",
    "Vinf/5.0000e-01",
    "Vinf/8.0000e-01"
]

vinf_labels = [
    r"V$_{\infty}$ = 0.1 V$_{esc}$",
    r"V$_{\infty}$ = 0.5 V$_{esc}$",
    r"V$_{\infty}$ = 0.8 V$_{esc}$"
]

# FIDUCIAL MODEL FOR EACH GRID

best_lines_grid = [
    "Mbh/1.0000e+07",
    "Rmin/1.5000e+01",
    "Vinf/1.0000e-01"
]

best_lines_labels = [
    r"M$_{BH}$ = 10$^7$ M$_{\odot}$",
    r"R$_{min}$ = 15 R$_{ISCO}$",
    r"V$_{\infty}$ = 0.1 V$_{esc}$"
]
