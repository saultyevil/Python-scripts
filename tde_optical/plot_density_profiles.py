#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import pandas as pd
from sys import exit
from platform import system
from sys import argv
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
import numpy as np
from PyPython import SpectrumUtils, WindUtils
from astropy.io import ascii
from PyPython import PythonUtils as Utils

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

MPROT = 1.672661e-24


def sightline_coords(x: np.ndarray, inclination: float) \
        -> np.ndarray:
    """
    Return the z coordinates for a given inclination angle and x coordinates.
    """

    return x * np.tan(np.pi / 2 - np.deg2rad(inclination))


def extract_density_profile(t: pd.DataFrame, inclination: float, density_type: str = "ne") \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract the requested density from the table t."""

    stride = np.max(t["j"]) + 1

    x = t["x"][::stride]
    z = sightline_coords(x, inclination)
    density = np.zeros_like(z)

    if len(x) != len(z):
        print("len(x) =", len(x))
        print("len(z) =", len(z))
        exit(1)

    for i in range(len(x)):
        j = 0
        while t["x"][j] < x[i]:
            j += 1
        k = 0
        tt = t[j:j+stride]
        while tt["z"][k] < z[i]:
            k += 1
        index = j + k
        try:
            label = density_type
            if density_type == "nh":
                label = "rho"
            density[i] = t[label][index]
            if density_type == "nh":
                density[i] /= MPROT
        except KeyError:
            print("The density {} is unknown or not implemented yet".format(density_type))
            exit(1)

    return z, z, density


def plot_density_profile_inclination(models: List[str], inclination: Union[str, float, int], labels: List[str],
                                     filename: str, density_type: str = "ne") \
        -> Tuple[plt.Figure, plt.Axes]:
    """Plot the density profiles along a specific inclination for the provided
    models."""

    density_type = density_type.lower()

    try:
        inclination = float(inclination)
    except ValueError as e:
        print(e)
        print("Unable to convert the provided inclination into a floating point number")
        exit(1)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    r = 0  # Scope related hack ;-)

    for i, m in enumerate(models):
        grid = ascii.read(m, format="basic", data_start=1)
        x, z, density = extract_density_profile(grid, inclination, density_type)
        r = np.sqrt(x ** 2 + z ** 2)
        ax.loglog(r[density != 0], density[density != 0], label=labels[i], linewidth=3)

    # ax.set_xlim(np.min(r[r != 0]), np.max(r))
    if density_type == "rho":
        ax.set_ylabel(r"Mass Density $\rho$ [g / cm$^{-3}$]", fontsize=15)
    elif density_type == "ne":
        ax.set_ylabel(r"Electron Number Density $n_{e}$ [cm$^{-3}$]", fontsize=15)
    else:
        ax.set_ylabel(r"Hydrogen Number Density $N_{H}$ [cm$^{-3}$]", fontsize=15)
    ax.set_xlabel(r"Cylindrical Radius $R$ [cm]", fontsize=15)
    ax.legend(fontsize=15)

    fig.tight_layout()

    plt.savefig("density/" + filename + "_i_{}".format(inclination) + ".png")
    plt.close()

    return fig, ax


if __name__ == "__main__":

    mbh_grid = [
        "/home/saultyevil/PySims/tde_optical/grid/round1/Mbh/1.0000e+06/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Mbh/1.0000e+07/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Mbh/1.0000e+08/tde_uv.master.txt"
    ]

    mbh_labels = [
        r"M$_{BH}$ = 10$^6$ M$_{\odot}$",
        r"M$_{BH}$ = 10$^7$ M$_{\odot}$",
        r"M$_{BH}$ = 10$^8$ M$_{\odot}$",
    ]

    rmin_grid = [
        "/home/saultyevil/PySims/tde_optical/grid/round1/Rmin/5.0000e+00/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Rmin/1.0000e+01/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Rmin/1.5000e+01/tde_uv.master.txt"
    ]

    rmin_labels = [
        r"R$_{min}$ = 5 R$_{ISCO}$",
        r"R$_{min}$ = 10 R$_{ISCO}$",
        r"R$_{min}$ = 15 R$_{ISCO}$",
    ]

    vinf_grid = [
        "/home/saultyevil/PySims/tde_optical/grid/round1/Vinf/1.0000e-01/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Vinf/5.0000e-01/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Vinf/8.0000e-01/tde_uv.master.txt"
    ]

    vinf_labels = [
        r"V$_{\infty}$ = 0.1 V$_{esc}$",
        r"V$_{\infty}$ = 0.5 V$_{esc}$",
        r"V$_{\infty}$ = 0.8 V$_{esc}$"
    ]

    # Mbh grid

    fig, ax = plot_density_profile_inclination(mbh_grid.copy(), "60", mbh_labels, "Mbh_ne", "ne")
    fig, ax = plot_density_profile_inclination(mbh_grid.copy(), "60", mbh_labels, "Mbh_Nh", "nh")
    fig, ax = plot_density_profile_inclination(mbh_grid.copy(), "60", mbh_labels, "Mbh_rho", "rho")

    # Rmin grid

    fig, ax = plot_density_profile_inclination(rmin_grid.copy(), "60", rmin_labels, "Rmin_ne", "ne")
    fig, ax = plot_density_profile_inclination(rmin_grid.copy(), "60", rmin_labels, "Rmin_Nh", "nh")
    fig, ax = plot_density_profile_inclination(rmin_grid.copy(), "60", rmin_labels, "Rmin_rho", "rho")

    # Vinf grid

    fig, ax = plot_density_profile_inclination(vinf_grid.copy(), "60", vinf_labels, "Vinf_ne", "ne")
    fig, ax = plot_density_profile_inclination(vinf_grid.copy(), "60", vinf_labels, "Vinf_Nh", "nh")
    fig, ax = plot_density_profile_inclination(vinf_grid.copy(), "60", vinf_labels, "Vinf_rho", "rho")

    # Plot the "best" lines

    files = [
        "/home/saultyevil/PySims/tde_optical/grid/round1/Mbh/1.0000e+07/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Rmin/1.5000e+01/tde_uv.master.txt",
        "/home/saultyevil/PySims/tde_optical/grid/round1/Vinf/1.0000e-01/tde_uv.master.txt"
    ]

    labels = [
        r"M$_{BH}$ = 10$^7$ M$_{\odot}$",
        r"R$_{min}$ = 15 R$_{ISCO}$",
        r"V$_{\infty}$ = 0.1 V$_{esc}$"
    ]

    fig, ax = plot_density_profile_inclination(files, "60", labels, "BestLines_ne", "ne")
    fig, ax = plot_density_profile_inclination(files, "60", labels, "BestLines_nh", "nh")
    fig, ax = plot_density_profile_inclination(files, "60", labels, "BestLines_rho", "rho")
