#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import pandas as pd
from sys import exit, argv
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
import numpy as np
from astropy.io import ascii
from path import *

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

MPROT = 1.672661e-24


def sightline_coords(x: np.ndarray, inclination: float):
    """
    Return the z coordinates for a given inclination angle and x coordinates.
    """

    return x * np.tan(np.pi / 2 - np.deg2rad(inclination))


def extract_density_profile(t: pd.DataFrame, inclination: float, density_type: str = "ne") \
        -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Extract the requested density from the table t."""

    stride = np.max(t["j"]) + 1
    x = np.array(t["x"][::stride])
    z = sightline_coords(x, inclination)
    density = np.zeros_like(z)

    assert (len(x) == len(z))

    index = 0
    prev_index = 0
    col = 0

    for i in range(len(x)):
        j = 0
        while t["x"][j] < x[i]:
            j += 1
            if j > len(t["x"]):
                j = -1
                break
        if j == -1:
            continue
        k = 0
        tt = t[j:j + stride]
        while tt["z"][k] < z[i]:
            k += 1
            if k > stride - 1:
                k = -1
                break
        if k == -1:
            continue
        prev_index = index
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

        # prev_r =

    return z, z, density


def plot_density_profile_inclination(directories: List[str], inclination: Union[str, float, int], labels: List[str],
                                     filename: str, density_type: str = "ne") \
        -> Tuple[plt.Figure, plt.Axes]:
    """Plot the density profiles along a specific inclination for the provided
    models."""

    density_type = density_type.lower()
    models = get_the_models(directories, "tde_uv.master.txt")

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    try:
        inclination = float(inclination)
    except ValueError as e:
        print(e)
        print("Unable to convert the provided inclination into a floating point number")
        return fig, ax

    for i, m in enumerate(models):
        print(m, inclination, density_type)
        grid = ascii.read(m, format="basic", data_start=1)

        x, z, density = extract_density_profile(grid, inclination, density_type)
        r = np.sqrt(x ** 2 + z ** 2)
        ax.loglog(r[density != 0], density[density != 0], label=labels[i], linewidth=3)
        ax.text(0.1, 0.1, "i = {}".format(inclination), transform=ax.transAxes)
       
    gridf = ascii.read("/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.master.txt",
                       format="basic", data_start=1)
    xf, zf, densityf = extract_density_profile(gridf, inclination, density_type)
    rf = np.sqrt(xf ** 2 + zf ** 2)
    ax.loglog(rf[densityf != 0], densityf[densityf != 0], label="UV Model", linewidth=3, color="k")

    if density_type == "rho":
        ax.set_ylabel(r"Mass Density $\rho$ [g / cm$^{-3}$]", fontsize=15)
    elif density_type == "ne":
        ax.set_ylabel(r"Electron Number Density $n_{e}$ [cm$^{-3}$]", fontsize=15)
    else:
        ax.set_ylabel(r"Hydrogen Number Density $N_{H}$ [cm$^{-3}$]", fontsize=15)
    ax.set_xlabel(r"Cylindrical Radius $R$ [cm]", fontsize=15)
    ax.legend(fontsize=15)

    fig.tight_layout()

    fig.savefig("density/" + filename + "_i{}".format(inclination) + ".pdf", dpi=300)
    fig.savefig("density/" + filename + "_i{}".format(inclination) + ".png", dpi=300)
    plt.close()

    return fig, ax


def main(argc: int, argv: List[str]) -> None:
    """
    Main function of the script

    Parameters
    ----------
    argc: int
        The number of command line arguments provided
    argv: List[str]
        The command line arguments provided
    """

    incls = ["20", "35", "60", "75"]

    # Mbh grid
    for i in incls:
        plot_density_profile_inclination(mbh_grid.copy(), i, mbh_labels, "ne_Mbh", "ne")
        plot_density_profile_inclination(mbh_grid.copy(), i, mbh_labels, "rho_Mbh", "rho")
        
        # Rmin grid
        
        plot_density_profile_inclination(rmin_grid.copy(), i, rmin_labels, "ne_Rmin", "ne")
        plot_density_profile_inclination(rmin_grid.copy(), i, rmin_labels, "rho_Rmin", "rho")
        
        # Vinf grid
        
        plot_density_profile_inclination(vinf_grid.copy(), i, vinf_labels, "ne_Vinf", "ne")
        plot_density_profile_inclination(vinf_grid.copy(), i, vinf_labels, "rho_Vinf", "rho")
        
        # Plot the "best" lines
        
        plot_density_profile_inclination(best_lines_grid.copy(), i, best_lines_labels, "ne_zBestLines", "ne")
        plot_density_profile_inclination(best_lines_grid.copy(), i, best_lines_labels, "rho_zBestLines", "rho")

    return


if __name__ == "__main__":
    main(len(argv), argv)
