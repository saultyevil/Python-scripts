#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import pandas as pd
from sys import exit, argv
from matplotlib import pyplot as plt
from typing import List, Tuple, Union
import numpy as np
from astropy.io import ascii
from path import *
from PyPython.WindUtils import extract_variable_along_sightline

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

MPROT = 1.672661e-24


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

        x, z, density = extract_variable_along_sightline(inclination, density_type, grid=grid, legacy=True)
        r = np.sqrt(x ** 2 + 0 ** 2)
        ax.loglog(r[density != 0], density[density != 0], label=labels[i], linewidth=3, alpha=0.75)
        ax.text(0.1, 0.1, "i = {}".format(inclination), transform=ax.transAxes)

    gridf = ascii.read("/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.master.txt",
                       format="basic", data_start=1)
    xf, zf, densityf = extract_variable_along_sightline(inclination, density_type, grid=gridf, legacy=True)
    rf = np.sqrt(xf ** 2 + 0 ** 2)
    ax.loglog(rf[densityf != 0], densityf[densityf != 0], label="UV Model", linewidth=2, color="k", alpha=0.75)

    if density_type == "rho":
        ax.set_ylabel(r"Mass Density $\rho$ [g / cm$^{-3}$]", fontsize=15)
    elif density_type == "ne":
        ax.set_ylabel(r"Electron Number Density $n_{e}$ [cm$^{-3}$]", fontsize=15)
    else:
        ax.set_ylabel(r"Ionization Parameter", fontsize=15)
    ax.set_xlabel(r"Cylindrical Radius $R$ [cm]", fontsize=15)
    ax.legend(fontsize=15)

    fig.tight_layout()

    fig.savefig("density/z" + filename + "_i{}".format(inclination) + ".pdf", dpi=300)
    fig.savefig("density/z" + filename + "_i{}".format(inclination) + ".png", dpi=300)
    plt.close()

    return fig, ax


def plot_density_profile(directories: List[str], inclinations: List[str], labels: List[str], filename: str,
                         density_type: str = "ne") \
        -> Tuple[plt.Figure, plt.Axes]:
    """Plot the density profiles along a specific inclination for the provided
    models."""

    density_type = density_type.lower()
    models = get_the_models(directories, "tde_uv.master.txt")

    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex="col", sharey="row")
    ax = np.reshape(ax, (4,))

    if len(inclinations) != 4:
        print("Expected four inclinations")
        return fig, ax

    for ii, inclination in enumerate(inclinations):

        try:
            inclination = float(inclination)
        except ValueError as e:
            print(e)
            print("Unable to convert the provided inclination into a floating point number")
            return fig, ax

        for i, m in enumerate(models):
            print(m, inclination, density_type)
            grid = ascii.read(m, format="basic", data_start=1)

            x, z, density = extract_variable_along_sightline(inclination, density_type, grid=grid, legacy=True)
            r = np.sqrt(x ** 2 + 0 ** 2)
            ax[ii].loglog(r[density != 0], density[density != 0], label=labels[i], linewidth=3)
            ax[ii].text(0.1, 0.1, "i = {}".format(inclination), transform=ax[ii].transAxes, fontsize=15, alpha=0.75)

        gridf = ascii.read("/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.master.txt",
                           format="basic", data_start=1)
        xf, zf, densityf = extract_variable_along_sightline(inclination, density_type, grid=gridf, legacy=True)
        rf = np.sqrt(xf ** 2 + 0 ** 2)
        ax[ii].loglog(rf[densityf != 0], densityf[densityf != 0], label="UV Model", linewidth=2, color="k",
                      alpha=0.75)

    ax[0].legend(fontsize=15)

    if density_type == "rho":
        fig.text(0.025, 0.5, r"Mass Density $\rho$ [g / cm$^{-3}$]", ha="center", va="center", rotation="vertical",
                 fontsize=15)
    elif density_type == "ne":
        fig.text(0.025, 0.5, r"Electron Number Density $n_{e}$ [cm$^{-3}$]", ha="center", va="center",
                 rotation="vertical", fontsize=15)
    else:
        fig.text(0.025, 0.5, r"Ionization Parameter", ha="center", va="center", rotation="vertical", fontsize=15)

    fig.text(0.5, 0.03, r"Cylindrical Radius $R$ [cm]", ha="center", va="center", rotation="horizontal", fontsize=15)

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0, wspace=0)

    fig.savefig("density/" + filename + ".pdf", dpi=300)
    fig.savefig("density/" + filename + ".png", dpi=300)
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

    incls = ["20", "60", "75", "89"]

    # Mbh grid

    plot_density_profile(mbh_grid.copy(), incls, mbh_labels, "ne_Mbh", "ne")
    plot_density_profile(mbh_grid.copy(), incls, mbh_labels, "rho_Mbh", "rho")
    plot_density_profile(mbh_grid.copy(), incls, mbh_labels, "ip_Mbh", "ip")

    plot_density_profile(mbh_fixed_grid.copy(), incls, mbh_labels, "ne_Mbh_fixed", "ne")
    plot_density_profile(mbh_fixed_grid.copy(), incls, mbh_labels, "rho_Mbh_fixed", "rho")
    plot_density_profile(mbh_fixed_grid.copy(), incls, mbh_labels, "ip_Mbh_fixed", "ip")

    # Rmin grid

    plot_density_profile(rmin_grid.copy(), incls, rmin_labels, "ne_Rmin", "ne")
    plot_density_profile(rmin_grid.copy(), incls, rmin_labels, "rho_Rmin", "rho")
    plot_density_profile(rmin_grid.copy(), incls, rmin_labels, "ip_Rmin", "ip")

    # Vinf grid

    plot_density_profile(vinf_grid.copy(), incls, vinf_labels, "ne_Vinf", "ne")
    plot_density_profile(vinf_grid.copy(), incls, vinf_labels, "rho_Vinf", "rho")
    plot_density_profile(vinf_grid.copy(), incls, vinf_labels, "ip_Vinf", "ip")

    # Plot the "best" lines

    plot_density_profile(best_lines_grid.copy(), incls, best_lines_labels, "ne_zBestLines", "ne")
    plot_density_profile(best_lines_grid.copy(), incls, best_lines_labels, "rho_zBestLines", "rho")
    plot_density_profile(best_lines_grid.copy(), incls, best_lines_labels, "ip_zBestLines", "ip")

    for i in incls:

        # Mbh grid

        plot_density_profile_inclination(mbh_grid.copy(), i, mbh_labels, "ne_Mbh", "ne")
        plot_density_profile_inclination(mbh_grid.copy(), i, mbh_labels, "rho_Mbh", "rho")
        plot_density_profile_inclination(mbh_grid.copy(), i, mbh_labels, "ip_Mbh", "ip")

        plot_density_profile_inclination(mbh_fixed_grid.copy(), i, mbh_labels, "ne_Mbh_fixed", "ne")
        plot_density_profile_inclination(mbh_fixed_grid.copy(), i, mbh_labels, "rho_Mbh_fixed", "rho")
        plot_density_profile_inclination(mbh_fixed_grid.copy(), i, mbh_labels, "ip_Mbh_fixed", "ip")

        # Rmin grid

        plot_density_profile_inclination(rmin_grid.copy(), i, rmin_labels, "ne_Rmin", "ne")
        plot_density_profile_inclination(rmin_grid.copy(), i, rmin_labels, "rho_Rmin", "rho")
        plot_density_profile_inclination(rmin_grid.copy(), i, rmin_labels, "ip_Rmin", "ip")

        # Vinf grid

        plot_density_profile_inclination(vinf_grid.copy(), i, vinf_labels, "ne_Vinf", "ne")
        plot_density_profile_inclination(vinf_grid.copy(), i, vinf_labels, "rho_Vinf", "rho")
        plot_density_profile_inclination(vinf_grid.copy(), i, rmin_labels, "ip_Vinf", "ip")

        # Plot the "best" lines

        plot_density_profile_inclination(best_lines_grid.copy(), i, best_lines_labels, "ne_zBestLines", "ne")
        plot_density_profile_inclination(best_lines_grid.copy(), i, best_lines_labels, "rho_zBestLines", "rho")
        plot_density_profile_inclination(best_lines_grid.copy(), i, best_lines_labels, "ip_zBestLines", "ip")

    return


if __name__ == "__main__":
    main(len(argv), argv)
