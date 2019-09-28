#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of the script is to create velocity plots for a Python simulation.
It's currently in a WIP state, and should likely be incorporated into
py_plot.py or another similar plot.
"""

import numpy as np
import py_plot_util as ppu
import py_plot as pp
from matplotlib import pyplot as plt
from typing import Union


def plot_velocity_magnitude(root: str, projection: str, ret_mag: bool = False) -> Union[None, np.ndarray]:
    """
    Plot the velocity magnitude exclusively - I may want to do this for some reason?

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    projection: str
        The coordinate projection of the simulation
    ret_mag: bool, optional
        If True, the velocity magnitude will be returned instead of plotted

    Returns
    -------
    v: np.ndarray
        The magnitude of the velocity for each grid cell
    """

    aprojections = ["rectilinear", "polar"]
    if projection not in aprojections:
        print("Projection {} unknown. Allowed projections are {}".format(projection, aprojections))

    xx, yx, vx = ppu.extract_wind_var(root, "v_x", "wind", coord=projection)
    xy, yy, vy = ppu.extract_wind_var(root, "v_y", "wind", coord=projection)
    xz, yz, vz = ppu.extract_wind_var(root, "v_z", "wind", coord=projection)
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

    if ret_mag:
        return v

    if projection == "polar":
        pp.polar_wind_plot(xx, yx, v, 0, "|v|", "wind", (1, 1))
    else:
        fig, ax = plt.subplots(1, 1, figsize=(12, 8), squeeze=False)
        fig, ax = pp.rectilinear_wind_plot(fig, ax, xx, yx, v, 0, 0, "|v|", "wind")

    plt.savefig("{}.velocitymag.png".format(root))
    plt.close()

    return


def plot_velocity_maps(root: str, projection: str):
    """
    Plot the 3 direction components of the wind as well as the magnitude of the
    velocity.

    Parameters
    ----------
    root: str
        The root name of the Python simulation
    projection: str
        The coordinate projection of the simulation
    """

    aprojections = ["rectilinear", "polar"]
    if projection not in aprojections:
        print("Projection {} unknown. Allowed projections are {}".format(projection, aprojections))

    xx, yx, vx = ppu.extract_wind_var(root, "v_x", "wind", coord=projection)
    xy, yy, vy = ppu.extract_wind_var(root, "v_y", "wind", coord=projection)
    xz, yz, vz = ppu.extract_wind_var(root, "v_z", "wind", coord=projection)
    v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    d = [vx, vy, vz, v]
    n = ["vx", "vy", "vz", "|v|"]

    nrows = 2
    ncols = 2
    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 8), squeeze=False)

    iidx = 0
    for i in range(nrows):
        for j in range(ncols):
            if projection == "polar":
                pp.polar_wind_plot(xx, yx, d[iidx], iidx, n[iidx], "wind", (2, 2))
            else:
                fig, ax = pp.rectilinear_wind_plot(fig, ax, xx, yx, d[iidx], i, j, n[iidx], "wind")
            iidx += 1

    theta = np.deg2rad(80)
    xc = np.linspace(xx[0, 0], xx[-1, 0], 50)
    zx = xc * np.tan(theta)
    ax[0, 0].plot(xc, zx, label="80 degree")
    ax[0, 0].set_xlim(1, 1e18)
    ax[0, 0].set_ylim(1, 1e8)

    if projection == "rectilinear":
        fig.tight_layout(rect=[0.02, 0.02, 0.98, 0.98])
    plt.savefig("{}.velocitymaps.png".format(root))
    plt.show()

    return


if __name__ == "__main__":
    plot_velocity_maps("tde_agn", "rectilinear")
    # plot_velocity_maps("tde_spherical", "polar")
    # velocity_magnitude("tde_agn", "rectilinear")
    # velocity_magnitude("tde_spherical", "polar")
