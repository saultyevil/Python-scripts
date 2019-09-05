#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from typing import Union


def flared_disk_z(r: Union[np.ndarray, float], rout: float, zmax: float, zscale: float) -> Union[np.ndarray, float]:
    """
    Returns the vertical height z of a flared disc given the parameters defined
    below.

    Parameters
    ----------
    r: np.ndarray, float

    rout: float

    zmax: float

    zscale: float

    Returns
    -------
    z: np.ndarray, float
        The vertical height of the flared disc as a position(s) r.
    """

    return rout * zmax * (r / rout) ** zscale


def doit(rmin: float, rmax: float, npoints: int = 500):
    """

    Returns
    -------

    """

    rout = 1e15
    zmax = 1e-2
    zscale = 1.5

    r = np.linspace(rmin, rmax, npoints)
    z = flared_disk_z(r, rout, zmax, zscale)

    title = "$R_{out}$ = " + str(rout) + " $Z_{max}$ = " + str(zmax * rout) + " $Z_{scale}$ = " + str(zscale)
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    ax.semilogy(r, z)
    ax.set_xlabel(r"Disk Radius, $R$")
    ax.set_ylabel(r"Disk Height, $z(R)$")
    ax.set_title(title)

    fig.tight_layout()
    plt.show()

    return


if __name__ == "__main__":
    doit(2.65e13, 1e15)