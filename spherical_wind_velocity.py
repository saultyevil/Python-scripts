#!/usr/bin/env python3

"""
Plot various velocity power laws for the Castor & Lamers prescription for
a stellar wind.
"""


import numpy as np
from constants import *
from matplotlib import pyplot as plt


def cl_velocity(r, v0, vinf, rstar, beta):
    """
    Calculate the velocity for a 1D spherical wind using the velocity law as
    given in Castor & Lamers 1979.

    Parameters
    ----------
    r: float
        The distance to calculate the velocity at
    v0: float
        The initial velocity at the bottom of the spherical wind
    vinf: float
        The terminal velocity at the edge of the wind
    rstar: float
        The radius of the central source (star)
    beta: float
        The acceleration exponents. Controls how rapidly the velocity
        is accelerated to the terminal velocity vinf.

    Returns
    -------
    v: float
        The velocity at a point r for a 1D spherical wind
    """

    return v0 + (vinf - v0) * (1 - rstar / r) ** beta


def plot_power_law():
    """
    Create a plot for various different values of beta for a 1D spherical wind.
    """

    plt.figure(figsize=(12, 8))

    betas = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 12.4, 15.0, 20.0]

    # vinf = 1
    mobj = 3e7 * MSOL
    rin = 2.65e13
    rout = 1e17
    v0 = 2.73e6
    vinf = 1.73e9

    nres = int(1e5)
    rgrid = np.linspace(rin, rout, nres)
    for beta in betas:
        v_r = cl_velocity(rgrid, v0, vinf, rin, beta)
        plt.semilogx(rgrid, v_r / vinf, label=r"$\beta$ = " + str(beta))

    # plt.xlim(1.0, rout)
    plt.ylim(0.0, 1.0)
    plt.xlabel(r"r/$R_{*}$", fontsize=14)
    plt.ylabel(r"v(r)/$v_{\infty}$", fontsize=14)
    plt.legend(loc="lower right")
    plt.title(r"$v_{\infty}$ = " + "{:2.1e} cm/s".format(vinf) + " = {:1.2f}c".format(vinf / C))

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return


if __name__ == "__main__":
    plot_power_law()
