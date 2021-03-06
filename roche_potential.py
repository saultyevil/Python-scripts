#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Create a plot of the Roche potential for a simple binary system.
"""

import numpy as np
from constants import PI, MSOL, G
from matplotlib import pyplot as plt


def calc_orbital_separation(
    period: float, m_star: float, m_seco: float
) -> float:
    """
    Calculate the oribtal separation a for a binary system given the orbital
    period and masses - uses Kepler's law.

    The period should be given in hours.

    Parameters
    ----------
    period: float
      The period of the binary system
    m_star: float
      The mass of the primary (most compact) object in the system
    m_seco: float
      The mass of the secondary object in the system

    Returns
    -------
    t: float
      The orbitical separation of the binary system in m.
    """

    period *= 3600
    tmp = period ** 2 * (m_star + m_seco) / (16 * PI ** 2)

    return tmp ** (1 / 3)


def roche_potential(
    m1: float, m2: float, a: float, omega: float, x: np.ndarray, xc: float, y: np.ndarray, z: np.ndarray
) -> np.ndarray:
    """
    Calculate the Roche potential.

    Parameters
    ----------
    m1: float
      The mass of primary object in the system.
    m2: float
      The mass of the secondary object in the system.
    a: float
      The orbital separation of the system.
    omega: float
      The
    x: np.array[float]
      The x-range to calculate the Roche potential over.
    xc: float
      Some value
    y: np.array[float]
      The y-range to calculate the Roche potential over.
    z: np.array[float]
      The z-range to calculate the Roche potential over.

    Returns
    -------
    phi: np.array[float]
      The Roche potential at a bunch of points r (x**2+y**2+z**2) etc
    """

    r1 = x ** 2 + y ** 2 + z ** 2
    r2 = (x - a) ** 2 + y ** 2 + z ** 2

    phi_1 = - G * m1 / r1
    phi_2 = - G * m2 / r2
    phi_3 = - 0.5 * omega ** 2 * ((x - xc) ** 2 + y ** 2)
    phi = phi_1 + phi_2 + phi_3

    return phi


def dimensionless_roche_potential(
    r1: float, r2: float, q: float, x: float, y: foat
) -> float:
    """
    Calculate the shape of the Roche potential independent of the size of the
    binary system.
    """

    phi_1 = 2 / (1 + q) / r1
    phi_2 = 2 * q / (1 + q) / r2
    phi_3 = (x + q / (1 + q)) ** 2 - y ** 2
    phi = phi_1 + phi_2 + phi_3

    return phi


def main():
    """
    Main function for plotting Roche potential as a function of orbial
    separation.
    """

    # Standard CV
    period = 5.57
    m_star = 0.80 * MSOL
    m_seco = 0.60 * MSOL

    # TDE
    # period = 131
    # m_star = 3e7 * MSOL
    # m_seco = 1.6 * MSOL

    # Calculate some system parameters
    q = m_star / m_seco
    a = calc_orbital_separation(period, m_star, m_seco)
    omega = G * (m_star + m_seco) / a ** 3
    xc = a * m_seco / (m_star + m_seco)

    # Now we can calculate the Roche potential for a given xrange
    l_xlim, u_xlim = 3 * np.array([-a, a])
    x_range = np.linspace(l_xlim, u_xlim, 500)
    y_range = 0
    z_range = 0
    phi = roche_potential(m_star, m_seco, a, omega, x_range, xc, y_range, z_range)

    plt.plot(x_range / a, phi)
    # plt.axhline(np.max(phi), 0, 1, linestyle="--", linewidth=1, color="k")
    plt.xlabel(r"$x / a$")
    plt.ylabel(r"Roche Potential $\Phi (x)$")
    plt.ylim(-3, 0)
    plt.xlim(l_xlim / a, u_xlim / a)
    plt.show()

    return


if __name__ == "__main__":
    main()
