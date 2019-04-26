#!/usr/bin/env python3


import numpy as np
from consts import *
from scipy.optimize import brentq
from matplotlib import pyplot as plt

verbose = False


class SV93wind:
    """
    Create a SV93 style wind object
    """

    def __init__(self, m_cobj, mdot_wind, r_min, r_max, theta_min, theta_max,
                 accel_length, accel_exp, v_inf, gamma):
        """
        Initialise the SV wind object with the parameters given.
        """

        self.v0 = 6e5  # 6 km/s in cm/s
        self.gamma = gamma
        self.mobj = m_cobj * MSOL
        self.mdot = mdot_wind * MSOL_PER_YEAR
        self.rmin = r_min
        self.rmax = r_max
        self.thetamin = np.deg2rad(theta_min)
        self.thetamax = np.deg2rad(theta_max)
        self.Rv = accel_length
        self.alpha = accel_exp
        self.vinf = v_inf

    def _find_theta(self, r0):
        """
        Determine the angle at which the wind emerges from at a special radius
        r from the disk surface
        """

        x = ((r0 - self.rmin) / (self.rmax - self.rmin)) ** self.gamma

        if r0 <= self.rmin:
            theta = np.arctan(np.tan(self.thetamax * r0 / self.rmin))
        elif r0 >= self.rmax:
            theta = self.thetamax
        else:
            theta = self.thetamin + (self.thetamax - self.thetamin) * x

        return theta

    def _r0_guess_func(self, r, x):
        """
        Note that r is a position along the disk
        """

        theta = self._find_theta(r)
        rho = np.sqrt(x[0] ** 2 + x[1] ** 2)
        rho_guess = r + np.tan(theta) * x[2]

        return rho_guess - rho  # We want to make this zero

    def _find_r0(self, x):
        """
        Determine r0 for a point in the x, y plane
        """

        # If the vector is in the x-y plane, then this is simple
        if x[2] == 0:
            return np.sqrt(x[0] ** 2 + x[1] ** 2)

        # For when the vector is not solely in the x-y plane
        rho_min = self.rmin + x[2] * np.tan(self.thetamin)
        rho_max = self.rmax + x[2] * np.tan(self.thetamax)
        rho = np.sqrt(x[0] ** 2 + x[1] ** 2)

        if rho <= rho_min:
            return self.rmin * rho / rho_min
        elif rho >= rho_max:
            return self.rmax * rho - rho_max
        else:
            return brentq(self._r0_guess_func, self.rmin, self.rmax, args=(x))

    def _vesc(self, r0):
        """
        Calculate the escape velocity at a point r0
        """

        return np.sqrt(2 * G * self.mobj / r0)

    def _poloidal_velocity(self, l, r0):
        """
        Calculate the polodial velocity for a polodial distance l along a wind
        stream line with fixed.
        """

        if l < 0:
            return self.v0

        tmp = (l / self.Rv) ** self.alpha
        v_term = self.vinf * self._vesc(r0)
        vl = self.v0 + (v_term - self.v0) * (tmp / (tmp + 1))

        return vl

    def _velocity_vector(self, x):
        """
        Determine the 3d velocity vector in cartesian coordinates
        """

        r0 = self._find_r0(x)
        theta = self._find_theta(r0)

        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        l = np.sqrt((r - r0) ** 2 + x[2] ** 2)
        vl = self._poloidal_velocity(l, r0)

        v = np.zeros(3)
        v[0] = vl * np.sin(theta)
        if r > 0:
            v[1] = np.sqrt(G * self.mobj * r0) / r
        else:
            v[1] = 0
        v[2] = np.abs(vl * np.cos(theta))

        return v


def plot_power_law():
    """
    Create a plot of the SV power law for various different values.
    """

    # Parameters for the SV wind system, see:
    #  - http://adsabs.harvard.edu/abs/1993ApJ...409..372S
    #  - http://adsabs.harvard.edu/abs/2014A%26A...561A..14K

    v_inf = 1
    gamma = 1
    accel_exp = 1.0
    m_cobj = 3e7                    # Msol
    mdot_wind = 2e-2                # Msol/yr
    Rstar = 2.65e13                 # cm
    r_min = 10 * Rstar              # cm
    r_max = 55 * Rstar              # cm
    theta_min = 70                  # degrees
    theta_max = 82                  # degrees
    accel_length = 5e16             # cm

    # Create the SV wind object for the parameters given above
    wind = SV93wind(m_cobj, mdot_wind, r_min, r_max, theta_min, theta_max, accel_length, accel_exp, v_inf, gamma)

    # This next part of the code will create a plot for different values of alpha
    # which controls how quickly the wind is accelerated. To do this, we will assume
    # the stream line originates from r_min
    n_resolution = 500
    wind_alphas = [0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4.0, 8.0, 10.0, 15.0]

    plt.figure(figsize=(12, 8))

    for al in wind_alphas:
        wind.alpha = al
        r0 = wind.rmin
        l = np.linspace(wind.rmin, 2 * wind.Rv, n_resolution)
        vl = np.zeros(n_resolution)
        for i in range(n_resolution):
            vl[i] = wind._poloidal_velocity(l[i], r0)
        vinf = wind.vinf * wind._vesc(r0)
        plt.plot(l / wind.Rv, vl / vinf, label=r"$\alpha$ = {}".format(al))

    plt.xlabel(r"$l$/$R_{v}$", fontsize=15)
    plt.ylabel(r"$v_{l}$/$v_{\infty}$", fontsize=15)
    plt.xlim(0, 2)
    plt.legend()
    plt.show()

    return


if __name__ == "__main__":
    plot_power_law()
