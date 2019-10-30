#!/usr/bin/env python3

"""
Plot the velocity law for a Schlosman and Vitello wind for a CV disk wind.
"""


import ss_disk
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
                 accel_length, accel_exp, v_inf, gamma, v0=6e5):
        """
        Initialise the SV wind object with the parameters given. Note that the
        default value for v0 is 6e5 cm/s.
        """

        self.v0 = v0
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

    def find_theta(self, r0):
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

    def r0_guess_func(self, r, x):
        """
        Note that r is a position along the disk
        """

        theta = self.find_theta(r)
        rho = np.sqrt(x[0] ** 2 + x[1] ** 2)
        rho_guess = r + np.tan(theta) * x[2]

        return rho_guess - rho  # We want to make this zero

    def find_r0(self, x):
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
            return brentq(self.r0_guess_func, self.rmin, self.rmax, args=(x))

    def vesc(self, r0):
        """
        Calculate the escape velocity at a point r0
        """

        return np.sqrt(2 * G * self.mobj / r0)

    def poloidal_velocity(self, l, r0):
        """
        Calculate the polodial velocity for a polodial distance l along a wind
        stream line with fixed.
        """

        if l < 0:
            return self.v0

        tmp = (l / self.Rv) ** self.alpha
        v_term = self.vinf * self.vesc(r0)
        vl = self.v0 + (v_term - self.v0) * (tmp / (tmp + 1))

        return vl

    def velocity_vector(self, x):
        """
        Determine the 3d velocity vector in cartesian coordinates
        """

        r0 = self.find_r0(x)

        if self.rmin > r0 or r0 > self.rmax:
            print("r0 outside rmin or rmax of wind")
            print("r0 = ", r0)
            exit(1)

        theta = self.find_theta(r0)
        # if self.thetamin < theta < self.thetamax:
        #     print("theta cannot be smaller than thetamin or larger than thetamax")
        #     print(theta)
        #     exit(1)

        r = np.sqrt(x[0] ** 2 + x[1] ** 2)
        l = np.sqrt((r - r0) ** 2 + x[2] ** 2)
        vl = self.poloidal_velocity(l, r0)

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
    r_star = 2.65e13                # cm
    r_min = 1 * r_star              # cm
    r_max = 30 * r_star             # cm
    theta_min = 70                  # degrees
    theta_max = 82                  # degrees
    accel_length = 5e16             # cm

    # Create the SV wind object for the parameters given above
    wind = SV93wind(m_cobj, mdot_wind, r_min, r_max, theta_min, theta_max, accel_length, accel_exp, v_inf, gamma)

    # This next part of the code will create a plot for different values of alpha
    # which controls how quickly the wind is accelerated. To do this, we will assume
    # the streamline originates from r_min where r0 is the launch radius of the
    # streamline
    n_resolution = 500
    wind_alphas = [0.1, 0.2, 0.4, 0.8, 1.0, 2.0, 4.0, 8.0, 10.0, 15.0]
    # wind_alphas = [0.5, 1.0, 2.0]

    print("Plotting SV93 polodial velocity power law for:")
    print(wind_alphas)

    plt.figure(figsize=(12, 8))

    r0 = 1.5 * wind.rmin
    if r0 < wind.rmin:
        print("r0 < wind.rmin")
        exit(1)
    if r0 > wind.rmax:
        print("r0 > wind.rmax")
        exit(1)

    v0_sound_speed = 0
    if v0_sound_speed:
        print("SV93 v0 set to {} sound speed at {}rmin".format(v0_sound_speed, r0 / wind.rmin))
        teff = ss_disk.t_eff(r0, wind.mobj, wind.mdot, r_star)
        wind.v0 = 1e6 * np.sqrt(teff / 1e4)  # Taken from Frank, King & Raine 1985
        wind.v0 *= v0_sound_speed
        print("teff = {:e}".format(teff))
    print("v0 = {:e}".format(wind.v0))
    rad_max = 5e17 # 2 * wind.Rv
    l = np.linspace(wind.rmin, rad_max, n_resolution)
    for al in wind_alphas:
        wind.alpha = al
        vl = np.zeros(n_resolution)
        for i in range(n_resolution):  # TODO: if I was clever earlier, this wouldn't need to be a loop
            vl[i] = wind.poloidal_velocity(l[i], r0)
        vinf = wind.vinf * wind.vesc(r0)
        plt.plot(l / wind.Rv, vl / vinf, label=r"$\alpha$ = {}".format(al))
    print("vinf = {:.3f} c".format(vinf / C))

    plt.xlabel(r"$l$/$R_{v}$", fontsize=15)
    plt.ylabel(r"$v_{l}$/$v_{\infty}$", fontsize=15)
    plt.title(r"$r_{0} = $" + str(r0 / wind.rmin) + r"$r_{min}$")
    plt.xlim(0, rad_max / wind.Rv)
    plt.ylim(0, 1)
    plt.legend()
    plt.show()

    return


if __name__ == "__main__":
    plot_power_law()
