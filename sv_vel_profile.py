#!/usr/bin/env python3


import numpy as np
from scipy.optimize import brentq
from astropy import constants


verbose = True


class SV93wind:
    """
    Create a SV93 style wind object
    """

    bigG = constants.G.cgs.value

    def __init__(self, m_cobj, mdot_wind, r_min, r_max, theta_min, theta_max,
                 accel_length, accel_exp, v_inf):
        self._v0 = 6e5   # 6 km/s in cm/s
        self._gamma = 1  # Place stream lines uniform between rmin-rmax
        self._mobj = m_cobj
        self._mdot = mdot_wind
        self._rmin = r_min
        self._rmax = r_max
        self._thetamin = theta_min
        self._thetamax = theta_max
        self._Rv = accel_length
        self._alpha = accel_exp
        self._vinf = v_inf

    def find_theta(self, r0):
        """
        Determine the angle at which the wind emerges from at a special radius
        r from the disk surface
        """

        x = ((r0 - self._rmin) / (self._rmax - self._rmin))**self._gamma

        if r0 <= self._rmin:
            theta = np.arctan(np.tan(self._thetamax * r0 / self._rmin))
        elif r0 >= self._rmax:
            theta = self._thetamax
        else:
            theta = self._thetamin + (self._thetamax - self._thetamin) * x

        return theta

    def r0_guess_func(self, r, x):
        """
        Note that r is a position along the disk
        """

        theta = self.find_theta(r)
        rho = np.sqrt(x[0]**2 + x[1]**2)
        rho_guess = r + np.tan(theta) * x[2]

        return rho_guess - rho  # We want to make this zero

    def find_r0(self, x):
        """
        Determine r0 for a point in the x, y plane
        """

        # If the vector is in the x-y plane, then this is simple
        if x[2] == 0:
            return np.sqrt(x[0]**2 + x[1]**2)

        # For when the vector is not solely in the x-y plane
        rho_min = self._rmin + x[2] * np.tan(self._thetamin)
        rho_max = self._rmax + x[2] * np.tan(self._thetamax)
        rho = np.sqrt(x[0]**2 + x[1]**2)

        if rho <= rho_min:
            return self._rmin * rho / rho_min
        elif rho >= rho_max:
            return self._rmax * rho - rho_max
        else:
            return brentq(self.r0_guess_func, self._rmin, self._rmax, args=(x))

    def poloidal_velocity(self, l, r0):
        """
        Calculate the polodial velocity for a polodial distance l along a wind
        stream line with fixed.
        """

        if l < 0:
            return self._v0

        tmp = (l / self._Rv)**self._alpha
        v_esc = np.sqrt (2 * self.bigG * self._mobj / r0)
        vl = self._v0 + (self._vinf * v_esc - self._v0) * (tmp / (tmp + 1))

        return vl

    def velocity(self, x):
        """
        Determine the 3d velocity vector in cartesian coordinates
        """

        r0 = self.find_r0(x)
        theta = self.find_theta(r0)

        r = np.sqrt(x[0]**2 + x[1]**2)
        l = np.sqrt((r - r0)**2 + x[2]**2)
        vl = self.poloidal_velocity(l, r0)

        v = np.zeros(3)
        v[0] = vl * np.sin(theta)
        if r > 0:
            v[1] = np.sqrt(self.bigG * self._mobj * r0) / r
        else:
            v[1] = 0
        v[2] = np.abs(vl * np.cos(theta))

        return v


def main():
    """
    Main controlling function
    """

    Rstar = 2.65e13
    Msun = constants.M_sun.cgs.value

    v_inf = 3
    accel_exp = 1.5
    m_cobj = 3e7 * Msun           # g
    mdot_wind = 6.305e23          # g/s
    r_min = 10 * Rstar             # cm
    r_max = 55 * Rstar            # cm
    theta_min = np.deg2rad(20)    # degrees
    theta_max = np.deg2rad(65)    # degrees
    accel_length = 5e16           # cm
    wind_radmax = 1e17            # cm

    # Create SV wind object
    sv = SV93wind(m_cobj, mdot_wind, r_min, r_max, theta_min, theta_max,
                  accel_length, accel_exp, v_inf)

    # Create a 3d cylindrical grid - note that there is no theta dependence as
    # axisymmetry is assumed
    n_resolution = 30
    rgrid = np.linspace(sv._rmin, wind_radmax, n_resolution)
    tgrid = np.zeros(n_resolution)
    zgrid = np.linspace(0, wind_radmax, n_resolution)
    r, theta, z = np.meshgrid(rgrid, tgrid, zgrid)

    x = np.zeros(3)
    for i in range(r.shape[0]):
        for j in range(theta.shape[0]):
            for k in range(theta.shape[0]):
                x[0] = r[i, j, k]
                x[1] = theta[i, j, k]
                x[2] = z[i, j, k]
                v = sv.velocity(x)
                print(v)

    return

if __name__ == "__main__":
    main()
