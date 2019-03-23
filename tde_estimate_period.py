#!/usr/bin/env python3

import numpy as np
from astropy import constants as c

G = c.G.value
Msun = c.M_sun.value
Rsun = c.R_sun.value


def KeplerPeriod(M, a3):
    return np.sqrt((4 * np.pi **2 * a3) / (G * M))


def DisruptionRadius(Rstar, Mbh, Mstar):
    return Rstar * (Mbh / Mstar) ** (1.0 / 3.0)


class BlackHole:
    def __init__(self, M):
        self._M = M * Msun


class SecondaryStar:
    def __init__(self, M, R):
        self._M = M * Msun
        self._R = R * Rsun

bh = BlackHole(3e7)
star = SecondaryStar(1.6, 1)
Rd = DisruptionRadius(star._R, bh._M, star._M)
period_seconds = KeplerPeriod(bh._M, Rd ** 3)
period_hours = period_seconds / 60
print("Disruption radius = {:e} cm".format(Rd * 100))
print("Period = {:.2f} hrs".format(period_hours))
