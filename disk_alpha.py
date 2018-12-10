#!/usr/bin/env python3

import numpy as np
from scipy.constants import pi
from matplotlib import pyplot as plt


def teff(R):
    mobj = 1.988435e40    # g: 0.8 msun
    mdot = 1.261057e23     # g/s: 1e-10 msun/yr
    rstar = 9e12       # cm
    sigma = 5.67e-5   # erg cm^-2 K^-4 s^-1
    G = 6.67e-8        # cm^3 g^-1 s^-2
    return ((3*G*mobj*mdot)/(8*pi*R**3*sigma)*(1-(rstar/R)**0.5))**0.25
    
Rstar = 9e12
N = 100
Rs = np.linspace(Rstar, N*Rstar, 100)
ticks = []
ticklabels =[]    
for i in range(1, N+1 ):
    st = "{}".format(i)+"$R_{*}$"
    ticks.append(i * Rstar)
    ticklabels.append(st)
plt.plot(Rs, teff(Rs))
plt.xticks(ticks[::10], ticklabels[::10])
plt.ylabel("T(R)")
plt.xlabel("R")
plt.show()
