#!/usr/bin/env python3

import py_plot
import numpy as np
from astropy import constants as const
from astropy.modeling import blackbody
from astropy import units as u
from matplotlib import pyplot as plt

def B_lamda(T_kelvin, lamda_cm):
    B1 = (2 * const.h.cgs.value * const.c.cgs.value**2) / lamda_cm ** 5
    hnu = (const.h.cgs.value * const.c.cgs.value) / (lamda_cm * T_kelvin *
                                                     const.k_B.cgs.value)
    B2 = 1 / (np.exp(hnu) - 1)

    return B1 * B2

blagspec = py_plot.load_blag_spec()
plt.semilogy(blagspec[:, 0], blagspec[:, 1], label="Blagorodnova et al. 2018")

# bbflux = B_lamda(T, lambda_range * 1e-8)*(4*np.pi*tde_dist**2)
# print(1/(4*np.pi*tde_dist**2), bbflux, bbflux* (1/(4*np.pi*tde_dist**2)),
#       bbflux* (1/(4*np.pi*tde_dist)))


# astropy monochromatic bb flux
#    Returns
#    -------
#    flux : `~astropy.units.Quantity`
#        Blackbody monochromatic flux in
#        :math:`erg \\; cm^{-2} s^{-1} \\mathring{A}^{-1} sr^{-1}`.

lambda_range = np.linspace(1100, 3300, 100)  # Angstrom
tde_dist = 1.079987153448e+27 # 358 Mpc in cm
T = 43300  # kelvin
bbflux = blackbody.blackbody_lambda(lambda_range, T) / (4*np.pi*tde_dist**2)

plt.semilogy(lambda_range, bbflux, label="BB = {} K".format(T))
plt.xlabel("Wavelength, $\AA$")
plt.ylabel("$F_{\lambda}$ (ergs/s/cm$^{2}$/$\AA$)")
plt.legend()
plt.savefig("bb_fit.pdf")
plt.show()
