#!/usr/bin/env python3

import numpy as np
from PyPython.Constants import *
from matplotlib import pyplot as plt
from typing import Union


def alpha_disc_effective_temperature(
    ri: Union[np.ndarray, float], rstar: float, mbh: float, mdot: float
) -> np.ndarray:
    """Standard alpha-disc effective temperature profile."""

    teff_4 = (3 * G * mbh * mdot) / (8 * PI * ri ** 3 * STEFAN_BOLTZMANN)
    teff_4 *= 1 - (rstar / ri) ** 0.5
    return teff_4 ** 0.25


def eddington_critical_disc_effective_temperature(
    ri: Union[np.ndarray, float], mbh: float, mdot: float, ledd: float, rg: float, risco: float
):
    """The effective temperature profile from Strubbe and Quataert 2009."""

    fnt = 1 - np.sqrt(risco / ri)
    teff_4 = (3 * G * mbh * mdot * fnt) / (8 * PI * ri ** 3 * STEFAN_BOLTZMANN)
    teff_4 *= (0.5 + (0.25 + 6 * fnt * (mdot * C ** 2 / ledd) ** 2 * (ri / rg) ** -2) ** 0.5) ** -1
    return teff_4 ** 0.25


mbh = 10 ** 6 * MSOL
mdot = 2.2e-2 * MSOL_PER_YEAR
rmin = 8.86e11
rmax = 1e15
rg = G * mbh / C ** 2
risco = 6 * rg
ledd = 4 * PI * G * mbh * C * MPROT / THOMPSON
radii = np.logspace(np.log10(rmin), np.log10(rmax), 3001)

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(
    radii, alpha_disc_effective_temperature(radii, rmin, mbh, mdot), label=r"$\alpha$-disc", linewidth=3
)
ax.plot(
    radii, eddington_critical_disc_effective_temperature(radii, mbh, mdot, ledd, rg, risco),
    label=r"Eddington critical disc", linewidth=3
)
ax.set_xlabel(r"Radius $R$")
ax.set_ylabel(r"Effective Temperature $T_{eff}$")
ax.legend()
ax.set_xscale("log")
ax.set_yscale("log")
plt.show()

output = np.zeros((3001, 2))
output[:, 0] = radii
output[:, 1] = eddington_critical_disc_effective_temperature(radii, mbh, mdot, ledd, rg, risco)
print(output)
np.savetxt("disc_temperature_profile.txt", output)
