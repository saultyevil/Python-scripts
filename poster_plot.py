#!/usr/bin/env python3

import py_util
import numpy as np
from matplotlib import  pyplot as plt


SMOOTH = 11
VERBOSE = False
spec_files = py_util.find_spec_files()
blag_spec = py_util.get_blagordnova_spec(SMOOTH, VERBOSE)

fig, ax = plt.subplots(1, 1, figsize=(12,8))
ax.semilogy(blag_spec[:, 0], blag_spec[:, 1], label="Blagorodnova et al. 2018 (in prep.)")
labels = ["AGN Style Wind at 75° inclination", "CV Style Wind at 62° inclination"]
j = 0
for file in spec_files:
    spec = py_util.read_spec_file(file, " ")

    idx = 0
    for i in range(spec.shape[1]):
        if spec[0, i].isdigit() and float(spec[0, i]) == 62:
            spec[0, i] = 62
            idx = i
            break
        elif spec[0, i].isdigit() and float(spec[0, i]) == 75:
            spec[0, i] = 75
            idx = i
            break
    flux = np.array(spec[1:, idx], dtype=float)
    smoothflux = py_util.smooth_spectra(flux, SMOOTH, VERBOSE)
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
    OBSERVE_DIST = 1.079987153448e+27
    Z = 0.07897
    default_dist = 100 * 3.08567758128e18  # 100 pc
    flux_dist = smoothflux * (default_dist ** 2 / OBSERVE_DIST ** 2)
    z_wav = wavelength * (Z + 1)

    ax.semilogy(z_wav, flux_dist, label=labels[j])
    j+=1

ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)",
              fontsize=17)
ax.tick_params(labelsize=14)
ax.set_xlim(1150, 2700)
ax.set_ylim(1e-17, 1e-14)
ax.legend(fontsize=14)
ax.set_title("Preliminary wind models", fontsize=17)
fig.tight_layout()
plt.savefig("poster_spec.png")
plt.show()