#!/usr/bin/env python3

import tde_util
import py_plot_util
import numpy as np
from matplotlib import pyplot as plt


SMOOTH = 11
VERBOSE = False
spec_files = py_plot_util.find_spec_files()
blag_spec = tde_util.iptf15af_spec(SMOOTH, VERBOSE)

fig, ax = plt.subplots(1, 1, figsize=(12,8))
Z = 0.07897
ax.semilogy(blag_spec[:, 0] / (Z + 1), blag_spec[:, 1], label="iPTF15af: Blagorodnova et al. 2019")
labels = ["AGN Style Wind at 75° inclination", "CV Style Wind at 62° inclination"]
j = 0
for file in spec_files:
    spec = py_plot_util.read_spec_file(file, " ")

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
    smoothflux = py_plot_util.smooth_1d_array(flux, SMOOTH, VERBOSE)
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
    OBSERVE_DIST = 1.079987153448e+27
    default_dist = 100 * 3.08567758128e18  # 100 pc
    flux_dist = smoothflux * (default_dist ** 2 / OBSERVE_DIST ** 2)
    z_wav = wavelength

    ax.semilogy(z_wav, flux_dist, label=labels[j])
    j+=1

ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=17)
ax.tick_params(labelsize=14)
ax.set_xlim(1150, 2500)
ax.set_ylim(1e-17, 1e-14)
ax.legend()
ax.set_title("Preliminary wind models", fontsize=17)
fig.tight_layout()
plt.savefig("poster_spec.png")
plt.show()