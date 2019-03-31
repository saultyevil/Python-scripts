#!/usr/bin/env python3

import py_plot_util
import numpy as np
from matplotlib import pyplot as plt

SMOOTH = 11
VERBOSE = False
OBSERVE_DIST = 1.079987153448e+27      # 350 Mpc
DEFAULT_DIST = 100 * 3.08567758128e18  # 100 pc

fig, ax = plt.subplots(1, 1, figsize=(12, 8))

# Read in and normalise the observed spectra
cenkspec = py_plot_util.get_ASSASN14li_spec(SMOOTH, VERBOSE)
blag_spec = py_plot_util.get_iPTF15af_spec(SMOOTH, VERBOSE)
cenk_max = cenkspec[:, 1].max()
cenk_norm = cenkspec[:, 1] / cenk_max
blag_max = blag_spec[:, 1].max()
blag_norm = blag_spec[:, 1] / blag_max

# Plot the observed spectra
ax.semilogy(blag_spec[:, 0] / (0.07897 + 1), blag_norm, label="iPTF15af: Blagordnova et al. (2019)")
ax.semilogy(cenkspec[:, 0] / (0.02058 + 1), cenk_norm * 20, label="ASASSN-14li: Cenko et al. (2016)")

# Load in the Python model - assume this is being called in the model directory
# spec_file = "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11/tde.spec"
spec_file = py_plot_util.find_spec_files()[0]
print(spec_file)
spec = py_plot_util.read_spec_file(spec_file, " ")

# Horrid hacky code for easier indexing :^) please ignore
angles = [28, 45, 62, 80]
idxs = []
for angle in angles:
    for i in range(spec.shape[1]):
        if spec[0, i].isdigit() and float(spec[0, i]) == angle:
            spec[0, i] = float(angle)
            idxs.append(i)

# Pull out the correct flux for the angle and scale appropriately
labels = ["Matom model 28°", "Matom model 45°", "Matom model 62°", "Matom model 80°"]  # the python interpreter may HATE this
offset = [0.1, 1, 10, 25]
for i in range(len(angles)):
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)
    flux = np.array(spec[1:, idxs[i]], dtype=float)
    smoothflux = py_plot_util.smooth_flux(flux, SMOOTH, VERBOSE)
    flux_dist = smoothflux * (DEFAULT_DIST ** 2 / OBSERVE_DIST ** 2)  # + 6e-16 *offset[i])
    flux_max = flux_dist.max()
    flux_norm = flux_dist / flux_max
    ax.semilogy(wavelength, flux_norm * offset[i], label=labels[i])

# Plot attributes and etc
ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
ax.set_ylabel(r"Normalised Flux $F_{\lambda}$ + const", fontsize=17)
ax.tick_params(labelsize=14)
ax.set_xlim(1150, 2800)
ax.legend(fontsize=14)
ax.set_title("Preliminary wind model", fontsize=17)

fig.tight_layout()
plt.savefig("tde_models.png")
plt.show()
