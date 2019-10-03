#!/usr/bin/env python3

"""
Extinction along the line of sight, AV = 0.093 mag
"""

from platform import system
from PyPython import SpectrumUtils
import numpy as np
from matplotlib import pyplot as plt


def mag2flux(mag, zero):
    return zero * 10**(-mag / 2.5)


smooth = 25
verbose = False
Z = 0.07897

# load in cv spec
spec_file = "/Users/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/tde_macro/tde.spec"
spec = SpectrumUtils.read_spec(spec_file, numpy=True)
wavelength = np.array(spec[1:, 1], dtype=float)

# flux indexs should be 11 (45 deg) or 12 (62 deg) -- rescale to iPTF15af distance
flux45 = SpectrumUtils.smooth_spectrum(np.array(spec[1:, 11], dtype=float), smooth)
flux45 *= 3.08567758128e20**2 / 1.079987153448e+27**2
flux62 = SpectrumUtils.smooth_spectrum(np.array(spec[1:, 12], dtype=float), smooth)
flux62 *= 3.08567758128e20**2 / 1.079987153448e+27**2

# load iPTF15af and put into rest frame
spec_dir = ""
sys = system()
if sys == "Linux":
    spec_dir = "/home/saultyevil/"
elif sys == "Darwin":
    spec_dir = "/Users/saultyevil/"
else:
    exit(1)
spec_dir += "PySims/tde/observed_spec/Blagorodnova_iPTF15af.dat"
blag = np.loadtxt(spec_dir)
blag[:, 1] = SpectrumUtils.smooth_spectrum(blag[:, 1], smooth)
blag[:, 0] /= (Z + 1)

mags = np.array([
    19.48,   # UVW1
    19.42,   # UVM2
    19.16,   # UVW2
    # 19.27,   # U
    # 18.81,   # B
    # 17.91,   # V
    20.09,   # g
    20.39,   # r
    20.31,   # i
])

# Taken from http://svo2.cab.inta-csic.es/svo/theory/fps/
zero_points = np.array([
    1.624e-8,
    2.193e-8,
    2.640e-8,
    # 8.880e-9,
    # 5.810e-9,
    # 3.730e-9,
    5.055e-9,
    2.904e-9,
    7.738e-9
])

# Taken from http://svo2.cab.inta-csic.es/svo/theory/fps/
central_waves = np.array([
    2604.57,
    2246.00,
    1941.22,
    # 3463.14,
    # 4371.22,
    # 5441.20,
    4700.33,
    6174.48,
    3684.29
])

labels = [
    "Swift (+52d): UVW1",
    "Swift (+52d): UWM2",
    "Swift (+52d): UVW2",
    # "Swift (+52d): U",
    # "Swift (+52d): B",
    # "Swift (+52d): V",
    "P60 (+53.9d): g",
    "P60 (+53.9d): r",
    "LasCumbres (+49.5d): I"
]

# convert magnitudes to flux
fluxs = mag2flux(mags, zero_points)
for i in range(len(fluxs)):
    print(labels[i], central_waves[i], fluxs[i])

# finally plot
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.semilogy(wavelength, flux45, label="Python model: 45 degrees")
ax.semilogy(wavelength, flux62, label="Python model: 62 degrees")
ax.semilogy(blag[:, 0], blag[:, 1], label="iPTF15af: Blagorodnova et al. (2019)")

for i in range(len(fluxs)):
    leg_label = labels[i] + r" ({:4.0f} $\AA$)".format(central_waves[i])
    ax.semilogy(central_waves[i], fluxs[i], "d", label=leg_label)

ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
ax.set_xlabel(r"Wavelength ($\AA$)")
ax.set_xlim(500, 10000)
ax.legend()
fig.tight_layout()
plt.savefig("/Users/saultyevil/Dropbox/iIPTF15af_python_sed_photometry.png")
plt.close()
