#!/usr/bin/env python3

import sys
import py_util
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import convolve, boxcar

try:
    angle = sys.argv[1]
except:
    print("No angle provided ;_;")
    sys.exit(1)
aaaa = py_util.find_spec_files()
specfile = aaaa[0]
print("Loading spec file {}".format(specfile))
spec = py_util.read_spec_file(specfile)
print(spec)
allowed = py_util.check_viewing_angle(angle, spec)
if not allowed:
    print("Invalid angle {}".format(angle))
    sys.exit(1)
idx = 0
for i in range(spec.shape[1]):
    if spec[0, i].isdigit() and float(spec[0, i]) == float(angle):
        spec[0, i] = float(spec[0, i])
        idx = i
        break
default_dist = 100 * 3.08567758128e18
tde_dist = 1.079987153448e+27
flux = np.array(spec[1:, idx], dtype=float) * (default_dist**2/tde_dist**2)
flux = np.reshape(flux, (len(flux), ))
smooth = 100
smoothflux = convolve(flux, boxcar(smooth) / float(smooth), mode="same")
print("\n\n---------------------------\n", smoothflux)
print("Peak flux = {} ergs/s".format(np.max(smoothflux)))
peakidx, = np.where(smoothflux == np.max(smoothflux))
wavpeak = float(spec[1+peakidx[0], 1]) * 1e-10
print("Peak lambda = {} m".format(wavpeak))
Tbb = 2.897772917e-3 / wavpeak
print("Tbb = {} K".format(Tbb))
plt.semilogy(np.array(spec[1:, 1],dtype=float), smoothflux)
plt.xlabel("lambda")
plt.ylabel("flux")
plt.savefig("bbpeak.pdf")
