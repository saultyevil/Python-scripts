#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import sys
from platform import system
import numpy as np
from matplotlib import pyplot as plt
from PyPython import SpectrumUtils

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

from consts import *
import tde_util as tu

SMOOTH = 5
VERBOSITY = False

iptf = tu.iptf15af_spec(SMOOTH, VERBOSITY)

KMS = 1e5

vel_blueshift = 5000 * KMS  # cm/s
f0 = C / (1550 * ANGSTROM)
f = (1 + vel_blueshift / C) * f0
wl = 2.98e8 / f / 1e-10
print(wl)

z = 0.07897
fig, ax = plt.subplots(1, 1, figsize=(12, 5))
ax.semilogy(iptf[:, 0] / (z + 1), iptf[:, 1], label="iptf15af")
ax.axvline(wl, color="k")
SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
plt.show()
