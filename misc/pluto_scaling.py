#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt

# small ratio grid on Iridis 5
ncores = [1, 2, 4, 8, 16, 32, 40, 80, 120, 160, 200]
times = [299, 151, 78, 41, 23, 15, 15, 33, 41, 43, 47]

plt.figure(figsize=(12, 8))
plt.plot(ncores, times[0] / np.array(times), "--D", label="Pluto")
plt.plot(ncores, ncores, "--", label="Ideal scaling")
plt.xlabel("Number of cores (40 cores per node)")
plt.ylabel("T$_{n}$ / T$_{1}$")
plt.xlim(0, 200)
plt.ylim(0, 25)
plt.legend()
plt.savefig("pluto_small_ratio_iridis5_scaling.png")
plt.show()
