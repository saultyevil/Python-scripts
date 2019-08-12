#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import py_plot_util as ppu
import py_plot as pp
from matplotlib import pyplot as plt

root = "tde_spherical_cno"

xx, yx, vx = ppu.get_wind_data(root, "v_x", "wind", coord="polar")
xy, yy, vy = ppu.get_wind_data(root, "v_y", "wind", coord="polar")
xz, yz, vz = ppu.get_wind_data(root, "v_z", "wind", coord="polar")

v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

pp.plot_polar_wind(xx, yx, v, 0, "v", "wind", (1, 1))
plt.show()
