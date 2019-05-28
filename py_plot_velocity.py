#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import py_plot_util

root = "tde.pf"

xx, yx, vx = py_plot_util.get_wind_data(root, "vx", "wind")
xy, yy, vy = py_plot_util.get_wind_data(root, "vy", "wind")
xz, yz, vz = py_plot_util.get_wind_data(root, "vz", "wind")

v = np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)

