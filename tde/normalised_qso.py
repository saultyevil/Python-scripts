#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import sys
from platform import system
import numpy as np

if system() == "Darwin":
    sys.path.append("/Users/saultyevil/Scripts")
else:
    sys.path.append("/home/saultyevil/Scripts")

import tde_util as tu

qso = tu.sdss_qso_spec(verbose=True)
twl = 1550
wl = qso[:, 0]
fl = qso[:, 1]
idx = 0
for i in range(len(wl)):
    idx = i
    if wl[i] > twl:
        break
tfl = fl[idx]
print(wl, fl)
print(idx, tfl)
