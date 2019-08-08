#!/usr/bin/env python

import tde_util as tu
import astropy.units as u
from astropy.table import Table
from dust_extinction.parameter_averages import F99
from matplotlib import pyplot as plt

"""
iPTF16fnl:
    Av = 0.22
    Rv = 3.1
"""

iptf16fnl = tu.iPTF16fnl_spec(5)
Av = 0.22
Rv = 3.1
Ebv = Av / Rv
print(Ebv)
ext = F99(Rv=Rv)

plt.semilogy(iptf16fnl[:, 0], iptf16fnl[:, 1] / ext.extinguish(iptf16fnl[:, 0] * u.angstrom, Ebv=Ebv))
plt.show()
