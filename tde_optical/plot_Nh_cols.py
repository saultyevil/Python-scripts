#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as plt
from path import *

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

col = []
grid = mbh_grid + rmin_grid + vinf_grid
name = grid.copy() + ["fiducial"]

lstyle = ["-", "--", "-."] * 3 + [":"]
col = ["C0"] * 3 + ["C1"] * 3 + ["C2"] * 3 + ["C3"]

fig, ax = plt.subplots(1, 1, figsize=(12, 6))

modelspecs = get_the_models(grid, "diag_tde_uv/tde_uv_0.diag") + \
	["/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/diag_tde_cv/tde_cv_0.diag"]

for i, m in enumerate(modelspecs):
	print(m)
	t = []
	with open(m, "r") as f:
		lines = f.readlines()
	for l in lines:
		if l.find("Hydrogen column density:") != -1:
			tmp = l.split()
			angle = tmp[0][1:3]
			if float(angle) < 30 or float(angle) > 75:
				continue
			rho = tmp[4]
			if rho != 0:
				t.append([angle, rho])
	t = np.array(t, dtype=float)
	t = np.unique(t, axis=0)
	ax.semilogy(t[:, 0], t[:, 1], label=name[i], linewidth=2, linestyle=lstyle[i], alpha=0.7, color=col[i])

ax.set_ylabel(r"Hydrogen Column Density N$_{H}$ [cm$^{-2}$]", fontsize=15)
ax.set_xlabel(r"Inclination angle [$^{\circ}$]", fontsize=15)
ax.set_xlim(30, 75)
ax.legend(ncol=2, fontsize=13)
fig.tight_layout()
fig.savefig("Nh_col.png")
fig.savefig("Nh_col.pdf")
plt.show()
