from matplotlib import pyplot as plt
import numpy as np
from sys import argv, exit
import pandas as pd

argc = len(argv)
if argc != 2:
    print("Require photon number")
    exit(1)
try:
    nphot = int(argv[1])
except ValueError:
    if argv[1] != "find":
        print("Unable to convert {} to int".format(argv[1]))
        exit(1)

photons = np.loadtxt("python.ext.txt", skiprows=1, dtype=str)
photons = photons[:, 1:-2].astype(float)
with open("python.ext.txt", "r") as f:
    headers = f.readline().replace(",", " ").split()
headers = headers[1:-1]
df = pd.DataFrame(photons, columns=headers)
if argv[1] == "find":
    phot_istat = df[df["p->istat"] == 3]
    print(phot_istat["p->np"].values.astype(int))
    exit(0)
phot_single = df[df["p->np"] == nphot]
print(phot_single)
final_istat = phot_single["p->istat"].values.astype(int)
final_istat = final_istat[-1]
if final_istat != 3:
    print("Photon wasn't killed due to hitting star")
    exit(1)
x = phot_single["p->x[0]"].values.astype(float)
z = phot_single["p->x[1]"].values.astype(float)


fig, ax = plt.subplots(1, 2, figsize=(15, 8))

ax[0].plot(x, z, "x-")
ax[0].plot(x[0], z[0], "*", label="start")
ax[0].plot(x[-1], z[-1], "D", label="end")
ax[0].set_xlabel("x")
ax[0].set_ylabel("z")
ax[0].axhline(y=0, color='k')
ax[0].axvline(x=0, color='k')
ax[0].legend()

ax[1].plot(x, z, "x-")
ax[1].plot(x[0], z[0], "*", label="start")
ax[1].plot(x[-1], z[-1], "D", label="end")
ax[1].set_xlabel("x")
ax[1].set_ylabel("z")
ax[1].axhline(y=0, color='k')
ax[1].axvline(x=0, color='k')
ax[1].set_xscale("symlog")
ax[1].set_yscale("symlog")
ax[1].legend()


plt.show()
