from PyPython import SpectrumUtils, Utils
from matplotlib import pyplot as plt
import tde_spectra

spec_files = SpectrumUtils.find_specs()
spec_files = sorted(spec_files)
nspecs = len(spec_files)
print(spec_files)

PARSEC = 3.086E18
DEFAULT_DIST = 100 * PARSEC
SMOOTH = 3
nrows = 2
ncols = 3
fig, ax = plt.subplots(nrows, ncols, figsize=(15, 8))

j = 0
k = 0

for i in range(nspecs):

    if j > nrows - 1:
        k += 1
        j = 0

    print(j, k)

    ipt = tde_spectra.iptf15af_spec(SMOOTH)

    spec = SpectrumUtils.read_spec(spec_files[i])
    wl = spec["Lambda"].values.astype(float)
    if spec_files[i].find("agn") != -1:
        fl = spec["75"].values.astype(float)
    else:
        fl = spec["60"].values.astype(float)
    fl = SpectrumUtils.smooth_spectrum(fl, SMOOTH)
    fl *= DEFAULT_DIST ** 2 / (350 * 1e6 * PARSEC) ** 2

    root, wd = Utils.split_root_directory(spec_files[i])


    ax[j, k].semilogy(ipt[:, 0] / (0.07897 + 1), ipt[:, 1])
    ax[j, k].semilogy(wl, fl, label=wd)
    ax[j, k].set_xlabel("wavelength")
    ax[j, k].set_ylabel("flux")
    ax[j, k].legend()
    ax[j, k].set_xlim(1000, 3000)
    ymax, ymin = SpectrumUtils.ylims(wl, fl, 1000, 3000, scale=5)
    ax[j, k].set_ylim(ymin, ymax)
    SpectrumUtils.plot_line_ids(ax[j, k], SpectrumUtils.common_lines())

    j += 1

fig.tight_layout()
plt.savefig("comparison.png")
plt.show()
