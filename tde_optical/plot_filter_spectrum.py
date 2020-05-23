#!/usr/bin/env python


from sys import argv, exit
import numpy as np
from matplotlib import pyplot as plt
from PyPython.SpectrumUtils import smooth, read_spec, spec_inclinations
from PyPython.SpectrumUtils import common_lines, plot_line_ids
from PyPython.Constants import C
from PyPython.Reverb import create_spectra_from_delay_dump


argc = len(argv)
if argc != 4:
    print("Expected three arguments: root directory norm")
    exit(1)

root = argv[1]
wd = argv[2]
norm = float(argv[3])

full_spectrum = read_spec("{}/{}.log_spec".format(wd, root))
fmin = np.min(full_spectrum["Freq."])
fmax = np.max(full_spectrum["Freq."])
nbins = len(full_spectrum["Freq."])

try:
    filtered_spectrum = np.loadtxt("{}/{}.delay_dump.spec".format(wd, root))
except IOError:
    filtered_spectrum = create_spectra_from_delay_dump(root, wd, fmin, fmax, nbins, norm, True, True)

smooth_amount = 5
inclinations = spec_inclinations(full_spectrum)

print("Plotting spectra for ", inclinations, "degrees")
for i, inclination in enumerate(inclinations):
    fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(C * 1e8 / filtered_spectrum[:-1, 0], smooth(filtered_spectrum[:-1, i + 1], smooth_amount),
            label="Filtered spectrum", linewidth=1.4, alpha=0.75)
    ax.plot(full_spectrum["Lambda"].values, smooth(full_spectrum[inclination].values, smooth_amount),
            label="Spectrum", linewidth=1.4, alpha=0.75)

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_xlabel(r"Wavelength [$\AA$]", fontsize=15)
    ax.set_ylabel(r"Flux at 100 Pc [ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", fontsize=15)
    ax.legend(fontsize=13)  # loc="lower left")
    ax = plot_line_ids(ax, common_lines(), logx=True)

    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    fig.savefig("{}/{}_i{}.delay_dump_spectrum.png".format(wd, root, inclination), dpi=300)
    fig.savefig("{}/{}_i{}.delay_dump_spectrum.pdf".format(wd, root, inclination), dpi=300)

    plt.close()
