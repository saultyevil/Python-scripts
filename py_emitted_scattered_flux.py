#!/usr/bin/env python3


import sys
import py_plot_util
import numpy as np
from matplotlib import pyplot as plt


def plot_emitted_scattered_flux():
    """
    Plot the emitted and scattered flux of a Python simulation
    """  
    
    try:
        SMOOTH = sys.argv[1]
    except:
        SMOOTH = 1
    
    spec_file = py_plot_util.find_spec_files()
    if len(spec_file) > 1:
        print("Can only work with 1 spec file at once")
        return
    spec_file = spec_file[0]  # as find_spec_files returns a list

    rootname, filepath = py_plot_util.get_root_name_and_path(spec_file)
    con_frac = py_plot_util.check_convergence(filepath, rootname)
    with open(spec_file, "r") as f:
        ver = []
        nlines = 2
        for i in range(nlines):
            ver.append(f.readline().replace("\n", "").split())

    headers = ["Created", "Scattered"]
    spec = py_plot_util.read_spec_file(spec_file, " ")
    wavelength = np.array(spec[1:, spec[0, :] == "Lambda"], dtype=float)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    for i in range(len(headers)):
        flux = np.array(spec[1:, spec[0, :] == headers[i]], dtype=float)
        smoothflux = py_plot_util.smooth_flux(flux, SMOOTH)
        ax.plot(wavelength, smoothflux, label=headers[i])
    ax.set_xlim(wavelength.min(), wavelength.max())
    ax.set_xlabel(r"Wavelength ($\AA$)", fontsize=17)
    ax.set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)",
                          fontsize=17)
    ax.set_title("Convergence fraction = {}\nPython Version = {}"
                 "\nGit commit = {}".format(con_frac, ver[0][-1], ver[1][-1][:7]))
    ax.legend()

    plt.savefig("{}_emitted_scattered_flux.png".format(rootname))

    return


if __name__ == "__main__":
    plot_emitted_scattered_flux()
