#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from PyPython import SpectrumUtils

PARSEC = 3.086E18
SMOOTH = 15
VERBOSITY = False


def plot_optical(mod_spec_path: str, opt_spec_path: str, obj_name: str, days: str):
    """
    Parameters
    ----------
    mod_spec_path: str
        The file path of the model spectrum.
    opt_spec_path:str
        The file path of the optical spectrum (real data).
    obj_name: str
        The name of the object the optical spectrum belongs to.
    days: str
        The time since outburst of the spectrum
    """

    nrows = 2
    ncols = 1
    wmin = 4000
    wmax = 7000

    fig, ax = plt.subplots(nrows, ncols, figsize=(12, 8), sharex="col")

    try:
        mspec = SpectrumUtils.read_spec(mod_spec_path)
    except IOError:
        print("Can't open model spectrum {}".format(mod_spec_path))
        return
    try:
        ospec = np.loadtxt(opt_spec_path)
    except IOError:
        print("Can't open optical spectrum {}".format(opt_spec_path))
        return

    if obj_name.lower() == "at2018zr":
        ospec[:, 1] *= 1e-15

    ax[0].semilogy(ospec[:, 0], SpectrumUtils.smooth_spectrum(ospec[:, 1], 3),
                   label=obj_name + r"$\Delta t$ = " + days)
    ax[0].set_xlim(wmin, wmax)
    # ax[0].set_ylim(1e-17, 5e-16)
    # ax[0].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)

    wl = mspec["Lambda"].values.astype(float)
    for inc in ["10", "28", "45", "62", "80"]:
        fl = mspec[inc].values.astype(float)
        if obj_name.lower() == "at2018zr":
            fl *= (100 * PARSEC) ** 2 / (337 * 1e6 * PARSEC) ** 2
        fl = SpectrumUtils.smooth_spectrum(fl, SMOOTH)
        ax[1].semilogy(wl, fl, label=r"$i = $" + inc + r"$^{\circ}$")

    ax[1].legend(loc="upper right")
    ax[1].set_xlim(wmin, wmax)
    ax[1].set_ylim(1e-17, 5e-16)
    ax[1].set_xlabel(r"Wavelength ($\AA$)", fontsize=15)
    # ax[1].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", fontsize=15)

    fig.text(0.025, 0.5, r"Flux $F_{\lambda}$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]", ha="center", va="center",
             rotation="vertical", fontsize=15)

    SpectrumUtils.plot_line_ids(ax[0], SpectrumUtils.common_lines(), fontsize=13, offset=40)
    SpectrumUtils.plot_line_ids(ax[1], SpectrumUtils.common_lines(), fontsize=13, offset=40)

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.subplots_adjust(hspace=0)

    plt.savefig("{}_optical_comparison.png".format(obj_name))
    plt.show()

    return


if __name__ == "__main__":
    plot_optical("/home/saultyevil/PySims/tde/paper_models/He_II/cv_star_model_original/tde.spec",
                 "/home/saultyevil/PySims/tde/observed_spec/at2018zr_50d_optical.dat",
                 "AT2018zr",
                 "50")
