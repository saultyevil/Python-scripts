#!/usr/bin/env python 3
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
from path import get_the_models, mbh_labels
from typing import List, Tuple
from PyPython.SpectrumUtils import smooth

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15

mbh_grid = [
    "Mbh/continuum/1.0000e+06",
    "Mbh/continuum/1.0000e+07",
    "Mbh/continuum/1.0000e+08"
]

alpha = 0.75
sm = 5


def plot(mbh_grid: List[str], labels: List[str], inclination: str) \
        -> Tuple[plt.Figure, plt.Axes]:
    """Plot the underlying continuum for the Mbh grid"""

    spectra = get_the_models(mbh_grid, "tde_uv.spec")

    fig, ax = plt.subplots(1, 1, figsize=(12, 5))

    for i, s in enumerate(spectra):
        ax.loglog(s["Lambda"], smooth(s[inclination], sm) * s["Lambda"], linewidth=2, alpha=alpha, label=labels[i])

    ax.legend()
    ax.set_xlabel(r"Wavelength $\lambda$ [$\AA$]", fontsize=15)
    ax.set_ylabel(r" $\lambda$ $F_{\lambda}$  at 100 pc [erg s$^{-1}$ cm$^{-2}$]", fontsize=15)

    fig.savefig("misc/Mbh_continuum_i{}.pdf".format(inclination))
    plt.close()

    return fig, ax


if __name__ == "__main__":

    incl = ["10", "60", "75"]

    for i in incl:
        plot(mbh_grid.copy(), mbh_labels, i)
