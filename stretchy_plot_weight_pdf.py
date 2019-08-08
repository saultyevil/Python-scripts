#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Plots the PDF of the weight due to a single stretched path.

Usage:
    [python] stretchy_plot_weight_pdf.py fname alpha nscats

        - fname: the filename of the data to read in
        - alpha: the value of the stretching parameter
        - nscats: the photons to extract based on number of scatters
"""


import sys
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


def weight_pdf(weight, alpha):
    """
    Evaluate Christian's magnum opus equation
    """

    g1 = alpha ** (1 - alpha / (alpha - 1)) / (1 - alpha)
    g2 = weight ** (-alpha / (alpha - 1) - 1)

    return g1 * g2


if __name__ == "__main__":

    argv = sys.argv
    argc = len(argv)
    if argc != 3:
        print(__doc__)
        sys.exit(1)
    else:
        fname = argv[1]
        alpha = float(argv[2])

    sim = pd.read_csv(fname, delim_whitespace=True)
    # try:
    #     phots = sim[sim["nscats"] == nscats]
    # except KeyError:
    #     print("Invalid table format, generated from stretchy.py pls")
    #     sys.exit(1)

    pweight = sim.values.astype(float)
    print(pweight)
    print("maxw", pweight.max())
    print("minw", pweight[pweight!=0].min())

    minw = pweight[pweight!=0].min()
    maxw = pweight.max()
    nweight = 100
    w = np.linspace(minw, maxw, nweight)
    g = weight_pdf(w, alpha)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    nbins = nweight


    # ax.hist(pweight, cumulative=True, density=True, bins=nbins)

    # ax.hist(pweight, nbins, cumulative=True, density=True)

    norm = np.ones_like(pweight) / float(len(pweight))
    logbins = np.logspace(np.log10(pweight[pweight != 0].min()), np.log10(pweight.max()), nbins)

    ax.hist(pweight, density=True, bins=nbins, cumulative=True)  # weights=norm
    # ax.set_xscale("log")

    # ax.set_xscale("log")

    plt.loglog(w, g)

    plt.show()
