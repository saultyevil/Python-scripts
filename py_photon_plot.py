#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to look at various things from the photon sample
in a Python simulation. To do this, we need to read in the diagnostic save_photon
file. Python generally is required to be run in the diagnostic
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import argparse as ap
from typing import Tuple


def get_input() -> int:
    """Get the input choices from the command line."""

    p = ap.ArgumentParser(description=__doc__)
    p.add_argument("wcycle", type=int, help="The ionisation cycle to extract photons from")
    # p.add_argument("comment", type=str, help="The photons with comment to extract")
    args = p.parse_args()

    return args.wcycle  # , args.comment


def read_photon_file(fname:str = "python.ext.txt"):
    """Read in the Python photon extra diagnostic file."""

    header = ["PHOTON", "wcycle", "np", "freq", "w", "x", "y", "z", "nx", "ny", "nz", "grid", "istat", "origin",
              "nres", "comment"]
    df = pd.read_csv(fname, delim_whitespace=True, names=header)
    df = df.drop("PHOTON", axis=1)

    return df


def extract_photons(df: pd.DataFrame, wcycle: int, comment: str):
    """Extract photons from a certain cycle and with a certain comment and
    change the data type to numeric."""

    if type(wcycle) is not int:
        wcycle = int(wcycle)

    backup = df.copy(deep=True)

    df = df[df["wcycle"] == wcycle]
    if len(df.index) == 0:
        print("No photons from cycle {}. Returning original.".format(wcycle))
        return backup
    df = df[df["comment"] == comment]
    if len(df.index) == 0:
        print("No photons with comment {}. Returning original.".format(comment))
        return backup
    df = df.drop(["wcycle", "comment"], axis=1)
    df = df.apply(pd.to_numeric)

    return df


def bin_photons(freq: np.ndarray, w: np.ndarray, nbins: int = 10000):
    """Bin the photon weights into frequency bins."""

    freq_min = freq.min()
    freq_max = freq.max()
    dfreq = (freq_max - freq_min) / nbins
    bins = np.linspace(freq_min, freq_max, nbins)
    hist = np.zeros((nbins, 2))
    hist[:, 0] = bins

    for i in range(len(freq)):
        k = int((freq[i] - freq_min) / dfreq)
        if k < 0:
            k = 0
        if k > nbins - 1:
            k = nbins - 1
        hist[k, 1] += w[i]

    return hist


def main():
    """Main function"""

    # wcycle, comment = get_input()
    # wcycle = get_input()

    ncycles = 10
    nrows = 2
    ncols = 5
    fig, ax = plt.subplots(nrows, ncols, figsize=(18, 6))

    j = 0
    k = 0

    for i in range(ncycles):
        wcycle = i

        photons = read_photon_file()
        photons_b4 = extract_photons(photons, wcycle, "beforeTransport")
        photons_af = extract_photons(photons, wcycle, "afterTransport")

        freq_b4 = photons_b4["freq"].values.astype(float)
        w_b4 = photons_b4["w"].values.astype(float)
        hist_b4 = bin_photons(freq_b4, w_b4)
        freq_af = photons_af["freq"].values.astype(float)
        w_af = photons_af["w"].values.astype(float)
        hist_af = bin_photons(freq_af, w_af)

        print("Cycle {}".format(wcycle))
        print("B4: Minimum frequency = {:1.3e} Hz".format(freq_b4.min()))
        print("AF: Minimum frequency = {:1.3e} Hz".format(freq_af.min()))
        print("B4: Maximum frequency = {:1.3e} Hz".format(freq_b4.max()))
        print("AF: Maximum frequency = {:1.3e} Hz".format(freq_af.max()))
        print("B4: Total photon weight = {:1.3e}".format(np.sum(w_b4)))
        print("AF: Total photon weight = {:1.3e}\n".format(np.sum(w_af)))

        if i > ncols - 1:
            j = 1
        else:
            j = 0

        if k > ncols - 1:
            k = 0

        ax[j, k].loglog(hist_b4[:, 0], np.cumsum(hist_b4[:, 1]), label="beforeTransport")
        ax[j, k].loglog(hist_af[:, 0], np.cumsum(hist_af[:, 1]), label="afterTransport")
        ax[j, k].set_xlabel(r"$\nu$")
        ax[j, k].set_ylabel(r"$W_{tot}(\nu^{*} < \nu)$")
        ax[j, k].set_title("Cycle {}".format(wcycle))

        j += 1
        k += 1

    ax[-1, -1].legend(loc="lower center")
    fig.tight_layout()
    plt.savefig("photon_hist.png")
    plt.show()

    return


if __name__ == "__main__":
    main()
