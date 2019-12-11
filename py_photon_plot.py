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
from typing import Tuple, List, Union


def get_input() -> Tuple[str, int, str]:
    """Get the input choices from the command line."""

    p = ap.ArgumentParser(description=__doc__)
    p.add_argument("mode", type=str, help="The mode for this script to run as")
    p.add_argument("wcycle", type=int, help="The ionisation cycle to extract photons from")
    p.add_argument("extract", type=str, help="The value to extract a photon with")
    args = p.parse_args()

    return args.mode, args.wcycle, args.extract


def read_photon_file(fname: str = "python.ext_0.txt", myheader: List[str] = None) -> pd.DataFrame:
    """Read in the Python photon extra diagnostic file."""

    header = ["PHOTON", "wcycle", "np", "freq", "w", "x", "y", "z", "nx", "ny", "nz", "grid", "istat", "origin",
              "nres", "nscat", "nrscat", "comment"]
    if myheader:
        header = myheader
    df = pd.read_csv(fname, delim_whitespace=True, names=header)
    df = df.drop("PHOTON", axis=1)

    return df


def extract_photons(df: pd.DataFrame, wcycle: int, column: str, extract: Union[str, int, float]) -> pd.DataFrame:
    """Extract photons from a certain cycle and with a certain comment and
    change the data type to numeric."""

    if type(wcycle) is not int:
        wcycle = int(wcycle)
    if type(column) is not str:
        column = str(column)

    backup = df.copy(deep=True)

    df = df[df["wcycle"] == wcycle]
    if len(df.index) == 0:
        print("No photons from cycle {}. Returning original.".format(wcycle))
        return backup
    df = df[df[column] == extract]
    if len(df.index) == 0:
        print("No photons with value {} = {}. Returning original.".format(column, extract))
        return backup
    df = df.drop(["wcycle", "comment"], axis=1)
    df = df.apply(pd.to_numeric)

    del backup  # Just in case we have any memory leaks because of Pandas

    return df


def bin_photon_weights_in_frequency(freq: np.ndarray, w: np.ndarray, nbins: int = 100) -> np.ndarray:
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


def plot_photon_frequency_distribution() -> None:
    """Plot the photons weights binned in frequency space"""

    ncycles = 10
    nrows = 2
    ncols = 5

    # header = ["PHOTON", "wcycle", "np", "freq", "w", "x", "y", "z", "nx", "ny", "nz", "grid", "istat", "origin",
    #           "nres", "comment"]

    for ii in range(2):

        k = 0
        fig, ax = plt.subplots(nrows, ncols, figsize=(18, 6))

        for i in range(ncycles):
            wcycle = i

            photons = read_photon_file()
            photons_b4 = extract_photons(photons, wcycle, "comment", "beforeTransport")
            photons_af = extract_photons(photons, wcycle, "comment", "afterTransport")

            freq_b4 = photons_b4["freq"].values.astype(float)
            w_b4 = photons_b4["w"].values.astype(float)
            hist_b4 = bin_photon_weights_in_frequency(freq_b4, w_b4)
            freq_af = photons_af["freq"].values.astype(float)
            w_af = photons_af["w"].values.astype(float)
            hist_af = bin_photon_weights_in_frequency(freq_af, w_af)

            print("\nCycle {}".format(wcycle))
            print("---------")
            print("B4: Minimum frequency = {:1.3e} Hz".format(freq_b4.min()))
            print("AF: Minimum frequency = {:1.3e} Hz\n".format(freq_af.min()))
            print("B4: Maximum frequency = {:1.3e} Hz".format(freq_b4.max()))
            print("AF: Maximum frequency = {:1.3e} Hz\n".format(freq_af.max()))
            print("B4: Total photon weight = {:1.3e}".format(np.sum(w_b4)))
            print("AF: Total photon weight = {:1.3e}\n".format(np.sum(w_af)))

            if i > ncols - 1:
                j = 1
            else:
                j = 0

            if k > ncols - 1:
                k = 0

            if ii == 0:
                ax[j, k].loglog(hist_b4[:, 0], hist_b4[:, 1], label="beforeTransport")
                ax[j, k].loglog(hist_af[:, 0], hist_af[:, 1], label="afterTransport")
            elif ii == 1:
                ax[j, k].loglog(hist_b4[:, 0], np.cumsum(hist_b4[:, 1]), label="beforeTransport")
                ax[j, k].loglog(hist_af[:, 0], np.cumsum(hist_af[:, 1]), label="afterTransport")

            ax[j, k].set_xlabel(r"$\nu$")
            if ii == 0:
                ax[j, k].set_ylabel(r"$W$")
            elif ii == 1:
                ax[j, k].set_ylabel(r"$W_{tot}(\nu^{*} < \nu)$")
            ax[j, k].set_title("Cycle {}".format(wcycle))

            j += 1
            k += 1

        ax[-1, -1].legend(loc="lower center")
        fig.tight_layout()
        if ii == 0:
            plt.savefig("photon_freq_hist.png")
        elif ii == 1:
            plt.savefig("photon_cumulative_freq_hist.png")
        plt.show()

    return


def plot_weight_hist():
    """Simply use matplotlib to create a histogram of the photon weights for
    each cycle."""

    ncycles = 10
    nrows = 2
    ncols = 5

    k = 0
    fig, ax = plt.subplots(nrows, ncols, figsize=(18, 6))

    for i in range(ncycles):
        wcycle = i

        photons = read_photon_file()

        photons_b4 = extract_photons(photons, wcycle, "comment", "beforeTransport")
        photons_af = extract_photons(photons, wcycle, "comment", "afterTransport")
        w_b4 = photons_b4["w"].values.astype(float)
        w_af = photons_af["w"].values.astype(float)

        if i > ncols - 1:
            j = 1
        else:
            j = 0

        if k > ncols - 1:
            k = 0

        nbins = 50
        ax[j, k].hist(w_b4, label="beforeTransport", alpha=0.5)
        ax[j, k].hist(w_af, label="afterTransport", alpha=0.5)
        ax[j, k].set_xscale("log")
        ax[j, k].set_title("W_bf = {:3.2e} : W_af = {:3.2e}".format(np.sum(w_b4), np.sum(w_af)))

        j += 1
        k += 1

    ax[-1, -1].legend(loc="lower center")
    fig.tight_layout()
    plt.savefig("photon_weight_hist.png")
    plt.show()

    return


def extract_scattered_photons(wcycle: int, nscats: str):
    """Extract a photon which has scattered nscats time from the ionisation
    cycle wcycle."""

    try:
        nscats = int(nscats)
    except ValueError:
        print("Unable to convert nscats = {} into int".format(nscats))
        return

    photons = read_photon_file()
    nphotons = len(photons.index)
    photons = extract_photons(photons, wcycle, "nscat", nscats)

    if len(photons.index) == nphotons:
        print("All photons returned - something has probably gone wrong!")
        return

    print(photons)

    return


def main():
    """Main function"""

    modes = ["frequency_distribution", "extract_scatters", "weight_hist"]
    mode, wcycle, extract = get_input()
    if mode == modes[0]:
        plot_photon_frequency_distribution()
    elif mode == modes[1]:
        extract_scattered_photons(wcycle, extract)
    elif mode == modes[2]:
        plot_weight_hist()
    else:
        print("Unknown mode: {}".format(mode))
        print("Known modes: ")
        for mode in modes:
            print(" - {}".format(mode))

    return


if __name__ == "__main__":
    main()
