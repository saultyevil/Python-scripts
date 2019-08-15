#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from datetime import datetime
import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams
import pandas as pd


def plot_weight_tau():
    """
    Plot the possible photon weight change as a function of tau_scat.
    """

    def new_weight(tau, alpha):
        """
        Return the value by which the original weight is modified by.
        """

        return np.exp(tau * (alpha - 1)) / alpha

    def path_alpha(tau_path):
        """
        Return the value of alpha for a given tau_path.
        """

        return 1 / (1 + tau_path)

    min_tau = 0
    max_tau = 1e3
    ntau = int(1e4)
    tau = np.linspace(min_tau, max_tau, ntau)

    # min_tau_path = 1e-7
    # max_tau_path = 1e7
    # ntau_path = 10  # int((max_tau - min_tau) + 2)
    # tau_path = np.linspace(min_tau_path, max_tau_path, ntau_path)
    alpha = [1.0, 0.9, 0.8, 0.5, 0.3, 0.1, 0.01, 0.0001, 1e-06, 1e-08, 1e-10]
    ntau_path = len(alpha)

    weight_factors = np.zeros((ntau_path, ntau))
    for i in range(ntau_path):
        # alpha = path_alpha(tau_path[i])

        weight_factors[i, :] = new_weight(tau, alpha[i])

    plt.figure(figsize=(12, 8))
    for i in range(ntau_path):
        plt.loglog(tau, weight_factors[i, :], label=r"$\alpha$ = " + "{:.2e}".format(alpha[i]))

    plt.legend()
    plt.xlabel(r"$\tau_{scat}$")
    plt.ylabel(r"$W(\tau_{scat})$")
    plt.savefig("path_stretching_weight.png")
    plt.close()

    return


def plot_weight_file(fname):
    """
    Plot the photon weights from a file generated by stretchy.py

    path_stretch_on      bool
                         open the pas file instead of true
    """

    with open(fname, "r") as f:
        data = []
        header = f.readline().split()
        for line in f:
            l = line.split()
            data.append(l)
    data = np.array(data, dtype=float)

    header = header[2]
    data = data[:, 2]
    # print(header)
    # print(data)
    # print(data[data != 0].min(), data.max())

    nbins = 20

    try:
        nplots = data.shape[1]
    except IndexError:
        nplots = 1

    size = (5, 12)
    if nplots == 1:
        size = (10, 5)

    fig, ax = plt.subplots(nplots, 1, figsize=size)
    for i in range(nplots):
        try:
            weights = np.ones_like(data[:, i])/float(len(data[:, i]))
            ax[i].hist(data[:, i], nbins,  weights=weights, log=True)
            ax[i].set_xlabel("{}".format(header[i]))
            ax[i].set_ylabel("Count")
        except IndexError:
            nbins = 200
            weights = np.ones_like(data)/float(len(data))
            logbins = np.logspace(np.log10(data[data!=0].min()), np.log10(data.max()), nbins)
            # print(logbins)
            ax.hist(data, bins=logbins,  weights=weights, log=True)
            ax.set_xscale("log")
            # ax.set_xlim(np.log10(data[data!=0].min()), np.log10(data.max()))
            ax.set_xlabel("{}".format(header[i]))
            ax.set_ylabel("Count")

    fig.tight_layout()
    plt.savefig("photon_hists.png")
    plt.close()

    return


def plot_tau_alpha_hists(tau_scat, tau_path, alpha, extra_name):
    """
    Plot histograms for the sampled tau_scats and the values of alpha and
    tau_path.

    pls ignore the ABSOLUTE state of this function
    """

    if type(tau_scat) != np.ndarray:
        tau_scat = np.array(tau_scat, dtype=float)
    if type(tau_path) != np.ndarray:
        tau_path = np.array(tau_path, dtype=float)
    if type(alpha) != np.ndarray:
        alpha = np.array(alpha, dtype=float)

    # print(tau_scat)
    # print(tau_path)
    # print(alpha)

    # print(alpha.min(), alpha.max())

    nplots = 3
    tmp = np.zeros((tau_scat.shape[0], nplots))
    tmp[:, 0] = tau_scat
    tmp[:, 1] = tau_path
    tmp[:, 2] = alpha
    name = ["tau_scat", "tau_path", "alpha"]
    nbins = 200

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    i = 0
    # print(name[i])
    data = tmp[:, i]
    weights = np.ones_like(data) / float(len(data))
    logbins = np.logspace(np.log10(data[data != 0].min()), np.log10(data.max()), nbins)
    ax.hist(data, bins=logbins, weights=weights, log=True)
    ax.set_xscale("log")
    ax.set_xlabel("{}".format(name[i]))
    fig.tight_layout()
    plt.savefig("{}_{}.png".format(name[i], extra_name))
    plt.close()
    fig.clf()
    i += 1

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    # print(name[i])
    data = tmp[:, i]
    # print(data, type(data), type(data[0]))
    weights = np.ones_like(data) / float(len(data))
    logbins = np.logspace(np.log10(data[data != 0].min()), np.log10(data.max()), nbins)
    ax.hist(data, bins=logbins, weights=weights, log=True)
    ax.set_xlabel("{}".format(name[i]))
    ax.set_xscale("log")
    fig.tight_layout()
    plt.savefig("{}_{}.png".format(name[i], extra_name))
    plt.close()
    fig.clf()
    i += 1

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    # print(name[i])
    data = tmp[:, i]
    weights = np.ones_like(data) / float(len(data))
    logbins = np.logspace(np.log10(data[data != 0].min()), np.log10(data.max()), nbins)
    ax.hist(data, log=True, weights=weights, bins=logbins)
    ax.set_xlabel("{}".format(name[i]))
    ax.set_xscale("log")
    fig.tight_layout()
    plt.savefig("{}_{}.png".format(name[i], extra_name))
    plt.close()
    fig.clf()
    i += 1

    return


def plot_hist_panels():
    """
    Plot a bunch of things
    """

    nrows = 3
    ncols = 4
    nbins = 250

    tau_scat_filename = "tau_scats_alpha_"
    tau_paths_filename = "tau_paths_alpha_"
    sim_alphas = ["1.0", "0.9", "0.8", "0.5", "0.3", "0.1", "0.01", "0.0001", "1e-06", "1e-08", "1e-10", "taupath"]
    sim_alphas.reverse()

    nplots = len(sim_alphas)

    fig1, axpath = plt.subplots(nrows, ncols, figsize=(20, 8), squeeze=False)

    print("Plotting tau paths\n\n")

    iplot = 0
    for i in range(nrows):
        for j in range(ncols):
            if iplot > nplots - 1:
                break
            print("Plotting {} tau_scat".format(iplot), sim_alphas[iplot])
            fnamepath = tau_paths_filename + sim_alphas[iplot] + ".txt"
            print("Reading data in: ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
            datapath = np.loadtxt(fnamepath)
            weightspath = np.ones_like(datapath / float(len(datapath)))
            logbinspath = np.logspace(np.log10(datapath[datapath != 0].min()), np.log10(datapath.max()), nbins)
            axpath[i, j].hist(datapath, log=True, weights=weightspath, bins=logbinspath)
            axpath[i, j].set_xscale("log")
            axpath[i, j].set_xlabel(r"$\tau_{path}$")
            axpath[i, j].set_title(r"$\alpha$ = " + sim_alphas[iplot])

            iplot += 1

    fig1.tight_layout()
    plt.savefig("tau_path_hists.png")
    plt.close()
    fig1.clf()

    fig2, axscat = plt.subplots(nrows, ncols, figsize=(20, 8), squeeze=False)

    print("\n\nPlotting tau scats\n\n")

    iplot = 0
    for i in range(nrows):
        for j in range(ncols):
            if iplot > nplots - 1:
                break
            print("Plotting {} tau_scat".format(iplot), sim_alphas[iplot])
            fnamescat = tau_scat_filename + sim_alphas[iplot] + ".txt"
            print("Reading data in: ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
            datascat = np.loadtxt(fnamescat)
            weightsscat = np.ones_like(datascat / float(len(datascat)))
            logbinsscat = np.logspace(np.log10(datascat[datascat != 0].min()), np.log10(datascat.max()), nbins)
            axscat[i, j].hist(datascat, log=True, weights=weightsscat, bins=logbinsscat)
            axscat[i, j].set_xscale("log")
            axscat[i, j].set_xlabel(r"$\tau_{scat}$")
            axscat[i, j].set_title(r"$\alpha$ = " + sim_alphas[iplot])

            iplot += 1

    fig2.tight_layout()
    plt.savefig("tau_scat_hists.png")
    plt.close()
    fig2.clf()
    #
    print("\n\nPlotting weights\n\n")

    fig3, axweight = plt.subplots(nrows, ncols, figsize=(20, 8), squeeze=False)

    iplot = 0
    for i in range(nrows):
        for j in range(ncols):
            if iplot > nplots - 1:
                break
            print("Plotting {} weights".format(iplot), sim_alphas[iplot])
            fnameweight = "photon_statistics_pasalpha_{}_.txt".format(sim_alphas[iplot])
            if sim_alphas[iplot] == 1.0:
                axweight[i, j].set_title("all weights = 1")
                iplot += 1
                continue
            print("Reading data in: ", datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
            dataweight = np.loadtxt(fnameweight, skiprows=1)
            dataweight = dataweight[:, 2]
            try:
                weightweight = np.ones_like(dataweight / float(len(dataweight)))
                logbinweight = np.logspace(np.log10(dataweight[dataweight != 0].min()), np.log10(dataweight.max()), nbins)
            except ValueError:
                print("ValueError: all weights are probably zero so skipping")
                axweight[i, j].set_title("sum(w) == 0")
                iplot += 1
                continue
            axweight[i, j].hist(dataweight, log=True, weights=weightweight, bins=logbinweight)
            axweight[i, j].set_xscale("log")
            axweight[i, j].set_xlabel(r"$W$")
            axweight[i, j].set_title(r"$\alpha$ = " + sim_alphas[iplot])

            iplot += 1

    fig3.tight_layout()
    plt.savefig("weights_hist.png")
    plt.close()
    fig3.clf()

    print("\n\nDone\n\n")

    return


def plot_nstretch_weight(fname):
    """
    Plot the number of stretched paths against nscat
    """

    data = np.loadtxt(fname, skiprows=1)
    nstretch = data[:, 1]
    weight = data[:, 2]

    print("npoints {:e}".format(len(weight)))

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    weight = np.log10(weight)
    print(weight)

    ax.scatter(nstretch, weight)
    ax.set_xlabel("Number of Stretched Paths")
    ax.set_ylabel("Weight")
    ax.set_yscale("log")
    ax.set_ylim(0, weight.max())

    plt.savefig("nstretch_weight_{}.png".format(fname))
    plt.show()

    return


def plot_intensity():
    """
    Plot many intensity things
    """

    fig, ax = plt.subplots(1,1, figsize=(12, 8))
    sim_alphas = ["1.0", "0.9", "0.8", "0.5", "0.3", "0.1", "0.01", "taupath"]
    iname = "intensity_alpha"
    tname = "theta_alpha"

    nplots = len(sim_alphas)
    for i in range(nplots):
        iiname = iname + "_{}.txt".format(sim_alphas[i])
        ttname = tname + "_{}.txt".format(sim_alphas[i])

        intensity = np.loadtxt(iiname)
        theta = np.loadtxt(ttname)
        ax.semilogy(np.cos(theta), intensity, label=r"$\alpha$ = " + sim_alphas[i])

    CHAND_SLAB_SOL = np.array(
            [[00.00, 1.26938], [18.19, 1.22945], [25.84, 1.18947], [31.79, 1.14943],
             [36.87, 1.10931], [41.41, 1.06911], [45.57, 1.02882], [49.46, 0.98842],
             [53.13, 0.94789], [56.63, 0.90722], [60.00, 0.86637], [63.26, 0.82530],
             [66.42, 0.78938], [69.51, 0.74234], [72.54, 0.70029], [75.52, 0.65770],
             [78.46, 0.61439], [81.37, 0.57001], [84.26, 0.52397], [87.13, 0.47490],
             [90.00, 0.41441]])

    ax.semilogy(np.cos(np.deg2rad(CHAND_SLAB_SOL[:, 0])), CHAND_SLAB_SOL[:, 1], "o-", label="Analytic Solution")

    ax.set_xlabel(r"$\mu = \cos(\theta)$")
    ax.set_ylabel(r"Intensity")
    ax.set_xlim(0, 1)
    ax.legend()

    fig.tight_layout()
    plt.savefig("intensity.png")

    plt.show()

    return


def plot_history_histograms(fname):
    """
    Plot the weight and tau scat history histograms.

    Parameters
    ----------
    fname: str
        The file name of the data to be plotted.

    """

    print("reading file: ", fname)
    d = pd.read_table(fname, delim_whitespace=True, header=None)
    d.columns = ["nphot", "nscats", "w"]

    nbins = 200
    nrows = 5
    ncols = 5
    fig, ax = plt.subplots(nrows, ncols, figsize=(16, 9))
    iscat = 1
    for i in range(nrows):
        for j in range(ncols):
            try:
                w = d[d["nscats"] == iscat]["w"]
            except KeyError:
                continue
            norm_weights = np.ones_like(w / float(len(w)))
            lbins = np.logspace(np.log10(w[w != 0].min()), np.log10(w.max()), nbins)
            ax[i, j].hist(w, log=True, weights=norm_weights, bins=lbins)
            ax[i, j].set_xscale("log")
            ax[i, j].set_title("nscats = {}".format(iscat))
            iscat += 1

    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.show()

    return


if __name__ == "__main__":

    if len(sys.argv) == 3:
        plot = sys.argv[1]
        arg = sys.argv[2]
    else:
        print("usage  : path_stretching_plots_weights.py plot alpha")
        print("choices: weights tau_weight panels tau_alpha_hists nstretch_weight intensity")
        sys.exit(1)

    if plot == "weights":
        plot_weight_file(arg)
    elif plot == "tau_weight":
        plot_weight_tau()
    elif plot == "panels":
        plot_hist_panels()
    elif plot == "tau_alpha_hists":
        alphas = "alphas_alpha_{}.txt".format(arg)
        taupaths = "tau_paths_alpha_{}.txt".format(arg)
        tauscat = "tau_scats_alpha_{}.txt".format(arg)
        now = datetime.now()
        datestring = now.strftime("%d/%m/%Y %H:%M:%S")
        print("reading in data: ", datestring)
        tau_scat = np.loadtxt(tauscat)
        tau_path = np.loadtxt(taupaths)
        arg = np.loadtxt(alphas)
        extra_name = "alpha_{}".format(arg)
        now = datetime.now()
        datestring = now.strftime("%d/%m/%Y %H:%M:%S")
        print("jesus christ it's finally read in:", datestring)
        plot_tau_alpha_hists(tau_scat, tau_path, arg, extra_name)
    elif plot == "nstretch_weight":
        plot_nstretch_weight(arg)
    elif plot == "intensity":
        plot_intensity()
    elif plot == "history_hist":
        plot_history_histograms(arg)
    else:
        print("you fucked up. choices: weights tau_weight panels tau_alpha_hists nstretch_weight intensity")
        sys.exit(1)