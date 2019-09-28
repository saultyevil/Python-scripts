#!/usr/bin/env python3

import py_plot_util
import numpy as np
from os import mkdir, path
from matplotlib import pyplot as plt

plot_show = False
FILE_TYPE = "png"
dirname = "sout"  # change to . to output into the current directory


def read_and_reshape_data(filename):
    """
    Read in the Snake grid output and reshape it into a 3d array of
    (ncycles, ncells, ncols).

    Parameters
    ----------
    None

    Returns
    -------
    sgrid: (ncycles, ncells, ncols) array of str
        The reshaped grid output from Snake. The columns are as follows:
            col 0: cell number
            col 1: spatial coordinate
            col 2: cell density
            col 3: cell rosseland mean opacity
            col 4: cell optical depth
            col 5: cumulative optical depth
            col 6: cell temperature
    """

    sgrid = np.array(py_plot_util.read_spec(filename), dtype=float)

    # Figure out the number of cells, cycles and cols from the dimensions
    ncells = int(sgrid[:, 0].max()) + 1
    ncycles = int(sgrid.shape[0] / ncells)
    ncols = sgrid.shape[1]

    return np.reshape(sgrid, (ncycles, ncells, ncols))


def plot_each_var_and_cycle(sgrid):
    """
    The function for plotting the cell conditions as a function of disk height.

    Parameters
    ----------
    sgrid: (ncycles, ncells, ncols) array of str
        The reshaped grid output from Snake

    Returns
    -------
    None
    """

    print("\n--------------------------------------------------------------------------------\n")
    print(" Plotting a 4x4 grid of the cell conditions")
    print("\n--------------------------------------------------------------------------------")

    nrows = 2
    ncols = 2
    ncycles = sgrid.shape[0]

    # Output the plots to a directory if desired
    if path.exists(dirname):
        print("\n Outputting plots into the directory '{}'\n".format(dirname))
    else:
        print("\n Creating directory '{}' for plots\n".format(dirname))
        mkdir(dirname)

    # Loop over the cycles and hopefully output successfully :-)
    for i in range(0, ncycles):
        print(" Plotting cycle {}".format(i))
        inds = [[2, 3], [5, 6]]
        ylabs = [[r"Density, $\rho$", r"Rosseland Opacity, $\kappa_{R}$"],
                 [r"Optical Depth To Escape, $\tau$", r"Temperature, $T$"]]
        fig, ax = plt.subplots(nrows, ncols, figsize=(15, 15))
        for j in range(ncols):
            for k in range(nrows):
                # Horrid hack so I don't have to see a matplotlib error which I don't know how to suppress
                if sgrid[i, :, inds[j][k]].min() == sgrid[i, :, inds[j][k]].max():
                    ax[j, k].set_ylim(0.1 * sgrid[i, :, inds[j][k]].min(), 10.0 * sgrid[i, :, inds[j][k]].min())
                # Use a log scale for some plots, I guess
                if ylabs[j][k] == r"Optical Depth To Escape, $\tau$" or ylabs[j][k] == r"Temperature, $T$":
                    ax[j, k].plot(sgrid[i, :, 1], sgrid[i, :, inds[j][k]])
                else:
                    ax[j, k].semilogy(sgrid[i, :, 1], sgrid[i, :, inds[j][k]])
                ax[j, k].set_xlim(sgrid[i, :, 1].min(), sgrid[i, :, 1].max())
                ax[j, k].set_xlabel("Height, $z$", fontsize=13)
                ax[j, k].set_ylabel(ylabs[j][k], fontsize=13)
        fig.suptitle("Snake: cycle {}".format(i), y=0.95, fontsize=25)
        plt.savefig("{}/cycle_{}.{}".format(dirname, i, FILE_TYPE))
        if plot_show:
            plt.show()
        else:
            plt.close()

    # # Plot above but only cycles 0, 1, 2 and the final
    # fig, ax = plt.subplots(nrows, ncols, figsize=(15, 15))
    # for i in [0, 1, 2, -1]:
    #     inds = [[2, 3], [5, 6]]
    #     ylabs = [[r"Density, $\rho$", r"Rosseland Opacity, $\kappa_{R}$"],
    #              [r"Optical Depth To Escape, $\tau$", r"Temperature, $T$"]]
    #     markers = ["-", "--", ":", "-."]
    #     for j in range(ncols):
    #         for k in range(nrows):
    #             # Horrid hack so I don't have to see a matplotlib error which I don't know how to suppress
    #             if sgrid[i, :, inds[j][k]].min() == sgrid[i, :, inds[j][k]].max():
    #                 ax[j, k].set_ylim(0.1 * sgrid[i, :, inds[j][k]].min(), 10.0 * sgrid[i, :, inds[j][k]].min())
    #             if j == 0 and k == 0:
    #                 ax[j, k].semilogy(sgrid[i, :, 1], sgrid[i, :, inds[j][k]], markers[i], label="Cycle {}".format(i))
    #             else:
    #                 ax[j, k].semilogy(sgrid[i, :, 1], sgrid[i, :, inds[j][k]], markers[i])
    #             ax[j, k].set_xlim(sgrid[i, :, 1].min(), sgrid[i, :, 1].max())
    #             ax[j, k].set_xlabel("Disk height, $z$", fontsize=13)
    #             ax[j, k].set_ylabel(ylabs[j][k], fontsize=13)
    # fig.legend()
    # fig.suptitle("Snake: cycles {}".format([0, 1, 2, -1]), y=0.95, fontsize=25)
    # plt.savefig("{}/cycle_comparison.{}".format(dirname, filetype))
    # if plot_show:
    #     plt.show()
    # else:
    #     plt.close()

    return


def plot_comparsion():
    """
    Quick function to plot a comparison of two simulation runs. Hard coding the file names for now...
    """

    sgridlow = read_and_reshape_data("sgrid_1e-8.out")
    sgridhigh = read_and_reshape_data("sgrid_1e-5.out")

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.semilogy(sgridlow[-1, :, 1], sgridlow[-1, :, 6], label=r"$\rho = 10^{-8}$ g/cm$^{3}$")
    ax.semilogy(sgridhigh[-1, :, 1], sgridhigh[-1, :, 6], label=r"$\rho = 10^{-5}$ g/cm$^{3}$")
    ax.set_xlim(sgridlow[-1, :, 1].min(), sgridlow[-1, :, 1].max())
    ax.set_xlabel("Height, $z$", fontsize=13)
    ax.set_ylabel("Temperature, $T$", fontsize=13)
    ax.legend()
    fig.suptitle(r"$\rho = 10^{-8}$ g/cm$^{3}$ Vs. $\rho = 10^{-5}$ g/cm$^{3}$")
    plt.savefig("temperature_comparison.{}".format(FILE_TYPE))
    if plot_show:
        plt.show()
    else:
        plt.close()

    return

def main():
    """
    The main steering function of the script.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """

    print("\n--------------------------------------------------------------------------------\n")
    print(" Snake plotting script")

    sgrid = read_and_reshape_data("sgrid.out")
    plot_each_var_and_cycle(sgrid)
    # plot_comparsion()

    print("\n--------------------------------------------------------------------------------\n")

    return sgrid


if __name__ == "__main__":
    sgrid = main()
