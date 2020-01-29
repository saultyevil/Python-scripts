#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The purpose of this script is to create a Spectral Energy Distribution (SED) for
a cell in Python. To do this, the simulation must have been run using the
matrix_pow ionisation solver.
"""


import argparse as ap
from PyPython import WindUtils
from PyPython import Utils
from PyPython.SpectrumUtils import smooth_spectrum
from subprocess import Popen, PIPE
import numpy as np
from matplotlib import pyplot as plt
from typing import List
from astropy import constants as consts


def get_spec_model(root: str, nx: int, nz: int, i: int, j: int, nbands: int = 4) -> List[str]:
    """
    Get the spectral model for a specific cell from py_wind

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    nx: int
        The number of grid cells in the x direction.
    nz: int
        The number of grid cells in the z direction.
    i: int
        The i-th index for the grid cell in question.
    j: int
        The j-th index for th grid cell in question.
    nbands: int [optional]
        The number of bands in the cell spectral model.

    Returns
    -------
    spectral_model_bands: List[str]
        A list containing the spectral model bands output from py_wind.
    """

    everything_output = py_wind(root, nx, nz, i, j).split("\n")
    for i in range(len(everything_output)):
        line = everything_output[i]
        if line == "Spectral model details:":
            break
    i += 1
    spectral_models_bands = everything_output[i:i + nbands]

    return spectral_models_bands


def py_wind(root: str, nx: int, nz: int, i: int, j: int):
    """
    Run py_wind to get the "everything" output to be able to parse out the
    spectral model bands.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    nx: int
        The number of grid cells in the x direction.
    nz: int
        The number of grid cells in the z direction.
    i: int
        The i-th index for the grid cell in question.
    j: int
        The j-th index for th grid cell in question.
    commands: List[str] [optional]
        Commands for py_wind to run

    Returns
    -------
    stdout: str
        The screen output from py_wind.
    """

    elem = WindUtils.get_wind_elem_number(nx, nz, i, j)
    cmds = np.array(["1", "e", str(elem)])
    np.savetxt("_tmpcmd.txt", cmds, fmt="%s")
    sh = Popen("Setup_Py_Dir; py_wind {} < _tmpcmd.txt".format(root), stdout=PIPE, stderr=PIPE, shell=True)
    stdout, stderr = sh.communicate()

    # if stderr:
    #     print(stderr.decode("utf-8"))

    Utils.remove_data_sym_links("./")

    return stdout.decode("utf-8")


def plot_cell_sed(model_bands: List[str], filename: str, icell: int, jcell: int, smooth: int = 30) -> None:
    """
    Create a plot of a cell SED given the model bands.

    Parameters
    ----------
    model_bands: List[str]
        A list containing the spectral model bands output from py_wind.
    filename: str
        The file name of the plot.
    """

    numin = []
    bandmin = []
    numax = []
    bandmax = []
    model = []
    pl_log_w = []
    pl_alpha = []
    exp_w = []
    exp_temp = []

    for line in model_bands:
        data = line.split()
        numin.append(float(data[1]))
        bandmin.append(float(data[3][:-1]))
        numax.append(float(data[5]))
        bandmax.append(float(data[7][:-1]))
        model.append(int(data[9]))
        pl_log_w.append(float(data[11]))
        pl_alpha.append(float(data[13]))
        exp_w.append(float(data[15]))
        exp_temp.append(float(data[17]))

    freq = []
    f_nu = []

    for i in range(len(numin)):
        if numax[i] > numin[i]:
            freq_temp = np.logspace(np.log10(numin[i]), np.log10(numax[i]), 101)
            for nu in freq_temp:
                freq.append(nu)
                if model[i] == 1:
                    f_nu.append(10 ** (pl_log_w[i] + np.log10(nu) * pl_alpha[i]))
                elif model[i] == 2:
                    f_nu.append(
                        exp_w[i] * np.exp((-1.0 * consts.h.cgs.value * nu) / (exp_temp[i] * consts.k_B.cgs.value)))
                else:
                    f_nu.append(0.0)
        else:
            freq.append(bandmin[i])
            freq.append(bandmax[i])
            f_nu.append(0.0)
            f_nu.append(0.0)

    xi = 0.0
    for i in range(len(freq) - 1):
        if freq[i] > 3.288e15:
            xi = xi + ((f_nu[i + 1] + f_nu[i]) / 2.0) * (freq[i + 1] - freq[i])

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    ax.loglog(freq, f_nu, "k-", linewidth=3)
    ax.loglog(freq, smooth_spectrum(np.array(f_nu, dtype=np.float64), smooth), "r-", label="smoothed")
    ax.set_xlim(np.min(freq), np.max(freq))
    ax.set_xlabel(r"Frequency, $\nu$")
    ax.set_ylabel(r"$J_{\nu}$ in cell (ergs s$^{-1}$ cm$^{-3}$ Sr$^{-1}$ Hz$^{-1}$)")
    ax.set_title("Cell SED i = {} j = {}".format(icell, jcell))
    ax.legend(loc="lower left")

    # Add axes labels for photon energy
    axx = ax.twiny()
    mn, mx = ax.get_xlim()
    hztoenergy = consts.h.value / consts.e.value
    axx.set_xlim(mn * hztoenergy, mx * hztoenergy)
    axx.set_xlabel("Energy, eV")
    axx.set_xscale("log")

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    fig.savefig(filename + '.png')
    plt.show(fig)

    return


def main() -> None:
    """
    Main function of the script. Parses arguments from the command line.
    """

    p = ap.ArgumentParser(description=__doc__)
    p.add_argument("root", help="The root name of the Python simulation")
    p.add_argument("nx", type=int, help="The number of cells in the i-direction")
    p.add_argument("nz", type=int, help="The number of cells in the j-direction")
    p.add_argument("i", type=int, help="The i index of the cell")
    p.add_argument("j", type=int, help="The j index of the cell")
    args = p.parse_args()

    model = get_spec_model(args.root, args.nx, args.nz, args.i, args.j)
    filename = "cell_i{}_j{}".format(args.i, args.j)
    plot_cell_sed(model, filename, args.i, args.j)


if __name__ == "__main__":
    main()
