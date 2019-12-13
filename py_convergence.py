#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The purpose of this script is to determine if a simulation has converged or not,
and to create any plots which are related to the convergence of a simulation.
"""

import numpy as np
from matplotlib import pyplot as plt
from PyPython import Simulation
from PyPython import Utils
from PyPython import Quotes
from typing import List


COL_WIDTH = 100


def plot_convergence(root: str, convergence: List[float], converging: List[float] = None, tr: List[float] = None,
                    te: List[float] = None, te_max: List[float] = None, hc: List[float] = None,  wd: str = "./"):
    """
    Create a detailed plot of the convergence of a Python simulation, including,
    if provided, a breakdown of the different convergence criteria. 

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    convergence: List[float]
        The convergence fraction of the simulation for each cycle.
    converging: List[float] [optional]
        The converging fraction of the simulation for each cycle.
    tr: List[float] [optional]
        The fraction of cells which have converged radiation temperature.
    te: List[float] [optional]
        The fraction of cells which have converged electron temperature.
    te_max: List[float] [optional]
        The fraction of cells which hit the electron temperature limit.
    hc: List[float] [optional]
        The fraction of cells which have converged heating and cooling rates.
    wd: str [optional]
        The directory containing the Python simulation.
    """

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ncycles = len(convergence)
    cycles = np.arange(1, ncycles + 1, 1)
    ax.set_xlim(1, ncycles)
    ax.set_ylim(0, 1)
    ax.set_xticks(cycles[::2])

    ax.plot(cycles, convergence, label="Convergence")
    if converging:
        ax.plot(cycles, converging, label="Converging")
    if tr:
        ax.plot(cycles, tr, "--", label="Radiation temperature")
    if te:
        ax.plot(cycles, te, "--", label="Electron temperature")
    if te_max:
        ax.plot(cycles, te_max, "--", label="Electron temperature max")
    if hc:
        ax.plot(cycles, hc, "--", label="Heating/Cooling")

    ax.legend()
    ax.set_xlabel("Cycle")
    ax.set_ylabel("Fraction of Cells Passed")
    ax.set_title("Final Convergence = {:4.2f}%".format(float(convergence[-1]) * 100))
    fig.tight_layout(rect=[0.015, 0.015, 0.985, 0.985])
    plt.savefig("{}/{}_convergence.png".format(wd, root))
    plt.close()

    return


def get_convergence(root: str, wd: str = "./") -> None:
    """
    Print out the convergence of a Python simulation and then create a detailed
    plot of the convergence and convergence break down of the simulation.

    Parameters
    ----------
    root: str
        The root name of the Python simulation.
    wd: str [optional]
        The directory containing the Python simulation.
    """

    convergence = Simulation.check_convergence(root, wd, return_per_cycle=True)
    converging = Simulation.check_convergence(root, wd, return_per_cycle=True, return_converging=True)
    tr, te, te_max, hc = Simulation.check_convergence_criteria(root, wd)

    ncycles = len(convergence)
    for i in range(ncycles):
        print("Cycle {:2d} / {:2d}: {:5.2f}% of cells converged and {:5.2f}% of cells are still converging"
              .format(i + 1, ncycles - 1, convergence[i] * 100, converging[i] * 100))
    print("")

    plot_convergence(root, convergence, converging, tr, te, te_max, hc, wd)

    return


def main():
    """
    Main function of the script.
    """

    print("-" * COL_WIDTH, "\n")

    Quotes.random_quote()

    pfs = Utils.find_parameter_files()
    for i in range(len(pfs)):
        root, wd = Utils.split_root_directory(pfs[i])
        print("-" * COL_WIDTH)
        print("\nGetting the convergence for {} in directory {}\n".format(root, wd[:-1]))
        get_convergence(root, wd)

    print("-" * COL_WIDTH)

    return


if __name__ == "__main__":
    main()
