#!/usr/bin/env python

"""
The purpose of this script is to create a comparison figure to showcase the
differences between the two provided models. This script assumes that the two
models are only slightly different. But, in theory, could be used to compare
two models with the same coordinate system. Of course, you can also hack the
script if you want to compare a cylindrical model to a polar/spherical model.

TODO: use argparse to get input
"""

import numpy as np
from sys import exit, argv
from matplotlib import pyplot as plt
from matplotlib import gridspec
from PyPython import SpectrumUtils
from PyPython import WindUtils
from PyPython import Simulation
from typing import Union
import tde_spectra


def wind_plot(fig: plt.Figure, ax: plt.Axes, x: np.ndarray, z: np.ndarray, var: np.ndarray,
              title: str, projection: str) -> Union[plt.Figure, plt.Axes]:

    if projection == "rectilinear":
        if title[:11] != "convergence":
            im = plt.pcolor(np.log10(x), np.log10(z), np.log10(var))
        else:
            im = plt.pcolor(np.log10(x), np.log10(z), var)
        plt.xlabel("x")
        plt.ylabel("z")
        plt.xlim(np.log10(x[1, 1]), np.log10(x[-1, -1]))
        plt.ylim(np.log10(z[1, 1]), np.log10(z[-1, -1]))
    elif projection == "polar":
        if title[:11] != "convergence":
            im = plt.pcolor(z, np.log10(x), np.log10(var))
        else:
            im = plt.pcolor(z, np.log10(x), var)
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_thetamin(0)
        ax.set_thetamax(90)
        ax.set_rlabel_position(90)
        ax.set_ylabel("Log[R]")
        rmin = x[1][0]
        rmax = x[-2][0]
        ax.set_rlim(np.log10(rmin), np.log10(rmax))
    else:
        print("Unknown projection {}".format(projection))
        exit(1)
    fig.colorbar(im)
    ax.set_title(title)

    return fig, ax


def create_plot(root1: str, model1_path: str, root2: str, model2_path: str, inclination: str, smooth: int,
                projection: str, wmin: float, wmax: float) -> None:

    pdims = (3, 4)
    fig = plt.figure(figsize=(25, 15))
    gridspec.GridSpec(pdims[0], pdims[1])

    model1_path += "/"
    model2_path += "/"

    # Plot the spectrum comparison
    try:
        model1_spec = SpectrumUtils.read_spec(model1_path + root1 + ".spec")
    except IOError:
        print("Can't open file {}".format(model1_path + root1 + ".spec"))
        return
    try:
        model2_spec = SpectrumUtils.read_spec(model2_path + root2 + ".spec")
    except IOError:
        print("Can't open file {}".format(model2_path + root2 + ".spec"))
        return

    model1_convergence = Simulation.check_convergence(root1, model1_path)
    model2_convergence = Simulation.check_convergence(root2, model2_path)

    model1_wl = model1_spec["Lambda"].values.astype(float)
    model1_fl = SpectrumUtils.smooth_spectrum(model1_spec[inclination[0]].values.astype(float), smooth)
    model2_wl = model2_spec["Lambda"].values.astype(float)
    model2_fl = SpectrumUtils.smooth_spectrum(model2_spec[inclination[0]].values.astype(float), smooth)

    ipt = tde_spectra.iptf15af_spec(smooth)

    ax = plt.subplot2grid(pdims, (0, 0), colspan=2)
    plt.semilogy(model1_wl, model1_fl, label="Model 1")
    plt.semilogy(model2_wl, model2_fl, "--", label="Model 2")
    plt.semilogy(ipt[:, 0], ipt[:, 1], label="iPTF15af")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(model1_wl, model1_fl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.title("Model 1 convergence = {}".format(model1_convergence))
    plt.legend(loc="lower right")

    model1_wl = model1_spec["Lambda"].values.astype(float)
    model1_fl = SpectrumUtils.smooth_spectrum(model1_spec[inclination[1]].values.astype(float), smooth)
    model2_wl = model2_spec["Lambda"].values.astype(float)
    model2_fl = SpectrumUtils.smooth_spectrum(model2_spec[inclination[1]].values.astype(float), smooth)
    ax = plt.subplot2grid(pdims, (0, 2), colspan=2)
    plt.semilogy(model1_wl, model1_fl, label="Model 1")
    plt.semilogy(model2_wl, model2_fl, "--", label="Model 2")
    plt.semilogy(ipt[:, 0], ipt[:, 1], label="iPTF15af")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(model1_wl, model1_fl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.title("Model 2 convergence = {}".format(model2_convergence))
    plt.legend(loc="lower right")

    model1_input = "{}/{}.0.master.txt".format(model1_path, root1)
    model2_input = "{}/{}.0.master.txt".format(model2_path, root2)

    # Plot electron temperature comparison
    x, z, te = WindUtils.extract_wind_var(root1, "t_e", "wind", model1_path, projection, input_file=model1_input)
    ax = plt.subplot2grid(pdims, (1, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) model 1", projection)
    x, z, te = WindUtils.extract_wind_var(root2, "t_e", "wind", model2_path, projection, input_file=model2_input)
    ax = plt.subplot2grid(pdims, (2, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) model 2", projection)

    # Plot Si IV
    x, z, ip = WindUtils.extract_wind_var(root1, "ip", "wind", model1_path, projection, input_file=model1_input)
    ax = plt.subplot2grid(pdims, (1, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) model 1", projection)
    x, z, ip = WindUtils.extract_wind_var(root2, "ip", "wind", model2_path, projection, input_file=model2_input)
    ax = plt.subplot2grid(pdims, (2, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) model 2", projection)

    # Plot C_i04 parameter comparison
    x, z, c4 = WindUtils.extract_wind_var(root1, "c4", "wind", model1_path, projection, input_file=model1_input)
    ax = plt.subplot2grid(pdims, (1, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) model 1", projection)
    x, z, c4 = WindUtils.extract_wind_var(root2, "c4", "wind", model2_path, projection, input_file=model2_input)
    ax = plt.subplot2grid(pdims, (2, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) model 2", projection)

    # Plot convergence
    x, z, c = WindUtils.extract_wind_var(root1, "converge", "wind", model1_path, projection, input_file=model1_input)
    ax = plt.subplot2grid(pdims, (1, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence model 1", projection)
    x, z, c = WindUtils.extract_wind_var(root2, "converge", "wind", model2_path, projection, input_file=model2_input)
    ax = plt.subplot2grid(pdims, (2, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence model 2", projection)

    fig.tight_layout()
    plt.savefig("{}_{}_comparison.png".format(root1, root2))
    plt.close()

    return


if __name__ == "__main__":

    wmin = 1000
    wmax = 3000
    smooth = 5

    if len(argv) != 5:
        print("Invalid run mode")
        print("py_compare_models.py root1 model1_path root2 model2_path")

    root1 = argv[1]
    model1 = argv[2]
    root2 = argv[3]
    model2 = argv[4]

    create_plot(root1, model1, root2, model2, ["60", "75"], smooth, "rectilinear", wmin, wmax)
