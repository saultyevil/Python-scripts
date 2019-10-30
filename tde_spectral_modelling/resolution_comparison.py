#!/usr/bin/env python

"""
Compare low and high resolution thingies fuck u
"""

import numpy as np
from sys import exit
from platform import system
from matplotlib import pyplot as plt
from matplotlib import gridspec
from PyPython import SpectrumUtils
from PyPython import WindUtils
from PyPython import Simulation
from typing import Union


if system() == "Darwin":
    home = "/Users"
else:
    home = "/home"
home += "/saultyevil/PySims/tde"
lo_res_home = home + "/low_resolution/paper_models/"
hi_res_home = home + "/paper_models/"


def wind_plot(fig: plt.Figure, ax: plt.Axes, x: np.ndarray, z: np.ndarray, var: np.ndarray,
              title: str, projection: str) -> Union[plt.Figure, plt.Axes]:
    """
    Create the wind plot...

    Parameters
    ----------
    fig: plt.figure
        The figure object
    ax: plt.Axes
        The plot object
    x: np.ndarray[float]
        The x or r coordinates
    z: np.ndarray[float]
        The z or theta coordinates
    var: np.array[float]
        The variable to plot
    title: str
        The title of the plot
    projection
        The projection of the plot

    Returns
    -------
    fig: plt.figure
        The figure object
    ax: plt.Axes
        The plot object
    """

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


def create_plot(root: str, wd_high_res: str, wd_low_res: str, inclination: str, smooth: int, projection: str,
                wmin: float, wmax: float, ofname: str) -> None:
    """
    Create comparison plots

    Parameters
    ----------
    root: str
        The root name of both simulations.
    wd_high_res: str
        The directory containing the high resolution simulation.
    wd_low_res: str
        The directory containing the low resolution simulation.
    inclination: str
        The inclination to use for the spectrum.
    smooth: int
        The amount of smoothing for the spectrum.
    projection: str
        The coordinate system of the simulation: rectilinear or polar.
    wmin: float
        Smallest wavelength to plot
    wmax: float
        Largest wavelength to plot
    ofname: str
        The name of the output file
    """

    pdims = (3, 4)
    fig = plt.figure(figsize=(25, 15))
    gridspec.GridSpec(pdims[0], pdims[1])

    wd_high_res += "/"
    wd_low_res += "/"

    # Plot the spectrum comparison
    try:
        hispec = SpectrumUtils.read_spec(wd_high_res + root + ".spec")
    except IOError:
        print("Can't open file {}".format(wd_high_res + root + ".spec"))
        return
    try:
        lospec = SpectrumUtils.read_spec(wd_low_res + root + ".spec")
    except IOError:
        print("Can't open file {}".format(wd_low_res + root + ".spec"))
        return
    
    hi_conv = Simulation.check_convergence(root, wd_high_res)
    lo_conv = Simulation.check_convergence(root, wd_low_res)

    hiwl = hispec["Lambda"].values.astype(float)
    hifl = SpectrumUtils.smooth_spectrum(hispec[inclination[0]].values.astype(float), smooth)
    lowl = lospec["Lambda"].values.astype(float)
    lofl = SpectrumUtils.smooth_spectrum(lospec[inclination[0]].values.astype(float), smooth)
    ax = plt.subplot2grid(pdims, (0, 0), colspan=2)
    plt.semilogy(hiwl, hifl, label="High resolution")
    plt.semilogy(lowl, lofl, "--", label="low resolution")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(hiwl, hifl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.title("High resolution convergence = {}".format(hi_conv))
    plt.legend(loc="lower right")

    hiwl = hispec["Lambda"].values.astype(float)
    hifl = SpectrumUtils.smooth_spectrum(hispec[inclination[1]].values.astype(float), smooth)
    lowl = lospec["Lambda"].values.astype(float)
    lofl = SpectrumUtils.smooth_spectrum(lospec[inclination[1]].values.astype(float), smooth)
    ax = plt.subplot2grid(pdims, (0, 2), colspan=2)
    plt.semilogy(hiwl, hifl, label="High resolution")
    plt.semilogy(lowl, lofl, "--", label="low resolution")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(hiwl, hifl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.title("Low resolution convergence = {}".format(lo_conv))
    plt.legend(loc="lower right")

    hi_input = "{}/{}.0.master.txt".format(wd_high_res, root)
    lo_input = "{}/{}.0.master.txt".format(wd_low_res, root)

    # Plot electron temperature comparison
    x, z, te = WindUtils.extract_wind_var(root, "t_e", "wind", wd_high_res, projection, input_file=hi_input)
    ax = plt.subplot2grid(pdims, (1, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) high res", projection)
    x, z, te = WindUtils.extract_wind_var(root, "t_e", "wind", wd_low_res, projection, input_file=lo_input)
    ax = plt.subplot2grid(pdims, (2, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) low res", projection)

    # Plot ionisation parameter comparison
    x, z, ip = WindUtils.extract_wind_var(root, "ip", "wind", wd_high_res, projection, input_file=hi_input)
    ax = plt.subplot2grid(pdims, (1, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) high res", projection)
    x, z, ip = WindUtils.extract_wind_var(root, "ip", "wind", wd_low_res, projection, input_file=lo_input)
    ax = plt.subplot2grid(pdims, (2, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) low res", projection)

    # Plot C_i04 parameter comparison
    x, z, c4 = WindUtils.extract_wind_var(root, "c4", "wind", wd_high_res, projection, input_file=hi_input)
    ax = plt.subplot2grid(pdims, (1, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) high res", projection)
    x, z, c4 = WindUtils.extract_wind_var(root, "c4", "wind", wd_low_res, projection, input_file=lo_input)
    ax = plt.subplot2grid(pdims, (2, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) low res", projection)

    # Plot convergence
    x, z, c = WindUtils.extract_wind_var(root, "converge", "wind", wd_high_res, projection, input_file=hi_input)
    ax = plt.subplot2grid(pdims, (1, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence high res", projection)
    x, z, c = WindUtils.extract_wind_var(root, "converge", "wind", wd_low_res, projection, input_file=lo_input)
    ax = plt.subplot2grid(pdims, (2, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence low res", projection)

    fig.tight_layout()
    plt.savefig(ofname + ".png")
    plt.close()

    return


if __name__ == "__main__":

    wmin = 1000
    wmax = 3000
    smooth = 5

    simulations = [
        "smooth/agn/solar", "smooth/agn/cno", "smooth/cv/solar", "smooth/cv/cno", "smooth/spherical/solar", 
        "smooth/spherical/cno", "clump/1e-1/agn/solar", "clump/1e-1/agn/cno", "clump/1e-1/cv/solar", 
        "clump/1e-1/cv/cno", "clump/1e-1/spherical/solar", "clump/1e-1/spherical/cno"
    ]
    roots = [
        "tde_agn", "tde_agn", "tde_cv", "tde_cv", "tde_spherical", "tde_spherical", "tde_agn", "tde_agn", "tde_cv",
        "tde_cv", "tde_spherical", "tde_spherical"
    ]
    inclinations = [
        ["75", "85"], ["75", "85"], ["60", "85"], ["60", "85"], ["60", "85"], ["60", "85"], ["75", "85"], 
        ["75", "85"], ["60", "85"], ["60", "85"], ["60", "85"], ["60", "85"]
    ]

    hires = []
    lores = []
    projections = []
    outnames = []
    for i in range(len(simulations)):
        hires.append(hi_res_home + simulations[i])
        lores.append(lo_res_home + simulations[i])
        if simulations[i].find("spherical") != -1:
            projections.append("polar")
        else:
            projections.append("rectilinear")
        if simulations[i].find("cno") != -1:
            abun = "cno"
        elif simulations[i].find("solar") != -1:
            abun = "solar"
        else:
            abun = i
        if simulations[i].find("clump") != -1:
            wind = "clump"
        elif simulations[i].find("smooth") != -1:
            wind = "smooth"
        else:
            wind = i
        outnames.append("{}_{}_{}".format(roots[i], wind, abun))

    for i in range(len(simulations)):
        create_plot(roots[i], hires[i], lores[i], inclinations[i], smooth, projections[i], wmin, wmax, outnames[i])
