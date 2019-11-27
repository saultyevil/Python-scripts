#!/usr/bin/env python

"""
Compare low and original thingies fuck u
"""

import numpy as np
from sys import exit
from platform import system
from matplotlib import pyplot as plt
from matplotlib import gridspec
from PyPython import SpectrumUtils
from PyPython import WindUtils
import matplotlib.colors as colors


if system() == "Darwin":
    home = "/Users"
else:
    home = "/home"
home += "/saultyevil/PySims/tde"

original_home = home + "/paper_models_matrix_pow/"
bug_fix_home = home + "/bug645_reflect/knox_fix"

original_simulations = [
    "clump/1e-1/agn/solar",
    "clump/1e-1/cv/solar",
    "clump/1e-1/cv/cno",
    "clump/1e-1/spherical/solar",
    "smooth/agn/solar",
    "smooth/cv/solar",
    "smooth/cv/cno",
    "smooth/spherical/solar"
]

bug_simulations = [
    "agn_clump_0.1",
    "cv_clump_0.1_solar",
    "cv_clump_0.1_cno",
    "spherical_clump_0.1",
    "agn_solar_1",
    "cv_solar_1",
    "cv_cno_1",
    "spherical_solar_1"
]


def wind_plot(fig, ax, x, z, var, title, projection):
    """Plot the wind"""
    if projection == "rectilinear":
        if title[:11] != "convergence":
            if title.find("log10(c4)") == -1:
                im = plt.pcolor(np.log10(x), np.log10(z), np.log10(var), vmin=0, vmax=5)
            else:
                im = plt.pcolor(np.log10(x), np.log10(z), np.log10(var), vmin=-10, vmax=0)
        else:
            im = plt.pcolor(np.log10(x), np.log10(z), var, vmin=0, vmax=3)
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


def create_plot(root, original_fname, bugfix_fname, inclination, smooth, projection, wmin, wmax, ofname):
    """Create comparison plots"""

    pdims = (3, 4)
    fig = plt.figure(figsize=(25, 15))
    gridspec.GridSpec(pdims[0], pdims[1])

    original_fname += "/"
    bugfix_fname += "/"

    # Plot the spectrum comparison
    try:
        original_spec = SpectrumUtils.read_spec(original_fname + root + ".spec")
    except IOError:
        print("Can't open file {}".format(original_fname + root + ".spec"))
        return
    try:
        bug_fix_spec = SpectrumUtils.read_spec(bugfix_fname + root + ".spec")
    except IOError:
        print("Can't open file {}".format(bugfix_fname + root + ".spec"))
        return

    orig_wl = original_spec["Lambda"].values.astype(float)
    orig_fl = SpectrumUtils.smooth_spectrum(original_spec[inclination[0]].values.astype(float), smooth)
    fix_wl = bug_fix_spec["Lambda"].values.astype(float)
    fix_fl = SpectrumUtils.smooth_spectrum(bug_fix_spec[inclination[0]].values.astype(float), smooth)
    ax = plt.subplot2grid(pdims, (0, 0), colspan=2)
    plt.semilogy(orig_wl, orig_fl, label="original")
    plt.semilogy(fix_wl, fix_fl, "--", label="bug fix")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(orig_wl, orig_fl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.legend(loc="lower right")

    orig_wl = original_spec["Lambda"].values.astype(float)
    orig_fl = SpectrumUtils.smooth_spectrum(original_spec[inclination[1]].values.astype(float), smooth)
    fix_wl = bug_fix_spec["Lambda"].values.astype(float)
    fix_fl = SpectrumUtils.smooth_spectrum(bug_fix_spec[inclination[1]].values.astype(float), smooth)
    ax = plt.subplot2grid(pdims, (0, 2), colspan=2)
    plt.semilogy(orig_wl, orig_fl, label="original")
    plt.semilogy(fix_wl, fix_fl, "--", label="bug fix")
    plt.xlabel(r"Wavelength [$\AA$]")
    plt.ylabel(r"Flux $F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    plt.xlim(wmin, wmax)
    yup, ylo = SpectrumUtils.ylims(orig_wl, orig_fl, wmin, wmax)
    plt.ylim(ylo, yup)
    SpectrumUtils.plot_line_ids(ax, SpectrumUtils.common_lines())
    plt.legend(loc="lower right")

    orig_input = "{}/{}.0.master.txt".format(original_fname, root)
    fix_input = "{}/{}.0.master.txt".format(bugfix_fname, root)

    # Plot electron temperature comparison
    x, z, te = WindUtils.extract_wind_var(root, "t_e", "wind", original_fname, projection, input_file=orig_input)
    ax = plt.subplot2grid(pdims, (1, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) original", projection)
    x, z, te = WindUtils.extract_wind_var(root, "t_e", "wind", bugfix_fname, projection, input_file=fix_input)
    ax = plt.subplot2grid(pdims, (2, 0), projection=projection)
    wind_plot(fig, ax, x, z, te, "log10(te) bug fix", projection)

    # Plot ionisation parameter comparison
    x, z, ip = WindUtils.extract_wind_var(root, "ip", "wind", original_fname, projection, input_file=orig_input)
    ax = plt.subplot2grid(pdims, (1, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) original", projection)
    x, z, ip = WindUtils.extract_wind_var(root, "ip", "wind", bugfix_fname, projection, input_file=fix_input)
    ax = plt.subplot2grid(pdims, (2, 1), projection=projection)
    wind_plot(fig, ax, x, z, ip, "log10(ip) bug fix", projection)

    # Plot C_i04 parameter comparison
    x, z, c4 = WindUtils.extract_wind_var(root, "c4", "wind", original_fname, projection, input_file=orig_input)
    ax = plt.subplot2grid(pdims, (1, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) original", projection)
    x, z, c4 = WindUtils.extract_wind_var(root, "c4", "wind", bugfix_fname, projection, input_file=fix_input)
    ax = plt.subplot2grid(pdims, (2, 2), projection=projection)
    wind_plot(fig, ax, x, z, c4, "log10(c4) bug fix", projection)

    # Plot convergence
    x, z, c = WindUtils.extract_wind_var(root, "converge", "wind", original_fname, projection, input_file=orig_input)
    ax = plt.subplot2grid(pdims, (1, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence original", projection)
    x, z, c = WindUtils.extract_wind_var(root, "converge", "wind", bugfix_fname, projection, input_file=fix_input)
    ax = plt.subplot2grid(pdims, (2, 3), projection=projection)
    wind_plot(fig, ax, x, z, c, "convergence bug fix", projection)

    fig.tight_layout()
    plt.savefig(ofname + ".png")
    plt.close()

    return


if __name__ == "__main__":

    wmin = 1000
    wmax = 3000
    smooth = 5

    nsims = 8
    roots = ["tde_agn", "tde_cv", "tde_cv", "tde_spherical", "tde_agn", "tde_cv", "tde_cv", "tde_spherical"]
    inclinations = [["75", "85"], ["60", "85"], ["60", "85"], ["60", "85"], ["75", "85"], ["60", "85"], ["60", "85"], ["60", "85"]]
    projections = ["rectilinear", "rectilinear", "rectilinear", "polar", "rectilinear", "rectilinear", "rectilinear", "polar"]
    outnames = ["agn_clump", "cv_solar_clump", "cv_cno_clump", "spherical_clump", "agn_smooth", "cv_solar_smooth", "cv_cno_smooth", "spherical_smooth"]

    for i in range(nsims):
        original_simulations[i] = original_home + "/" + original_simulations[i]
        bug_simulations[i] = bug_fix_home + "/" + bug_simulations[i]
        create_plot(roots[i], original_simulations[i], bug_simulations[i], inclinations[i], smooth, projections[i], wmin, wmax, outnames[i])
