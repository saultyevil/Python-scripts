#!/usr/bin/env python3


from typing import List, Tuple
import numpy as np
from PyPython import WindUtils
from matplotlib import pyplot as plt
from consts import *
from copy import copy

SHOW_PLOT = True

plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15


def sightline_coords(x: np.ndarray, theta: float):
    """
    Return the vertical coordinates for a sightline given the x coordinates
    and the inclination of the sightline.

    Parameters
    ----------
    x: np.ndarray[float]
        The x-coordinates of the sightline
    theta: float
        The opening angle of the sightline

    Returns
    -------
    z: np.ndarray[float]
        The z-coordinates of the sightline
    """

    return x * np.tan(np.pi / 2 - theta)


EPS = 1e-10


def renorm(a, scalar):
    """Renormalise the vector because we must do so"""

    x = np.dot(a, a)

    if x < EPS:
        print("renorm: Cannot renormalize a vector of length 0\n")
        print("renorm: %e %e %e\n", a[0], a[1], a[2])
        print("renorm returning -1\n")
        return -1

    x = scalar / np.sqrt(x)
    a[0] *= x
    a[1] *= x
    a[2] *= x

    return a


def project_from_xyz_cyl(a, b):
    """Attempt to convert something from cartesian to cylindrical coordinates."""

    result = np.zeros(3)
    n_rho = np.zeros(3)
    n_z = np.zeros(3)

    n_rho[0] = a[0]
    n_rho[1] = a[1]
    n_rho[2] = 0
    rc = renorm(n_rho, 1.0)
    if type(rc) == int:
        return rc

    n_z[0] = n_z[1] = 0
    n_z[2] = 1

    n_phi = np.cross(n_z, n_rho)

    result[0] = np.dot(b, n_rho)
    result[1] = np.dot(b, n_phi)
    result[2] = b[2]

    return result


M_bh = 3e7 * MSOL


def rotational_velocity(r0, r):
    """Return the value of the rotational velocity at a given point r in km/s."""

    vk = np.sqrt(G * M_bh / r0)

    return vk * r0 / r


def wind_plot(root: str, output_name: str, wind_names: List[str], wind_types: List[str], wd: str = "./",
              input_file: str = None, projection: str = "rectilinear", scale: str = "loglog",
              line_of_sights: bool = None, plot_indices: bool = False, ndims: str = "2d",
              subplot_dims: Tuple[int, int] = None, fig_size: Tuple[int, int] = (10, 15), plot_title: str = None,
              filetype: str = "png", verbose: bool = False) -> None:
    """
    Creates a figure containing 2D plots of each provided wind variable. For
    each variable in wind_names, there should be the corresponding variable
    type which can either be "wind" or "ion". The subplot dimensions should
    multiply together to be larger than the number of variables provided to
    the function.
    """

    n = wind_plot.__name__
    allowed_projections = ["rectilinear", "polar"]

    if ndims.lower() != "2d":
        raise NotImplementedError("{}: only understand ndims 2d at the moment".format(n))

    # Check to make sure each wind variable has a type corresponding type
    if len(wind_names) != len(wind_types):
        print("{}: wind_names and wind_types should be the same length".format(n))
        return

    # Check to ensure the subplot dimensions are provided and valid
    if subplot_dims is None:
        subplot_dims = (2, 2)
    if subplot_dims[0] * subplot_dims[1] < len(wind_names):
        print("{}: not enough subplot panels to plot all the provided wind variables".format(n))
        return

    # Check to see if the projection is understood and set up the figure and axes objects
    if projection not in allowed_projections:
        print("{}: projection {} not allowed, allowed values: rectilinear or polar".format(n, projection))
        return

    fig, ax = plt.subplots(subplot_dims[0], subplot_dims[1], figsize=fig_size, squeeze=False,
                           sharex="col", sharey="row")

    incls = ["10", "60", "75"]
    lstyle = ["k-", "k--", "k-.", "k:"]
    namezzz = [r"$\log_{10}$(Electron Temperature) [K]", r"$\log_{10}$(Electron Density) [cm$^{-3}$]",
               r"$\log_{10}$(Carbon IV Density) [cm$^{-3}$]", r"$\log_{10}$(Polodial Velocity) [km s$^{-1}$]"]

    if plot_indices:
        scale = "linlin"

    index = 0
    for i in range(subplot_dims[0]):
        for j in range(subplot_dims[1]):
            if index > len(wind_names) - 1:
                break

            wind_name = wind_names[index]
            wind_type = wind_types[index]

            if wind_name != "v_l":
                try:
                    x, z, w = WindUtils.get_wind_variable(root, wind_name, wind_type, wd, projection,
                                                         input_file=input_file, return_indices=plot_indices)
                except Exception:
                    print("Exception occured >:(. Can't plot {} for some reason".format(wind_name))
                    index += 1
                    continue

            r0 = 2.65e13

            with np.errstate(divide="ignore"):
                if wind_name == "v_x":
                    vk = copy(w)
                    dims = w.shape
                    for ii in range(dims[0]):
                        for jj in range(dims[1]):
                            if type(vk[ii, jj]) == np.float64:
                                vk[ii, jj] = rotational_velocity(r0, x[ii, jj]) / 1e5
                    im = ax[i, j].pcolor(x, z, np.log10(vk))
                elif wind_name == "v_l" or wind_name == "v_x":
                    xvx, zvx, wvx = WindUtils.get_wind_variable(root, "v_x", wind_type, wd, projection,
                                                               input_file=input_file, return_indices=plot_indices)
                    xvy, zvy, wvy = WindUtils.get_wind_variable(root, "v_y", wind_type, wd, projection,
                                                               input_file=input_file, return_indices=plot_indices)
                    xvz, zvz, wvz = WindUtils.get_wind_variable(root, "v_z", wind_type, wd, projection,
                                                               input_file=input_file, return_indices=plot_indices)
                    x = copy(xvx)
                    z = copy(zvx)
                    n1, n2 = wvx.shape
                    vl = np.zeros_like(wvx)
                    for ii in range(n1):
                        for jj in range(n2):
                            r = [xvx[ii, jj], 0, zvx[ii, jj]]
                            if type(wvx[ii, jj]) != np.float64 or type(wvy[ii, jj]) != np.float64 \
                                    or type(wvz[ii, jj]) != np.float64:
                                vl[ii, jj] = 0
                            else:
                                v = [wvx[ii, jj], wvy[ii, jj], wvz[ii, jj]]
                                v_cyl = project_from_xyz_cyl(r, v)
                                if type(v_cyl) == int:
                                    continue
                                vl[ii, jj] = np.sqrt(v_cyl[0] ** 2 + v_cyl[2] ** 2) / 1e5
                    if wind_name == "v_l":
                        im = ax[i, j].pcolor(xvx, zvx, np.log10(vl))
                    elif wind_name == "v_x":
                        im = ax[i, j].pcolor(xvx, zvx, np.log10(v_cyl[1] / 1e5))
                else:
                    if wind_name == "c4":
                        x, z, w = WindUtils.get_wind_variable(root, "i04", "ion", wd, projection,
                                                             input_file="/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar/tde_cv.0.C.txt",
                                                             return_indices=plot_indices)
                        im = ax[i, j].pcolor(x, z, np.log10(w)) # , vmin=-10, vmax=0)
                    else:
                        im = ax[i, j].pcolor(x, z, np.log10(w))

            print(wind_name, wind_type)

            if i == 0 and j == 0:
                for k in range(len(incls)):
                    xsight = np.linspace(0, np.max(x), int(1e5))
                    zsight = sightline_coords(xsight, np.deg2rad(float(incls[k])))
                    ax[i, j].plot(xsight, zsight, lstyle[k], label=incls[k] + r"$^{\circ}$ sightline")

            fig.colorbar(im, ax=ax[i, j])  # , orientation="horizontal")

            ax[i, j].set_xlim(np.min(x[x != 0]), np.max(x))
            ax[i, j].set_ylim(np.min(z[z != 0]), np.max(z))

            if scale == "loglog" or scale == "logx":
                ax[i, j].set_xscale("log")
            if scale == "loglog" or scale == "logy":
                ax[i, j].set_yscale("log")

            if i == 0 and j == 0:
                ax[i, j].legend(loc="lower right")

            ax[i, j].text(0.03, 0.93, namezzz[index], ha="left", va="center", rotation="horizontal", fontsize=15,
                          transform=ax[i, j].transAxes)

            index += 1

    fig.tight_layout(rect=[0.03, 0.03, 0.97, 0.97])
    # fig.subplots_adjust(hspace=0)

    # fig.subplots_adjust(right=0.825, wspace=0, hspace=0)
    # cax = fig.add_axes([0.85, 0.06, 0.035, 0.91])
    # fig.colorbar(im, cax=cax)

    fig.text(0.5, 0.02, r"x [cm]", ha="center", va="center", rotation="horizontal", fontsize=17)
    fig.text(0.025, 0.5, r"z [cm]", ha="center", va="center", rotation="vertical", fontsize=17)

    if plot_title:
        fig.suptitle(plot_title)
    # fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig("clumpy_wind_properties.png".format(wd, output_name, filetype))

    if SHOW_PLOT:
        plt.show(block=True)
    else:
        plt.close()

    return


root = "tde_cv"
output_name = "wind"
projection = "rectilinear"
path = "/home/saultyevil/PySims/tde_uv/models/clump/1e-1/cv/solar"
inclination_angles = None
DIMS = "2d"

wind_names = ["t_e", "ne", "c4", "v_l"]
wind_types = ["wind"] * len(wind_names)
wind_plot(root, output_name + "_wind", wind_names, wind_types, path, projection=projection,
          line_of_sights=inclination_angles, ndims=DIMS, verbose=True, subplot_dims=(2, 2),
          fig_size=(17, 15))
