#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import py_util
import numpy as np
from matplotlib import pyplot as plt


SMOOTH = 11
VERBOSE = False
PLT_SHOW = False
PLT_BLAG = False


def plot_3_inclinations(files, titles, suptitle, inds, labels, ncols, nrows, outname, ylims):
    """
    Plot things
    """

    print(suptitle)

    fig, ax = plt.subplots(nrows, ncols, figsize=(8, 12), squeeze=False)

    if PLT_BLAG:
        blag_spec = py_util.get_blagordnova_spec(SMOOTH, VERBOSE)

    file = 0
    for i in range(nrows):
        for j in range(ncols):
            if file >= len(files):
                break

            if PLT_BLAG:
                Z = 0.07897
                ax[i, j].semilogy(blag_spec[:, 0] / (Z + 1), blag_spec[:, 1])

            spec = py_util.read_spec_file(files[file])
            wavelength = np.array(spec[1:, 1], dtype=float)
            flux1 = np.array(spec[1:, inds[0]], dtype=float)
            flux1 = py_util.smooth_spectra(flux1, SMOOTH, VERBOSE)
            fmax = np.max(flux1)
            flux1 /= fmax

            flux2 = np.array(spec[1:, inds[1]], dtype=float)
            flux2 = py_util.smooth_spectra(flux2, SMOOTH, VERBOSE)
            fmax = np.max(flux2)
            flux2 /= fmax
            flux2 *= 100


            flux3 = np.array(spec[1:, inds[2]], dtype=float)
            flux3 = py_util.smooth_spectra(flux3, SMOOTH, VERBOSE)
            fmax = np.max(flux3)
            flux3 /= fmax
            flux3 *= 10000

            if PLT_BLAG:
                const = (100 * 3.08567758128e18)**2 / 1.079987153448e+27**2
                flux1 *= const
                flux2 *= const
                flux3 *= const

            ax[i, j].semilogy(wavelength, flux1, label=labels[0])
            ax[i, j].semilogy(wavelength, flux2, label=labels[1])
            ax[i, j].semilogy(wavelength, flux3, label=labels[2])
            ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")

            if ylims == "cv_macro":
                ax[i, j].set_xlim(500, 3000)

            if PLT_BLAG:
                ax[i, j].set_ylim(1e-17, 1e-13)
            else:
                ax[i, j].set_ylim(1e-4, 1e6)

            # elif ylims == "agn":
            #     ax[i, j].set_ylim(1e-5, 1)
            # elif ylims == "cv":
            #     ax[i, j].set_ylim(1e-4, 1e-1)

            ax[i, j].set_title(titles[file])
            ax[i, j].legend(loc="lower right")
            file += 1

    fig.suptitle(suptitle)
    fig.tight_layout()
    if PLT_BLAG:
        outname = "blag_" + outname
    plt.savefig(outname)

    print("Plotted {} specs".format(file))
    if PLT_SHOW:
        plt.show()
    else:
        plt.close()

    return

def plot_agn_grid():
    """
    Plot grid of AGN wind models
    """

    inds = [
        12,
        13,
        14
    ]

    legend = [
        "70 deg",
        "75 deg",
        "80 deg"
    ]

    #
    # mdot_wind / mdot_acc = 1
    #

    files = [
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_0/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_0/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_2/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_2/mratio_1/1-30R/tde.spec",

    ]

    titles = [
        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = 2",

    ]

    outname = "agn_plots_mratio1.pdf"
    suptitle = "Mdot_wind / Mdot_acc = 1.0"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 2, outname, "agn")

    files = [
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_-3/m_ratio1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_-2/m_ratio1/1-30R/tde.spec",
    ]

    titles = [
        "Mdot_acc = 2e-1 mdot_exp = -3",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 2, outname, "agn")

    #
    # mdot_wind / mdot_acc = 0.1
    #

    files = [
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_0/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_0/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_2/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_2/mratio_0.1/1-30R/tde.spec",
        # "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-1/Mdot_exp_-3/mratio_0.1/1-30R/tde.spec",
        # the above never appears to have been run
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/disk_mdot_2e-3/Mdot_exp_-2/mratio_0.1/1-30R/tde.spec",
    ]

    titles = [
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 2",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        # "Mdot_acc = 2e-1 mdot_exp = -3",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    suptitle = "Mdot_wind / Mdot_acc = 0.1"
    outname = "agn_plots_mratio0.1.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")


    #
    # disk mdot = 2e-2
    #

    files = [
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/mdot_2e-2_bigger_disk/mratio_0.1/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/mdot_2e-2_bigger_disk/mratio_0.5/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/AGN_Model/mdot_2e-2_bigger_disk/mratio_1/tde.spec",
    ]

    titles = [
        "mdot_wind / mdot_acc = 0.1",
        "mdot_wind / mdot_acc = 0.5",
        "mdot_wind / mdot_acc = 1.0",
    ]

    suptitle = "aaa"
    outname = "agn_plots_mdot_acc_2e-2.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 1, 3, outname, "agn")

    return


def plot_cv_grid():
    """
    Plot grid of CV wind models
    """

    inds = [
        11,
        12,
        13,
    ]

    legend = [
        "45 deg",
        "62 deg",
        "80 deg",
    ]

    #
    # Mdot_exp
    #

    files = [

        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_1/tde.spec",
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_2/tde.spec",
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_3/tde.spec",
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_-1/tde.spec",

    ]

    titles = [

        "mdot exp = 1",
        "mdot exp = 2",
        "mdot exp = 3",
        "mdot exp = -1",

    ]

    suptitle = "mdot_acc = 2e-2, mdot_wind / mdot_acc = 1, rmin-rmax 1-30Rstar, disk 1e15"
    outname = "cv_plots_mdot_exp.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 2, outname, "cv")

    files = [
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_-2/tde.spec",
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_-3/tde.spec",
        "/home/saultyevil/PySims/TDE/Plots/CV_mdot/exp_0/tde.spec",
    ]

    titles = [
        "mdot exp = -2",
        "mdot exp = -3",
        "mdot exp = 0",
    ]

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 2, outname, "cv")

    #
    # m_ratios
    #

    files = [
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/CV_Model/disk_mdot_2e-2/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/CV_Model/disk_mdot_2e-2/mratio_0.5/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/CV_Model/disk_mdot_2e-2/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/Wind_1Rstar_models/CV_Model/disk_mdot_2e-2/mratio_1.5/tde.spec",
    ]

    titles = [
        "mdot_wind / mdot_acc = 0.1",
        "mdot_wind / mdot_acc = 0.5",
        "mdot_wind / mdot_acc = 1.0",
        "mdot_wind / mdot_acc = 1.5",
    ]

    suptitle = "mdot_acc = 2e-2, rmin-rmax 1-30R, disk 1e15"
    outname = "cv_plots_mratio.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 1, 4, outname, "cv")

    #
    # sim11 macro
    #

    files = [
        "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/tde_macro/tde.spec",
        "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/sim11_V0_sound_speed/tde.spec",
        "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/tde_macro_clump_0.1/tde.spec",
        "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/tde_macro_clump_0.01/tde.spec",
        "/home/saultyevil/PySims/TDE/Star_model/macro_models/sim11_more/sim11_no_censrc/tde.spec",
    ]

    titles = [
        "baseline",
        "sv.v0 sound speed",
        "clump 0.1",
        "clump 0.01",
        "no central source",
    ]

    suptitle = "mdot_acc = 2e-2, mdot_wind / mdot_acc = 0.25, rmin-rmax 1-18R, disk 5e14"
    outname = "cv_plot_macro_11.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv_macro")

    #
    # Star models
    #

    files = [
        "/home/saultyevil/PySims/TDE/Star_model/simple_atom_models/sim6/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/Star_model/simple_atom_models/sim7/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/Star_model/simple_atom_models/sim8/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/Star_model/simple_atom_models/sim8_maxscat/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/Star_model/simple_atom_models/sim11/tdeCV.spec",
    ]

    titles = [
        "mrat = 1, 5-18R, disk 5e14",
        "mrat = 1, 1-18R, disk 5e14",
        "mrat = 0.25, 1-18R, disk 5e14",
        "mrat = 0.25, 1-18R, disk 5e14, maxscat",
        "mrat = 0.25,  1-18R, disk 5e14, wind 1e17, vinf 1",
    ]

    suptitle = "mdot_acc = 4e-2"
    outname = "cv_plot_star_models.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv_macro")

    return


def main():
    """
    Main controlling function
    """

    plot_agn_grid()
    plot_cv_grid()

    return

if __name__ == "__main__":
    main()