#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import tde_util
import py_plot_util
import numpy as np
from matplotlib import pyplot as plt


SMOOTH = 15
VERBOSE = False
PLT_SHOW = False


def plot_3_inclinations(files, titles, suptitle, inds, labels, ncols, nrows, outname, ylims):
    """
    Plot things
    """

    # print(files, len(files))

    print("Plotting " + suptitle)

    if ncols == 1:
        figsi = (4, 12)
    else:
        figsi = (8, 12)

    fig, ax = plt.subplots(nrows, ncols, figsize=figsi, squeeze=False)

    blag = tde_util.iptf15af_spec(SMOOTH, VERBOSE)
    blag[:, 0] /= (0.07897 + 1)
    fmax = np.max(blag[:, 1])
    blag[:, 1] /= fmax
    
    file = 0
    for i in range(nrows):
        for j in range(ncols):
            # print(file)
            if file >= len(files):
                break
            # print(files[file])
            spec = py_plot_util.read_spec_file(files[file])
            wavelength = np.array(spec[1:, 1], dtype=float)
            
            ax[i, j].semilogy(blag[:, 0], blag[:, 1]*100, label="iPTF15af")

            # Read in the fluxes and normalise and scale appropriately
            # inclination angle 1
            flux1 = np.array(spec[1:, inds[0]], dtype=float)
            flux1 = py_plot_util.smooth_1d_array(flux1, SMOOTH, VERBOSE)
            fmax = np.max(flux1)
            flux1 /= fmax
            flux1 *= 10000
            # inclination angle 2
            flux2 = np.array(spec[1:, inds[1]], dtype=float)
            flux2 = py_plot_util.smooth_1d_array(flux2, SMOOTH, VERBOSE)
            fmax = np.max(flux2)
            flux2 /= fmax
            flux2 *= 100  # move up axis
            # inclination angle 3
            flux3 = np.array(spec[1:, inds[2]], dtype=float)
            flux3 = py_plot_util.smooth_1d_array(flux3, SMOOTH, VERBOSE)
            fmax = np.max(flux3)
            flux3 /= fmax
            flux3 /= 1  # move up axis MORE

            # Finally plot...
            ax[i, j].semilogy(wavelength, flux1, label=labels[0])
            ax[i, j].semilogy(wavelength, flux2, label=labels[1])
            ax[i, j].semilogy(wavelength, flux3, label=labels[2])
            ax[i, j].set_xlabel(r"Wavelength ($\AA$)")
            ax[i, j].set_ylabel(r"$F_{\lambda}$ (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")

            # This may need changing...
            ax[i, j].set_ylim(1e-2, 1e5)
            
            # Set ylims due to inconsistency with spectrum cycles
            # if ylims == "cv_macro":
            #     ax[i, j].set_xlim(500, 3000)
            ax[i, j].set_xlim(800, 3000)

            ax[i, j].set_title(titles[file], fontsize=10)
            ax[i, j].legend(loc="lower right")
            file += 1

    fig.suptitle(suptitle)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(outname)

    print("Plotted {} specs".format(file))
    if PLT_SHOW:
        plt.show()
    else:
        plt.close()

    return


def plot_agn_grid_1_30():
    """
    Plot grid of AGN wind models
    """

    inds = [
        11,
        13,
        15
    ]

    legend = [
        "60 deg",
        "75 deg",
        "85 deg"
    ]

    #
    # mdot_wind / mdot_acc = 1
    #

    files = [
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_0/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_0/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_2/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_2/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_-3/m_ratio1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_-2/m_ratio1/1-30R/tde.spec",
    ]

    titles = [
        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = 2",
        "Mdot_acc = 2e-1 mdot_exp = -3",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    outname = "agn_plots_mratio1.pdf"
    suptitle = "Mdot_wind / Mdot_acc = 1.0, disk 1e15, 1-30R, 70-82degs"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")

    #
    # mdot_wind / mdot_acc = 0.1
    #

    files = [
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_0/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_0/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_2/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_2/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_-2/mratio_0.1/1-30R/tde.spec",
    ]

    titles = [

        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    suptitle = "Mdot_wind / Mdot_acc = 0.1, disk 1e15, 1-30R, 70-82degs"
    outname = "agn_plots_mratio0.1.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")


    #
    # disk mdot = 2e-2
    #

    files = [
        "/home/saultyevil/PySims/TDE/agn_model/mdot_2e-2_bigger_disk/mratio_0.1/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/mdot_2e-2_bigger_disk/mratio_0.5/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/mdot_2e-2_bigger_disk/mratio_1/tde.spec",
    ]

    titles = [
        "mdot_wind / mdot_acc = 0.1",
        "mdot_wind / mdot_acc = 0.5",
        "mdot_wind / mdot_acc = 1.0",
    ]

    suptitle = "mdot_acc 2e-2 disk 1e15 rmin-rmax 1-40R, 70-82degs"
    outname = "agn_plots_mdot_acc_2e-2.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")

    return


def plot_agn_grid_1_20():
    """
    Plot grid of AGN wind models
    """

    inds = [
        11,  # 70 degrees
        13,  # 75 degrees
        15   # 80 degrees
    ]

    legend = [
        "60 deg",
        "75 deg",
        "85 deg"
    ]

    #
    # mdot_wind / mdot_acc = 1
    #

    files = [
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_0/mratio_1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_0/mratio_1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_2/mratio_1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_2/mratio_1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_-3/m_ratio1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_-2/m_ratio1/1-20R/tde.spec",
    ]

    titles = [
        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = 2",
        "Mdot_acc = 2e-1 mdot_exp = -3",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    outname = "agn_plots_mratio1_1_20.pdf"
    suptitle = "Mdot_wind / Mdot_acc = 1.0, disk 1e15, 1-20R, 70-82degs"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")

    #
    # mdot_wind / mdot_acc = 0.1
    #

    files = [
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_0/mratio_0.1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_0/mratio_0.1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_2/mratio_0.1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-1/Mdot_exp_2/mratio_0.1/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/agn_model/disk_mdot_2e-3/Mdot_exp_-2/mratio_0.1/1-20R/tde.spec",
    ]

    titles = [
        "Mdot_acc = 2e-3 mdot_exp = 0",
        "Mdot_acc = 2e-1 mdot_exp = 0",
        "Mdot_acc = 2e-3 mdot_exp = 2",
        "Mdot_acc = 2e-1 mdot_exp = 2",
        "Mdot_acc = 2e-3 mdot_exp = -2",
    ]

    suptitle = "Mdot_wind / Mdot_acc = 0.1, disk 1e15, 1-20R, 70-82degs"
    outname = "agn_plots_mratio0.1_1_20.pdf"

    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")

    return


def plot_agn_grid_runs():
    files = [
        "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/agn_grid_simple/1-20Rwind/tde.spec",
        "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/agn_grid_simple/1-30Rwind/tde.spec",
        "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/agn_grid_simple/1-40Rwind/tde.spec",
        "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/agn_grid_simple/1-20Rwind_wind_1e16/tde.spec",
        "/home/saultyevil/Dropbox/DiskWinds/PySims/TDE/agn_grid_simple/1-40R_mdot_acc_2e-3/tde.spec",
    ]

    titles = [
        "1-20Rwind: mdot_acc 2e-1 mrat 0.1 wind radmax 1e17",
        "1-30Rwind: mdot_acc 2e-1 mrat 0.1",
        "1-40Rwind: mdot_acc 2e-1 mrat 0.1",
        "1-20Rwind_wind_1e16: mdot_acc 2e-1 mrat 0.1 wind radmax 1e16",
        "1-40R_mdot_acc_2e-3: mratio 0.1 wind radmax 1e16",
    ]

    inds = [
        11,
        13,
        15
    ]

    legend = [
        "60 deg",
        "75 deg",
        "85 deg"
    ]

    outname = "agn_gird_plots.pdf"
    suptitle = ""
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "agn")

    return


def plot_cv_grid():
    """
    Plot grid of CV wind models
    """

    inds = [
        10,
        12,
        13,
    ]

    legend = [
        "28 deg",
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
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

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

    outname = "cv_plots_mdot_exp_more.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

    #
    # m_ratios
    #

    files = [
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/mratio_0.1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/mratio_0.5/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/mratio_1/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/mratio_1.5/tde.spec",
    ]

    titles = [
        "mdot_wind / mdot_acc = 0.1",
        "mdot_wind / mdot_acc = 0.5",
        "mdot_wind / mdot_acc = 1.0",
        "mdot_wind / mdot_acc = 1.5",
    ]

    suptitle = "mdot_acc = 2e-2, rmin-rmax 1-30R, disk 1e15"
    outname = "cv_plots_mratio.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

    #
    # wind size
    #

    files = [
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/wind_size/1-20R/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/wind_size/1-30R/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/wind_size/1-50R/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/disk_mdot_2e-2/wind_size/1-70R/tde.spec"
    ]
    
    titles = [
        "1-20R",
        "1-30R",
        "1-50R",
        "1-70R"
    ]

    suptitle = "mdot_acc 2e-3, mwind/macc = 0.1, disk 1e15, NOTE wind size may be bullshit"
    outname = "cv_plot_windsize.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

    files = [
        "/home/saultyevil/PySims/TDE/cv_model/tweak/co_mass_2e7/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/tweak/disk_mdot_2e-1/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/tweak/Rdisk/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_model/tweak/Rdisk2/tde.spec"
        
    ]
    
    titles = [
        "macc = 2e-2, mwin/macc=1, 1-30R, disk 1e15",
        "macc=2e-1, mwin/macc=1, 1-30R, disk 1e15",
        "macc=6e-2, mwin/macc=0.33, 1-30, disk 1e14",
        "macc=3e-2, mwin/macc=1.33, 1-30, disk 4e14"
    ]
    
    suptitle = "Some 'tweaks'"
    outname = "cv_plot_tweaks.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

    return


def star_models():
    """
    Plot grid of star which are basicaly CV wind models
    """
   
    inds = [
       10,
       12,
       13,
   ]

    legend = [
       "28 deg",
       "62 deg",
       "80 deg",
   ]
   
    #
    # sim11 macro
    #

    files = [
        "/home/saultyevil/PySims/TDE/cv_star_model/macro_models/sim11_more/tde_macro/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/macro_models/sim11_more/sim11_V0_sound_speed/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/macro_models/sim11_more/tde_macro_clump_0.1/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/macro_models/sim11_more/tde_macro_clump_0.01/tde.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/macro_models/sim11_more/sim11_no_censrc/tde.spec",
    ]

    titles = [
        "baseline",
        "sv.v0 sound speed",
        "clump 0.1",
        "clump 0.01",
        "no cen src (dimmer when not in arb units)",
    ]

    suptitle = "mdot_acc = 2e-2, mdot_wind / mdot_acc = 0.25, rmin-rmax 1-18R, disk 5e14"
    outname = "cv_plot_macro_11.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv_macro")

    #
    # simple atom models
    #

    files = [
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim6/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim7/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim8/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim8_maxscat/tdeCV.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim10/tdeStar.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim11/tdeCV.spec",
    ]

    titles = [
        "mrat = 1, 5-18R, disk 5e14",
        "mrat = 1, 1-18R, disk 5e14",
        "mrat = 0.25, 1-18R, disk 5e14",
        "mrat = 0.25, 1-18R, disk 5e14, maxscat",
        "mrat = 0.25, 7-55R, disk 1e14, co.radius event horizon",
        "mrat = 0.25,  1-18R, disk 5e14, wind 1e17, vinf 1",
    ]

    suptitle = "mdot_acc = 4e-2"
    outname = "cv_plot_star_models.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")
    
    #
    # simple atom banding test
    #
    
    files = [
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim12/tdeStar.spec",
        "/home/saultyevil/PySims/TDE/cv_star_model/simple_atom_models/sim12_agn_banding/tdeStar.spec",
    ]
    
    titles = [
        "cv banding",
        "agn banding",
    ]
    
    suptitle = "co.rad 8.86e12, Mwind = 0.25Macc, 1055R,1, macc = 4e-2"
    outname = "cv_plot_banding_test.pdf"
    plot_3_inclinations(files, titles, suptitle, inds, legend, 2, 3, outname, "cv")

    return


def main():
    """
    Main controlling function
    """

    plot_agn_grid_1_30()
    plot_agn_grid_1_20()
    plot_agn_grid_runs()  # new stuff :-)
    plot_cv_grid()
    star_models()

    return


if __name__ == "__main__":
    main()
