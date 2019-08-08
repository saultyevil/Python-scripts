#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import tde_util
import numpy as np
import pandas as pd
import py_plot_util
import tde_spec_plot
from consts import *
from matplotlib import pyplot as plt


SMOOTH = 15
WMIN = 100
WMAX = 6000
TDE_OBJ = "iPTF15af"
DEFAULT_DIST = 100 * PARSEC
FTYPE = "png"
VERBOSE = False


grid_angles = ["20", "40", "60", "62", "70", "75", "80", "85", "89"]
COMP_MODEL_SPEC_PATH = "zzz_orig/tde.spec"

try:
    comp_model_spec = py_plot_util.read_spec_file(COMP_MODEL_SPEC_PATH, " ", pandas_table=True)
    comp_model_angles = py_plot_util.spec_inclinations_pandas(comp_model_spec)
except IOError:
    print("Can't open base agn model. Is something wrong?")
    comp_model_spec = np.zeros((10,2))
    comp_model_spec = grid_angles.copy()


def plot_comparison(grid_name, root, i_to_plot, subplots, sim_dirs, tde_obj="iPTF15af", wmin=WMIN, wmax=WMAX,
                    smooth=SMOOTH):
    """
    Create a comparison plot of a bunch of spec files.
    
    Parameter
    ---------
    grid_name   The grid parameter 
    root        The common root name of the simulations to plot
    subplots    The shape of the subplot panels (tuple)
    sim_dirs    A list of directories for the grid simulation
    
    Returns
    -------
    None
    """

    nrows, ncols = subplots
    fig, ax = plt.subplots(nrows, ncols, figsize=(32, 16), squeeze=False)
    print("Plotting {} grid,\n\t{},\nfor angles,\n\t{}".format(grid_name, sim_dirs, i_to_plot))
    
    # hack to make i_to_plot a list
    if type(i_to_plot) == float or type(i_to_plot) == int or type(i_to_plot) == str:
        i_to_plot = [i_to_plot]
    
    tde_spec_plot.TDE_OBJ = tde_obj  # hacky fix
    tde, dist, reference = tde_spec_plot.get_tde_spectrum()
    line_ids = py_plot_util.common_lines()

    index = 0 
    for i in range(nrows):
        for j in range(ncols):
            
            if index > len(i_to_plot) - 1:
                break
            
            # Plot TDE spectrum for comparison
            ax[i, j].semilogy(tde[:, 0], tde[:, 1], "b", label=TDE_OBJ)
            
            # Plot AGN base model - if the inclination angle is in the spec file
            if i_to_plot[index] in comp_model_angles:
                base_lambda = comp_model_spec["Lambda"].values.astype(float)
                base_flux = comp_model_spec[i_to_plot[index]].values.astype(float)
                base_flux = py_plot_util.smooth_1d_array(base_flux, 30)
                base_flux *= DEFAULT ** 2 / dist ** 2
                ax[i, j].semilogy(base_lambda, base_flux, "k", label="zzz_orig")
                
            ymin = 1e-14
            ymax = 1e-13
            for d in sim_dirs:
                spec_filename = d + "/" + root + ".spec"
                try:
                    spectrum = py_plot_util.read_spec_file(spec_filename, pandas_table=True)
                except IOError:
                    print("Can't open spec {}, check pls".format(spec_filename))
                inclinations = py_plot_util.spec_inclinations_pandas(spectrum)
                if i_to_plot[index] not in inclinations:
                    continue
                wavelength = spectrum["Lambda"].values.astype(float)
                flux = spectrum[i_to_plot[index]].values.astype(float)
                flux = py_plot_util.smooth_1d_array(flux, smooth)
                flux *= DEFAULT_DIST ** 2 / dist ** 2
                ax[i, j].semilogy(wavelength, flux, label=spec_filename)
                
                tmp_max, tmp_min = py_plot_util.define_ylims(wavelength, flux, wmin, wmax)
                if tmp_max > ymax:
                    ymax = tmp_max
                if tmp_min < ymin:
                    ymin = tmp_min
                
            ax[i, j].legend(loc="upper right")
            ax[i, j].set_ylim(ymin, ymax)
            ax[i, j].set_xlim(wmin, wmax)
            ax[i, j].set_title("i = " + i_to_plot[index])

            py_plot_util.plot_line_ids(ax[i, j], line_ids)
            
            # increment index counter
            index += 1
            
    fig.suptitle(grid_name)
            
    grid_name = grid_name.replace("/", "_")
        
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    plt.savefig(grid_name+".pdf")
    plt.show()
            
    return
    

if __name__ == "__main__":
    print("Don't run this :-)")
