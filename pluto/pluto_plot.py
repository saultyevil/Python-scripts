#!/usr/bin/env python

"""
Provided by Nick Higginbottom
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import pluto_v_hydro_sub as vhs
import pluto_lin_plot_cont as lpc
import pickle
from astropy import constants as ac
import pyPLUTO as pp

from time import time
from scipy.optimize import brentq

# Set up some basic parameters of the system

Mbh = 7.0 * ac.M_sun.cgs.value
T_x = 4.0 * 1.4e7
mu = 0.6
gamma = 5. / 3.

# Compute the compton radius

Ric = (ac.G.cgs.value * Mbh * ac.m_p.cgs.value * mu / (ac.k_B.cgs.value * (T_x / 4.0)))

# First get scaling factors from the definitions file

inp = open('definitions.h', 'ro')
for line in inp.readlines():
    data = line.split()
    if len(data) > 1:
        if data[1] == 'UNIT_DENSITY':
            UNIT_DENSITY = float(data[2])
        elif data[1] == 'UNIT_LENGTH':
            UNIT_LENGTH = float(data[2])
        elif data[1] == 'UNIT_VELOCITY':
            UNIT_VELOCITY = float(data[2])
inp.close()

# Compute deived scaling factors

UNIT_MASS = (UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH)
UNIT_ACCELERATION = (UNIT_VELOCITY * UNIT_VELOCITY / UNIT_LENGTH)
UNIT_FORCE = (UNIT_MASS * UNIT_ACCELERATION)
UNIT_TIME = (UNIT_LENGTH / UNIT_VELOCITY)
UNIT_PRESSURE = (UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY)

# Compute the number that transforms from pressure to temperature

KELVIN = mu * UNIT_VELOCITY * UNIT_VELOCITY * ac.m_p.cgs / ac.k_B.cgs

# open the actual data file

try:
    fname = int(sys.argv[1])
except:
    print("No filename given")
    sys.exit(1)

# set a file name for the streamline data - this takes ages to compute, so we want to store it

streamfile = str(fname) + 'stream1.pk1'

# See if the user input a theta bin to start the streamlines from - if not default to -5

if len(sys.argv) > 2:
    istart = int(sys.argv[2])
else:
    istart = -5  # Changed originally from itheta_disk

# Load the data -- uses pyPLUTO (pp) see http://plutocode.ph.unito.it/

D = pp.pload(fname)

# Now convert to zeus format so the old scripts work

data = vhs.pluto_to_zeus(D)

# Get density and pressure in code units

density = data["Data"]["DENSITY"]["data"]
pressure = data["Data"]["PRESSURE"]["data"]

###Compute the speed of sound everywhere in the grid.

g_gamma = (5. / 3.)
c = np.sqrt(g_gamma * pressure / density) * UNIT_VELOCITY

# Load or compute temperatures

try:
    temperature = np.transpose(D.T)
except:
    print("No temperature data - computing")
    temperature = pressure * KELVIN / density

nd = density / (1.43 * ac.m_p.value)

# Convert velocities to cgs - needed for streamline calculations

data["Data"]["1-VELOCITY"]["data"] = data["Data"]["1-VELOCITY"]["data"] * UNIT_VELOCITY
data["Data"]["2-VELOCITY"]["data"] = data["Data"]["2-VELOCITY"]["data"] * UNIT_VELOCITY
data["Data"]["3-VELOCITY"]["data"] = data["Data"]["3-VELOCITY"]["data"] * UNIT_VELOCITY

# And load into handy locally named variables

v_r = data["Data"]["1-VELOCITY"]["data"]
v_t = data["Data"]["2-VELOCITY"]["data"]
v_p = data["Data"]["3-VELOCITY"]["data"]

# Obtain geometric variables

r = data["Data"]["1-VELOCITY"]["x2"] * UNIT_LENGTH
theta = data["Data"]["1-VELOCITY"]["x1"]
rmax = max(r)
rmin = min(r)

# This is needed for the scripts to work

data["Coord_sys"] = 'spol'

# Get the time

try:
    time = D.SimTime * UNIT_TIME
except:
    time = "no time info"

# See if the user input a maximum radius to start streamlines from - if not default

if len(sys.argv) > 3:
    radmax = float(sys.argv[3])
else:
    radmax = 0.95 * rmax

# Test for sanity

if radmax < rmin:
    print("User input radmax is less than rmin, resetting")
    radmax = rmax

# Compute a set of streamline roots

wroot = np.logspace(np.log10(1.1 * rmin), np.log10(radmax), 15)

# Before using this script, we need to ensure the radial infor for this dictionary is in cgs

data["Data"]["1-VELOCITY"]["x2"] = data["Data"]["1-VELOCITY"]["x2"] * UNIT_LENGTH

# Try to compute streamlines - it will only execute if there is no currently computed streamline file

try:
    savefile = open(streamfile, 'rb')
except:
    vhs.multi_stream(data, wroot=wroot, picklefile=streamfile, area_flag="n", tstart=istart, npoints=10000)
    savefile = open(streamfile, 'rb')

# Open the streamline file - it is either existing, or will have been saved above

stream_data = pickle.load(savefile)
streamw = stream_data["sw"]
streamz = stream_data["sz"]
savefile.close()

# We will be sending a sub-dictionary to the plotter - ensure the values for thias dictionary are in cgs

data["Data"]["DENSITY"]["data"] = data["Data"]["DENSITY"]["data"] * UNIT_DENSITY
data["Data"]["DENSITY"]["x2"] = data["Data"]["DENSITY"]["x2"] * UNIT_LENGTH

# Set up the dictionary to send to the plotter

data1 = data["Data"]["DENSITY"]
# data1["fig_size"]=[3,3] #Figure size in inches - keep it square!!
data1["long_name"] = ""  # Figure name
data1["linlog"] = "log"  # log to plot quantity in log space, lin for linear
data1["dist_scale"] = 1.  # Scaling factor for distances, e.g. cab use Ric to plot in compton radii
data1["xyname"] = "/cm"  # Name for x and z axes - ideally used to show units, either cm or /Ric
data1["rmax"] = rmax  # Maximum radius to plot - doesnt work that well!
data1["cmap"] = 'inferno'  # Color map for contours
data1["fill_cont"] = "fill"
data1["contour_label"] = r"Density $\rm{(g~cm^{-3})}$"  # Label for contour scale
data1["data_scale"] = np.linspace(-18, -12, 701)  # levels for countors
data1["data_ticks"] = np.linspace(-18, -12, 7)  # Levels for contour scale
fig = lpc.polar_contour(data1)  # Call the script - resulting figure is in fig

# We now want to plot contours of mach number of the top

# First we obtain the poloidal velocity

vpol = np.sqrt(v_r ** 2 + v_t ** 2)

# And we divide this by the sound speed to get M - mach number everywhere

M = vpol / c

# Now we overplot the Mach contours


CS = fig.axes[2].contour(np.pi / 2.0 - theta, r, np.transpose(M), np.linspace(1, 5, 5), colors='k', hold='on')

# And now we overplot the streamlines

for i in range(len(streamw)):
    fig.axes[3].plot(streamw[i], streamz[i], '0.5')

# We now want to compute the velocity vectors, to plot with arrows, we need to go from r,theta to w,z

# Make some empty arrays

v_w = np.zeros(np.shape(v_r))
v_z = np.zeros(np.shape(v_r))
mach_number = []

# Populate the arrays by simply computing v_w and v_z from v_r and v_t

for i in range(len(r)):
    for j in range(len(theta)):
        v_z[j][i] = (v_r[j][i] * np.cos(theta[j]) - v_t[j][i] * np.sin(theta[j]))
        v_w[j][i] = (v_r[j][i] * np.sin(theta[j]) + v_t[j][i] * np.cos(theta[j]))

# We now interpolate on our r,theta grid to get a recatnagular grid of arrows

# Make our choices of grid - change with caution!

vw_grid = np.linspace(rmax * np.sin(theta[0]), rmax * 0.95, 20)
vz_grid = np.linspace(rmax * np.cos(theta[-2]), rmax * 0.95, 15)

# Make arrays to fit the grid
vw_vals = np.zeros([len(vw_grid), len(vz_grid)])
vz_vals = np.zeros([len(vw_grid), len(vz_grid)])

# Interpolate on our r-theta grid to get the rectuangular grid

for i in range(len(vw_grid)):
    for j in range(len(vz_grid)):
        r_test = np.sqrt(vw_grid[i] ** 2 + vz_grid[j] ** 2)
        theta_test = np.arctan(vw_grid[i] / vz_grid[j])
        if r_test < max(r) and r_test > min(r) and theta_test > min(theta) and theta_test < max(theta):
            vw_vals[i][j] = (vhs.interpolate(theta, r, theta_test, r_test, v_w))
            vz_vals[i][j] = (vhs.interpolate(theta, r, theta_test, r_test, v_z))
        else:
            vw_vals[i][j] = 0.0
            vz_vals[i][j] = 0.0

# Work out the maximum velocity in the grid to scale the arrows

# vel=np.sqrt(vw_vals*vw_vals+vz_vals*vz_vals)
# vmax=10**(float(int(np.log10(np.max(vel)))+1))

vmax = 1e8  # Or just set it by hand

vw_vals[-5][-2] = vmax  # This sets one point to the maximum velocity, to make a label

# Plot the arrows

fig.axes[3].quiver(vw_grid, vz_grid, np.transpose(vw_vals), np.transpose(vz_vals), units='width', scale=vmax,
                   scale_units='inches')

# Plot the label - rather complicated this

if vmax < (100. * 1000.):
    fig.axes[3].text(0.8 * rmax, (1.72 / 2.0) * rmax, r"$\rm{" + str(vmax) + "~cm~s^{-1}}$")
else:
    fig.axes[3].text(0.8 * rmax, (1.72 / 2.0) * rmax, r"$\rm{" + str(vmax / 100. / 1000.) + "~km~s^{-1}}$")

fig.axes[3].text((1.6 / 2.0) * rmax, rmax, r"Time=" + str(np.float(time)) + " s")
title_name = str(fname)

# Add a title

fig.axes[3].text((1 / 2.0) * rmax, 1.05 * rmax, title_name)

# Save the figure

plt.savefig(str(fname) + '_dens_stream.png')

# Close the figure

plt.close(fig)

# We now plot a temperature plot


data1["data"] = temperature
# data1["fig_size"]=[3,3]
data1["long_name"] = ""
data1["linlog"] = "log"
data1["dist_scale"] = 1.0
data1["xyname"] = "/cm"
data1["rmax"] = rmax
data1["fill_cont"] = "fill"
data1["contour_label"] = r"Temperature $\rm{(K)}$"
# data1["data_scale"]=np.linspace(4.5,7.5,501)
data1["data_scale"] = [2, 3, 4, 5, 6, 7, 8]
# data1["data_ticks"]=np.linspace(4.5,7.5,5)
data1["data_ticks"] = [2, 3, 4, 5, 6, 7, 8]
fig = lpc.polar_contour(data1)

plt.savefig(str(fname) + '_temperature.png')
plt.close(fig)
