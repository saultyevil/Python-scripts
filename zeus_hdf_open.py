#!/usr/bin/env python

"""
Most of the code was provided by the wonderful Nick Higginbottom, who can be found on GitHub here,
	https://github.com/Higginbottom/
"""


import sys
import numpy as np
import pluto_hdf_utils
from astropy import constants as c


def write_z_v_density(fname, z, density, Ric):
	with open(fname, "w") as f:
		f.write("# R = {:2.4f} Ric\n".format(Ric))
		f.write("# z          rho\n")
		for i in range(len(z)):
			f.write("{:^+7.5e} {:^+7.5e}\n".format(z[i], density[i]))
	return


try:
	fname=sys.argv[1]
except:
	print("No filename provided")
	sys.exit(1)

# Read in data using pyhdf
data=pluto_hdf_utils.get_hdf_data(fname)

if len(sys.argv) > 3:
	itheta_disk=int(sys.argv[3])
else:
	itheta_disk = -1
if len(sys.argv) > 4:
	iradmax = int(sys.argv[4])
else:
	iradmax = -1

try:
	print("HDF file produced at t=",data["Time"])
except:
	print("No time information")

# Add the data into their own arrays
density=data["Data"]["DENSITY"]["data"]
nd=density/(1.43*c.m_p.value)
energy=data["Data"]["TOTAL ENERGY"]["data"]
v_r=data["Data"]["1-VELOCITY"]["data"]
v_t=data["Data"]["2-VELOCITY"]["data"]
v_p=data["Data"]["3-VELOCITY"]["data"]
line_c=data["Data"]["LINE_C"]["data"]
xray_h=data["Data"]["XRAY_H"]["data"]
comp_c=data["Data"]["COMP_C"]["data"]
comp_h=data["Data"]["COMP_H"]["data"]
brem_c=data["Data"]["BREM_C"]["data"]
adiab_hc=data["Data"]["ADIABATIC"]["data"]
line_c_pre=data["Data"]["LINE_PRE"]["data"]
xray_h_pre=data["Data"]["XRAY_PRE"]["data"]
comp_c_pre=data["Data"]["COMP_C_PRE"]["data"]
comp_h_pre=data["Data"]["COMP_H_PRE"]["data"]
brem_c_pre=data["Data"]["BREM_PRE"]["data"]
theta=(data["Data"]["2-VELOCITY"]["x1"])
r=(data["Data"]["2-VELOCITY"]["x2"])
zeus_density=data["Data"]["DENSITY"]["data"]

# Various constants and numbers from the data
t_x=5.6e7
R_ic=c.G.cgs.value*7.0*c.M_sun.cgs.value*0.6*c.m_p.cgs.value/(c.k_B.cgs.value*t_x/4.0)
zeus_nh=zeus_density/(1.43*c.m_p.cgs).value
zeus_ne=zeus_nh*1.21	
temp=(2.0/3.0)*energy
temperature=temp/((density/(c.m_p.cgs.value*0.6))*c.k_B.cgs.value)
pressure=energy*(5./3.-1.)
xi=data["Data"]["XI"]["data"]
nt=len(theta)
nr=len(r)

try:
	rindex = int(sys.argv[2])
except:
	print("No r index given")
	sys.exit(1)

print("------------------")
print("Ric = {:e}".format(R_ic))
print("index = {}".format(rindex))
print("r[index] = {:e}".format(r[rindex]))
print("r/Ric = {}".format(r[rindex] / R_ic))
print("n theta = {}".format(nt))
print("n r = {}".format(nr))
print("density.shape = {}".format(density.shape))
print("temperature.shape = {}".format(temperature.shape))
print("------------------")

zidx = 0
nlim = 37
z = np.zeros(nlim)
rho = np.zeros(nlim)
print("\nR = {:2.4} Ric".format(r[rindex] / R_ic))
print("# z          rho\n")
for i in range(nt-nlim, nt):
	ztemp = z[zidx] = r[rindex] * np.cos(theta[i])
	dtemp = rho[zidx] = density[i][rindex]
	T = temperature[i][rindex]
	zidx += 1
	print("{:^+7.5e} {:^+7.5e} {:^+7.5e}".format(ztemp, dtemp, T))
write_z_v_density("grid.out", z,  rho, r[rindex] / R_ic)
