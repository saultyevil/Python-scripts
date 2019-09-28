#!/usr/bin/env python

"""
Provided by Nick Higginbottom. Thanks!
"""


import sys
import numpy as np
from astropy import constants as c
import pyPLUTO as pp

GAMMA= 5.0 / 3.0


#First get scaling factorz from the definitions file

inp=open('definitions.h','ro')
for line in inp.readlines():
	data=line.split()
	if len(data)>1:
		if data[1]=='UNIT_DENSITY':
			UNIT_DENSITY=float(data[2])
		elif data[1]=='UNIT_LENGTH':
			UNIT_LENGTH=float(data[2])
		elif data[1]=='UNIT_VELOCITY':
			UNIT_VELOCITY=float(data[2])

#Compute deived scaling factors

UNIT_MASS=(UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
UNIT_ACCELERATION=(UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
UNIT_FORCE=(UNIT_MASS*UNIT_ACCELERATION)
UNIT_TIME=(UNIT_LENGTH/UNIT_VELOCITY)
UNIT_PRESSURE=(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)

#Compute the numeber that transforms from pressure to temperature

KELVIN=UNIT_VELOCITY*UNIT_VELOCITY*c.m_p.cgs/c.k_B.cgs


inp.close()


#open the actual data file

fname=int(sys.argv[1])
D=pp.pload(fname)


#Get the pressure and density

density=np.transpose(D.rho)*UNIT_DENSITY
pressure=np.transpose(D.prs)*UNIT_PRESSURE

#Compute the internal energy from the pressure

energy=pressure/(GAMMA - 1.)

#Compute/get number densities

#nd=density/(1.43*c.m_p.value)
try:
	ne=np.transpose(D.ne)
	nh=np.transpose(D.nh)
except:
	print("No ne or nh fields, using 1.43 as scaling to nh")
	nh=density/(1.43*c.m_p.cgs).value
	ne=nh*1.21

#Get the velocities

v_r=np.transpose(D.vx1)*UNIT_VELOCITY
v_t=np.transpose(D.vx2)*UNIT_VELOCITY
v_p=np.transpose(D.vx3)*UNIT_VELOCITY

#And compute the speed

v=np.sqrt(v_r**2+v_t**2)

#Get the cooling rates (if here)

try:
	line_c=np.transpose(D.line_c)
	xray_h=np.transpose(D.xray_h)
	comp_c=np.transpose(D.comp_c)
	comp_h=np.transpose(D.comp_h)
	brem_c=np.transpose(D.brem_c)
except:
	print("No cooling rate info")

try:
	line_c_pre=np.transpose(D.line_c_pre)
	xray_h_pre=np.transpose(D.xray_h_pre)
	comp_c_pre=np.transpose(D.comp_c_pre)
	comp_h_pre=np.transpose(D.comp_h_pre)
	brem_c_pre=np.transpose(D.brem_c_pre)
except:
	print("No cooling rate prefactors")


#Get optcially thin ionization parameter if here

try:
	xi=np.transpose(D.XI)
except:
	print("No ionization parameter")
	
#Get temperature - if present or calculate it
	
try:
	temperature=np.transpose(D.T)
except:
	print("No temperature data - computing")
	temperature=pressure/UNIT_PRESSURE*KELVIN*0.6/(density/UNIT_DENSITY)

#Get the geometric quantities

r=D.x1*UNIT_LENGTH
theta=D.x2
