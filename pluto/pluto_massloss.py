#!/usr/bin/env python

import csv, sys, os, array, warnings,io
import numpy as np
from astropy import constants as c
from matplotlib import pyplot as plt
import pyPLUTO as pp
import pluto_v_hydro_sub as vhs
import subprocess






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


UNIT_MASS=(UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
UNIT_ACCELERATION=(UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
UNIT_FORCE=(UNIT_MASS*UNIT_ACCELERATION)
UNIT_TIME=(UNIT_LENGTH/UNIT_VELOCITY)
UNIT_PRESSURE=(UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)
KELVIN=UNIT_VELOCITY*UNIT_VELOCITY*c.m_p.cgs/c.k_B.cgs




if len(sys.argv)>1:
    istart=int(sys.argv[1])
else:
    istart=1

istart=np.max([1,istart])


if len(sys.argv)>2:
    itheta_disk=int(sys.argv[2])
else:
    itheta_disk=-1


ml=[]


cmdline="tail -1 dbl.out"
proc=subprocess.Popen(cmdline,shell=True,stdout=subprocess.PIPE) #This mess gets the last hdffile
ifile=int(proc.stdout.read().split()[0])

time = []
massloss = []
for ii in range(istart,ifile+1):
    sys.stdout = open(os.devnull, 'w')
    fname=ii
    try:
        D=pp.pload(ii)
        error=0
    except:
        error=1
    data=vhs.pluto_to_zeus(D)
    sys.stdout = sys.__stdout__

    density=data["Data"]["DENSITY"]["data"]
    pressure=data["Data"]["PRESSURE"]["data"]
    energy=data["Data"]["TOTAL ENERGY"]["data"]*UNIT_FORCE*UNIT_LENGTH

    temperature=pressure*KELVIN*0.6/(density)

    density=density*UNIT_DENSITY

    nd=density/(1.43*c.m_p.value)
    v_r=data["Data"]["1-VELOCITY"]["data"]*UNIT_VELOCITY
    v_t=data["Data"]["2-VELOCITY"]["data"]*UNIT_VELOCITY

    v=np.sqrt(v_r**2+v_t**2)
    theta=(data["Data"]["2-VELOCITY"]["x1"])
    r=(data["Data"]["2-VELOCITY"]["x2"])*UNIT_LENGTH

    nt=len(theta)
    nr=len(r)

    t_x=5.6e7


    R_ic=c.G.cgs.value*7.0*c.M_sun.cgs.value*0.6*c.m_p.cgs.value/(c.k_B.cgs.value*t_x/4.0)






    theta_ratio=(theta[2]-theta[1])/(theta[1]-theta[0])
    dtheta=[]
    dtheta.append((theta[1]-theta[0])/(0.5*(1.0+theta_ratio)))
    theta_min=theta[0]-0.5*dtheta[-1]
    theta_max=theta_min+dtheta[-1]*(1.0-theta_ratio**len(theta))/(1.0-theta_ratio)
    for i in range(len(theta)-1):
        dtheta.append(theta_ratio*dtheta[-1])

    r_ratio=(r[2]-r[1])/(r[1]-r[0])
    dr=[]
    dr.append((r[1]-r[0])/(0.5*(1.0+r_ratio)))
    rmin=r[0]-0.5*dr[-1]
    rmax=rmin+dr[-1]*(1.0-r_ratio**len(r))/(1.0-r_ratio)
    for i in range(len(r)-1):
        dr.append(dr[0]*r_ratio**float(i+1))

    outer_mass_loss=density[:,-1]*v_r[:,-1]
    cum_outer_mass_loss=[]
    cum_outer_mass_loss2=[]
    cum_outer_mass_loss.append(outer_mass_loss[0]*np.sin(theta[0])*(dtheta[0]))
    cum_outer_mass_loss2.append(outer_mass_loss[-1]*np.sin(theta[-1])*(dtheta[-1]))
    for i in range(len(outer_mass_loss)-1):
        cum_outer_mass_loss.append(cum_outer_mass_loss[-1]+(outer_mass_loss[i+1])*np.sin(theta[i+1])*(dtheta[i+1]))
    cum_outer_mass_loss=np.array(cum_outer_mass_loss)*rmax*rmax*4.0*np.pi #4pi means we are looking at the mass outflow from the whole sphere


    rho=density[itheta_disk]
    vz=(v_r[itheta_disk]*np.cos(theta[itheta_disk])-v_t[itheta_disk]*np.sin(theta[itheta_disk]))
    disk_mass_loss=(rho*vz)
    disk_mass_loss_r=disk_mass_loss*r
    cum_disk_mass_loss=[]
    cum_disk_mass_loss2=[]
    nr=len(disk_mass_loss)
    cum_disk_mass_loss.append(disk_mass_loss[0]*2.0*np.pi*(r[0]*dr[0]))
    cum_disk_mass_loss2.append(disk_mass_loss[-1]*2.0*np.pi*(r[-1]*dr[-1]))
    for j in range(len(disk_mass_loss)-1):
        cum_disk_mass_loss.append(cum_disk_mass_loss[-1]+disk_mass_loss[j+1]*2.0*np.pi*(r[j+1]*dr[j+1]))
        cum_disk_mass_loss2.append(cum_disk_mass_loss2[-1]+disk_mass_loss[nr-2-j]*2.0*np.pi*(r[nr-2-j]*dr[nr-2-j]))



    cum_disk_mass_loss=2.0*np.array(cum_disk_mass_loss) #Both sides of the disk
    cum_disk_mass_loss2=2.0*np.array(cum_disk_mass_loss2) #Both sides of the disk



    if error==0:
        # print fname,D.SimTime*UNIT_TIME,cum_outer_mass_loss[itheta_disk]
        time.append(D.SimTime*UNIT_TIME)
        # aa = 7.02066390388e-09 # * UNIT_MASS # / UNIT_TIME
        massloss.append(cum_outer_mass_loss[itheta_disk])

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.plot(time, massloss/4e17)
ax.set_ylabel("Mass loss rate", fontsize=14)
ax.set_xlabel("Time (s)", fontsize=14)
plt.savefig("masslss_time.png")
plt.show()



#	fig=plt.figure(figsize=(12,8))
#	ax=fig.add_subplot(111)

#	levels=np.linspace(0.5,2.0,100)
#	rho=data["Data"]["DENSITY"]["data"][-2]
#	r=data["Data"]["DENSITY"]["x2"]

#	ax.semilogy(r,rho)
#	plt.axvline(x=9.2321292e+10)
#	plt.ylim([1e-14,1e-11])
#	fig.savefig("dens_"+str(j).zfill(3)+".png")
#	j=j+1
#	plt.close(fig)
