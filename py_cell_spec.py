#!/usr/bin/env python

import matplotlib.pyplot as plt

import csv, sys, os, array, warnings
import numpy as np

from astropy import constants as consts
from astropy import units as u


inp=open(sys.argv[1],'r')


numin=[]
bandmin=[]
numax=[]
bandmax=[]
model=[]
pl_log_w=[]
pl_alpha=[]
exp_w=[]
exp_temp=[]



for line in inp.readlines():
	data=line.split()
	numin.append(float(data[1]))
	bandmin.append(float(data[3][:-1]))
	numax.append(float(data[5]))
	bandmax.append(float(data[7][:-1]))
	model.append(int(data[9]))
	pl_log_w.append(float(data[11]))
	pl_alpha.append(float(data[13]))
	exp_w.append(float(data[15]))
	exp_temp.append(float(data[17]))
	

freq=[]
f_nu=[]

for i in range(len(numin)):
	if numax[i]>numin[i]:
		freq_temp=np.logspace(np.log10(numin[i]),np.log10(numax[i]),101)
		for nu in freq_temp:
			freq.append(nu)
			if model[i]==1:
				f_nu.append(10**(pl_log_w[i]+np.log10(nu)*pl_alpha[i]))
			elif model[i]==2:
				f_nu.append(exp_w[i]*np.exp((-1.0*consts.h.cgs.value*nu)/(exp_temp[i]*consts.k_B.cgs.value)))
			else:
				f_nu.append(0.0)
	else:
		freq.append(bandmin[i])
		freq.append(bandmax[i])
		f_nu.append(0.0)
		f_nu.append(0.0)
					
				
xi=0.0				
for i in range(len(freq)-1):				
	if freq[i]>3.288e15:
		xi=xi+((f_nu[i+1]+f_nu[i])/2.0)*(freq[i+1]-freq[i])	
		
			
			
fig=plt.figure(figsize=(12,8))
ax=fig.add_subplot(111)
ax.loglog(freq,f_nu)
#ax.set_ylim([1e-30,1e-3])
#plt.xlim([1e14,1e19])
#plt.ylim([1e-14,1e-6])

ax.set_title(sys.argv[1])

ax.set_xlabel(r"$\rm{Frequency(Hz)}$")
ax.set_ylabel(r"$\rm{J_{\nu}~in~cell(ergs~s^{-1}~cm^{-3}~Sr^-1~Hz^-1)}$")
fig.savefig(sys.argv[1]+'.png')
plt.close(fig)





