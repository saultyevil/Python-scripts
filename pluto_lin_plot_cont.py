#!/usr/bin/env python -i

"""
Provided by Nick Higginbottom
"""

import csv, sys, os, array, warnings, subprocess
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from matplotlib import rc
import matplotlib.colors as col
import matplotlib.cm as cm
# Force matplotlib to not use any Xwindows backend.
mpl.use('Agg')
import matplotlib.pyplot as plt


#polar_contour plots an input dictionary as a polar plot, with rectilinear axes
#data_inp should be a dictionary with the data in "data", r in "x2" and theta in "x1".
#"long_name" contains the text to be used as the figure title, and if it is something the
#code understands, it will automatically put a legend on the contour. 
#If the dictionary contains a member "linlog" then this controls whether we want to 
#plot the log of the quantity, otherwise it defaults to log (normally sensible)
#if the disctionary contains a member "dist_scale" then this is used to scale the r and z
#axes, otherwise it is set to 1

def polar_contour(data_inp):

#The data as read in from a zeus hdf file is not the same way round as python expects!

	try:
		data=np.transpose(data_inp["data"])
	except:
		print("There is no data to plot, make sure the dictionary has 'data' in it")


	try:
		if data_inp["linlog"]=="lin" or data_inp["linlog"]=="log":
			linlog=data_inp["linlog"]
		else:
			print("I dont understand the switch ",data_inp["linlog"]," should be 'lin' or 'log'")
			print("Defaulting to log")
			linlog="log"
	except:
		linlog="log"

#Obtain the length scale from the dictionary, or default to 1

	try:
		dist_scale=data_inp["dist_scale"]
	except:
		dist_scale=1.0

#Obtain the length scale from the dictionary, or default to 1

	try:
		fig_size=data_inp["fig_size"]
	except:
		fig_size=(7,7.5)

	try:
		font_size=data_inp["font_size"]
	except:
		font_size=11.2


#Obtain the scale for the data from the dictionary, or do my best

	try:
		scale=data_inp["data_scale"]
	except:
		if linlog=="log":
			
#			imin=int(np.log10(np.percentile(data,5.)))			
#			imax=int(np.log10(np.percentile(data,95.)))
			try:
				imin=int(np.log10(1.1*np.min(data)))
			except:
				print("minimum of data=",np.min(data))
				imin=0.0
			imax=int(np.log10(0.9*np.max(data)))+1
			scale=np.linspace(imin,imax,(imax-imin)*100+1)
		elif linlog=="lin":
			min_sign=abs(np.percentile(data,5.))/np.percentile(data,5.)
			max_sign=abs(np.percentile(data,95.))/np.percentile(data,95.)
			imin=min_sign*10.**float(int(np.log10(min_sign*np.percentile(data,5.))))
			imax=max_sign*10.**float(int(np.log10(max_sign*np.percentile(data,95.))))
			scale=np.linspace(imin,imax,101)	

#Obtain the tick range for the data from the dictionary, or try my best

	try: 
		ticks=data_inp["data_ticks"]
	except: 

		if linlog=="log":
#			imin=int(np.log10(np.percentile(data,5.)))			
#			imax=int(np.log10(np.percentile(data,95.)))
			try:
				imin=int(np.log10(np.min(data)))
			except:
				print("minimum of data=",np.min(data))
				imin=0.0
			imax=int(np.log10(np.max(data)))+1
			ticks=np.linspace(imin,imax,(imax-imin)+1)
		elif linlog=="lin":
			min_sign=abs(np.percentile(data,5.))/np.percentile(data,5.)
			max_sign=abs(np.percentile(data,95.))/np.percentile(data,95.)
			imin=min_sign*10.**float(int(np.log10(min_sign*np.percentile(data,5.))))
			imax=max_sign*10.**float(int(np.log10(max_sign*np.percentile(data,95.))))
			ticks=np.linspace(imin,imax,11)


#Extract the r and theta corrdinates, I hope they were set correctly....

	try:	
		r=data_inp["x2"]
		theta=data_inp["x1"]
	except:
		print("One or other of the coordinates is empty, make sure the dict has 'x1' and 'x2'")


#Try to work out if we want a filled contour plot, or lines, or both - the default is both

	try:
		fill_cont=data_inp["fill_cont"]
	except:
		fill_cont="both"
#		print fill_cont

	try:
		cticks=data_inp["contour_ticks"]
	except:
		cticks=ticks


	try:
		rmax=data_inp["rmax"]
	except:
		rmax=max(r)

		
	try:
		xylabel=data_inp["xyname"]
	except:
		xylabel="(cm)"

	try:
		c_label=data_inp["contour_label"]
	except:
		if  data_inp["long_name"][0:7]=="DENSITY" or data_inp["long_name"][0:7]=="Density":
			c_label=r"$\rm{\log~\rho~(g~cm^{-3})}$"
		elif data_inp["long_name"][0:6]=="Number":
			c_label=r"$\rm{\log~n_{H}~(cm^{-3})}$"
		elif data_inp["long_name"][0:4]=="TEMP" or data_inp["long_name"][0:4]=="Temp":
			c_label=r"$\rm{\log~T~(K)}$"
		elif data_inp["long_name"][0:4]=="MACH" or data_inp["long_name"][0:4]=="Mach":
			c_label=r"Dimensionless Mach number"
		else:
			c_label=r"NO LABEL SUPPLIED - SET WITH contour label"
		
	try:
		cmap=data_inp["cmap"]
	except:
		cmap="plasma"



	theta=(np.pi/2)-theta
	
	
	rc('font',**{'family':'serif','serif':['Times']})
	rc('font',size=font_size)
	rc('text', usetex=True)



	#Attempt to set a sensible scale  

	
	fig=plt.figure(figsize=fig_size)
	
	#These lines put the colour bar in a good location
	if fill_cont=="fill" or fill_cont=="both":
		ax0=fig.add_axes([0.1,0.0,0.8,0.7],projection='polar',frameon=False)
		
		ax0.get_xaxis().set_visible(False)		
		ax0.get_yaxis().set_visible(False)
		cax=ax0.contourf([1,1],[1,1],np.zeros([2,2]),scale,extend='both',cmap=cmap)
		ax0.set_ylim(0,1e10)
		cbar=plt.colorbar(cax,orientation='horizontal') 
		cbar.set_ticks(ticks)
		cbar.set_ticklabels(ticks)
		cbar.ax.set_xlabel(c_label)
		
	
	#These lines plot the actual data
	
	ax1=fig.add_axes([-0.7,-0.55,1.6,1.493],projection='polar',frameon=False)
#	ax1=fig.add_axes([-0.7,-0.55,2.0,2.0],projection='polar',frameon=False)

	
	if fill_cont=="fill" or fill_cont=="both":
		if linlog=="log":
			ax1.contourf(theta, r/dist_scale, np.log10(data),scale,extend='both',cmap=cmap)
		elif linlog=="lin":
			ax1.contourf(theta, r/dist_scale, data,scale,extend='both',cmap=cmap)
		else:
			print("I dont understand the switch ",linlog," should be 'lin' or 'log'")
#		ax1.set_rmax(rmax)
	if fill_cont=="cont" or fill_cont=="both":
		if linlog=="log":
			print(ticks)
			CS=ax1.contour(theta, r/dist_scale, np.log10(data),cticks,colors='k',hold='on')
		elif linlog=="lin":
			CS=ax1.contour(theta, r/dist_scale, data,cticks,colors='k',hold='on')
		else:
			print("I dont understand the switch ",linlog," should be 'lin' or 'log'")
		for c in CS.collections:
    			c.set_linestyle('solid')
#		plt.clabel(CS, inline=1, fontsize=10)



	ax1.set_rlim((0,rmax/dist_scale))
	ax1.get_xaxis().set_visible(False)
	ax1.get_yaxis().set_visible(False)
	
	#Finally, these lines apply a recangular set of axes over the top, to allow scales
	
	ax2=fig.add_axes([0.1,0.196,0.8,0.746],xlim=[0,rmax/dist_scale],ylim=[0,rmax/dist_scale])
	ax2.set_xlabel(r"$\rm{\omega "+xylabel+"}$")
	ax2.set_ylabel(r"$\rm{z"+xylabel+"}$")

	ax2.set_title(data_inp["long_name"])
	ax2.patch.set_facecolor('none')
#We return the figure, so the user can modify it or write it out as required

	return(fig)

