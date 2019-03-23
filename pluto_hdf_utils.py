#!/usr/bin/env python

"""
This code was provided by the wonderful Nick Higginbottom, who can be found on GitHub here,
	https://github.com/Higginbottom/
"""


from pyhdf import SD
import numpy as np


'''The following routine reads in hdf data from the file fname
   It takes a filename as input argument, and returns a dictionary
   The upper level of the dictionay contains
   Filename - the filename originally supplied
   Coord_sys - the coordinate system - NB, it only knows about spol at the moment
   Time - the time stamp of the simulation
   N_Data - the number of different data types
   Dims - the number of dimensions
   Data - This is in turn a dictionary that contains the different data
	This dictionary contains infividual dictionaries, each of which has
		data - an array of the data
		x1 - the first corrdinate (usually theta)
		x2 - the second coordinate (usually r)'''


def get_hdf_data(fname):


	hdf = SD.SD(fname)
	info= hdf.datasets()

	#Lets see what is inside

	hdf.info()
	data_sets=[]

	for name in sorted(info.keys()):
		if name[0:4]=="Data":
			sds=hdf.select(name)
			long_name=sds.attributes()["long_name"]
			for i in range(len(sds.attributes()["long_name"])):
				if long_name[i:i+2]=="AT":
					junk=i-1
			short_name=long_name[:junk]
			data_sets.append([name,long_name,short_name])

#Get the time from the last long name 
	if long_name[junk+9:] != "********":
		time=float(long_name[junk+9:])
	else:
		time=0.0

#Get the dimensions from the last data set

	dims=len((sds.info()[2]))

#Get the coordinate system from the last data set

	coord_sys=sds.attributes()["coordsys"]
	if coord_sys=='spherical polar':
		coord_sys='spol'
	else:
		print("get_hdf_data: I don't understand coordinate system ",coord_sys)
		exit()

#Now we know which of the datasets contain real data, we can extract all the data

	
	alldat={}
	alldat["Filename"]=fname	
	alldat["Coord_sys"]=coord_sys
	alldat["Time"]=time
	alldat["Data"]={}
	alldat["Data_names"]=np.array(data_sets)[:,2]
	alldat["N_data"]=len(data_sets)
	alldat["Dims"]=dims

#Loop over all the data sets in the hdf file - name each of the resulting dictionaries with the short name

	for i in range (len(data_sets)):
		sds=hdf.select(data_sets[i][0])
		data = sds.get()
		c1=info[data_sets[i][0]][0][0]
		c2=info[data_sets[i][0]][0][1]
		sds=hdf.select(c1)
		x1=sds.get()
		sds=hdf.select(c2)
		x2=sds.get()
		alldat["Data"][data_sets[i][2]]={}
		alldat["Data"][data_sets[i][2]]["long_name"]=data_sets[i][1]
		alldat["Data"][data_sets[i][2]]["data"]=data
		alldat["Data"][data_sets[i][2]]["x1"]=x1
		alldat["Data"][data_sets[i][2]]["x2"]=x2
		alldat["Data"][data_sets[i][2]]["x1_name"]=c1
		alldat["Data"][data_sets[i][2]]["x2_name"]=c2
	
	return(alldat)



'''
This subroutine reads in ASCII style zeus data. It takes three arguments, the filenama
of the data, and that of the r and theta files that give the grid cell coordnates.
The routine retyurns a dictionary that is identical to that computed from an HDF
file
'''




def get_zeus_data(fname,r_file,theta_file):

	if theta_file=='':
		dim=1
		print ("I think we are dealing with a 1D data file")
	else:
		dim=2
		print ("I think we are dealing with a 2D data file")



	inp=open(fname,"r")
	ir=[]
	itheta=[]
	data=[]

	raw_names=inp.readline()

	for line in inp.readlines():
		data_temp=line.split()
		print (len(data_temp))
		ir.append(int(data_temp[0]))
		itheta.append(int(data_temp[1]))
		temp=[]
		for i in range(len(data_temp)-2):
			temp.append(float(data_temp[i+2]))
		data.append(temp)
	inp.close()
	
	r=[]
	theta=[]
	
	if dim==2:
		inp=open(theta_file,"r")
		for line in inp.readlines():
			data_temp=line.split()
			try:
				if int(data_temp[0]) >=np.min(itheta) and int(data_temp[0]) <= np.max(itheta):
					theta.append(float(data_temp[2]))
				if int(data_temp[0]) == np.max(itheta):
					break
			except:
				print ("Something wrong with theta data file ",theta_file)
		inp.close()
	else:
		theta.append(np.pi/2.0) #We give theta a length of 1 so the reshape command still works.
	
	inp=open(r_file,"r")
	for line in inp.readlines():
		data_temp=line.split()
		try:
			if int(data_temp[0]) >=np.min(ir) and int(data_temp[0]) <= np.max(ir):
				r.append(float(data_temp[2]))
			if int(data_temp[0]) == np.max(ir):
				break
		except:
			print ("Something wrong with r data file ",r_file)
	inp.close()

	data_sets=["DENSITY","1-VELOCITY","2-VELOCITY","3-VELOCITY","TOTAL ENERGY"]
	
	alldat={}
	alldat["Filename"]=fname	
	alldat["Coord_sys"]="spol"
	alldat["Time"]=0
	alldat["Data"]={}
	alldat["N_data"]=5
	alldat["Dims"]=dim
	alldat["Data_names"]=data_sets
	

	


	for i in range (len(data_sets)):
		alldat["Data"][data_sets[i]]={}
		alldat["Data"][data_sets[i]]["long_name"]=data_sets[i]
		alldat["Data"][data_sets[i]]["data"]=np.reshape(np.array(data)[:,i],(len(theta),len(r)))
		alldat["Data"][data_sets[i]]["x1"]=np.array(theta)
		alldat["Data"][data_sets[i]]["x2"]=np.array(r)
		alldat["Data"][data_sets[i]]["x1_name"]=theta_file
		alldat["Data"][data_sets[i]]["x2_name"]=r_file
	
	return(alldat)



