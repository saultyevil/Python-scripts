#!/usr/bin/env python

"""
Provided by Nick Higginbottom
"""

from scipy.interpolate import griddata
from astropy import constants as c
import numpy as np
import pickle

#########################################
#
#streamline_2d returns an array of x,y coordinates from a dictionary containing x1,x2,v1,v2 and 
#the start point, which should be in cartesian coordinates x,y
#rmax is an optional value which defines the outer edge that we want to trak point to.
#tdisk is also optional, it defined where the disk starts, and if a streamline ends up here we truncate it.
#
#   v1 is the radial velocity
#   v2 is the velocity in the theta direction
#   x1 is theta
#   x2 is r
#   start is a 2 component array with the w,z starting point for the streamline
#   coord_sys at the moment only accepts spol. If any other systems are required, this will need coding up
#   rmax is an optional argument, it is the outer radius to which we want to track the streamline
#   tdisk is the grid cell number which defines the surface of the disk. -2 is the penultimate cell in the theta array
#   npoints is the number of points we will track the streamline for before abandoning it.
#
#######################################

def streamline_2d(v1,v2,x1,x2,start,coord_sys,rmax=-1,tdisk=-2,npoints=10000,opt='cart'):

#At the moment, this only works for spol coordinates. 
	if coord_sys!='spol':
		print("Error: Streamline 2d - unknown coordinate type ",coord_sys)


# If rmax is not supplied, or is greater than the maximum radius in the simulation, set it to just inside the outer radius.
	if rmax==-1 or rmax>=max(x2):
		rmax=max(x2)*0.999

# Set up two arrays to hold the r, theta coordinates of the streamline
	r=[]
	theta=[]	

#First, compute the theta coordinate of the start of the streamline
#start[1] is the y coordinate. We need to deal with the situation where a user asks for y=0, this will give an infinity 
	if start[1]==0.0:
		tstart=np.pi/2.0   #Deal with potential infinity
	else:
		tstart=np.arctan(start[0]/start[1])

#Check whether the user has sent a root location that is inside the disk. This is, by default, 2 cells from the end
	if tstart>x1[tdisk]:
		tstart=x1[tdisk]*0.9999999
		print("Warning: Streamline_2d - requested start point in disk - resetting to just above disk")
	
#Append the root location as the first point in the theta array
	theta.append(tstart)

#Now compute the starting radius
	rstart=np.sqrt(start[1]**2+start[0]**2)\

#Check that we are in the computational domain
	if rstart<x2[0]:
		rstart=x2[0]*1.00001
	elif rstart>rmax:
		rstart=rmax

#Append the root location as the first point in the r array
	r.append(rstart)



# Set ap a grid of length scales for the calculation - this is used to ensure that we choose sensible time steps for each cell.
#We will divide the length scale in a cell by the maximum velocity in the cell to give a time to cross the cell. We then divide by 10.
#We calculate the length of each side of the wedge shaped cell, then set the length scale to the smallest length

#Set up an empty array
	lengths=np.ones([len(x1),len(x2)])
#Loop over the r and theta arrays
	for j in range (len(x1)):
		for k in range(len(x2)):
			if j>0:	
				dl1=np.abs(x2[k]*np.cos(x1[j])-x2[k]*np.cos(x1[j-1]))
			else:
				dl1=np.abs(x2[k]*np.cos(x1[j])-x2[k]*np.cos(x1[j+1]))
			if k>0:
				dl2=np.abs(x2[k]-x2[k-1])
			else:
				dl2=np.abs(x2[k]-x2[k+1])

			lengths[j][k]=np.min([dl1,dl2])


	push=1.01   #This is a variable to try to ensure we dont end up exactly on a boundary
	status="tracking" #This flag is set to tracking while we are following a streamline. Any condition that means we want to stop following should result in this flag being changed to something else, and the loop being exited at the end.
	i=0
	while status=="tracking":	

#Compute the two velocity components by interpolation 
		tempv1=interpolate(x1,x2,theta[-1],r[-1],v1)
		tempv2=interpolate(x1,x2,theta[-1],r[-1],v2)

#Find out which cell we are in - this is just used to get the relvant length scale
		j,k,frac1,frac2=get_index(x1,x2,theta[-1],r[-1])

#We step across a cell my moving lengths given by the velocity x a time. This next line computes a time step for
#the current cell, simply by taking the smallest length, and dividing by the fastest velocity, and then by 10.	
		dt=(lengths[j][k]/np.max(np.abs([tempv1,tempv2])))/10.0

#Compute the new theta and r coordinates, just multiply the relevant velocities by the time step.
#We will check that these two 
		tnew=theta[-1]+(np.arctan(tempv2/r[-1])*dt)
		rnew=r[-1]+tempv1*dt


#We are now going to do a little testing. We want to make sure that if we enter a new cell, we dont jump too 
#far across that new cell.
#We begin by setting two new times, that we will use in a minute. Set the to the age of the universe as default
#NB, this was the age of the universe (to 1SF) as of 23rd March 2015. 
		r_time=4e17
		theta_time=4e17

#The following tests see if we will go up or down a cell in r or theta, and sets a time that will just about take us
#into the new cell. This means we won't race across a cell
		if tnew < x1[j] and theta[-1] != x1[j] :
			theta_time=push*(abs((x1[j]-theta[-1])/(np.arctan(tempv2/r[-1]))))
		if tnew > x1[j+1] and theta[-1] != x1[j+1]:
			theta_time=push*(abs((x1[j+1]-theta[-1])/(np.arctan(tempv2/r[-1]))))
		if rnew > x2[k+1]:
			r_time=push*(abs((x2[k+1]-r[-1])/tempv1))
		if rnew < x2[k]:
			r_time=push*(abs((x2[k]-r[-1])/tempv1))

#Compute new, new values of r and theta
		tnew=theta[-1]+(np.arctan(tempv2/r[-1])*min(dt,r_time,theta_time))
		rnew=r[-1]+tempv1*min(dt,r_time,theta_time)

#Now we tests to see if we have reached any of the limits that would cause us to stop tracking the streamline
#The message that the code outputs is pretty self explanatory 
		if tnew > x1[tdisk] :
			print("At upper theta limit (hit disk) after ",i," steps")
			if i>0:
				status="hit_disk"
			else:
				status="never_launched"
		elif tnew < x1[0] :
			print("At lower theta limit after ",i," steps")	
			status="hit_pole"
		elif rnew > rmax: 
			print("At outer radial boundary after ",i," steps")
			status="reached_outer"
		elif rnew < x2[0]:
			print("At inner radial boundary after ",i," steps")
			status="reached_inner"
		elif tnew==theta[-1] and rnew==r[-1]:
			print("Stalled after ",i," steps")
			status="stalled"
		elif i>npoints:
			print("Tracked for ",i," steps, aborting")
			status="too_many_steps"
		else:	
#If there are no stop conditions, append the new r and theta coordinates onto the arrays
			r.append(rnew)
			theta.append(tnew)
			i=i+1

#We want to return cartesian coordinates so loop over the r,theta coordinates and turn them back.
	x=[]
	z=[]
	for i in range(len(r)):
		x.append(r[i]*np.sin(theta[i]))
		z.append(r[i]*np.cos(theta[i]))
	if opt=='cart':
		return(np.array(x),np.array(z),status)
	elif opt=='vt':
		return(np.array(r),np.array(theta),status)

############################################################
#
#Multi stream produces a set of streamlines, defaulting to start at the second theta bin. This can be output as a
# pickled binary file, if a filename is specified, or returns the two arrays of streamx and streamw. This routine also
# returns the area along the streamlines, projected to the xy plane.
#The user supplies a dictionary which should contain r and theta coordinates, and vr, vtheta velocities.
#If the data is imported by the routines contained in input_sub, then the data will be correctly formatted.
#
#  wroot is an array of the w coordinates of streamlines
#  zroot is optional, if empty then the code computes z coordinates assuming the second theta bin is the start
#			if not, then it should contain the same number of coordinates as wroot. 
#  picklefile is the name of a file to write the streamlines out to, if empty, no picklefile is written
#  tstart is optional, it represents the number of the theta array that represents the cell above a disk
#  area_flag says whether the user wishes to compte the area of a stremline parallel to the w plane. This requires 
#			two extra streamlines, so triples the execution time
#  rmax is the maximum radius that we want to track a streamline out to
#  npoints, this is the number of points that the streamline can reach before we abandon it. 
#
##################################################################

def multi_stream (data,wroot,zroot=(),picklefile="",tstart=-2,area_flag="n",rmax=-1,npoints=100000,verbose='yes'):

# Extract the arrays we will be using from the dictionary
	v_r=data["Data"]["1-VELOCITY"]["data"]
	v_theta=data["Data"]["2-VELOCITY"]["data"]
	theta=data["Data"]["1-VELOCITY"]["x1"]
	r=data["Data"]["1-VELOCITY"]["x2"]
	
	print(r[0])

# If the user didnt define rmax, just set it to the outer edge
	if rmax==-1:
		rmax=max(r)

#We will be numbering the streamlines, the first is number zero
	nstream=0

#Set up empty arrays to hold all the streamlines

#First the coordinates
	streamw=[]
	streamz=[]
#Now the area of the streamline
	area_stream=[]
#We save the original roots
	roots=[]
#A flag to say wether we calculated areas
	area_flag1=[]
#The area can't necessarily be calculated for all z points, this tells us which z point first has an area
	area_min_z=[]
#The end condition of each streamline (hit disk, escaped, stalled etc)
	status=[]

#Do a loop over all requested roorts
	for i in range(len(wroot)):
#If we have zroots for all wroots, then use them
		if len(zroot) == len(wroot):
			start =(wroot[i],zroot[i])
#Otherwise calculate the zroot from the disk surface
		else:	
			start =(wroot[i],wroot[i]*(np.tan(np.pi/2.0-theta[tstart])))
#We try to generate a stream, this try command will crash if any of the commands before except fail. Obviously,
#this is most likely to be the call to streamline2d.
		try:
#Just to tell the user where we are - it can take a while to do this...
			if verbose=='yes':
				print("Multi stream: computing stream ",i," of ",len(wroot))
#This is the call to compute the streamline
			w,z,stat=streamline_2d(v_r,v_theta,theta,r,start,data["Coord_sys"],rmax=rmax,tdisk=tstart,npoints=npoints)
#Check if we want an area for the streamline
			if area_flag=="y" and len(w)>1:
#If so, remind the user what we are doing, and call stremline_area
				if verbose=='yes':
					print("Multi stream: computing stream area ",i," of ",len(wroot))
				area,min_z=streamline_area(v_r,v_theta,theta,r,start,data["Coord_sys"],w=w,z=z,tdisk=tstart,npoints=npoints)
#If not, fill the area array with dummy data
			else:
				area=[-1]
#If something went wrong, we hold the error code and report it, hopefully the user can work out what went wrong
		except Exception as e:
			if verbose=='yes':
				print("Error computing stream:",str(e))
			if len(w)==0:
				if verbose=='yes':
					print("Error multi stream - cannot compute streamlines for w0=",wroot[i])
				w=[]
				z=[]
				stat="error"
#This would be odd, it means that there was an error along the streamline, well, lets tell the user but keep the data
			else:
				if verbose=='yes':
					print("For some reason the streamline from w0=",wroot[i]," was truncated. Status=",stat) 
#If we got a nice streamline, tell the user how much work it was to do it, then copy the data to all of the 
#Arrays. NB, these end up as lists of arays/	
		if len(w)>1:
			if verbose=='yes':
				print("Stream ",i," has length ",len(w))
			roots.append(start)
			streamw.append(np.array(w))
			streamz.append(np.array(z))
			status.append(stat)
#If we asked for area, and got area, append all that data, set the flags accordingly.
			if area_flag=="y" and len(area)>0 and area[0]!=-1:
				area_stream.append(np.array(area))
				area_flag1.append(1)
				area_min_z.append(min_z)
			else:
				area_stream.append(np.array([-1]))
				area_flag1.append(0)
				area_min_z.append(-1)
		else:
			if verbose=='yes':
				print("Stream ",i," has length 0")
			roots.append(start)
			streamw.append(np.array(w))
			streamz.append(np.array(z))
			status.append(stat)
#Move onto the next streamline
		nstream=nstream+1

#Turn the data into arrays - easier to address in future.
	streamw=np.array(streamw)
	streamz=np.array(streamz)
	area_stream=np.array(area_stream)
	area_flag1=np.array(area_flag1)
	area_min_z=np.array(area_min_z)



#We will save this data - it can take a while to generate for lots of streamlines. first, we neatly package it in
#a dictionary.
	

	stream_data={'sw':streamw,'sz':streamz,'area':area_stream,'area_flag':area_flag1,'area_min_z':area_min_z}
	stream_data["roots"]=np.array(roots)
	stream_data["type"]='stream'
	stream_data["nstream"]=nstream
	stream_data["status"]=status
#If we have been given a filename, create the pickle file and save it. We dont test wether the file exists here, 
#so be careful to check in the calling routine. This one is dumb, it just overwrites.
	if picklefile != "":
		savefile=open(picklefile,'wb')
		pickle.dump(stream_data,savefile)
		savefile.close()
		return()
#If there wasnt a savefile - return the dictionary
	else:
		return (stream_data)

###########################################################
#
#Streamline_area is a routine to compute the area of a streamline.
#THe user supplies v theta, v r, r and theta, and possibly an already computed
#streamlime, then the code moves a small distance in and out in w, and computes two
#new streamlines. The area of this tube, parallel to the w plane is then computed.
#
#   v1 is the radial velocity
#   v2 is the velocity in the theta direction
#   x1 is theta
#   x2 is r
#   start is a 2 component array with the w,z starting point for the streamline
#   coord_sys at the moment only accepts spol. If any other systems are required, this will need coding up
#   w and z are optional, and contain a precomputed streamline. If not present, then the code
#		just calculates it, this is just a little timesaveer
#   tdisk is the grid cell number which defines the surface of the disk. -2 is the penultimate cell in the theta array
#   npoints is the number of points we will track the streamline for before abandoning it.
#
###########################################################

def streamline_area(v1,v2,x1,x2,start,coord_sys,w=[-1],z=[-1],tdisk=-2,npoints=100000):

#dr is the distance in and out that we generate new streamlines from. This is 1% of the current w root.
#I can imagnine that this hardwired number might well not suit all cases....

	dr=start[0]/100.0

#If we have not been given a stremline, call streamline_2d to generate on from the root
	if w[0]==[-1]:
		x,y=streamline_2d(v1,v2,x1,x2,(start[0],start[1]),coord_sys,tdisk=tdisk,npoints=npoints)
	else:
		x=w
		y=z
#Try to generate the inner streamline if we fail for some reason, report the error to the user and give up.
	try:
		xin,yin,status=streamline_2d(v1,v2,x1,x2,(start[0]-dr,start[1]/10.0),coord_sys,tdisk=tdisk,npoints=npoints)
	except Exception as e:
		print("Error: streamline_area computing inner streamline:",str(e))
		area=[-1.0]
		return(area)

#Try to generate the outer streamline if we fail for some reason, report the error to the user and give up.
	try:
		xout,yout,status=streamline_2d(v1,v2,x1,x2,(start[0]+dr,start[1]/10.0),coord_sys,tdisk=tdisk,npoints=npoints)
	except Exception as e:
		print("Error: streamline_area computing outer streamline:",str(e))
		area=[-1.0]
		return(area)

#This method of computing the streamline area only works if the streamlines are going in one direction. 
#So we do a quick check that all three lines do just that. This is pretty brute force, and will give up
#if you have a streamline that doubles back on itself. But is that a sensible streamline to analyze anyway??
#We just set a flag, and set it to zero if any points are smaller than the previous one!
	sort=1
	for i in range(len(yin)-1):
		if yin[i+1]<yin[i]:
			sort=0 
	for i in range(len(y)-1):
		if y[i+1]<y[i]:
			sort=0 
	for i in range(len(yout)-1):
		if yout[i+1]<yout[i]:
			sort=0 

#If our flag is still equal to 1, we are good. So lets work out the area
	if sort==1:
		ymin=max(y[0],yin[0],yout[0])   #This line finds the lowest y coordinate in all three streamlines. This is the first y  coordinate that we can get a sensible area for.
		area=[]
		imax=np.searchsorted(y,min(y[-1],yin[-1],yout[-1]))  #This line finds the last point in the streamline that has a corresponding y coordinate in both the inner and outer streamline

		for i in range(imax):
			if y[i]>ymin:
#This finds the point along the inner streamline with a y axis coordinate just above the y axis point of the central line
				i_in=np.searchsorted(yin,y[i])
#And this find the point in the outer line
				i_out=np.searchsorted(yout,y[i])
#This finds how far between the two inner points we need to be to exactly match the y axis point of the central line
				fracin=(y[i]-yin[i_in-1])/(yin[i_in]-yin[i_in-1])
#And for the outer
				fracout=(y[i]-yout[i_out-1])/(yout[i_out]-yout[i_out-1])
#So we work out the r coords of those inner and outer y points, interpolated
				rin=xin[i_in-1]+fracin*(xin[i_in]-xin[i_in-1])
				rout=xout[i_out-1]+fracout*(xout[i_out]-xout[i_out-1])
#And this is the area of that annulus
				area.append(4*np.pi*x[i]*(rout-rin))
#If the streamline doesnt always go up - we give up. This might be a bit wimpy - but needs a whole new way of computing area
	else:
		print("Area cannot be computed - chaotic streamlines")
		area=[-1]
		ymin=-1
		
#Return the area array, and the minimum value of y.

	return area,ymin

####################################################
#
# interpolate simply performs a 2D linear interpolation to find the value of data, at point
#    x1_coord, x2_coord, on a grid of x1, x2.
#  data is an array of the data
#  x1 is normally theta, x2 is r - this code is written assuming spol. Might need rewriting or
# extending if we want to use more coordinate types
#
###########################################3


def interpolate(x1,x2,x1_coord,x2_coord,data):

#First we call get_index to find the number of the grid that our coords are in. This is
#broken out into a seperate routine, because sometimes we just want to know the index, without
#interpolating.

	if len(x1)==1:  #Theta array is only 1 long, so we have 1D array
		j,k,frac1,frac2=get_index(x1,x2,x1_coord,x2_coord)
		d1=data[0][k]
		d2=data[0][k+1]
		d3=data[0][k]
		d4=data[0][k+1]	
	else:        # We have 2D data
		j,k,frac1,frac2=get_index(x1,x2,x1_coord,x2_coord)
		d1=data[j+1][k]
		d2=data[j+1][k+1]
		d3=data[j][k]
		d4=data[j][k+1]



#frac1 and frac2 are the fractional distance in theta and r. We now interpolarte first in
#theta, at the upper and lower r boundaries, then we interpolate those two points in r
#to obtain the final result.

	upper=d1+frac2*(d2-d1)
	lower=d3+frac2*(d4-d3)
	ans=lower+frac1*(upper-lower)

	return (ans)



###############################################################
#
# logterp simply performs a 2D logarithmic interpolation to find the value of data, at point
#    x1_coord, x2_coord, on a grid of x1, x2.
#  data is an array of the data
#  x1 is normally theta, x2 is r - this code is written assuming spol. Might need rewriting or
# extending if we want to use more coordinate types
#I'm not 100% sure that this work - use with caution
#
###########################################3


def logterp(x1,x2,x1_coord,x2_coord,data):


#	print "Warning: Logterp - the creator, NSH is not sure this works properly."

	j,k,frac1,frac2=get_index_log(x1,x2,x1_coord,x2_coord)
	if data[j+1][k]<0:
		sign=-1.0
	else:
		sign=1.0

	d1=np.log10(sign*data[j+1][k])
	d2=np.log10(sign*data[j+1][k+1])
	d3=np.log10(sign*data[j][k])
	d4=np.log10(sign*data[j][k+1])


	upper=d1+frac2*(d2-d1)
	lower=d3+frac2*(d4-d3)
	ans=sign*10.0**(lower+frac1*(upper-lower))

	return (ans)




####################################################
#
#     col_dens computes the column density back to the origin from a given point in polar coordinates.
#       most routines here take cartesian coords, but it makes more sense to use r theta, since all you
#      need to do is find theta, then run back in r
#   x1 should be theta
#   x2 should be r. 
#   x1_coord is the required theta coord we are going to - 
#   x2_coord is the required r coord
#   density is a 2d array of errrrrrm, density.
#
#################################################


def col_dens (x1,x2,x1_coord,x2_coord,density):

#First, generate the radial grid from the r dimensions. This gives us a list of radial widths of all the cells.

	dr=grid_gen(x2)[0]

#first find the cell we are in. Frac1 gives us the required interpolation distance in theta, frac2 is r
	
	j,k,frac1,frac2=get_index(x1,x2,x1_coord,x2_coord)

#now generate a 1D density by averaging the two adjacent density strips

	rho=density[j]+frac1*(density[j+1]-density[j])

#We start off by adding half the width of the innermost cell

	col_dens=rho[0]*0.5*dr[0]

#We now loop over cell 1 to cell k. Each time adding half the previous cell, and half the current cell

	for i in range(1,k+1):
		col_dens=col_dens+rho[i-1]*0.5*dr[i-1]+rho[i]*0.5*dr[i]

#Finally, we add on the last bit. This might be half the current cell, plus a bit, or a fraction of the current cell.

	if x2_coord-x2[k] > dr[k]:
		col_dens=col_dens+(rho[k]*0.5*dr[k])
		col_dens=col_dens+(rho[k+1]*frac2*(0.5*dr[k+1]))
	else:
		col_dens=col_dens+(rho[k]*frac2*0.5*dr[k])


	return (col_dens)


#################################################################
#
# Get index is a simple routine which just gives j,k where x1_coord 
#  lies between x1[j] and x1[j+1] and the same for k. It is used in interpolate, logterp
#   and any other place where we want to know which cell a coordinate is in and how far
#   it is from the boundaries in the two coordinates.
#    It was written and tested for spherical polar coordinates, but is probably general
#
#################################################################  

def get_index(x1,x2,x1_coord,x2_coord):




	if len(x1)==1: #Deals with 1D case where theta is not a variable
		j=0
		frac1=0.0
	else:      #Find the theta cell by just stepping up through the cells to find a bracket
		for j in range(len(x1)-1):
			if x1_coord>=x1[j] and x1_coord<=x1[j+1]:
				frac1=(x1_coord-x1[j])/(x1[j+1]-x1[j])
				break
	for k in range(len(x2)-1):  #Find the r cell in the same simple way
		if x2_coord>=x2[k] and x2_coord<x2[k+1]:
			frac2=(x2_coord-x2[k])/(x2[k+1]-x2[k])
			break

# These next lines were written to deal with problems when a point was essentially on the inner
# our outer boundary. This caused problems with this simple routine, and so the following lines
# just set it onto the boundary by brute force if it is close.

	if x1_coord-x1[0]<x1[0]/1e6: #There is some numerical problem....
		j=0
		frac1=0.0
	if x1[-1]-x1_coord<x1[-1]/1e6:
		j=len(x1)-2
		frac1=1.0
	if x2_coord-x2[0]<x2[0]/1e6:
		k=0
		frac2=0.0
	if x2[-1]-x2_coord<x1[-1]/1e6:
		k=len(x2)-2
		frac2=1.0




	return (j,k,frac1,frac2)


#################################################################
#
# get_index_log is a simple routine which just gives j,k where x1_coord 
#  lies between x1[j] and x1[j+1] and the same for k. It is used in interpolate, logterp
#   and any other place where we want to know which cell a coordinate is in and how far
#   it is from the boundaries in the two coordinates.
#    It was written and tested for spherical polar coordinates, but is probably general
#   this differs from get_index by interpolating in log space, this is used by logterp.
#
#################################################################


def get_index_log(x1,x2,x1_coord,x2_coord):


	

	for j in range(len(x1)-1):
		if x1_coord>=x1[j] and x1_coord<=x1[j+1]:
			frac1=(np.log10(x1_coord)-np.log10(x1[j]))/(np.log10(x1[j+1])-np.log10(x1[j]))
			break
	for k in range(len(x2)-1):
		if x2_coord>=x2[k] and x2_coord<x2[k+1]:
			frac2=(np.log10(x2_coord)-np.log10(x2[k]))/(np.log10(x2[k+1])-np.log10(x2[k]))
			break

	return (j,k,frac1,frac2)


#################################################################
#
# hc_rate computes the heating/cooling rate from PSK00, Blondin 94 etc, 
# returned as an array, so element 0 is the sum, but all the comonents are included 
#lf, cf, bf, and pf are optional terms that allow the user to artificially change the coefficient 
#of the various heating and cooling terms. See higginbottom and proga 15
#
# xi is the ionization parameter (not log xi)
# t is the temperature of the gas
# tx is the xray temperature
# T_l is the line temperature - defaults to 1.3e5 if it isn't included
# if the coefficients are omitted, then they default to 1
#
####################################################

def hc_rate(xi,t,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,T_l=1.3e5):
	g_comp=cf*(8.9e-36*xi*(tx-4.0*t))				#Compton heating and cooling
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))		#Xray heating (photoionization/recomb)
	lb=-1*bf*(3.3e-27*t**0.5)					#Bremstrahlung
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5+1e-24)	#Line
	hcrate=g_comp+gx+lb+ll

	return (np.array([hcrate,g_comp,gx,lb,ll]))


####################################################################
#
# hc_rate2 is just like hc_rate but splits up the heating and cooling components of X-ray and compton.
# I'm sure this could be done more cleanly
# It also included bremstrahlung heating. If the density is set to other than the default of
# 0.0 then this is computed and returned.
#
##########################################################################

def hc_rate2(xi,t,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,T_l=1.3e5,n_dens=0.0):
#	print xi,t,tx
	comp_heat=cf*(8.9e-36*xi*(tx))					#Compton heatinng
	comp_cool=-1.0*cf*(8.9e-36*xi*(4.0*t))				#Compton cooling
	pi_heat=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0))			#Xray heating
	rec_cool=pf*-1.0*(1.5e-21*xi**0.25*t**(-0.5)*(t/tx))		#Xray cooling
	lb=-1*bf*(3.3e-27*t**0.5)					#Bremstrahlung cooling
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5) #Bremstrahlung heating
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5+1e-24)	#Line
	hcrate=comp_heat+comp_cool+pi_heat+rec_cool+lb+ll
	return (np.array([hcrate,comp_heat,comp_cool,pi_heat,rec_cool,brem_heat,lb,ll]))

def hc_rate3(xi,t,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,T_l=1.3e5,n_dens=0.0):
#	print xi,t,tx
	comp_heat=cf*(8.9e-36*xi*(tx))					#Compton heatinng
	comp_cool=-1.0*cf*(8.9e-36*xi*(4.0*t))				#Compton cooling
	pi_heat=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))			#Xray heating
	lb=-1*bf*(3.3e-27*t**0.5)					#Bremstrahlung cooling
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5) #Bremstrahlung heating
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5+1e-24)	#Line
	hcrate=comp_heat+comp_cool+pi_heat+lb+ll
	return (np.array([hcrate,comp_heat,comp_cool,pi_heat,brem_heat,lb,ll]))


def hc_rate4(xi,t,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,T_l=1.3e5,n_dens=0.0):
#	print xi,t,tx
	comp_heat=cf*(8.9e-36*xi*(tx))					#Compton heatinng
	comp_cool=-1.0*cf*(8.9e-36*xi*(4.0*t))				#Compton cooling
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)
	ll=-1*lf*(1e-16*np.exp(-T_l/t)*xi**-0.5*t**-1.0)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)
	#	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)
	hcrate=comp_heat+comp_cool+gx+lb+ll
	return (np.array([hcrate,comp_heat,comp_cool,gx,lb,ll]))


#########################################################################
#
# hc_rate1 is a brentq friendly version of the routine, which has t as the first argument, and only returns the rate.
# it can be used to compute an equilibrium temperature using brentq or similar root finding algorithm
#
########################################################################


def hc_rate1(t,xi,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,const=0.0,n_dens=0.0,T_l=1.3e5):
	g_comp=cf*(8.9e-36*xi*(tx-4.0*t))
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5)
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5+1e-24)
	hcrate=g_comp+gx+lb+ll+const+brem_heat

	return (hcrate)
	
	#########################################################################
	#
	# hc_rate1a is a brentq friendly version of the routine, which has t as the first argument, and only returns the rate.
	# it can be used to compute an equilibrium temperature using brentq or similar root finding algorithm. This is modified
	# over hc_rate1a to take account of a misunderstanding in how the prefactor should be applied to line cooling
	#
	########################################################################


	
def hc_rate1a(t,xi,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,const=0.0,n_dens=0.0,T_l=1.3e5):
	g_comp=cf*(8.9e-36*xi*(tx-4.0*t))
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5)
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5)+1e-24
	hcrate=g_comp+gx+lb+ll+const+brem_heat

	return (hcrate)




	#########################################################################
	#
	# hc_rate1b is a brentq friendly version of the routine, which has t as the first argument, and only returns the rate.
	# it can be used to compute an equilibrium temperature using brentq or similar root finding algorithm. This is modified
	# over hc_rate1a to take account of a misunderstanding in how the prefactor should be applied to line cooling
	#
	########################################################################


	
def hc_rate1b(t,xi,tx,lf=1.0,chf=1.0,ccf=1.0,bf=1.0,pf=1.0,const=0.0,n_dens=0.0,T_l=1.3e5):
	g_comp=chf*(8.9e-36*xi*tx)
	l_comp=ccf*(-1*8.9e-36*xi*4.0*t)
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5)
	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5)+1e-24
	hcrate=g_comp+l_comp+gx+lb+ll+const+brem_heat

	return (hcrate)




	#########################################################################
	#
	# hc_rate1c is a brentq friendly version of the routine, which has t as the first argument, and only returns the rate.
	# has wierd new line code
	#
	########################################################################


	
def hc_rate1c(t,xi,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,const=0.0,n_dens=0.0,T_l=1.3e5):
	g_comp=cf*(8.9e-36*xi*(tx-4.0*t))
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5)
	ll=-1*lf*(1e-16*np.exp(-T_l/t)*xi**-0.5*t**-1.0)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)
#	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)
	hcrate=g_comp+gx+lb+ll+const+brem_heat

	return (hcrate)

	#########################################################################
	#
	# hc_rate1d is a brentq friendly version of the routine, which has t as the first argument, and only returns the rate.
	# has wierd new line code and also a variable ne_vs_nh 
	#
	########################################################################


	
	
def hc_rate1d(t,xi,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0,const=0.0,n_dens=0.0,T_l=1.3e5):
	
	x1=10**(-51.59417133+12.27740153*np.log10(t))
	x2=10**(-3.80749689+0.86092628*np.log10(t))
	ne=np.min([x1,x2,1.21])
		
	
	
	g_comp=cf*(8.9e-36*xi*(tx-4.0*t))*ne
	gx=pf*(1.5e-21*xi**0.25*t**(-0.5)*(1.0-t/tx))
	lb=-1*bf*(3.3e-27*t**0.5)*ne
	brem_heat=3.3e-27*xi*n_dens/(4.0*np.pi)/(c.sigma_sb.cgs.value)*t**(-3.5)*ne
	ll=-1*lf*(1e-16*np.exp(-T_l/t)*xi**-0.5*t**-1.0)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)*ne
#	ll=-1*lf*(1.7e-18*np.exp(-T_l/t)*xi**-1.0*t**-0.5)+min(1e-24,5e-27*np.sqrt(t),1.5e-17/t)
	hcrate=g_comp+gx+lb+ll+const+brem_heat

	return (hcrate)	
	
	x1=10**(-51.59417133+12.27740153*log10(temp))
	x2=10**(-3.80749689+0.86092628*log10(temp))
	
	
	

###############################################################################
#
#
#   A routine to find the value of T where L=0 for a given xi. 
#   The lf,cf,bf,pf parameters allow the user to apply a scaling factor to the heating and 
#  cooling processes for testing purposes. Leave them out and all are set to 1.0.
#   This is largely superceeded in most of my codes by the use of brentq with 
#   hc_rate1 as the function
#
#################################################################################

def T_equil(xi,tx,lf=1.0,cf=1.0,bf=1.0,pf=1.0):

	icalls=0
	roots=[]
	tol=1.0


	t_array=np.logspace(3,10,10)
	for i in range(len(t_array)-1):
		icalls=icalls+1
		if hc_rate(xi,t_array[i],tx,lf=lf,cf=cf,bf=bf,pf=pf)[0]*hc_rate(xi,t_array[i+1],tx)[0] < 0.0:
			roots.append(i)

#	print "There are ",len(roots)," roots"
#	for i in range(len(roots)):
#		print "Root ",i+1," is between t=",t_array[roots[i]]," and ",t_array[roots[i]+1]

	if len(roots) ==0:
		print("There is no equilibrium for xi=",xi)
		return 0.0
	try: tupper=t_array[roots[0]+1]
	except:
		print(roots)

	tlower=t_array[roots[0]]
	gap=tupper-tlower
	
	while gap>tol:
		tmiddle=(tupper+tlower)/2.0
		if hc_rate(xi,tupper,tx,lf=lf,cf=cf,bf=bf,pf=pf)[0]*hc_rate(xi,tmiddle,tx,cf=cf,bf=bf,pf=pf)[0] < 0.0:
			icalls=icalls+2
			tlower=tmiddle
		elif hc_rate(xi,tlower,tx,lf=lf,cf=cf,bf=bf,pf=pf)[0]*hc_rate(xi,tmiddle,tx,cf=cf,bf=bf,pf=pf)[0] < 0.0:
			icalls=icalls+2
			tupper=tmiddle
		else:
			return ans
		gap=tupper-tlower
		
#	print "Extied loop with tupper=",tupper," and tlower=",tlower
#	print icalls
	ans=(tupper+tlower)/2.0


	return ans




##############################################################
#
# rec_grid produces a rectangular gridded data set from r, theta, data inputs. 
# this is to allow plotting using a regular countour plot in python
#
#######################################################

def rec_grid(r,theta,data):

	w=[]
	z=[]
	data1=[]


	for i in range(len(theta)):
		for j in range(len(r)):
			w.append(r[j]*np.sin(theta[i]))
			z.append(r[j]*np.cos(theta[i]))
			data1.append(data[i][j])
	w=np.array(w)
	z=np.array(z)
	data1=(np.array(data1))
	wi=np.logspace(np.log10(np.min(w)),np.log10(np.max(w)),1000)
	zi=np.logspace(np.log10(np.min(z)),np.log10(np.max(z)),1000)
	datai=griddata((w,z),data1,(wi[None,:],zi[:,None]),method='linear')
	print(len(datai),len(datai[0]))
	for i in range(len(datai)):
		for j in range(len(datai[0])):
			if (np.sqrt(wi[j]*wi[j]+zi[i]*zi[i]) < np.min(r)) or (np.sqrt(wi[j]*wi[j]+zi[i]*zi[i]) > np.max(r)) or np.arctan(wi[j]/zi[i]) > np.max(theta) or np.arctan(wi[j]/zi[i]) < np.min(theta):
				datai[i][j]=-999				
	data2=np.ma.masked_equal(datai,-999)
	w=wi	
	z=zi


	return (w,z,data2,data1)


##########################################################
#
#This subroutine calculates the gradient for a 2d array data.
#It returns two arrays, one for the r direction, and one for theta, 
#  both are in change in data per cm  
#
#############################################################


def pot_grad(r,theta,data):

	ddt,ddr=np.gradient(data)
#We now have two array, each is the change per grid, we now just divide by the actual distance.

	for i in range(len(theta)):
		for j in range(len(r))	:
			if j==0:
				ddr[i][j]=ddr[i][j]/(r[1]-r[0])
			elif j==len(r)-1:
				ddr[i][j]=ddr[i][j]/(r[-1]-r[-2])
			else:
				ddr[i][j]=ddr[i][j]/((r[j+1]-r[j-1])/2.0)
			if i==0:
				ddt[i][j]=ddt[i][j]/(r[j]*(theta[1]-theta[0]))
			elif i==len(theta)-1:
				ddt[i][j]=ddt[i][j]/(r[j]*(theta[-1]-theta[-2]))
			else:
				ddt[i][j]=ddt[i][j]/(r[j]*(theta[i+1]-theta[i-1])/2.0)

	return ddt,ddr


#########################################################
#
# outer_root_gen generates streamline roots on the outer edge of a simulation based upon 
# density and mass loss rate. These roots can then be used to produce streamlines 
# by sending the roots to multistream, with the relvant options. 
# the ignore parameter allows one to disregard the last theta bin if this is the disk...
# It produces 20 roots by default, but this can be overridden by setting the npoints 
# option
#
#  Inputs:
#		data - a dictionary containing the data to be used
#  		npoints (optional) - the number of roots to be generated
#		ignore (optional) the number of cells close to the disk to ignore
#
#############################################################


def outer_root_gen(data,npoints=20,ignore=0):

	v_r=data["Data"]["1-VELOCITY"]["data"]
	r=data["Data"]["1-VELOCITY"]["x2"]
	theta=data["Data"]["1-VELOCITY"]["x1"]
	density=data["Data"]["DENSITY"]["data"]


	outer_density=[]
	cum_density=[]
	mass_flux=[]
	cum_mass_flux=[]
	vr=[]

	if ignore==0:
		tmax=np.max(theta)
	else:
		tmax=np.max(theta[:-ignore])

	
	#Generate a fine mesh of theta points
	theta1=np.linspace(np.min(theta),tmax,10001)

	#In general, starting streamlines exactly at the outer surface causes problems
	r_test=np.max(r)*0.999
	dtheta=theta1[1]-theta1[0]
	#Generate arrays of all the cumlative quantities, this is pretty cheap, so lets just make them all
	for i in range(len(theta1)):
		theta_test=theta1[i]
		#Density
		d1=interpolate(theta,r,theta_test,r_test,density)
		#Radial velocity
		vr.append(interpolate(theta,r,theta_test,r_test,v_r))
		#This is the flux of mass moving out of the simulation at this point
		mass_flux.append(vr[-1]*d1)
		outer_density.append(d1)
		if i==0: #This is the first point
			cum_density.append(d1)
			cum_mass_flux.append(mass_flux[-1]*np.sin(theta1[0])*dtheta )
		else:
			cum_density.append(d1+cum_density[-1])
			cum_mass_flux.append(mass_flux[-1]*np.sin(theta1[i])*dtheta+cum_mass_flux[-1])

	cum_mass_flux=np.array(cum_mass_flux)*max(r)*max(r)*4.0*np.pi
	#The total mass flux, and density is just the last point in the respective arrays
	rho_tot=cum_density[-1]
	mass_flux_tot=cum_mass_flux[-1]

	#The steps in density and outgoing mass flux
	drho=rho_tot/npoints
	dmf=mass_flux_tot/npoints
	
	#Empty arrays for the roots, based in density at the outer edge
	wroot_rho=[]
	zroot_rho=[]

	#Empty arrays for the roots based on outgoing mass flux
	wroot_mf=[]
	zroot_mf=[]


	j=1.0
	k=1.0
	for i in range(len(theta1)):
		if cum_density[i]>drho*j:
			j=j+1.0
			wroot_rho.append(r_test*np.sin(theta1[i]))
			zroot_rho.append(r_test*np.cos(theta1[i]))
		if cum_mass_flux[i]>dmf*k:
			k=k+1.0
			wroot_mf.append(r_test*np.sin(theta1[i]))
			zroot_mf.append(r_test*np.cos(theta1[i]))

	roots={}
	#we return the roots
	roots["wroot_rho"]=wroot_rho
	roots["zroot_rho"]=zroot_rho
	roots["wroot_mf"]=wroot_mf
	roots["zroot_mf"]=zroot_mf
	#We also return the data used to generate the roots (why not??)
	roots["theta_grid"]=theta1
	roots["outer_rho"]=outer_density
	roots["cum_rho"]=cum_density
	roots["mass_flux"]=mass_flux
	roots["v_r"]=vr
	roots["cum_mass_flux"]=cum_mass_flux
	

	return roots
			

###########################################
#
# 1d_streamline produces a streamline type file but for a 1D case - it is just a set of r and areas, assuming area goes as r**2
#
#  Inputs
#		data - a dictionary containing the datafile
#		picklefile - the name of the picklefile in which to store the points describing the streamline
#		gridding (optional) - the number of points to generate between each actual radius - allows finer grids
#
#######################################################

def oned_streamline(data,picklefile,gridding=1):
	
	r=data["Data"]["1-VELOCITY"]["x2"]

	streamw=[]
	streamz=[]
	area_stream=[]
	
	for i in range(len(r)-1):
		for j in range(gridding):
			streamw.append(r[i]+(r[i+1]-r[i])/float(gridding)*float(j))
			streamz.append(0.0)
			if i==0:
				area_stream.append(1.0)
			else:
				area_stream.append(streamw[-1]**2/r[0]**2)

	streamw=np.array(streamw)
	streamz=np.array(streamz)
	area_stream=np.array(area_stream)




	stream_data={'sw':streamw,'sz':streamz,'area':area_stream,'area_flag':[1],'nstream':1}
	stream_data["type"]='oned'
	savefile=open(picklefile,'wb')
	pickle.dump(stream_data,savefile)
	savefile.close()
	return()


###################################################################
#
#  Sightline generates a set of data points which define a series of 
#  sightlines along which one wishes to make plots of physical parameters.
#  
#  Inputs
#		data - a dictionary containing the data
#		picklefile (optional)- the name of the picklefile that we are going to store the 
#					sightline data in. If it 
#		angles - an array containing the angles we want to compute data for
#		gridding (optional) - number of points between each radial coordinate
#
######################################################################

def sightline_gen(data,picklefile,angles,gridding=1):

	theta=data["Data"]["1-VELOCITY"]["x1"]
	r=data["Data"]["1-VELOCITY"]["x2"]


	roots=[]
	streamw=[]
	streamz=[]
	area_flag1=[]
	area=[]	
	area_min_z=[]


	for i in range(len(angles)):
		w=[]
		z=[]
		for j in range(len(r)):
			w.append(r[j]*np.sin(np.radians(angles[i])))
			z.append(r[j]*np.cos(np.radians(angles[i])))
		streamw.append(w)
		streamz.append(z)
		area_flag1.append(0)
		area.append(0.0)
		area_min_z.append(-1.0)
		roots.append([w[0],z[0]])


	if picklefile != "":
		stream_data={'sw':streamw,'sz':streamz,'area_flag':area_flag1,'area':area,'area_min_z':area_min_z}
		stream_data["type"]="sight" 
		stream_data["roots"]=np.array(roots)
		stream_data["nstream"]=len(angles)
		stream_data["angles"]=np.array(angles)
		savefile=open(picklefile,'wb')
		pickle.dump(stream_data,savefile)
		savefile.close()
		return()

	else:
		return (streamw,streamz,area_stream,area_flag1)

###############################################
#
#  Smooth Array - a simple smoothing algorithm 
#
#    Inputs:
#		array - an array of data
#		bin (optional) - the bin size to smmoth over
#
#############################################



def smooth_array(array,bin=20):
		smooth=[]
		for n in range(int(bin/2),len(array)-int(bin/2)):
			temp=0.0
			for i in range(bin):
				temp=temp+array[n-int(bin/2)+i]
			smooth.append(temp/float(bin))
		
		return np.array(smooth)

##################################################
#
# const_ang_mom_stream produces a streamline from a point based on equation 2.1 in Icke (1981)
# It takes as inputs just the starting radius and the maximum radius that we want to track the streamline to
#
#################################################

def const_ang_mom_stream(R,rmax):

	r_array=np.logspace(np.log10(R),np.log10(rmax),1001)
	 
	t1=(r_array/R)**(2.0/3.0)
	t1=np.sqrt(t1-1.0)	

	h_array=t1*r_array

	return r_array,h_array


##################################################################
#
# grid_gen takes a set of coordinates, and works out the grid. 
# It returns the grid spacing, min, max, and the ratio of each
# grid point to the previous one
#
###############################################################
		
def grid_gen(x):

	x_ratio=(x[2]-x[1])/(x[1]-x[0])
	dx=[]
	dx.append((x[1]-x[0])/(0.5*(1.0+x_ratio)))
	xmin=x[0]-0.5*dx[-1]
	xmax=xmin+dx[-1]*(1.0-x_ratio**len(x))/(1.0-x_ratio)
	for i in range(len(x)-1):
		dx.append(x_ratio*dx[-1])
	
	return np.array(dx),xmin,xmax,x_ratio
	
##################################################################
#
# athena_to_zeus takes an athena data structure and converts it
# to be just like a zeus structure, so all the scripts work
#
###############################################################	
	
def athena_to_zeus(data,j=0,gamma=5./3.,type="hdf"):
	
	
	r=[]
	theta=[]

	


	if type=="hdf":
		v_r=data["vel1"][j]
		v_theta=data["vel2"][j]
		v_phi=data["vel3"][j]
		density=data["rho"][j]
		pressure=data["press"][j]
		for i in range(len(data["x1f"])-1):
			r.append((data["x1f"][i]+data["x1f"][i+1])/2.0)

		for i in range(len(data["x2f"])-1):
			theta.append((data["x2f"][i]+data["x2f"][i+1])/2.0)
	elif type=="vtk":
		v_r=data[3]['vel'][0,:,:,0]
		v_theta=data[3]['vel'][0,:,:,1]
		v_phi=data[3]['vel'][0,:,:,2]
		density=data[3]["rho"][0]
		pressure=data[3]["press"][0]
		for i in range(len(data[0])-1):
			r.append((data[0][i]+data[0][i+1])/2.)
		for i in range(len(data[1])-1):
			theta.append((data[1][i]+data[1][i+1])/2.0)
	else:
		print("Unknown type of file")
		
		
	r=np.array(r)
	theta=np.array(theta)
	energy=pressure/(gamma-1.)
	
	
	
	data1={}
	data1["Data"]={}
	data1["Data"]["DENSITY"]={}
	data1["Data"]["3-VELOCITY"]={}
	data1["Data"]["2-VELOCITY"]={}
	data1["Data"]["1-VELOCITY"]={}
	data1["Data"]["TOTAL ENERGY"]={}
	

	data1["Data"]["DENSITY"]["data"]=density
	data1["Data"]["DENSITY"]["x2"]=r
	data1["Data"]["DENSITY"]["x1"]=theta

	data1["Data"]["3-VELOCITY"]["data"]=v_phi
	data1["Data"]["3-VELOCITY"]["x2"]=r
	data1["Data"]["3-VELOCITY"]["x1"]=theta


	data1["Data"]["2-VELOCITY"]["data"]=v_theta
	data1["Data"]["2-VELOCITY"]["x2"]=r
	data1["Data"]["2-VELOCITY"]["x1"]=theta
	
	
	data1["Data"]["1-VELOCITY"]["data"]=v_r
	data1["Data"]["1-VELOCITY"]["x2"]=r
	data1["Data"]["1-VELOCITY"]["x1"]=theta

	data1["Coord_sys"]="spol"
	
	data1["Data"]["DENSITY"]["long_name"]="Athena file"
	
	data1["Data"]["TOTAL ENERGY"]["data"]=energy
	data1["Data"]["TOTAL ENERGY"]["x2"]=r
	data1["Data"]["TOTAL ENERGY"]["x1"]=theta
	
	
	return (data1)


##################################################################
#
# pluto_to_zeus takes an pluto dbl data structure and converts it
# to be just like a zeus structure, so all the scripts work
#
###############################################################	


def pluto_to_zeus(data,j=0,gamma=5./3.,type="hdf"):
	
	
	r=[]
	theta=[]
	v_r=data.vx1
	v_theta=data.vx2
	v_phi=data.vx3
	density=data.rho
	r=data.x1
	theta=data.x2
	

	
	
	data1={}
	data1["Data"]={}
	data1["Data"]["DENSITY"]={}
	data1["Data"]["3-VELOCITY"]={}
	data1["Data"]["2-VELOCITY"]={}
	data1["Data"]["1-VELOCITY"]={}
	data1["Data"]["TOTAL ENERGY"]={}
	data1["Data"]["PRESSURE"]={}
	data1["Data"]["LINE_C"]={}
	data1["Data"]["XRAY_H"]={}
	data1["Data"]["COMP_C"]={}
	data1["Data"]["COMP_H"]={}
	data1["Data"]["BREM_C"]={}
	
	data1["Data"]["LINE_C_PRE"]={}
	data1["Data"]["XRAY_H_PRE"]={}
	data1["Data"]["COMP_C_PRE"]={}
	data1["Data"]["COMP_H_PRE"]={}
	data1["Data"]["BREM_C_PRE"]={}
	
	data1["Data"]["T"]={}
	data1["Data"]["XI"]={}
	data1["Data"]["ne"]={}
	data1["Data"]["nh"]={}
	
	try:
		T=data.T
		xi=data.XI
		comp_c=data.comp_c
		comp_h=data.comp_h
		line_c=data.line_c
		brem_c=data.brem_c
		xray_h=data.xray_h		
		data1["Data"]["LINE_C"]["data"]=np.transpose(line_c)
		data1["Data"]["XRAY_H"]["data"]=np.transpose(xray_h)
		data1["Data"]["COMP_C"]["data"]=np.transpose(comp_c)
		data1["Data"]["COMP_H"]["data"]=np.transpose(comp_h)
		data1["Data"]["BREM_C"]["data"]=np.transpose(brem_c)
		data1["Data"]["T"]["data"]=np.transpose(T)
		data1["Data"]["XI"]["data"]=np.transpose(xi)
	except:
		print("No cooling data")
		data1["Data"]["LINE_C"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["XRAY_H"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["COMP_C"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["COMP_H"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["BREM_C"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["T"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["XI"]["data"]=np.zeros(np.shape(np.transpose(density)))
		
		
	try:
		data1["Data"]["LINE_C_PRE"]["data"]=np.transpose(data.line_c_pre)
		data1["Data"]["XRAY_H_PRE"]["data"]=np.transpose(data.xray_h_pre)
		data1["Data"]["COMP_C_PRE"]["data"]=np.transpose(data.comp_c_pre)
		data1["Data"]["COMP_H_PRE"]["data"]=np.transpose(data.comp_h_pre)
		data1["Data"]["BREM_C_PRE"]["data"]=np.transpose(data.brem_c_pre)	
	except:
		print("No prefactor data")
		data1["Data"]["LINE_C_PRE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["XRAY_H_PRE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["COMP_C_PRE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["COMP_H_PRE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["BREM_C_PRE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		
	try:
		ne=data.ne
		nh=data.nh		
		data1["Data"]["ne"]["data"]=np.transpose(ne)
		data1["Data"]["nh"]["data"]=np.transpose(nh)
	except:
		data1["Data"]["ne"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["nh"]["data"]=np.zeros(np.shape(np.transpose(density)))
	
	
	
	try:
		pressure=data.prs
		energy=pressure/(gamma-1.)
		data1["Data"]["PRESSURE"]["data"]=np.transpose(pressure)
		data1["Data"]["TOTAL ENERGY"]["data"]=np.transpose(energy)
	except:
		print ("No pressure, ISOTHERMAL")
		data1["Data"]["PRESSURE"]["data"]=np.zeros(np.shape(np.transpose(density)))
		data1["Data"]["TOTAL ENERGY"]["data"]=np.zeros(np.shape(np.transpose(density)))
		
		
	data1["Data"]["DENSITY"]["data"]=np.transpose(density)
	data1["Data"]["DENSITY"]["x2"]=r
	data1["Data"]["DENSITY"]["x1"]=theta

	data1["Data"]["3-VELOCITY"]["data"]=np.transpose(v_phi)
	data1["Data"]["3-VELOCITY"]["x2"]=r
	data1["Data"]["3-VELOCITY"]["x1"]=theta


	data1["Data"]["2-VELOCITY"]["data"]=np.transpose(v_theta)
	data1["Data"]["2-VELOCITY"]["x2"]=r
	data1["Data"]["2-VELOCITY"]["x1"]=theta
	
	
	data1["Data"]["1-VELOCITY"]["data"]=np.transpose(v_r)
	data1["Data"]["1-VELOCITY"]["x2"]=r
	data1["Data"]["1-VELOCITY"]["x1"]=theta


	
	return (data1)


