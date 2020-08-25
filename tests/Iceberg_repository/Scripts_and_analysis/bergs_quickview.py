#!/usr/bin/env python

def bergs_yearday(bergs):
	bergs=bergs_read(bergs,'year')
	bergs=bergs_read(bergs,'day')
	for b in xrange (0,len(bergs.berg)):
		bergs.berg[b].yearday=bergs.berg[b].year+(bergs.berg[b].day/372)
		#Note I have not included all the warnings about dimensions already existing.
        return bergs


def bergs_read(bergs,varname):
	nc=bergs.nc
	bergs_list=bergs.berg
	var = nc.variables[varname][:]
	for b in xrange (0,len(bergs.berg)):
		k=[]
		js=bergs.berg[b].js
		je=bergs.berg[b].je
		for l in xrange (0,len(js)):
			temp=range(js[l],je[l]+1)
			for i in xrange (0, len(temp)):
				k.append(temp[i])
			setattr(bergs.berg[b],varname,var[k])
        return bergs


def bergs_hash(year, day, x, y, m):
	ilat=(90.+y)/180
        ilon=(360 % x)/360;
        iyear=np.floor(year);
        iday=(day/372);
        im=math.log(m,10)-7.9;
        h=ilon+1e3*ilat+1e6*iday+1e9*iyear+1e12*im;
        return h



def bergs_open(filename):

	#Creating a classes
	class berg_handle:
		def __init__(self,berg,nc):
			self.berg=berg
			self.nc=nc

	class berg_type:
		def __init__(self,hash):
			self.hash=hash
			self.year0=[]
			self.day0=[]
			self.lon0=[]
			self.lat0=[]
			self.mass0=[]
			self.js=[]
			self.je=[]
			self.width=[]
			self.lon=[]
			self.lat=[]
			self.year=[]
			self.day=[]
			self.uvel=[]
			self.vvel=[]
			self.uo=[]
			self.vo=[]
			self.ui=[]
			self.ua=[]
			self.vu=[]
			self.mass_of_bits=[]
			self.heat_density=[]
			self.thickness=[]
			self.length=[]
			self.ssh_x=[]
			self.ssh_y=[]
			self.sst=[]
			self.cn=[]
			self.hi=[]
			self.yearday=[]
			self.iceberg_num=[]

	#Importing the data
	nc = Dataset(filename, mode='r')
	lat = nc.variables['lat'][:]
	N=len(lat)
	js=np.where(lat>91)[0]
	je=np.concatenate([js[1:]-1, np.array([N-1])])
	lon = nc.variables['lon'][:]
	year = nc.variables['year'][:]
	day = nc.variables['day'][:]
	#mass = nc.variables['mass'][:]
	iceberg_num = nc.variables['iceberg_num'][:]

	N=len(js)
	h=np.zeros(N)
	# Making a list of identifying hashes
	k=0
	for j in js+1:
		#h[k]=bergs_hash(year[j], day[j],lon[j],lat[j],mass[j])
		h[k]=iceberg_num[j]
		k=k+1

	#Creates a sorted array of unique elements of h
	h=sorted(list(set(h)))

	#Create a list of bergs, named by their hash values, h
	berg=[]
	for i in xrange (0,len(h)):
		berg.append(berg_type(h[i]))

	for i in xrange(0, len(js)):
		j=js[i]+1
		#l=np.where(h==bergs_hash(year[j], day[j],lon[j],lat[j],mass[j]))[0]
		l=np.where(h==iceberg_num[j])[0]
		berg[l].year0.append(year[j])
		berg[l].day0.append(day[j])
		berg[l].lon0.append(lon[j])
		berg[l].lat0.append(lat[j])
		#berg[l].mass0.append(mass[j])
		berg[l].js.append(js[i]+2)
		berg[l].je.append(je[i])

	#Sort each segment
	for i in xrange (0,len(h)):
		yd=year[berg[i].js]+(day[berg[i].js]/373)
		#Sort according to descending yd, and then pull out indicies for js, and je.
		
		berg[i].js=[berg[i].js[n] for n,_ in sorted(list(enumerate(yd)),reverse=True)]
		berg[i].je=[berg[i].je[n] for n,_ in sorted(list(enumerate(yd)),reverse=True)]

	bergs_handle=berg_handle(berg,nc)

	return bergs_handle

###############################################################################################
#################################  Beginning of Script  #######################################
###############################################################################################
from netCDF4 import Dataset
import numpy as np
#from PIL import *
import math
import os
import matplotlib
matplotlib.use("GTKAgg")
import matplotlib.pyplot as plt
from pylab import *
import pdb
import argparse
#This combination allows you to access the varibale imputted from the command line.
import sys

#Clear screen
#os.system('clear')

#This is a command which puts you into testing mode so that you can troubleshoot (like a breakpoint in matlab)
#pdb.set_trace()

#If no files are given at runtime, then use the files selected here:
if len(sys.argv)==1:
	Number_of_files=2  # How many files would you like to view?

	#Choosing which files to use
	#Weddel Sea Version
	#fileroot='/home/aas/Iceberg_Project/Weddell_Sea/'
	#dirname1='one_iceberg_testing_Verlet';
	#dirname2='one_iceberg_testing_RK';
	
	#filename1=fileroot+dirname1+'/iceberg_trajectories.nc'
	#filename2=fileroot+dirname2+'/iceberg_trajectories.nc'

	#Other files
	#filename1='/home/aas/Iceberg_Project/Test_data1/Quick_Global_RK15/iceberg_trajectories.nc';
	#filename1= '/home/aas/Iceberg_Project/Weddell_Sea/one_iceberg_testing/iceberg_trajectories.nc'
	filename1= '../Rolling/iceberg_trajectories.nc'

else:
	Number_of_files=len(sys.argv)-1
	
#Imput parameters
berg_num=296190 #index of ic

#Flags and on/off switches
plot_timeseries=0
plot_trajectory=1
plot_bergsize=1
fix_horizontal_scale=0
plot_all_bergs=0  #plots all the icebergs and not just the one given by n
plot_initial_size=0
plot_size_at_start_and_finish_only=0
short_form=1 #only looks up the lat / lon

#Creating and filling up the bergs
for loop in xrange(0,Number_of_files):
	if len(sys.argv)==1:
		if loop==0: filename=filename1
		if loop==1: filename=filename2
	else:
		filename=sys.argv[loop+1]
	print filename	

	berg=[]; del berg
	bergs=bergs_open(filename)
	bergs=bergs_read(bergs,'width')
	bergs=bergs_read(bergs,'length')
	bergs=bergs_read(bergs,'lat')
	bergs=bergs_read(bergs,'lon')
	bergs=bergs_yearday(bergs)
	
	N=len(bergs.berg)
	if plot_all_bergs==0:
		N=1
	for n in xrange (0, N):
		if plot_all_bergs==0:
			n=berg_num

		#Loading variables
		lon=bergs.berg[n].lon[:]
		lat=bergs.berg[n].lat[:]
		width=bergs.berg[n].width[:]
		length=bergs.berg[n].length[:]
		yearday=bergs.berg[n].yearday[:]

		#Finding the index of the initial time
		initial_time_ind=sorted(list(enumerate(yearday)),key=lambda x: x[0])[-1][0]
		final_time_ind=sorted(list(enumerate(yearday)),key=lambda x: x[0])[0][0]

###########################   Begining plotting  #################################

		if plot_trajectory==1:
			#plt.subplot(2,1,1)
			if loop==0:
				if n==1:
					plt.plot(lon,lat,'ko-')
				if n!=1:
					plt.plot(lon,lat,'go-')
			else:
				plt.plot(lon,lat,'bo-')
			plt.xlabel('longitude (deg)')
			plt.ylabel('latitude (deg)')
			plt.title('Iceberg trajectory')
			plt.grid(True)
			#savefig("test.png")

			if plot_bergsize==1:
				Radius_earth=6378.135*1000;
				circ_ind=np.linspace(0,2*pi,100);
				for k in xrange (0 ,len(lat)):
    					L_eff=sqrt(((width[k]*length[k])/pi));
					plt.plot(lon[k],lat[k],'bo',linewidth=5)
				        #d_lat=(L_eff/Radius_earth)*(180/pi);
					#d_lon=L_eff/(111.320*(10**(3))*cos(lat[k]*pi/180));
				        d_lat=(L_eff/Radius_earth)*(180/pi);
					d_lon=(L_eff/(Radius_earth*cos(lat[k]*pi/180)))*(180/pi);
                                        if plot_size_at_start_and_finish_only==0:
					        plt.plot(lon[k]+(d_lon*cos(circ_ind)),lat[k]+(d_lat*sin(circ_ind)),'b');
					if k==initial_time_ind and plot_initial_size==1:
						plt.plot(lon[k]+(d_lon*cos(circ_ind)),lat[k]+(d_lat*sin(circ_ind)),'r');
					if k==final_time_ind:
						plt.plot(lon[k]+(d_lon*cos(circ_ind)),lat[k]+(d_lat*sin(circ_ind)),'r');
				
			#Setting axes to have equal distance
				if fix_horizontal_scale==1:
					y_scale=(max(lat)-min(lat))*(Radius_earth*pi/180)
					x_scale=(max(lon)-min(lon))*(111.320*(10**(3))*cos(np.mean(lat)/180*pi))
					d=0.6*max(x_scale,y_scale)
					d_lat=d*(180/(Radius_earth*pi))
					d_lon=d/(111.320*(10**(3))*cos(np.mean(lat)*pi/180))
					plt.xlim([((min(lon)+max(lon))/2)-d_lon,((min(lon)+max(lon))/2)+d_lon])
					plt.ylim([((min(lat)+max(lat))/2)-d_lat,((min(lat)+max(lat))/2)+d_lat])
	
			plt.plot(lon[initial_time_ind],lat[initial_time_ind],'ro')
plt.show()

if plot_timeseries==1:
	n=0
	#plt.subplot(2,1,2)
	plt.plot(yearday,lon,'ko-')
	plt.xlabel('yearday')
	plt.ylabel('lon (deg)')
	plt.grid(True)
	plt.savefig("test7.png")
	plt.show()


print 'Script complete'
