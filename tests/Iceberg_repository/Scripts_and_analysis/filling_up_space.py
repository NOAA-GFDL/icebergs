#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("GTKAgg")
from pylab import *
from scipy import interpolate
#import matplotlib.pyplot as plt
import argparse
import pdb
import netCDF4 as nc
from hexagon_area import Divide_hexagon_into_4_quadrants_old
from hexagon_area import Hexagon_into_quadrants_using_triangles



def spread_mass_to_ocean(Nx,Ny,i,j,mass_on_ocean,x,y,Area,Mass,element_type,grid_area,static_berg,lat1, lat2, lat3, lat4, lon1, lon2, lon3, l):
		#Using physical space

		#Initialize weights for each cell	
		yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0. 
		yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

		if element_type=='square':
		#if True:
			L = min(( (np.sqrt(Area)  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim length of the square)
			xL=min(0.5, max(0., 0.5-(x/L)))
			xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
			xC=max(0., 1.-(xL+xR))
			yD=min(0.5, max(0., 0.5-(y/L)))
			yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
			yC=max(0., 1.-(yD+yU))

			yDxL=yD*xL#*grd%msk[i-1,j-1]
			yDxC=yD*xC#*grd%msk[i  ,j-1]
			yDxR=yD*xR#*grd%msk[i+1,j-1]
			yCxL=yC*xL#*grd%msk[i-1,j  ]
			yCxR=yC*xR#*grd%msk[i+1,j  ]
			yUxL=yU*xL#*grd%msk[i-1,j+1]
			yUxC=yU*xC#*grd%msk[i  ,j+1]
			yUxR=yU*xR#*grd%msk(i+1,j+1]
			yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

		if element_type=='hexagon':
		#if False:#element_type=='hexagon':
			H = min(( (np.sqrt(Area/(2*sqrt(3))) ;  #No longer no dimensionalizing, but rather using actual distances
			S=(2/np.sqrt(3))*H #Side of the hexagon
			grid_dist=1e15 #This should be changed to the minimum grid spacinge.
			if S>0.5*(grid_dist):
				print 'Elements must be smaller than a whole gridcell', 'i.e.: S= ' , S , '>=0.5',i,j
				halt

			#Subtracting the position of the nearest corner from x,y
			#[origin_x, origin_y]=find_closest_corner(x,y,
			origin_x=1 ; origin_y=1
			if x<0.5:
				origin_x=0
			if y<0.5:
				origin_y=0

			x0=(x-origin_x) #Position of the hexagon center, relative to origin at the nearest vertex
			y0=(y-origin_y)
			
			#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
			(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,0.)
			if min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <0:
				print 'Yolo'
				print x0,y0,H
				#print min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4))
				print Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4

			Area_Q1=Area_Q1/Area_hex
			Area_Q2=Area_Q2/Area_hex
			Area_Q3=Area_Q3/Area_hex
			Area_Q4=Area_Q4/Area_hex

			#Now, you decide which quadrant belongs to which mass on ocean cell.
			#Top right vertex
			if x>=0.5 and y>= 0.5:
				yUxR=Area_Q1
				yUxC=Area_Q2
				yCxC=Area_Q3
				yCxR=Area_Q4
			
			#Top left vertex
			if x<0.5 and y>= 0.5:
				yUxC=Area_Q1
				yUxL=Area_Q2
				yCxL=Area_Q3
				yCxC=Area_Q4
			
			#Bottom left vertex
			if x<0.5 and y< 0.5:
				yCxC=Area_Q1
				yCxL=Area_Q2
				yDxL=Area_Q3
				yDxC=Area_Q4

			#Bottom right vertex
			if x>=0.5 and y< 0.5:
				yCxR=Area_Q1
				yCxC=Area_Q2
				yDxC=Area_Q3
				yDxR=Area_Q4

			#if Sector<4 and Sector>-1:
			#	print Sector
				
					
			#Check that this is true
			if abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>0.001:
				print 'All the mass is not being used!!!'
				#print W1 , W2 , W3 , W4 , W5 , W6
				#print Area_Upper, Area_Lower, Area_right, Area_left
				print 'Areas: ',Area_hex,Area_hex*Area_Q1, Area_hex*Area_Q2, Area_hex*Area_Q3, Area_hex*Area_Q4
				print 'x0=',x0, 'y0=',y0, 'H=', H
				print 'Total area= ',(Area_Q1+Area_Q2+Area_Q3+Area_Q4)#, Sector
				

		#Accounting for masked points 
		a=1. ; b=1. ; c=1. ; d=1.
		if i==0:
			a=0.;
		if j==0: 
			b=0.;
		if i==Nx-1:
			c=0;
		if j==Ny-1:
			d=0;

		fraction_used= ((yDxL*a*b) + (yDxC*b) + (yDxR*b*c) +(yCxL*a) + (yCxR*c) + (yUxL*a*d) + (yUxC*d) + (yUxR*c*d) + (yCxC))
		if static_berg>0.5:
			fraction_used=1.

		mass_on_ocean[i,j,1]=mass_on_ocean[i,j,1]+(a*b*yDxL*Mass/fraction_used)
		mass_on_ocean[i,j,2]=mass_on_ocean[i,j,2]+(b*yDxC*Mass/fraction_used)
		mass_on_ocean[i,j,3]=mass_on_ocean[i,j,3]+(b*c*yDxR*Mass/fraction_used)
		mass_on_ocean[i,j,4]=mass_on_ocean[i,j,4]+(a*yCxL*Mass/fraction_used)
		mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+(yCxC*Mass/fraction_used)
		mass_on_ocean[i,j,6]=mass_on_ocean[i,j,6]+(c*yCxR*Mass/fraction_used)
		mass_on_ocean[i,j,7]=mass_on_ocean[i,j,7]+(a*d*yUxL*Mass/fraction_used)
		mass_on_ocean[i,j,8]=mass_on_ocean[i,j,8]+(d*yUxC*Mass/fraction_used)
		mass_on_ocean[i,j,9]=mass_on_ocean[i,j,9]+(c*d*yUxR*Mass/fraction_used)


		return mass_on_ocean


def Create_icebergs(Radius, R_earth):
	dx_berg=[]
	dy_berg=[]
	lat=[]
	lon=[]
	berg_count=0
	
	#Flags
	adjust_lat_ref=True


	#Developing some elements in cartesian coordinates:
	lat_init=20.0 ; lon_init=0.0 #Starting point


	for i in range(40):
		x_start=Radius
		y_start=(Radius)
		x_val=x_start+(2*i*Radius)
		for j in range(40):
			y_val=y_start+(2*j*Radius)  ;# x_val=x_start + (np.sqrt(3)*Radius*i)

			berg_count=berg_count+1
			dx_berg.append(x_val)
			dy_berg.append(y_val)

	Number_of_bergs=berg_count

	#Convert_to_lat_lon==True:
	dlat_dy=(180/np.pi)*(1/R_earth)
	#dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
	#dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
	for i in range(Number_of_bergs):
		#Finding latittude
		lat.append(lat_init+(dlat_dy*dy_berg[i]))
		#Finding longitude
		if adjust_lat_ref==True:
			lat_ref=lat_init+(dlat_dy*dy_berg[i])
			#print lat_ref
		dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
		lon.append(lon_init+(dlon_dx*dx_berg[i] ))

	return [Number_of_bergs, lat, lon] 
	
def Create_lat_lon_grid(lat_min,lat_max,lon_min,lon_max, N_lat, N_lon):
	grd_lat=np.linspace(lat_min,lat_max,N_lat)
	grd_lon=np.linspace(lon_min,lon_max,N_lon)
	field=np.zeros((N_lat,N_lon))

	return [grd_lat, grd_lon, field]
	
def Projecting_to_local_cartesian(lat0, lon0, lat_corner, lon_corner, R_earth):
	#Creates a local tangent plan around an iceberg at lat0, lon0
	d_lat=(lat_corner-lat0)
	d_lon=(lon_corner-lon0)

	dy=R_earth*(np.pi/180.0)*d_lat
	dx=R_earth*(np.pi/180.0)*cos(lat0*(np.pi/180))*d_lon
	#dx=R_earth*(np.pi/180.0)*cos(lat_corner*(np.pi/180))*d_lon

	return dy,dx


def find_corners_of_cell(lat_berg,lon_berg,grd_lat,grd_lon):
	j=np.sum(grd_lat<lat_berg)
	i=np.sum(grd_lon<lon_berg)
	print lat_berg,lon_berg
	print grd_lat[j-1],  grd_lat[j]
	print grd_lon[i-1],  grd_lon[i]
	i1=i-1 ; j1=j-1
	i2=i ; j2=j-1
	i3=i   ; j3=j
	i4=i-1   ; j4=j

	return [j1,j2,j3,j4,i1,i2,i3,i4]

def plot_grid_cell_and_point(lat_berg, lon_berg, lat1, lat2, lat3, lat4, lon1, lon2, lon3, lon4):
	plt.plot(lon_berg,lat_berg,'o')
	#plt.plot(np.array([lon1, lon2]), np.array([lat1, lat2]))
	#plt.plot(np.array([lon2, lon3]), np.array([lat2, lat3]))
	plt.plot(np.array([lon1, lon2, lon3, lon4, lon1]), np.array([lat1, lat2, lat3, lat4, lat1]))
	plt.plot(lon1,lat1,'*r')
	plt.plot(lon2,lat2,'*c')
	plt.plot(lon4,lat3,'*g')
	plt.plot(lon4,lat4,'*b')
	eps_lat=(lat3-lat2)/2.
	eps_lon=(lon2-lon1)/2.
	#print 'lon:',  lon1,lon2, lon3, lon4
	#print 'lat:', lat1,lat2, lat3, lat4
	print eps_lat, eps_lon
	plt.ylim([lat2-eps_lat, lat3+eps_lat])
	plt.xlim([lon1-eps_lon, lon2+eps_lon])


####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#Parameters
	R_earth=6360000.0
	Radius=5000.0 #In meters

	#Lat Lon Parameters
	lat_min=20.0 ;	lat_max=50.0 ;	N_lat=5
	lon_min=0.0  ;	lon_max=40.0 ;	N_lon=5
	
	
	



	#Defining Lat Lon
	(grd_lat, grd_lon, field)= Create_lat_lon_grid(lat_min,lat_max,lon_min,lon_max, N_lat, N_lon)

	#Creating Icebergs
	#(Number_of_bergs, lat_berg, lon_berg)=Create_icebergs(Radius, R_earth)



	lat_berg=32.0
	lon_berg=22.2

	#Finding grid index of cell of berg:
	j1,j2,j3,j4,i1,i2,i3,i4	=find_corners_of_cell(lat_berg,lon_berg,grd_lat,grd_lon)
	lat1=grd_lat[j1] ; lat2=grd_lat[j2]; lat3=grd_lat[j3]; lat4=grd_lat[j4]
	lon1=grd_lon[i1] ; lon2=grd_lon[i2]; lon3=grd_lon[i3]; lon4=grd_lon[i4]
	i_berg=i1 ;j_berg=j1 #Lower left point


	#Converting points to tangent plane
	(y_berg,x_berg)= Projecting_to_local_cartesian(lat_berg, lon_berg, lat_berg, lon_berg, R_earth)
	(y1,x1)= Projecting_to_local_cartesian(lat_berg, lon_berg, lat1, lon1, R_earth)
	(y2,x2)= Projecting_to_local_cartesian(lat_berg, lon_berg, lat2, lon2, R_earth)
	(y3,x3)= Projecting_to_local_cartesian(lat_berg, lon_berg, lat3, lon3, R_earth)
	(y4,x4)= Projecting_to_local_cartesian(lat_berg, lon_berg, lat4, lon4, R_earth)

	#Normalization of the shape in tangent plane
	#x_cell=(x_berg-x[i_val])/sqrt(grid_area)+0.5



	#Plotting everything
	plt.figure(figsize=(16,8))
	plt.subplot(1,2,1)
	plot_grid_cell_and_point(lat_berg,lon_berg, lat1, lat2, lat3, lat4, lon1, lon2, lon3, lon4)
	plt.xlabel('lon') ; plt.ylabel('lat') ; 
	plt.subplot(1,2,2)
	plot_grid_cell_and_point(y_berg,x_berg, y1, y2, y3, y4, x1, x2, x3, x4)
	plt.xlabel('x') ; plt.ylabel('y') ; 
	print x1,x2,x3,x4
	print y1,y2,y3,y4
		





	#Plotting berg positions
	#plt.pcolor(grd_lon,grd_lat,field)
	#plt.colorbar
	#plt.scatter(lon_berg, lat_berg,c='red')

	plt.show()
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














