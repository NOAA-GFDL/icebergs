#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("GTKAgg")
from pylab import *
#import matplotlib.pyplot as plt
import pdb
import netCDF4 as nc
from hexagon_area import Divide_hexagon_into_4_quadrants_old
from hexagon_area import Hexagon_into_quadrants_using_triangles


def load_ISOMIP_ice_geometry(filename):
	with nc.Dataset(filename) as file:
		ocean_mask = file.variables['openOceanMask'][:,:]
		upperSurface = file.variables['upperSurface'][:,:]
		lowerSurface = file.variables['lowerSurface'][:,:]
		x = file.variables['x'][:]
		y = file.variables['y'][:]

	ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
	h_ice=upperSurface-lowerSurface #The ice thickness
	M=ice_mask.shape

	return [x,y,ice_mask,h_ice]

def Switch_x_and_y_directions(lat,lon, A, B, C):
        print 'Switching x and y directions!!!'
        lat_new=[]
        lon_new=[]
        for i in range(len(lat)):
                lat_new.append(lon[i])
                lon_new.append(lat[i])
        M=A.shape
        A_new=np.zeros([M[1],M[0]])
        B_new=np.zeros([M[1],M[0]])
        C_new=np.zeros([M[1],M[0]])
        for i in range(M[0]):
                for j in range(M[1]):
                        A_new[j,i]=A[i,j]
                        B_new[j,i]=B[i,j]
                        C_new[j,i]=C[i,j]

        return  [lat_new,lon_new,A_new,B_new, C_new]


def create_clipped_icethickness_file(h_ice,area,mass,grid_area,gravity,New_ice_thickness_filename):
	#Creating clipped file
        [ny, nx]= h_ice.shape ; 
	
        g=Dataset(New_ice_thickness_filename,'w') # w if for creating a file

        g.createDimension('nx',nx)
        g.createDimension('ny',ny)

        thick_h=g.createVariable('thick','f8',('ny','nx'))
        area_h=g.createVariable('area','f8',('ny','nx'))
        p_surf_h=g.createVariable('p_surf','f8',('ny','nx'))

        thick_h.units = 'm'
        thick_h.standard_name = 'ice shelf thickness (clipped)'
        area_h.units = 'm2'
        area_h.standard_name = 'ice shelf area'
        p_surf_h.units = 'Pa'
        p_surf_h.standard_name = 'surface pressure due to ice shelf'

	p_surf=(gravity*mass)/grid_area
        g.variables['thick'][:]=h_ice
        g.variables['area'][:]=area
        g.variables['p_surf'][:]=p_surf
        print 'Creating clipped ice file: ' , New_ice_thickness_filename
        
        g.sync()
        g.close()


def load_ISOMIP_reduced_ice_geometry(ice_filename,topog_filename):
	with nc.Dataset(ice_filename) as file:
		h_ice = file.variables['thick'][:,:]
		area =  file.variables['area'][:,:]
		M=h_ice.shape

	y=np.linspace(1000,79000,M[0],endpoint=True)
	#x=np.linspace(321000,799000,M[1],endpoint=True)
	x=np.linspace(1000,479000,M[1],endpoint=True)
	ice_mask=h_ice>0.
	
	return [x,y,ice_mask,h_ice]

####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#Flags
	Switch_x_and_y_to_rotate_90=False


	ISOMIP_ice_geometry_filename='input_files/Isomip_ice_geometry.nc'
	ISOMIP_reduced_ice_geometry_filename='input_files/isomip_ice_shelf1.nc'
	ISOMIP_reduced_ice_geometry_filename='output_files/isomip_ice_shelf1_clipped.nc'
	Weddell_ice_geometry_filename='input_files/Bedmap2_gridded_subset_Weddell_Sea_Region.nc'
	ISOMIP_topography_filename='input_files/Isomip_topog.nc'
	Output_filename='output_files/isomip_ice_shelf1_new.nc'

	(x,y,ice_mask,h_ice)=load_ISOMIP_reduced_ice_geometry(ISOMIP_reduced_ice_geometry_filename,ISOMIP_topography_filename)

	#Parameters
	rho_ice=850.
	gravity=9.8


	grid_area=(x[1]-x[0])*(x[1]-x[0])
	New_area=h_ice*0.+grid_area
	#New_area=ice_mask*grid_area
	z_extra=0.00000000001
	h_ice=h_ice+z_extra
	#h_ice=(h_ice>0.)*0.00000000000001
	New_mass=New_area*h_ice*rho_ice
	h_ice=(0*h_ice)+z_extra
	if Switch_x_and_y_to_rotate_90 is True:
		lat=[];lon=[]
		[lat,lon, h_ice,New_area,New_mass]=Switch_x_and_y_directions(lat,lon, h_ice,New_area,New_mass)
		x_temp=x
		x=y
		y=x_temp
	create_clipped_icethickness_file(h_ice,New_area,New_mass,grid_area,gravity,Output_filename)
	


	#Plotting the data
	plot_data=(h_ice)
        vmax=np.max(abs(plot_data))
        cNorm = mpl.colors.Normalize(vmin=0., vmax=vmax)
        plt.pcolormesh(x,y,plot_data,cmap='bwr',norm=cNorm)
        plt.title('Ice geometry')
	plt.colorbar()
	

	plt.show()
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














