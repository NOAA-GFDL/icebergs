#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import numpy as np
import os
import numpy
import netCDF4 as nc
#import matplotlib
#import matplotlib.pyplot as plt
#from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys

#Clear screen
def main():

	os.system('clear')
	input_geometry_filename='input_files/Isomip_ocean_geometry.nc'

	#Flags
	new_filename='output_files/one_file.nc'

	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file


	with nc.Dataset(input_geometry_filename) as file:
		

		Depth = file.variables['D'][:,:]

	M= Depth.shape
	ny=M[0]
	nx=M[1]
	print nx,ny
	
	#Creating the topog file
	#g=Dataset(new_filename,'w') # w if for creating a file
	g=nc.Dataset(new_filename,'w',format='NETCDF3_CLASSIC')
	

	lat=g.createDimension('lat',ny)
	lon=g.createDimension('lon',nx)
	time=g.createDimension('time', None)

	lat=g.createVariable('lat','f4',('lat'))
	lon=g.createVariable('lon','f4',('lon'))
	time=g.createVariable('time','f4',('time'))
	one_var=g.createVariable('one_var','f4',('lat','lon'))

	time.long_name = "time" ;
	time.units = "days since 1900-01-01 00:00:00" ;
	time.cartesian_axis = "T" ;
	time.calendar_type = "NOLEAP" ;
	time.calendar = "NOLEAP" ;
	time.bounds = "time_bounds" ;
	time.modulo = " " ;

	lon.units = "degrees_E" ;
	lon.cartesian_axis = "X" ;
	#lon.edges = "xb" ;

	#lat.long_name = "latitude" ;
	lat.units = "degrees_N" ;
	lat.cartesian_axis = "Y" ;
	#lat.edges = "yb" ;

	one_var.long_name="unnamed_variable"
	one_var.units="m/s"
	

	
	g.variables['one_var'][:]=1.
	g.variables['lon'][:]=np.linspace(0.5,239.5,nx,endpoint=True)
	g.variables['lat'][:]=np.linspace(0.5,39.5,ny,endpoint=True)
	g.variables['time'][0]=1.

	g.sync()
	g.close()





	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
