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
	input_filename='input_files/isomip_ice_shelf1.nc'
	new_filename='output_files/flat_isomip_ice_shelf1.nc'

	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file

	rho_ice=850.

	with nc.Dataset(input_filename) as file:
		thick = file.variables['thick'][:,:]
		area = file.variables['area'][:,:]
		height = file.variables['height'][:,:]

	M= thick.shape
	ny=M[0]
	nx=M[1]
	print nx,ny
	print thick.shape

	#Setting thickness to a prescribed value
	prescribed_thickness=1.
	thick[np.where(thick>0)]=prescribed_thickness

	
	#Creating the file
	g=Dataset(new_filename,'w') # w if for creating a file

	yt=g.createDimension('yt',ny)
	xt=g.createDimension('xt',nx)

	thick_h=g.createVariable('thick','f4',('yt','xt'))
	g.variables['thick'][:]=thick>0
	thick_h=g.createVariable('area','f4',('yt','xt'))
	g.variables['area'][:]=area
	height_h=g.createVariable('height','f4',('yt','xt'))
	g.variables['height'][:]=height

	g.sync()
	g.close()





	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
