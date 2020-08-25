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
	use_flat_bottom=True
	use_side_walls=True

	#Parameters
	constant_depth=750.

	if use_flat_bottom==True:
		new_filename='output_files/topog_flat'
	else:
		new_filename='output_files/topog'
	if use_side_walls==True:
		new_filename=new_filename+'_boundaries'
	new_filename=new_filename+'.nc'


	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file


	with nc.Dataset(input_geometry_filename) as file:
		Depth = file.variables['D'][:,:]

	M= Depth.shape
	ny=M[0]
	nx=M[1]
	print nx,ny

	if use_flat_bottom==True:
		Depth[:,:]=(Depth[:,:]*0.0)+constant_depth
	if use_side_walls==True:
		Depth[0,:]=0.0
		Depth[ny-1,:]=0.0
		Depth[:,0]=0.0
		Depth[:,nx-1]=0.0
	
		Depth[ny-2,:]=0.0
		Depth[ny-3,:]=0.0
		Depth[ny-4,:]=0.0
		Depth[ny-5,:]=0.0
		Depth[ny-6,:]=0.0

	#Creating the topog file
	g=Dataset(new_filename,'w') # w if for creating a file

	yt=g.createDimension('yt',ny)
	xt=g.createDimension('xt',nx)

	depth_h=g.createVariable('depth','f4',('yt','xt'))
	g.variables['depth'][:]=Depth

	print 'Writing file: ' , new_filename
	
	g.sync()
	g.close()





	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
