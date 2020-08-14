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
	input_filename='input_files/Isomip_ice_geometry.nc'
	new_filename='output_files/ISOMIP_MIB.nc'

	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file

	rho_ice=850.

	with nc.Dataset(input_filename) as file:
		ocean_mask = file.variables['openOceanMask'][:,:]
		upperSurface = file.variables['upperSurface'][:,:]
		lowerSurface = file.variables['lowerSurface'][:,:]
		h_ice=upperSurface-lowerSurface
		mass=h_ice*rho_ice
		x = file.variables['x'][:]
		y = file.variables['y'][:]


	M= mass.shape
	ny=M[0]
	nx=M[1]
	print nx,ny
	print mass.shape

	#subsampling data onto a grid half the size
	MIB=np.zeros((ny/2,nx/2))
	print MIB.shape
	for i in range(0,nx,2):
		for j in range(0,ny,2):
			MIB[j/2,i/2]=mass[j,i]

	
	#Creating the topog file
	g=Dataset(new_filename,'w') # w if for creating a file

	time=g.createDimension('time',1)
	yt=g.createDimension('yt',ny)
	xt=g.createDimension('xt',nx)

	mass_h=g.createVariable('mass','f4',('time','yt','xt'))
	g.variables['mass'][:]=mass

	#Creating subsampled version
	yh=g.createDimension('yh',ny/2)
	xh=g.createDimension('xh',nx/2)
	MIB_h=g.createVariable('MIB','f4',('time','yh','xh'))
	g.variables['MIB'][:]=MIB


	g.sync()
	g.close()





	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
