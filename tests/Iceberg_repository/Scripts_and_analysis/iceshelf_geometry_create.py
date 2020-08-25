#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import numpy as np
import os
import numpy
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys

#Clear screen
def main():

	os.system('clear')
	input_geometry_filename='input_files/Isomip_ocean_geometry.nc'
	input_iceshelf_filename='input_files/Ocean1_3D_no_calving.nc'

	new_filename='output_files/Ocean1_3D_no_calving_trimmed.nc'


	#Parameters
	rho_water=1028
	rho_ice=918

	#f=Dataset(input_geometry_filename,'r')
	#g=Dataset(new_filename,'w') # w if for creating a file


	with nc.Dataset(input_geometry_filename) as file:
		Depth = file.variables['D'][:,:]
	with nc.Dataset(input_iceshelf_filename) as file:
		thick = file.variables['thick'][:,:]
	with nc.Dataset(input_iceshelf_filename) as file:
		area= file.variables['area'][:,:]



	M=thick.shape
	thick_orig=np.zeros((M[0],M[1]))
	for i in range(M[0]):
		for j in range(M[01]):
			thick_orig[i,j]=thick[i,j]
	
	#Making ice shelf symetric
	M=thick.shape
	thick_tmp=thick
	print M
	for i in range(M[0]/2):
		for j in range(M[1]):
			i_sym=(M[0]-1-i)
			if j==0:
				print i, i_sym
			thick_tmp[i,j]=((thick[i,j]+thick[i_sym,j])/2.)
			thick_tmp[i_sym,j]=((thick[i,j]+thick[i_sym,j])/2.)
	thick=thick_tmp

	#Making draft increase going inwards
	M=thick.shape
	thick_tmp=thick
	print M
	for k in range(5):
		#for i in range(M[0]/2):
		#for i in range(M[0]/8):
		for i in range(5):
			#for j in range(M[1]):
			for j in range(118):
				i_sym=(M[0]-1-i)
				if (thick_tmp[i,j]<thick_tmp[i+1,j]):
					thick_tmp[i,j]=thick_tmp[i+1,j]
					#thick_tmp[i,j]=((thick_tmp[i+1,j] +thick_tmp[i+2,j])/2)
				thick_tmp[i_sym,j]=thick_tmp[i,j]

		thick=thick_tmp
	



	draft=thick*(rho_water/rho_ice)

	ocean_thick=Depth-draft
	ocean_thick[np.where(ocean_thick<0)]=0.

	#ocean_thick=draft
	
	subplot(4,1,1)
	cNorm = mpl.colors.Normalize(vmin=0, vmax=10)
	cmap='jet'
	plt.pcolormesh(ocean_thick,norm=cNorm,cmap=cmap)
	plt.colorbar()
	
	subplot(4,1,2)
	cNorm = mpl.colors.Normalize(vmin=0, vmax=1000)
	cmap='jet'
	plt.pcolormesh(thick_orig,norm=cNorm,cmap=cmap)
	plt.colorbar()


	subplot(4,1,3)
	cNorm = mpl.colors.Normalize(vmin=0, vmax=1000)
	cmap='jet'
	plt.pcolormesh(thick,norm=cNorm,cmap=cmap)
	plt.colorbar()

	subplot(4,1,4)
	cNorm = mpl.colors.Normalize(vmin=0, vmax=10)
	cmap='jet'
	plt.pcolormesh(thick-thick_orig,norm=cNorm,cmap=cmap)
	plt.colorbar()



	ny=M[0]
	nx=M[1]
	print nx,ny


	#Creating the topog file
	g=Dataset(new_filename,'w') # w if for creating a file

	yt=g.createDimension('yt',ny)
	xt=g.createDimension('xt',nx)

	thick_h=g.createVariable('thick','f4',('yt','xt'))
	g.variables['thick'][:]=thick
	area_h=g.createVariable('area','f4',('yt','xt'))
	g.variables['area'][:]=area

	print 'Writing file: ' , new_filename
	
	g.sync()
	g.close()



	plt.show()


	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
