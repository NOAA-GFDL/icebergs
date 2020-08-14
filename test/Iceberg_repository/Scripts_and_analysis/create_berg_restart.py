#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import numpy as np
import os
#import matplotlib
#import matplotlib.pyplot as plt
#from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys

#Clear screen
os.system('clear')

new_filename = 'output_files/icebergs_test.res.nc'

f = Dataset(new_filename,'w', format='NETCDF3_CLASSIC')


#Creating the calving file
#f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')

if True:
	if True:

		i=f.createDimension('i', None)
		
		i=f.createVariable('i','i')
		#i=f.createVariable('i','i',('i'))

		lon=f.createVariable('lon','d',('i'))
		lon.long_name = "longitude" ;
		lon.units = "degrees_E" ;
		lon.checksum = "               0" ;

		lat=f.createVariable('lat','d',('i'))
		lat.long_name = "latitude" ;
		lat.units = "degrees_N" ;
		lat.checksum = "               0" ;

		uvel=f.createVariable('uvel','d',('i'))
		uvel.long_name = "zonal velocity" ;
		uvel.units = "m/s" ;
		uvel.checksum = "               0" ;

		vvel=f.createVariable('vvel','d',('i'))
		vvel.long_name = "meridional velocity" ;
		vvel.units = "m/s" ;
		vvel.checksum = "               0" ;

		mass=f.createVariable('mass','d',('i'))
		mass.long_name = "mass" ;
		mass.units = "kg" ;
		mass.checksum = "               0" ;

		axn=f.createVariable('axn','d',('i'))
		axn.long_name = "explicit zonal acceleration" ;
		axn.units = "m/s^2" ;
		axn.checksum = "               0" ;

		ayn=f.createVariable('ayn','d',('i'))
		ayn.long_name = "explicit meridional acceleration" ;
		ayn.units = "m/s^2" ;
		ayn.checksum = "               0" ;

		bxn=f.createVariable('bxn','d',('i'))
		bxn.long_name = "inplicit zonal acceleration" ;
		bxn.units = "m/s^2" ;
		bxn.checksum = "               0" ;

		byn=f.createVariable('byn','d',('i'))
		byn.long_name = "implicit meridional acceleration" ;
		byn.units = "m/s^2" ;
		byn.checksum = "               0" ;

		ine=f.createVariable('ine','i',('i'))
		ine.long_name = "i index" ;
		ine.units = "none" ;
		ine.packing = 0 ;
		ine.checksum = "               0" ;

		jne=f.createVariable('jne','i',('i'))
		jne.long_name = "j index" ;
		jne.units = "none" ;
		jne.packing = 0 ;

		thickness=f.createVariable('thickness','d',('i'))
		thickness.long_name = "thickness" ;
		thickness.units = "m" ;
		thickness.checksum = "               0" ;

		width=f.createVariable('width','d',('i'))
		width.long_name = "width" ;
		width.units = "m" ;
		width.checksum = "               0" ;

		length=f.createVariable('length','d',('i'))
		length.long_name = "length" ;
		length.units = "m" ;
		length.checksum = "               0" ;

		start_lon=f.createVariable('start_lon','d',('i'))
		start_lon.long_name = "longitude of calving location" ;
		start_lon.units = "degrees_E" ;
		start_lon.checksum = "               0" ;

		start_lat=f.createVariable('start_lat','d',('i'))
		start_lat.long_name = "latitude of calving location" ;
		start_lat.units = "degrees_N" ;
		start_lat.checksum = "               0" ;

		start_year=f.createVariable('start_year','i',('i'))
		start_year.long_name = "calendar year of calving event" ;
		start_year.units = "years" ;
		start_year.packing = 0 ;
		start_year.checksum = "               0" ;

		iceberg_num=f.createVariable('iceberg_num','i',('i'))
		iceberg_num.long_name = "identification of the iceberg" ;
		iceberg_num.units = "dimensionless" ;
		iceberg_num.packing = 0 ;
		iceberg_num.checksum = "               0" ;

		start_day=f.createVariable('start_day','d',('i'))
		start_day.long_name = "year day of calving event" ;
		start_day.units = "days" ;
		start_day.checksum = "               0" ;

		start_mass=f.createVariable('start_mass','d',('i'))
		start_mass.long_name = "initial mass of calving berg" ;
		start_mass.units = "kg" ;
		start_mass.checksum = "               0" ;

		mass_scaling=f.createVariable('mass_scaling','d',('i'))
		mass_scaling.long_name = "scaling factor for mass of calving berg" ;
		mass_scaling.units = "none" ;
		mass_scaling.checksum = "               0" ;

		mass_of_bits=f.createVariable('mass_of_bits','d',('i'))
		mass_of_bits.long_name = "mass of bergy bits" ;
		mass_of_bits.units = "kg" ;
		mass_of_bits.checksum = "               0" ;

		heat_density=f.createVariable('heat_density','d',('i'))
		heat_density.long_name = "heat density" ;
		heat_density.units = "J/kg" ;
		heat_density.checksum = "               0" ;

		halo_berg=f.createVariable('halo_berg','d',('i'))
		halo_berg.long_name = "halo_berg" ;
		halo_berg.units = "dimensionless" ;
		halo_berg.checksum = "               0" ;

		static_berg=f.createVariable('static_berg','d',('i'))
		static_berg.long_name = "static_berg" ;
		static_berg.units = "dimensionless" ;
		static_berg.checksum = "               0" ;

print 'Created restart file: ' ,new_filename
print 'Script complete'

f.sync()
f.close()





print 'Script complete'
