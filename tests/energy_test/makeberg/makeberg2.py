#!/usr/bin/env python
import numpy as np
import math
from netCDF4 import Dataset
from pylab import *
#import pdb
import netCDF4 as nc

#
# Initialize 2 iceberg elements, which can later be bonded in the fortran code if needed
#

def Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg,uvel,vvel):

	print 'Writing iceberg restart files, with ' , Number_of_bergs  , 'icebergs..'
	# To copy the global attributes of the netCDF file

	#Input and output files
	#Create Empty restart file. This is later read so that the attributes can be used.
	Empty_restart_filename='output_files/Empty_icebergs.res.nc'
	create_empty_iceberg_restart_file(Empty_restart_filename)
	#Empty_restart_filename='input_files/icebergs.res.nc'

	#Read empty restart file
	f=Dataset(Empty_restart_filename,'r') # r is for read only
	#Write a new restart file
	g=Dataset('output_files/' + Ice_geometry_source + '_icebergs.res.nc','w', format='NETCDF3_CLASSIC') # w if for creating a file

	for attname in f.ncattrs():
		    setattr(g,attname,getattr(f,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in f.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		g.createDimension(dimname,Number_of_bergs)

	# To copy the variables of the netCDF file

	for varname,ncvar in f.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		#var[:] = ncvar[0]  #I commented out this line because it was causing errors. I'm not sure if it is needed.

		if varname=='i':
			var[:]=Number_of_bergs

		if varname=='iceberg_num':
			for j in range(Number_of_bergs):
				#var[j]=j+1
				var[j]=iceberg_num[j]

		if varname=='uvel_old' or varname=='vvel_old' or varname=='axn' or varname=='ayn'\
		or varname=='bxn' or varname=='byn' or  varname=='halo_berg' or varname=='heat_density' or varname=='lon_old' or varname=='lat_old' \
		or varname=='mass_of_bits' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' \
		or varname=='start_lat' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' or varname=='lat_old':\
			var[:]=0

                if varname=='uvel':
                        var[:]=uvel

                if varname=='vvel':
                        var[:]=vvel

		if varname=='mass_scaling':
			var[:]=mass_scaling

		if varname=='thickness':
			for j in range(Number_of_bergs):
				var[j]=thickness[j]

		if varname=='mass':
			for j in range(Number_of_bergs):
				var[j]=mass[j]

		if varname=='width'  or varname=='length':
			for j in range(Number_of_bergs):
				var[j]=width[j]

		if varname=='lon':
			for j in range(Number_of_bergs):
				var[j]=lon[j]

		if varname=='lat':
			for j in range(Number_of_bergs):
				var[j]=lat[j]

		if varname=='static_berg':
			for j in range(Number_of_bergs):
				var[j]=static_berg[j]


	f.close()
	g.close()


def create_empty_iceberg_restart_file(Empty_restart_filename):

	f = Dataset(Empty_restart_filename,'w', format='NETCDF3_CLASSIC')

	i=f.createDimension('i', None)
	lon=f.createVariable('i','i')

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
	jne.checksum = "               0" ;

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

	f.sync()
	f.close()



#-----------------------------------#
#               Main                #
#-----------------------------------#

grdres=1000.0
xc=4500.0
yc=15001.0
R_frac=0.45
radius=(np.sqrt(3)/2.)*(R_frac*grdres) #S is < 0.5 grid res
element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*radius**2)
width=np.sqrt(element_area)
thickness=300.0
rhoi=850.0
mass=thickness*rhoi*element_area

berg_count=2
berg_x=[xc-radius,xc+radius]
berg_y=[yc,yc]
berg_h=[thickness,thickness]
berg_width=[width,width]
berg_mass=[mass,mass]
berg_mass_scaling=[1,1]
berg_id=[1,2]
berg_static=[0,0]
berg_uvel=[0,0.2]
berg_vvel=[0,0]

print('Number of bergs',berg_count)

#Create iceberg restart file
Ice_geometry_source='Generic'
Create_iceberg_restart_file(berg_count,berg_x,berg_y,berg_h,berg_width,berg_mass,berg_mass_scaling,berg_id,Ice_geometry_source,berg_static,berg_uvel,berg_vvel)
