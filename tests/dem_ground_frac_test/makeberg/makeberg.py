#!/usr/bin/env python
import numpy as np
import math
from netCDF4 import Dataset
from pylab import *
#import pdb
import netCDF4 as nc

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

grdxmin=0; grdxmax=225000e3
grdymin=0; grdymax=225000e3
grdres=5000.0
R_frac=0.45
radius=1.5e3
thickness1=200.0
thickness2=200.0
rho_ice=850.0

nbergs=1 #number of conglomerates

#center coords of each (rectangular) conglomerate berg (CB=conglomerate berg)
CBxc=np.array([50000]); CByc=np.array([50000])
#side lengths
CBxl=np.array([15000]); CByl=np.array([35000])


CByl=CByl.astype(int); CBxl=CBxl.astype(int)

#--- xmax, xmin, ymax, ymin for each CB --
CBxmin=CBxc-(0.5*CBxl); CBxmax=CBxc+(0.5*CBxl)
CBymin=CByc-(0.5*CByl); CBymax=CByc+(0.5*CByl)
#this shouldn't happen:
CBxmin[CBxmin<grdxmin]=grdxmin; CBxmax[CBxmax>grdxmax]=grdxmax
CBymin[CBxmin<grdymin]=grdymin; CBymax[CBymax>grdymax]=grdymax
CBxc=0.5*(CBxmin+CBxmax); CByc=0.5*(CBymin+CBymax)

#--- min and max thicknesses for each CB ---
h1=thickness1
h2=thickness2
CBhmax=np.array([h1]); CBhmin=np.array([h2])

#radii
CBrad=np.array([radius])
print('radii',CBrad)

berg_x=[]; berg_y=[]
berg_id=[]; berg_static=[]
berg_width=[]; berg_bonds=[]
berg_h=[]; berg_mass_scaling=[]
berg_mass=[]
berg_uvel=[]; berg_vvel=[]

berg_count=0

for i in range(nbergs):
        x_start=CBxmin[i]+(CBrad[i]*2./np.sqrt(3))

        if x_start>CBxmax[i]:
                x_start=CBxmax[i]
        y_start0=CBymin[i]+CBrad[i]
        if y_start0>CBymax[i]:
                y_start0=CBymax[i]
        element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*(CBrad[i])**2)
        #H (thickness) is a linear function of a berg element's position from the
        #center of the CB. H=Hmax at the center of the CB and H=Hmin at a corner of the
        #CB, where the dist of a corner from the center is:
        cdistb=np.sqrt((CBxmin[i]-CBxc[i])**2+(CBymin[i]-CByc[i])**2)
        j=0
        x_val=x_start
        offset=0.0

        uvel=0.1
        vvel=0.0

        while x_val<=CBxmax[i] and x_val>=CBxmin[i]:
                y_start=y_start0+((j%2)*CBrad[i])+offset
                k=0
                y_val=y_start
                while y_val<=(CBymax[i]+offset):
                        berg_count=berg_count+1
                        berg_id.append(berg_count)
                        berg_x.append(x_val)
                        berg_y.append(y_val)
                        berg_width.append(sqrt(element_area))
                        #dist of berg elem from center of CB
                        bdistc=np.sqrt((x_val-CBxc[i])**2+(y_val-CByc[i])**2)
                        bh=CBhmin[i]*bdistc/cdistb + CBhmax[i]*(1-bdistc/cdistb)

                        berg_h.append(bh) #thickness
                        berg_mass_scaling.append(1)
                        berg_mass.append(bh*rho_ice*element_area)
                        berg_static.append(0)
                        berg_uvel.append(uvel)
                        berg_vvel.append(vvel)

                        k=k+1
                        y_val=y_start+(2*k*CBrad[i])
                j=j+1

                x_val=x_start+(np.sqrt(3)*CBrad[i]*j)

print('Number of bergs',berg_count)

#Create iceberg restart file
Ice_geometry_source='Generic'
Create_iceberg_restart_file(berg_count,berg_x,berg_y,berg_h,berg_width,berg_mass,berg_mass_scaling,berg_id,Ice_geometry_source,berg_static,berg_uvel,berg_vvel)
