#!/usr/bin/env python

from netCDF4 import Dataset  
import numpy as np 
import math
#import os
from pylab import *
#import pdb
import netCDF4 as nc

nx = 20
ny = 20
New_ice_thickness_filename='testberg.nc'
grid_area = 10.e3*10.e3 #10 km grid res
h_ice = 300.0

New_ice_thickness_filename='testtabularberg.nc'
g=Dataset(New_ice_thickness_filename,'w') # w if for creating a file

g.createDimension('nx',nx)
g.createDimension('ny',ny)
thick_h=g.createVariable('thick','f8',('ny','nx'))
area_h=g.createVariable('area','f8',('ny','nx'))
lat_h=g.createVariable('lat','f8','ny')
lon_h=g.createVariable('lon','f8','nx')
#lat_h=g.createVariable('lat','f8',('ny','nx'))
#lon_h=g.createVariable('lon','f8',('ny','nx'))

thick_h.units = 'm'
thick_h.standard_name = 'ice shelf thickness'
area_h.units = 'm2'
area_h.standard_name = 'ice shelf area'
lon_h.unit = 'm'
lon_h.standard_name = 'longitude'
lat_h.unit = 'm'
lat_h.standard_name = 'latitude'



#g.variables['thick'][:]=h_ice
g.variables['area'][:,:]=grid_area

#create a circular tabular berg centered on the following coords (km)
cx = 60.e3
cy = 60.e3

for i in range(nx):

        tx = float(i) * 10.e3
        g.variables['lon'][i] = tx
        
        for j in range(ny):
               
                ty = float(j) * 10.e3
                dist = sqrt((tx-cx)*(tx-cx) + (ty-cy)*(ty-cy))
                #g.variables['lon'][i,j] = tx
                #g.variables['lat'][i,j] = ty
        
                if (dist < 50.e3):
                        g.variables['thick'][i,j]=h_ice
                else:
                        g.variables['thick'][i,j]=0.0

for j in range(ny):
        g.variables['lat'][j] = float(j) * 10.e3

print('Creating file: ' , New_ice_thickness_filename)

g.sync()
g.close()

        
