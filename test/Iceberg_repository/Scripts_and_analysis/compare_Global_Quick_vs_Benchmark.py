!I still need to write this code comparing the files. 


from midas.rectgrid import *
from PIL import *
import netCDF4 as nc
import numpy as np
import os


sgrid_path = '/archive/gold/datasets/OM4_025/mosaic.v20140610.unpacked/ocean_hgrid.nc'
topog_path = '/archive/gold/datasets/OM4_025/mosaic.v20140610.unpacked/ocean_topog.nc'
tideamp_path = '/archive/gold/datasets/OM4_025/INPUT/tidal_amplitude.v20140616.nc'
seawifs_path = '/archive/gold/datasets/OM4_025/INPUT/seawifs_1998-2006_smoothed_2X.v20140629.nc'

sgrid=supergrid(file=sgrid_path,cyclic_x=True,tripolar_n=True)
grid=quadmesh(supergrid=sgrid)

location='Weddell'
xextent={'Weddell':(-85.,20.2)}
yextent={'Weddell':(grid.lath[0],-60.)}

subreg=grid.geo_region(x=xextent[location],y=yextent[location])
h2=nc.Dataset(topog_path).variables['h2'][subreg['y'],subreg['x']]
depth=nc.Dataset(topog_path).variables['depth'][subreg['y'],subreg['x']]
tideamp=nc.Dataset(tideamp_path).variables['tideamp'][subreg['y'],subreg['x']]


fnam='topog.nc'
f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')

ny,nx=depth.shape

f.createDimension('ny',ny)
f.createDimension('nx',nx)
f.createDimension('ntiles',1)
f.createVariable('depth','f8',('ny','nx'))
f.createVariable('h2','f8',('ny','nx'))
f.variables['depth'][:]=depth
f.variables['h2'][:]=h2
f.sync()
f.close()


fnam='tideamp.nc'
f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')
ny,nx=tideamp.shape
f.createDimension('ny',ny)
f.createDimension('nx',nx)
f.createVariable('tideamp','f4',('ny','nx'))
f.variables['tideamp'][:]=tideamp
f.sync()
f.close()


fnam='ocean_hgrid.nc'
f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')
f.createDimension('ny',ny*2)
f.createDimension('nx',nx*2)
f.createDimension('nyp',ny*2+1)
f.createDimension('nxp',nx*2+1)
f.createDimension('string',255)
f.createVariable('y','f8',('nyp','nxp'))
f.createVariable('x','f8',('nyp','nxp'))
f.createVariable('dy','f8',('ny','nxp'))
f.createVariable('dx','f8',('nyp','nx'))
f.createVariable('area','f8',('ny','nx'))
f.createVariable('angle_dx','f8',('nyp','nxp'))
f.createVariable('tile','S1',('string'))

x_read=subreg['x']
y_read=subreg['y']


if x_read[-1] < x_read[0]:
    ib=np.where(np.roll(x_read,shift=-1)<x_read)[0][0]
    tmp = np.arange(x_read[0]*2,x_read[ib]*2+2)
    x_read=np.concatenate((tmp,np.arange(x_read[ib+1]*2,x_read[-1]*2+3)))
    y_read = np.arange(y_read[0]*2,y_read[-1]*2+3)    
else:
    x_read = np.arange(x_read[0]*2,x_read[-1]*2+3)
    y_read = np.arange(y_read[0]*2,y_read[-1]*2+3)


f.variables['y'].units='degrees'
f.variables['x'].units='degrees'
f.variables['dy'].units='meters'
f.variables['dx'].units='meters'
f.variables['area'].units='m2'
f.variables['angle_dx'].units='degrees'

yn=sgrid.y[y_read,:]
yn=yn[:,x_read]
f.variables['y'][:]=yn
xn=sgrid.x[y_read,:]
xn=xn[:,x_read]
f.variables['x'][:]=xn
dyn=sgrid.dy[y_read[:-1],:]
dyn=dyn[:,x_read]
f.variables['dy'][:]=dyn
dxn=sgrid.dx[y_read,:]
dxn=dxn[:,x_read[:-1]]
f.variables['dx'][:]=dxn
arean=sgrid.area[y_read[:-1],:]
arean=arean[:,x_read[:-1]]
f.variables['area'][:]=arean
ang=sgrid.angle_dx[y_read,:]
ang=ang[:,x_read]
f.variables['angle_dx'][:]=ang
f.variables['tile'][0] = 't'  ## This is stupid
f.variables['tile'][1] = 'i'
f.variables['tile'][2] = 'l'
f.variables['tile'][3] = 'e'
f.variables['tile'][4] = '1'
f.sync()
f.close()





