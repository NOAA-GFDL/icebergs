#from midas.rectgrid import *
#from PIL import *
import netCDF4 as nc
import numpy as np
import os
import matplotlib as ml
import matplotlib.pyplot as plt

#These lines could be used to help select a subset of another data set. Not used for now.
#sgrid_path = '/archive/gold/datasets/OM4_025/mosaic.v20140610.unpacked/ocean_hgrid.nc'
#topog_path = '/archive/gold/datasets/OM4_025/mosaic.v20140610.unpacked/ocean_topog.nc'
#calving_path = '/home/Alon.Stern/Iceberg_Project/Weddell_Sea/subsampling_calving/example_calving.nc'
#sgrid=supergrid(file=sgrid_path,cyclic_x=True,tripolar_n=True)
#grid=quadmesh(supergrid=sgrid)
#location='Weddell'
#xextent={'Weddell':(-85.,20.2)}
#yextent={'Weddell':(grid.lath[0],-60.)}
#subreg=grid.geo_region(x=xextent[location],y=yextent[location])
#depth=nc.Dataset(topog_path).variables['depth'][subreg['y'],subreg['x']]
#ny,nx=depth.shape
#calving=nc.Dataset(calving_path).variables['CALVING']
#calving_time=nc.Dataset(calving_path).variables['time']



#Size of the domain
nx=420
ny=200

R=3.
x_pos_ind_1=200
y_pos_ind_1=180

x_pos_ind_2=400
y_pos_ind_2=180

Max_C=0.0058 #Maximum comes from the other calving flux file.

#Defining the amount of calving
calving=np.zeros((1,nx,ny))
C=np.zeros((nx,ny))
for i in xrange(0,nx):
    for j in xrange(1,ny):
        r_sq= (i-x_pos_ind_1)**2+(j-y_pos_ind_1)**2
        if r_sq  < 3*(R**2):
            C[i,j]=Max_C*np.exp(-(r_sq/(R**2))) 
        r_sq= (i-x_pos_ind_2)**2+(j-y_pos_ind_2)**2 
        if r_sq  < 3*(R**2):
            C[i,j]=Max_C*np.exp(-(r_sq/(R**2)))
   # calving[0,i,j]=C[i,j]



#Creating the calving file
fnam='calving_new.nc'
f=nc.Dataset(fnam,'w',format='NETCDF3_CLASSIC')
time=f.createDimension('time', None)
yt=f.createDimension('yt',ny)
xt=f.createDimension('xt',nx)

yt=f.createVariable('yt','f4',('yt'))
xt=f.createVariable('xt','f4',('xt'))
CALVING=f.createVariable('CALVING','f4',('time','yt','xt'))
time=f.createVariable('time','f4',('time'))

CALVING.long_name = "frozen runoff" ;
CALVING.units = "kg/(m^2*s)" ;
CALVING.missing_value =-10000;# -1.e+34f ;
CALVING.FillValue = -10000;#-1.e+34f ;
CALVING.cell_methods = "time: mean" ;
CALVING.time_avg_info = "average_T1,average_T2,average_DT" ;

time.long_name = "time" ;
time.units = "days since 1900-01-01 00:00:00" ;
time.cartesian_axis = "T" ;
time.calendar_type = "NOLEAP" ;
time.calendar = "NOLEAP" ;
time.bounds = "time_bounds" ;
time.modulo = " " ;

xt.long_name = "longitude" ;
xt.units = "degrees_E" ;
xt.cartesian_axis = "X" ;
xt.edges = "xb" ;

yt.long_name = "latitude" ;
yt.units = "degrees_N" ;
yt.cartesian_axis = "Y" ;
yt.edges = "yb" ;

f.filename = "19090101.ice_month.nc" ;
f.title = "CM2G_LM3_1990_DYNCODE_pi-control_C2_25sep2012_test" ;
f.grid_type = "regular" ;
f.grid_tile = "N/A" ;
f.history = "Mon Apr 22 14:33:03 2013: ncks -v time,xt,yt,CALVING 19090101.ice_month.nc example_calving.nc" ;
f.NCO = "4.0.3" ;


f.variables['time'][:]=0.
#f.variables['yt'][:]=yn
#f.variables['xt'][:]=xy
f.variables['CALVING'][:]=calving

print f.dimensions
for dimobj in f.dimensions.values():
    print dimobj
f.sync()
f.close()



H = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])  # added some commas and array creation code
#H=np.zeros((nx,ny))
H=C
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
ax.set_title('colorMap')
plt.imshow(H)
ax.set_aspect('equal')
cax = fig.add_axes([0.12, 0.1, 0.78, 0.8])
cax.get_xaxis().set_visible(False)
cax.get_yaxis().set_visible(False)
cax.patch.set_alpha(0)
cax.set_frame_on(False)
plt.colorbar(orientation='vertical')
plt.show()
