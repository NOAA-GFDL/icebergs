
# coding: utf-8

# In[2]:

import scipy.io.netcdf
import numpy
get_ipython().magic('pylab inline')


# Model grid will have ni x nj h-point grid cells

# In[3]:

nj,ni = 40,240


# A mosaic file has a supergrid of an integer refinement (e.g. 2*ni x 2*nj cells)

# In[4]:

snj, sni = 2*nj, 2*ni


# Uniform resolution of the ocean model grid will be dxo x dyo

# In[5]:

dxo = 2.e3
dyo = 2.e3


# Using a supergrid refinement of 2, the grid lengths, area and coordinates are:

# In[6]:

dx = (0.5 * dxo) * numpy.ones((snj+1,sni))
dy = (0.5 * dyo) * numpy.ones((snj,sni+1))
area = (0.25 * (dxo * dyo)) * numpy.ones((snj,sni))
x = numpy.zeros((snj+1,sni+1))
x[:,1:] = numpy.cumsum( dx, axis=-1 )
y = numpy.zeros((snj+1,sni+1))
y[1:,:] = numpy.cumsum( dy, axis=-2 )


# In[7]:

# A helper function to use when writing strings in a netcdf file
def set_string(variable, value):
    """Sets "variable" to "value" padded with blanks where
    "variable" is a netcdf variable object and "value" is a string."""
    variable[:] = '\000' * variable.shape[0]
    variable[:len(value)] = value


# Create a netcdf3 horizontal ocean-grid file...

# In[8]:

rg = scipy.io.netcdf_file('ocean_hgrid.nc','w',version=2)


# In[9]:

rg.createDimension('nx',sni)
rg.createDimension('ny',snj)
rg.createDimension('nxp',sni+1)
rg.createDimension('nyp',snj+1)
rg.createDimension('string',255)


# In[10]:

dx_h = rg.createVariable('dx','f8',('nyp','nx',))
dx_h.units = 'm'
dy_h = rg.createVariable('dy','f8',('ny','nxp',))
dy_h.units = 'm'
area_h = rg.createVariable('area','f8',('ny','nx',))
area_h.units = 'm2'
x_h = rg.createVariable('x','f8',('nyp','nxp',))
x_h.units = 'm'
y_h = rg.createVariable('y','f8',('nyp','nxp',))
y_h.units = 'm'
tile = rg.createVariable('tile','c',('string',))


# In[11]:

dx_h[:,:] = dx[:,:]
dy_h[:,:] = dy[:,:]
area_h[:,:] = area[:,:]
x_h[:,:] = x[:,:]
y_h[:,:] = y[:,:]
tile[:] = '\000' * 255
tile[:5] = 'tile1'
set_string(tile,'tile1')


# In[12]:

rg.close()


# Create a netcdf3 ocean-mask file...

# In[13]:

rg = scipy.io.netcdf_file('ocean_mask.nc','w',version=2)
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
mask = rg.createVariable('mask','f4',('ny','nx',))
mask[:,:] = 1.
rg.close()


# Create a mosaic description for the above grid file...

# In[14]:

rg = scipy.io.netcdf_file('ocean_mosaic.nc','w',version=2)
rg.createDimension('ntiles',1)
rg.createDimension('string',255)
mosaic = rg.createVariable('mosaic','c',('string',))
mosaic.standard_name = 'grid_mosaic_spec'
mosaic.children = 'contacts'
mosaic.grid_descriptor = ''
gridlocation = rg.createVariable('gridlocation','c',('string',))
gridlocation.standard_name = 'grid_file_location'
gridfiles = rg.createVariable('gridfiles','c',('ntiles','string',))
gridtiles = rg.createVariable('gridtiles','c',('ntiles','string',))
rg.grid_version = '0.2'
# Fill in data
mosaic[:] = '\000' * 255
mosaic[:12] = 'ocean_mosaic'
gridlocation[:] = '\000' * 255
gridlocation[:2] = './'
gridfiles[:] = '\000' * 255
gridfiles[0,:14] = 'ocean_hgrid.nc'
gridtiles[:] = '\000' * 255
gridtiles[0,:5] = 'tile1'
rg.close()


# In[15]:

get_ipython().system('pwd')


# In[16]:

# Fake a ocean_topog.nc file
rg = scipy.io.netcdf_file('ocean_topog.nc','w',version=2)
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
rg.createDimension('ntiles',1)
depth = rg.createVariable('depth','f4',('ny','nx',))
depth[:,:] = 1000.
depth[-1,-1] =0
rg.close()


# In[17]:

#fre_nctools/tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ocean_topog.nc


# In[18]:

rg = scipy.io.netcdf_file('atmos_mosaic_tile1Xland_mosaic_tile1_BLAH.nc','w',version=2) # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',1)
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
tile1_cell[:] = [ni,nj]
tile2_cell[:] = [ni,nj]
xgrid_area[:] = dxo * dyo
tile1_distance[:] = 0.
tile2_distance[:] = 0.
rg.close()


# In[19]:

rg = scipy.io.netcdf_file('atmos_mosaic_tile1Xocean_mosaic_tile1_BLAH.nc','w',version=2) # atmos_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-1) # -1 is for a single land point
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
tile1_cell[:] = [ni,nj]
tile2_cell[:] = [ni,nj]
xgrid_area[:] = dxo * dyo
for j in range(nj):
    for i in range(ni):
        if (i,j) != (ni-1,nj-1):
            tile1_cell[i+ni*j] = [i,j]
            tile2_cell[i+ni*j] = [i,j]
            tile1_distance[i+ni*j] = [0,0]
            tile2_distance[i+ni*j] = [0,0]
            xgrid_area[i+ni*j] = dxo * dyo
rg.close()


# In[20]:

rg = scipy.io.netcdf_file('land_mosaic_tile1Xocean_mosaic_tile1_BLAH.nc','w',version=2) # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-1) # -1 is for a single land point
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','f8',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','f8',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
tile1_cell[:] = [ni,nj]
tile2_cell[:] = [ni,nj]
xgrid_area[:] = dxo * dyo
for j in range(nj):
    for i in range(ni):
        if (i,j) != (ni-1,nj-1):
            tile1_cell[i+ni*j] = [i,j]
            tile2_cell[i+ni*j] = [i,j]
            tile1_distance[i+ni*j] = [0,0]
            tile2_distance[i+ni*j] = [0,0]
            xgrid_area[i+ni*j] = dxo * dyo
rg.close()


# In[24]:

rg = scipy.io.netcdf_file('mosaic_BLAH.nc','w',version=2) # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('nfile_aXo',1) # -1 is for a single land point
rg.createDimension('nfile_aXl',1) # -1 is for a single land point
rg.createDimension('nfile_lXo',1) # -1 is for a single land point

atm_mosaic_dir = rg.createVariable('atm_mosaic_dir','c',('string',))
atm_mosaic_dir.standard_name = 'directory_storing_atmosphere_mosaic'

atm_mosaic_file = rg.createVariable('atm_mosaic_file','c',('string',))
atm_mosaic_file.standard_name = 'atmosphere_mosaic_file_name'

atm_mosaic = rg.createVariable('atm_mosaic','c',('string',))
atm_mosaic.standard_name = 'atmosphere_mosaic_name'

lnd_mosaic_dir = rg.createVariable('lnd_mosaic_dir','c',('string',))
lnd_mosaic_dir.standard_name = 'directory_storing_land_mosaic'

lnd_mosaic_file = rg.createVariable('lnd_mosaic_file','c',('string',))
lnd_mosaic_file.standard_name = 'land_mosaic_file_name'

lnd_mosaic = rg.createVariable('lnd_mosaic','c',('string',))
lnd_mosaic.standard_name = 'land_mosaic_name'

ocn_mosaic_dir = rg.createVariable('ocn_mosaic_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_mosaic'

ocn_mosaic_file = rg.createVariable('ocn_mosaic_file','c',('string',))
ocn_mosaic_file.standard_name = 'ocean_mosaic_file_name'

ocn_mosaic = rg.createVariable('ocn_mosaic','c',('string',))
ocn_mosaic.standard_name = 'ocean_mosaic_name'


ocn_topog_dir = rg.createVariable('ocn_topog_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_topog'

ocn_topog_file = rg.createVariable('ocn_topog_file','c',('string',))
ocn_topog_file.standard_name = 'ocean_topog_file_name'

aXo_file = rg.createVariable('aXo_file','c',('string',))
aXo_file.standard_name = 'atmXocn_exchange_grid_file'

aXl_file = rg.createVariable('aXl_file','c',('string',))
aXl_file.standard_name = 'atmXlnd_exchange_grid_file'

lXo_file = rg.createVariable('lXo_file','c',('string',))
lXo_file.standard_name = 'lndXocn_exchange_grid_file'

#Global attributes
rg.grid_version = '0.2'
rg.code_version = "$Name:  $"
rg.history = "/net2/aja/workspace/MOM6-examples/ice_ocean_SIS2/OM4_025/preprocessing/fre_nctools/tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --ocean_topog ocean_topog.nc"
# Fill in data
atm_mosaic_dir[:] = '\000' * 255
atm_mosaic_dir[:2] = './'

atm_mosaic_file[:] = '\000' * 255
atm_mosaic_file[:15] = 'ocean_mosaic.nc'

atm_mosaic[:] = '\000' * 255
atm_mosaic[:12] = 'atmos_mosaic'

lnd_mosaic_dir[:] = '\000' * 255
lnd_mosaic_dir[:2] = './'

lnd_mosaic_file[:] = '\000' * 255
lnd_mosaic_file[:15] = 'ocean_mosaic.nc'

lnd_mosaic[:] = '\000' * 255
lnd_mosaic[:11] = 'land_mosaic'

ocn_mosaic_dir[:] = '\000' * 255
ocn_mosaic_dir[:2] = './'

ocn_mosaic_file[:] = '\000' * 255
ocn_mosaic_file[:15] = 'ocean_mosaic.nc'

ocn_mosaic[:] = '\000' * 255
ocn_mosaic[:12] = 'ocean_mosaic'

ocn_topog_dir[:] = '\000' * 255
ocn_topog_dir[:2] = './'

ocn_topog_file[:] = '\000' * 255
ocn_topog_file[:14] = 'ocean_topog.nc'

aXo_file[:] = '\000' * 255
aXo_file[:40] = 'atmos_mosaic_tile1Xocean_mosaic_tile1.nc'

aXl_file[:] = '\000' * 255
aXl_file[:39] = 'atmos_mosaic_tile1Xland_mosaic_tile1.nc'

lXo_file[:] = '\000' * 255
lXo_file[:39] = 'land_mosaic_tile1Xocean_mosaic_tile1.nc'

rg.close()



# In[ ]:



