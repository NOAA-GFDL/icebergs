echo 'Beginning setup...'



#Loading python libraries (its possible that not all of these are needed)
module unload PrgEnv-intel/5.2.8
module load PrgEnv-gnu/5.2.82
module load python_scipy/0.17.0

module load python
module load python_netcdf4
module load python_numpy
module load python_pygtk
module load python_matplotlib
module load python_scipy/0.17.0



git submodule init
git submodule update

cd ISOMIP_Calving/python_scripts
mkdir -p output_files

#Creating mosaic files using script
echo 'Creating mosaic files...'
./Create_mosaic.py -coupled_run -nx=240 -ny=40 -nz=36  -L=80.0 -W=480.0 
echo 'Creating mosaic files complete'

#Copying mosaic files to INPUT directory  (Files created by python script Create_mosaic)
cp output_files/atmos_mosaic_tile1Xland_mosaic_tile1.nc ../INPUT/
cp output_files/atmos_mosaic_tile1Xocean_mosaic_tile1.nc ../INPUT/
cp output_files/land_mosaic_tile1Xocean_mosaic_tile1.nc ../INPUT/
cp output_files/grid_spec.nc ../INPUT/
cp output_files/ocean_hgrid.nc ../INPUT/
cp output_files/ocean_mask.nc ../INPUT/
cp output_files/ocean_mosaic.nc ../INPUT/
cp output_files/ocean_topog.nc ../INPUT/
cp output_files/topog.nc ../INPUT/

#Creating iceberg and bond restart files, and ice shelf file
echo 'Creating iceberg restart files...'
./initialize_bergs_in_pattern.py -Make_icebergs_non_static_later=False -Run_plotting_subroutine=False
echo 'Creating iceberg restart files complete'

#Copying iceberg restarts and iceshelf file to INPUT directory  (Files created by python script initialize_bergs_in_pattern_Tech.py)
cp output_files/ISOMIP_icebergs.res.nc ../INPUT/icebergs.res.nc
cp output_files/ISOMIP_bonds_iceberg.res.nc ../INPUT/bonds_iceberg.res.nc
cp output_files/ISOMIP_Ice_Shelf_clipped.nc ../INPUT/isomip_ice_shelf1_clipped.nc

#Copying 3D_LAYER_WARM_TPY_IC.nc to INPUT
cp input_files/3D_LAYER_WARM_TPY_IC.nc ../INPUT/
cd ..
echo 'Setup complete!'
