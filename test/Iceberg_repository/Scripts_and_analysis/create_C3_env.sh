#This script creates a C3 environment

mkdir build/intel
cat << EOFA > build/intel/env
module unload PrgEnv-pgi
module unload PrgEnv-pathscale
module unload PrgEnv-intel
module unload PrgEnv-gnu
module unload PrgEnv-cray
module load PrgEnv-intel
module unload netcdf
module load cray-netcdf
EOFA

