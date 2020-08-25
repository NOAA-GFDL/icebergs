#!/bin/csh -fx
#PBS -N Tech_Bergs
#PBS -l walltime=16:00:00
#PBS -l size=128
#PBS -S /bin/tcsh
#PBS -r n
#PBS -m ae
#PBS -j oe
#PBS -E
#PBS -A gfdl_o

#set compiler = intel


module load totalview
env
pwd
time aprun -n 128 ../../build/intel/ice_ocean_SIS2/repro/MOM6
#date

#time aprun ...

#date
