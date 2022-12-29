#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("tkagg")
from pylab import *
from scipy import interpolate
#import matplotlib.pyplot as plt
import argparse
import pdb
import netCDF4 as nc


import sys
pypath = '.'
for dir_name in os.listdir(pypath):
        dir_path = os.path.join(pypath, dir_name)
        if os.path.isdir(dir_path):
                sys.path.insert(0, dir_path)

from hexagon_area import Divide_hexagon_into_4_quadrants_old
from hexagon_area import Hexagon_into_quadrants_using_triangles


def parseCommandLine():
        """
        Parse the command line positional and optional arguments.
        This is the highest level procedure invoked from the very end of the script.
        """

        parser = argparse.ArgumentParser(description=
        '''
        Generate files for iceberg and bond restart files, matching given ice thickness fields.
        ''',
        epilog='Written by Alex Huth, 2021, based on Iceberg_repository code from Alon Stern.')

        #Adding an extra boolian type argument
        #p = argparse.ArgumentParser()
        parser.register('type','bool',str2bool) # add type keyword to registries
        #p.add_argument('-b',type='bool')  # do not use 'type=bool'

        #Flags
        parser.add_argument('-save_restart_files', type='bool', default=True,
                          help=''' Writes an iceberg restart file,  icebergs.res ''')

        parser.add_argument('-save_new_h_ice_file', type='bool', default=True,
                          help=''' Writes a new gridded ice thickness file, which is what the ice thickness after the berg mass is interpolated back onto the grid''')

        parser.add_argument('-Convert_to_lat_lon', type='bool', default=False,
                          help=''' Converts the iceberg positions to lat / lon coordinates  ''')

        parser.add_argument('-input_is_cartesian', type='bool', default=True,
                          help=''' The positions in the ice thickness input file is given in cartesian coordiantes''')
        parser.add_argument('-lat_ref_grid_to_cartesian', type=float, default=-9999,
                          help='''Reference latitude for converting to cartesian grid   ''')

        parser.add_argument('-A68', type='bool', default=False,
                          help=''' Use velocities and angular velocities for A68 experiment''')

        parser.add_argument('-displace_y', type=float, default=0.0,
                          help='''displace the y coords of the berg (meters)   ''')

        parser.add_argument('-ang_vel_prescribed', type=float, default=0.0,
                          help='''Prescibed angular velocity   ''')

        parser.add_argument('-dx', type=float, default=250.0,
                          help='''dx for regrid   ''')

        parser.add_argument('-icethres', type=float, default=0.5,
                          help='''threshold for whether a grid cell is ice or not   ''')

        #Iceberg setup flags
        parser.add_argument('-only_choose_one_berg', type='bool', default=False,
                          help=''' When true, only one iceberg (with number chosen_berg_num) is written to icebergs.res. This is for debugging  ''')

        parser.add_argument('-chosen_berg_num', type=int, default=1,
                          help='''When only_choose_one_berg=True, then only the iceberg with this number is written to the icebergs.res file  ''')

        parser.add_argument('-scale_the_grid_to_lat_lon', type='bool', default=False,
                          help=''' This option scales dimensions by a fixed amount. Option should be removed   ''')

        parser.add_argument('-adjust_lat_ref', type='bool', default=True,
                          help=''' This option is only used when Convert_to_lat_lon is true. If true, then the reference latitudes depends on iceberg position.\
                                          This still needs to be tested using lat lon test cases.''')

        parser.add_argument('-set_all_thicknesses_to_one', type='bool', default=False,
                          help=''' Sets all non-zero ice thickness to Th_prescibed (default 1), useful for debugging.  ''')

        parser.add_argument('-Th_prescribed', type=float, default=1.0,
                          help='''Prescibed ice thickness used when the flag set_all_thicknesses_to_one=True   ''')

        parser.add_argument('-uvel', type=float, default=0.0,
                          help='''Prescibed ice u velocity (m/s)''')

        parser.add_argument('-vvel', type=float, default=0.0,
                          help='''Prescibed ice v velocity (m/s)''')

        parser.add_argument('-Interpolate_from_four_corners', type='bool', default=True,
                          help=''' The thickness of an iceberg is calcaulated from the 4 corners of the gridded thickness. When flag is false, the near thickness is used  ''')

        parser.add_argument('-Switch_x_and_y_to_rotate_90', type='bool', default=False,
                          help=''' Rotates the whole domain by 90 degrees by switching lat and lon of icebergs.  ''')

        parser.add_argument('-set_all_domain_to_ice', type='bool', default=False,
                          help=''' Sets the ice thickness = Th_prescibed (default 1) throughout the entire domain. Applied only when flag set_all_thicknesses_to_one=True   ''')

        parser.add_argument('-Use_default_radius', type='bool', default=False,
                          help=''' Calculates an appropriate ice element radius based on the ocean grid spacing.   ''')

        parser.add_argument('-Remove_stationary_bergs', type='bool', default=False,
                          help=''' All stationary icebergs are removed.  ''')

        #Static vs non-static flags (These are needed for calving away from Static shelf)
        parser.add_argument('-set_all_bergs_static_by_default', type='bool', default=True,
                        help='''Bergs are static unless instructed otherwise (e.g.: after applying calving)   ''')

        parser.add_argument('-set_some_bergs_static_by_default', type='bool', default=False,
                          help=''' Bergs in part of the domain are set to static by default (see subroutine for details)  ''')

        parser.add_argument('-Make_icebergs_non_static_later', type='bool', default=False,
                        help=''' Run subroutine which makes some icebergs non-static later (eg: icebergs around calved tabular icebergs).  Note that bonds are \
                                        only created for non_static, so by making bergs non_static later, it means they can move, but are not bonded.  ''')

        #Mass spreadng flags
        parser.add_argument('-Switch_regridding_element_type', type='bool', default=False,
                          help=''' When true, element shape is switched from hexagon to square (and visa versa) before regridding. Used in debugging to see if having \
                          the wrond elemnt type makes a bid difference. Maybe this option should be removed(?)''')

        parser.add_argument('-Fill_in_the_boundaries', type='bool', default=True,
                          help=''' Adds extra (static)  ice elements at the edges of the domain when the ice shelf extends right up until the boundaries of the domain.\
                                          Only used with hexagonal ice elements (for now).''')

        parser.add_argument('-regrid_icebergs_onto_grid', type='bool', default=True,
                          help=''' After icebergs have been defined, the iceberg mass is interpolated back onto the ocean grid  ''')

        parser.add_argument('-element_type', type=str, default='hexagon',
                        help='''Shape of element used in mass spreading. Options are: 'hexagon' or 'square'   ''')


        #Plotting flags - Only affect how the ice elements are plotted and does not affect files created by this script
        parser.add_argument('-Run_plotting_subroutine', type='bool', default=True,
                          help=''' Subrouting is run plots icebergs and ice thickness. All other plotting flags require this to be true.  ''')

        parser.add_argument('-plot_circles', type='bool', default=False,
                          help=''' Circles are plotted to scale which show how large each ice element is. (Useful when converting to lat lon)  ''')

        parser.add_argument('-plot_ice_mask', type='bool', default=False,
                          help=''' The ice mask (showing where ice shelf is present) is plotted under iceberg positions. Option may be removed.  ''')

        parser.add_argument('-plot_ice_thickness', type='bool', default=True,
                          help=''' The ice thickness is plotted under iceberg positions   ''')

        parser.add_argument('-plot_icebergs_positions', type='bool', default=True,
                          help='''  Positions of ice elements are plotted  ''')

        parser.add_argument('-plot_h_ice_new', type='bool', default=True,
                          help='''  The newly interpolated ice thickness is plotted under iceberg positions ''')

        parser.add_argument('-plot_bonds', type='bool', default=False,
                          help=''' Bonds are plotted between the ice elements   ''')

        #Bond related flags
        parser.add_argument('-Create_icebergs_bonds', type='bool', default=True,
                          help=''' Bonds are created between icebergs  ''')

        parser.add_argument('-break_some_bonds', type='bool', default=True,
                          help=''' After the bonds have been set up, some are broken. This simulates a iceberg / ice shelf calving.\
                                          Note that this flag is also used to make a calving event even if bonds are not used since \
                                          the "calving" process also sets some icebergs to static and others to non-static''')

        parser.add_argument('-Allow_bonds_for_static_iceberg', type='bool', default=False,
                          help=''' If True then bonds are allowed for ice elements which are static.   ''')

        parser.add_argument('-Allow_bonds_with_boundary_bergs', type='bool', default=False,
                          help=''' When true, bonds are allowed to form with elements close to the boundary. \
                                          When false, bonds are allowed for bergs closer with i,j < N_bergs_before_bd.\
                                          This option is a little confusing and might need to be removed (>)''')

        #Experimental setup flags
        parser.add_argument('-Ice_geometry_source', type=str, default='ISOMIP',
                        help=''' Name of experiment which is being run. Options are 'ISOMIP', 'Generic' and 'Weddell'\
                                        These experiments override some of the flags above and have standary input files  (see below).\
                                        This will be made more general later so that ice thickness files are passed into the script with standarized formate through driver.\
                                        ''')

        parser.add_argument('-ISOMIP_ice_geometry_filename', type=str, default='input_files/Isomip_ice_geometry.nc',
                        help=''' Ice thickness from the ISOMIP experiment, given on ice grid, rather than ocean grid.  ''')

        parser.add_argument('-ISOMIP_reduced_ice_geometry_filename', type=str, default='input_files/Ocean1_3D_no_calving_trimmed.nc',
                        help=''' Ice thickness from ISOMIP experiment, given on ocean grid. Later all input files will be made to have this format\
                                        (or to specify if they are on the ice grid, what format the ocean should have)''')

        parser.add_argument('-Weddell_ice_geometry_filename', type=str, default='input_files/Bedmap2_gridded_subset_Weddell_Sea_Region.nc',
                        help='''  Ice thickness from the ice shelves in the Weddell Sea, given on ice grid, rather than ocean grid.  ''')

        parser.add_argument('-Generic_ice_geometry_filename', type=str, \
                        default='/lustre/f1/unswept/Alon.Stern/MOM6-examples_Alon/ice_ocean_SIS2/Drifting_tabular/python_scripts/output_files/Ice_shelf_file.nc',
                        help=''' Ice thickness used in a generic setup.   ''')


        parser.add_argument('-ISOMIP_reduced', type='bool', default=True,
                        help=''' Flag which specifies that we are using the reduced ISOMIP file (on the ocean grid), rather than geomety on the ice grid. \
                                        Reduced uses 2X2 grid, not reduced uses 1X1 grid \
                                        This should be made more general later.''')

        #Parameters
        parser.add_argument('-Radius', type=float, default='850.0',
                        help='''Radius of ice element (m). This is overriden if Use_default_radius=True\
                                        Note that Hexagon only valid if the Radius, S< half gridcell  (about 0.85 using 2km grid)''')

        parser.add_argument('-rho_ice', type=float, default=850.0,
                        help=''' Density of icebergs (kg/m^3)  ''')

        parser.add_argument('-gravity', type=float, default=9.8,
                        help=''' Gravitational acceleration  (m/s^2)''')

        parser.add_argument('-mass_scaling', type=float, default=1.0,
                        help=''' Number of icebergs represented by one ice element  ''')

        #parser.add_argument('-R_earth', type=float, default=6360000.,
        #                 help=''' Radius of the earch (m) - used in Lat/Lon conversions.   ''')

        parser.add_argument('-R_earth', type=float, default=6378000.,
                        help=''' Radius of the earch (m) - used in Lat/Lon conversions.   ''')

        parser.add_argument('-buffer_number', type=int, default=0,
                        help=''' Amount of points from the boundaries where ice thickness is set to zero for debugging. This should be removed.  ''')

        parser.add_argument('-IA_scaling', type=float, default=1.,
                        help='''  A scaling parameter which allows the interactive radius be different from radius for testing the bonds.\
                                        (This has not been used in a long while, and perhaps should be removed)''')

        optCmdLineArgs = parser.parse_args()
        return optCmdLineArgs

def str2bool(string):
        if string.lower() in  ("yes", "true", "t", "1"):
                Value=True
        elif string.lower() in ("no", "false", "f", "0"):
                Value=False
        else:
                print( '**********************************************************************')
                print( 'The input variable ' ,str(string) ,  ' is not suitable for boolean conversion, using default')
                print( '**********************************************************************')

                Value=None
                return

        return Value



def Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,ang_vel_in,uvel_out,vvel_out,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg):

        print( 'Writing iceberg restart files, with ' , Number_of_bergs  , 'icebergs..')
        # To copy the global attributes of the netCDF file

        #Input and output files
        #Create Empty restart file. This is later read so that the attributes can be used.
        Empty_restart_filename='output_files/Empty_icebergs.res.nc'
        create_empty_iceberg_restart_file(Empty_restart_filename)
        #Empty_restart_filename='input_files/icebergs.res.nc'

        #Read empty restart file
        f=Dataset(Empty_restart_filename,'r') # r is for read only
        #Write a new restart file
        print('Ice_geometry_source',Ice_geometry_source)
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

                if varname=='ang_vel':
                        for j in range(Number_of_bergs):
                                var[j]=ang_vel_in #1.0.#ang_vel_out[j]

                if varname=='uvel':
                        for j in range(Number_of_bergs):
                                var[j]=uvel_out[j]

                if varname=='vvel':
                        for j in range(Number_of_bergs):
                                var[j]=vvel_out[j]

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



def Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source):
        #Creating the bond restart file

        print( 'Writing bond restart files with  ', Number_of_bonds, ' Bonds')
        # To copy the global attributes of the netCDF file

        #Input and output files
        #Create Empty restart file. This is later read so that the attributes can be used.
        Empty_bond_restart_filename='output_files/Empty_bonds_icebergs.res.nc'
        create_empty_bond_restart_file(Empty_bond_restart_filename)
        #Empty_bond_restart_filename='input_files/bonds_iceberg.res.nc'

        h=Dataset(Empty_bond_restart_filename,'r') # r is for read only
        q=Dataset('output_files/' + Ice_geometry_source + '_bonds_iceberg.res.nc','w', format='NETCDF3_CLASSIC') # w if for creating a file

        for attname in h.ncattrs():
                    setattr(q,attname,getattr(h,attname))


        # To copy the dimension of the netCDF file
        for dimname,dim in h.dimensions.iteritems():
                # if you want to make changes in the dimensions of the new file
                # you should add your own conditions here before the creation of the dimension.
                #g.createDimension(dimname,len(dim))
                q.createDimension(dimname,Number_of_bonds)

        # To copy the variables of the netCDF file

        for varname,ncvar in h.variables.iteritems():
                # if you want to make changes in the variables of the new file
                # you should add your own conditions here before the creation of the variable.
                var = q.createVariable(varname,ncvar.dtype,ncvar.dimensions)
                #Proceed to copy the variable attributes
                for attname in ncvar.ncattrs():
                        setattr(var,attname,getattr(ncvar,attname))
                #Finally copy the variable data to the new created variable
                #var[:] = ncvar[0]
                var[:] = 0.

                if varname=='i':
                        var[:]=Number_of_bonds

                if varname=='first_berg_num':
                        for j in range(Number_of_bonds):
                                var[j]=first_berg_num[j]

                if varname=='first_berg_ine':
                        for j in range(Number_of_bonds):
                                var[j]=first_berg_ine[j]

                if varname=='first_berg_jne':
                        for j in range(Number_of_bonds):
                                var[j]=first_berg_jne[j]

                if varname=='other_berg_num':
                        for j in range(Number_of_bonds):
                                var[j]=other_berg_num[j]

                if varname=='other_berg_ine':
                        for j in range(Number_of_bonds):
                                var[j]=other_berg_ine[j]

                if varname=='other_berg_jne':
                        for j in range(Number_of_bonds):
                                var[j]=other_berg_jne[j]

        h.close()
        q.close()


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

        ang_vel=f.createVariable('ang_vel','d',('i'))
        ang_vel.long_name = "angular velocity" ;
        ang_vel.units = "rad/s" ;
        ang_vel.checksum = "               0" ;

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

def create_empty_bond_restart_file(Empty_bond_restart_filename):

        f = Dataset(Empty_bond_restart_filename,'w', format='NETCDF3_CLASSIC')

        i=f.createDimension('i', None)
        i=f.createVariable('i','i')

        first_berg_ine=f.createVariable('first_berg_ine','i',('i'))
        first_berg_ine.long_name = "iceberg ine of first berg in bond" ;
        first_berg_ine.units = "dimensionless" ;
        first_berg_ine.packing = 0 ;
        first_berg_ine.checksum = "               0" ;

        first_berg_jne=f.createVariable('first_berg_jne','i',('i'))
        first_berg_jne.long_name = "iceberg jne of first berg in bond" ;
        first_berg_jne.units = "dimensionless" ;
        first_berg_jne.packing = 0 ;
        first_berg_jne.checksum = "               0" ;

        first_berg_num=f.createVariable('first_berg_num','i',('i'))
        first_berg_num.long_name = "iceberg id first berg in bond" ;
        first_berg_num.units = "dimensionless" ;
        first_berg_num.packing = 0 ;
        first_berg_num.checksum = "               0" ;

        other_berg_ine=f.createVariable('other_berg_ine','i',('i'))
        other_berg_ine.long_name = "iceberg ine of second berg in bond" ;
        other_berg_ine.units = "dimensionless" ;
        other_berg_ine.packing = 0 ;
        other_berg_ine.checksum = "               0" ;

        other_berg_jne=f.createVariable('other_berg_jne','i',('i'))
        other_berg_jne.long_name = "iceberg jne of second berg in bond" ;
        other_berg_jne.units = "dimensionless" ;
        other_berg_jne.packing = 0 ;
        other_berg_jne.checksum = "                0" ;

        other_berg_num=f.createVariable('other_berg_num','i',('i'))
        other_berg_num.long_name = "iceberg id second berg in bond" ;
        other_berg_num.units = "dimensionless" ;
        other_berg_num.packing = 0 ;
        other_berg_num.checksum = "               0" ;

        f.sync()
        f.close()

def Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,rho_ice,Radius,h_ice,A68,uvel_in,vvel_in,x,y,\
                width,Interpolate_from_four_corners,element_area,element_type,static_berg):
        thickness=[]
        mass=[]
        uvel_out=[]
        vvel_out=[]
        ang_vel_out=[]
        ny,nx=h_ice.shape
        dx=x[1]-x[0]
        grid_area=dx*dx
        for berg_count in range(Number_of_bergs):
                x_val=dx_berg[berg_count] ; y_val=dy_berg[berg_count]
                i_val=int(floor(x_val/dx))
                j_val=int(floor(y_val/dx))

                #Interpolate thickness from 4 corners - possibly do this later, but I worry about when you are between a shelf and non shelf piece.
                if Interpolate_from_four_corners==True:

                        x_cell=(x_val-x[i_val])/sqrt(grid_area)+0.5
                        y_cell=(y_val-y[j_val])/sqrt(grid_area)+0.5

                        mass_on_ocean=np.zeros([nx,ny,10])   #Setting up matrix to spread mass to ocean.  Note that I have used 10 points  so that I can ignore 0 and match with python numbering
                        mass_val=1.
                        mass_on_ocean=spread_mass_to_ocean(nx,ny,i_val,j_val,mass_on_ocean,x_cell,y_cell,element_area,mass_val,element_type,grid_area,static_berg[berg_count])
                        Th=0.
                        h=np.array([0,0,0,0,0,0,0,0,0,0])
                        if i_val>0 and j_val>0:
                                h[1]=h_ice[j_val-1,i_val-1]
                        if j_val>0:
                                h[2]=h_ice[j_val-1,i_val]
                        if j_val>0 and i_val<(nx-1):
                                h[3]=h_ice[j_val-1,i_val+1]
                        if i_val>0:
                                h[4]=h_ice[j_val,i_val-1]
                        if True:
                                h[5]=h_ice[j_val,i_val]
                        if i_val<(nx-1):
                                h[6]=h_ice[j_val,i_val+1]
                        if i_val>0 and j_val<(ny-1):
                                h[7]=h_ice[j_val+1,i_val-1]
                        if j_val<(ny-1):
                                h[8]=h_ice[j_val+1,i_val]
                        if i_val<(nx-1) and j_val<(ny-1):
                                h[9]=h_ice[j_val+1,i_val+1]
                        for k in range(1,10):
                                Th=Th+(mass_on_ocean[i_val,j_val,k]*h[k])
                                #print i_val,j_val,k,mass_on_ocean[i_val,j_val,k]
                        #if abs(Th-1)>0.0000001:
                        #       print 'Thickness',Th-1,i_val,j_val
                        #       test=0.
                        #       for k in range(1,10):
                        #               print k,mass_on_ocean[i_val,j_val,k], h[k]
                        #               test=test+(mass_on_ocean[i_val,j_val,k])
                        #       print test
                        #       halt


                Th=h_ice[j_val,i_val]
                thickness.append(Th)

                if (A68):
                        uvel_out.append(uvel_in)
                        vvel_out.append(vvel_in)
                        ang_vel_out.append(0)
                else:
                        uvel_out.append(uvel_in)
                        vvel_out.append(vvel_in)
                        ang_vel_out.append(0)
                #Check that all thicknesses are positive
                if Th<=0:
                        print ('Thickness is less than or equal to zero,0!',i_val, j_val,x_val,y_val)
                        halt
                #mass.append(Th*rho_ice*element_area)
                mass.append(Th*rho_ice*(width[berg_count])**2)
                #if element_type=='square':
                #       mass.append(Th*rho_ice*((2*Radius)**2)) # For square elements
                #else:
                #       mass.append(Th*rho_ice*np.pi*(Radius**2))  #For circular elements
        return [thickness, mass, uvel_out, vvel_out, ang_vel_out]


def check_if_it_is_in_domain(x_val,y_val,x_min,x_max,y_min,y_max,R_earth,lat_init,adjust_lat_ref,dx,dy):
        point_is_in_domain=True

        if (x_val >= (x_max-x_min+(dx))) or (x_val<= 0) or  (y_val >= (y_max-y_min+(dy))) or (y_val <= 0):
                point_is_in_domain=False

        return point_is_in_domain

def check_if_it_is_ice(x_val,y_val,ice_mask,dx,icethres):
        i_val=int(floor(x_val/dx))
        j_val=int(floor(y_val/dx))
        #print('jval',j_val,'ival',i_val)


############# ##### ####### #########

        if ice_mask[j_val,i_val]>icethres: #0.5:#99:
                it_is_ice=True
        else:
                it_is_ice=False
        return it_is_ice

def calculate_element_area(element_type,Radius):
        if element_type=='square':
                element_area=(2*Radius)**2
        elif element_type=='hexagon':
                element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*(Radius)**2) #Area of hexagon around circle (used for packing)
                #Another derivation uses innner hexagon with two more triangles added, which is a 1/6 of the hexagon area each (two since there are 6, shared by 3 each)
                #element_area=(4./3.)*H_i, where H_i=(3.*np.sqrt(3.)/2.)*((Radius)**2)  is the area of the inner hexagon (with sides equal to the radius)

        return element_area

def add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,\
                iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice):
        tol=0.0000000000001
        if New_thickness[j,i]>1.+tol:
                print( 'The new thickness is too big!!!', i,j,  New_thickness[j,i])
                halt
        if ((Nbh_thickness_value-New_thickness[j,i])>tol):
                y_val=y_shift[j]
                x_val=x_shift[i]
                if Nbh_thickness_value>0. and (abs(Nbh_thickness_value-1)<tol) and (Nbh_mass_value>New_mass[j,i]):
                        berg_count=berg_count+1
                        dx_berg.append(x_val)
                        dy_berg.append(y_val)
                        iceberg_num.append(berg_count)
                        width_val=np.sqrt((Nbh_thickness_value-New_thickness[j,i])*grid_area)
                        if Nbh_thickness_value-New_thickness[j,i]<0:
                                print( 'Stop', Nbh_thickness_value, New_thickness[j,i],i,j)
                        mass_val=Nbh_mass_value-New_mass[j,i]
                        thickness_val=mass_val/(width_val*width_val*rho_ice)
                        thickness.append(thickness_val)
                        width.append(width_val)
                        mass.append(mass_val)
                        static_berg.append(1.)
        return [berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]


def add_extra_bergs_on_boundary(dx,dy,x,y,Number_of_bergs,element_area,rho_ice,X_min,X_max,Y_min,Y_max,R_earth,lat_init,adjust_lat_ref,ice_mask,\
                iceberg_num,width,dx_berg,dy_berg,h_ice,element_type,static_berg,thickness,mass):
        eps=0.000001
        grid_area=dx*dy
        Nx=len(x)  ; Ny=len(y)
        x_shift=(x-np.min(x))+(dx/2) ; y_shift=(y-np.min(y))+(dy/2)

        thickness_temp= [1. for i in range(Number_of_bergs)] ;
        mass_temp= [element_area*rho_ice*1. for i in range(Number_of_bergs)] ;
        New_area=regrid_iceberg_thickness(dy_berg,dx_berg,Number_of_bergs,thickness_temp,mass_temp,h_ice,x_shift,y_shift,rho_ice,element_type,static_berg,plot_outcome=False)
        New_thickness=(New_area)/(rho_ice*grid_area)  #Should be equal to one everywhere where there is ice shelf.
        New_mass=regrid_iceberg_thickness(dy_berg,dx_berg,Number_of_bergs,thickness,mass,h_ice,x_shift,y_shift,rho_ice,element_type,static_berg,plot_outcome=False)
        berg_count=Number_of_bergs


        #Vertical boundary
        for i in np.array([0,Nx-1]):
                #for j in range(Ny):
                for j in range(1,Ny-1):
                        Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
                        if i==0:
                                Nbh_thickness_value=New_thickness[j,1]  #Neighbouring shelf value
                                #x_val=0+(eps*dx)
                        if i==Nx-1:
                                Nbh_thickness_value=New_thickness[j,Nx-2]  #Neighbouring shelf value
                                #x_val=x_shift[i]+(dx/2)-eps
                        [berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
                                        Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)

        ##Horizontal boundary
        for j in np.array([0,Ny-1]):
                #for i in range(1,Nx-1):
                for i in range(1,Nx-1):
                        Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
                        if j==0:
                                Nbh_thickness_value=New_thickness[1,i]  #Neighbouring shelf value
                                #y_val=0+(eps*dy)
                        if j==Ny-1:
                                Nbh_thickness_value=New_thickness[Ny-2,i]  #Neighbouring shelf value
                                #y_val=y_shift[j]+(dy/2)-eps
                        #Nbh_thickness_value=np.mean(New_thickness[range(1,Ny-1),i])  #Neighbouring shelf value
                        [berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
                                        Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)

        #Doing the corners seperately
        for j in np.array([0,Ny-1]):
                #for i in range(1,Nx-1):
                for i in np.array([0,Nx-1]):
                        Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
                        if j==0 and i==0:
                                Nbh_thickness_value=New_thickness[1,1]  #Neighbouring shelf value
                        if j==0 and i==Nx-1:
                                Nbh_thickness_value=New_thickness[1,Nx-2]  #Neighbouring shelf value
                        if j==Ny-1 and i==0:
                                Nbh_thickness_value=New_thickness[Ny-2,1]  #Neighbouring shelf value
                        if j==Ny-1 and i==Nx-1:
                                Nbh_thickness_value=New_thickness[Ny-2,Nx-2]  #Neighbouring shelf value
                        [berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
                                        Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)


        Number_of_bergs=berg_count

        return [dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg,thickness,mass]


def remove_stationary_bergs(dx_berg, dy_berg,iceberg_num,width,Number_of_bergs,static_berg):
        dx_berg_new=[]
        dy_berg_new=[]
        iceberg_num_new=[]
        width_new=[]
        static_berg_new=[]
        for i in range(Number_of_bergs):
                if static_berg[i]<0.5:
                        dx_berg_new.append(dx_berg[i])
                        dy_berg_new.append(dy_berg[i])
                        iceberg_num_new.append(iceberg_num[i])
                        width_new.append(width[i])
                        static_berg_new.append(static_berg[i])

        Number_of_bergs_new=len(iceberg_num_new)
        return [dx_berg_new, dy_berg_new,iceberg_num_new, width_new,Number_of_bergs_new,static_berg_new]

def Create_icebergs(lon_init,lat_init,Radius,R_earth, x, y,ice_mask,h_ice,A68,uvel_in,vvel_in,Convert_to_lat_lon,rho_ice,\
                element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,Interpolate_from_four_corners,\
                Fill_in_the_boundaries,set_all_bergs_static_by_default,break_some_bonds,Remove_stationary_bergs,set_some_bergs_static_by_default,icethres):
        print( 'Starting to create icebergs...')
        dx_berg=[]  #x distance in cartesian of berg from lon_init
        dy_berg=[]  #y distance in cartesian of berg from lat_init
        dlon_berg=[] #x distance in lon of berg from lon_init
        #dlat_berg=[] #y distance in lat of berg from lat_init
        lon=[] #Longitude of iceberg
        lat=[] #Latitude of iceberg
        iceberg_num=[] #ID of iceberg

        element_area=calculate_element_area(element_type,Radius)
        #width=np.sqrt(element_area)
        width=[]

        #dx=x[1]-x[0]  ;dy=y[1]-y[0]
        #X_min=np.min(x)-dx; X_max=np.max(x)+dx  #In lat lon or cartesian
        #Y_min=np.min(y)-dy; Y_max=np.max(y)+dy  #In lat lon or cartesian
        X_min=np.min(x); X_max=np.max(x)  #In lat lon or cartesian
        Y_min=np.min(y); Y_max=np.max(y)  #In lat lon or cartesian

        #N=2*int(ceil((R_max)/(2*Radius))+2)

        if element_type=='square':
                N=int(ceil((X_max-X_min)/(Radius)))
                M=int(ceil((Y_max-Y_min)/(Radius)))
        else:
                N=2*int(ceil((X_max-X_min)/(Radius)))
                M=2*int(ceil((Y_max-Y_min)/(Radius)))

        #N=10
        #M=6
        #MP4
        dx=x[1]-x[0]
        dy=y[1]-y[0]
        Lx=X_max+dx
        Ly=Y_max+dx

        berg_count=0
        #for j in range(N):
        for i in range(N):
                if element_type=='square':
                        x_start=Radius#+Radius/2
                        y_start=Radius#+Radius/2
                        x_val=x_start+(2*i*Radius)

                #Hexagonal
                else:
                        y_start=(Radius)+(((i%2)*Radius))
                        x_start=((2/sqrt(3))*Radius)
                        x_val=x_start + (np.sqrt(3)*Radius*i)

                for j in range(M):
                #for i in range(M):
                        #x_val=x_start+(2*i*Radius)  ; y_val=y_start
                        y_val=y_start+(2*j*Radius)  ;# x_val=x_start + (np.sqrt(3)*Radius*i)
                        if check_if_it_is_in_domain(x_val,y_val,X_min,X_max,Y_min,Y_max,R_earth,lat_init,adjust_lat_ref,dx,dy):
                        #if True:
                                #R_val=np.sqrt(((x_val-x0)**2) +((y_val-y0)**2))
                                if check_if_it_is_ice(x_val,y_val,ice_mask,dx,icethres):
                                        #Don't allow points closer than R from the boundary (these are sorted out later)
                                        #if abs(y_val-Ly)<Radius or  ((abs(x_val-Lx)<((2/sqrt(3))*Radius)) and element_type=='hexagon') or ((abs(x_val-Lx)<Radius) and element_type=='square'):
                                        berg_count=berg_count+1
                                        iceberg_num.append(berg_count)
                                        dx_berg.append(x_val)
                                        dy_berg.append(y_val)
                                        width.append(np.sqrt(element_area))
        #print 'dx_berg',dx_berg
        #print 'dy_berg',dy_berg



        Number_of_bergs=berg_count
        print( 'Icebergs created. Number of bergs = ', Number_of_bergs)

        #Deciding if icebergs are static or not
        static_berg = [0. for i in dx_berg]
        if set_all_bergs_static_by_default==True:
                static_berg = [1. for i in static_berg]
        if set_some_bergs_static_by_default is True:
                for i in range(len(dx_berg)):
                        if dx_berg[i]>Lx/2:
                                static_berg[i]=1.
        if break_some_bonds==True:
                #static_berg =Change_static_berg_after_calving(Number_of_bergs,lat,lon, static_berg)
                static_berg =Change_static_berg_after_calving(Number_of_bergs,dy_berg,dx_berg, static_berg)
                static_berg=static_berg[0]

        if Remove_stationary_bergs is True:
                Fill_in_the_boundaries=False
                [dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg]=remove_stationary_bergs(dx_berg, dy_berg,iceberg_num,\
                                width,Number_of_bergs,static_berg)
                print( 'After purge  Number of bergs = ', Number_of_bergs)

        #Defining the thickness of the icebergs
        (thickness, mass, uvel_out, vvel_out,ang_vel_out)=Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,rho_ice,Radius, h_ice,A68,uvel_in,vvel_in,x,y,\
                        width,Interpolate_from_four_corners,element_area,element_type,static_berg)

        N_bergs_before_bd=Number_of_bergs
        if Fill_in_the_boundaries==True and element_type=='hexagon':
                [dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg,thickness,mass]= add_extra_bergs_on_boundary(dx,dy,x,y,Number_of_bergs,element_area,rho_ice,X_min,X_max,Y_min,\
                                Y_max,R_earth,lat_init,adjust_lat_ref,ice_mask,iceberg_num,width,dx_berg,dy_berg,\
                                h_ice,element_type,static_berg,thickness,mass)
                print( 'Number of icebergs after accounting for boundaries = ', Number_of_bergs)



        if Convert_to_lat_lon==True:
                #Defining lon lat positions:
                #dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
                #dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
                for i in range(Number_of_bergs):
                        #Finding latittude
                        dlat_dy=(180/np.pi)*(1/R_earth)

                        lat.append(lat_init+(dlat_dy*dy_berg[i]))

                        #Finding longitude
                        if adjust_lat_ref==True:
                                lat_ref=lat_init+(dlat_dy*dy_berg[i])
                                #print lat_ref
                        dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
                        lon.append(lon_init+(dlon_dx*dx_berg[i] ))

                        #dlon_berg.append(dlon_dx*dx_berg[i])
                        #lon.append(lon_init+dlon_berg[i])
        else:
                if scale_the_grid_to_lat_lon==True:
                        Scale_up=1./2000.
                        Radius=Radius*Scale_up
                        dx_berg = [(i*Scale_up) for i in dx_berg] ; dy_berg = [i*Scale_up for i in dy_berg]
                        x=x*Scale_up ; y=y*Scale_up
                        dx=dx*Scale_up ;dy=dy*Scale_up
                x=(x-np.min(x))+(dx/2) ; y=(y-np.min(y))+(dy/2)
                lon=dx_berg  ; lat=dy_berg

        #Note that static_berg calculations used to be here, after the conversion. I have moved them. I hope that this does not affect answers.

        return (Number_of_bergs,lon,lat,iceberg_num,dx_berg,dy_berg,thickness,ang_vel_out,uvel_out,vvel_out,mass,width,x,y,Radius,static_berg,N_bergs_before_bd)



def Create_calving_event(lat1,lon1,lat2,lon2):
        [R_calve, Calve_lon, Calve_lat]=get_calving_parameters()

        #Circular calving from ice front
        R1=np.sqrt((lon1-Calve_lon)**2+ (lat1-Calve_lat)**2)
        R2=np.sqrt((lon2-Calve_lon)**2+ (lat2-Calve_lat)**2)

        #Iceberg splitting down the midde (length wise)
        #R1=lat1
        #R2=lat2
        #R_calve=500*1000
        #Iceberg splitting down the midde (width wise)
        R1=lon1
        R2=lon2
        R_calve=250*1000

        bond_broken=False
        if ((R1 < R_calve)*(R2>R_calve) )  > 0.5 :
                bond_broken=True
        if ((R2 < R_calve)*(R1>R_calve) )  > 0.5 :
                bond_broken=True

        return bond_broken

def get_calving_parameters():
        #Calve_lat =20
        #Calve_lon=160
        #R_calve=15
        Calve_lat =20.2*2000
        Calve_lon=165*2000
        R_calve=12*2000
        #R_calve=15.*2000.
        #R_calve=1.*2000.

        return [R_calve, Calve_lon, Calve_lat]

def Change_static_berg_after_calving(Number_of_bergs,lat,lon, static_berg):
        [R_calve, Calve_lon, Calve_lat]=get_calving_parameters()
        count=0.
        for i in range(Number_of_bergs):
                R1=np.sqrt((lon[i]-Calve_lon)**2+ (lat[i]-Calve_lat)**2)
                #Making calved icebergs not be static.
                if R1<R_calve:
                        count=count+1.
                        #print 'An iceberg is now static!'
                        static_berg[i]=0.
        print( 'Amount of unstatic icebergs after calving = ', count)

        return [static_berg]

def Change_static_berg_after_bonding(Number_of_bergs,lat,lon, static_berg):
        #For now, it just uses a bigger circle.
        multiplier=1.2
        [R_calve, Calve_lon, Calve_lat]=get_calving_parameters()
        count=0.
        for i in range(Number_of_bergs):
                R1=np.sqrt((lon[i]-Calve_lon)**2+ (lat[i]-Calve_lat)**2)
                #Making calved icebergs not be static.
                if R1<(R_calve*multiplier):
                        count=count+1.
                        #print 'An iceberg is now static!'
                        static_berg[i]=0.
        print ('Amount of unstatic icebergs after bonding and restatic= ', count)

        return [static_berg]

def find_max_number_of_bonds(first_berg_num,Number_of_bonds):
        count=1
        best_count=1
        sorted_bonds=np.sort(first_berg_num)
        for k in range(1,Number_of_bonds):
                previous_berg=sorted_bonds[k-1]
                #print k, sorted_bonds[k]
                if sorted_bonds[k] == sorted_bonds[k-1]:
                        count=count+1
                        best_count=max(count,best_count)
                else:
                        count=1
        Max_number_of_bonds=best_count
        return Max_number_of_bonds

def Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds,static_berg, Allow_bonds_for_static_iceberg\
                ,N_bergs_before_bd,Allow_bonds_with_boundary_bergs):
        print( 'Starting to create bonds...')
        #Defining Bonds:
        Bond=np.zeros((Number_of_bergs, Number_of_bergs))
        bond_broken=False
        first_berg_num=[]  # Initializing bond list first berg
        first_berg_ine=[]  # Initializing bond list first berg
        first_berg_jne=[]  # Initializing bond list first berg
        first_berg_lat=[]  # Initializing bond list first berg
        first_berg_lon=[]  # Initializing bond list first berg
        other_berg_num=[]  # Initializing bond list other berg
        other_berg_ine=[]  # Initializing bond list other berg
        other_berg_jne=[]  # Initializing bond list other berg
        other_berg_lat=[]  # Initializing bond list other berg
        other_berg_lon=[]  # Initializing bond list other berg
        bond_count=0
        for i in range(Number_of_bergs):
                if (static_berg[i]<0.5) or (Allow_bonds_for_static_iceberg is True):
                        #for j in range(Number_of_bergs):
                        for j in range(i):
                                if (static_berg[j]<0.5) or (Allow_bonds_for_static_iceberg is True):
                                        if i!=j:
                                                if ((i <N_bergs_before_bd) and (j <N_bergs_before_bd))  or (Allow_bonds_with_boundary_bergs is True):
                                                        R_dist=np.sqrt(((dx_berg[i]-dx_berg[j])**2) + ((dy_berg[i]-dy_berg[j])**2))
                                                        if break_some_bonds==True:
                                                                bond_broken=Create_calving_event(lat[i],lon[i],lat[j],lon[j])
                                                        if R_dist < (2.01*Radius) and (bond_broken==False):
                                                                bond_count=bond_count+2
                                                                #Connect bond in the first direction
                                                                first_berg_num.append(iceberg_num[i]);  first_berg_ine.append(999);
                                                                other_berg_num.append(iceberg_num[j]);  other_berg_ine.append(999);
                                                                first_berg_jne.append(999);first_berg_lat.append(lat[i]); first_berg_lon.append(lon[i])
                                                                other_berg_jne.append(999);other_berg_lat.append(lat[j]); other_berg_lon.append(lon[j])
                                                                #Connect bond in the other direction
                                                                first_berg_num.append(iceberg_num[j]);  first_berg_ine.append(999);
                                                                other_berg_num.append(iceberg_num[i]);  other_berg_ine.append(999);
                                                                first_berg_jne.append(999);first_berg_lat.append(lat[j]); first_berg_lon.append(lon[j])
                                                                other_berg_jne.append(999);other_berg_lat.append(lat[i]); other_berg_lon.append(lon[i])
        Number_of_bonds=bond_count
        Max_number_of_bonds=find_max_number_of_bonds(first_berg_num,Number_of_bonds)
        print( 'Number of bonds created = ' , Number_of_bonds)
        print( 'Maximum number of bonds = ' ,  Max_number_of_bonds)

        return [ Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon]

def load_ISOMIP_ice_geometry(filename,buffer_number):
        with nc.Dataset(filename) as file:
                ocean_mask = file.variables['openOceanMask'][:,:]
                upperSurface = file.variables['upperSurface'][:,:]
                lowerSurface = file.variables['lowerSurface'][:,:]
                x = file.variables['x'][:]
                y = file.variables['y'][:]

        ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
        h_ice=upperSurface-lowerSurface #The ice thickness
        M=ice_mask.shape

        #Setting the boundaries to non-ice
        A=np.arange(0,buffer_number)
        B=np.arange(M[0]-buffer_number,M[0])
        C=np.arange(M[1]-buffer_number,M[1])
        ice_mask[A,:]=0; ice_mask[B,:]=0
        ice_mask[:,A]=0; ice_mask[:,C]=0

        return [x,y,ice_mask,h_ice]


def create_clipped_icethickness_file(h_ice,area,mass,grid_area,gravity,Ice_geometry_source):
        #Creating clipped file
        [ny, nx]= h_ice.shape ;

        New_ice_thickness_filename='output_files/'+ Ice_geometry_source +'_Ice_Shelf_clipped.nc'
        g=Dataset(New_ice_thickness_filename,'w') # w if for creating a file

        g.createDimension('nx',nx)
        g.createDimension('ny',ny)

        thick_h=g.createVariable('thick','f8',('ny','nx'))
        area_h=g.createVariable('area','f8',('ny','nx'))
        p_surf_h=g.createVariable('p_surf','f8',('ny','nx'))

        thick_h.units = 'm'
        thick_h.standard_name = 'ice shelf thickness (clipped)'
        area_h.units = 'm2'
        area_h.standard_name = 'ice shelf area'
        p_surf_h.units = 'Pa'
        p_surf_h.standard_name = 'surface pressure due to ice shelf'

        p_surf=(gravity*mass)/grid_area
        g.variables['thick'][:]=h_ice
        g.variables['area'][:]=area
        g.variables['p_surf'][:]=p_surf
        print ('Creating clipped ice file: ' , New_ice_thickness_filename)

        g.sync()
        g.close()


def load_ISOMIP_reduced_ice_geometry(ice_filename,buffer_number):
        with nc.Dataset(ice_filename) as file:
                h_ice = file.variables['thick'][:,:]
                area =  file.variables['area'][:,:]
                M=h_ice.shape

        y=np.linspace(1000,79000,M[0],endpoint=True)
        #x=np.linspace(321000,799000,M[1],endpoint=True)
        x=np.linspace(1000,479000,M[1],endpoint=True)
        ice_mask=h_ice>0.


        count=0. #MP1
        #Setting the boundaries to non-ice
        A=np.arange(0,buffer_number)
        B=np.arange(M[0]-buffer_number,M[0])
        C=np.arange(M[1]-buffer_number,M[1])
        #if buffer_number>0:
        ice_mask[A,:]=0; ice_mask[B,:]=0
        ice_mask[:,A]=0; ice_mask[:,C]=0

        return [x,y,ice_mask,h_ice]

def load_generic_ice_geometry(filename,y_disp):
        print( 'Loading from ' , filename)
        with nc.Dataset(filename) as file:
                #h_ice = file.variables['thick'][:]
                #area = file.variables['area'][:,:]
                h_ice = file.variables['b_id'][:]
                x = file.variables['lon'][:]
                y = file.variables['lat'][:]+y_disp

                ice_mask=h_ice>0
                ice_mask[np.where(ice_mask<0)]=0.
                #ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
                h_ice[np.where(np.isnan(h_ice))]=0.
                ice_mask[np.where(np.isnan(ice_mask))]=0.
                dx=x[1]-x[0]
                #lon_init=x[0]-(dx/2.)+360
                #lat_init=y[0]-(dx/2.)
                lon_init=x[0]+360
                lat_init=y[0]

                #dy=y[1]-y[0]
                lon_init=x[0]-(dx/2.)+360
                lat_init=y[0]-(dx/2.)
                print( 'Bottom corner = [ ', lon_init,' , ' , lat_init, ' dx = ' , dx)

        return [x,y,ice_mask,h_ice,lon_init, lat_init]

def load_Weddel_ice_geometry(filename):
        with nc.Dataset(filename) as file:
                h_ice = file.variables['thickness'][:,:]
                ice_mask = file.variables['icemask_grounded_and_shelves'][:,:]
                x = file.variables['lon'][:]
                y = file.variables['lat'][:]

                ice_mask[np.where(ice_mask<-1)]=0.
                #ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
                h_ice[np.where(np.isnan(h_ice))]=0.
                ice_mask[np.where(np.isnan(ice_mask))]=0.

        return [x,y,ice_mask,h_ice]

def Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num,static_berg):
        print ('You have chosen to choose just one icebergs!!!')
        for k in range(len(lat)):
                if iceberg_num[k]==chosen_berg_num:
                        berg_ind=k
        lon_temp=lon[berg_ind]; lon=[] ; lon.append(lon_temp)
        lat_temp=lat[berg_ind]; lat=[] ; lat.append(lat_temp)
        thickness_temp=thickness[berg_ind]; thickness=[] ; thickness.append(thickness_temp)
        mass_temp=mass[berg_ind]; mass=[] ; mass.append(mass_temp)
        width_temp=width[berg_ind]; width=[] ; width.append(width_temp)
        static_berg_temp=static_berg[berg_ind]; static_berg=[] ; static_berg.append(static_berg_temp)
        iceberg_num_temp=iceberg_num[berg_ind]; iceberg_num=[] ; iceberg_num.append(iceberg_num_temp)
        Number_of_bergs=1
        return [Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num]

def plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon, \
                plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions,static_berg,h_ice_new,plot_h_ice_new):
        print ('Starting to plot...')
        Radius=Radius*IA_scaling
        circ_ind=np.linspace(0,2*pi,100);

        if plot_ice_mask==True:
                print(shape(x),shape(y),shape(ice_mask))
                cNorm = mpl.colors.Normalize(vmin=0., vmax=1.)
                plt.pcolor(x,y,ice_mask,norm=cNorm)
                plt.title('Ice mask')

        elif plot_h_ice_new is True: # Ice thickness after interpolation
                cNorm = mpl.colors.Normalize(vmin=np.min(h_ice_new), vmax=np.max(h_ice_new))
                plt.pcolor(x,y,h_ice_new , norm=cNorm)
                plt.title('New thickness')

        elif plot_ice_thickness==True:
                cNorm = mpl.colors.Normalize(vmin=0., vmax=1000.)
                plt.pcolor(x,y,h_ice,cmap='jet',norm=cNorm)
                plt.title('Old thickness')

        if plot_icebergs_positions==True:
                variable=thickness
                #variable=static_berg
                plt.scatter(lon, lat,color='yellow')
                cNorm = mpl.colors.Normalize(vmin=0., vmax=np.max(variable))
                plt.scatter(lon, lat,c=variable,norm=cNorm,cmap='jet',s=150)
                #cNorm = mpl.colors.Normalize(vmin=-1, vmax=1.)
                #plt.scatter(lon, lat,c=static_berg,norm=cNorm,cmap='jet',s=150)
                #plt.title('Iceberg positions')

        #plt.plot(lon, lat,'bo-',linewidth=5)
        if plot_circles==True:
                for k in range(Number_of_bergs):
                        if Convert_to_lat_lon==True:
                                dR_lat=(Radius/R_earth)*(180/np.pi)
                                dR_lon=(Radius/R_earth)*(180/np.pi)*(1/np.cos(lat[k]*np.pi/180))
                                plt.plot(lon[k]+(dR_lon*cos(circ_ind)),lat[k]+(dR_lat*sin(circ_ind)),'c');
                        else:
                                plt.plot(lon[k]+(Radius*cos(circ_ind)),lat[k]+(Radius*sin(circ_ind)),'c');



        #plt.plot(lon, lat,'bo-')
        if Convert_to_lat_lon==True:
                plt.xlabel('longitude (deg)')
                plt.ylabel('latitude (deg)')
        else:
                plt.xlabel('longitude (m)')
                plt.ylabel('latitude (m)')
        plt.colorbar()
        plt.grid(True)


def plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds):

        for k in range(Number_of_bonds):
                x_bond=[]
                y_bond=[]
                x_bond.append(first_berg_lon[k])
                x_bond.append(other_berg_lon[k])
                y_bond.append(first_berg_lat[k])
                y_bond.append(other_berg_lat[k])
                plt.plot(x_bond, y_bond,'r',linewidth=5)


def spread_mass_to_ocean(Nx,Ny,i,j,mass_on_ocean,x,y,Area,Mass,element_type,grid_area,static_berg):
                #Note that the x,y coming into this routine are the position within a cell (from 0 to 1), with 0.5,0.5 being in the center of the cell.

                #Initialize weights for each cell
                yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0.
                yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

                if element_type=='square':
                #if True:
                        L = min(( (np.sqrt(Area)  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim length of the square)
                        xL=min(0.5, max(0., 0.5-(x/L)))
                        xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
                        xC=max(0., 1.-(xL+xR))
                        yD=min(0.5, max(0., 0.5-(y/L)))
                        yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
                        yC=max(0., 1.-(yD+yU))

                        yDxL=yD*xL#*grd%msk[i-1,j-1]
                        yDxC=yD*xC#*grd%msk[i  ,j-1]
                        yDxR=yD*xR#*grd%msk[i+1,j-1]
                        yCxL=yC*xL#*grd%msk[i-1,j  ]
                        yCxR=yC*xR#*grd%msk[i+1,j  ]
                        yUxL=yU*xL#*grd%msk[i-1,j+1]
                        yUxC=yU*xC#*grd%msk[i  ,j+1]
                        yUxR=yU*xR#*grd%msk(i+1,j+1]
                        yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

                if element_type=='hexagon':
                #if False:#element_type=='hexagon':
                        H = min(( (np.sqrt(Area/(2*sqrt(3)))  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim Apothen of the hexagon)
                        S=(2/np.sqrt(3))*H #Side of the hexagon
                        if S>0.5:
                                print( 'Elements must be smaller than a whole gridcell', 'i.e.: S= ' , S , '>=0.5',i,j)
                                halt

                        #Subtracting the position of the nearest corner from x,y
                        origin_x=1 ; origin_y=1
                        if x<0.5:
                                origin_x=0
                        if y<0.5:
                                origin_y=0

                        x0=(x-origin_x) #Position of the hexagon center, relative to origin at the nearest vertex
                        y0=(y-origin_y)

                        #(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
                        (Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,0.)
                        if min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <0:
                                print( 'Yolo')
                                print( x0,y0,H)
                                #print min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4))
                                print( Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

                        Area_Q1=Area_Q1/Area_hex
                        Area_Q2=Area_Q2/Area_hex
                        Area_Q3=Area_Q3/Area_hex
                        Area_Q4=Area_Q4/Area_hex


                        #Now, you decide which quadrant belongs to which mass on ocean cell.
                        #Top right vertex
                        if x>=0.5 and y>= 0.5:
                                yUxR=Area_Q1
                                yUxC=Area_Q2
                                yCxC=Area_Q3
                                yCxR=Area_Q4

                        #Top left vertex
                        if x<0.5 and y>= 0.5:
                                yUxC=Area_Q1
                                yUxL=Area_Q2
                                yCxL=Area_Q3
                                yCxC=Area_Q4

                        #Bottom left vertex
                        if x<0.5 and y< 0.5:
                                yCxC=Area_Q1
                                yCxL=Area_Q2
                                yDxL=Area_Q3
                                yDxC=Area_Q4

                        #Bottom right vertex
                        if x>=0.5 and y< 0.5:
                                yCxR=Area_Q1
                                yCxC=Area_Q2
                                yDxC=Area_Q3
                                yDxR=Area_Q4

                        #if Sector<4 and Sector>-1:
                        #       print Sector


                        #Check that this is true
                        if abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>0.001:
                                print ('All the mass is not being used!!!')
                                #print W1 , W2 , W3 , W4 , W5 , W6
                                #print Area_Upper, Area_Lower, Area_right, Area_left
                                print( 'Areas: ',Area_hex,Area_hex*Area_Q1, Area_hex*Area_Q2, Area_hex*Area_Q3, Area_hex*Area_Q4)
                                print( 'x0=',x0, 'y0=',y0, 'H=', H)
                                print ('Total area= ',(Area_Q1+Area_Q2+Area_Q3+Area_Q4))#, Sector


                #Accounting for masked points
                a=1. ; b=1. ; c=1. ; d=1.
                if i==0:
                        a=0.;
                if j==0:
                        b=0.;
                if i==Nx-1:
                        c=0;
                if j==Ny-1:
                        d=0;

                fraction_used= ((yDxL*a*b) + (yDxC*b) + (yDxR*b*c) +(yCxL*a) + (yCxR*c) + (yUxL*a*d) + (yUxC*d) + (yUxR*c*d) + (yCxC))
                if static_berg>0.5:
                        fraction_used=1.

                mass_on_ocean[i,j,1]=mass_on_ocean[i,j,1]+(a*b*yDxL*Mass/fraction_used)
                mass_on_ocean[i,j,2]=mass_on_ocean[i,j,2]+(b*yDxC*Mass/fraction_used)
                mass_on_ocean[i,j,3]=mass_on_ocean[i,j,3]+(b*c*yDxR*Mass/fraction_used)
                mass_on_ocean[i,j,4]=mass_on_ocean[i,j,4]+(a*yCxL*Mass/fraction_used)
                mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+(yCxC*Mass/fraction_used)
                mass_on_ocean[i,j,6]=mass_on_ocean[i,j,6]+(c*yCxR*Mass/fraction_used)
                mass_on_ocean[i,j,7]=mass_on_ocean[i,j,7]+(a*d*yUxL*Mass/fraction_used)
                mass_on_ocean[i,j,8]=mass_on_ocean[i,j,8]+(d*yUxC*Mass/fraction_used)
                mass_on_ocean[i,j,9]=mass_on_ocean[i,j,9]+(c*d*yUxR*Mass/fraction_used)


                return mass_on_ocean


def calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, Lx):

        Lx_2=Lx/2.

        alpha=x2-x1
        delta=y2-y1
        beta=x4-x1
        epsilon=y4-y1
        gamma=(x3-x1)-(alpha+beta)
        kappa=(y3-y1)-(delta+epsilon)

        a=(kappa*beta-gamma*epsilon)
        dx=np.mod(x-(x1-Lx_2),Lx)+(x1-Lx_2)-x1
        dy=y-y1
        b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
        c=(alpha*dy-delta*dx)

        #print 'alpha,beta,gamma',alpha,beta,gamma
        #print 'delta,epsilon,kappa',delta,epsilon,kappa
        #print 'dx,dy',dx,dy
        #print 'A,B,C,',a,b,c

        if (abs(a)>1.e-12):
                d=0.25*(b**2)-a*c
                if (d>=0.):
                        yy1=-(0.5*b+sqrt(d))/a
                        yy2=-(0.5*b-sqrt(d))/a
                        if (abs(yy1-0.5).lt.abs(yy2-0.5)):
                                yj=yy1;
                        else:
                                yj=yy2;
                else:
                        print ('ERROR', halt)
        else:
                if (b!=0.):
                        yj=-c/b
                else:
                        yj=0.

        a=(alpha+gamma*yj)
        b=(delta+kappa*yj)
        if (a!=0.):
                xi=(dx-beta*yj)/a
        elif (b!=0.):
                xi=(dy-epsilon*yj)/b
        else:
                c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
                if (c!=0.):
                        xi=(epsilon*dx-beta*dy)/c
                else:
                        print( 'ERROR', halt2)

        return xi,yj



def regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type,static_berg,plot_outcome=True):
        Nx=len(x)  ; Ny=len(y)
        dx=x[1]-x[0]    ;dy=y[1]-y[0]
        New_mass=h_ice*0 #Initializing regrided field to zero
        Area = [(mass[i]/(rho_ice*thickness[i])) for i in range(Number_of_bergs)] ;
        grid_area=dx*dy #Assuming a regular grid
        Orig_mass=h_ice*rho_ice*grid_area

        mass_on_ocean=np.zeros([Nx,Ny,10])   #Setting up matrix to spread mass to ocean.  Note that I have used 10 points  so that I can ignore 0 and match with python numbering
        #mass_on_ocean=np.zeros([Nx,Ny])   #Setting up matrix to spread mass to ocean.
        #Note: You may want to import the grid that the ocean model actually sees.
        for berg_count in range(Number_of_bergs):
                x_val=lon[berg_count]
                y_val=lat[berg_count]
                Area_val=Area[berg_count]
                mass_val=mass[berg_count]

                j_val=int(floor(y_val/dy))
                i_val=int(floor(x_val/dx))
                x_cell=(x_val-x[i_val])/sqrt(grid_area)+0.5
                y_cell=(y_val-y[j_val])/sqrt(grid_area)+0.5

                #if plot_outcome==True:
                #       print 'YOLO0', x_val-123626.65521,y_val
                #       print 'YOLO1', x[i_val],y[j_val]
                #       print 'YOLO2', x_cell,y_cell

                #Using the same interpolation as in the model
                #x1=x[i_val]-(dx/2.)   ; y1=y[j_val]-(dx/2.)
                #x2=x[i_val]+(dx/2.) ; y2=y[j_val]-(dx/2.)
                #x3=x[i_val]+(dx/2.) ; y3=y[j_val]+(dx/2.)
                #x4=x[i_val]-(dx/2.) ; y4=y[j_val]+(dx/2.)
                #[xi,yj]=calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x_val, y_val, Lx=100000000000000.)
                #[x_cell,y_cell]=calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x_val, y_val, Lx=100000000000000.)
                #print 'YOLO:x1,x2,x3,x4',x1,x2,x3,x4
                #print 'YOLO:y1,y2,y3,y4',y1,y2,y3,y4
                #print 'YOLO:xi,yj,x_cell,y_cell', xi,yj, x_cell,y_cell

                mass_on_ocean=spread_mass_to_ocean(Nx,Ny,i_val,j_val,mass_on_ocean,x_cell,y_cell,Area_val,mass_val,element_type,grid_area,static_berg[berg_count])

        #Adding mass_onto_ocean
        for i in range(1,Nx-1):
                for j in range(1,Ny-1):
                        New_mass[j,i]=mass_on_ocean[i,j,5]  \
                                        +  ( ( (mass_on_ocean[i-1,j-1,9]+mass_on_ocean[i+1,j+1,1])   \
                                        +  (mass_on_ocean[i+1,j-1,7]+mass_on_ocean[i-1,j+1,3]) ) \
                                        +   ( (mass_on_ocean[i-1,j  ,6]+mass_on_ocean[i+1,j  ,4])   \
                                        +  (mass_on_ocean[i  ,j-1,8]+mass_on_ocean[i  ,j+1,2]) ) )

        #Adding mass to boundary cells
        for j in range(1,Ny-1):
                i=0
                New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i+1,j+1,1] +  mass_on_ocean[i+1,j-1,7] \
                        +  mass_on_ocean[i+1,j  ,4]     +  mass_on_ocean[i  ,j-1,8] + mass_on_ocean[i  ,j+1,2]
                i=Nx-1
                New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i-1,j-1,9] + mass_on_ocean[i-1,j+1,3] \
                                        +   mass_on_ocean[i-1,j ,6] + mass_on_ocean[i ,j-1,8]+mass_on_ocean[i ,j+1,2]

        for i in range(1,Nx-1):
                j=0
                New_mass[j,i]=mass_on_ocean[i,j,5] +   mass_on_ocean[i+1,j+1,1] +  mass_on_ocean[i-1,j+1,3]  \
                                +  mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i+1,j  ,4] +  mass_on_ocean[i  ,j+1,2]

                j=Ny-1
                New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i-1,j-1,9]  +  mass_on_ocean[i+1,j-1,7] + \
                                mass_on_ocean[i-1,j  ,6]+mass_on_ocean[i+1,j  ,4] +  mass_on_ocean[i  ,j-1,8]

        #Corners
        i=0 ;j=0
        New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i+1,j+1,1] + mass_on_ocean[i+1,j  ,4] + mass_on_ocean[i  ,j+1,2]
        i=0 ;j=Ny-1
        New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i+1,j-1,7] + mass_on_ocean[i+1,j  ,4] + mass_on_ocean[i  ,j-1,8]
        i=Nx-1; j=0
        New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i-1,j+1,3] + mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i  ,j+1,2]
        i=Nx-1; j=Ny-1
        New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i-1,j-1,9] + mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i  ,j-1,8]


        if plot_outcome==True: #MP2
                #print 'Warning: the interpolation scheme is not the same as the one used in the model'

                plot_data=(Orig_mass-New_mass)/Orig_mass
                cNorm = mpl.colors.Normalize(vmin=-.1, vmax=.1)
                vmax=np.max(abs(plot_data))
                plt.pcolormesh(x,y,plot_data,cmap='bwr',norm=cNorm)
                plt.title('Error between bergs and shelf')

                #plot_data=(New_mass/New_thickness/rho_ice)
                #vmax=np.max(abs(plot_data))
                #cNorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
                #plt.pcolormesh(x,y,plot_data,cmap='jet',norm=cNorm)
                #plt.title('h_ice')

                #plot_data=(New_mass-Orig_mass)/(rho_ice*grid_area)
                #cNorm = mpl.colors.Normalize(vmin=-0.00000000001, vmax=0.00000000001)
                #vmax=np.max(abs(plot_data))
                #plt.pcolormesh(x,y,plot_data,cmap='bwr',norm=cNorm)

                plt.colorbar()
                #plt.xlim([min(x)-,max(x)])
                #plt.ylim([min(y),max(y)])
                #plt.show()

        return New_mass

def Switch_x_and_y_directions(lat,lon, A, B, C):
        print( 'Switching x and y directions!!!')
        lat_new=[]
        lon_new=[]
        for i in range(len(lat)):
                lat_new.append(lon[i])
                lon_new.append(lat[i])
        M=A.shape
        A_new=np.zeros([M[1],M[0]])
        B_new=np.zeros([M[1],M[0]])
        C_new=np.zeros([M[1],M[0]])
        for i in range(M[0]):
                for j in range(M[1]):
                        A_new[j,i]=A[i,j]
                        B_new[j,i]=B[i,j]
                        C_new[j,i]=C[i,j]

        return  [lat_new,lon_new,A_new,B_new, C_new]

def interpolate_field_onto_new_grid(x,y,x_new,y_new,field):
        xx, yy = np.meshgrid(x, y)
        #f = interpolate.interp2d(x, y, field, kind='cubic')
        f = interpolate.interp2d(x, y, field, kind='linear')
        field_new = f(x_new, y_new)
        return field_new

def convert_input_to_catesian_coordinates(lon,lat,ice_mask,h_ice,R_earth,dx,lat_ref):
        #Parameters
        thickness_threshhold=0.01

        lon_max=np.max(lon) ; lon_min=np.min(lon)
        lat_max=np.max(lat) ; lat_min=np.min(lat)
        #lat_ref=np.mean(lat)

        #55.8 one more element of support near finger
        #lat_ref=-55.75
        print('lat_ref',lat_ref)
        #llon,llat=np.meshgrid(lon,lat)
        #lat_ref=np.mean(llat[np.where(h_ice>0.01)])
        #print('lat_ref new',lat_ref)

        #Converting to cartesian grid
        x=(np.pi/180)*R_earth*np.cos(lat_ref*np.pi/180.)*lon
        y=(np.pi/180)*R_earth*lat

        # x=R_earth*np.cos(lat*np.pi/180)*np.cos(lon*np.pi/180)
        # y=R_earth*np.cos(lon*np.pi/180)*np.cos(lat*np.pi/180)

        print('minx',np.min(x))
        print('maxx',np.max(x))
        print('miny',np.min(y))
        print('maxy',np.max(y))
        #x=x-np.min(x)
        #y=y-np.min(y)

        #Creating regular grid
        x_new=np.linspace(np.min(x), np.max(x), num=len(x))
        y_new=np.linspace(np.min(y),np.max(y), num=len(y))
        x_new=np.arange(np.min(x), np.max(x), step=dx)
        y_new=np.arange(np.min(y), np.max(y), step=dx)

        h_ice_new=interpolate_field_onto_new_grid(x,y,x_new,y_new,h_ice)
        h_ice_new[np.where(h_ice_new<thickness_threshhold)]=0.
        ice_mask_new=interpolate_field_onto_new_grid(x,y,x_new,y_new,ice_mask)
        #ice_mask_new=np.round(ice_mask_new)
        #print ice_mask_new[50,:]

        return [x_new,y_new,ice_mask_new,h_ice_new]


####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main(args):

        #Flags
        save_restart_files=args.save_restart_files
        save_new_h_ice_file=args.save_new_h_ice_file
        Convert_to_lat_lon=args.Convert_to_lat_lon
        input_is_cartesian=args.input_is_cartesian
        A68=args.A68

        #Iceberg setup flags
        only_choose_one_berg=args.only_choose_one_berg  ; chosen_berg_num=args.chosen_berg_num
        scale_the_grid_to_lat_lon=args.scale_the_grid_to_lat_lon
        adjust_lat_ref=args.adjust_lat_ref
        set_all_thicknesses_to_one=args.set_all_thicknesses_to_one   ; Th_prescribed=args.Th_prescribed
        uvel_in=args.uvel; vvel_in=args.vvel
        Interpolate_from_four_corners=args.Interpolate_from_four_corners
        Switch_x_and_y_to_rotate_90=args.Switch_x_and_y_to_rotate_90
        set_all_domain_to_ice=args.set_all_domain_to_ice
        Use_default_radius=args.Use_default_radius

        #Static vs non-static flags
        set_all_bergs_static_by_default=args.set_all_bergs_static_by_default
        set_some_bergs_static_by_default=args.set_some_bergs_static_by_default
        Make_icebergs_non_static_later=args.Make_icebergs_non_static_later

        #Mass spreadng flags
        element_type=args.element_type
        Switch_regridding_element_type=args.Switch_regridding_element_type
        Fill_in_the_boundaries=args.Fill_in_the_boundaries
        regrid_icebergs_onto_grid=args.regrid_icebergs_onto_grid

        #Plotting flags
        Run_plotting_subroutine=args.Run_plotting_subroutine
        plot_circles=args.plot_circles
        plot_ice_mask=args.plot_ice_mask
        plot_ice_thickness=args.plot_ice_thickness
        plot_icebergs_positions=args.plot_icebergs_positions
        plot_h_ice_new=args.plot_h_ice_new
        plot_bonds=args.plot_bonds

        #Bond related flags
        Create_icebergs_bonds=args.Create_icebergs_bonds
        break_some_bonds=args.break_some_bonds
        Remove_stationary_bergs=args.Remove_stationary_bergs
        Allow_bonds_for_static_iceberg=args.Allow_bonds_for_static_iceberg
        Allow_bonds_with_boundary_bergs=args.Allow_bonds_with_boundary_bergs


        #Which experiment
        Ice_geometry_source=args.Ice_geometry_source
        ISOMIP_reduced=args.ISOMIP_reduced
        ISOMIP_ice_geometry_filename=args.ISOMIP_ice_geometry_filename
        ISOMIP_reduced_ice_geometry_filename=args.ISOMIP_reduced_ice_geometry_filename
        Weddell_ice_geometry_filename=args.Weddell_ice_geometry_filename
        Generic_ice_geometry_filename=args.Generic_ice_geometry_filename


        #Parameters
        Radius=args.Radius
        rho_ice=args.rho_ice
        gravity=args.gravity
        mass_scaling=args.mass_scaling
        R_earth=args.R_earth
        print('R_earth',R_earth)
        buffer_number=args.buffer_number
        ang_vel_in=args.ang_vel_prescribed
        y_disp=args.displace_y
        lat_ref=args.lat_ref_grid_to_cartesian
        dx_in=args.dx
        icethres=args.icethres

        #Let interactive radius be different from radius for testing the model:
        IA_scaling=args.IA_scaling

        #Applying some special case flags
        if Ice_geometry_source=='ISOMIP':
                Convert_to_lat_lon=False       ; input_is_cartesian=True ;
        if Ice_geometry_source=='Generic':
                Convert_to_lat_lon=True       ; input_is_cartesian=False ; set_all_bergs_static_by_default=False ; Fill_in_the_boundaries=False ; Use_default_radius=False
        if Ice_geometry_source=='Weddell':
                Convert_to_lat_lon=False        ; input_is_cartesian=False

        #Applying Warnings
        if element_type=='square' and Interpolate_from_four_corners==True:
                print ('Square packing with R dividing dx, works best with interpolation off.')

        print( 'Element type= ', element_type, ';  Switched_regridding= ', Switch_regridding_element_type)
        if break_some_bonds==True:
                print ('Some bonds are being broken')

        #####################################################################################
        #####################################################################################
        Radius=Radius/IA_scaling

        if  Ice_geometry_source=='ISOMIP':
                if ISOMIP_reduced==False:
                        lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
                        (x,y,ice_mask,h_ice)=load_ISOMIP_ice_geometry(ISOMIP_ice_geometry_filename,buffer_number)
                        ice_filename=ISOMIP_ice_geometry_filename
                else:
                        lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
                        (x,y,ice_mask,h_ice)=load_ISOMIP_reduced_ice_geometry(ISOMIP_reduced_ice_geometry_filename,buffer_number)
                        ice_filename=ISOMIP_reduced_ice_geometry_filename

        if  Ice_geometry_source=='Weddell':
                lon_init=-32.9  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
                (x,y,ice_mask,h_ice)=load_Weddel_ice_geometry(Weddell_ice_geometry_filename)
                ice_filename=Weddell_ice_geometry_filename
                Radius=0.5*10000.

        if  Ice_geometry_source=='Generic':
                (x,y,ice_mask,h_ice, lon_init, lat_init)=load_generic_ice_geometry(Generic_ice_geometry_filename,y_disp)

        if input_is_cartesian is False:
                dx=dx_in #3500. #250. #1000.0 #2000.0 #4000.#1000.
                if (lat_ref==-9999):
                        lat_ref=np.mean(y)
                lat_grid=y.copy() ; lon_grid=x.copy()+360; ice_mask_grid=h_ice.copy()
                #print('ice_mask_grid',ice_mask_grid)
                (x,y,ice_mask,h_ice)=convert_input_to_catesian_coordinates(x,y,ice_mask,h_ice,R_earth,dx,lat_ref)

        if set_all_thicknesses_to_one==True:
                print( 'All thicknesses being set equal to 1')
                if set_all_domain_to_ice is True:
                        ice_mask[:,:]=1.
                h_ice[np.where(ice_mask>0.5)]=Th_prescribed

        if Use_default_radius is True:
                dx=x[1]-x[0]  #This is used as the ocean grid size
                R_frac=0.45
                if element_type=='square':
                        Radius=R_frac*dx  #(less than half the grid size)
                if element_type=='hexagon':
                        Radius= (np.sqrt(3)/2.) *(R_frac*dx)   #(S is < half the grid size)
        print ('Radius set to R= ' , Radius)

        #Define the positions,thickness, mass,  of the icebergs

        (Number_of_bergs,lon,lat,iceberg_num,dx_berg, dy_berg,thickness,ang_vel_out,uvel_out,vvel_out,mass,width,x,y,Radius,static_berg,N_bergs_before_bd)= Create_icebergs(lon_init,lat_init,\
                                                                                                                                                                            Radius,R_earth, x, y,ice_mask,h_ice,A68,uvel_in,vvel_in,Convert_to_lat_lon,rho_ice,element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,\
                                                                                                                                                                            Interpolate_from_four_corners,Fill_in_the_boundaries, set_all_bergs_static_by_default,break_some_bonds,Remove_stationary_bergs,set_some_bergs_static_by_default,icethres)

        print ('Maximum thickness:', np.max(thickness), np.min(thickness))
        #Define the positions of the iceberg bonds
        if Create_icebergs_bonds==True:
                (Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon)=\
                                Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds,static_berg,Allow_bonds_for_static_iceberg\
                                ,N_bergs_before_bd,Allow_bonds_with_boundary_bergs)

        if Make_icebergs_non_static_later is True:
                Change_static_berg_after_bonding(Number_of_bergs,lat,lon, static_berg)


        temp_mass=[mass[i]/thickness[i] for i in range(len(thickness))]
        if only_choose_one_berg==True:
                (Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num)= Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num,static_berg)


        if regrid_icebergs_onto_grid==True:
                if scale_the_grid_to_lat_lon==True:
                        print ('Regridding should be run with scale_the_grid_to_lat_lon off')
                else:
                        if Switch_regridding_element_type==True:
                                if element_type=='square':
                                        element_type='hexagon'
                                else:
                                        element_type='square'
                        temp_thickness=[1. for i in thickness]
                        temp_mass=[mass[i]/(thickness[i]*rho_ice) for i in range(len(thickness))]
                        temp_rho_ice=1.
                        New_area=regrid_iceberg_thickness(lat,lon,Number_of_bergs,temp_thickness,temp_mass,h_ice,x,y,temp_rho_ice,element_type,static_berg,plot_outcome=False)
                        New_mass=regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type,static_berg,plot_outcome=False)
                if save_new_h_ice_file==True:
                        grid_area=(x[1]-x[0])*(x[1]-x[0])
                        h_ice_new=New_mass*0
                        M=h_ice.shape
                        for j in range(M[0]):
                                for i in range(M[1]):
                                        if New_area[j,i]>0:
                                                h_ice_new[j,i]=New_mass[j,i]/(rho_ice*New_area[j,i])

        if Switch_x_and_y_to_rotate_90 is True:
                [lat,lon, h_ice_new,New_area,New_mass]=Switch_x_and_y_directions(lat,lon, h_ice_new,New_area,New_mass)


        if save_restart_files==True:

                #Creating iceberg restart file
                Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,ang_vel_in,uvel_out,vvel_out,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg)

                #Creating bond restart file
                if Create_icebergs_bonds==True:
                        Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source)

                if save_new_h_ice_file==True:
                        create_clipped_icethickness_file(h_ice_new,New_area,New_mass,grid_area,gravity,Ice_geometry_source)

        # Plotting the positions and bonds of the newly formed formation
        if Run_plotting_subroutine==True:
                if Convert_to_lat_lon==True:
                        plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon,\
                                                   plot_circles,h_ice,ice_mask_grid,lon_grid,lat_grid,plot_ice_mask,\
                                                   plot_ice_thickness,thickness,plot_icebergs_positions,static_berg,\
                                                   h_ice_new,plot_h_ice_new)

                if (plot_bonds==True) and (Create_icebergs_bonds):
                        plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds)

        field=New_area/grid_area -1.
        M=field.shape
        if Run_plotting_subroutine:
                plt.show()
        print ('Script complete')



if __name__ == '__main__':
        optCmdLineArgs= parseCommandLine()
        main(optCmdLineArgs)
        #sys.exit(main())
