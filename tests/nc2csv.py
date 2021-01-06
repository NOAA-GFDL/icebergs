#!/usr/bin/env python

from netCDF4 import Dataset
import pdb
import netCDF4 as nc
import numpy as np
import csv
import math
import os
import argparse

#CONVERT BERG TRAJECTORIES FROM NETCDF TO CSV FILE FOR VISUALIZATION IN PARAVIEW
#returns a series of csv files (one for each time step)
#To view csv in paraview:
#1. load all csv files in the series
#2. Filters --> Table To Points
#3. In the TableToPoints filter, assign x for X Column, and y for Y Column. Check 2D Points.


def parseCommandLine():
    parser = argparse.ArgumentParser(description=
                                     '''Convert iceberg trajectories netcdf file to csv file.''',
                                     epilog='Written by Alex Huth, 2020')
    parser.register('type','bool',str2bool) # add type keyword to registries
    parser.add_argument('-fin', type=str, default='iceberg_trajectories.nc',
                        help=''' provide netcdf filename''')
    parser.add_argument('-fout', type=str, default='./csvresults/iceberg_trajectories',
                        help=''' provide filename prefix for csv''')
    parser.add_argument('-short_traj', type='bool', default=False,
		        help=''' indicate whether save_short_traj was used ''')
    parser.add_argument('-mts', type='bool', default=True,
       	                help=''' indicate whether mts was used ''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args ):

    print('reading file')

    fin=args.fin
    fout=args.fout
    short_traj=args.short_traj
    mts=args.mts
    if short_traj==True:
        mts=False

    with nc.Dataset(fin) as file:
        lon = file.variables['lon'][:]
        lat = file.variables['lat'][:]
        year = file.variables['year'][:]
        day = file.variables['day'][:]
        id_cnt = file.variables['id_cnt'][:]
        id_ij = file.variables['id_ij'][:]
        if short_traj==True:
            fields = ['x','y','time','id']
        else:
            uvel = file.variables['uvel'][:]
            vvel = file.variables['vvel'][:]
            # uo = file.variables['uo'][:]
            # vo = file.variables['vo'][:]
            # ui = file.variables['ui'][:]
            # vi = file.variables['vi'][:]
            # ua = file.variables['ua'][:]
            # va = file.variables['va'][:]
            pmass = file.variables['mass'][:]
            # mass_of_bits= file.variables['mass_of_bits'][:]
            # heat_density= file.variables['heat_density'][:]
            thickness = file.variables['thickness'][:]
            width = file.variables['width'][:]
            length = file.variables['length'][:]
            # ssh_x = file.variables['ssh_x'][:]
            # ssh_y = file.variables['ssh_y'][:]
            # sst = file.variables['sst'][:]
            # sss = file.variables['sss'][:]
            # cn = file.variables['cn'][:]
            # hi = file.variables['hi'][:]
            axn = file.variables['axn'][:]
            ayn = file.variables['ayn'][:]
            bxn = file.variables['bxn'][:]
            byn = file.variables['byn'][:]
            #halo_berg = file.variables['halo_berg'][:]
            ocean_depth  = file.variables['od'][:]
            if mts==True:
                nbonds_id = file.variables['n_bonds'][:]
                #uvprev = file.variables['uvel_prev'][:]
                #vvprev = file.variables['vvel_prev'][:]
                axnf = file.variables['axn_fast'][:]
                aynf = file.variables['ayn_fast'][:]
                bxnf = file.variables['bxn_fast'][:]
                bynf = file.variables['byn_fast'][:]
                tot_axn = file.variables['tot_accel_x'][:]
                tot_ayn = file.variables['tot_accel_y'][:]
                fields = ['x','y','time','id','uv','vv','mass','H','W','L','axn','ayn','bxn','byn','ocean_depth','axnf','aynf','bxnf','bynf','nbonds_id','accelx','accely']
            else:
                fields = ['x','y','time','id','uv','vv','mass','H','W','L','axn','ayn','bxn','byn','ocean_depth']

    idmid=int('2',8)**32
    ptime=year+day/365.15
    pidin=id_cnt*idmid+id_ij
    utime=np.unique(ptime)

    for i in range(len(utime)):
        x = lon[ptime == utime[i]]
        y = lat[ptime == utime[i]]
        time = ptime[ptime == utime[i]]
        pid = pidin[ptime == utime[i]]
        if short_traj==False:
            uv = uvel[ptime == utime[i]]
            vv = vvel[ptime == utime[i]]
            mass = pmass[ptime == utime[i]]
            H = thickness[ptime == utime[i]]
            W = width[ptime == utime[i]]
            L = length[ptime == utime[i]]
            ax = axn[ptime == utime[i]]
            ay = ayn[ptime == utime[i]]
            bx = bxn[ptime == utime[i]]
            by = byn[ptime == utime[i]]
            od = ocean_depth[ptime == utime[i]]
            if mts==True:
                nb = nbonds_id[ptime == utime[i]]
                axf = axnf[ptime == utime[i]]
                ayf = aynf[ptime == utime[i]]
                bxf = bxnf[ptime == utime[i]]
                byf = bynf[ptime == utime[i]]
                tot_ax = tot_axn[ptime == utime[i]]
                tot_ay = tot_ayn[ptime == utime[i]]
        filename = fout+str(i)+".csv"
        with open(filename, 'w') as file:
            writer = csv.writer(file)
            writer.writerow(fields)
            if short_traj==True:
                for j in range(len(x)):
                    writer.writerow([x[j],y[j],time[j],pid[j]])
            else:
                if mts==True:
                    for j in range(len(x)):
                        writer.writerow([x[j],y[j],time[j],pid[j],uv[j],vv[j],mass[j],H[j],W[j],L[j],ax[j],ay[j],bx[j],by[j],od[j],axf[j],ayf[j],bxf[j],byf[j],nb[j],tot_ax[j],tot_ay[j]])
                else:
                    for j in range(len(x)):
                        writer.writerow([x[j],y[j],time[j],pid[j],uv[j],vv[j],mass[j],H[j],W[j],L[j],ax[j],ay[j],bx[j],by[j],od[j]])

print('Script complete')

def str2bool(string):
	if string.lower() in  ("yes", "true", "t", "1"):
		Value=True
	elif string.lower() in ("no", "false", "f", "0"):
		Value=False
	else:
		print '**********************************************************************'
		print 'The input variable ' ,str(string) ,  ' is not suitable for boolean conversion, using default'
		print '**********************************************************************'

		Value=None
		return

	return Value

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
