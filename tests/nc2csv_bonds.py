#!/usr/bin/env python

from netCDF4 import Dataset
import pdb
import netCDF4 as nc
import numpy as np
import csv
import math
import os
import argparse

#CONVERT BOND TRAJECTORIES FROM NETCDF TO CSV FILE FOR VISUALIZATION IN PARAVIEW
#returns a series of csv files (one for each time step)
#To view csv in paraview:
#1. load all csv files in the series
#2. Filters --> Table To Points
#3. In the TableToPoints filter, assign x for X Column, and y for Y Column. Check 2D Points.
#4. Use calculator to define a vector field for visualizing the bond:
#   bond_vect = n1*iHat+n2*iHat
#5. Make a line glyph w/ orientation given by bond_vect components and magnitude by length


def parseCommandLine():
    parser = argparse.ArgumentParser(description=
                                     '''Convert iceberg bond trajectories netcdf file to csv file.''',
                                     epilog='Written by Alex Huth, 2020')
    parser.register('type','bool',str2bool) # add type keyword to registries
    parser.add_argument('-fin', type=str, default='bond_trajectories.nc',
                        help=''' provide netcdf filename''')
    parser.add_argument('-fout', type=str, default='./csvresults/bond_trajectories',
                        help=''' provide filename prefix for csv''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs

def main(args ):

    print('reading file')

    fin=args.fin
    fout=args.fout

    with nc.Dataset(fin) as file:
        lon = file.variables['lon'][:]
        lat = file.variables['lat'][:]
        year = file.variables['year'][:]
        day = file.variables['day'][:]
        length = file.variables['length'][:]
        n1 = file.variables['n1'][:]
        n2 = file.variables['n2'][:]
        #rot = file.variables['rotation'][:]
        #rrot = file.variables['rel_rotation'][:]
        #ns = file.variables['n_strain'][:]
        #nsr = file.variables['n_strain_rate'][:]
        ns = file.variables['nstress'][:]
        ss = file.variables['sstress'][:]
        id_cnt1 = file.variables['id_cnt1'][:]
        id_ij1 = file.variables['id_ij1'][:]
        id_cnt2 = file.variables['id_cnt2'][:]
        id_ij2 = file.variables['id_ij2'][:]
        #fields = ['x','y','time','len','n1','n2','rot','rrot','ns','nsr','berg_id1','berg_id2']
        fields = ['x','y','time','len','n1','n2','ns','ss','berg_id1','berg_id2']
    idmid=int('2',8)**32
    btime=year+day/365.15
    bidin1=id_cnt1*idmid+id_ij1
    bidin2=id_cnt2*idmid+id_ij2

    # uid=np.unique(bidin1)
    # uid2=np.unique(bidin2)
    # print('uid',uid)
    # print('uid2',uid2)

    #there are originally 2 copies of each bond (one copy for each particle in the bond)
    #Remove the duplicates:
    biddiff=bidin2-bidin1

    lon = lon[biddiff>0]
    lat = lat[biddiff>0]
    btime = btime[biddiff>0]
    length = length[biddiff>0]
    n1 = n1[biddiff>0]
    n2 = n2[biddiff>0]
    # rot=rot[biddiff>0]
    # rrot=rrot[biddiff>0]
    # ns=ns[biddiff>0]
    # nsr=nsr[biddiff>0]
    ns=ns[biddiff>0]
    ss=ss[biddiff>0]
    bidin1 = bidin1[biddiff>0]
    bidin2 = bidin2[biddiff>0]

    utime=np.unique(btime)

    for i in range(len(utime)):
        x = lon[btime == utime[i]]
        y = lat[btime == utime[i]]
        time = btime[btime == utime[i]]
        blen = length[btime == utime[i]]
        bn1 = n1[btime == utime[i]]
        bn2 = n2[btime == utime[i]]
        # brot = rot[btime == utime[i]]
        # brrot = rrot[btime == utime[i]]
        # bns = ns[btime == utime[i]]
        # bnsr = nsr[btime == utime[i]]
        bns = ns[btime == utime[i]]
        bss = ss[btime == utime[i]]
        bid1 = bidin1[btime == utime[i]]
        bid2 = bidin2[btime == utime[i]]

        filename = fout+str(i)+".csv"
        with open(filename, 'w') as file:
            writer = csv.writer(file)
            writer.writerow(fields)
            for j in range(len(x)):
               # writer.writerow([x[j],y[j],time[j],blen[j],bn1[j],bn2[j],brot[j],brrot[j],bns[j],bnsr[j],bid1[j],bid2[j]])
                writer.writerow([x[j],y[j],time[j],blen[j],bn1[j],bn2[j],bns[j],bss[j],bid1[j],bid2[j]])
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
