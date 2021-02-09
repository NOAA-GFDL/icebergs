#!/usr/bin/env python3

from netCDF4 import Dataset
import pdb
import netCDF4 as nc
import numpy as np
import matplotlib
import math
import os
matplotlib.use("tkagg")
from pylab import *
import matplotlib.pyplot as plt
import argparse

def parseCommandLine():
    parser = argparse.ArgumentParser(description=
    '''Plot energy from bond trajectories.''',
    epilog='Written by Alex Huth, 2020')
    parser.add_argument('-fname', type=str, default='bond_trajectories.nc',
                    help=''' provide filename to plot''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname

    with nc.Dataset(filename) as file:
        rot    = file.variables['rotation'][:]
        rrot    = file.variables['rel_rotation'][:]
        ns  = file.variables['n_strain'][:]
        nsr   = file.variables['n_strain_rate'][:]
        day   = file.variables['day'][:]

    print('rrotmin','rrotmax',np.min(nsr),np.max(nsr))

    ud = np.unique(day)
    dt = ud[1]-ud[0]

    print('dt',dt)

    fig = plt.figure(1)
    ax = plt.subplot(111)
    
    for j in range(len(ud)):

        #ENERGY STUFF
        rot1 = rot[day == ud[j]]
        rrot1 = rrot[day == ud[j]]
        ns1 = ns[day == ud[j]]
        nsr1 = nsr[day == ud[j]]
        ud1 = day[day == ud[j]]


        ax.plot(ud1, rot1,  color='k', label='rotation')
        ax.plot(ud1, rrot1, color='c', label='relative rot', linestyle='dashed')
        #ax.plot(ud1, ns1, color='b', label='normal strain', linestyle='dashed')
        ax.plot(ud1, nsr1, color='g', label='normal strain rate', linestyle='dashed')


    #ax.legend(loc="upper right")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)
    plt.xlabel('time (days)')
    plt.ylabel('strains')

    plt.show()

print('Script complete')

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
