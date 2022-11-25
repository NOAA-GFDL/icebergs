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
    '''Plot energy from iceberg trajectories.''',
    epilog='Written by Alex Huth, 2020')
    parser.add_argument('-fname', type=str, default='iceberg_trajectories.nc',
                    help=''' provide filename to plot''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname

    with nc.Dataset(filename) as file:
        Ee    = file.variables['Ee'][:]
        Ed    = file.variables['Ed'][:]
        Eext  = file.variables['Eext'][:]
        Eec   = file.variables['Ee_contact'][:]
        Edc   = file.variables['Ed_contact'][:]
        Efrac = file.variables['Efrac'][:]
        day   = file.variables['day'][:]
        vx    = file.variables['uvel_prev'][:]
        vy    = file.variables['vvel_prev'][:]
        mass  = file.variables['mass'][:]
        x    = file.variables['lon'][:]
        y    = file.variables['lat'][:]
        vx2    = file.variables['uvel'][:]
        vy2    = file.variables['vvel'][:]

    crop=False
    tc=1.0
    if crop:
        Ee = Ee[day<=tc]
        Ed = Ed[day<=tc]
        Eext=Eext[day<=tc]
        Eec=Eec[day<=tc]
        Edc=Edc[day<=tc]
        Efrac=Efrac[day<=tc]
        vx=vx[day<=tc]
        vy=vy[day<=tc]
        mass=mass[day<=tc]
        x=x[day<=tc]
        y=y[day<=tc]
        vx2=vx2[day<=tc]
        vy2=vy2[day<=tc]       
        day=day[day<=tc]



    # vx=np.copy(vx2); vy=np.copy(vy2)

    ud = np.unique(day)
    step0 = 0
    dt = ud[1]-ud[0]

    print('dt',dt)

    tot_Ee = np.zeros(len(ud)-step0)
    tot_Ed = np.zeros(len(ud)-step0)
    tot_Eext = np.zeros(len(ud)-step0)
    tot_Eec = np.zeros(len(ud)-step0)
    tot_Edc = np.zeros(len(ud)-step0)
    tot_Efrac = np.zeros(len(ud)-step0)
    tot_Ek = np.zeros(len(ud)-step0)
    tot_E  = np.zeros(len(ud)-step0)
    tot_momentum_x = np.zeros(len(ud)-step0)
    tot_momentum_y = np.zeros(len(ud)-step0)
    tot_momentum = np.zeros(len(ud)-step0)
    tot_momentum_x2 = np.zeros(len(ud)-step0)
    tot_momentum_y2 = np.zeros(len(ud)-step0)
    tot_momentum2 = np.zeros(len(ud)-step0)

    # tot_vel_diff = np.zeros(len(ud)-step0)


    #regular version
    mom_x0 = mass[day == ud[step0]]*vx[day == ud[step0]]
    mom_y0 = mass[day == ud[step0]]*vy[day == ud[step0]]
    mom_x0_tot = np.sum(mom_x0); mom_y0_tot = np.sum(mom_y0)

    #'uncorrected' version
    mom_x02 = mom_x0; mom_y02=mom_y0
    mom_x0_tot2=mom_x0_tot; mom_y0_tot2=mom_y0_tot

    #if needed, 'correct' the regular version by changing the reference frame
    #so that mom_x0_tot>0 and mom_y0_tot>0
    dvx=0.0; dvy=0.0
    if (mom_x0_tot==0):
        print('changing x reference frame!')
        dvx=0.01#e-8
        mom_x0 = mass[day == ud[step0]]*(vx[day == ud[step0]]+dvx)
        mom_x0_tot = np.sum(mom_x0)
    if (mom_y0_tot==0):
        print('changing y reference frame!')
        dvy=0.01#e-8
        mom_y0 = mass[day == ud[step0]]*(vy[day == ud[step0]]+dvy)
        mom_y0_tot = np.sum(mom_y0)


    for i in range(len(ud)-step0):

        j=i+step0

        #ENERGY STUFF
        Ee1 = Ee[day == ud[j]]
        Ed1 = Ed[day == ud[j]]
        Eec1 = Eec[day == ud[j]]
        Edc1 = Edc[day == ud[j]]
        Efrac1 = Efrac[day == ud[j]]
        Eext1 = Eext[day == ud[j]]
        Ek1 = 0.5*mass[day == ud[j]]*(vx[day == ud[j]]**2+vy[day == ud[j]]**2)

        tot_Ee[i] = np.sum(Ee1)
        tot_Ed[i] = np.sum(Ed1)
        tot_Eec[i] = np.sum(Eec1)
        tot_Edc[i] = np.sum(Edc1)

        tot_Efrac[i] = np.sum(Efrac1)
        tot_Eext[i] = np.sum(Eext1)
        tot_Ek[i] = np.sum(Ek1)

        tot_E[i] = tot_Ee[i]+tot_Ed[i]+tot_Eec[i]+tot_Edc[i]+tot_Efrac[i]+tot_Eext[i]+tot_Ek[i]

        #MOMENTUM STUFF

        #--reg version, corrected--
        mom_x = mass[day == ud[j]] * (vx[day == ud[j]]+dvx)
        mom_y = mass[day == ud[j]] * (vy[day == ud[j]]+dvy)
        tot_momentum_x[i] = (np.sum(mom_x)-mom_x0_tot)/mom_x0_tot
        tot_momentum_y[i] = (np.sum(mom_y)-mom_y0_tot)/mom_y0_tot

        #--reg version, uncorrected--
        mom_x2 = mass[day == ud[j]] * (vx[day == ud[j]])
        mom_y2 = mass[day == ud[j]] * (vy[day == ud[j]])
        tot_momentum_x2[i] = (np.sum(mom_x2)-mom_x0_tot2)/mom_x0_tot2
        tot_momentum_y2[i] = (np.sum(mom_y2)-mom_y0_tot2)/mom_y0_tot2

        # ua =  vx[day == ud[j]]; va =  vy[day == ud[j]]
        # ub = vx2[day == ud[j]]; vb = vy2[day == ud[j]]
        # tot_vel_diff[i] = np.sum(sqrt((ua-ub)**2 + (va-vb)**2))


    ud=ud[0:len(ud)-step0]

    fig = plt.figure(1)
    ax = plt.subplot(111)
    ax.plot(ud, tot_E,  color='k', label='total')
    ax.plot(ud, tot_Eext, color='c', label='external', linestyle='dashed')
    ax.plot(ud, tot_Ee, color='b', label='elastic (bonded)', linestyle='dashed')
    ax.plot(ud, tot_Ed, color='g', label='damping (bonded)', linestyle='dashed')
    ax.plot(ud, tot_Eec, color='b', label='elastic (collision)', linestyle='dotted')
    ax.plot(ud, tot_Edc, color='g', label='damping (collision)', linestyle='dotted')
    ax.plot(ud, tot_Ek, color='r', label='kinetic', linestyle='dotted')
    ax.plot(ud, tot_Efrac, color='m', label='dissipated fracture', linestyle='dotted')

    print('min/max tot_E',np.min(tot_E),np.max(tot_E))

    #ax.legend(loc="upper right")
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)
    plt.xlabel('time (days)')
    plt.ylabel('energy (J)')
    #axes=plt.gca()
    #axes.set_ylim([-1.e-14,1.e-14])

    fig2 = plt.figure(2)
    ax = plt.subplot(111)
    ax.plot(ud[step0:], tot_momentum_x[step0:],  color='k', label='x momentum')
    ax.plot(ud[step0:], tot_momentum_x2[step0:],  color='y', label='x momentum 2', linestyle='dashed')
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)


    fig3 = plt.figure(3)
    ax = plt.subplot(111)
    ax.plot(ud[step0:], tot_momentum_y[step0:],  color='k', label='y momentum')
    if (abs(mom_y0_tot2)>0):
        ax.plot(ud[step0:], tot_momentum_y2[step0:],  color='y', label='y momentum 2', linestyle='dashed')
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)


    # fig4 = plt.figure(4)
    # ax = plt.subplot(111)
    # ax.plot(ud[step0:], tot_vel_diff[step0:],  color='k', label='vel diff')
    # chartBox = ax.get_position()
    # ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
    # ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)


    plt.show()

print('Script complete')

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
