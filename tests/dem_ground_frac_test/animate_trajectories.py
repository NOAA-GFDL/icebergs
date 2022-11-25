#!/usr/bin/env python

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
from matplotlib.animation import FuncAnimation
import argparse

def parseCommandLine():
    parser = argparse.ArgumentParser(description=
    '''Generate animation of iceberg trajectories.''',
    epilog='Written by Alex Huth, 2020')
    parser.add_argument('-fname', type=str, default='iceberg_trajectories.nc',
                    help=''' provide filename to plot''')
    parser.add_argument('-s', type=int, default='14000',
                        help='''plotted particle size''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname
    psize=args.s
    #filename = 'iceberg_trajectories.nc'
    field = 'lon'
    field2 = 'lat'

    with nc.Dataset(filename) as file:
        x = file.variables['lon'][:]/1.e3
        y = file.variables['lat'][:]/1.e3
        day = file.variables['day'][:]
        length = file.variables['length'][:]
        width = file.variables['width'][:]
        thick = file.variables['thickness'][:]
        od = file.variables['od'][:]

    ud = np.unique(day)
    t = ud[0]

    radius = length*width*(1./(2*np.sqrt(3)))

    #for determining if grounded:
    rho_bergs=850.
    rho_seawater=1025.
    h_to_ground=20 #200.
    draught=(rho_bergs/rho_seawater)*thick
    if (h_to_ground==0):
        groundfrac=draught-od
    else:
        groundfrac=1.-(od-draught)/h_to_ground
    groundfrac[groundfrac<0.]=0.
    groundfrac[groundfrac>1.]=1.

    #groundfrac[thick<200.]=0.
    #groundfrac[thick>=200.]=1.


    # frame info
    num_frames = len(ud)
    movie_len = 5 #seconds
    #frame_len = 1000.0*movie_len/num_frames/10000000
    frame_len = movie_len/num_frames/1000

    xmin = 40
    xmax = 100
    ymin = 20 #xmin
    ymax = 80 #xmax

    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))

        scat.set_offsets(data)

        t = ud[i]
        time_text.set_text('time = %.1f days' % t )

        r1 = radius[day == ud[i]]
        scat.set_sizes(r1/psize)

        hstat = groundfrac[day == ud[i]]
        cstring=[]
        for j in range(len(hstat)):
            if (hstat[j]==0):
                cstring.append("none")
            else:
                cstring.append("b")
        scat.set_color(cstring)
        scat.set_edgecolor('r')

        return scat,time_text

    def init():
        scat.set_offsets([])
        return scat,

    def onClick(event):
        global anim_running
        if anim_running:
            ani.event_source.stop()
            anim_running = False
        else:
            ani.event_source.start()
            anim_running = True

    # Now we can do the plotting!
    f = plt.figure(figsize=(5,5))
    f.tight_layout()
    ax1 = plt.subplot(111,xlim=(xmin, xmax), ylim=(ymin, ymax))
    scat = ax1.scatter([],[],marker='o')#,edgecolor='red')
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

    # Change major ticks to show every 20.
    ax1.xaxis.set_major_locator(MultipleLocator(5))
    ax1.yaxis.set_major_locator(MultipleLocator(5))

    ax1.xaxis.set_minor_locator(MultipleLocator(1))
    ax1.yaxis.set_minor_locator(MultipleLocator(1))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax1.grid(which='major', color='#CCCCCC', linestyle=':')
    ax1.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('y (km)')
    ax1.set_title('Iceberg trajectory')

    f.canvas.mpl_connect('button_press_event', onClick)

    ani = FuncAnimation(
        f,animate,init_func=init,frames=num_frames,
        interval=frame_len,blit=True,repeat=True)

    print('time',t)

    plt.show()
    #animation.save("iceberg_traj_animation.mp4")


print('Script complete')

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
