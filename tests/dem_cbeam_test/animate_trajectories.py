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
    parser.add_argument('-s', type=int, default='100000',
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
        hs = file.variables['thickness'][:]


    radius = length*width*(1./(2*np.sqrt(3)))

    #print('day',day)
    crop=False
    tc=1.25
    if crop:
        x=x[day<=tc]
        y=y[day<=tc]
        day=day[day<=tc]

    ud = np.unique(day)
    #print('ud',ud)
    t = ud[0]

    # frame info
    num_frames = len(ud)
    movie_len = 4.0 #seconds
    frame_len = 1000.0*movie_len/num_frames

    xshift=101
    yshift=(151.e3+5000.0)/1.e3

    xmin = 90-xshift
    xmax = 260-xshift
    ymin = 0-yshift
    ymax = 170-yshift

    x=x-xshift
    y=y-yshift

    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        # hstat = hs[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))
        hstat= day[day == ud[i]]

        ud_max=np.max(ud)

        cstring=[]

        for j in range(len(hstat)):
            if (hstat[j]<ud_max):
                cstring.append("r")
            else:
                cstring.append("b")

        scat.set_offsets(data)
        scat.set_color(cstring)
        scat.set_edgecolor('k')
        #scat2.set_offsets(data)

        t = ud[i]
        time_text.set_text('time = %.1f days' % t )

        r1 = radius[day == ud[i]]

        scat.set_sizes(r1/psize)
        return scat,time_text #,scat2

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
    #scat = ax1.scatter([],[],marker='o',facecolor='w',s=psize,edgecolor='red')
    scat = ax1.scatter([],[],marker='o')#,facecolor='w',edgecolor='red')
    #scat2 = ax1.scatter([],[],marker='.',color='k',s=5)
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

    # Change major ticks to show every 20.
    ax1.xaxis.set_major_locator(MultipleLocator(50))
    ax1.yaxis.set_major_locator(MultipleLocator(50))

    ax1.xaxis.set_minor_locator(MultipleLocator(10))
    ax1.yaxis.set_minor_locator(MultipleLocator(10))

    # Turn grid on for both major and minor ticks and style minor slightly
    # differently.
    ax1.grid(which='major', color='#CCCCCC', linestyle=':')
    ax1.grid(which='minor', color='#CCCCCC', linestyle=':')

    ax1.set_xlabel('x (km)')
    ax1.set_ylabel('displacement (km)')
    ax1.set_title('Iceberg trajectory')


    #analytical solution
    thick=1
    xa=np.linspace(0,150000,100)
    P=-(1.5e10)#/3. #Newtons
    l=29*5000 #total beam length (m) between center of start and end particles (x-dir)
    h=3.*5000 #total beam height (m) between center of start and end particles (y-dir)
    AI= thick*(h**3)/12. #second  moment of area
    YM=1.e9 #young's modulus (Pa)
    w= P*(xa**2.)*(3.*l-xa)/(6.*YM*AI) #deflection

    xa=(xa)/1.e3
    w=w/1.e3
    ax1.plot(xa,w)

    f.canvas.mpl_connect('button_press_event', onClick)
    ani = FuncAnimation(
        f,animate,init_func=init,frames=num_frames,
        interval=frame_len,blit=True,repeat=False)

    print('time',t)

    plt.show()
    #animation.save("iceberg_traj_animation.mp4")


print('Script complete')

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
