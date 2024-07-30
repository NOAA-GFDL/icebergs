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
    parser.add_argument('-s', type=float, default='0.0025',
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
        x = file.variables['lon'][:] #/1.e3
        y = file.variables['lat'][:] #/1.e3
        day = file.variables['day'][:]*24*60*60
        length = file.variables['length'][:]
        width = file.variables['width'][:]
        hs = file.variables['thickness'][:]


    radius = 0.5*sqrt(length*width) #*(1./(2*np.sqrt(3)))

    #print('day',day)
    crop=False
    tc=2e-5#1.25
    if crop:
        x=x[day<=tc]
        y=y[day<=tc]
        radius=radius[day<=tc]
        day=day[day<=tc]

    ud = np.unique(day)
    print('ud',ud)
    t = ud[0]

    # frame info
    num_frames = len(ud)
    movie_len = 4.0 #seconds
    frame_len = 1000.0*movie_len/num_frames

    d=0.5 #5000.0
    r=d/2.
    xshift=101.e3
    yshift=(151.e3+2*r) #/1.e3
    x=x-xshift
    y=y-yshift

    xmin=np.min(x)-1
    xmax=np.max(x)+1
    ymin=np.min(y)-1
    ymax=np.max(y)+1


    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))
        hstat= day[day == ud[i]]

        ud_max=np.max(ud)

        cstring=[]

        for j in range(len(hstat)):
            cstring.append("none")

        scat.set_offsets(data)
        scat.set_color(cstring)
        scat.set_edgecolor('k')

        t = ud[i]
        time_text.set_text('time = %.1f sec' % t )

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
    scat = ax1.scatter([],[],marker='o')
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

    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('displacement (m)')
    ax1.set_title('Iceberg trajectory')

    #analytical solution
    xa=np.arange(0,29)
    xa=xa*0.5
    print('xa',xa)
    l=np.max(xa)
    print('l',l)
    P=-(1.5e5) #Newtons
    YM=1.e9 #young's
    thick=1
    h=0.5
    AI= thick*(h**3)/12. #second  moment of area
    w1=-P*xa*(4.*xa*xa-3.*l*l)/(48.*YM*AI)
    w2=P*(xa-l)*(l*l-8.*l*xa+4*xa*xa)/(48.*YM*AI)
    w=w1
    w[xa>0.5*l]=w2[xa>0.5*l] #deflection
    print('w',w)

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
