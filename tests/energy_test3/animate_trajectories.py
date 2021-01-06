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
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname
    #filename = 'iceberg_trajectories.nc'
    field = 'lon'
    field2 = 'lat'

    with nc.Dataset(filename) as file:
        x = file.variables['lon'][:]/1.e3
        y = file.variables['lat'][:]/1.e3
        day = file.variables['day'][:]

    crop=False
    tc=1.25
    if crop:
        x=x[day<=tc]
        y=y[day<=tc]
        day=day[day<=tc]

    ud = np.unique(day)
    t = ud[0]

    # frame info
    num_frames = len(ud)
    movie_len = 10.0 #seconds
    frame_len = 1000.0*movie_len/num_frames

    xmin = 0
    xmax = 30 #45 #30
    ymin = 7.5 #0 #7.5
    ymax = 37.5 #45 #37.5

    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))

        scat.set_offsets(data)

        t = ud[i]
        time_text.set_text('time = %.1f days' % t )
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
    scat = ax1.scatter([],[],marker='o',facecolor='w',s=40,edgecolor='red')
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
