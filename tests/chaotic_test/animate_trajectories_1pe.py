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
    parser.add_argument('-fname', type=str, default='iceberg_trajectories.nc.0000',
                    help=''' provide filename to plot''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname
    field = 'lon'
    field2 = 'lat'

    with nc.Dataset(filename) as file:
        x = file.variables['lon'][:]/1.e4
        y = file.variables['lat'][:]/1.e4
        day = file.variables['day'][:]
        hs = file.variables['halo_berg'][:]
        idij = file.variables['id_ij'][:]
        idcnt = file.variables['id_cnt'][:]

    ud = np.unique(day)
    t = ud[0]

    print('unique hs',np.unique(hs))

    # frame info
    num_frames = len(ud)
    movie_len = 5.0 #seconds
    frame_len = 1000.0*movie_len/num_frames

    xmin = 25
    xmax = 65
    ymin = 40
    ymax = 70

    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        hstat = hs[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))
        idd = idij[day == ud[i]]
        idc = idcnt[day == ud[i]]

        cstring=[]
        for j in range(len(hstat)):

            if (hstat[j]<0):
                cstring.append("y")
            elif (hstat[j]==0):
                cstring.append("b")
            elif (hstat[j]==1):
                cstring.append("r")
            elif (hstat[j]==2):
                cstring.append("g")
            elif (hstat[j]==10):
                cstring.append("c")
            else:
                cstring.append("k")

        scat.set_offsets(data)
        scat.set_color(cstring)
        scat.set_edgecolor('none')

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
    scat = ax1.scatter([],[],marker='o',s=6)
    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

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
        interval=frame_len,blit=True,repeat=False)

    print('time',t)

    jsd=49
    jsc=50
    jec=60
    jed=61.5
    isd=39
    isc=40
    iec=50
    ied=51.5

    plt.plot([isc,isc],[jsc,jec],'k-',lw=1)
    plt.plot([iec,iec],[jsc,jec],'k-',lw=1)
    plt.plot([isc,iec],[jsc,jsc],'k-',lw=1)
    plt.plot([isc,iec],[jec,jec],'k-',lw=1)

    plt.plot([isd,isd],[jsd,jed],'b:',lw=1)
    plt.plot([ied,ied],[jsd,jed],'b:',lw=1)
    plt.plot([isd,ied],[jsd,jsd],'b:',lw=1)
    plt.plot([isd,ied],[jed,jed],'b:',lw=1)

    plt.show()
    #animation.save("iceberg_traj_animation.mp4")


print('Script complete')

if __name__ == '__main__':
    optCmdLineArgs=	parseCommandLine()
    anim_running = True
    main(optCmdLineArgs)
