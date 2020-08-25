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


print('reading file')

filename = 'iceberg_trajectories.nc'
field = 'lon'
field2 = 'lat'

with nc.Dataset(filename) as file:
    x = file.variables['lon'][:]/1.e3
    y = file.variables['lat'][:]/1.e3
    day = file.variables['day'][:]

ud = np.unique(day)
t = ud[0]

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

# frame info
num_frames = len(ud)
movie_len = 5.0 #seconds
frame_len = 1000.0*movie_len/num_frames

smax = max(max(x),max(y),20)
xmin = 0.0
xmax = smax
ymin = 0.0
ymax = smax


# Now we can do the plotting!
f = plt.figure(figsize=(5,5))
f.tight_layout()
ax1 = plt.subplot(111,xlim=(xmin, xmax), ylim=(ymin, ymax))
scat = ax1.scatter([],[],marker='o',facecolor='w',s=110,edgecolor='red')
time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)

ax1.set_xticks(np.arange(0, xmax+1, 1)) 
ax1.set_yticks(np.arange(0, ymax+1, 1)) 
ax1.set_xlabel('x (km)')
ax1.set_ylabel('y (km)')
ax1.set_title('Iceberg trajectory')

anim_running = True
f.canvas.mpl_connect('button_press_event', onClick)
 
ani = FuncAnimation(
    f,animate,init_func=init,frames=num_frames,
    interval=frame_len,blit=True)

print('time',t)

plt.grid()
plt.show()

#animation.save("iceberg_traj_animation.mp4")



#old:
#mid = int(np.ceil(len(ud)/2))
#t1 = ud[0]
#t2 = ud[mid]
#t3 = ud[len(ud)-1]
#print('plotting')
#plt.figure()
#figsize=(15,10))
#plt.scatter(x[day == t1],y[day == t1],color = 'blue', marker='*')
#plt.scatter(x[day == t2],y[day == t2],color = 'green',marker='o')
#plt.scatter(x[day == t3],y[day == t3],color = 'red',marker='*')
#plt.xlabel('x-extent (km)')
#plt.ylabel('y-extent (km)')
#plt.show()


print('Script complete')







