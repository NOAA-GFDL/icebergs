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
import xarray as xr
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import cartopy.crs as ccrs
import copy
from matplotlib.figure import figaspect
from matplotlib import ticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import matplotlib.ticker as mticker
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

from pyproj import Transformer
#from scipy import stats
#import pandas as pd
import matplotlib.ticker as plticker


from matplotlib.gridspec import GridSpec
import matplotlib.path as mpath
import matplotlib.dates as mdates
import cartopy.mpl.ticker as ctk

def parseCommandLine():
    t0=9+13./24
    parser = argparse.ArgumentParser(description=
    '''Generate animation of iceberg trajectories.''',
    epilog='Written by Alex Huth, 2020')
    parser.add_argument('-fname', type=str, default='iceberg_trajectories.nc',
                    help=''' provide filename to plot''')
    parser.add_argument('-s', type=int, default='400',
                        help='''plotted particle size''')
    parser.add_argument('-t0', type=float, default=t0,
                        help=''' initial time (days)''')
    parser.add_argument('-plotvar', type=bool, default=False,
                        help=''' plot a var''')
    parser.add_argument('-zoom', type=bool, default=False,
                        help=''' zoom to A68 grounding point''')
    optCmdLineArgs = parser.parse_args()
    return optCmdLineArgs



def main(args):

    print('reading file')

    filename=args.fname
    psize=args.s
    t0=args.t0
    pvar=args.plotvar
    zoom=args.zoom
    #if (zoom):
    #    psize=75
    #filename = 'iceberg_trajectories.nc'
    field = 'lon'
    field2 = 'lat'

    with nc.Dataset(filename) as file:
        x = file.variables['lon'][:]-360 #/1.e3
        y = file.variables['lat'][:] #/1.e3
        day = file.variables['day'][:]
        length = file.variables['length'][:]
        width = file.variables['width'][:]
        thick = file.variables['thickness'][:]
        od = file.variables['od'][:]
        nbonds = file.variables['n_bonds'][:]
        #uvel = file.variables['uvel_prev'][:]
        #vvel = file.variables['vvel_prev'][:]

    print('min day',day.min())
    day=day+t0

    thevar=nbonds
    # thres=17.5
    # x=x[day<thres]
    # y=y[day<thres]
    # length=length[day<thres]
    # width=width[day<thres]
    # thick=thick[day<thres]
    # od=od[day<thres]
    # day=day[day<thres]

    ud = np.unique(day)
    t = ud[0]

    #print('uinit',uvel[day==t])
    #print('vinit',vvel[day==t])

    # print('lon',x)
    # print('lat',y)
    radius = np.sqrt(length*width/4)#*(1./(2*np.sqrt(3)))

    #for determining if grounded:
    rho_bergs=850.
    rho_seawater=1025.
    h_to_ground=30.
    draught=(rho_bergs/rho_seawater)*thick
    groundfrac=1.-(od-draught)/h_to_ground
    groundfrac[groundfrac<0.]=0.
    groundfrac[groundfrac>1.]=1.

    groundfrac2=groundfrac*0
    groundfrac2[draught>od]=1.

    #groundfrac[thick<200.]=0.
    #groundfrac[thick>=200.]=1.


    # frame info
    num_frames = len(ud)
    movie_len = 5 #seconds
    #frame_len = 1000.0*movie_len/num_frames/10000000
    frame_len = movie_len/num_frames/1000

    if (zoom):
        xmin = -39 #38.5
        xmax = -36
        ymin = -57 #-56.5
        ymax = -54
    else:
        xmin = -41
        xmax = -36
        ymin = -58
        ymax = -54

    anim_running = True

    def animate(i):

        x1 = x[day == ud[i]]
        y1 = y[day == ud[i]]
        data = np.hstack((x1[:,np.newaxis],y1[:,np.newaxis]))

        scat.set_offsets(data)

        #t = ud[i]+14+13./24
        #t = ud[i]+16#9+13./24
        t=ud[i]
        time_text.set_text('date = Dec %.3f' % t )

        r1 = radius[day == ud[i]]
        scat.set_sizes(r1/psize)

        hstat = groundfrac[day == ud[i]]
        hstat2 = groundfrac2[day == ud[i]]
        cstring=[]
        cstring2=[]

        if pvar==True:
        #var1=od[day==ud[i]]
            var1=thevar[day==ud[i]]
            scat.set_color(cmap(norm(var1)))
        else:
            for j in range(len(hstat)):
                if (hstat2[j]==0):
                    if (hstat[j]==0):
                        cstring.append("none")
                        cstring2.append('r')
                    else:
                        cstring.append("b")
                        cstring2.append('b')
                else:
                    cstring.append("c")
                    cstring2.append('c')
            scat.set_color(cstring)
            scat.set_edgecolor(cstring2)

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


    minlon = xmin
    maxlon = xmax
    midlon = (minlon+maxlon)/2

    minlat = ymin
    maxlat = ymax
    midlat = (ymin+ymax)/2

    extent=[xmin, xmax, ymin, ymax]
    extent2=[xmin-1, xmax+1, ymin-1, ymax+1]
    extent3=[xmin-1, ymin-1, xmax+1, ymax+1]

    subplot_kws=dict(projection=ccrs.NearsidePerspective(central_longitude=midlon, central_latitude=midlat))

    bottom_x=np.vstack(np.linspace(minlon,maxlon,3)); bottom_y=bottom_x*0+minlat
    top_x=np.vstack(np.linspace(maxlon,minlon,3)); top_y=top_x*0+maxlat
    lons=np.vstack((bottom_x,top_x)); lats=np.vstack((bottom_y,top_y))
    boundary_path=np.hstack((lons,lats))
    boundary_path = mpath.Path(boundary_path)

    old=False

    if pvar==True:
        old=True

    # Now we can do the plotting!
    f = plt.figure(figsize=(8,8))
    f.tight_layout()

    if old:
        ax1 = plt.subplot(111,xlim=(xmin, xmax), ylim=(ymin, ymax))
        scat = ax1.scatter([],[],marker='o')#,edgecolor='red')
    else:
        ax1 = plt.axes(projection=ccrs.NearsidePerspective(central_longitude=midlon, central_latitude=midlat))
        scat = ax1.scatter([],[],transform=ccrs.PlateCarree(),marker='o')#,edgecolor='red')

    if pvar==True:
        div = make_axes_locatable(ax1)
        cax = div.append_axes('right', size='5%', pad='5%')
        cmap = plt.cm.get_cmap('RdYlBu')
        norm = plt.Normalize(thevar.min(), thevar.max())
        #norm = plt.Normalize(165.7, 165.8)
        cb1 = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap,
                                               norm=norm,
                                               orientation='vertical')


    time_text = ax1.text(0.02, 0.95, '', transform=ax1.transAxes)


    #plt.scatter(323.4-360,-56.6,marker='o')
    #plt.scatter(322.52-360,-55.39,marker='o')

    xt=np.array([-37.32708,-37.05308,-37.13542,-37.32708])
    yt=np.array([-55.01042,-55.01042,-54.93958,-55.01042])
    plt.plot(xt,yt)


    if (pvar==False):
        dir='./data'
        oscar=xr.load_dataset(f'{dir}/a68_experiment_ocean_surf_vel_oscar_dec2020_HOURLY_ll_p125.nc')
        #oscar.isel(time=400).plot.quiver(x="lon",y="lat",u="uo",v="vo",hue='vel',vmax=0.2, scale=8)
        vv=oscar.isel(time=400)
        norm=plt.cm.colors.Normalize(vmin=0,vmax=0.35)
        if (old==False):
            ax1.quiver(vv.lon.values,vv.lat.values,vv.uo.values,vv.vo.values,vv.vel.values,scale=8,cmap='viridis',norm=norm,transform=ccrs.PlateCarree())
        else:
            ax1.quiver(vv.lon,vv.lat,vv.uo,vv.vo,vv.vel,scale=8,cmap='viridis',norm=norm)
        f.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'),ax=ax1, shrink=0.75)
        #330 is dec 14th at 5PM. 408 is dec 16th at 5 PM


    if old==True:
        ax1.set_aspect('equal', adjustable='box')

        # Change major ticks to show every 20.
        ax1.xaxis.set_major_locator(MultipleLocator(5))
        ax1.yaxis.set_major_locator(MultipleLocator(5))

        ax1.xaxis.set_minor_locator(MultipleLocator(1))
        ax1.yaxis.set_minor_locator(MultipleLocator(1))

        # Turn grid on for both major and minor ticks and style minor slightly
        # differently.
        ax1.grid(which='major', color='#CCCCCC', linestyle=':')
        ax1.grid(which='minor', color='#CCCCCC', linestyle=':')

        ax1.set_xlabel('lon')
        ax1.set_ylabel('lat')
        ax1.set_title('Iceberg trajectory')
    else:
        xticks=list(range(minlon,maxlon+1)); yticks=list(range(minlat,maxlat+1))
        dir_labels=True
        proj_to_data = ccrs.PlateCarree()._as_mpl_transform(ax1) - ax1.transData
        bound_in_target = proj_to_data.transform_path(boundary_path)
        ax1.set_extent(extent2, crs=ccrs.PlateCarree())
        ax1.set_boundary(bound_in_target)
        g=ax1.gridlines(color='black',alpha=0.5,linestyle='--',draw_labels=True,x_inline=False, y_inline=False,xlocs=xticks,ylocs=yticks)
        g.top_labels=False; g.right_labels=False; g.rotate_labels=False
        #g.xformatter=LongitudeFormatter(direction_label=dir_labels); g.yformatter=LatitudeFormatter(direction_label=dir_labels)
        g.xlabel_style = {'size': 12, 'color': 'gray'}; g.ylabel_style = {'size': 12, 'color': 'gray'}


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
