#!/usr/bin/env python

###############################################################################################
from netCDF4 import Dataset
import sys
import argparse
import numpy as np
import os
import numpy
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
from pylab import *
#import pdb
#import argparse
#This combination allows you to access the varibale imputted from the command line.
#import sys


def load_timeseries_data_from_file(filename,field):
	print 'Loading data from file: ', filename	
	with nc.Dataset(filename) as file:
			data = file.variables[field][:]
			Time = file.variables['time'][:]
			x = file.variables['xq'][:]
			y = file.variables['yh'][:]
	
	print 'Data size:' ,data.shape
	
	if len(data.shape)==3:
		data=np.squeeze(np.sum(data,axis=1))  #Mean over first variable
		data=np.squeeze(np.mean(data,axis=1))  #Mean over first variable
	if len(data.shape)==4:
			#data=np.squeeze(np.mean(data,axis=0))  #Mean over first variable
			data=np.squeeze(np.sum(data,axis=1))  #Mean over first variable
			data=np.squeeze(np.mean(data,axis=1))  #Mean over first variable
			data=np.squeeze(np.mean(data,axis=1))  #Mean over first variable

	return [data, Time,x,y]

def plot_time_series_data(data,Time,field): 
	print 'Starting to plot...'	
	
	plt.plot(Time,data)
	plt.title(field)



#Clear screen
def main():

	#os.system('clear')
	parser = argparse.ArgumentParser()
	parser.add_argument('--file', default=None, help='The input data file in NetCDF format.')
	parser.add_argument('--field', default='u', help='The field to be viewed')
	args = parser.parse_args()

	filename=args.file
	field=args.field

	(time_series, Time) =load_timeseries_data_from_file(filename,field,x,y)

	plot_time_series_data(time_series,Time, field)




	plt.show()


	print 'Script complete'



if __name__ == '__main__':
        main()
        #sys.exit(main())
