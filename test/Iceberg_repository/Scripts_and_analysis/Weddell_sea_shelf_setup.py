#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("GTKAgg")
from pylab import *
#import matplotlib.pyplot as plt
import pdb
import netCDF4 as nc
from scipy import interpolate
from hexagon_area import Divide_hexagon_into_4_quadrants_old
from hexagon_area import Hexagon_into_quadrants_using_triangles


def Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg):
	
	print 'Writing iceberg restart files...'
	# To copy the global attributes of the netCDF file  

	#Input and output files	
	f=Dataset('input_files/icebergs.res.nc','r') # r is for read only
	g=Dataset('output_files/' + Ice_geometry_source + '_icebergs.res.nc','w') # w if for creating a file

	for attname in f.ncattrs():
		    setattr(g,attname,getattr(f,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in f.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		g.createDimension(dimname,Number_of_bergs)

	# To copy the variables of the netCDF file

	for varname,ncvar in f.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = g.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():  
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		#var[:] = ncvar[0]  #I commented out this line because it was causing errors. I'm not sure if it is needed.
		
		if varname=='i':
			var[:]=Number_of_bergs
		
		if varname=='iceberg_num':
			for j in range(Number_of_bergs):
				#var[j]=j+1
				var[j]=iceberg_num[j]

		if varname=='uvel' or varname=='vvel' or varname=='uvel_old' or varname=='vvel_old' or varname=='axn' or varname=='ayn'\
		or varname=='bxn' or varname=='byn' or  varname=='halo_berg' or varname=='heat_density' or varname=='lon_old' or varname=='lat_old' \
		or varname=='mass_of_bits' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' \
		or varname=='start_lat' or varname=='start_mass' or  varname=='start_day' or varname=='start_year' or varname=='start_lon' or varname=='lat_old':\
			var[:]=0

		if varname=='mass_scaling':
			var[:]=mass_scaling

		if varname=='thickness':
			for j in range(Number_of_bergs):
				var[j]=thickness[j]
		
		if varname=='mass':
			for j in range(Number_of_bergs):
				var[j]=mass[j]

		if varname=='width'  or varname=='length':
			for j in range(Number_of_bergs):
				var[j]=width[j]

		if varname=='lon':
			for j in range(Number_of_bergs):
				var[j]=lon[j]

		if varname=='lat':
			for j in range(Number_of_bergs):
				var[j]=lat[j]

		if varname=='static_berg':
			for j in range(Number_of_bergs):
				var[j]=static_berg[j]


	f.close()
	g.close()

	
def Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source):
	#Creating the bond restart file

	print 'Writing bond restart files...'
	# To copy the global attributes of the netCDF file  

	#Input and output files	
	h=Dataset('input_files/bonds_iceberg.res.nc','r') # r is for read only
	q=Dataset('output_files/' + Ice_geometry_source + '_bonds_iceberg.res.nc','w') # w if for creating a file

	for attname in h.ncattrs():
		    setattr(q,attname,getattr(h,attname))


	# To copy the dimension of the netCDF file
	for dimname,dim in h.dimensions.iteritems():
		# if you want to make changes in the dimensions of the new file
		# you should add your own conditions here before the creation of the dimension.
		#g.createDimension(dimname,len(dim))
		q.createDimension(dimname,Number_of_bonds)

	# To copy the variables of the netCDF file

	for varname,ncvar in h.variables.iteritems():
		# if you want to make changes in the variables of the new file
		# you should add your own conditions here before the creation of the variable.
		var = q.createVariable(varname,ncvar.dtype,ncvar.dimensions)
		#Proceed to copy the variable attributes
		for attname in ncvar.ncattrs():  
			setattr(var,attname,getattr(ncvar,attname))
		#Finally copy the variable data to the new created variable
		#var[:] = ncvar[0]
		var[:] = 0.

		if varname=='i':
			var[:]=Number_of_bonds

		if varname=='first_berg_num':
			for j in range(Number_of_bonds):
				var[j]=first_berg_num[j]

		if varname=='first_berg_ine':
			for j in range(Number_of_bonds):
				var[j]=first_berg_ine[j]

		if varname=='first_berg_jne':
			for j in range(Number_of_bonds):
				var[j]=first_berg_jne[j]

		if varname=='other_berg_num':
			for j in range(Number_of_bonds):
				var[j]=other_berg_num[j]

		if varname=='other_berg_ine':
			for j in range(Number_of_bonds):
				var[j]=other_berg_ine[j]

		if varname=='other_berg_jne':
			for j in range(Number_of_bonds):
				var[j]=other_berg_jne[j]

	h.close()
	q.close()



def Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius,x_ind_vec,y_ind_vec,h_ice,x,y,\
		width,Interpolate_from_four_corners,element_area,element_type,Find_thickness_using_vector,static_berg):
	thickness=[]
	mass=[]
	ny,nx=h_ice.shape
	dx=x[1]-x[0]
	grid_area=dx*dx
	for berg_count in range(Number_of_bergs):
		x_val=dx_berg[berg_count] ; y_val=dy_berg[berg_count]
		R_ind=(abs(xi-x_val)+(abs(yi-y_val))).argmin()
		i_val=floor(x_val/dx)
		j_val=floor(y_val/dx)

		#Interpolate thickness from 4 corners - possibly do this later, but I worry about when you are between a shelf and non shelf piece.
		if Interpolate_from_four_corners==True:

			x_cell=(x_val-x[i_val])/sqrt(grid_area)+0.5
			y_cell=(y_val-y[j_val])/sqrt(grid_area)+0.5

			mass_on_ocean=np.zeros([nx,ny,10])   #Setting up matrix to spread mass to ocean.  Note that I have used 10 points  so that I can ignore 0 and match with python numbering
			mass_val=1.
			mass_on_ocean=spread_mass_to_ocean(nx,ny,i_val,j_val,mass_on_ocean,x_cell,y_cell,element_area,mass_val,element_type,grid_area,static_berg[berg_count])
			Th=0.
			h=np.array([0,0,0,0,0,0,0,0,0,0])
			if i_val>0 and j_val>0:
				h[1]=h_ice[j_val-1,i_val-1]
			if j_val>0:
				h[2]=h_ice[j_val-1,i_val]
			if j_val>0 and i_val<(nx-1):
				h[3]=h_ice[j_val-1,i_val+1]
			if i_val>0:
				h[4]=h_ice[j_val,i_val-1]
			if True:
				h[5]=h_ice[j_val,i_val]
			if i_val<(nx-1):
				h[6]=h_ice[j_val,i_val+1]
			if i_val>0 and j_val<(ny-1):
				h[7]=h_ice[j_val+1,i_val-1]
			if j_val<(ny-1):
				h[8]=h_ice[j_val+1,i_val]
			if i_val<(nx-1) and j_val<(ny-1):
				h[9]=h_ice[j_val+1,i_val+1]
			for k in range(1,10):
				Th=Th+(mass_on_ocean[i_val,j_val,k]*h[k])
				#print i_val,j_val,k,mass_on_ocean[i_val,j_val,k]
			#if abs(Th-1)>0.0000001:
			#	print 'Thickness',Th-1,i_val,j_val
			#	test=0.
			#	for k in range(1,10):
			#		print k,mass_on_ocean[i_val,j_val,k], h[k]
			#		test=test+(mass_on_ocean[i_val,j_val,k])
			#	print test
			#	halt
				

		elif Find_thickness_using_vector==True:
			Th=h_ice_vec[R_ind][0]
		else:#MP3
			Th=h_ice[j_val,i_val]
		thickness.append(Th)
		#Check that all thicknesses are positive
		if Th<=0:
			print 'Thickness is less than or equal to zero,0!',i_val, j_val,x_val,y_val
			halt
		#mass.append(Th*rho_ice*element_area) 
		mass.append(Th*rho_ice*(width[berg_count])**2) 
		#if element_type=='square': 
		#	mass.append(Th*rho_ice*((2*Radius)**2)) # For square elements
		#else:
		#	mass.append(Th*rho_ice*np.pi*(Radius**2))  #For circular elements
	return [thickness, mass]

def Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref,adjust_lat_ref):
	N=len(x)  ; M = len(y)
	x_min=np.min(x)  ;y_min=np.min(y)
	count=-1
	#R=np.zeros([N*M,1])
	xi=np.zeros([N*M,1])
	yi=np.zeros([N*M,1])
	ice_mask_vec=np.zeros([N*M,1])
	h_ice_vec=np.zeros([N*M,1])
	x_ind_vec=np.zeros([N*M,1])
	y_ind_vec=np.zeros([N*M,1])

	#print N,M
	for i in range(N):
		for j in range(M):
			count=count+1
			if input_is_cartesian==True:
				#R[count]=sqrt(((x[i]-x_min)**2) + ((y[j]-y_min)**2))
				xi[count]=x[i]-x_min
				yi[count]=y[j]-y_min
			else:
				#x and y are given in lon and lat
				dlon=x[i]-x_min
				dlat=y[j]-y_min
				dy_dlat=(np.pi/180)*(R_earth)			
				if adjust_lat_ref==True:
					lat_ref=y[j]

				dx_dlon=(np.pi/180)*(R_earth)*(np.cos(lat_ref*np.pi/180))
				dx=dlon*dx_dlon
				dy=dlat*dy_dlat
				#R[count]=sqrt((dx**2)+(dy**2))
				xi[count]=dx
				yi[count]=dy
			ice_mask_vec[count,0]=ice_mask[j,i]
			h_ice_vec[count,0]=h_ice[j,i]
			x_ind_vec[count,0]=i
			y_ind_vec[count,0]=j
			

	ice_mask_vec[np.where(np.isnan(ice_mask_vec))]=0
	#cNorm = mpl.colors.Normalize(vmin=-1., vmax=5.)
	#plt.pcolor(x,y,ice_mask,norm=cNorm,cmap='jet')
	#plt.scatter(x_temp, y_temp,c=ice_mask_vec[:,0],cmap='jet',norm=cNorm)
	plt.show()

	return [xi,yi, ice_mask_vec,h_ice_vec,x_ind_vec,y_ind_vec]


def check_if_it_is_in_domain(x_val,y_val,x_min,x_max,y_min,y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,dx,dy):
	point_is_in_domain=True
	if input_is_cartesian==False:
		dlat_dy=(180/np.pi)*(1/R_earth)
		if adjust_lat_ref==True:
			lat_ref=lat_init+(dlat_dy*y_val)
		else:
			lat_ref=y_max 
		dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180))
		x_val=(x_val*dlon_dx)
		y_val=(y_val*dlat_dy)
	
	if (x_val >= (x_max-x_min+(dx))) or (x_val<= 0) or  (y_val >= (y_max-y_min+(dy))) or (y_val <= 0):
		point_is_in_domain=False
	#print 'x', x_val,  x_max
	#print 'y', y_val,  y_max, floor(y_val/dx), floor(y_max/dx)
	
	return point_is_in_domain

def check_if_it_is_ice(x_val,y_val,xi,yi,ice_mask_vec,ice_mask,input_is_cartesian,Find_thickness_using_vector,dx):
	#R_val=np.sqrt(((x_val)**2) +((y_val)**2))
	#R_ind=(abs(R-R_val)).argmin()
	if Find_thickness_using_vector==True:
		R_ind= (abs(x_val-xi) + abs(yi-y_val)).argmin()
		if ice_mask_vec[R_ind]==1:
			it_is_ice=True
		else:
			it_is_ice=False

	else:
		i_val=floor(x_val/dx)
		j_val=floor(y_val/dx)
		#print i_val,j_val,x_val,y_val
		if ice_mask[j_val,i_val]==1.:
			it_is_ice=True
		else:
			it_is_ice=False

	return it_is_ice

def calculate_element_area(element_type,Radius):
	if element_type=='square':
		element_area=(2*Radius)**2
	elif element_type=='hexagon':
		element_area=(3.*np.sqrt(3.)/2.)*((4./3.)*(Radius)**2) #Area of hexagon around circle (used for packing)  
		#Another derivation uses innner hexagon with two more triangles added, which is a 1/6 of the hexagon area each (two since there are 6, shared by 3 each)
		#element_area=(4./3.)*H_i, where H_i=(3.*np.sqrt(3.)/2.)*((Radius)**2)  is the area of the inner hexagon (with sides equal to the radius)

	return element_area

def add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,\
		iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice):
	tol=0.0000000000001
	if New_thickness[j,i]>1.+tol:
		print 'The new thickness is too big!!!', i,j,  New_thickness[j,i]
		halt
	if ((Nbh_thickness_value-New_thickness[j,i])>tol):
		y_val=y_shift[j]
		x_val=x_shift[i]
		if Nbh_thickness_value>0. and (abs(Nbh_thickness_value-1)<tol) and (Nbh_mass_value>New_mass[j,i]):
			berg_count=berg_count+1
			dx_berg.append(x_val)
			dy_berg.append(y_val)
			iceberg_num.append(berg_count)
			width_val=np.sqrt((Nbh_thickness_value-New_thickness[j,i])*grid_area)
			if Nbh_thickness_value-New_thickness[j,i]<0:
				print 'Stop', Nbh_thickness_value, New_thickness[j,i],i,j
			mass_val=Nbh_mass_value-New_mass[j,i]
			thickness_val=mass_val/(width_val*width_val*rho_ice)
			thickness.append(thickness_val)
			width.append(width_val)
			mass.append(mass_val)
			static_berg.append(1.)
	return [berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]


def add_extra_bergs_on_boundary(dx,dy,x,y,Number_of_bergs,element_area,rho_ice,X_min,X_max,Y_min,Y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,xi,yi,ice_mask_vec,ice_mask,\
		iceberg_num,width,dx_berg,dy_berg,h_ice,element_type,Find_thickness_using_vector,static_berg,thickness,mass):
	eps=0.000001
	grid_area=dx*dy
	Nx=len(x)  ; Ny=len(y)
	x_shift=(x-np.min(x))+(dx/2) ; y_shift=(y-np.min(y))+(dy/2)

	thickness_temp= [1. for i in range(Number_of_bergs)] ;
	mass_temp= [element_area*rho_ice*1. for i in range(Number_of_bergs)] ;
	New_area=regrid_iceberg_thickness(dy_berg,dx_berg,Number_of_bergs,thickness_temp,mass_temp,h_ice,x_shift,y_shift,rho_ice,element_type,static_berg,plot_outcome=False)
	New_thickness=(New_area)/(rho_ice*grid_area)  #Should be equal to one everywhere where there is ice shelf.
	New_mass=regrid_iceberg_thickness(dy_berg,dx_berg,Number_of_bergs,thickness,mass,h_ice,x_shift,y_shift,rho_ice,element_type,static_berg,plot_outcome=False)
	berg_count=Number_of_bergs
	
	
	#Vertical boundary
	for i in np.array([0,Nx-1]):
		#for j in range(Ny):
		for j in range(1,Ny-1):
			Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
			if i==0:
				Nbh_thickness_value=New_thickness[j,1]  #Neighbouring shelf value
				#x_val=0+(eps*dx)
			if i==Nx-1:
				Nbh_thickness_value=New_thickness[j,Nx-2]  #Neighbouring shelf value
				#x_val=x_shift[i]+(dx/2)-eps
			[berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
					Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)

	##Horizontal boundary
	for j in np.array([0,Ny-1]):
		#for i in range(1,Nx-1):
		for i in range(1,Nx-1):
			Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
			if j==0:
				Nbh_thickness_value=New_thickness[1,i]  #Neighbouring shelf value
				#y_val=0+(eps*dy)
			if j==Ny-1:
				Nbh_thickness_value=New_thickness[Ny-2,i]  #Neighbouring shelf value
				#y_val=y_shift[j]+(dy/2)-eps
			#Nbh_thickness_value=np.mean(New_thickness[range(1,Ny-1),i])  #Neighbouring shelf value
			[berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
					Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)

	#Doing the corners seperately
	for j in np.array([0,Ny-1]):
		#for i in range(1,Nx-1):
		for i in np.array([0,Nx-1]):
			Nbh_mass_value=h_ice[j,i]*rho_ice*grid_area  #Neighbouring shelf value
			if j==0 and i==0:
				Nbh_thickness_value=New_thickness[1,1]  #Neighbouring shelf value
			if j==0 and i==Nx-1:
				Nbh_thickness_value=New_thickness[1,Nx-2]  #Neighbouring shelf value
			if j==Ny-1 and i==0:
				Nbh_thickness_value=New_thickness[Ny-2,1]  #Neighbouring shelf value
			if j==Ny-1 and i==Nx-1:
				Nbh_thickness_value=New_thickness[Ny-2,Nx-2]  #Neighbouring shelf value
			[berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg]=add_and_extra_boundary_berg(i,j,New_thickness,y_shift,x_shift,\
					Nbh_thickness_value,Nbh_mass_value,New_mass,berg_count,dx_berg,dy_berg,iceberg_num,thickness,width,mass,static_berg,grid_area,rho_ice)


	Number_of_bergs=berg_count
	
	return [dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg,thickness,mass]


def remove_stationary_bergs(dx_berg, dy_berg,iceberg_num,width,Number_of_bergs,static_berg):
	dx_berg_new=[]
	dy_berg_new=[]
	iceberg_num_new=[]
	width_new=[]
	static_berg_new=[]
	for i in range(Number_of_bergs):
		if static_berg[i]<0.5:
			dx_berg_new.append(dx_berg[i])
			dy_berg_new.append(dy_berg[i])
			iceberg_num_new.append(iceberg_num[i])
			width_new.append(width[i])
			static_berg_new.append(static_berg[i])

	Number_of_bergs_new=len(iceberg_num_new)
	return [dx_berg_new, dy_berg_new,iceberg_num_new, width_new,Number_of_bergs_new,static_berg_new]

def Create_icebergs(lon_init,lat_init,Radius,R_earth, x, y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian,\
		element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,Interpolate_from_four_corners,\
		Fill_in_the_boundaries,set_all_bergs_static_by_default,break_some_bonds,Find_thickness_using_vector,Remove_stationary_bergs):
	print 'Starting to create icebergs...'
	dx_berg=[]  #x distance in cartesian of berg from lon_init
	dy_berg=[]  #y distance in cartesian of berg from lat_init
	dlon_berg=[] #x distance in lon of berg from lon_init
	#dlat_berg=[] #y distance in lat of berg from lat_init
	lon=[] #Longitude of iceberg
	lat=[] #Latitude of iceberg
	iceberg_num=[] #ID of iceberg
	
	element_area=calculate_element_area(element_type,Radius)
	#width=np.sqrt(element_area)
	width=[]
	
	#Create a vector of distances from a minimum point.
	(xi,yi,ice_mask_vec,h_ice_vec,x_ind_vec,y_ind_vec)=Create_distance_and_ice_mask_vector(x,y,ice_mask,h_ice,input_is_cartesian,R_earth,lat_ref,adjust_lat_ref)
	
	#dx=x[1]-x[0]  ;dy=y[1]-y[0]
	#X_min=np.min(x)-dx; X_max=np.max(x)+dx  #In lat lon or cartesian
	#Y_min=np.min(y)-dy; Y_max=np.max(y)+dy  #In lat lon or cartesian
	X_min=np.min(x); X_max=np.max(x)  #In lat lon or cartesian
	Y_min=np.min(y); Y_max=np.max(y)  #In lat lon or cartesian

	#N=2*int(ceil((R_max)/(2*Radius))+2)
	
	if element_type=='square':
		#N=int(ceil((np.max(xi))/(2*Radius))+2) +2
		#M=int(ceil((np.max(yi))/(2*Radius))+2) +2
		N=int(ceil((X_max-X_min)/(Radius)))
		M=int(ceil((Y_max-Y_min)/(Radius)))
	else:
		#N=int(ceil((np.max(xi))/(Radius/2))+2) +2
		#M=int(ceil((np.max(yi))/((1/sqrt(3))*Radius/2))+2) +2
		N=2*int(ceil((X_max-X_min)/(Radius)))
		M=2*int(ceil((Y_max-Y_min)/(Radius)))
	#N=10
	#M=6
	#MP4

 	dx=x[1]-x[0]	
 	dy=y[1]-y[0]
	Lx=X_max+dx
	Ly=Y_max+dx

	berg_count=0
	#for j in range(N):
	for i in range(N):
		if element_type=='square':
			x_start=Radius
			y_start=(Radius)
			x_val=x_start+(2*i*Radius)

		#Hexagonal	
		else:   
			y_start=(Radius)+(((i%2)*Radius))
			x_start=((2/sqrt(3))*Radius)
			x_val=x_start + (np.sqrt(3)*Radius*i)

		for j in range(M):
		#for i in range(M):
			#x_val=x_start+(2*i*Radius)  ; y_val=y_start
			y_val=y_start+(2*j*Radius)  ;# x_val=x_start + (np.sqrt(3)*Radius*i)
			if check_if_it_is_in_domain(x_val,y_val,X_min,X_max,Y_min,Y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,dx,dy):
			#if True:
				#R_val=np.sqrt(((x_val-x0)**2) +((y_val-y0)**2))
				if check_if_it_is_ice(x_val,y_val,xi,yi,ice_mask_vec,ice_mask,input_is_cartesian,Find_thickness_using_vector,dx):
					#Don't allow points closer than R from the boundary (these are sorted out later)
					#if abs(y_val-Ly)<Radius or  ((abs(x_val-Lx)<((2/sqrt(3))*Radius)) and element_type=='hexagon') or ((abs(x_val-Lx)<Radius) and element_type=='square'):
					berg_count=berg_count+1
					iceberg_num.append(berg_count)
					dx_berg.append(x_val)
					dy_berg.append(y_val)
					width.append(np.sqrt(element_area))

	Number_of_bergs=berg_count
	print 'Icebergs created. Number of bergs = ', Number_of_bergs
	
	#Deciding if icebergs are static or not
	static_berg = [0. for i in dx_berg]
	if set_all_bergs_static_by_default==True:
		static_berg = [1. for i in static_berg]
	if break_some_bonds==True:
		#static_berg =Change_static_berg_after_calving(Number_of_bergs,lat,lon, static_berg)
		static_berg =Change_static_berg_after_calving(Number_of_bergs,dy_berg,dx_berg, static_berg)
		static_berg=static_berg[0]
		
	if Remove_stationary_bergs is True:
		Fill_in_the_boundaries=False
		[dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg]=remove_stationary_bergs(dx_berg, dy_berg,iceberg_num,\
				width,Number_of_bergs,static_berg)
		print 'After purge  Number of bergs = ', Number_of_bergs


	#Defining the thickness of the icebergs
	(thickness, mass)=Define_iceberg_thickness_and_mass(Number_of_bergs,dx_berg,dy_berg,h_ice_vec,xi,yi,rho_ice,Radius,x_ind_vec,y_ind_vec, h_ice,x,y,\
			width,Interpolate_from_four_corners,element_area,element_type,Find_thickness_using_vector,static_berg)	

	N_bergs_before_bd=Number_of_bergs
	if Fill_in_the_boundaries==True and element_type=='hexagon':
	 	[dx_berg, dy_berg,iceberg_num, width,Number_of_bergs,static_berg,thickness,mass]= add_extra_bergs_on_boundary(dx,dy,x,y,Number_of_bergs,element_area,rho_ice,X_min,X_max,Y_min,\
				Y_max,input_is_cartesian,R_earth,lat_init,adjust_lat_ref,xi,yi,ice_mask_vec,ice_mask,iceberg_num,width,dx_berg,dy_berg,\
				h_ice,element_type,Find_thickness_using_vector,static_berg,thickness,mass)
		print 'Number of icebergs after accounting for boundaries = ', Number_of_bergs
	


	if Convert_to_lat_lon==True:
		#Defining lon lat positions:
		#dlat_berg=(180/np.pi)*(1/R_earth)*dy_berg
		#dlon_berg=(180/np.pi)*(1/R_earth)*(1/cos(lat_init*np.pi/180))*dx_berg
		for i in range(Number_of_bergs):
			#Finding latittude
			dlat_dy=(180/np.pi)*(1/R_earth)

			lat.append(lat_init+(dlat_dy*dy_berg[i]))

			#Finding longitude
			if adjust_lat_ref==True:
				lat_ref=lat_init+(dlat_dy*dy_berg[i])
			dlon_dx=(180/np.pi)*(1/R_earth)*(1/np.cos(lat_ref*np.pi/180)) #Note that this depends on the latitude of the iceberg. Could approx this with lat_init.
			lon.append(lon_init+(dlon_dx*dx_berg[i] ))

			#dlon_berg.append(dlon_dx*dx_berg[i])
			#lon.append(lon_init+dlon_berg[i])
	else:
		if scale_the_grid_to_lat_lon==True:
			Scale_up=1./2000.
			Radius=Radius*Scale_up
			dx_berg = [(i*Scale_up) for i in dx_berg] ; dy_berg = [i*Scale_up for i in dy_berg] 
			x=x*Scale_up ; y=y*Scale_up
			dx=dx*Scale_up ;dy=dy*Scale_up
		x=(x-np.min(x))+(dx/2) ; y=(y-np.min(y))+(dy/2)
		lon=dx_berg  ; lat=dy_berg
	
	#Note that static_berg calculations used to be here, after the conversion. I have moved them. I hope that this does not affect answers.

	return (Number_of_bergs,lon,lat,iceberg_num,dx_berg,dy_berg,thickness, mass,width,x,y,Radius,static_berg,N_bergs_before_bd)



def Create_calving_event(lat1,lon1,lat2,lon2):
	[R_calve, Calve_lon, Calve_lat]=get_calving_parameters()
	
	R1=np.sqrt((lon1-Calve_lon)**2+ (lat1-Calve_lat)**2)
	R2=np.sqrt((lon2-Calve_lon)**2+ (lat2-Calve_lat)**2)

	bond_broken=False
	if ((R1 < R_calve)*(R2>R_calve) )  > 0.5 :
		bond_broken=True
	if ((R2 < R_calve)*(R1>R_calve) )  > 0.5 :
		bond_broken=True

	return bond_broken

def get_calving_parameters():
	#Calve_lat =20  
	#Calve_lon=160
	#R_calve=15
	Calve_lat =20.2*2000  
	Calve_lon=160*2000
	R_calve=15*2000
	#R_calve=1.*2000.
	
	return [R_calve, Calve_lon, Calve_lat]

def Change_static_berg_after_calving(Number_of_bergs,lat,lon, static_berg):
	[R_calve, Calve_lon, Calve_lat]=get_calving_parameters()
	count=0.
	for i in range(Number_of_bergs):
		R1=np.sqrt((lon[i]-Calve_lon)**2+ (lat[i]-Calve_lat)**2)
		#Making calved icebergs not be static.
		if R1<R_calve:
			count=count+1.
			#print 'An iceberg is now static!'
			static_berg[i]=0.
	print 'Amount of unstatic icebergs after calving = ', count

	return [static_berg]

def find_max_number_of_bonds(first_berg_num,Number_of_bonds):
	count=1
	best_count=1
	sorted_bonds=np.sort(first_berg_num)
	for k in range(1,Number_of_bonds):
		previous_berg=sorted_bonds[k-1]
		print k, sorted_bonds[k]
		if sorted_bonds[k] == sorted_bonds[k-1]:
			count=count+1
			best_count=max(count,best_count)
		else:
			count=1
	Max_number_of_bonds=best_count
	return Max_number_of_bonds

def Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds,static_berg, Allow_bonds_for_static_iceberg\
		,N_bergs_before_bd,Allow_bonds_with_boundary_bergs):
	print 'Starting to create bonds...'
	#Defining Bonds:
	Bond=np.zeros((Number_of_bergs, Number_of_bergs))
	bond_broken=False
	first_berg_num=[]  # Initializing bond list first berg
	first_berg_ine=[]  # Initializing bond list first berg
	first_berg_jne=[]  # Initializing bond list first berg
	first_berg_lat=[]  # Initializing bond list first berg
	first_berg_lon=[]  # Initializing bond list first berg
	other_berg_num=[]  # Initializing bond list other berg
	other_berg_ine=[]  # Initializing bond list other berg
	other_berg_jne=[]  # Initializing bond list other berg
	other_berg_lat=[]  # Initializing bond list other berg
	other_berg_lon=[]  # Initializing bond list other berg
	bond_count=0
	for i in range(Number_of_bergs):
		if (static_berg[i]<0.5) or (Allow_bonds_for_static_iceberg is True):
			#for j in range(Number_of_bergs):
			for j in range(i):
				if (static_berg[j]<0.5) or (Allow_bonds_for_static_iceberg is True):
					if i!=j:
						if ((i <N_bergs_before_bd) and (j <N_bergs_before_bd))  or (Allow_bonds_with_boundary_bergs is True):
							R_dist=np.sqrt(((dx_berg[i]-dx_berg[j])**2) + ((dy_berg[i]-dy_berg[j])**2))
							if break_some_bonds==True:
								bond_broken=Create_calving_event(lat[i],lon[i],lat[j],lon[j])
							if R_dist < (2.01*Radius) and (bond_broken==False):
								bond_count=bond_count+2
								#Connect bond in the first direction
								first_berg_num.append(iceberg_num[i]);	first_berg_ine.append(999); 
								other_berg_num.append(iceberg_num[j]); 	other_berg_ine.append(999);
								first_berg_jne.append(999);first_berg_lat.append(lat[i]); first_berg_lon.append(lon[i])
								other_berg_jne.append(999);other_berg_lat.append(lat[j]); other_berg_lon.append(lon[j])
								#Connect bond in the other direction
								first_berg_num.append(iceberg_num[j]);	first_berg_ine.append(999);
								other_berg_num.append(iceberg_num[i]); 	other_berg_ine.append(999);
								first_berg_jne.append(999);first_berg_lat.append(lat[j]); first_berg_lon.append(lon[j])
								other_berg_jne.append(999);other_berg_lat.append(lat[i]); other_berg_lon.append(lon[i])
	Number_of_bonds=bond_count
	Max_number_of_bonds=find_max_number_of_bonds(first_berg_num,Number_of_bonds)
	print 'Number of bonds created = ' , Number_of_bonds
	print 'Maximum number of bonds = ' ,  Max_number_of_bonds

	return [ Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon]

def load_ISOMIP_ice_geometry(filename,buffer_number):
	with nc.Dataset(filename) as file:
		ocean_mask = file.variables['openOceanMask'][:,:]
		upperSurface = file.variables['upperSurface'][:,:]
		lowerSurface = file.variables['lowerSurface'][:,:]
		x = file.variables['x'][:]
		y = file.variables['y'][:]

	ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
	h_ice=upperSurface-lowerSurface #The ice thickness
	M=ice_mask.shape

	#Setting the boundaries to non-ice
	A=np.arange(0,buffer_number)
	B=np.arange(M[0]-buffer_number,M[0])
	C=np.arange(M[1]-buffer_number,M[1])
	ice_mask[A,:]=0; ice_mask[B,:]=0
	ice_mask[:,A]=0; ice_mask[:,C]=0
	
	return [x,y,ice_mask,h_ice]


def create_clipped_icethickness_file(h_ice,area,mass,grid_area,gravity):
	#Creating clipped file
        [ny, nx]= h_ice.shape ; 
	
	New_ice_thickness_filename='output_files/isomip_ice_shelf1_clipped.nc'
        g=Dataset(New_ice_thickness_filename,'w') # w if for creating a file

        g.createDimension('nx',nx)
        g.createDimension('ny',ny)

        thick_h=g.createVariable('thick','f8',('ny','nx'))
        area_h=g.createVariable('area','f8',('ny','nx'))
        p_surf_h=g.createVariable('p_surf','f8',('ny','nx'))

        thick_h.units = 'm'
        thick_h.standard_name = 'ice shelf thickness (clipped)'
        area_h.units = 'm2'
        area_h.standard_name = 'ice shelf area'
        p_surf_h.units = 'Pa'
        p_surf_h.standard_name = 'surface pressure due to ice shelf'

	p_surf=(gravity*mass)/grid_area
        g.variables['thick'][:]=h_ice
        g.variables['area'][:]=area
        g.variables['p_surf'][:]=p_surf
        print 'Creating clipped ice file: ' , New_ice_thickness_filename
        
        g.sync()
        g.close()


def load_ISOMIP_reduced_ice_geometry(ice_filename,buffer_number,topog_filename):
	with nc.Dataset(ice_filename) as file:
		h_ice = file.variables['thick'][:,:]
		area =  file.variables['area'][:,:]
		M=h_ice.shape

	y=np.linspace(1000,79000,M[0],endpoint=True)
	#x=np.linspace(321000,799000,M[1],endpoint=True)
	x=np.linspace(1000,479000,M[1],endpoint=True)
	ice_mask=h_ice>0.


	count=0. #MP1
	#Setting the boundaries to non-ice
	A=np.arange(0,buffer_number)
	B=np.arange(M[0]-buffer_number,M[0])
	C=np.arange(M[1]-buffer_number,M[1])
	#if buffer_number>0:
	ice_mask[A,:]=0; ice_mask[B,:]=0
	ice_mask[:,A]=0; ice_mask[:,C]=0
	
	return [x,y,ice_mask,h_ice]

def load_Weddel_ice_geometry(filename):
	with nc.Dataset(filename) as file:
		h_ice = file.variables['thickness'][:,:]
		ice_mask = file.variables['icemask_grounded_and_shelves'][:,:]
		x = file.variables['lon'][:]
		y = file.variables['lat'][:]
		ice_mask[np.where(ice_mask<-1)]=0.
		
		h_ice[np.where(np.isnan(h_ice))]=0.
		ice_mask[np.where(np.isnan(ice_mask))]=0.

		#ice_mask=1-ocean_mask #one if it ice, zero if it is ocean
	
	return [x,y,ice_mask,h_ice]

def Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num,static_berg):
	print 'You have chosen to choose just one icebergs!!!'
	for k in range(len(lat)):
		if iceberg_num[k]==chosen_berg_num:
			berg_ind=k
	lon_temp=lon[berg_ind]; lon=[] ; lon.append(lon_temp)
	lat_temp=lat[berg_ind]; lat=[] ; lat.append(lat_temp)
	thickness_temp=thickness[berg_ind]; thickness=[] ; thickness.append(thickness_temp)
	mass_temp=mass[berg_ind]; mass=[] ; mass.append(mass_temp)
	width_temp=width[berg_ind]; width=[] ; width.append(width_temp)
	static_berg_temp=static_berg[berg_ind]; static_berg=[] ; static_berg.append(static_berg_temp)
	iceberg_num_temp=iceberg_num[berg_ind]; iceberg_num=[] ; iceberg_num.append(iceberg_num_temp)
	Number_of_bergs=1
	return [Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num]

def plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon, \
		plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions,static_berg):
	print 'Starting to plot...'	
	Radius=Radius*IA_scaling
	circ_ind=np.linspace(0,2*pi,100);

	if plot_ice_mask==True:
		cNorm = mpl.colors.Normalize(vmin=0., vmax=1.)
		plt.pcolormesh(x,y,ice_mask,norm=cNorm)
	elif plot_ice_thickness==True:
		cNorm = mpl.colors.Normalize(vmin=0., vmax=1000.)
		plt.pcolor(x,y,h_ice,cmap='jet',norm=cNorm)


	if plot_icebergs_positions==True:
		#plt.scatter(lon, lat,color='yellow')
		#cNorm = mpl.colors.Normalize(vmin=0., vmax=1000.)
		plt.scatter(lon, lat,c=thickness,norm=cNorm,cmap='jet',s=150)
		cNorm = mpl.colors.Normalize(vmin=-1, vmax=1.)
		#plt.scatter(lon, lat,c=static_berg,norm=cNorm,cmap='jet',s=150)

	#plt.plot(lon, lat,'bo-',linewidth=5)
	if plot_circles==True:
		for k in range(Number_of_bergs):
			if Convert_to_lat_lon==True:
				dR_lat=(Radius/R_earth)*(180/np.pi)
				dR_lon=(Radius/R_earth)*(180/np.pi)*(1/np.cos(lat[k]*np.pi/180))
				plt.plot(lon[k]+(dR_lon*cos(circ_ind)),lat[k]+(dR_lat*sin(circ_ind)),'b');
			else:
				plt.plot(lon[k]+(Radius*cos(circ_ind)),lat[k]+(Radius*sin(circ_ind)),'b');



	#plt.plot(lon, lat,'bo-')
	if Convert_to_lat_lon==True:
		plt.xlabel('longitude (deg)')
		plt.ylabel('latitude (deg)')
	else:
		plt.xlabel('longitude (m)')
		plt.ylabel('latitude (m)')
	plt.title('Iceberg initial positions')
	plt.grid(True)
 	plt.colorbar()


def plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds):

	for k in range(Number_of_bonds):
		x_bond=[]
		y_bond=[]
		x_bond.append(first_berg_lon[k])
		x_bond.append(other_berg_lon[k])
		y_bond.append(first_berg_lat[k])
		y_bond.append(other_berg_lat[k])
		plt.plot(x_bond, y_bond,'r',linewidth=5)


def spread_mass_to_ocean(Nx,Ny,i,j,mass_on_ocean,x,y,Area,Mass,element_type,grid_area,static_berg):
		#Note that the x,y coming into this routine are the position within a cell (from 0 to 1), with 0.5,0.5 being in the center of the cell.

		#Initialize weights for each cell	
		yDxL=0.  ; yDxC=0. ; yDxR=0. ; yCxL=0. ; yCxR=0. 
		yUxL=0.  ; yUxC=0. ; yUxR=0. ; yCxC=1.

		if element_type=='square':
		#if True:
			L = min(( (np.sqrt(Area)  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim length of the square)
			xL=min(0.5, max(0., 0.5-(x/L)))
			xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
			xC=max(0., 1.-(xL+xR))
			yD=min(0.5, max(0., 0.5-(y/L)))
			yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
			yC=max(0., 1.-(yD+yU))

			yDxL=yD*xL#*grd%msk[i-1,j-1]
			yDxC=yD*xC#*grd%msk[i  ,j-1]
			yDxR=yD*xR#*grd%msk[i+1,j-1]
			yCxL=yC*xL#*grd%msk[i-1,j  ]
			yCxR=yC*xR#*grd%msk[i+1,j  ]
			yUxL=yU*xL#*grd%msk[i-1,j+1]
			yUxC=yU*xC#*grd%msk[i  ,j+1]
			yUxR=yU*xR#*grd%msk(i+1,j+1]
			yCxC=1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )

		if element_type=='hexagon':
		#if False:#element_type=='hexagon':
			H = min(( (np.sqrt(Area/(2*sqrt(3)))  / np.sqrt(grid_area))),1) ;  #Non dimensionalize element length by grid area. (This gives the non-dim Apothen of the hexagon)
			S=(2/np.sqrt(3))*H #Side of the hexagon
			if S>0.5:
				print 'Elements must be smaller than a whole gridcell', 'i.e.: S= ' , S , '>=0.5',i,j
				halt

			#Subtracting the position of the nearest corner from x,y
			origin_x=1 ; origin_y=1
			if x<0.5:
				origin_x=0
			if y<0.5:
				origin_y=0

			x0=(x-origin_x) #Position of the hexagon center, relative to origin at the nearest vertex
			y0=(y-origin_y)
			
			#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
			(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,0.)
			if min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4)) <0:
				print 'Yolo'
				print x0,y0,H
				#print min(min(Area_Q1,Area_Q2),min(Area_Q3, Area_Q4))
				print Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4

			Area_Q1=Area_Q1/Area_hex
			Area_Q2=Area_Q2/Area_hex
			Area_Q3=Area_Q3/Area_hex
			Area_Q4=Area_Q4/Area_hex


			#Now, you decide which quadrant belongs to which mass on ocean cell.
			#Top right vertex
			if x>=0.5 and y>= 0.5:
				yUxR=Area_Q1
				yUxC=Area_Q2
				yCxC=Area_Q3
				yCxR=Area_Q4
			
			#Top left vertex
			if x<0.5 and y>= 0.5:
				yUxC=Area_Q1
				yUxL=Area_Q2
				yCxL=Area_Q3
				yCxC=Area_Q4
			
			#Bottom left vertex
			if x<0.5 and y< 0.5:
				yCxC=Area_Q1
				yCxL=Area_Q2
				yDxL=Area_Q3
				yDxC=Area_Q4

			#Bottom right vertex
			if x>=0.5 and y< 0.5:
				yCxR=Area_Q1
				yCxC=Area_Q2
				yDxC=Area_Q3
				yDxR=Area_Q4

			#if Sector<4 and Sector>-1:
			#	print Sector
				
					
			#Check that this is true
			if abs(yCxC-(1.-( ((yDxL+yUxR)+(yDxR+yUxL)) + ((yCxL+yCxR)+(yDxC+yUxC)) )))>0.001:
				print 'All the mass is not being used!!!'
				#print W1 , W2 , W3 , W4 , W5 , W6
				#print Area_Upper, Area_Lower, Area_right, Area_left
				print 'Areas: ',Area_hex,Area_hex*Area_Q1, Area_hex*Area_Q2, Area_hex*Area_Q3, Area_hex*Area_Q4
				print 'x0=',x0, 'y0=',y0, 'H=', H
				print 'Total area= ',(Area_Q1+Area_Q2+Area_Q3+Area_Q4)#, Sector
				

		#Accounting for masked points 
		a=1. ; b=1. ; c=1. ; d=1.
		if i==0:
			a=0.;
		if j==0: 
			b=0.;
		if i==Nx-1:
			c=0;
		if j==Ny-1:
			d=0;

		fraction_used= ((yDxL*a*b) + (yDxC*b) + (yDxR*b*c) +(yCxL*a) + (yCxR*c) + (yUxL*a*d) + (yUxC*d) + (yUxR*c*d) + (yCxC))
		if static_berg>0.5:
			fraction_used=1.

		mass_on_ocean[i,j,1]=mass_on_ocean[i,j,1]+(a*b*yDxL*Mass/fraction_used)
		mass_on_ocean[i,j,2]=mass_on_ocean[i,j,2]+(b*yDxC*Mass/fraction_used)
		mass_on_ocean[i,j,3]=mass_on_ocean[i,j,3]+(b*c*yDxR*Mass/fraction_used)
		mass_on_ocean[i,j,4]=mass_on_ocean[i,j,4]+(a*yCxL*Mass/fraction_used)
		mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+(yCxC*Mass/fraction_used)
		mass_on_ocean[i,j,6]=mass_on_ocean[i,j,6]+(c*yCxR*Mass/fraction_used)
		mass_on_ocean[i,j,7]=mass_on_ocean[i,j,7]+(a*d*yUxL*Mass/fraction_used)
		mass_on_ocean[i,j,8]=mass_on_ocean[i,j,8]+(d*yUxC*Mass/fraction_used)
		mass_on_ocean[i,j,9]=mass_on_ocean[i,j,9]+(c*d*yUxR*Mass/fraction_used)


		return mass_on_ocean


def calc_xiyj(x1, x2, x3, x4, y1, y2, y3, y4, x, y, Lx):

	Lx_2=Lx/2.
  
	alpha=x2-x1
	delta=y2-y1
	beta=x4-x1
	epsilon=y4-y1
	gamma=(x3-x1)-(alpha+beta)
	kappa=(y3-y1)-(delta+epsilon)

	a=(kappa*beta-gamma*epsilon)
	dx=np.mod(x-(x1-Lx_2),Lx)+(x1-Lx_2)-x1
	dy=y-y1
	b=(delta*beta-alpha*epsilon)-(kappa*dx-gamma*dy)
	c=(alpha*dy-delta*dx)
	
	#print 'alpha,beta,gamma',alpha,beta,gamma
	#print 'delta,epsilon,kappa',delta,epsilon,kappa
	#print 'dx,dy',dx,dy
	#print 'A,B,C,',a,b,c

	if (abs(a)>1.e-12):
		d=0.25*(b**2)-a*c
		if (d>=0.):
			yy1=-(0.5*b+sqrt(d))/a
			yy2=-(0.5*b-sqrt(d))/a
			if (abs(yy1-0.5).lt.abs(yy2-0.5)):
				yj=yy1;
			else:
				yj=yy2;
    		else:
			print 'ERROR', halt
	else:
		if (b!=0.):
			yj=-c/b
		else:
			yj=0.

	a=(alpha+gamma*yj)
	b=(delta+kappa*yj)
	if (a!=0.):
		xi=(dx-beta*yj)/a
	elif (b!=0.):
		xi=(dy-epsilon*yj)/b
	else:
		c=(epsilon*alpha-beta*delta)+(epsilon*gamma-beta*kappa)*yj
		if (c!=0.):
			xi=(epsilon*dx-beta*dy)/c
		else:
			print 'ERROR', halt2

	return xi,yj



def regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type,static_berg,plot_outcome=True):
	Nx=len(x)  ; Ny=len(y)
	dx=x[1]-x[0]    ;dy=y[1]-y[0]
	New_mass=h_ice*0 #Initializing regrided field to zero
	Area = [(mass[i]/(rho_ice*thickness[i])) for i in range(Number_of_bergs)] ;
	grid_area=dx*dy #Assuming a regular grid
	Orig_mass=h_ice*rho_ice*grid_area
	
	mass_on_ocean=np.zeros([Nx,Ny,10])   #Setting up matrix to spread mass to ocean.  Note that I have used 10 points  so that I can ignore 0 and match with python numbering
	#mass_on_ocean=np.zeros([Nx,Ny])   #Setting up matrix to spread mass to ocean.
	#Note: You may want to import the grid that the ocean model actually sees.
	for berg_count in range(Number_of_bergs):  
		x_val=lon[berg_count]  
		y_val=lat[berg_count]
		Area_val=Area[berg_count]
		mass_val=mass[berg_count]

		j_val=floor(y_val/dy)
		i_val=floor(x_val/dx)
		x_cell=(x_val-x[i_val])/sqrt(grid_area)+0.5
		y_cell=(y_val-y[j_val])/sqrt(grid_area)+0.5
		mass_on_ocean=spread_mass_to_ocean(Nx,Ny,i_val,j_val,mass_on_ocean,x_cell,y_cell,Area_val,mass_val,element_type,grid_area,static_berg[berg_count])

	#Adding mass_onto_ocean
	for i in range(1,Nx-1):
		for j in range(1,Ny-1):
			New_mass[j,i]=mass_on_ocean[i,j,5]  \
					+  ( ( (mass_on_ocean[i-1,j-1,9]+mass_on_ocean[i+1,j+1,1])   \
					+  (mass_on_ocean[i+1,j-1,7]+mass_on_ocean[i-1,j+1,3]) ) \
					+   ( (mass_on_ocean[i-1,j  ,6]+mass_on_ocean[i+1,j  ,4])   \
					+  (mass_on_ocean[i  ,j-1,8]+mass_on_ocean[i  ,j+1,2]) ) )
	
	#Adding mass to boundary cells
	for j in range(1,Ny-1):
		i=0
		New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i+1,j+1,1] +  mass_on_ocean[i+1,j-1,7] \
			+  mass_on_ocean[i+1,j  ,4] 	+  mass_on_ocean[i  ,j-1,8] + mass_on_ocean[i  ,j+1,2] 
		i=Nx-1	
		New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i-1,j-1,9] + mass_on_ocean[i-1,j+1,3] \
					+   mass_on_ocean[i-1,j ,6] + mass_on_ocean[i ,j-1,8]+mass_on_ocean[i ,j+1,2]
	
	for i in range(1,Nx-1):
		j=0
		New_mass[j,i]=mass_on_ocean[i,j,5] +   mass_on_ocean[i+1,j+1,1] +  mass_on_ocean[i-1,j+1,3]  \
				+  mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i+1,j  ,4] +  mass_on_ocean[i  ,j+1,2]

		j=Ny-1	
		New_mass[j,i]=mass_on_ocean[i,j,5] +  mass_on_ocean[i-1,j-1,9] 	+  mass_on_ocean[i+1,j-1,7] + \
				mass_on_ocean[i-1,j  ,6]+mass_on_ocean[i+1,j  ,4] +  mass_on_ocean[i  ,j-1,8]

	#Corners
	i=0 ;j=0
	New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i+1,j+1,1] + mass_on_ocean[i+1,j  ,4] + mass_on_ocean[i  ,j+1,2]  
	i=0 ;j=Ny-1
	New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i+1,j-1,7] + mass_on_ocean[i+1,j  ,4] + mass_on_ocean[i  ,j-1,8]  
	i=Nx-1; j=0
	New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i-1,j+1,3] + mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i  ,j+1,2]
	i=Nx-1; j=Ny-1	
	New_mass[j,i]=mass_on_ocean[i,j,5] + mass_on_ocean[i-1,j-1,9] + mass_on_ocean[i-1,j  ,6] + mass_on_ocean[i  ,j-1,8]


	if plot_outcome==True: #MP2
		#print 'Warning: the interpolation scheme is not the same as the one used in the model'	
		
		#plot_data=(Orig_mass-New_mass)/Orig_mass
		#cNorm = mpl.colors.Normalize(vmin=-.1, vmax=.1)
		#vmax=np.max(abs(plot_data))
		#plt.pcolormesh(x,y,plot_data,cmap='bwr',norm=cNorm)
		#plt.title('Error between bergs and shelf')

		#plot_data=h_ice
		plot_data=New_mass
		vmax=np.max(abs(plot_data))
		cNorm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
		plt.pcolormesh(x,y,plot_data,cmap='jet',norm=cNorm)
		plt.title('h_ice')

		#plot_data=(New_mass-Orig_mass)/(rho_ice*grid_area)
		#cNorm = mpl.colors.Normalize(vmin=-0.00000000001, vmax=0.00000000001)
		#vmax=np.max(abs(plot_data))
		#plt.pcolormesh(x,y,plot_data,cmap='bwr',norm=cNorm)

		plt.colorbar()
		#plt.xlim([min(x)-,max(x)])
		#plt.ylim([min(y),max(y)])
		#plt.show()

	return New_mass

def Switch_x_and_y_directions(lat,lon, A, B, C):
	print 'Switching x and y directions!!!'
	lat_new=[]
	lon_new=[]
	for i in range(len(lat)):
		lat_new.append(lon[i])
		lon_new.append(lat[i])
	M=A.shape
	A_new=np.zeros([M[1],M[0]])
	B_new=np.zeros([M[1],M[0]])
	C_new=np.zeros([M[1],M[0]])
	for i in range(M[0]):
		for j in range(M[1]):
			A_new[j,i]=A[i,j]
			B_new[j,i]=B[i,j]
			C_new[j,i]=C[i,j]

	return	[lat_new,lon_new,A_new,B_new, C_new]

def interpolate_field_onto_new_grid(x,y,x_new,y_new,field):
	xx, yy = np.meshgrid(x, y)
	#f = interpolate.interp2d(x, y, field, kind='cubic')
	f = interpolate.interp2d(x, y, field, kind='linear')
	field_new = f(x_new, y_new)
	return field_new

def convert_input_to_catesian_coordinates(lon,lat,ice_mask,h_ice,R_earth,dx):
	#Parameters
	thickness_threshhold=0.01

	lon_max=np.max(lon) ; lon_min=np.min(lon)
	lat_max=np.max(lat) ; lat_min=np.min(lat)
	lat_ref=np.mean(lat)

	#Converting to cartesian grid
        x=(np.pi/180)*R_earth*np.cos(lat_ref*np.pi/180.)*lon
	y=(np.pi/180)*R_earth*lat
	x=x-np.min(x)
	y=y-np.min(y)
	
	#Creating regular grid
	x_new=np.linspace(np.min(x), np.max(x), num=len(x))
	y_new=np.linspace(np.min(y),np.max(y), num=len(y))
	x_new=np.arange(np.min(x), np.max(x), step=dx)
	y_new=np.arange(np.min(y), np.max(y), step=dx)

	h_ice_new=interpolate_field_onto_new_grid(x,y,x_new,y_new,h_ice)
	h_ice_new[np.where(h_ice_new<thickness_threshhold)]=0.
	ice_mask_new=interpolate_field_onto_new_grid(x,y,x_new,y_new,ice_mask)
	#ice_mask_new=np.round(ice_mask_new)
	#print ice_mask_new[50,:]
	
	return [x_new,y_new,ice_mask_new,h_ice_new]

####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#Flags
	save_restart_files=False
	save_new_h_ice_file=True
	Convert_to_lat_lon=False
	input_is_cartesian=True

	#Iceberg setup flags
	only_choose_one_berg=False  ; chosen_berg_num=4564
	scale_the_grid_to_lat_lon=False  #Remember to change this back when providing fields for Isomip
	adjust_lat_ref=True
	set_all_thicknesses_to_one=False   ; Th_prescribed=1.
	Interpolate_from_four_corners=False
	set_all_bergs_static_by_default=True
	Find_thickness_using_vector=False
	Switch_x_and_y_to_rotate_90=False
	set_all_domain_to_ice=False

	#Mass spreadng flags
	Switch_regridding_element_type=False
	Fill_in_the_boundaries=False
	regrid_icebergs_onto_grid=True

	#Plotting flags
	plotting_bergs_over_grid=True
	plot_circles=False
	plot_ice_mask=True
	plot_ice_thickness=True
	plot_icebergs_positions=True

	#Bond related flags
	Create_icebergs_bonds=False
	plot_bonds=True
	break_some_bonds=False
	Remove_stationary_bergs=False 
	Allow_bonds_for_static_iceberg=False
	Allow_bonds_with_boundary_bergs=False

	element_type='square' 
	#element_type='hexagon'

	#Which experiment
	#Ice_geometry_source='ISOMIP'   ; Convert_to_lat_lon=False       ; input_is_cartesian=True ; 
	Ice_geometry_source='Weddell'     ; Convert_to_lat_lon=False        ; input_is_cartesian=False

	ISOMIP_reduced=True  #Reduced uses 2X2 grid, not reduced uses 1X1 grid

	ISOMIP_ice_geometry_filename='input_files/Isomip_ice_geometry.nc'
	ISOMIP_reduced_ice_geometry_filename='input_files/isomip_ice_shelf1.nc'
	Weddell_ice_geometry_filename='input_files/Bedmap2_gridded_subset_Weddell_Sea_Region.nc'
	ISOMIP_topography_filename='input_files/Isomip_topog.nc'


	#Parameters
	thickness=100.
	#Radius=0.25*1000
	#Radius=sqrt(3)/2.*1000
	#Radius=1.*1000
	Radius=0.5*10000.  #Hexagon only valid for S<half gridcell.  (about 0.85 using 2km)
	print 'Radius = ', Radius
	rho_ice=850.
	gravity=9.8
	mass_scaling=1.
	R_earth=6360.*1000.
	buffer_number=0   #Number of buffer points on the sides of the domain

	#Let interactive radius be different from radius for testing the model:
	IA_scaling=1.#(1./2.)
	Radius=Radius/IA_scaling

	#Here we create the lons and lats for a tabular iceberg
	#N= 5  # Number of rows in iceberg
	#M= 4   # Number of columns in iceberg

	if element_type=='square' and Interpolate_from_four_corners==True:
		print 'Square packing with R dividing dx, works best with interpolation off.'

	print 'Element type= ', element_type, ';  Switched_regridding= ', Switch_regridding_element_type
	if break_some_bonds==True:
		print 'Some bonds are being broken'

	#####################################################################################
	#####################################################################################

	if  Ice_geometry_source=='ISOMIP':
		if ISOMIP_reduced==False:
			lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
			(x,y,ice_mask,h_ice)=load_ISOMIP_ice_geometry(ISOMIP_ice_geometry_filename,buffer_number)
			ice_filename=ISOMIP_ice_geometry_filename
		else:
			lon_init=0  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
			(x,y,ice_mask,h_ice)=load_ISOMIP_reduced_ice_geometry(ISOMIP_reduced_ice_geometry_filename,buffer_number,ISOMIP_topography_filename)
			ice_filename=ISOMIP_reduced_ice_geometry_filename
		topog_filename=ISOMIP_topography_filename


	if  Ice_geometry_source=='Weddell':
		lon_init=-32.9  ; lat_init=-70.  #latitude  of bottom left corner of iceberg
		(x,y,ice_mask,h_ice)=load_Weddel_ice_geometry(Weddell_ice_geometry_filename)
		ice_filename=Weddell_ice_geometry_filename
	
	if input_is_cartesian is False:
		dx=10000.
		lat_grid=y ; lon_grid=x
		(x,y,ice_mask,h_ice)=convert_input_to_catesian_coordinates(x,y,ice_mask,h_ice,R_earth,dx)
		input_is_cartesian=True


	if set_all_thicknesses_to_one==True:
		print 'All thicknesses being set equal to 1'
		if set_all_domain_to_ice is True:
			ice_mask[:,:]=1.
		h_ice[np.where(ice_mask>0.5)]=Th_prescribed

	
	lat_ref=np.max(y)
	if (Convert_to_lat_lon==True)  and (input_is_cartesian==True):
		lat_ref=lat_init
		x=x-np.min(x)
		y=y-np.min(y)
		x=lon_init+((x/R_earth)*(180./np.pi)*(1/np.cos(lat_ref*np.pi/180.)))
		y=lat_init+((y/R_earth)*(180./np.pi))
		input_is_cartesian=False

	#Define the positions,thickness, mass,  of the icebergs
	(Number_of_bergs,lon,lat,iceberg_num,dx_berg, dy_berg,thickness,mass,width,x,y,Radius,static_berg,N_bergs_before_bd)= Create_icebergs(lon_init,lat_init,\
			Radius,R_earth, x, y,ice_mask,h_ice,Convert_to_lat_lon,rho_ice,input_is_cartesian,element_type,scale_the_grid_to_lat_lon,lat_ref,adjust_lat_ref,\
			Interpolate_from_four_corners,Fill_in_the_boundaries, set_all_bergs_static_by_default,break_some_bonds,Find_thickness_using_vector,Remove_stationary_bergs)

	print 'Maximum thickness:', np.max(thickness), np.min(thickness)
	#Define the positions of the iceberg bonds
	if Create_icebergs_bonds==True:
		(Number_of_bonds, first_berg_num,first_berg_ine,first_berg_jne,first_berg_lat,first_berg_lon, other_berg_num,other_berg_ine, other_berg_jne,other_berg_lat,other_berg_lon)=\
				Define_iceberg_bonds(Number_of_bergs,iceberg_num,lat,lon,dx_berg, dy_berg,Radius,break_some_bonds,static_berg,Allow_bonds_for_static_iceberg\
				,N_bergs_before_bd,Allow_bonds_with_boundary_bergs)

	temp_mass=[mass[i]/thickness[i] for i in range(len(thickness))]
	if only_choose_one_berg==True:
		(Number_of_bergs,lon,lat,thickness,width,mass,iceberg_num)= Select_just_one_berg(lon,lat,thickness,width,mass,iceberg_num,chosen_berg_num,static_berg)


	if regrid_icebergs_onto_grid==True:
		if scale_the_grid_to_lat_lon==True:
			print 'Regridding should be run with scale_the_grid_to_lat_lon off'
		else:
			if Switch_regridding_element_type==True:
				if element_type=='square':
					element_type='hexagon'
				else:
					element_type='square'
			temp_thickness=[1. for i in thickness]
			temp_mass=[mass[i]/(thickness[i]*rho_ice) for i in range(len(thickness))]
			temp_rho_ice=1.
			New_area=regrid_iceberg_thickness(lat,lon,Number_of_bergs,temp_thickness,temp_mass,h_ice,x,y,temp_rho_ice,element_type,static_berg,plot_outcome=False)
			New_mass=regrid_iceberg_thickness(lat,lon,Number_of_bergs,thickness,mass,h_ice,x,y,rho_ice,element_type,static_berg,plot_outcome=False)
		if save_new_h_ice_file==True:
			grid_area=(x[1]-x[0])*(x[1]-x[0])
			h_ice_new=New_mass*0
			M=h_ice.shape
			for j in range(M[0]):
				for i in range(M[1]):
					if New_area[j,i]>0:
						h_ice_new[j,i]=New_mass[j,i]/(rho_ice*New_area[j,i])  

	if Switch_x_and_y_to_rotate_90 is True:
		[lat,lon, h_ice_new,New_area,New_mass]=Switch_x_and_y_directions(lat,lon, h_ice_new,New_area,New_mass)
			

	if save_restart_files==True:	
		#Creating iceberg restart file
		Create_iceberg_restart_file(Number_of_bergs, lon,lat,thickness,width,mass,mass_scaling,iceberg_num,Ice_geometry_source,static_berg)
		
		#Creating bond restart file
		if Create_icebergs_bonds==True:
			Create_bond_restart_file(Number_of_bonds,first_berg_num,first_berg_ine,first_berg_jne,other_berg_ine,other_berg_jne,iceberg_num,other_berg_num,Ice_geometry_source)
			
		if save_new_h_ice_file==True:
			create_clipped_icethickness_file(h_ice_new,New_area,New_mass,grid_area,gravity)

	# Plotting the positions and bonds of the newly formed formation
	if plotting_bergs_over_grid==True:
		plotting_iceberg_positions(lat,lon,Number_of_bergs,R_earth,Radius,IA_scaling,Convert_to_lat_lon,\
				plot_circles,h_ice,ice_mask,x,y,plot_ice_mask,plot_ice_thickness,thickness,plot_icebergs_positions,static_berg)
		if (plot_bonds==True) and (Create_icebergs_bonds):
			plotting_iceberg_bonds(first_berg_lat,first_berg_lon,other_berg_lat,other_berg_lon,Number_of_bonds)


	#field=New_mass/(rho_ice*grid_area)
	#field=New_area/grid_area -1. 
	#M=field.shape	
	#print 'field',(field[M[0]-1,:])
	#print 'field',(field[M[0]-2,:])
	#print 'h_ice',h_ice_new[:,:]
	#print 'h_ice-1',h_ice_new[:,:]-1.
	plt.show()
	print 'Script complete'



if __name__ == '__main__':
	main()
	#sys.exit(main())














