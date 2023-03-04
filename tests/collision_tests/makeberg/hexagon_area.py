#!/usr/bin/env python

#First import the netcdf4 library
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import numpy as np  # http://code.google.com/p/netcdf4-python/
import matplotlib
import math
import os
matplotlib.use("tkagg")
from pylab import *
#import matplotlib.pyplot as plt
import pdb
import netCDF4 as nc



#######################

def Divide_hexagon_into_4_quadrants_old(x0,y0,H):
	S=(2/sqrt(3))*H
	Area_hex=3*(np.sqrt(3)/2)*(S**2)  #Area of the hexagon (should be equal to Area/grid_area, since it is not dim)	 - check this
	#print 'y0,H',y0, H
	#Defining boundaries of hexagon, and if statements to see which side of the boundary you are on
	W1= False ; W2=False ;W3=False ;W4=False ;W5=False ;W6=False ;H1=False ;V1=False ;V2=False
	W1=(-y0<-sqrt(3)*(-x0) + ((sqrt(3)*(S))));#upper right
	W2=(-y0<(H));#Top
	W3=(-y0<sqrt(3)*(-x0) + ((sqrt(3)*(S))));#Upper left
	W4=(-y0<sqrt(3)*(-x0) + (-(sqrt(3)*(S))));#Lower right
	W5=(-y0>-H);#Bottom
	W6=(-y0<-sqrt(3)*(-x0) + (-(sqrt(3)*(S))));#Lower left


	T1=(-x0<S) #Right
	T2=(-y0<(H));#Top
	T3= (-x0>-S) #Left
	T4=(-y0>-H);#Bottom

	H1=(y0<0);
	V1=(x0<-(S/2));
	V2=(x0<(S/2));

	#Deciding if the origin is within the hexagon
	#print W1 , W2 , W3 , W4 , W5 , W6
	#if In_hex:
	#	print In_hex

	#Calculating the area of the top and bottom half of the hexagon, 2 Cases for the majority above and below the y0=0 line
	#(and two more for the hexagon totally above and below the y0=0 line)
	if abs(y0)<H:
		Trapesium=((sqrt(3)*H)-(abs(y0)/sqrt(3)))*(H-abs(y0));
		if y0>=0:
			Area_Lower=Trapesium;
			Area_Upper=Area_hex-Trapesium;
		else:
			Area_Upper=Trapesium;
			Area_Lower=Area_hex-Trapesium;
	else:
		if y0>=0:
			Area_Lower=0.;
			Area_Upper=Area_hex;
		else:
			Area_Lower=Area_hex;
			Area_Upper=0.;

	#Calcularing Left and Right area of the hexagon, about the x0=0 line, 3 cases:
	#(and two more for when the hexagon is totally to the left or right of the x0=0 line)
	if abs(x0)<S:
		if abs(x0)<S/2:
			Rectangle=(abs(x0)*2*H);
			Big_side =(Area_hex/2) +Rectangle;
			Small_side=Area_hex-Big_side;
		else:
			Triangle=(sqrt(3))*((S-abs(x0))**2);
			Small_side=Triangle;
			Big_side=Area_hex-Small_side;
		if x0>=0.:
			Area_right=Big_side;
			Area_left=Small_side;
		else:
			Area_right=Small_side;
			Area_left=Big_side;
	else:
		if x0>=0.:
			Area_right=Area_hex;
			Area_left=0.;
		else:
			Area_right=0.;
			Area_left=Area_hex;

	In_hex= W1 & W2 & W3 & W4 & W5 & W6;
	In_hex_box=T1 & T2 & T3 & T4
	Area_Q1=0.; Area_Q2=0. ; Area_Q3=0.; Area_Q4=0.;
	Sector=0
	#if In_hex==False: #Then the hexagon is completely contained in the middle cell
	if In_hex_box==False: #Then the hexagon is completely contained in the middle cell
		Sector=-1
		#mass_on_ocean[i,j,5]=mass_on_ocean[i,j,5]+Mass
		if min(Area_Upper,Area_Lower)==0.:
			Sector=-2
			if Area_Upper==0.:
				Area_Q3=Area_left;
				Area_Q4=Area_right;
			if Area_Lower==0.:
				Area_Q1=Area_right;
				Area_Q2=Area_left;

		elif min(Area_right,Area_left)==0.:
			Sector=-3
			if Area_right==0.:
				Area_Q2=Area_Upper;
				Area_Q3=Area_Lower;

			if Area_left==0.:
				Area_Q1=Area_Upper;
				Area_Q4=Area_Lower;


		#yCxC=1.
		#print 'out of hex'

	else:
		#Determine which sector within the hexagon you are in. (sectors 1 to 6 go counter clockwise starting with top right)
		if (H1==True):	#Bottom half
			if V1:
				#if W1==False:
				if ((y0+(sqrt(3)*(x0+S)))<=0.):
					Sector=1;
				else:
					Sector=2;
			elif (V1==False) & (V2==True):
				Sector=3;
			else:
				#if (W3==True):
				if ((y0-(sqrt(3)*(x0-S)))>=0.):
					Sector=4;
				else:
					Sector=5;
		else:  #Bottom half
			if V1:
				#if W6==False:
				if ((y0 -(sqrt(3)*(x0+S)))>=0.):
					Sector=10;
				else:
					Sector=9;
			elif (V1==False) & (V2==True):
				Sector=8;
			else:
				#if (W4==True):
				if ((y0+(sqrt(3)*(x0-S)))<=0.):
					Sector=7;
				else:
					Sector=6;

		#print Sector

		#If the hexagon is in Sector 1,3,4 or 6, then the intersetion of the hexagon and the corresponding sector forms a baby triangle
		#If the hexagon is in Sector 2,5 then the intersetion of the hexagon and the corresponding sector forms a baby trapesoid
		if Sector==2 or Sector==4  or Sector==7 or Sector==9:
			Baby_triangle=(1/(2*sqrt(3)))*((-abs(y0)+(sqrt(3)*(S-abs(x0))))**2);
		else:
			#Baby_trap= (H-abs(y0)) * ((-H-abs(y0)+(2*sqrt(3)*(S-abs(x0))))/(2*sqrt(3)));
			Baby_trap=(H-abs(y0))*((S-abs(x0) - ((H+abs(y0))/(2*sqrt(3))))) ;

		#Finally, we assign the correct areas in each quadrant (Q1,Q2,Q3,Q4), depending on which sector you are in.
		C1=0.;C2=0.;C3=0.;C4=0.;

		#Corner cases
		if Sector==2:
			Area_Q1=Baby_triangle;
			Area_Q2=Area_Upper-Area_Q1
			Area_Q3=Area_left-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==4:
			Area_Q2=Baby_triangle;
			Area_Q1=Area_Upper-Area_Q2
			Area_Q3=Area_left-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==7:
			Area_Q3=Baby_triangle;
			Area_Q2=Area_left-Area_Q3
			Area_Q1=Area_Upper-Area_Q2
			Area_Q4=Area_right-Area_Q1

		if Sector==9:
			Area_Q4=Baby_triangle;
			Area_Q1=Area_right-Area_Q4
			Area_Q2=Area_Upper-Area_Q1
			Area_Q3=Area_left-Area_Q2

		#Center cases
		if Sector==3:
			if x0<=0.:
				Area_Q1=Baby_trap;
				Area_Q2=Area_Upper-Area_Q1;
				Area_Q3=Area_left-Area_Q2;
				Area_Q4=Area_right-Area_Q1;
			else:
				Area_Q2=Baby_trap;
				Area_Q1=Area_Upper-Area_Q2;
				Area_Q3=Area_left-Area_Q2;
				Area_Q4=Area_right-Area_Q1;

		if Sector==8:
			if x0<=0.:
				Area_Q4=Baby_trap;
				Area_Q3=Area_Lower-Area_Q4;
				Area_Q1=Area_right-Area_Q4;
				Area_Q2=Area_Upper-Area_Q1;
			else:
				Area_Q3=Baby_trap;
				Area_Q4=Area_Lower-Area_Q3;
				Area_Q1=Area_right-Area_Q4;
				Area_Q2=Area_Upper-Area_Q1;

		#Outside triangle cases:
		if Sector==1:
			Area_Q1=0.;
			Area_Q2=Area_Upper;
			Area_Q4=Area_right;
			Area_Q3=Area_left-Area_Q2;

		if Sector==5:
			Area_Q2=0.;
			Area_Q1=Area_Upper;
			Area_Q3=Area_left;
			Area_Q4=Area_Lower-Area_Q3;

		if Sector==6:
			Area_Q3=0.;
			Area_Q2=Area_left;
			Area_Q4=Area_Lower;
			Area_Q1=Area_right-Area_Q4;

		if Sector==10:
			Area_Q4=0.;
			Area_Q3=Area_Lower;
			Area_Q1=Area_right;
			Area_Q2=Area_left-Area_Q3;

		#print x0,y0,Sector
	return [Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4]



def Hexagon_into_quadrants_using_triangles_across_lines(x0,y0,H,theta,	Px, Py, Qx, Qy, Rx, Ry, Sx, Sy ):
	#In this routine, the hexagon is divided across lines PQ, RS.
	#Length of side of Hexagon
	S=(2/sqrt(3))*H;

	#Finding positions of corners
	C1x=S		; C1y=0.  #Corner 1 (right)
	C2x=H/sqrt(3.)	; C2y=H;  #Corner 2 (top right)
	C3x=-H/sqrt(3.) ; C3y=H;  #Corner 3 (top left)
	C4x=-S		; C4y=0.; #Corner 4 (left)
	C5x=-H/sqrt(3.) ; C5y=-H; #Corner 5 (bottom right)
	C6x=H/sqrt(3.)	; C6y=-H; #Corner 6 (bottom right)

	#Finding positions of corners
	[C1x,C1y]=rotate_and_translate(C1x,C1y,theta,x0,y0)
	[C2x,C2y]=rotate_and_translate(C2x,C2y,theta,x0,y0)
	[C3x,C3y]=rotate_and_translate(C3x,C3y,theta,x0,y0)
	[C4x,C4y]=rotate_and_translate(C4x,C4y,theta,x0,y0)
	[C5x,C5y]=rotate_and_translate(C5x,C5y,theta,x0,y0)
	[C6x,C6y]=rotate_and_translate(C6x,C6y,theta,x0,y0)

	#Area of Hexagon is the sum of the triangles
	[T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4]=Triangle_divided_into_four_parts(x0,y0,C1x,C1y,C2x,C2y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T012
	[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_parts(x0,y0,C2x,C2y,C3x,C3y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	[T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4]=Triangle_divided_into_four_parts(x0,y0,C3x,C3y,C4x,C4y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T034
	[T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4]=Triangle_divided_into_four_parts(x0,y0,C4x,C4y,C5x,C5y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	[T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4]=Triangle_divided_into_four_parts(x0,y0,C5x,C5y,C6x,C6y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	[T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4]=Triangle_divided_into_four_parts(x0,y0,C6x,C6y,C1x,C1y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023

	#Area of Hexagon is the sum of the triangles
	#print '1:', Triangle_divided_into_four_parts(x0,y0,C1x,C1y,C2x,C2y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T012
	#print '1b:', Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y); #T012
	#print '2', Triangle_divided_into_four_parts(x0,y0,C2x,C2y,C3x,C3y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	#print '2b',Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y); #T023
	#print '3:',Triangle_divided_into_four_parts(x0,y0,C3x,C3y,C4x,C4y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T034
	#print '3b',Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y); #T034
	#print '4',Triangle_divided_into_four_parts(x0,y0,C4x,C4y,C5x,C5y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	#print '4b:',Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y); #T023
	#print '5', Triangle_divided_into_four_parts(x0,y0,C5x,C5y,C6x,C6y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	#print '5b',Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y); #T023
	#print '6:',Triangle_divided_into_four_parts(x0,y0,C6x,C6y,C1x,C1y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
	#print '6b',Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y); #T023

	#Summing up
	Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area;
	Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1;
	Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2;
	Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3;
	Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4;
	#Area_Q4=Area_hex-(Area_Q1+Area_Q2+Area_Q3)




	#if  (abs(x0-(-0.186672393129))<0.00001)  and (abs(y0-(-0.436171))<0.00001):
	#	print 'Triangle 1,xo,y0',x0,y0
	#	print 'Triangle 1H,S',H,S
	#	print 'Triangle 1,Corners',C1x,C1y,C2x,C2y
	#	print 'Triangle 1, Area full', T12_Area
	#	print 'Triangle 1, Area T1', T12_Q1
	#	print 'Triangle 1, Area T2', T12_Q2
	#	print 'Triangle 1, Area T3', T12_Q3
	#	print 'Triangle 1, Area T4', T12_Q4
	#	#print 'Triangle 2',x0,y0,C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4
	#	#print 'Triangle 3',x0,y0,C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4
	#	#print 'Triangle 4',x0,y0,C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4
	#	#print 'Triangle 5',x0,y0,C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4
	#	#print 'Triangle 6',x0,y0,C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4

	Area_Q1=max(Area_Q1,0.);
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.);


	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if (abs(Error)>0.01):
		print  ('diamonds, hex error, H,x0,y0, Error', H, x0 , y0, Error)
		print ('diamonds, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,	 Area_Q4)
		print ('Hexagon error is large!!')



	#Adjust Areas so that the error is zero by subtracting the error from the largest sector.
	if  (((Area_Q1>=Area_Q2) and (Area_Q1>=Area_Q3)) and (Area_Q1>=Area_Q4)):
		#print 'fix1',Error
		Area_Q1=Area_Q1+Error
	elif  (((Area_Q2>=Area_Q1) and (Area_Q2>=Area_Q3)) and (Area_Q2>=Area_Q4)):
		#print 'fix2', Error
		Area_Q2=Area_Q2+Error
	elif  (((Area_Q3>=Area_Q1) and (Area_Q3>=Area_Q2)) and (Area_Q3>=Area_Q4)):
		#print 'fix3',Error
		Area_Q3=Area_Q3+Error
	elif  (((Area_Q4>=Area_Q1) and (Area_Q4>=Area_Q2)) and (Area_Q4>=Area_Q3)):
		#print 'fix4',Error
		Area_Q4=Area_Q4+Error
	else:
		print ('There is some thing wrong with this hexagon. Something very wrong')

	#Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if ((abs(Error)>0.01)):
		print ('The hexagon error is still too large!!!', Error)



	return [Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4]

def Hexagon_into_quadrants_using_triangles(x0,y0,H,theta):

	#Length of side of Hexagon
	S=(2/sqrt(3))*H;

	#Finding positions of corners
	C1x=S		; C1y=0.  #Corner 1 (right)
	C2x=H/sqrt(3.)	; C2y=H;  #Corner 2 (top right)
	C3x=-H/sqrt(3.) ; C3y=H;  #Corner 3 (top left)
	C4x=-S		; C4y=0.; #Corner 4 (left)
	C5x=-H/sqrt(3.) ; C5y=-H; #Corner 5 (bottom right)
	C6x=H/sqrt(3.)	; C6y=-H; #Corner 6 (bottom right)

	#Finding positions of corners
	[C1x,C1y]=rotate_and_translate(C1x,C1y,theta,x0,y0)
	[C2x,C2y]=rotate_and_translate(C2x,C2y,theta,x0,y0)
	[C3x,C3y]=rotate_and_translate(C3x,C3y,theta,x0,y0)
	[C4x,C4y]=rotate_and_translate(C4x,C4y,theta,x0,y0)
	[C5x,C5y]=rotate_and_translate(C5x,C5y,theta,x0,y0)
	[C6x,C6y]=rotate_and_translate(C6x,C6y,theta,x0,y0)

	#Area of Hexagon is the sum of the triangles
	[T12_Area,T12_Q1,T12_Q2,T12_Q3,T12_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C1x,C1y,C2x,C2y); #T012
	[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C2x,C2y,C3x,C3y); #T023
	[T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C3x,C3y,C4x,C4y); #T034
	[T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C4x,C4y,C5x,C5y); #T023
	[T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y); #T023
	[T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4]=Triangle_divided_into_four_quadrants(x0,y0,C6x,C6y,C1x,C1y); #T023

	#Summing up
	Area_hex=T12_Area+T23_Area+T34_Area+T45_Area+T56_Area+T61_Area;
	Area_Q1=T12_Q1+T23_Q1+T34_Q1+T45_Q1+T56_Q1+T61_Q1;
	Area_Q2=T12_Q2+T23_Q2+T34_Q2+T45_Q2+T56_Q2+T61_Q2;
	Area_Q3=T12_Q3+T23_Q3+T34_Q3+T45_Q3+T56_Q3+T61_Q3;
	Area_Q4=T12_Q4+T23_Q4+T34_Q4+T45_Q4+T56_Q4+T61_Q4;
	#Area_Q4=Area_hex-(Area_Q1+Area_Q2+Area_Q3)




	#if  (abs(x0-(-0.186672393129))<0.00001)  and (abs(y0-(-0.436171))<0.00001):
	#	print 'Triangle 1,xo,y0',x0,y0
	#	print 'Triangle 1H,S',H,S
	#	print 'Triangle 1,Corners',C1x,C1y,C2x,C2y
	#	print 'Triangle 1, Area full', T12_Area
	#	print 'Triangle 1, Area T1', T12_Q1
	#	print 'Triangle 1, Area T2', T12_Q2
	#	print 'Triangle 1, Area T3', T12_Q3
	#	print 'Triangle 1, Area T4', T12_Q4
	#	#print 'Triangle 2',x0,y0,C2x,C2y,C3x,C3y,T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4
	#	#print 'Triangle 3',x0,y0,C3x,C3y,C4x,C4y,T34_Area,T34_Q1,T34_Q2,T34_Q3,T34_Q4
	#	#print 'Triangle 4',x0,y0,C4x,C4y,C5x,C5y,T45_Area,T45_Q1,T45_Q2,T45_Q3,T45_Q4
	#	#print 'Triangle 5',x0,y0,C5x,C5y,C6x,C6y,T56_Area,T56_Q1,T56_Q2,T56_Q3,T56_Q4
	#	#print 'Triangle 6',x0,y0,C6x,C6y,C1x,C1y,T61_Area,T61_Q1,T61_Q2,T61_Q3,T61_Q4

	Area_Q1=max(Area_Q1,0.);
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.);


	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if (abs(Error)>0.01):
		print ( 'diamonds, hex error, H,x0,y0, Error', H, x0 , y0, Error)
		print ('diamonds, hex error, Areas',Area_hex, (Area_Q1+Area_Q2 + Area_Q3+Area_Q4), Area_Q1,  Area_Q2 , Area_Q3,	 Area_Q4)
		print ('Hexagon error is large!!')



	#Adjust Areas so that the error is zero by subtracting the error from the largest sector.
	if  (((Area_Q1>=Area_Q2) and (Area_Q1>=Area_Q3)) and (Area_Q1>=Area_Q4)):
		#print 'fix1',Error
		Area_Q1=Area_Q1+Error
	elif  (((Area_Q2>=Area_Q1) and (Area_Q2>=Area_Q3)) and (Area_Q2>=Area_Q4)):
		#print 'fix2', Error
		Area_Q2=Area_Q2+Error
	elif  (((Area_Q3>=Area_Q1) and (Area_Q3>=Area_Q2)) and (Area_Q3>=Area_Q4)):
		#print 'fix3',Error
		Area_Q3=Area_Q3+Error
	elif  (((Area_Q4>=Area_Q1) and (Area_Q4>=Area_Q2)) and (Area_Q4>=Area_Q3)):
		#print 'fix4',Error
		Area_Q4=Area_Q4+Error
	else:
		print ('There is some thing wrong with this hexagon. Something very wrong')

	#Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	Error=Area_hex-(Area_Q1+Area_Q2+Area_Q3+Area_Q4)
	if ((abs(Error)>0.01)):
		print ('The hexagon error is still too large!!!', Error)



	return [Area_hex ,Area_Q1, Area_Q2, Area_Q3, Area_Q4]


def find_intersect_using_metric(Ax,Ay, Bx, By, Cx, Cy):
	#Line AB
	[AB_px, AB_py]=intercept_of_a_line(Ax,Ay,Bx,By,'x'); #x_intercept
	[AB_qx, AB_qy]=intercept_of_a_line(Ax,Ay,Bx,By,'y'); #y_intercept
	M_AB=point_in_interval_metric(Ax,Ay,Bx,By,AB_px,AB_py, AB_qx, AB_qy)
	#Line AC
	[AC_px, AC_py]=intercept_of_a_line(Ax,Ay,Cx,Cy,'x'); #x_intercept
	[AC_qx, AC_qy]=intercept_of_a_line(Ax,Ay,Cx,Cy,'y'); #y_intercept
	M_AC=point_in_interval_metric(Ax,Ay,Cx,Cy,AC_px, AC_py,AC_qx, AC_qy)
	#Line AB
	[BC_px, BC_py]=intercept_of_a_line(Bx,By,Cx,Cy,'x'); #x_intercept
	[BC_qx, BC_qy]=intercept_of_a_line(Bx,By,Cx,Cy,'y'); #y_intercept
	M_BC=point_in_interval_metric(Bx,By,Cx,Cy,BC_px,BC_py,BC_qx, BC_qy)

	if M_AB<=min(M_BC, M_AC):
		return AB_px, AB_py, AB_qx, AB_qy
	elif M_AC<=min(M_BC, M_AB):
		return AC_px, AC_py, AC_qx, AC_qy
	elif M_BC<=min(M_AC, M_AB):
		return BC_px, BC_py, BC_qx, BC_qy
	else:
		print ('You should not get here')
		halt

def find_intersect_using_metric_and_lines(Ax,Ay, Bx, By, Cx, Cy,Px, Py, Qx, Qy, Rx, Ry, Sx, Sy):
	#Line AB
	[AB_px, AB_py]=intercept_of_two_lines(Ax,Ay,Bx,By,Px, Py, Qx, Qy); #PQ_intercept
	[AB_qx, AB_qy]=intercept_of_two_lines(Ax,Ay,Bx,By,Rx, Ry, Sx, Sy); #RS_intercept
	M_AB=point_in_interval_metric(Ax,Ay,Bx,By,AB_px,AB_py, AB_qx, AB_qy)
	#Line AC
	[AC_px, AC_py]=intercept_of_two_lines(Ax,Ay,Cx,Cy,Px, Py, Qx, Qy); #PQ_intercept
	[AC_qx, AC_qy]=intercept_of_two_lines(Ax,Ay,Cx,Cy, Rx, Ry, Sx, Sy); #RS_intercept
	M_AC=point_in_interval_metric(Ax,Ay,Cx,Cy,AC_px, AC_py,AC_qx, AC_qy)
	#Line AB
	[BC_px, BC_py]=intercept_of_two_lines(Bx,By,Cx,Cy,Px, Py, Qx, Qy); #PQ_intercept
	[BC_qx, BC_qy]=intercept_of_two_lines(Bx,By,Cx,Cy, Rx, Ry, Sx, Sy); #RS_intercept
	M_BC=point_in_interval_metric(Bx,By,Cx,Cy,BC_px,BC_py,BC_qx, BC_qy)

	if M_AB<=min(M_BC, M_AC):
		return AB_px, AB_py, AB_qx, AB_qy
	elif M_AC<=min(M_BC, M_AB):
		return AC_px, AC_py, AC_qx, AC_qy
	elif M_BC<=min(M_AC, M_AB):
		return BC_px, BC_py, BC_qx, BC_qy
	else:
		print ('You should not get here')
		halt


def point_in_interval_metric(Ax,Ay,Bx,By,px,py,qx,qy):
	metric=0.0
	#print 'Ax,Ay,Bx,By',Ax,Ay,Bx,By,px,py,qx,qy
	#print 'px,py,qx,qy',px,py,qx,qy

	#Finds a metric for how close a point is to being inside an interval. If it is inside, then metric is equal to zero.
	Mx1= max(px - max(Ax,Bx),0.)  #Zero is px <= max(Ax,Bx)
	Mx2= max( min(Ax,Bx)-px ,0.)  #Zero is px <= max(Ax,Bx)
	My1= max(py - max(Ay,By),0.)  #Zero is px <= max(Ax,Bx)
	My2= max( min(Ay,By)-py ,0.)  #Zero is px <= max(Ax,Bx)
	p_metric=abs(Mx1)+abs(Mx2)+abs(My1)+abs(My2)
	#print 'Mx1, Mx2, My1, My2', Mx1, Mx2, My1, My2

	#Finds a metric for how close a point is to being inside an interval. If it is inside, then metric is equal to zero.
	Mx1= max(qx - max(Ax,Bx),0.)  #Zero is qx <= max(Ax,Bx)
	Mx2= max( min(Ax,Bx)-qx ,0.)  #Zero is qx <= max(Ax,Bx)
	My1= max(qy - max(Ay,By),0.)  #Zero is qx <= max(Ax,Bx)
	My2= max( min(Ay,By)-qy ,0.)  #Zero is qx <= max(Ax,Bx)
	q_metric=abs(Mx1)+abs(Mx2)+abs(My1)+abs(My2)

	metric=p_metric+q_metric
	return metric


def Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy):
	Area_key_quadrant=0.0 #Initializing
	Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	#if Area_triangle==0.0:
	#	Area_triangle=0.0 ; Area_Q1=0.0 ; Area_Q2=0.0; Area_Q3=0.0; Area_Q4= 0.0
	#	return [Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4]


	#Calculating area across axes
	[Area_Upper ,Area_Lower]=divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'x');
	[Area_Right ,Area_Left]=divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,'y');

	#Decide if the origin is in the triangle
	if point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.):	#Then you have to divide area 4 ways.
		#Find a line in the triangle that cuts both axes in/on the trianlge
		[px, py]=intercept_of_a_line(Ax,Ay,Bx,By,'x'); #x_intercept
		[qx, qy]=intercept_of_a_line(Ax,Ay,Bx,By,'y'); #y_intercept
		if (point_in_interval(Ax,Ay,Bx,By,px,py) & point_in_interval(Ax,Ay,Bx,By,qx,qy))==False:
			[px, py]=intercept_of_a_line(Ax,Ay,Cx,Cy,'x'); #x_intercept
			[qx, qy]=intercept_of_a_line(Ax,Ay,Cx,Cy,'y'); #y_intercept
			if (point_in_interval(Ax,Ay,Cx,Cy,px,py) & point_in_interval(Ax,Ay,Cx,Cy,qx,qy))==False:
				[px, py]=intercept_of_a_line(Bx,By,Cx,Cy,'x'); #x_intercept
				[qx, qy]=intercept_of_a_line(Bx,By,Cx,Cy,'y'); #y_intercept
				if (point_in_interval(Bx,By,Cx,Cy,px,py) & point_in_interval(Bx,By,Cx,Cy,qx,qy))==False:
					#point_in_interval_metric
					[px, py, qx, qy]= find_intersect_using_metric(Ax,Ay, Bx, By, Cx, Cy)
					#print 'Houston, we have a problem'
					#plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy)
					#halt

		#Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
		Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.);
		if Area_key_quadrant>0.0:
			if px>=0. and qy>=0.: #First quadrant
				Key_quadrant=1;
			elif px<0. and qy>=0.:	#Second quadrant
				Key_quadrant=2;
			elif px<0. and qy<0.:
				Key_quadrant=3;	 #Third quadrant
			else:
				Key_quadrant=4;	 #Forth quadrant

	#if (point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.) is False) or (Area_key_quadrant==0):
	if Area_key_quadrant==0:
	#else:	#Then at least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
		#print 'Mother...'
		#print 'Ax, Ay',Ax,Ay
		#print 'Bx, By',Bx,By
		#print 'Cx, Cy',Cx,Cy
		Area_key_quadrant=0;
		if (((Ax>0. and Ay>0.) or (Bx>0. and By>0.) or (Cx>0. and Cy>0.))==False)  and (Area_Upper+Area_Right<=Area_triangle):
		#if (((Ax>=0. and Ay>=0.) or (Bx>=0. and By>=0.) or (Cx>=0. and Cy>=0.))==False)  and (Area_Upper+Area_Right<=Area_triangle):
			#No points land in this quadrant and triangle does not cross the quadrant
			Key_quadrant=1;
		elif  (((Ax<0. and Ay>0) or (Bx<0. and By>0.) or (Cx<0. and Cy>0.))==False)  and (Area_Upper+Area_Left<=Area_triangle):
			#elif  (((Ax<0. and Ay>=0) or (Bx<0. and By>=0.) or (Cx<0. and Cy>=0.))==False)	 and (Area_Upper+Area_Left<=Area_triangle):
			Key_quadrant=2;
		elif (((Ax<0. and Ay<0.) or (Bx<0. and By<0.) or (Cx<0. and Cy<0.))==False) & (Area_Lower+Area_Left<=Area_triangle):
			#elif (((Ax<0. and Ay<0.) or (Bx<0. and By<0.) or (Cx<0. and Cy<0.))==False) & (Area_Lower+Area_Left<=Area_triangle):
			Key_quadrant=3;
		else:
			Key_quadrant=4;

	#Assign values to quadrants
	if Key_quadrant==1:
		Area_Q1=Area_key_quadrant;
		Area_Q2=Area_Upper-Area_Q1;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==2:
		Area_Q2=Area_key_quadrant;
		Area_Q1=Area_Upper-Area_Q2;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==3:
		Area_Q3=Area_key_quadrant;
		Area_Q2=Area_Left-Area_Q3;
		Area_Q1=Area_Upper-Area_Q2;
		#Area_Q4=Area_Right-Area_Q1;
		Area_Q4=Area_triangle-Area_Q1-Area_Q2-Area_Q3;
	elif Key_quadrant==4:
		Area_Q4=Area_key_quadrant;
		Area_Q1=Area_Right-Area_Q4;
		Area_Q2=Area_Upper-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	else:
	    print( 'Help, I need somebody, help!')
	    halt

	Area_Q1=max(Area_Q1,0.);
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.);

	Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
	#Marker 1
	if Error>1e-15:
	#if True:
		print ('The triangles are not accurate enough. This is a problem!')
		print ('Triangle corners: ' ,Ax,Ay,Bx,By,Cx,Cy)
		print ('Error',Error)
		print ('Key_quadrant', Key_quadrant )
		print ('Area Key_quadrant', Area_key_quadrant)
		print ('Upper, Lower',Area_Upper ,Area_Lower)
		print ('Left, Right ', Area_Right ,Area_Left)
		print ('Point in triangle', point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
		#plt.plot(np.array([Ax,Bx,Cx,Ax]),np.array([Ay, By, Cy,Ay]))
		#plt.plot(np.array([-0.5,0.5]),np.array([0.,0.]),'k')
		#plt.plot(np.array([0.,0.]),np.array([-0.5,0.5]),'k')
		#plt.plot(0,0,'*')
		#plt.show()
		#halt
		#return

	return [Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4]

def find_quadrant_using_lines(x0 , y0 , Px, Py, Qx, Qy, Rx, Ry, Sx, Sy): #Finding "quadrant" of point using lines PQ, RS as axes.  (This assumes PQ is like the x axis and RS like the y axis)
	#above_PQ=Point_above_line(x0,y0, Px, Py, Qx, Qy)
	#above_RS=Point_above_line(x0,y0, Rx, Ry, Sx, Sy)
	above_PQ=Above_on_below_line(x0,y0, Px, Py, Qx, Qy)
	above_RS=Above_on_below_line(x0,y0, Rx, Ry, Sx, Sy)
	if (above_PQ==1) and (above_RS==1):
		quadrant=1
	elif (above_PQ==1) and (above_RS==-1):
		quadrant=2
	elif (above_PQ==-1) and (above_RS==-1):
		quadrant=3
	elif (above_PQ==-1) and (above_RS==1):
		quadrant=4
	else:
		#print 'Point not in a quadrant!'
		quadrant=0
		#halt

	return quadrant

def Triangle_divided_into_four_parts(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy	):#Divides the triangle into 4 parts, across two intersecting lines. PQ, RS intersect at origin (for now)
	Area_key_quadrant=0.0
	Area_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	#if Area_triangle==0.0:
	#	Area_triangle=0.0 ; Area_Q1=0.0 ; Area_Q2=0.0; Area_Q3=0.0; Area_Q4= 0.0
	#	return [Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4]


	#Calculating area across axes
	[Area_Upper ,Area_Lower]=divding_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy);
	[Area_Right ,Area_Left] =divding_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Rx, Ry, Sx, Sy);

	#Decide if the origin is in the triangle
	if point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.):	#Then you have to divide area 4 ways.
		#Find a line in the triangle that cuts both axes in/on the trianlge
		[px, py]=intercept_of_two_lines(Ax,Ay,Bx,By,Px, Py, Qx, Qy); #PQ_intercept
		[qx, qy]=intercept_of_two_lines(Ax,Ay,Bx,By,Rx, Ry, Sx, Sy); #RS_intercept
		if (point_in_interval(Ax,Ay,Bx,By,px,py) & point_in_interval(Ax,Ay,Bx,By,qx,qy))==False:
			[px, py]=intercept_of_two_lines(Ax,Ay,Cx,Cy,Px, Py, Qx, Qy); #PQ_intercept
			[qx, qy]=intercept_of_two_lines(Ax,Ay,Cx,Cy,Rx, Ry, Sx, Sy); #RS_intercept
			if (point_in_interval(Ax,Ay,Cx,Cy,px,py) & point_in_interval(Ax,Ay,Cx,Cy,qx,qy))==False:
				[px, py]=intercept_of_two_lines(Bx,By,Cx,Cy,Px, Py, Qx, Qy); #PQ_intercept
				[qx, qy]=intercept_of_two_lines(Bx,By,Cx,Cy,Rx, Ry, Sx, Sy); #RS_intercept
				if (point_in_interval(Bx,By,Cx,Cy,px,py) & point_in_interval(Bx,By,Cx,Cy,qx,qy))==False:
					#print 'Houston, we have a problem'
					[px, py, qx, qy]= find_intersect_using_metric_and_lines(Ax,Ay, Bx, By, Cx, Cy,Px, Py, Qx, Qy, Rx, Ry, Sx, Sy)
					#plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy)
					#halt

		#Assigning quadrants. Key_quadrant is the quadrant with the baby triangle in it.
		Area_key_quadrant=Area_of_triangle(px,py,qx,qy,0.,0.);
		if px>=0. and qy>=0.: #First quadrant  (I think the zeros here are because the PQ, RS intersect at zero)
			Key_quadrant=1;
		elif px<0. and qy>=0.:	#Second quadrant
			Key_quadrant=2;
		elif px<0. and qy<0.:
			Key_quadrant=3;	 #Third quadrant
		else:
			Key_quadrant=4;	 #Forth quadrant

	if Area_key_quadrant==0:
		#else:	#Then at least one quadrant is empty, and this can be used to find the areas in the other quadrant.  Assigning quadrants. Key_quadrant is the empty quadrant.
		#print 'Mother...'
		#print 'Ax, Ay',Ax,Ay
		#print 'Bx, By',Bx,By
		#print 'Cx, Cy',Cx,Cy
		Area_key_quadrant=0;
		#Finding which quadrant the triangle points fall in.
		A_quad = find_quadrant_using_lines(Ax,Ay,Px, Py, Qx, Qy, Rx, Ry, Sx, Sy)
		B_quad = find_quadrant_using_lines(Bx,By,Px, Py, Qx, Qy, Rx, Ry, Sx, Sy)
		C_quad = find_quadrant_using_lines(Cx,Cy,Px, Py, Qx, Qy, Rx, Ry, Sx, Sy)
		#print 'A_quad, B_quad, C_quad',A_quad, B_quad, C_quad

		if (((A_quad==1) or (B_quad==1) or (C_quad==1))==False)	 and (Area_Upper+Area_Right<=Area_triangle):
			#No points land in this quadrant and triangle does not cross the quadrant
			Key_quadrant=1;
		elif  (((A_quad==2) or (B_quad==2) or (C_quad==2))==False)  and (Area_Upper+Area_Left<=Area_triangle):
			Key_quadrant=2;
		elif (((A_quad==3) or (B_quad==3) or (C_quad==3))==False) & (Area_Lower+Area_Left<=Area_triangle):
			Key_quadrant=3;
		else:
			Key_quadrant=4;

		#print Key_quadrant
	#Assign values to quadrants
	if Key_quadrant==1:
		Area_Q1=Area_key_quadrant;
		Area_Q2=Area_Upper-Area_Q1;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==2:
		Area_Q2=Area_key_quadrant;
		Area_Q1=Area_Upper-Area_Q2;
		Area_Q4=Area_Right-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	elif Key_quadrant==3:
		Area_Q3=Area_key_quadrant;
		Area_Q2=Area_Left-Area_Q3;
		Area_Q1=Area_Upper-Area_Q2;
		#Area_Q4=Area_Right-Area_Q1;
		Area_Q4=Area_triangle-Area_Q1-Area_Q2-Area_Q3;
	elif Key_quadrant==4:
		Area_Q4=Area_key_quadrant;
		Area_Q1=Area_Right-Area_Q4;
		Area_Q2=Area_Upper-Area_Q1;
		#Area_Q3=Area_Left-Area_Q2;
		Area_Q3=Area_triangle-Area_Q1-Area_Q2-Area_Q4;
	else:
	    print ('Help, I need somebody, help!')
	    halt

	Area_Q1=max(Area_Q1,0.);
	Area_Q2=max(Area_Q2,0.);
	Area_Q3=max(Area_Q3,0.);
	Area_Q4=max(Area_Q4,0.);

	Error=abs(Area_Q1+Area_Q2+Area_Q3+Area_Q4-Area_triangle)
	if Error>1e-15:
	#if True:
		print ('The triangles are not accurate enough. This is a problem!')
		print ('Triangle corners: ' ,Ax,Ay,Bx,By,Cx,Cy)
		print ('Error',Error)
		print ('Key_quadrant', Key_quadrant )
		print ('Area Key_quadrant', Area_key_quadrant )
		print ('Upper, Lower',Area_Upper ,Area_Lower)
		print ('Left, Right ', Area_Right ,Area_Left)
		print ('Point in triangle', point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))
		#return

	return [Area_triangle, Area_Q1, Area_Q2 ,Area_Q3 ,Area_Q4]

def divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1):
	if axes1=='x':	#Use the y-coordinates for if statements to see which side of the line you are on
		A0=Ay; B0=By;  C0=Cy;
	if axes1=='y':	#Use the y-coordinates for if statements to see which side of the line you are on
		A0=Ax; B0=Bx;  C0=Cx;

	#print 'A0,B0,C0', A0,B0,C0
	A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	if B0*C0>0.: #then B and C are on the same side	 (and non-zero)
		#print 'HHHH1'
		if A0*B0>=0.: #then all three on the the same side (if it equals zero, then A0=0 and the otehrs are not)
			if (A0>0.)  or	(A0==0. and  B0>0.):
				Area_positive= A_triangle;
				Area_negative= 0.;
			else:
				Area_positive= 0.;
				Area_negative= A_triangle;
		else:  #A is on the opposite side to B and C
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1);

	elif B0*C0<0.: #then B and C are on the opposite sides
		#print 'HHHH2'
		if A0*B0>=0.:  #then C is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Cx,Cy,Bx,By,Ax,Ay,axes1);
		else: #then B is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Bx,By,Cx,Cy,Ax,Ay,axes1);

	else:  #This is the case when either B or C is equal to zero (or both), A0 could be zero too.
		#print 'HHHH3'
		if (A0==0. and B0==0. and C0==0.):
			Area_positive= 0.;
			Area_negative= 0.;
		elif (A0*B0<0.)	 or  (A0*C0<0.):    #A, B are on opposite sides, and C is zero.	 OR  A, C are on opposite sides, and B is zero.
			[Area_positive, Area_negative]=Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1);
			#print 'HHHH6'
		elif (A0*B0>0.) or (A0*C0>0.) or (abs(A0)>0. and (B0==0.) and (C0==0.)):
			#print 'HHHH7'
			if (A0>0.):
				Area_positive= A_triangle;
				Area_negative= 0.;
			else:
				Area_positive= 0.;
				Area_negative= A_triangle;

		elif A0==0.:  #(Also, one of B,C is zero too)
			#print 'HHHH8'
			if B0>0. or C0>0.:
				Area_positive= A_triangle;
				Area_negative= 0.;
			elif B0<0. or C0<0.:
				Area_positive= 0.;
				Area_negative= A_triangle;
			else:
				print ('You should not get here1')
				halt
		else:
			print ('You should not get here2')
			print (Ax,Ay,Bx,By,Cx,Cy)
			halt

	return [Area_positive, Area_negative]

def divding_triangle_across_line2(Ax,Ay,Bx,By,Cx,Cy,  Px, Py, Qx, Qy ):
	A0= Above_on_below_line(Ax,Ay, Px, Py, Qx, Qy)
	B0= Above_on_below_line(Bx,By, Px, Py, Qx, Qy)
	C0= Above_on_below_line(Cx,Cy, Px, Py, Qx, Qy)

	#print 'A0,B0,C0', A0,B0,C0
	#print Px, Py, Qx, Qy

	A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	if B0*C0>0.: #then B and C are on the same side	 (and non-zero)
		#print 'GGG1'
		if A0*B0>=0.: #then all three on the the same side (if it equals zero, then A0=0 and the otehrs are not)
			if (A0>0.)  or	(A0==0. and  B0>0.):
				Area_positive= A_triangle;
				Area_negative= 0.;
			else:
				Area_positive= 0.;
				Area_negative= A_triangle;
		else:  #A is on the opposite side to B and C
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy );

	elif B0*C0<0.: #then B and C are on the opposite sides
		#print 'GGG2'
		if A0*B0>=0.:  #then C is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Cx,Cy,Bx,By,Ax,Ay, Px, Py, Qx, Qy );
		else: #then B is all alone
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Bx,By,Cx,Cy,Ax,Ay, Px, Py, Qx, Qy );

	else:  #This is the case when either B or C is equal to zero (or both), A0 could be zero too.
		#print 'GGG3'
		if (A0==0. and B0==0. and C0==0.):
			Area_positive= 0.;
			Area_negative= 0.;
		elif (A0*B0<0.)	 or  (A0*C0<0.):    #A, B are on opposite sides, and C is zero.	 OR  A, C are on opposite sides, and B is zero.
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy);
			#print 'GGG5', Area_positive, Area_negative
		elif (A0*B0>0.) or (A0*C0>0.) or (abs(A0)>0. and (B0==0.) and (C0==0.)):
			if (A0>0.):
				Area_positive= A_triangle;
				Area_negative= 0.;
			else:
				Area_positive= 0.;
				Area_negative= A_triangle;

		elif A0==0.:  #(Also, one of B,C is zero too)
			if B0>0. or C0>0.:
				Area_positive= A_triangle;
				Area_negative= 0.;
			elif B0<0. or C0<0.:
				Area_positive= 0.;
				Area_negative= A_triangle;
			else:
				print ('You should not get here1')
				halt
		else:
			print ('You should not get here2')
			print (Ax,Ay,Bx,By,Cx,Cy)
			halt

	return [Area_positive, Area_negative]


def divding_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy,   Px, Py, Qx, Qy ):	 #Triangle is ABC, line is PQ
	A0= Above_on_below_line(Ax,Ay, Px, Py, Qx, Qy)
	B0= Above_on_below_line(Bx,By, Px, Py, Qx, Qy)
	C0= Above_on_below_line(Cx,Cy, Px, Py, Qx, Qy)
	#print 'BBB0', A0, B0, C0

	A_triangle=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	# B and C are on the same side of the line
	#if (Point_above_line(Bx,By, Px, Py, Qx, Qy)==Point_above_line(Cx,Cy, Px, Py, Qx, Qy)) :
	#if ((B0==C0) and (B0!=0))  or ((B0==0 or C0==0) and np.max(B0,C0)>0.5)	  :
	if ((B0==C0) and (B0!=0))  or ((B0==0 or C0==0))   : #B and C on the same side or B or C is zero
		#print 'BBB1'
		[Area_positive, Area_negative]=Area_of_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy);
	else: #then B and C are on the opposite sides
		#if Point_above_line(Ax,Ay, Px, Py, Qx, Qy)==Point_above_line(Bx,By, Px, Py, Qx, Qy):  #then C is all alone
		if A0*B0>=0.:  #then C is all alone
			#print 'BBB2'
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Cx,Cy,Bx,By,Ax,Ay, Px, Py, Qx, Qy);
		else: #then B is all alone
			#print 'BBB3'
			[Area_positive, Area_negative]=Area_of_triangle_across_line(Bx,By,Cx,Cy,Ax,Ay, Px, Py, Qx, Qy);

	return [Area_positive, Area_negative]


def check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy,repeat_calculation=False):

	tol=0.00000000000000;
	dxc = qx - Ax;
	dyc = qy - Ay;

	dxl = Bx - Ax;
	dyl = By - Ay;

	cross = dxc * dyl - dyc * dxl;

	if abs(cross)<=tol:
		point_is_on_line=True
	else:
		point_is_on_line=False

	if repeat_calculation is False:
		if point_is_on_line != check_if_point_is_on_the_line(Bx,By,Ax,Ay,qx,qy, True):
			point_is_on_line=True

	return point_is_on_line

def intercept_of_a_line(Ax,Ay,Bx,By,axes1):
	No_intercept_val=100000000000.; #Huge value used to make sure that the intercept is outside the triange in the parralel case.
	#No_intercept_val=np.NaN;

	if axes1=='x': #x intercept
		if (Ay==By)==False:
			x0=Ax -(((Ax-Bx)/(Ay-By))*Ay);
			y0=0.;
		else:
			x0=No_intercept_val;
			y0=No_intercept_val;

	if axes1=='y': #y intercept
		if (Ax==Bx)==False:
			x0=0.;
			y0=-(((Ay-By)/(Ax-Bx))*Ax)+Ay;
		else:
			x0=No_intercept_val;
			y0=No_intercept_val;

	return [x0, y0]

def intercept_of_two_lines(Ax,Ay,Bx,By,Px,Py,Qx,Qy):
	No_intercept_val=100000000000.; #Huge value used to make sure that the intercept is outside the triange in the parralel case.
	#No_intercept_val=np.NaN;

	#One of the lines is actually a point:
	if ((Ax==Bx) and (Ay==By)) or ((Px==Qx) and (Py==Qy)):
		x0=No_intercept_val; y0=No_intercept_val;
		halt6
	else:
		if (Ax==Bx):
			(x0 , y0) =intercept_of_a_line(Px-Ax,Py,Qx-Ax,Qy,'y')
			x0=Ax
		#Here we check the other three directions to make code agree to machine precision with old version. The code could skip straight to the else.
		elif (Px==Qx):
			(x0 , y0) =intercept_of_a_line(Ax-Px,Ay,Bx-Px,By,'y')
			x0=Px
		elif (Ay==By):
			(x0 , y0) =intercept_of_a_line(Px,Py-Ay,Qx,Qy-Ay,'x')
			y0=Ay
		elif (Py==Qy):
			(x0 , y0) =intercept_of_a_line(Ax,Ay-Py,Bx,By-Py,'x')
			y0=Py
		else:
			m1=(Ay-By)/(Ax-Bx)
			if (Px==Qx):
				print ('Stop2')
				(x0 , y0) =intercept_of_a_line(Ax-Px,Ay,Bx-Px,By,'y')
				x0=Px
			else:
				m2=(Py-Qy)/(Px-Qx)
				if m1==m2:
					print( 'Stop3')
					x0=No_intercept_val;
					y0=No_intercept_val;
					halt2
				else:
					print( 'Stop4')
					x0=(1/(m1-m2))*( (m1*Ax -Ay) -(m2*Px - Py))
					y0=m1*(x0-Ax)+Ay

	return [x0, y0]


def Area_of_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axis1):#In this fuction, A is the point on one side of the axis, and B,C are on the opposite sides
	#print 'Area_of_triangle_across_axes' ,Ax,Ay,Bx,By,Cx,Cy
	A_triangle2=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);
	[pABx, pABy]=intercept_of_a_line(Ax,Ay,Bx,By,axis1);
	[pACx, pACy]=intercept_of_a_line(Ax,Ay,Cx,Cy,axis1);
	if axis1=='x':
		A0=Ay; #Value used for if statements (deciding up/down vs left/right)
	if axis1=='y':
		A0=Ax; #Value used for if statements (deciding up/down vs left/right)

	#print 'QQQ', Ax,Ay,pABx,pABy,pACx,pACy
	A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy);
	if (A0>=0.):
		Area_positive= A_half_triangle;
		Area_negative= A_triangle2-A_half_triangle;

	else:
		Area_positive= A_triangle2-A_half_triangle;
		Area_negative= A_half_triangle;

	return [Area_positive, Area_negative]

def Above_on_below_line(x,y,Px, Py, Qx, Qy): #If the point is above the line, then returns True, otherwise false
	if Px==Qx:
		if x>Px:
			Point_above_line=1.0
		elif x<Px:
			Point_above_line=-1.0
		else:
			Point_above_line=0.0
	else:
		m=(Py-Qy)/(Px-Qx)
		if y > m*(x-Px) + Py:
			Point_above_line=1.0
		elif y < m*(x-Px) + Py:
			Point_above_line=-1.0
		else:
			Point_above_line=0.0

	return Point_above_line

def Point_above_line(x,y,Px, Py, Qx, Qy): #If the point is above the line, then returns True, otherwise false
	if Px==Qx:
		if x>=Px:
			Point_above_line=True
		else:
			Point_above_line=False
	else:
		m=(Py-Qy)/(Px-Qx)
		if y >= m*(x-Px) + Py:
			Point_above_line=True
		else:
			Point_above_line=False

	return Point_above_line


def Area_of_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px,Py, Qx, Qy, debug=False):
	#In this fuction, B and C are the same side of the line. (I think this works even if B and/or C are on the line too)
	#A is unknown. PQ is the line.
	#print 'RRR Area_of_triangle_across_line' ,Ax,Ay,Bx,By,Cx,Cy

	if (Ax==Bx and Ay==By) or (Ax==Cx and Ay==Cy) or (Cx==Bx and Cy==By):
		#print 'AAA1'
		Area_positive=0.0
		Area_negative=0.0
	else:
		#print 'AAA2'
		A0= Above_on_below_line(Ax,Ay, Px, Py, Qx, Qy)
		B0= Above_on_below_line(Bx,By, Px, Py, Qx, Qy)
		C0= Above_on_below_line(Cx,Cy, Px, Py, Qx, Qy)
		A_triangle2=Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy);

		if check_if_point_is_on_the_line(Px,Py,Qx,Qy, Ax,Ay):  #Is point A on the line?
			#print 'AAA3'
			#if Point_above_line(Bx,By, Px,Py, Qx, Qy) or Point_above_line(Cx,Cy, Px,Py, Qx, Qy):
			if B0==1 or C0==1:
				Area_positive= A_triangle2;
				Area_negative= 0.0;
			else:
				Area_positive= 0.0;
				Area_negative= A_triangle2;
		else: #If A is not on the line
			#If A is on the same side of the line as B and C
			#All on the same side of the line
			if ((A0==B0 and A0==C0)	 or (A0==B0 and C0==0) or (A0==C0 and B0==0)):
			#if Point_above_line(Ax,Ay, Px,Py, Qx, Qy) == Point_above_line( Bx,By, Px,Py, Qx, Qy)  and (Point_above_line(Ax,Ay,Px,Py,Qx,Qy)==Point_above_line(Cx,Cy,Px,Py,Qx,Qy)):
				#print 'AAA4'
				#if Point_above_line(Bx,By, Px,Py, Qx, Qy) or Point_above_line(Cx,Cy, Px,Py, Qx, Qy):
				if A0==1 or C0==1:
					Area_positive= A_triangle2;
					Area_negative= 0.0;
				else:
					Area_positive= 0.0;
					Area_negative= A_triangle2;
			elif (C0==0 and B0==0 and (abs(A0)>0)):
				if Point_above_line(Ax,Ay, Px,Py, Qx, Qy):
					Area_positive= A_triangle2;
					Area_negative= 0.0;
				else:
					Area_positive= 0.0;
					Area_negative= A_triangle2;

			else:  #If A is not on the same side of the line as B and C
				#print 'AAA5'
				[pABx, pABy]=intercept_of_two_lines(Ax,Ay,Bx,By, Px, Py, Qx, Qy);
				[pACx, pACy]=intercept_of_two_lines(Ax,Ay,Cx,Cy, Px, Py, Qx, Qy);
				#print 'Tri2', Ax,Ay,pABx,pABy,pACx,pACy

				A_half_triangle=Area_of_triangle(Ax,Ay,pABx,pABy,pACx,pACy);
				if Point_above_line(Ax,Ay , Px,Py,Qx,Qy):
					Area_positive= A_half_triangle;
					Area_negative= A_triangle2-A_half_triangle;
				else:
					Area_positive= A_triangle2-A_half_triangle;
					Area_negative= A_half_triangle;

	return [Area_positive, Area_negative]

def Area_of_triangle(Ax,Ay,Bx,By,Cx,Cy):
	Area	=   abs(    0.5*((Ax*(By-Cy))+(Bx*(Cy-Ay))+(Cx*(Ay-By))) );
	return Area


def point_in_interval(Ax,Ay,Bx,By,px,py):
	point_is_in_interval=False
	#print 'NNN1'
	if ((px <= max(Ax,Bx)) and (px >= min(Ax,Bx))):
		#print 'NNN2', 'py=', py,  'max=',max(Ay,By), 'XXX', (py - max(Ay,By))
		if ((py <= max(Ay,By)) and (py >= min(Ay,By))):
			##print 'NNN3'
			point_is_in_interval=True
	return point_is_in_interval

def rotate_and_translate(px,py,theta,x0,y0):
	#This function takes a point px,py, and rotates it clockwise around the origin by theta degrees, and then translates by (x0,y0)
	 px_temp = (cos(theta*pi/180.)*px) + (sin(theta*pi/180.)*py)
	 py_temp = (-sin(theta*pi/180.)*px) + (cos(theta*pi/180.)*py)
	 px= px_temp +x0
	 py= py_temp +y0
	 return [px,py]


def point_in_triangle_old(Ax,Ay,Bx,By,Cx,Cy,qx,qy):
	#if False:
	if (Ax==qx and Ay==qy) or (Bx==qx and By==qy) or (Cx==qx and Cy==qy):  #Exclude the pathelogical case
			point_is_in_triangle = 0.;
	else:
		if (check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) or (check_if_point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy)) or (check_if_point_is_on_the_line(Bx,By,Cx,Cy,qx,qy))):
			point_is_in_triangle = 0;
		else:
			#Compute vectors
			v0x=Cx-Ax;
			v1x=Bx-Ax;
			v2x=qx-Ax;
			v0y=Cy-Ay;
			v1y=By-Ay;
			v2y=qy-Ay;

			#%Compute dot products
			dot00 = (v0x*v0x)+(v0y*v0y);
			dot01 = (v0x*v1x)+(v0y*v1y);
			dot02 = (v0x*v2x)+(v0y*v2y);
			dot11 = (v1x*v1x)+(v1y*v1y);
			dot12 = (v1x*v2x)+(v1y*v2y);

			#Compute barycentric coordinates
			invDenom= 1 / ((dot00 * dot11) - (dot01*dot01));
			u=((dot11*dot02)-(dot01*dot12))*invDenom;
			v=((dot00*dot12)-(dot01*dot02))*invDenom;

			point_is_in_triangle = (((u)>=0) & ((v)>=0) & ((u+v)<(1)));
	return point_is_in_triangle



def point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,qx,qy):
	point_is_in_triangle = 0;

	if (Ax==qx and Ay==qy) or (Bx==qx and By==qy) or (Cx==qx and Cy==qy):  #Exclude the pathelogical case
			point_is_in_triangle = 0.;
	else:
		if (check_if_point_is_on_the_line(Ax,Ay,Bx,By,qx,qy) or (check_if_point_is_on_the_line(Ax,Ay,Cx,Cy,qx,qy)) or (check_if_point_is_on_the_line(Bx,By,Cx,Cy,qx,qy))):
			point_is_in_triangle = 0;
		else:
			l0=(qx-Ax)*(By-Ay)-(qy-Ay)*(Bx-Ax)
			l1=(qx-Bx)*(Cy-By)-(qy-By)*(Cx-Bx)
			l2=(qx-Cx)*(Ay-Cy)-(qy-Cy)*(Ax-Cx)

			p0=np.sign( l0);
			p1=np.sign( l1);
			p2=np.sign( l2);
			if (l0 == 0.):
				p0=0.
			if (l1==0.):
				p1=0.
			if (l2==0.):
				p2=0.

			if ( (abs(p0)+abs(p2))+(abs(p1)) == abs((p0+p2)+(p1)) ):
				point_is_in_triangle = 1;
	return point_is_in_triangle


def roundoff(x,sig_fig):
	x_roundoff=round(x*(10**(sig_fig)))/(10**(sig_fig))
	#x_roundoff=(FLOAT (INT(x * (10.**sig_fig) + 0.5)) / (10.**sig_fig))
	return x_roundoff



def square_spreading_calculation(x,y,L):
	xL=min(0.5, max(0., 0.5-(x/L)))
	xR=min(0.5, max(0., (x/L)+(0.5-(1/L) )))
	xC=max(0., 1.-(xL+xR))
	yD=min(0.5, max(0., 0.5-(y/L)))
	yU=min(0.5, max(0., (y/L)+(0.5-(1/L) )))
	yC=max(0., 1.-(yD+yU))

def plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy, Px=None, Py=None, Qx=None, Qy=None, Rx=None, Ry=None, Sx=None, Sy=None):
	if Px is None:
		Px=-1. ; Py=0. ; Qx=1. ; Qy=0. #x-axis
		Rx=0. ; Ry=-1. ; Sx=0. ; Sy=1. #y-axis
	plt.plot(np.array([Px, Qx]), np.array([Py, Qy]))
	plt.plot(np.array([Rx, Sx]), np.array([Ry, Sy]))
	plt.plot(np.array([Ax, Bx, Cx, Ax]), np.array([Ay, By, Cy, Ay]))
	plt.plot(np.array([Ax, Bx, Cx, Ax]), np.array([Ay, By, Cy, Ay]),'8')
	text(Ax,Ay,'A',fontsize=20)
	text(Bx,By,'B',fontsize=20)
	text(Cx,Cy,'C',fontsize=20)
	plt.show()

####################################################################################################################################################
##########################################################  Main Program   #########################################################################
####################################################################################################################################################

def main():

	#test_case='hexagon'
	#test_case='triangle'
	#test_case='square'
	#test_case='point_in_triangle'
	#test_case='line_intercept'
	#test_case='Point_above_line'
	#test_case='Triangle_across_line'
	#test_case='Triangle_into_fours'
	test_case='many_hexagons'

	if test_case=='line_intercept':
		Ax=1. ; Ay=1.5 ; Bx=1. ; By=-3.4
		Px=0.2; Py=-1.0; Qx=0.2; Qy=2.2
		(x0,y0)=intercept_of_two_lines(Ax,Ay,Bx,By,Px,Py,Qx,Qy)
		plt.plot( np.array([Ax, Bx])  , np.array([Ay, By]))
		plt.plot(np.array([Px, Qx]), np.array([Py, Qy]))
		plt.plot(x0,y0,'*')
		#plt.plot(Ax,Ay,'*')
		#plt.plot(Bx,By,'*')
		#plt.plot(Px,Py,'*')
		#plt.plot(Qx,Qy,'*')

		print ('x0, y0' , x0 , y0)
		plt.show()


	if test_case=='Triangle_across_line':
		#y-axis
		Px=0.0 ; Py=-1.0 ;Qx=0.0 ;Qy=1.0
		axes1='y'

		#x-axis
		#Px=-1.0 ; Py=0.0 ;Qx=1.0 ;Qy=0.0
		#axes1='x'

		N=11 #Choose an odd number
		count=0
		for Ax in np.linspace(-0.5 , 0.5 , N):
			for Ay in np.linspace(-0.5 , 0.5 , N):
				for Bx in np.linspace(-0.5 , 0.5 , N):
					for By in np.linspace(-0.5 , 0.5 , N):
						for Cx in np.linspace(-0.5 , 0.5 , N):
							for Cy in np.linspace(-0.5 , 0.5 , N):
								count=count+1
								if count>0:
									x=divding_triangle_across_axes(Ax,Ay,Bx,By,Cx,Cy,axes1)
									y=divding_triangle_across_line2(Ax,Ay,Bx,By,Cx,Cy, Px,Py, Qx, Qy)
									#y=divding_triangle_across_line(Ax,Ay,Bx,By,Cx,Cy, Px,Py, Qx, Qy)
									#if abs(x[0]-y[0])>0 or	 abs(x[1]-y[1])>0 or  abs(x[0]-z[0])>0 or  abs(x[1]-z[1])>0:
									if abs(x[0]-y[0])>0 or	abs(x[1]-y[1])>0:
										print ('Broken!!!: ', 'Ax, Ay =', Ax, Ay, 'Bx, By =' , Bx, By, 'Cx, Cy =' ,Cx, Cy)
										print ('count', count)
										print (x)
										print( y)
										plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy)
										halt

		print ('Test Successful')



	if test_case=='Triangle_into_fours':
		Px=1. ; Py=0. ; Qx=2. ; Qy=0. #x-axis
		Rx=0. ; Ry=1. ; Sx=0. ; Sy=2. #y-axis

		count=0
		N=11 #Use an odd number
		for Ax in np.linspace(-0.5 , 0.5 , N):
			for Ay in np.linspace(-0.5 , 0.5 , N):
				for Bx in np.linspace(-0.5 , 0.5 , N):
					for By in np.linspace(-0.5 , 0.5 , N):
						for Cx in np.linspace(-0.5 , 0.5 , N):
							for Cy in np.linspace(-0.5 , 0.5 , N):
								count=count+1
								#print count
								if count>0:
									x=Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy)
									y=Triangle_divided_into_four_parts(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy )
									Error=abs(x[1]+x[2]+x[3]+x[4]-x[0])
									#if Error>1e-15:
									if abs(x[0]-y[0])>0 or	abs(x[1]-y[1])>0 or abs(x[2]-y[2])>0 or	 abs(x[3]-y[3])>0 or  abs(x[4]-y[4])>0 :
										print ('Broken!!!: ', 'Ax, Ay =', Ax, Ay, 'Bx, By =' , Bx, By, 'Cx, Cy =' ,Cx, Cy)
										print( 'count=', count)
										#print 'Error = ', Error
										print( x)
										#print y
										plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy)
										halt

		print ('Triangle_into_fours Test Successful')

		#[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy);
		#[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]= Triangle_divided_into_four_parts(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  )

	if test_case=='Point_above_line':
		x0=0.1 ;y0=0.3 ; Px=0.4; Py=-1.0; Qx=0.2; Qy=2.2
		plt.plot(np.array([Px, Qx]), np.array([Py, Qy]))
		plt.plot(x0,y0,'*')
		print (Point_above_line(x0,y0,Px, Py, Qx, Qy))
		plt.show()




	if test_case=='many_hexagons':
		Px=-1. ; Py=0. ; Qx=1. ; Qy=0. #x-axis
		Rx=0. ; Ry=-1. ; Sx=0. ; Sy=1. #y-axis
		N=101
		theta=0.0
		count=0
		for x0 in np.linspace(-0.5 , 0.5 , N):
			for y0 in np.linspace(-0.5 , 0.5 , N):
				for H in np.linspace(0.0 , 0.5 , N):
					count=count+1
					#print count
					if count>0:
						C5x=-H/sqrt(3.) ; C5y=-H; #Corner 5 (bottom right)
						C6x=H/sqrt(3.)	; C6y=-H; #Corner 6 (bottom right)
						[C5x,C5y]=rotate_and_translate(C5x,C5y,theta,x0,y0)
						[C6x,C6y]=rotate_and_translate(C6x,C6y,theta,x0,y0)
						#x = Triangle_divided_into_four_parts(x0,y0,C5x,C5y,C6x,C6y, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  ); #T023
						#y = Triangle_divided_into_four_quadrants(x0,y0,C5x,C5y,C6x,C6y); #T023
						#print x
						#print y
						#Ax=x0 ;  Ay=y0
						#Bx=C5x ; By=C5y ;
						#Cx=C6x ; Cy=C6y ;
						#plot_axes_and_triangel(Ax,Ay, Bx, By, Cx, Cy, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy)
						#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)=Hexagon_into_quadrants_using_triangles_across_lines(x0,y0,H,theta,  Px, Py, Qx, Qy, Rx, Ry, Sx, Sy )
						#(Area_hexb, Area_Q1b, Area_Q2b, Area_Q3b, Area_Q4b)= Hexagon_into_quadrants_using_triangles(x0,y0,H,theta=0)
						x=Hexagon_into_quadrants_using_triangles_across_lines(x0,y0,H,theta,  Px, Py, Qx, Qy, Rx, Ry, Sx, Sy )
						y= Hexagon_into_quadrants_using_triangles(x0,y0,H,theta=0)
						eps=1e-15
						#eps=0.0
						#if abs(Area_hex-Area_hexb)>eps or  abs(Area_Q1-Area_Q1b)>eps or abs(Area_Q2-Area_Q2b)>eps \
						#		or  abs(Area_Q3-Area_Q3b)>eps or  abs(Area_Q4-Area_Q4b)>eps :
						if abs(x[0]-y[0])>0 or	abs(x[1]-y[1])>0 or abs(x[2]-y[2])>0 or	 abs(x[3]-y[3])>0 or  abs(x[4]-y[4])>0 :
							print ('Broken!!!: ', x0, y0, H)
							print (Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4 )
							print (Area_hexb, Area_Q1b, Area_Q2b, Area_Q3b, Area_Q4b)
							print ('count=', count)
							halt

		print ('Many_Hexagons Test Successful')


	if test_case=='hexagon':
		H=1.
		S=2.*H/sqrt(3.)
		x0=0.
		y0=0.
		#(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Divide_hexagon_into_4_quadrants_old(x0,y0,H)
		#print x0,y0,H
		#print 'First version', Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4
		(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,theta=0)
		print ('Triangle version1',Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)
		(Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)= Hexagon_into_quadrants_using_triangles(x0,y0,H,theta=90.)
		print ('Triangle version2',Area_hex, Area_Q1, Area_Q2, Area_Q3, Area_Q4)

		#print 'Analysis'
		#print 'sum of quadrants: ', (Area_Q1+Area_Q2+Area_Q3+Area_Q4)
		print ('Scaled areas: ')
		print( (Area_Q1+Area_Q2+Area_Q3+Area_Q4)/Area_hex, Area_Q1/Area_hex, Area_Q2/Area_hex, Area_Q3/Area_hex, Area_Q4/Area_hex)


	elif test_case=='square':

		x0=0.6999999999999
		y0=0.3
		L=0.6
		square_spreading_calculation(x0,y0,L)

	elif test_case=='triangle':
		Ax=-0.5
		Ay=-0.3 -((8./9.)/10.)
		Bx=-0.277777777778
		By= -0.5

		Cx=0.5
		Cy=0.3 +((8./9.)/10.)

		[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]=Triangle_divided_into_four_quadrants(Ax,Ay,Bx,By,Cx,Cy); #T023

		#Px=1. ; Py=0. ; Qx=2. ; Qy=0.
		#Rx=0. ; Ry=1. ; Sx=0. ; Sy=2.
		#[T23_Area,T23_Q1,T23_Q2,T23_Q3,T23_Q4]= Triangle_divided_into_four_parts(Ax,Ay,Bx,By,Cx,Cy, Px, Py, Qx, Qy, Rx, Ry, Sx, Sy  )

		print ('Triangle 2: ',Ax,Ay,Bx,By,Cx,Cy)
		print ('Full',T23_Q1+T23_Q2+T23_Q3+T23_Q4)
		print ('Q1, Q2,Q3,Q4',T23_Q1,T23_Q2,T23_Q3,T23_Q4)
		print ('Full triangle area',T23_Area)
		print ('Error',(T23_Q1+T23_Q2+T23_Q3+T23_Q4-T23_Area))



	elif test_case=='point_in_triangle':
		Ax= -2.695732526092343E-012
		Ay=0.204344508198090
		Bx=-2.695750202346321E-012
		By= -8.433062639672301E-002
		Cx=0.249999999997304
		Cy=6.000694090068343E-002
		print ('Point in triangle', point_in_triangle(Ax,Ay,Bx,By,Cx,Cy,0.,0.))


	print ('Script complete')

if __name__ == '__main__':
	main()
	#sys.exit(main())
