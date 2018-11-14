# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:16:45 2018

@author: Maria
"""

import numpy as np
from numpy.linalg import inv



def nan_borehole(x,y,sc,a): 
##The function sets all the coordinates inside the borehole (r<1) to nan 
##Inputs: x- x coordinates
        # y- y coordinattes
        # sc- sizes of matrix
        # a- radius of the borehole (unit, m)
##Outputs: newx- x coodinates with the borehole coordinates=nan
         # newy- y coodinates with the borehole coordinates=nan
    for i in range(0,sc):
        for j in range(0,sc):
            r=(((x[i,j])**2)+((y[i,j])**2))/a
            if r<a:          
                x[i,j]=np.nan
                y[i,j]=np.nan
          
            elif x[i,j]==0: 
            #####when x=0 or y=0 is unstable so change to nan
                x[i,j]=np.nan
            elif y[i,j]==0:
                y[i,j]=np.nan               
    newx=x
    newy=y
    return newx,newy      


    
def beta(A): ##The matrix of reduced strain coefficients, equation (12)
##The function calculates the  reduced strain coefficents
##Inputs: A- the compliance matrix     
##Outputs: Beta =matrix of reduced strain coefficients   
    beta=np.zeros((6,6))
    for i in range(0,6):
        for j in range(0,6):
            beta[i,j]=(A[i,j]-((A[i,2]*A[j,2])/A[2,2]))
    return beta

def miu(Beta):## Eqn. (C.2)
##The function calculates the roots of equation (C.2)
##Inputs: Beta =matrix of reduced strain coefficients     
##Outputs: mius - the roots of equation (C.2)      
    weili2=np.poly1d([Beta[4,4],-2*Beta[3,4],Beta[3,3]])
    weili3=[Beta[0,4], -(Beta[0,3]+Beta[4,5]), Beta[1,4]+Beta[3,5], -Beta[1,3]]
    weili4=np.poly1d([Beta[0,0],-2*Beta[0,5],((2*Beta[0,1])+Beta[5,5]),-2*Beta[1,5],Beta[1,1]])
    ##For general case
    sixorder=np.polymul(weili4,weili2)-np.polymul(weili3, weili3) # Eqn.(C.1)
    rootsmiu=np.roots(sixorder)
    mius=rootsmiu[np.where((rootsmiu).imag>=0)]
    return mius ## solutions of the sixorder equation (C.2)

def lsoln (Beta,mius): 
##The function solves the functions in equation (C.2) and considers the special cases
##Inputs: Beta =matrix of reduced strain coefficients     
        # mius - the roots of equation (C.2)   
##Outputs: l2,l3,l4- functions in (C.2)
        #  mius - the roots of equation (C.2) recalculated for the special cases
    weili2=np.poly1d([Beta[4,4],-2*Beta[3,4],Beta[3,3]])
    weili3=[Beta[0,4], -(Beta[0,3]+Beta[4,5]), Beta[1,4]+Beta[3,5], -Beta[1,3]]
    weili4=np.poly1d([Beta[0,0],-2*Beta[0,5],((2*Beta[0,1])+Beta[5,5]),-2*Beta[1,5],Beta[1,1]])
    ## Eqn. (12)
    l4=np.polyval(weili4,mius)
    l3=np.polyval(weili3,mius)
    l2=np.polyval(weili2,mius)
    
    if ( np.all(np.abs(l3)< 10**(-6))): ##for special cases
        rootsmiu1=np.roots(weili4)
        roots1=rootsmiu1[np.where((rootsmiu1).imag>=0)]
        rootsmiu2=np.roots(weili2)
        roots2=rootsmiu2[np.where((rootsmiu2).imag>=0)]
        mius=np.zeros((3,1),dtype=complex)
        mius[0]=roots1[0]
        mius[1]=roots1[1]
        mius[2]=roots2[0]
        
    return l2,l3,l4,mius

def zeta (z,a,miu,laplace,sc):#Eqn. (D.6)
##The function solves the equation (D.2) and chooses the correct sign for plotting
    ##the analytical solution
##Inputs: z -complex variable mapping function
        # a- radius of the borehole (unit, m)
        # miu - the root of equation (C.2) 
        #laplace - variable from Eqn. (20)
        # sc- sizes of matrix
##Outputs:zeta- variable in equation (D.2) with the correct sign    
    eta01=(((z/a)**2)-1-(miu**2))**(0.5) ## numeritor of equation (D.6)
    eta02=-((((z/a)**2)-1-(miu**2))**(0.5))
    eta1=( (z/a)+ eta01)/(1-miu*1j)## equation (D.6)
    eta2=( (z/a)+ eta02)/(1-miu*1j)
    zeta1=2*laplace*eta1*eta01## Eqn.(19) which is the dominator of Eqn.(18)
    zeta2=2*laplace*eta2*eta02
    zeta=zeta1  
    ## Take the sign of the zeta with the largest absolute value
    ## This is referred to section 3.2 "Ambiguities in sign" in the manuscript.
    zeta[np.where(np.abs(zeta1)<(np.abs(zeta2)))]=zeta2[np.where(np.abs(zeta1)<(np.abs(zeta2)))]
    return zeta


def borehole_wall(theta2,geophys,R,a):
##The function grabs the data 180 degrees around the borehole wall 
##Inputs: theta2- ## Angle around borehole in degrees
        # geophys- the stress parameter being plotted around borehole
        # x2- x coordinates of the geophys matrix
        # y2- y coordinates of the geophys matrix
##Outputs: wall2- 180 degree theta2 coodinates around the borehole wall
         # geophys2- values of the geophys parameter 180 degrees around the borehole wall
    wall=theta2[np.where(R/a<=np.nanmin(R/a)+0.0086)]
    geophys=geophys[np.where(R/a<=np.nanmin(R/a)+0.0086)]
    wall2=np.array(zip(*sorted(zip(wall[np.where(wall>0)],geophys[np.where(wall>0)])))[0])
    geophys2=np.array(zip(*sorted(zip(wall[np.where(wall>0)],geophys[np.where(wall>0)])))[1])
    return wall2,geophys2

def borehole_wall360(theta2,geophys,R,a):
##The function grabs the data 360 degrees around the borehole wall 
##Inputs: theta2- ## Angle around borehole in degrees
        # geophys- the stress parameter being plotted around borehole
        # x2- x coordinates of the geophys matrix
        # y2- y coordinates of the geophys matrix
##Outputs: wall2- 360 degree theta2 coodinates around the borehole wall
         # geophys2- values of the geophys parameter 360 degrees around the borehole wall
    wall=theta2[np.where(R/a<=np.nanmin(R/a)+0.005)]
    geophys=geophys[np.where(R/a<=np.nanmin(R/a)+0.005)]
    wall2=np.array(zip(*sorted(zip(wall,geophys)))[0])
    geophys2=np.array(zip(*sorted(zip(wall,geophys)))[1])
    
    return wall2,geophys2

