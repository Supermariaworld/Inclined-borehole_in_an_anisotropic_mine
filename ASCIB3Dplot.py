# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 17:17:38 2018

@author: Maria
"""
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from scipy.interpolate import griddata

def plottin(geophys,theta,R,normrad,title,fgnm):
## The function plots the stress variable, geophys, 360 around the borehole in 
    ## in polar coordinates a radial distance, normrad, away from the borehole wall
## Inputs: geophys- the stress parameter being plotted around borehole
        # theta-angle around borehole in radians 
        # R-radial distance away from borehole 
        # normrad- distance into far field measured by normalized borehole radius 
        #title- title of the plot
        #fgnm- figure number
## Outputs: plot of one stress variable, geophys, in polar coordinates 
    mgeophys = ma.masked_invalid(geophys)
    mtheta = ma.masked_invalid(theta)
    mR = ma.masked_invalid(R)
    fig=plt.figure(fgnm)
    ax =fig.add_subplot(111,polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rticks(np.arange(1, normrad+1, 1))
    plt.pcolormesh(mtheta,mR,mgeophys, cmap='gist_ncar')
    plt.title(title)
    plt.axis('scaled')
    plt.axis([0,360,0,normrad])
    ax.grid(True)
    plt.xlabel('r/a')
    cbar=plt.colorbar()
    cbar.set_label('MPa',rotation=270)
    plt.show()
    
def plottin2(geophys,x2,y2,normrad,title,fgnm):
## The function plots the stress variable, geophys, 360 around the borehole in 
## in cartesian coordinates a radial distance, normrad, away from the borehole wall
## Inputs: geophys- the stress parameter being plotted around borehole
        # x2- x coordinates of the geophys matrix
        # y2- y coordinates of the geophys matrix
        # normrad- distance into far field measured by normalized borehole radius 
        #title- title of the plot
        #fgnm- figure number
## Outputs: plot of one stress variable, geophys, in cartesian coordinates 
    mgeophys = ma.masked_invalid(geophys)
    plt.figure(fgnm)
    plt.pcolormesh(x2,y2,mgeophys,cmap='gist_ncar')
    plt.title(title)
    plt.axis('scaled')
    plt.axis([-normrad,normrad,-normrad,normrad])
    plt.ylabel('r/a')
    plt.xlabel('r/a')
    cbar=plt.colorbar()
    cbar.set_label('MPa',rotation=270)
    plt.show()
    
def plottin3(radial,hoop,axial,taotz,taorz,taort,normrad,title):
## The function plots all of the stress variables, 360 around the borehole wall
## Inputs: radial,hoop,axial,taotz,taorz,taort - all directions of stress
        # normrad- distance into far field measured by normalized borehole radius 
        #title- title of the plot
## Outputs: plot of the stress at the borehole wall from -180 to 180 degrees
    plt.figure(9)
    plt.clf
    plt.plot(normrad,radial,'k', label='Radial')
    plt.plot(normrad,hoop,'blue',label='Hoop')
    plt.plot(normrad,axial,'r', label='Axial')
    plt.plot(normrad,taotz,'g',label='Tau_theta_z')
    plt.plot(normrad,taorz,'g--',label='Tau_r_z')
    plt.plot(normrad,taort,'k--',label='Tau_r_theta')
    plt.title(title)
    plt.axis([-180,180,-20,round(np.nanmax(hoop)+10,-1)])
    plt.xticks(np.arange(-180, 180+1, 30))
    plt.yticks(np.arange(-30, round(np.nanmax(hoop)+10,-1)+1, 20))
    plt.xlabel('Theta')
    plt.ylabel('Stress in MPa')
#    plt.legend(loc='upper right', shadow=True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid('on')
    plt.show()
    
def plottin4(radial,hoop,axial,taotz,taorz,taort,normrad,title, nradial,nhoop,naxial,ntaotz,ntaorz,ntaort):
    radial=np.round(radial,-1)
    hoop=np.round(hoop,2)
    axial=np.round(axial,2)
    taotz=np.round(taotz,1)
    taorz=np.round(taorz,1)
    taort=np.round(taort,1)
    
    nradial=np.round(nradial,-1)
    nhoop=np.round(nhoop,2)
    naxial=np.round(naxial,2)
    ntaotz=np.round(ntaotz,1)
    ntaorz=np.round(ntaorz,1)
    ntaort=np.round(ntaort,1) 
    
    pradial=abs((abs(abs(radial)-abs(nradial))/(abs(abs(radial)+abs(nradial))/2.0)) )*100
    phoop=abs((abs(abs(hoop)-abs(nhoop))/(abs(abs(hoop)+abs(nhoop))/2.0)) )*100
    paxial=abs((abs(abs(axial)-abs(naxial))/(abs(abs(axial)+abs(naxial))/2.0)) )*100
    
    ptaotz=abs((abs(abs(taotz)-abs(ntaotz))/(abs(abs(taotz)+abs(ntaotz))/2.0)) )*100
    ptaorz=abs((abs(abs(taorz)-abs(ntaorz))/(abs(abs(taorz)+abs(ntaorz))/2.0)) )*100
    ptaort=abs((abs(abs(taort)-abs(ntaort))/(abs(abs(taort)+abs(ntaort))/2.0)) )*100
    
    plt.figure(9)
    plt.clf
    plt.plot(normrad,pradial,'k', label='Radial')
    plt.plot(normrad,phoop,'blue',label='Hoop')
    plt.plot(normrad,paxial,'r', label='Axial')
#    plt.plot(normrad,ptaotz,'g',label='Tao_theta_z')
#    plt.plot(normrad,ptaorz,'g--',label='Tao_r_z')
#    plt.plot(normrad,ptaort,'k--',label='Tao_r_theta')
    plt.title(title)
    plt.xlim([-180,180])
    plt.xticks(np.arange(-180, 180+1, 30))
#    plt.yticks(np.arange(-20, round(np.nanmax(hoop)+10,-1)+1, 10))
    plt.xlabel('Theta')
    plt.ylabel('Percent Difference')
#    plt.legend(loc='upper right', shadow=True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid('on')
    plt.show()
    
def plottin4b(radial,hoop,axial,taotz,taorz,taort,normrad,title, nradial,nhoop,naxial,ntaotz,ntaorz,ntaort):
    
    radial=np.round(radial)
    hoop=np.round(hoop,2)
    axial=np.round(axial,2)
    taotz=np.round(taotz,4)
    taorz=np.round(taorz)
    taort=np.round(taort)
    
    nradial=np.round(nradial)
    nhoop=np.round(nhoop,2)
    naxial=np.round(naxial,2)
    ntaotz=np.round(ntaotz,4)
    ntaorz=np.round(ntaorz)
    ntaort=np.round(ntaort) 
    
    pradial=abs((abs(abs(radial)-abs(nradial))/(abs(abs(radial)+abs(nradial))/2.0)) )*100
    phoop=abs((abs(abs(hoop)-abs(nhoop))/(abs(abs(hoop)+abs(nhoop))/2.0)) )*100
    paxial=abs((abs(abs(axial)-abs(naxial))/(abs(abs(axial)+abs(naxial))/2.0)) )*100
    
    ptaotz=abs((abs(abs(taotz)-abs(ntaotz))/(abs(abs(taotz)+abs(ntaotz))/2.0)) )*100
    ptaorz=abs((abs(abs(taorz)-abs(ntaorz))/(abs(abs(taorz)+abs(ntaorz))/2.0)) )*100
    ptaort=abs((abs(abs(taort)-abs(ntaort))/(abs(abs(taort)+abs(ntaort))/2.0)) )*100
    
    plt.figure()
    plt.clf
#    plt.plot(normrad,pradial,'k', label='Radial')
#    plt.plot(normrad,phoop,'blue',label='Hoop')
#    plt.plot(normrad,paxial,'r', label='Axial')
    plt.plot(normrad,ptaotz,'g',label='Tau_theta_z')
    plt.plot(normrad,ptaorz,'g--',label='Tau_r_z')
    plt.plot(normrad,ptaort,'k--',label='Tau_r_theta')
    plt.title(title)
    plt.xlim([-180,180])
    plt.xticks(np.arange(-180, 180+1, 30))
#    plt.yticks(np.arange(-20, round(np.nanmax(hoop)+10,-1)+1, 10))
    plt.xlabel('Theta')
    plt.ylabel('Percent Difference')
#    plt.legend(loc='upper right', shadow=True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid('on')
    plt.show()    
  
def plottin4c(radial,hoop,axial,taotz,taorz,taort,normrad,title, nradial,nhoop,naxial,ntaotz,ntaorz,ntaort):
    radial=np.round(radial,2)
    hoop=np.round(hoop,2)
    axial=np.round(axial,2)
    taotz=np.round(taotz,1)
    taorz=np.round(taorz)
    taort=np.round(taort)
    
    nradial=np.round(nradial,2)
    nhoop=np.round(nhoop,2)
    naxial=np.round(naxial,2)
    ntaotz=np.round(ntaotz,1)
    ntaorz=np.round(ntaorz)
    ntaort=np.round(ntaort) 
    
    
    pradial=abs(radial)-abs(nradial)
    phoop=abs(hoop)-abs(nhoop)
    paxial=abs(axial)-abs(naxial)
    
    ptaotz=abs(taotz)-abs(ntaotz)
    ptaorz=abs(taorz)-abs(ntaorz)
    ptaort=abs(taort)-abs(ntaort)
    
    plt.figure()
    plt.clf
    plt.plot(normrad,pradial,'k', label='Radial')
    plt.plot(normrad,phoop,'blue',label='Hoop')
    plt.plot(normrad,paxial,'r', label='Axial')
#    plt.plot(normrad,ptaotz,'g',label='Tao_theta_z')
#    plt.plot(normrad,ptaorz,'g--',label='Tao_r_z')
#    plt.plot(normrad,ptaort,'k--',label='Tao_r_theta')
    plt.title(title)
    plt.xlim([-180,180])
    plt.xticks(np.arange(-180, 180+1, 30))
#    plt.yticks(np.arange(-20, round(np.nanmax(hoop)+10,-1)+1, 10))
    plt.xlabel('Theta')
    plt.ylabel('Stress Magnitude Difference (MPa)')
#    plt.legend(loc='upper right', shadow=True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid('on')
    plt.show() 
    
def plottin4d(radial,hoop,axial,taotz,taorz,taort,normrad,title, nradial,nhoop,naxial,ntaotz,ntaorz,ntaort):
    radial=np.round(radial,2)
    hoop=np.round(hoop,2)
    axial=np.round(axial,2)
    taotz=np.round(taotz,4)
    taorz=np.round(taorz,-1)
    taort=np.round(taort,-1)
    
    nradial=np.round(nradial,2)
    nhoop=np.round(nhoop,2)
    naxial=np.round(naxial,2)
    ntaotz=np.round(ntaotz,4)
    ntaorz=np.round(ntaorz,-1)
    ntaort=np.round(ntaort,-1) 
    
    pradial=abs(radial)-abs(nradial)
    phoop=abs(hoop)-abs(nhoop)
    paxial=abs(axial)-abs(naxial)
    
    ptaotz=abs(taotz)-abs(ntaotz)
    ptaorz=abs(taorz)-abs(ntaorz)
    ptaort=abs(taort)-abs(ntaort)
    
    plt.figure()
    plt.clf
#    plt.plot(normrad,pradial,'k', label='Radial')
#    plt.plot(normrad,phoop,'blue',label='Hoop')
#    plt.plot(normrad,paxial,'r', label='Axial')
    plt.plot(normrad,ptaotz,'g',label='Tau_theta_z')
    plt.plot(normrad,ptaorz,'g--',label='Tau_r_z')
    plt.plot(normrad,ptaort,'k--',label='Tau_r_theta')
    plt.title(title)
    plt.xlim([-180,180])
    plt.xticks(np.arange(-180, 180+1, 30))
#    plt.yticks(np.arange(-20, round(np.nanmax(hoop)+10,-1)+1, 10))
    plt.xlabel('Theta')
    plt.ylabel('Stress Magnitude Difference (MPa)')
#    plt.legend(loc='upper right', shadow=True)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.grid('on')
    plt.show()    
      
    
    
def plottin5(geophys,theta,normrad,sc):    
## The function creates a 3D the stress variable, geophys, 360 around the borehole 
## as a function of both theta and radius away from borehole
## Inputs: geophys- the stress parameter being plotted around borehole
        # theta-angle around borehole in radians 
        # normrad- distance into far field measured by normalized borehole radius 
        # sc- sizes of matrix
## Outputs: 3D plot of stress variable as a function of both theta and radius
    y=normrad
    x=theta*180/np.pi
    v=geophys
    y_flat=y.ravel()
    x_flat=x.ravel()
    v_flat=v.ravel()
    points=[]
    z_out=[]
    for i in range(len(y_flat)):
        if not np.isnan(y_flat[i]):
            points.append([x_flat[i],y_flat[i]])
            z_out.append(v_flat[i])
        
    xi=np.linspace(-180,180,sc)
    yi=np.linspace(0,11,sc)
    X,Y=np.meshgrid(xi,yi)
    zi=griddata(points,z_out,(xi[None,:],yi[:,None]),method='cubic')
    fig=plt.figure(5)
    ax = p3.Axes3D(fig)
    ax.plot_surface(X,Y,zi)
#    ax.scatter(x[::10,::10],y[::10,::10],z[::10,::10])#rstride=15,cstride=15)     
#    ax.set_xlim ([1,359])   
#    ax.set_ylim([1,10]) 
    ax.set_xlabel('Theta')
    ax.set_ylabel('R/a')
    ax.set_zlabel('Stress, MPa')
    plt.show()