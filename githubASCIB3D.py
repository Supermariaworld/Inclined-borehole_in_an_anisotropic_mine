# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 14:16:44 2018

@author: Maria Tibbo
"""
####### VALUES FOR MINE-BY TUNNEL
import numpy as np
from ASCIB3Dlibrary import*
from ASCIB3Dplot import*
from ASCIB3Drotation import*
#from noriteASCIB3D import nradialwall,nhoopwall,naxialwall,ntaotzwall,ntaorzwall,ntaortwall
############## Isotropic ###############
inputs=0
p0=2630
#Young=20
#Poisson=0.079
### Mine by
#Young=56.8
#Poisson=0.08781430017398313
Young=65
Poisson=0.25
zzC11=Young*(1-Poisson)/(1+Poisson)/(1- 2*Poisson)
zzC12=Young*Poisson/(1+Poisson)/(1-2*Poisson)
zzC44=(zzC11-zzC12)/2
c11=zzC11; c22=zzC11; c33=zzC11
c44=zzC44; c55=zzC44; c66=zzC44
c12=zzC12; c13=zzC12; c21=zzC12
c31=zzC12; c23=zzC12; c32=zzC12

###### BG Granite from Vermont ##### 
#### Due to microcracks

#inputs=0
#c11=132.7; c22=141.67; c33=156.17
#c44=20.59; c55=17.67; c66=15.78
#c12=7.41; c13=3.93; c23=6.43
#p0=2648
#Vp0=3540; Vs0=2835


######### Orthorombic ##########  Oliver Gaede, 2012
#inputs=0
#c11=35.6; c22=22.15; c33=16.97
#c44=6; c55=7; c66=8
#c12=8.2; c13=7.9; c23=3.85
#p0=2535
#Vp0=2746; Vs0=1764

######## VTI and TTI ######## Oliver Gaede, 2012
 ##Vertical axis of symmetry and tilted axis transverse isotropy
#inputs=0
#c11=45.2; c33=28; c44=7.05; c66=14.4; c13=19.67
#c12=c11-(2*c66); c22=c11; c23=c13; c55=c44
#p0=2535
#Vp0=3340; Vs0=1675

################################################
C1=np.matrix([[c11,c12,c13,0,0,0],
 [c12,c22,c23,0,0,0],
 [c13,c23,c33,0,0,0],
 [0,0,0,c44,0,0],
 [0,0,0,0,c55,0],
 [0,0,0,0,0,c66]])

### If using velocties making sure the inputs is set to 1
# inputs=1
if inputs==1:
    C1=np.matrix(C1vel)
    
######### Model Inputs ###########    
a=1.0 # radius of the borehole, most useful equal to 1 for normalization
taoxy0=0 ## Shear stresses (typically zero)
taoxz0=0
taoyz0=0
#NRS 17 sig1=54.9 sig 2= 41.2 sig 3=34.3
#sigmax=54.9 ##magnitude of maximum horizontal stress, MPa
#sigmay=41.2 ##magnitude of minimum horizontal stress, MPa
#sigmaz=34.3 ## magnitude of vertical stress, MPa
pw=0  ## magnitude of wellbore pressure, MPa

###Mine by
sigmax=60
sigmay=45
sigmaz=11
stress1=np.matrix([[sigmax,taoxy0,taoxz0], ## Far-field stress state
    [taoxy0,sigmay,taoyz0],
    [taoxz0,taoyz0,sigmaz]])
    
sc=1000 ## Size of matrix sc X sc
normrad=4## Distance into far field measured by normalized borehole radius 

#### Orientation of the borehole in degrees
#350 azi and 34 inclinaion for NRS 17
borehole1deg=54## azimuth alphazb
borehole2deg=90## inclination  alphayb
## Convert to Radians
borehole1=borehole1deg*np.pi/180 ## azimuth
borehole2=borehole2deg*np.pi/180 ## inclination
#NRS Collar Dip -56.0
	#Collar Az 350

############### Angle Option 1: material and far-field stress reference frames ###################

## Orientation of far-field stress frame defined referred to the NED (North(x)-East(y)-Downward(z)) coordinate system
orientS=[37,3,0]# orientS=[strikeS,dipS,rakeS] 
## vector containing the strike, dip, and rake of the stress tensor in degrees
## strikeS- strike of the sigmax-sigmay plane; dipS- dip angle of the sigmax-sigmay plane
## rakeS- the angle of line rake indicating the direction of sigmax within sigmax-sigmay plane

## Orientation of medium symmetry frame defined referred to the NED (North(x)-East(y)-Downward(z)) coordinate system
orientM=[0,0,0]##orientM=[strikeM,dipM,rakeM]
## vector containing the strike, dip, and rake of the stiffness matrix in degrees. 
## strikeM -strike of the foliation; dipM -dip angle of the foliation; rakeM -the angle of line rake indicating the direction of lineation
## Convert the geology angles into Euler rotation angles (See Fig. 3 in section 2.1)
alphazm1= orientM[0]*np.pi/180 - np.pi/2.0# The rotation angle of the medium symmetry frame around D-axis
alphaym= orientM[1]*np.pi/180# The rotation angle of the medium symmetry frame around y'-axis
alphazm2= np.pi/2.0 - orientM[2]*np.pi/180# The rotation angle of the medium symmetry frame around zm-axis
alphazs1= orientS[0]*np.pi/180 - np.pi/2.0# The rotation angle of the far-field stress frame around D-axis
alphays= orientS[1]*np.pi/180# The rotation angle of the far-field stress frame around y''-axis
alphazs2= np.pi/2.0 - orientS[2]*np.pi/180# The rotation angle of the far-field stress frame around zs-axis


############### Angle Option 2: material and far-field stress reference frames ###################
## Angle of the principle stress coordinate system from NEV in degrees
## NRS σ1 plunging 3° and trending 37° and σ3 trending 131°.
##NRS for 17 is azi= 37 inclin=10
principle1deg=145## azimuth alphazs
principle2deg=13 ## inclination  alphays
### This is the anisotropy of the medium angles in degress

symmetry1deg=0 ## azimuth
symmetry2deg=0 ## inclination
## Convert to Radians 
principle1=principle1deg*np.pi/180 ## azimuth
principle2=principle2deg*np.pi/180 ## inclination 
symmetry1=symmetry1deg*np.pi/180 ## azimuth
symmetry2=symmetry2deg*np.pi/180 ## inclination 

########################## Rotations of References Frames #############################
### Rotate the elastic stiffnesses from the medium symmetric coordinate system to borehole coordinate system
### Rotate the far-field stress from the in situ stress coordinate system to borehole coordinate system

########### option=1 or 2 select depending on what angle option is wanted
option=2
#### Reference Frame Rotation for angle option 1 ####
if option==1:
    stress=StressRotationEuler(stress1, alphazs1, alphays, alphazs2, borehole1, borehole2)
    C=StiffnessRotationEuler(C1, alphazm1, alphaym, alphazm2, borehole1, borehole2)
    
#### Reference Frame Rotation for angle option 2 ####
if option==2:
    stress=StressRotation3D(stress1, principle1, principle2, borehole1, borehole2)
    C=StiffnessRotation3D(C1, symmetry1, symmetry2, borehole1, borehole2)


############################# Main Solution #####################################
#### Once all reference frames are rotated into the borehole reference frame 
sigma1=stress[0,0]
sigma2=stress[1,1]
sigma3=stress[2,2]
sigma4=stress[1,2]
sigma5=stress[0,2]
sigma6=stress[0,1]

#################### Plotting Parameters #########################
X3=np.linspace(-normrad,normrad,sc)
Y1=np.linspace(-normrad,normrad,sc)
oldx2,oldy2=np.meshgrid(X3,Y1)
         
x2,y2=nan_borehole(oldx2,oldy2,sc,a) 
R=np.sqrt((x2**2)+(y2**2)) #Radius away from the borehole
theta=np.arctan2(y2,x2) ## Angle around borehole in radians
theta2=theta*180/np.pi## Angle around borehole in degrees
######################################################################


A=C**(-1) ## Compliance Matrix, Eqn. (13)
Beta=beta(A)  ## the matrix of reduced strain coefficients, Eqn. (12)

########## Roots and are complex or imaginary. 
mius=miu(Beta) ## General Case Mius, Eqn. (C.1)

l2,l3,l4,mius=lsoln(Beta,mius) ## linear differential operators,Eqn. (C.2),
## and new mius if it is a special case
miu1=mius[0]; miu2=mius[1]; miu3=mius[2] 

lamda=[-l3[0]/l2[0], -l3[1]/l2[1], -l3[2]/l4[2]]## lamdas, Eqn. (16) 
lamda1=lamda[0]; lamda2=lamda[1]; lamda3=lamda[2]

laplace=(miu2-miu1+(lamda2*lamda3*(miu1-miu3))+(lamda1*lamda3*(miu3-miu2))) ## Eqn. (20)

z1=x2+miu1*y2 ##Equation (15)
z2=x2+miu2*y2
z3=x2+miu3*y2

zeta1=zeta(z1,a,miu1,laplace,sc)## Eqn. (D.6)
zeta2=zeta(z2,a,miu2,laplace,sc)
zeta3=zeta(z3,a,miu3,laplace,sc)

#################### The analytical function, Appendix D #####################
phi1=1j*sigma6-sigma2+pw 
phi2=sigma6-1j*sigma1+1j*pw
phi3=sigma4-1j*sigma5

## Eq. (18)
dphi1=-(phi1*(miu2-lamda2*lamda3*miu3)+phi2*(lamda2*lamda3-1)+\
phi3*lamda3*(miu3-miu2 ))/zeta1

dphi2=-(-phi1*(miu1-lamda1*lamda3*miu3)-phi2*(lamda1*lamda3-1)+\
phi3*lamda3*(miu1-miu3))/zeta2

dphi3=-(phi1*(miu1*lamda2-miu2*lamda1)+phi2*(lamda1-lamda2)+\
phi3*(miu2-miu1))/zeta3
        
################################ Calculating Stresses #################################
## Eqn.(17)
istrbore1=2*( miu1**2*dphi1+ miu2**2*dphi2+ lamda3*miu3**2*dphi3).real
istrbore2=2*(dphi1+ dphi2+ lamda3*dphi3).real
istrbore4=-2*(lamda1*dphi1+ lamda2*dphi2+ dphi3).real
istrbore5=2*(lamda1*miu1*dphi1+ lamda2*miu2*dphi2+ miu3*dphi3).real
istrbore6=-2*(miu1*dphi1+ miu2*dphi2+ lamda3*miu3*dphi3).real
## Eqn.(B6)
istrbore3=-((A[2,0]*istrbore1)+(A[2,1]*istrbore2)+(A[2,3]*istrbore4)+\
            (A[2,4]*istrbore5) +(A[2,5]*istrbore6))/A[2,2]


### Cartesian Stresses Eqn.(21)
sigmax=sigma1+istrbore1
sigmay=sigma2+istrbore2
sigmaz=sigma3+istrbore3
taoyz=sigma4+istrbore4
taoxz=sigma5+istrbore5
taoxy=sigma6+istrbore6

### Cylindrical Stresses
##  Eqn.(22) used to convert the stress into cylindrical coordinate
hoop=(((np.sin(theta))**2)*(sigmax))+(((np.cos(theta))**2)*(sigmay))-((np.sin(theta*2))*(taoxy))
radial=(sigmax*((np.cos(theta))**2))+(((np.sin(theta))**2)*sigmay)+(np.sin(2*theta)*taoxy)
axial=sigmaz
taotz=((np.cos(theta))*taoyz)-((np.sin(theta))*taoxz)
taorz=((np.sin(theta))*taoyz)+((np.cos(theta))*taoxz)
taort=((sigmax)*((-0.5)*(np.sin(2*theta))))+((sigmay)*((0.5)*(np.sin(2*theta))))+((taoxy)*((np.cos(2*theta))))

#################################### Plotting Figures ####################################


### Cylindrical Stress ###
##
plottin(hoop,theta,R,normrad,'Hoop',1)
plottin(radial,theta,R,normrad,'Radial',2)
plottin(axial,theta,R,normrad,'Axial',3)
plottin(taotz,theta,R,normrad,'Tautz',4)
plottin(taorz,theta,R,normrad,'Taurz',5)
plottin(taort,theta,R,normrad,'Taurt',6)


#### Cartesian Stress ### 
##
#plottin2(sigmax,x2,y2,normrad,'Sigma X',7)
#plottin2(sigmay,x2,y2,normrad,'Sigma Y',8)
#plottin2(taoyz,x2,y2,normrad,'Taoyz',9)
#plottin2(taoxz,x2,y2,normrad,'Taoxz',10)
#plottin2(taoxy,x2,y2,normrad,'Taoxy',11)
#plottin2(sigmaz,x2,y2,normrad,'Sigma Z',12)

### Plotting the borehole wall stresses ###
### Solving for just the wall of the borehole ###
#
wall,hoopwall=borehole_wall360(theta2,hoop,R,a)
wall,radialwall=borehole_wall360(theta2,radial,R,a)
wall,axialwall=borehole_wall360(theta2,sigmaz,R,a)
wall,taotzwall=borehole_wall360(theta2,taotz,R,a)
wall,taorzwall=borehole_wall360(theta2,taorz,R,a)
wall,taortwall=borehole_wall360(theta2,taort,R,a)
plottin3(radialwall,hoopwall,axialwall,taotzwall,taorzwall,taortwall,wall,'Stress around Mine-By Tunnel Wall')
