# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 18:31:27 2018

@author: Maria
"""
import numpy as np
from numpy.linalg import inv

####################Velocityto stiffness Example #####################
###Li, 2017
#density=2.767 ## the density of the media
#vp90=4.76727 ## 90 degree P-wave phase velocty(parallel to the bedding)
#vp0=3.47572 ##0 degree P-wave phase velocty(vertical to the bedding)
#vsh0=2.63688 ## 90 degree SH-wave phase velocty(parallel to the bedding)
#vsh90=3.10092 ## 0 degree SH-3wave phase velocty(vertical to the bedding)
#vp45=4.06808 ## 45 degree P-wave phase velocty
#inputs=1
#C1vel=Velocitytostiffness(vp0, vsh0, vp90, vsh90, vp45,density)


#######################################################################

def StrRMatrixY(x): 
##The function calculates the matrix for rotating the stress tensor around y-axis.
##Inputs: x- rotation angle.
##Outputs: M- the matrix used to rotate the stress tensor around y-axis (E-axis in the NEV coordinates).
    M=np.matrix([[np.cos(x), 0, -np.sin(x)],
    [0, 1, 0],
    [np.sin(x), 0, np.cos(x)]]) 
    return M
        

def StrRMatrixZ(x):
##The function calculates the matrix for rotating the stress tensor around z-axis.
##Inputs: x- rotation angle.
##Outputs: M- the matrix used to rotate the stress tensor around z-axis (V-axis in the NEV coordinates).
    M=np.matrix([[np.cos(x), np.sin(x), 0],
    [-np.sin(x), np.cos(x), 0],
    [0, 0, 1]])
    return M

def BondY(alpha):
##The function calculates the matrix for rotating the elastic stiffness matrix around y-axis.
##Inputs: alpha- rotation angle.
##Outputs: M- the matrix used to rotate elastic stiffness matrix around y-axis (E-axis in the NEV coordinates).

    M=np.zeros([6,6])
    M[0,0]=np.cos(alpha)**2
    M[1,1]=1
    M[2,2]=np.cos(alpha)**2
    M[3,3]=np.cos(alpha)
    M[4,4]=np.cos(2*alpha)
    M[5,5]=np.cos(alpha)
    M[0,2]=np.sin(alpha)**2
    M[2,0]=np.sin(alpha)**2
    M[0,4]=-np.sin(2*alpha)
    M[4,0]=(1./2.)*np.sin(2*alpha)
    M[2,4]=np.sin(2*alpha)
    M[4,2]=-(1./2.)*np.sin(2*alpha)
    M[3,5]=np.sin(alpha)
    M[5,3]=-np.sin(alpha)
    return np.matrix(M)


def BondZ(alpha):
##The function calculates the matrix for rotating the elastic stiffness matrix around z-axis.
##Inputs:alpha- rotation angle.
##Outputs:M- the matrix used to rotate elastic stiffness matrix around z-axis (V-axis in the NEV coordinates).
    M=np.zeros([6,6])
    M[0,0]=np.cos(alpha)**2
    M[1,1]=np.cos(alpha)**2
    M[2,2]=1
    M[3,3]=np.cos(alpha)
    M[4,4]=np.cos(alpha)
    M[5,5]=np.cos(2*alpha)
    M[0,1]=np.sin(alpha)**2
    M[1,0]=np.sin(alpha)**2
    M[0,5]=np.sin(2*alpha)
    M[5,0]=-(1./2.)*np.sin(2*alpha)
    M[1,5]=-np.sin(2*alpha)
    M[5,1]=(1./2.)*np.sin(2*alpha)
    M[3,4]=-np.sin(alpha)
    M[4,3]=np.sin(alpha)
    return np.matrix(M)
    
def StressRotation3D(stress, principle1, principle2, borehole1, borehole2):
##The function rotates the stress tensor from the principle stress coordinate system to the borehole coodinate system
##Inputs:
    ##stress- the stress tensor in the principle stress frame
    ##principle1 and principle2 are the azimuth and tilt angle of the principle stress frame referred to NEV frame, respectively
    ##borehole1 and borehole2 are the azimuth and tilt angle of the borehole frame referred to NEV frame, respectively.
##Outputs:
    ##istress- the stress tensor in the borehole coordinate frame.
    M1=StrRMatrixZ(principle1)
    M2=StrRMatrixY(principle2)
    M3=StrRMatrixZ(borehole1)
    M4=StrRMatrixY(borehole2)
    M1t=M1.transpose()
    M2t=M2.transpose()
    M3t=M3.transpose()
    M4t=M4.transpose()
    istress=M4*M3*(inv(M1))*(inv(M2))*stress*(inv(M2t))*(inv(M1t))*M3t*M4t ##equation (6)
    return istress

def StressRotationEuler(stress, alphazs1, alphays, alphazs2, borehole1, borehole2):
##The function rotates the stress tensor from the principle stress coordinate system to the borehole coodinate system
#Inputs:
    ##stress- the stress tensor in the far-field stress frame
    ##alphazs1 The rotation angle of the far-field stress frame around D-axis
    ##alphays The rotation angle of the far-field stress frame around y'-axis
    ##alphazs2 The rotation angle of the far-field stress frame around zm-axis
    ##borehole1 and borehole2 are the azimuth and tilt angle of the borehole frame referred to NEV frame, respectively.
##Outputs:
    ##istress- the stress tensor in the borehole coordinate frame.
    N1=StrRMatrixZ(-alphazs1)# Function StrRMatrixZ is the rotation of stress around Z-axis
    N2=StrRMatrixY(-alphays)# Function StrRMatrixY is the rotation of stress around Y-axis
    N3=StrRMatrixZ(-alphazs2)
    N4=StrRMatrixZ(borehole1)
    N5=StrRMatrixY(borehole2)
    N1t=N1.transpose()
    N2t=N2.transpose()
    N3t=N3.transpose()
    N4t=N4.transpose()
    N5t=N5.transpose()
    istress=N5*N4*N1*N2*N3*stress*N3t*N2t*N1t*N4t*N5t ##equation (6)
    return istress

def StiffnessRotation3D(C1, symmetry1, symmetry2, borehole1, borehole2):
##The function rotates the elastic stiffness matrix from the mmedium symmetry coordinate system to the borehole coordinate system
##Inputs:
    ##C1- the elastic stiffness matrix in the media symmetric frame
    ##symmetry1 and symmetry2 are the azimuth and tilt angle of the medium symmetry frame referred to NEV frame, respectively.
    ##borehole1 and borehole2 are the azimuth and tilt angle of the borehole frame referred to NEV frame, respectively.
#Outputs:
    ##C- the elastic stiffness matrix in the borehole coordinate frame.
    M1=BondZ(symmetry1)
    M2=BondY(symmetry2)
    M3=BondZ(borehole1)
    M4=BondY(borehole2)
    M1t=M1.transpose()
    M2t=M2.transpose()
    M3t=M3.transpose()
    M4t=M4.transpose()
    C=M4*M3*(inv(M1))*(inv(M2))*C1*(inv(M2t))*(inv(M1t))*M3t*M4t ##equation (6)
    return C

def StiffnessRotationEuler(C1, alphazm1, alphaym, alphazm2, borehole1, borehole2):
## The function rotates the elastic stiffness matrix from the medium symmetry coordinate frame to the borehole coordinate frame
## Inputs:
    ##C1- the elastic stiffness matrix in the media symmetric frame
    ##alphazm1 The rotation angle of the medium symmetry frame around D-axis
    ##alphaym The rotation angle of the medium symmetry frame around y'-axis
    ##alphazm2 The rotation angle of the medium symmetry frame around zm-axis
    ##borehole1 borehole2 are the azimuth and inclination of the borehole.
    ##All the angles are input referred to the NED frame.
#Outputs:
    #C- the elastic stiffness matrix in the borehole coordinate frame.
    M1=BondZ(-alphazm1)# Function BondZ is the Bond rotation around Z-axis
    M2=BondY(-alphaym) # Function BondY is the Bond rotation around Y-axis
    M3=BondZ(-alphazm2)
    M4=BondZ(borehole1)
    M5=BondY(borehole2)
    M1t=M1.transpose()
    M2t=M2.transpose()
    M3t=M3.transpose()
    M4t=M4.transpose()
    M5t=M5.transpose()
    C=M5*M4*M1*M2*M3*C1*M3t*M2t*M1t*M4t*M5t ##equation (6)
    return C

def Velocitytostiffness(vp0, vsh0, vp90, vsh90, vp45,density):
##The function generates elastic stiffness matrix from the phase velocities.
##Inputs:
    ##density- the density of the media;
    ##vp0- the P-wave phase velocity in the direction vertical to the foliation(0 degree);
    ##vsh0- the SH-wave phase velocity in the direction vertical to the foliation(0 degree);
    ##vp90- the P-wave phase velocity in the direction parallel to the foliation(90 degree);
    ##vsh90- the SH-wave phase velocity in the direction parallel to the foliation(90 degree);
    ##vsh45- the P-wave phase velocity in the 45 degree.
##Outputs:
    #C1-elastic stiffness matrix
#The function also provides the anisotropic coefficients (ita, delta, and gama). Please add them to the function output if you need.
    C1=np.zeros((6,6))
    C1[0,0]=density*vp90**(2)
    C1[1,1]=density*vp90**(2)
    C1[2,2]=density*vp0**(2)
    C1[3,3]=density*vsh0**(2)
    C1[4,4]=density*vsh0**(2)
    C1[5,5]=density*vsh90**(2)
    C1[0,1]=C1[0,0]-2*C1[5,5]
    C1[1,0]=C1[0,1]
    C1[0,2]=-C1[3,3]+np.sqrt( (C1[0,0]+C1[3,3]-2*density*vp45**2)*(C1[2,2]+C1[3,3]-2*density*vp45**2) )
    C1[2,0]=C1[0,2]
    C1[2,1]=C1[0,2]
    C1[1,2]=C1[0,2]
    ita=(C1[0,0]-C1[2,2])/(2.0 *C1[2,2])
    delta=( (C1[0,2]+C1[3,3])**2.0-(C1[2,2]-C1[3,3])**2.0 )/( 2*C1[2,2]*(C1[2,2]-C1[3,3]))
    gamma=(C1[5,5]-C1[3,3])/2.0/C1[3,3]
    return C1,ita,delta,gamma

def thompsenstiffness(vp0, vsh0, ita,delta,gamma,density):
##The function generates elastic stiffness matrix from the phase velocities.
##Inputs:
    ##density- the density of the media;
    ##vp0- the P-wave phase velocity in the direction vertical to the foliation(0 degree);
    ##vsh0- the SH-wave phase velocity in the direction vertical to the foliation(0 degree);
    ##vp90- the P-wave phase velocity in the direction parallel to the foliation(90 degree);
    ##vsh90- the SH-wave phase velocity in the direction parallel to the foliation(90 degree);
    ##vsh45- the P-wave phase velocity in the 45 degree.
##Outputs:
    #C1-elastic stiffness matrix
#The function also provides the anisotropic coefficients (ita, delta, and gama). Please add them to the function output if you need.

    c33=(density*(vp0**(2.0)))
    c44=(density*(vsh0**(2.0)))
    c11=((ita*2.0*c33)+c33)
    c66=((gamma*2.0*c44)+c44)
    c13=((((delta*2.0*c33*(c33-c44))+((c33-c44)**(2.0)))**(1.0/2.0))-c44)
    return c11,c33,c44,c66,c13