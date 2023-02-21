# -*- coding: utf-8 -*-
"""
This selection of functions solves the zonal mean two-layer shallow water
equations on a beta channel (see Vallis book) for use in global circulation

For details on in and out see the driving code in the same directory
(driving_two_layer.py)

@author: andrewcharlton-perez
"""
import numpy as np
import matplotlib.pyplot as plt


# main solving function
def solve(ndays=200, tau=40, R=10, K1=5E5):

# setup model
# staggered grid (pv, eta on full points, u, flux, resid vel. on half points)
#######################

# channel parameters (keep fixed)
 global a, y, dy, ym, f0, ny, nym, p, amp 
 width=90   # width of channel in degrees
 chan_cent=45
 omega=7.292E-5
 a=6370E3
 f0=2*(7.292E-5)*np.sin(np.radians(chan_cent))    # channel centred at 45 N

# setup the grid
 width_y=a*np.tan(np.radians(width/2))
 y=np.linspace(-width_y,width_y,num=width)
 ny=len(y)
 ym=(y[0:-1]+y[1:])/2.   # mid=points for staggered grid
 nym=len(ym)
 dy=y[1]-y[0]
 p=0.5*np.pi/width_y
 
# mean depth of layers
 global H, H2, H1, gp
 H=10000  # total depth
 H2=5000  # bottom layer depth
 H1=H-H2  # top layer depth
 gp,rho1,rho2=red_grav(H1,H2)
 
##############################
# parameters
# convert relaxation rate to seconds
 tau=tau*86400. 
 # slope of background height is for constant shear of 30ms-1
 slope=30*f0/gp
 eta_e=slope*ym
  
# friction parameter
 R=1/(R*86400.)
  
# diffusion parameter for bottom layer from kinematic constraint using i.c.s
 K2=(((H1*gp*(p**2))/(f0**2))+1)*K1

###############################
# time parameters
 end_time=ndays*86400
 steps_per_day=32
 dt=86400/steps_per_day  # 16 timesteps per day
 nt=ndays*steps_per_day
 
# store output once per day
 step_store=np.arange(0,steps_per_day*(ndays+1),steps_per_day)

# storage arrays
 u_st=np.zeros((len(y),2,ndays+1))
 q_st=np.zeros((len(ym),2,ndays+1))
 flux_st=np.zeros((len(y),2,ndays+1))
 vR_st=np.zeros((len(y),2,ndays+1))
 eta_st=np.zeros((len(ym),ndays+1))

#############################
# initial conditions 

# zonal winds (peak at 40)
 amp=40
 u1=amp*np.cos(p*y)
 u2=np.zeros((len(y)))
 u_st[:,:,0]=np.vstack((u1,u2)).transpose()
# interface
 eta=((f0*amp)/(gp*p))*np.sin(p*ym)
 eta_st[:,0]=eta


# PV 
 q1=((amp*p)+(((f0**2)*amp)/(H1*p*gp)))*np.sin(p*ym)
 q2=-(((f0**2)*amp)/(H2*p*gp))*np.sin(p*ym)
 q_st[:,:,0]=np.vstack((q1,q2)).transpose()

# PV flux 
 flux1=np.hstack((0,-K1*np.diff(q1)/dy,0))
 flux2=np.hstack((0,-K2*np.diff(q2)/dy,0))
 flux_st[:,:,0]=np.vstack((flux1,flux2)).transpose()


# Residual velocity
 S=(eta_e-eta)/tau
 vR1,vR2=invert_resid(np.diff(S)/dy,   \
                     flux1,flux2,u2,R)
                                         
 vR_st[:,:,0]=np.vstack((vR1,vR2)).transpose()
  
 

###############################################################

# now solve prognostic equations in PV using FTCS scheme
# boundary condition is no flux of PV at boundaries 
 
 for t in np.arange(1,nt+1):
# first time-stepping PV and mass      
    
# layer 1 (top) - solve for interior points     
# boundary conditions of no eddy pv flux at boundary    
    q1_past=q1
    trans=-np.diff(flux1)/dy
    rad=(f0/H1)*S  
    q1=q1_past+dt*(trans+rad)


# bottom layer - solve for interior points
    q2_past=q2    
    trans=-np.diff(flux2)/dy
    rad=(f0/H2)*S
    fric=R*np.diff(u2)/dy    
    q2=q2_past+dt*(trans-rad+fric)

# solve for eta change (FTCS) - bc is no residual velocity at boundary
    eta_past=eta
    resid_div=np.diff(vR1)/dy
    eta=eta_past+dt*(S+(H1*resid_div))
 
# inversion steps to get new state
# wind
    u1,u2=int_u2(q1,q2,eta)    
    
# PV flux
    flux1=np.hstack((0,-K1*np.diff(q1)/dy,0))
    flux2=np.hstack((0,-K2*np.diff(q2)/dy,0))     
        
# radiative term
    S=(eta_e-eta)/tau    
            
# residual velocity
    vR1,vR2=invert_resid(np.diff(S)/dy,   \
                         flux1,flux2,u2,R)
    
       
# if on a whole number day store  
    if t in step_store:
        d=np.int(t/steps_per_day)
        u_st[:,:,d]=np.vstack((u1,u2)).transpose()
        q_st[:,:,d]=np.vstack((q1,q2)).transpose()
        eta_st[:,d]=eta
        vR_st[:,:,d]=np.vstack((vR1,vR2)).transpose()
        flux_st[:,:,d]=np.vstack((flux1,flux2)).transpose()
        
    
#############
# return steps

 yp=np.degrees(np.arctan(y/a)) 
 ymp=np.degrees(np.arctan(ym/a))
 tp=(step_store*dt)/86400.


# return all of the outputs but also the constants used
 return(yp,ymp,tp,u_st,q_st,eta_st,flux_st,vR_st,   \
        R,tau,eta_e,K1,K2,f0)


############################################

# other functions

def red_grav(H1,H2):
# based on input layer heights calculate the density at the mid-point of each level
# US standard atmosphere
# this gives reduced gravity (gp)
 T2=288.15-((6.5E-3)*(H2/2.))
 P2=101325*(1-2.25569E-5 * (H2/2.))**5.25616
 RHO2=P2/(287.04*T2)
 T1=288.15-((6.5E-3)*(H2+H1/2.))
 P1=101325*(1-2.25569E-5 * (H2+H1/2.))**5.25616
 RHO1=P1/(287.04*T1)
 gp=9.80655*(RHO2-RHO1)/RHO1
 return(gp,RHO1,RHO2)

# use thermal wind balance to calc. interface height for interior points
def thermal_wind(u_top,u_bot,dy):
 sh=(f0/gp)*(u_top-u_bot)
 return((sh[0:-1]+sh[1:])*dy)

# invert to get the residual velocity
def invert_resid(dSdy,v1q1,v2q2,u2,R):
 from numpy.linalg import inv   
 F=gp*dSdy-f0*(v1q1[1:-1]-v2q2[1:-1])-f0*R*u2[1:-1]

# tri-diagonal matrix
 c=((f0**2)*H/H2)+2*H1*gp/(dy**2)
 l=-H1*gp/(dy**2)
 M=np.diag(np.repeat(c,len(y)-2),k=0) \
  +np.diag(np.repeat(l,len(y)-3),k=-1) \
  +np.diag(np.repeat(l,len(y)-3),k=1)
  
 vR1=np.dot(inv(M),F)

 vR2=(-H1*vR1)/H2
  
# b.c.'s residual velocity zero at boundaries
 vR1=np.hstack((0,vR1,0))
 vR2=np.hstack((0,vR2,0))  
  
 return(vR1,vR2)
 
########## 
 
# invert to get lower layer velocity 
# (interior points and bcs of no flow)
def int_u2(q1,q2,eta):
 from numpy.linalg import inv    

# solve for u 
 F1=dy*((f0*eta/H1)-q1)
 F2=dy*(-q2-f0*eta/H2)

 u1=np.zeros((len(y)))
 u2=np.zeros((len(y)))
 for i in np.arange(1,len(y)):
     u1[i]=u1[i-1]+F1[i-1]
     u2[i]=u2[i-1]+F2[i-1]    

 return(u1,u2)


#############################################






