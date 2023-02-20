# -*- coding: utf-8 -*-
"""
This module has various ways of solving the Held-Hou problem
Each is a function which can be run with default parameters or by setting each 
in turn

The best way to use it is to run from a driving code which you have example of 
in the exercises directory or:

import held_hou
held_hou.dry(H=20000)

Created on Mon Nov 16 05:59:47 2015

@author: andrewcharlton-perez
"""

def dry(H=10000, a=6370E3, omega=7.292E-5, delta_theta=50,   \
                 theta_0=300, g=9.81, tau=10*86400, delta_thetav=30):   

 """
 This function solves the simple zonally symmetric, dry Held-Hou model

 Inputs are:

 H - scale height (m)
 a - radius of planet (m)
 omega - rotation rate of planet (s-1)
 delta_theta - difference in temp, equator to pole (K)
 theta_0 - base temperature at mid-point (K)
 g - acceleration due to gravity (ms-2)
 delta_thetav - temperature difference between layers (K)


 Default values are for the present day Earth can run with no options for
 standard case or change each from default

 Outputs are:

 latitude
 theta_e
 theta_m u (which are all vectors) and

 edge of cell (m)
 edge of cell in latitude
 difference in potential temp
 at equator
 scale for v (ms-1)
 scale for w (ms-1)
 """



 import numpy as np

# calculate independent variable y based on planet radius
 phi=np.linspace(-40,40,400)
 y=a*np.radians(phi)

# solve for unknowns (Y and theta_m0-theta_e0)
 Y=np.sqrt((5*delta_theta*g*H)/(3*(omega**2)*theta_0))
 delthe_eq=(5*(delta_theta**2)*g*H)/(18*(a**2)*(omega**2)*theta_0)

# potential temperature profiles
 theta_e0=theta_0+delta_theta/3.
 theta_m0=theta_e0-delthe_eq
 
 theta_e=theta_e0-(delta_theta/(a**2))*y**2
 theta_m=theta_m0-(((omega**2)*theta_0)/(2*(a**2)*g*H))*y**4

# cell edge in degrees
 PHI=np.degrees(np.arcsin(Y/a))

# calculating velocities
# zonal
 u=(omega*y**2)/a
# outside cell velocity is in themral wind balance
 u[np.where(np.abs(y) > Y)]=(g*H*delta_theta)/(omega*theta_0*a)

# vertical
 w=(-float(H)/(delta_thetav*tau))*(theta_m0-theta_e0)

# meridional
 v=Y*(w/H)

 return(phi,theta_e,theta_m,u,Y,PHI,delthe_eq,v,w)

####################################


def off(heating=6,Yw_fac=-1.2,Ys_fac=1.2,Y1_fac=1.2,delthe_fac=1.2,
        delta_theta=50):

 """
 This function iteratively solves the off equatorial version of the Held-Hou model
 The standard Earth values (as above) are used as constants here

 The single optional value used is the position of off-equatorial heating in degrees (other optional
 values used above are tuning constants for the optimization

 outputs are latitude, edge of winter cell, edge of summer cell, mid-point of
 cells (position of max upwelling), theta_m1 (the constant for the theta_m
 function), theta_e, theta_m, the vertical and horizontal velocity scales
 and then two messages from the iteration routine
 to check convergence, the first is True is converged and triggers a warning
 message if it didn't

 """



 import numpy as np
 from numpy import cos,sin
 from scipy.optimize import root
 import matplotlib.pyplot as plt

 

# use standard Earth values as constants but can be changed
 H=15000 
 a=6370.E3 
 omega=7.292E-5 
 delta_theta=delta_theta
 theta_0=300 
 g=9.81 
 tau=10*86400 
 delta_thetav=30
 delta_h=1./6.

# independent variable phi 
 y=np.arange(-4000,4010,10)*1E3
 phi=np.degrees(np.arcsin(y/a))
 Y0=a*sin(np.radians(heating))
 
# centred case (for checking of code and constraints)
 Y=np.sqrt((5*delta_theta*g*H)/(3*(omega**2)*theta_0))
 delthe_eq=(5*(delta_theta**2)*g*H)/(18*(a**2)*(omega**2)*theta_0) 

 # define functional form of theta_m and theta_e
 theta_e0=theta_0+(delta_theta/3.)
 theta_m0=theta_e0-delthe_eq
 M=(((omega**2)*theta_0)/(2*(a**2)*g*H))
 E=delta_theta/(a**2)
 
# integrals for centred case
# int_e_cen=(theta_e0*Y)-(E/3)*Y**3
# int_m_cen=(theta_m0*Y)-(M/5)*Y**5

 
 def theta_e(y,Y0):
     return(theta_e0-E*((y-Y0)**2))

 def theta_m(y,Y1,theta_m1):
     return(theta_m1-M*((y**2)-(Y1**2))**2)

# define constraints function
 def constraints(p):
  Yw,Ys,Y1,the_diff=p
         
# each constraint in turn
# edges      
  f1=-E*(Y0 - Ys)**2 + M*(Y1**2 - Ys**2)**2 + the_diff
  f2=-E*(Y0 - Yw)**2 + M*(Y1**2 - Yw**2)**2 + the_diff
 
  int_m_s=-M*(8*Y1**5 - 15*Y1**4*Ys + 10*Y1**2*Ys**3 - 3*Ys**5)/15      
  int_m_w=-M*(8*Y1**5 - 15*Y1**4*Yw + 10*Y1**2*Yw**3 - 3*Yw**5)/15       
 
  int_e_s=-E*(3*Y0**2*Y1 - 3*Y0**2*Ys - 3*Y0*Y1**2 + 3*Y0*Ys**2 + Y1**3 - Ys**3)/3
  int_e_w=-E*(3*Y0**2*Y1 - 3*Y0**2*Yw - 3*Y0*Y1**2 + 3*Y0*Yw**2 + Y1**3 - Yw**3)/3
  
# integral constraints      
  f3=the_diff*(Ys-Y1)-(int_e_s-int_m_s)
  f4=the_diff*(Yw-Y1)-(int_e_w-int_m_w)
          
  return(f1,f2,f3,f4)

 
### Iterating section here 
 
# set up starting points based on optional input mutlipliers of symmetric case
 sp=(Yw_fac*Y,Ys_fac*Y,Y1_fac*Y0,delthe_fac*delthe_eq)

# iterate
 sol=root(constraints,sp,method='hybr')  
 
 
# parse outputs of iteration here (but be sure to check message pass in case 
# no convergence) 
 Yw,Ys,Y1,the_diff=sol.x

# Yw,Ys,Y1,the_diff=sp
 PHIw=np.degrees(np.arcsin(Yw/a))
 PHIs=np.degrees(np.arcsin(Ys/a))
 PHI1=np.degrees(np.arcsin(Y1/a))  
 theta_m1=theta_e0-the_diff
 

# vertical velocity at centre of cell
 w=(-float(H)/(delta_thetav*tau))*(theta_m1-theta_e(Y1,Y0))

# meridional velocity in each cell
 vw=(Y1-Yw)*(w/H) 
 vs=(Y1-Ys)*(w/H)
 
 
 PHI0=heating
 return(phi,PHIw,PHIs,PHI1,PHI0,  \
        theta_m1,theta_e(y,Y0),theta_m(y,Y1,theta_m1),  \
        w,vw,vs,  \
        sol.success,sol.message)



###################################


def moist(width=0.3, H=10000, a=6370E3, omega=7.292E-5, delta_theta=50,   \
          theta_0=300, g=9.81, tau=10*86400, delta_thetav=30):

 """"
 This function iteratively solves for the mosit version of the Held-Hou model
 There is only one input paramter in the experiment (although the constants can
 be changed above) which is a scale for the width of the moist heating as a
 fraction of the width of the dry cell.
 
 Outputs are latitude, latitude of cell edge, latitude of edge of heating
 theta_m, theta_e, theta_estar, scale for v and scale for w
 """
     

 import numpy as np
 from scipy.optimize import root
 import math
 import matplotlib.pyplot as plt

# independent variable (y) in m
 y=np.arange(-4000,4010,10)*1E3
# for plotting convert this to latitude
 phi=np.degrees(np.arcsin(y/a))

# symmetric case
 Y=np.sqrt((5*delta_theta*g*H)/(3*(omega**2)*theta_0))
 delthe_eq=(5*(delta_theta**2)*g*H)/(18*(a**2)*(omega**2)*theta_0)
 theta_e0=theta_0+delta_theta/3.
 theta_m0=theta_e0-delthe_eq
 
# parameters of equation 
 M=(((omega**2)*theta_0)/(2*(a**2)*g*H))
 E=(delta_theta/(a**2)) 
 
# define functional form of potential temperatures
 def theta_estar(y,amp):
     return(theta_estary+amp*(np.exp(-(y**2)/(w_m**2))))
     
 def theta_e(y):
     return(theta_e0-E*(y**2))

 def theta_m(y):
     return(theta_m0-M*(y**4))    

# fixed parameters of moist solution
# width is a fraction of the dry case
 w_m=width*Y
# outside of heating region theta_estar is flat and equal to edge value            
 theta_estary=theta_m(Y)           
     
# define constraint function just constraining amplitude based on width
 def constraint(amp):
# match_area
  int_theta_m=(theta_m0*Y)-(M*Y**5)/5 
  int_theta_estar=theta_estary*Y + \
    np.sqrt(np.pi)*amp*math.erf(Y*np.sqrt(w_m**(-2)))/(2*np.sqrt(w_m**(-2))) \
 
  return(int_theta_estar-int_theta_m)
 
# find solution         
 sol=root(constraint,10,method='hybr')   
 amp=sol.x

# cell edge in degrees
 PHI=np.degrees(np.arcsin(Y/a))

# find width of heating in degrees
 zc=np.ravel(np.where((theta_estar(y[0:len(y)-1],amp) > theta_m(y[0:len(y)-1]))  \
        & (theta_estar(y[1:len(y)],amp) < theta_m(y[1:len(y)]))))
 PHIw=np.degrees(np.arcsin(y[zc[1]]/a))

# equilibrium potential temperature at equator is equal to amplitude
 theta_estar_eq=amp+theta_estary

 
# vertical velocity
 w=(-float(H)/(delta_thetav*tau))*(theta_m0-theta_estar_eq)

# meridional velocity
 v=Y*(w/H)

 return(phi,PHI,PHIw,theta_m(y),theta_e(y),theta_estar(y,amp),v,w)
  
  
##################  
 


#%%

#%%
