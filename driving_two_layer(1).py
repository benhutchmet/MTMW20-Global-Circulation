# -*- coding: utf-8 -*-
"""

This shows you how to use the function to run the rather complex code needed to
solve the two-layer shallow water model for the mid-lats

Inputs are:
Time of integration - you may need to change this for some parameter settings
to get to a convergent state
default: ndays=200

Timescale for radiative damping back to equilibrium state in days
default: tau=40 

Timescale for friction to damp away jet in bottom layer in days
default: R=10

Diffusision coefficient (a squared length scale) in m^-2
default: K1=5E5

The model has a very large number of parameters and if you like you can
change them in the function

Outputs are:
ym - mid-point grid in degrees latitude for PV and the interface height
y- grid in degrees latitude for zonal wind, PV flux and residual vel.
time - in days for the values which are saved (not all steps are 1 per day)
u - zonal wind (all outputs have the structure (latitude, level, time))
q - PV
eta - interface height (only one value)
flux - PV flux
resid - residual velocity

and then some parameters used in the model (as scaled within)
R - friction coefficient
tau - damping coefficient
eta_e - equilibrium position of interface
K1 - diffusion coefficient in layer 1
K2 - diffusion coefficient in layer 2
f0 - coriolis parameter at the centre of the channel

Created on Sun Nov 29 21:07:16 2015

@author: andrewcharlton-perez
"""

import numpy as np
import matplotlib.pyplot as plt
import two_layer_prog

y,ym,time,u,q,eta,flux,resid,  \
           R,tau,eta_e,K1,K2,f0=two_layer_prog.solve()

# example plots showing evolution of fields in each layer
# two types of plot are generated

# plot diagnostics to show state

for l in range(2):
  fig=plt.figure(figsize=(7,9))
  plt.suptitle('Level '+str(l+1))
  plt.subplot(321)
  plt.contourf(time,ym,q[:,l,:].squeeze()/f0,cmap=plt.get_cmap('RdBu_r'))
  plt.colorbar()
  plt.title('PV (scaled f$_{0}$)')
  plt.yticks([-45,-25,0,25,45])
 
  plt.subplot(322)
  plt.contourf(time,ym,eta[:,:].squeeze(),cmap=plt.get_cmap('RdBu_r'))
  plt.colorbar()
  plt.title('Interface height') 
  plt.yticks([-45,-25,0,25,45])
 
  plt.subplot(324)
  mp=np.ceil(np.max(u[:,l,:]))
  levs=np.linspace(-mp,mp,10)
  plt.contourf(time,y,u[:,l,:].squeeze(),cmap=plt.get_cmap('RdBu_r'),levels=levs)
  plt.colorbar()
  plt.title('Zon. wind') 
  plt.yticks([-45,-25,0,25,45])
 
  plt.subplot(323)
  plt.contourf(time,y,resid[:,l,:].squeeze(),cmap=plt.get_cmap('RdBu_r'))
  plt.colorbar()
  plt.title('Residual velocity') 
  plt.yticks([-45,-25,0,25,45])
 
  plt.subplot(325)
  plt.contourf(time,y,flux[:,l,:].squeeze()/f0,cmap=plt.get_cmap('RdBu_r'))
  plt.colorbar()
  plt.title('PV Flux (scaled f$_{0}$)')  
  plt.yticks([-45,-25,0,25,45])
 
  plt.tight_layout() 
  plt.savefig('two_layer_evolution_'+str(l+1)+'.png')