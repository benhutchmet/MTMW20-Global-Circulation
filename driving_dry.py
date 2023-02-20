# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 08:59:29 2015

@author: andrewcharlton-perez
"""

import numpy as np
import matplotlib.pyplot as plt
import held_hou


# this shows you how to run the Held-Hou function with standard options
# outputs are latitude, theta_e, theta_m u (which are all vectors) and
# edge of cell in m, edge of cell in latitude, difference in potential temp
# at equator, scale for v, scale for w

# you can change the following standard values in the options (put e.g. H=8000)
# H - depth
# a - radius of planet 
# omega - rotation rate of planet 
# delta_theta - scale of eq to pole temp. gradient
# theta_0 - reference pot. temp at mid-point 
# g - gravitational constant 
# tau - radiative damping timescale 
# delta_thetav - temperature difference between layers

# this line is all you need to run the model
phi,theta_e,theta_m,u,Y,PHI,delthe_eq,v,w=held_hou.dry()


##########

# This is an example of how you can plot this data
# make plot and add numbers for velocity scales etc. 
# get value at edge for scaling axis and then add another 20 points
edge_ind=np.argmin(np.abs(phi-PHI))+5

# plotting steps
fig=plt.figure(figsize=(10,5))
# theta values
ax=fig.add_subplot(1,2,1)
ax.plot(phi,theta_m,'k')
ax.plot(phi,theta_e,'k--')

# options
ax.set_xlim([-phi[edge_ind],phi[edge_ind]])
ax.set_ylim([theta_e[edge_ind],np.max(theta_e)])
ax.set_ylabel("Potential temperature / K")
ax.set_xlabel('Latitude')

# calculated values
mub=theta_e[edge_ind]
mut=np.max(theta_e)
delmu=mut-mub

ax.text(0,mub+.1*delmu,'Y $\sim$ '+"{0:.2g}".format(Y/1000)+' km',ha='center') 
ax.text(0,mub+.15*delmu,  \
  r'$\theta_{e0}$ - $\theta_{m0}$ $\sim$ '+"{0:.2g}".format(delthe_eq)+' K', \
  ha='center') 

#######

# u plot
ax=fig.add_subplot(1,2,2)
ax.plot(phi,u,'k')

# options 
ax.set_xlim([-phi[edge_ind],phi[edge_ind]])
ax.set_ylim([np.min(u),np.max(u)])
ax.set_ylabel('Zonal wind / ms$^{-1}$')
ax.set_xlabel('Latitude')


# calculated values
mu=np.max(u)
ax.text(0,mu-.1*mu,'v $\sim$ '+"{0:.2g}".format(v)+' ms$^{-1}$',ha='center') 
ax.text(0,mu-.15*mu,'w $\sim$ '+"{0:.2g}".format(w)+' ms$^{-1}$',ha='center') 

# this will show the figure on screen to save it hit save
plt.show()


