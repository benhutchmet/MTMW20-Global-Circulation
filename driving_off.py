# -*- coding: utf-8 -*-
"""

 This shows you how to run the Held-Hou function for the off equatorial case

 Setup for standard Earth parameters which you can change in the code if you
 wish

 Otherwise inputs are:
 heating - position of maximum heating in degrees

 and then inputs to do with the starting values (which are all relative to the
 symmetric case), this is quite sensitive since iterative you can change each
 or all of these for more extreme cases to see if you can get convergence
 there is a message if not. There is one factor for each input with standard
 values of 1.2 (20% more than the symmetric case)

 Yw_fac=-1.2
 Ys_fac=1.2
 Y1_fac=1.2
 delthe_fac=1.2


 outputs are latitude, edge of winter cell, edge of summer cell, mid-point of
 cells (position of max upwelling), theta_m1 (the constant for the theta_m
 function), theta_e, theta_m, the vertical and horizontal velocity scales
 and then two messages from the iteration routine
 to check convergence, the first is True is converged and triggers a warning
 message if it didn't

Created on Mon Nov 16 08:59:29 2015

@author: andrewcharlton-perez
"""

import numpy as np
import matplotlib.pyplot as plt
import held_hou


# this shows you how to run the Held-Hou off equatorial function
# this line is all you need to run the model
# phi,PHIw,PHIs,PHI1,PHI0,  \
#     theta_m1,theta_e,theta_m,  \
#       vert_sc,mer_w_sc,mer_s_sc,  \
#         conv,message=held_hou.off(heating=4)
#
# print('the value for phi0 is, ',PHI0)
# print('the value for phi1 is',PHI1)
#
#
# # flag for if solution hasn't converged so print warning message
# if conv != True:
#     print(message)
#
# ##########
#
# # This is an example of how you can plot this data
# # make plot and add numbers for velocity scales etc.
# # get value at edge for scaling axis and then add another 20 points
# edge_ind=np.argmin(np.abs(phi-PHIw))-40
# if edge_ind < 0: edge_ind=0
#
# # plotting steps
# fig=plt.figure(figsize=(10,5))
# # theta values
# ax=fig.add_subplot(1,1,1)
# ax.plot(phi,theta_m,'k')
# ax.plot(phi,theta_e,'k--')
#
# # calculated values
# mub=theta_e[edge_ind]
# mut=np.max(theta_e)
# delmu=mut-mub
#
# # options
# ax.set_xlim([-40,40])
# ax.set_ylim([theta_e[edge_ind],np.max(theta_e)+2.5])
# ax.set_ylabel("Potential temperature / K")
# ax.set_xlabel('Latitude')
#
# # plot lines to show max heating and position of cells
# ax.plot([PHI0,PHI0],[theta_e[edge_ind],np.max(theta_e)+2.5],'b:')
# ax.text(PHI0-1,mub+.45*delmu, r'$\phi_0$ $\sim$ '+"{0:.2g}".format(PHI0)+'$^\circ$',ha='right',color='b')
# ax.plot([PHI1,PHI1],[theta_e[edge_ind],np.max(theta_e)+2.5],'g--')
# ax.text(PHI1+1,mub+.45*delmu, r'$\phi_1$ $\sim$ '+"{0:.2g}".format(PHI1)+'$^\circ$',ha='left',color='g')
#
# # print cell edges
# ax.text(-25,mub+.1*delmu,r'$\phi_w$ $\sim$ '+"{0:.2f}".format(PHIw)+'$^\circ$',ha='center')
# ax.text(-25,mub+.15*delmu, r'$\phi_s$ $\sim$ '+"{0:.2g}".format(PHIs)+'$^\circ$',ha='center')
#
# # print velocity scale
# ax.text(-10,mub+.1*delmu,r'$w$ $\sim$ '+"{0:.2g}".format(vert_sc)+' $ms^{-1}$',ha='center')
# ax.text(-10,mub+.15*delmu, r'$v_w$ $\sim$ '+"{0:.2g}".format(mer_w_sc)+' $ms^{-1}$',ha='center')
# ax.text(-10,mub+.2*delmu, r'$v_s$ $\sim$ '+"{0:.2g}".format(mer_s_sc)+' $ms^{-1}$',ha='center')

#plt.show()

# save the figure
#fig.savefig('off_equatorial_ACP.png',dpi=300)

# ACP says to explore boundary conditions (min max solar heating)
# and the position of the global deserts (would this change?)
# and the position of the subtropical jet (would this change?)
# try to break the model - are the values it produces realistic?
# e.g. you woudn't get a phi0 of 70 degrees up in the midlatitudes
# what happens if you change the values of the constants? - sort of
# sensitivity analysis

# perform a simple test by running the function with 10 values for heating
# position and then plotting the results
# this is the array of heating positions
heating=np.linspace(0,6,50) # 8.59 degrees = 0.15 radians

# initialize the arrays for the output, each should have 2 dimensions
# one for the delta_theta = 50 and one for the delta_theta = 98
# phi is an array of the latitude of the heating position with size (801,10)
phi=np.zeros((801,len(heating)))
PHIw=np.zeros(len(heating))
PHIs=np.zeros(len(heating))
PHI1=np.zeros(len(heating))
PHI0=np.zeros(len(heating))
theta_m1=np.zeros(len(heating))
# theta m and theta e are arrays with size (801,10)
theta_e=np.zeros((801,len(heating)))
theta_m=np.zeros((801,len(heating)))
vert_sc=np.zeros(len(heating))
mer_w_sc=np.zeros(len(heating))
mer_s_sc=np.zeros(len(heating))

# initialize a list for the boolean conv flag
conv=np.zeros(len(heating),dtype=bool)
# initialize a list for the message
message=np.zeros(len(heating),dtype='str')

# initialize a second set of arrays for the delta_theta = 98 runs
phi98=np.zeros((801,len(heating)))
PHIw98=np.zeros(len(heating))
PHIs98=np.zeros(len(heating))
PHI198=np.zeros(len(heating))
PHI098=np.zeros(len(heating))
theta_m198=np.zeros(len(heating))
theta_e98=np.zeros((801,len(heating)))
theta_m98=np.zeros((801,len(heating)))
vert_sc98=np.zeros(len(heating))
mer_w_sc98=np.zeros(len(heating))
mer_s_sc98=np.zeros(len(heating))
conv98=np.zeros(len(heating),dtype=bool)
message98=np.zeros(len(heating),dtype='str')

# loop through the heating positions
for i in range(len(heating)):
    # initialize the output
    output=held_hou.off(heating=heating[i],delta_theta=48,theta_0=300)
    # unpack the output
    phi[:,i],PHIw[i],PHIs[i],PHI1[i],PHI0[i],theta_m1[i],theta_e[:,i],theta_m[:,i],vert_sc[i],mer_w_sc[i],mer_s_sc[i],conv[i],message[i]=output


    # flag for if solution hasn't converged so print warning message
    if conv[i] != True:
        print(message[i])

# loop through the 98 heating positions
for i in range(len(heating)):
    # initialize the output
    output98=held_hou.off(heating=heating[i],delta_theta=97,theta_0=300)
    # unpack the output
    phi98[:,i],PHIw98[i],PHIs98[i],PHI198[i],PHI098[i],theta_m198[i],theta_e98[:,i],theta_m98[:,i],vert_sc98[i],mer_w_sc98[i],mer_s_sc98[i],conv98[i],message98[i]=output98

    # flag for if solution hasn't converged so print warning message
    if conv98[i] != True:
        print(message98[i])

# plot the results
fig=plt.figure(figsize=(10,5))
ax=fig.add_subplot(1,1,1)
ax.plot(heating,PHIw*-1,'b',label=r'$\phi_w$ $\delta_\theta=1/3$')
ax.plot(heating,PHIs,'r',label=r'$\phi_s$ $\delta_\theta=1/3$')
ax.plot(heating,PHI1,'k',label=r'$\phi_1$ $\delta_\theta=1/3$')
# now plot the 98 delta_theta results
ax.plot(heating,PHIw98*-1,'b--',label=r'$\phi_w$ $\delta_\theta=1/6$')
ax.plot(heating,PHIs98,'r--',label=r'$\phi_s$ $\delta_\theta=1/6$')
ax.plot(heating,PHI198,'k--',label=r'$\phi_1$ $\delta_\theta=1/6$')
# text containing the value of theta_0
ax.text(5.4,0,r'$\theta_0=$'+str(300)+'K',fontsize=14)
ax.set_xlabel('$\phi_0$ (degrees)')
ax.set_ylabel('Latitude (degrees)')
ax.legend()
plt.show()

# save the figure
fig.savefig('heating_position_lindzen_hou_fig4.png')

#%%

