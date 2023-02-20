# -*- coding: utf-8 -*-
"""

This is driving code for the moist case
There are fewer options than for the standard dry case, but you can change the 
parameters in the function itself (where they are hard coded)

# outputs are latitude, latitude of cell edge, latitude of edge of heating
# theta_m, theta_e, theta_estar, scale for v and scale for w

# there is only one option to change which is a scale for the width of the moist
# heating in units of the fraction of the width of the dry cell

Created on Mon Nov 16 08:59:29 2015

@author: andrewcharlton-perez
"""

import numpy as np
import matplotlib.pyplot as plt
import held_hou


# this shows you how to run the Held-Hou moist function 

# set up ten values for the width between 0.1 and 5.0
widths=np.linspace(0.1,5.0,10)


# this line is all you need to run the model
phi,PHI,PHIw,theta_m,theta_e,theta_estar,v,w=held_hou.moist(width=0.3)

# run the model held_hou.moist for each width
# and store the results in a list
# this is a list of lists
# each list is the output from the model
# the first list is the output for width=0.1
# the second list is the output for width=0.2
# etc.
# the first element of each list is the latitude
# the second element is the latitude of the edge of the cell
# the third element is the latitude of the edge of the heating
# the fourth element is theta_m
# the fifth element is theta_e
# the sixth element is theta_estar
# the seventh element is the scale for v
# the eighth element is the scale for w
# the ninth element is the convergence flag
# the tenth element is the message
# the eleventh element is the width
# the twelfth element is the heating
# the thirteenth element is the latitude of the edge of the dry cell
# the fourteenth element is the latitude of the edge of the dry heating
# the fifteenth element is the latitude of the edge of the dry heating
# the sixteenth element is the latitude of the edge of the dry heating
# the seventeenth element is the latitude of the edge of the dry heating
# the eighteenth element is the latitude of the edge of the dry heating
# the nineteenth element is the latitude of the edge of the dry heating
# the twentieth element is the latitude of the edge of the dry heating
# the twenty-first element is the latitude of the edge of the dry heating
# the twenty-second element is the latitude of the edge of the dry heating
# set up the lists
# phi=[]
# PHI=[]
# PHIw=[]
# theta_m=[]
# theta_e=[]
# theta_estar=[]
# v=[]
# w=[]
# convergence=[]
#
# # run the model for each width
# v1,v2,v3,v4,v5,v6,v7,v8,v9=held_hou.moist(width=widths[0])
# phi.append(v1)
# PHI.append(v2)
# PHIw.append(v3)
#
#
# # plot the results
# fig=plt.figure(figsize=(10,5))
# # theta values
# ax=fig.add_subplot(1,1,1)
# for i in range(len(widths)):
#     ax.plot(phi[i],theta_m[i])
# ax.set_xlim([-phi[edge_ind],phi[edge_ind]])
# ax.set_ylim([theta_e[edge_ind],np.max(theta_estar)])
# ax.set_ylabel("Potential temperature / K")
# ax.set_xlabel('Latitude')
#
#
#
# # plot the results
# fig=plt.figure(figsize=(10,5))
# # theta values
# ax=fig.add_subplot(1,1,1)
# ax.plot(phi,theta_m,'k')
# ax.plot(phi,theta_e,'k--')
# ax.plot(phi,theta_estar,'r')
#
# # show the plot
# plt.show()


# ##########
#
# This is an example of how you can plot this data
# make plot and add numbers for velocity scales etc.
# get value at edge for scaling axis and then add another 20 points
edge_ind=np.argmin(np.abs(phi-PHI))+20
heat_ind=np.argmin(np.abs(phi-PHIw))

# plotting steps
fig=plt.figure(figsize=(10,5))
# theta values
ax=fig.add_subplot(1,1,1)
ax.plot(phi,theta_m,'k')
ax.plot(phi,theta_e,'k--')
ax.plot(phi,theta_estar,'r')


# options
ax.set_xlim([-phi[edge_ind],phi[edge_ind]])
ax.set_ylim([theta_e[edge_ind],np.max(theta_estar)])
ax.set_ylabel("Potential temperature / K")
ax.set_xlabel('Latitude')

# calculated values

mub=np.max(theta_estar)
delmu=np.max(theta_estar)-theta_estar[edge_ind]

ax.text(-15,mub-.1*delmu,r'$\phi(Y)$ $\sim$ '+"{0:.2g}".format(PHI)+'$^\circ$',ha='center')
ax.text(-15,mub-.15*delmu, r'$\phi(heat)$ $\sim$ '"{0:.2g}".format(PHIw)+'$^\circ$',ha='center')
ax.text(-15,mub-.2*delmu,'v $\sim$ '+"{0:.2g}".format(v[0])+' ms$^{-1}$',ha='center')
ax.text(-15,mub-.25*delmu,'w $\sim$ '+"{0:.2g}".format(w[0])+' ms$^{-1}$',ha='center')

ax.grid()

#######

plt.show()



