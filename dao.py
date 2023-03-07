"""
# delayed action oscillator code
# replicates the model of Suarez and Schopf
# Written by Cristina Perez, 2002 (c)
# Modified by Andrew Charlton-Perez, 2008
# Translated Matlab to Python, Andrew Charlton-Perez, 2015
"""


# model equation
def daodt(T,a,Td,beta):
    return(T - T**3 - a*Td + beta)

# solver
def calc_dao(alpha=0.6,N=20000,trans=400,beta=0):

 import numpy as np

# rescaling for hidden parameters
 h=0.01
 delta = trans/50.
 delta = np.int(delta/h)

# rescale beta (input as K per year)
 beta = beta*(50./365.)

# storage arrays for Temperature and rate of change
 T=np.zeros((delta+N+1))
 dt=np.zeros((delta+N+1))

# initial conditions (need to set up an oscillator)
 T[0:delta] = -2.+0.1*np.sin(2.*np.pi*h*np.arange(0,delta))
 dt[0:delta] = 0.1*2.*np.pi*np.cos(2.*np.pi*h*np.arange(0,delta))


# 4th order Runge-Kutta time stepping
 for i in range(delta,delta+N):
    Td1 = T[i-delta]
    k1 = daodt(T[i], alpha, Td1, beta)
    Td2 = T[i-int(delta/2)]
    k2 = daodt(T[i]+0.5*h*k1, alpha, Td2, beta)
    Td3 = Td2
    k3 = daodt(T[i]+0.5*h*k2, alpha, Td3, beta)
    Td4 = T[i]
    k4 = daodt(T[i]+h*k3, alpha, Td4, beta)

    T[i+1] = (h/6.)*(k1+(2.*k2)+(2.*k3)+k4)+T[i]
    dt[i] = (h/6.)*(k1+(2.*k2)+(2.*k3)+k4)


# calculate and scale time
 time = h*np.arange(0.,delta+N+1)
# scale to intrinsic coupled mode timescale
 time = (time*50.)/365.


# cut out the initialisation
 T=T[2*delta:-1]
 dt=dt[2*delta:-1]
 time=time[2*delta:-1]

 return(time,T,dt)
