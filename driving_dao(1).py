# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 07:42:44 2016


This code is used to drive the dao model (you need to make sure you have both 
files present)

There are three options here to change (you could explore the space by looping)
over this function for example
alpha - efficieny of wave damping (between 0 and 1)
trans - transit time for wave to cross Pacific in days 
N - length of integration (in units of half a day)
if you run with no options you get alpha=0.6, trans=400, N=20000


There is one additional parameters used in the MT4YA problem following 
Bouttle et al. (2006; doi: 10.1119/1.2358155)
beta - a crude way of representing warming SST as a constant increase over time
with units of degrees per year

This is set to zero as standard 


@author: andrewcharlton-perez
"""

import matplotlib.pylab as plt
import dao
from matplotlib import mlab


# drive calculation and plot output
time,T,dt=dao.calc_dao()


# example plot showing T, a phase portrait and a power spectrum
plt.figure(figsize=(7,10))
plt.subplot(311)
plt.plot(time,T, 'b-')
plt.xlabel('Time / years')
plt.ylabel('T anomaly / K')
plt.title('Temperature')
plt.grid()

plt.subplot(312)
plt.plot(T,dt, 'm+')
plt.xlabel('T anomaly / K')
plt.ylabel('dT/dt')
plt.title('Phase portrait')
plt.tight_layout()


plt.subplot(313)


# computer the power spectrum using Welch's method
time_step=(time[1]-time[0])
power, freq = mlab.psd(T, NFFT=10001, Fs=1. / time_step,
                      window=mlab.window_hanning, noverlap=2000)

period=1./freq


plt.plot(period,power, 'g')
plt.xlabel('Frequency / years')
plt.ylabel('Power')
plt.xlim([0,10])
plt.title('Power spectrum')
plt.tight_layout()

plt.show()

