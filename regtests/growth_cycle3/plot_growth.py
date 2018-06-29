#!/usr/bin/env python
"""
Compute the cross spectral density of two signals
"""
import numpy as np
import matplotlib.pyplot as plt
from StringIO import StringIO

obs   = np.loadtxt('120s_phillim_dia.site',skiprows=1)
test1 = np.loadtxt('120s_phillim_dia.site',skiprows=1)
#test2 = np.loadtxt('strang600nolim.site',skiprows=1)   
#test3 = np.loadtxt('nostrang600nolim.site',skiprows=1)
test4 = np.loadtxt('120s_phillim_dia.site',skiprows=1)

nobs = len(obs[:,1])
print nobs
xobs = obs[:,0] # time 
err_test = xobs
err_10_2s = xobs
err_test = 0
err_10_2s = 0
dt = 600./3600.
i = 0 
xobs[i] = 0
obs1 = obs[:,1]
t1 = test1[:,1]
#t2 = test2[:,1]
#t3 = test3[:,1]
t4 = test4[:,1]

err10_2s = (t1 - obs1)/t1 * 100.
err_test = (t4 - obs1)/t4 * 100.

while i < nobs-1:
  i = i + 1
  xobs[i] = xobs[i-1] + dt 
  #print i, xobs[i], obs1[i]

xmin = xobs.min()
xmax = xobs.max()
ymin = 0#err.min() 
ymax = 10#err.max() 

#fig = plt.subplot(121) 
plt.plot(xobs , obs1, color='green', linestyle='solid', linewidth=2, label='120s FBI Philim')
plt.plot(xobs , t1, color='blue', linestyle='solid', linewidth=2, label='nolim 10. sec')
#plt.plot(xobs , t2, color='red', linestyle='solid', linewidth=2, label='strang 600s nolim')
#plt.plot(xobs , t3, color='black', linestyle='solid', linewidth=2, label='no strang 600s nolim')
#plt.plot(xobs , t4, color='cyan', linestyle='solid', linewidth=2, label='test')
err10_2s[0] = 0.
err_test[0] = 0.

plt.plot(xobs , abs(err10_2s), color='red', linestyle='dashed', linewidth=1, label='err_10_2')
plt.plot(xobs , abs(err_test), color='black', linestyle='dashed', linewidth=1, label='err_test')
plt.axis([xmin, xmax, ymin, ymax])
plt.xlabel('Frequency [Hz]')
plt.ylabel('Energy Density [m2/Hz]')
plt.grid(True)
plt.legend(loc=0,ncol=1)
plt.show()


