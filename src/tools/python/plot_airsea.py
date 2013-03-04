"""
hexbin is an axes method or pyplot function that is essentially
a pcolor of a 2-D histogram with hexagonal cells.  It can be
much more informative than a scatter plot; in the first subplot
below, try substituting 'scatter' for 'hexbin'.
"""

import numpy as np
import scipy as sc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from StringIO import StringIO

#n = 100000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
#X = np.loadtxt('airsea1.dat')  # data in two columns
X = np.loadtxt('airsea.dat')  # data in two columns
windx   = X[:,0]
windy   = X[:,1]
wind    = X[:,2]
tauw    = X[:,3]
tauhf   = X[:,4]
taut    = X[:,5]
z0      = X[:,6]
ustar   = X[:,7]
alpha   = X[:,8]
cd      = X[:,9]
u1028   = X[:,2]/28.
charnock = 9.81*X[:,6]/(X[:,7]*X[:,7])


plt.subplot(231)
xmin = wind.min()
xmax = wind.max()
ymin = ustar.min()
ymax = ustar.max()
plt.hexbin(wind,ustar,bins='log', cmap=cm.jet)
#plt.plot(wind,ustar,'o',color='red')
#plt.plot(wind,u1028,'x',color='black')
plt.axis([xmin, xmax, ymin, ymax])
#cb = plt.colorbar()
#cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('Friction Velocity [m/s]')

plt.subplot(232)
xmin = wind.min()
xmax = wind.max()
ymin = z0.min()
ymax = z0.max()
plt.hexbin(wind,z0,bins='log', cmap=cm.jet)
#plt.plot(wind,z0,'o',color='red')
plt.axis([xmin, xmax, ymin, ymax])
#cb = plt.colorbar()
#cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('Roughness length [m]')

plt.subplot(233)
xmin = wind.min()
xmax = wind.max()
ymin = cd.min()
ymax = cd.max()
plt.hexbin(wind,cd,bins='log', cmap=cm.jet)
#plt.plot(wind,cd,'o',color='red')
plt.axis([xmin, xmax, ymin, ymax])
cb = plt.colorbar()
cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('Drag coefficent [-]')

plt.subplot(234)
xmin = wind.min()
xmax = wind.max()
ymin = alpha.min()
ymax = alpha.max()
plt.hexbin(wind,alpha,bins='log', cmap=cm.jet)
#plt.plot(wind,alpha,'o',color='red')
plt.axis([xmin, xmax, ymin, ymax])
#cb = plt.colorbar()
#cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('Charnock Coefficent [-]')

plt.subplot(235)
xmin = wind.min()
xmax = wind.max()
ymin = tauhf.min()
ymax = tauhf.max()
plt.hexbin(wind,tauhf,bins='log', cmap=cm.jet)
#plt.plot(wind,alpha,'o',color='red')
plt.axis([xmin, xmax, ymin, ymax])
#cb = plt.colorbar()
#cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('High Freq. Stress [-]')

plt.subplot(236)
xmin = wind.min()
xmax = wind.max()
ymin = tauw.min()
ymax = tauw.max()
plt.hexbin(wind,tauw,bins='log', cmap=cm.jet)
#plt.plot(wind,alpha,'o',color='red')
plt.axis([xmin, xmax, ymin, ymax])
cb = plt.colorbar()
cb.set_label('log(counts)')
plt.xlabel('Wind Velocity [m/s]')
plt.ylabel('Prognostic Stress [-]')


plt.show()
