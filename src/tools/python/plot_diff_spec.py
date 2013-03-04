#!/usr/bin/env python
"""
Illustrate simple contour plotting, contours on an image with
a colorbar for the contours, and labelled contours.

See also contour_image.py.
"""
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.colors import LogNorm

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

sp = np.loadtxt('44013_diffspec.dat')

IT = sp[:,0]
X = sp[:,1] 
Y = sp[:,2] * 10 
e_obs = sp[:,3]
e_sim = sp[:,4]
triang = tri.Triangulation(X, Y)
Z = sp[:,0] * 0.
LOGZ = sp[:,0] * 0.
a = 0
#while a < len(X):
#  if np.abs(e_obs[a]) > 0.:
#    Z[a] = (e_obs[a]-e_sim[a])/e_obs[a]*100.
#  a = a + 1
#  print a
Z = e_obs-e_sim
N = len(Z)
max_Z = np.max(Z)
min_Z = np.min(Z)
#Z[N-1] = -max_Z
#Z[N-2] = -min_Z

maxval = 5. 
minval = -maxval 

while a < len(X):
  if Z[a] > maxval:
    Z[a] = maxval 
  if Z[a] < minval:
    Z[a] = minval
  #if Z[a] > 0.:
  #  LOGZ[a] = np.log(Z[a])
  #if Z[a] < 0.:
  #  LOGZ[a] = -np.log(-Z[a])
  #print a, Z[a], LOGZ[a]
  a = a + 1

#plt.figure()
plt.subplot(211)
plt.gca().set_aspect('equal')
plt.tricontourf(triang, Z, 199, cmap=plt.cm.RdBu)
plt.colorbar()
#plt.tricontour(triang, Z, V)
plt.title('Contour plot of Delaunay triangulation')

plt.show()

