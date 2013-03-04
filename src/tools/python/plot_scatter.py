"""
hexbin is an axes method or pyplot function that is essentially
a pcolor of a 2-D histogram with hexagonal cells.  It can be
much more informative than a scatter plot; in the first subplot
below, try substituting 'scatter' for 'hexbin'.
"""
from pylab import *
import numpy as np
import scipy as sc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from StringIO import StringIO

#n = 100000
#x = np.random.standard_normal(n)
#y = 2.0 + 3.0 * x + 4.0 * np.random.standard_normal(n)
X = np.loadtxt('scatter.dat')  # data in two columns
time    = X[:,0]
hs_sp   = X[:,1]
hs_o    = X[:,2]
tm01_sp = X[:,3]
tm01_o  = X[:,4]
tm02_sp = X[:,5]
tm02_o  = X[:,6]

ymin_hs = hs_sp.min()
ymax_hs = hs_sp.max()
xmin_hs = hs_o.min()
xmax_hs = hs_o.max()

ymin_tm01 = tm01_sp.min()
ymax_tm01 = tm01_sp.max()
xmin_tm01 = tm01_o.min()
xmax_tm01 = tm01_o.max()

ymin_tm02 = tm02_sp.min()
ymax_tm02 = tm02_sp.max()
xmin_tm02 = tm02_o.min()
xmax_tm02 = tm02_o.max()

plt.subplot(231)
polycoeffs1 = np.polyfit(hs_sp, hs_o, 1)
print polycoeffs1
yfit1 = sc.polyval(polycoeffs1, hs_sp)
def fit(x):
    return (x)
xfit = array ( [amin(hs_sp), amax(hs_o) ] )
plt.hexbin(hs_o,hs_sp,bins='log',cmap=cm.jet,facecolor='white')
#plt.plot(hs_sp,hs_o,color='black',markerfacecolor='red')
#plt.plot(hs_sp,hs_o,'o')
plt.plot(xfit,fit(xfit),color='red',linestyle='dashed',markerfacecolor='red',linewidth=2,label='Best Fit')
plt.plot(hs_sp,yfit1,color='black',linestyle='dashed',linewidth=2,label='Fit to Data')
plt.axis([xmin_hs,xmax_hs,ymin_hs,ymax_hs])
plt.title("Scatter of u10")

sum_o   = hs_o.sum()
sum_sp  = hs_sp.sum()
mean_o  = sum_o/hs_o.size
mean_sp = sum_sp/hs_sp.size

tmp_crnms  = ((hs_o-mean_o)-(hs_sp-mean_sp))**2.
tmp_rms    = (hs_o-hs_sp)**2.
tmp2_crnms = tmp_crnms.sum()/hs_o.size
tmp2_rms   = tmp_rms.sum()/hs_o.size
cnrms      = np.sqrt(tmp2_crnms)
rms        = np.sqrt(tmp2_rms)
sci_cnrms  = cnrms/mean_o
sci_rms    = rms/mean_o
print cnrms
print sci_cnrms
print sci_rms
cb = plt.colorbar()
cb.set_label('log(counts)')
plt.xlabel('Hs WWM  [m]')
plt.ylabel('Hs Buoy [m]')

#plt.subplot(231)
#plt.hexbin(hs_sp,hs_o,bins='log',cmap=cm.jet,facecolor='blue')
#cb = plt.colorbar()
#cb.set_label('log(counts)')

plt.subplot(232)
polycoeffs1 = np.polyfit(tm01_sp, tm01_o, 1)
print polycoeffs1
yfit1 = sc.polyval(polycoeffs1, tm01_sp)
def fit(x):
    return (x)
xfit = array ( [amin(tm01_sp), amax(tm01_sp) ] )
plt.hexbin(tm01_o,tm01_sp,bins='log',cmap=cm.jet,facecolor='white')
#plt.plot(hs_sp,hs_o,color='black',markerfacecolor='red')
#plt.plot(tm01_sp,tm01_o,'o')
plt.plot(xfit,fit(xfit),color='red',linestyle='dashed',markerfacecolor='red',linewidth=2,label='Best Fit')
plt.plot(tm01_sp,yfit1,color='black',linestyle='dashed',linewidth=2,label='Fit to Data')
plt.axis([xmin_tm01,xmax_tm01,ymin_tm01,ymax_tm01])
plt.title("Scatter of u10")
cb = plt.colorbar()
cb.set_label('log(counts)')
plt.xlabel('Tm02 WWM  [s]')
plt.ylabel('Tm02 Buoy [s]')
diff = tm01_o - tm02_sp
diff2 = diff**2.
sum_o = tm01_o.sum()
sum_sp = tm01_sp.sum()
mean_o = sum_o/tm01_o.size
mean_sp = sum_sp/tm01_sp.size
tmp_crnms = ((tm01_o-mean_o)-(tm01_sp-mean_sp))**2.
tmp2_crnms = tmp_crnms.sum()*1./tm01_o.size
cnrms = np.sqrt(tmp2_crnms)
print cnrms
sci = cnrms/mean_o
print sci

plt.subplot(233)
polycoeffs1 = np.polyfit(tm02_sp, tm02_o, 1)
print polycoeffs1
yfit1 = sc.polyval(polycoeffs1, tm02_sp)
def fit(x):
    return (x)
xfit = array ( [amin(tm02_sp), amax(tm02_sp) ] )
plt.hexbin(tm02_o,tm02_sp,bins='log',cmap=cm.jet,facecolor='white')
#plt.plot(hs_sp,hs_o,color='black',markerfacecolor='red')
#plt.plot(tm02_sp,tm02_o,'o')
plt.plot(xfit,fit(xfit),color='red',linestyle='dashed',markerfacecolor='red',linewidth=2,label='Best Fit')
plt.plot(tm02_sp,yfit1,color='black',linestyle='dashed',linewidth=2,label='Fit to Data')
plt.axis([xmin_tm02,xmax_tm02,ymin_tm02,ymax_tm02])
plt.title("Scatter of u10")
cb = plt.colorbar()
cb.set_label('log(counts)')
plt.xlabel('Tm02 WWM  [s]')
plt.ylabel('Tm02 Buoy [s]')
diff = tm02_o - tm02_sp
diff2 = diff**2.
sum_o = tm02_o.sum()
sum_sp = tm02_sp.sum()
mean_o = sum_o/tm02_o.size
mean_sp = sum_sp/tm02_sp.size
tmp_crnms = ((tm02_o-mean_o)-(tm02_sp-mean_sp))**2.
tmp2_crnms = tmp_crnms.sum()*1./tm02_o.size
cnrms = np.sqrt(tmp2_crnms)
print cnrms
sci = cnrms/mean_o
print sci

sum_o = tm02_o.sum()
sum_sp = tm02_o.sum()
mean_o = sum_o/tm02_o.size
mean_sp = sum_sp/tm02_o.size
tmp_crnms = ((tm02_o-mean_o)-(tm02_o-mean_o))**2.
tmp2_crnms = tmp_crnms.sum()*1./tm02_o.size
cnrms = np.sqrt(tmp2_crnms)
print cnrms
sci = cnrms/mean_o
print sci


legend(loc=2)
plt.show()
