from pylab import *
import numpy as np
import scipy as sc
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from StringIO import StringIO

FNAME = np.chararray((20), itemsize=30)

FNAME[0]='44008_time_series.dat'
FNAME[1]='44013_time_series.dat'
FNAME[2]='44029_time_series.dat'

for i in range(0,3):
   #print i,FNAME[i]
   X = np.loadtxt(FNAME[i])  # data in two columns
#extract
   time    = X[:,0]
   hs_sp   = X[:,1]
   hs_o    = X[:,2]
   tm01_sp = X[:,3]
   tm01_o  = X[:,4]
   tm02_sp = X[:,5]
   tm02_o  = X[:,6]
   tp_sp = X[:,7]
   tp_o  = X[:,8]
#min/max
   xmin_hs = 0.#hs_o.min()
   xmax_hs = 10.#hs_o.max() + 1
   ymin_hs = 0.#hs_o.min()
   ymax_hs = 10.#hs_o.max() + 1
   xmin_tm01 = 0.#tm01_o.min()
   xmax_tm01 = 10.#tm01_o.max() + 1
   ymin_tm01 = 0.#tm01_o.min()
   ymax_tm01 = 10.#tm01_o.max() + 1
   xmin_tm02 = 0.#tm02_o.min()
   xmax_tm02 = 10.#tm02_o.max() + 1
   ymin_tm02 = 0.#tm02_o.min()
   ymax_tm02 = 10.#tm02_o.max() + 1
   xmin_tp = 0.#tp_o.min()
   xmax_tp = 15.#tp_o.max() + 1
   ymin_tp = 0.#tp_o.min()
   ymax_tp = 15.#tp_o.max() + 1

#time
   xmin_time = time.min() - time.min()
   xmax_time = time.max() - time.min()
   time = time - time.min()
#plot
   plt.subplot(321)
   plt.axis([xmin_time,xmax_time,ymin_hs,ymax_hs])
   plt.plot(time,hs_o,'o',color='red')
   plt.plot(time,hs_sp,'*',color='blue')
   plt.title("Time series")
   plt.ylabel('Hs Buoy [m]')
   plt.subplot(322)
   polycoeffs1 = np.polyfit(hs_o,hs_sp,1)
   yfit1 = sc.polyval(polycoeffs1, hs_sp)
   def fit(x):
     return (x)
   xfit = array ( [amin(hs_sp), amax(hs_o) ] )
   plt.plot(hs_o,hs_sp,'o')
   plt.axis([xmin_hs,xmax_hs,ymin_hs,ymax_hs])
   plt.title("Scatter")
   sum_o      = hs_o.sum()
   sum_sp     = hs_sp.sum()
   mean_o     = sum_o/hs_o.size
   mean_sp    = sum_sp/hs_sp.size
   tmp_crnms  = ((hs_o-mean_o)-(hs_sp-mean_sp))**2.
   tmp_rms    = (hs_o-hs_sp)**2.
   tmp2_crnms = tmp_crnms.sum()/hs_o.size
   tmp2_rms   = tmp_rms.sum()/hs_o.size
   cnrms      = np.sqrt(tmp2_crnms)
   rms        = np.sqrt(tmp2_rms)
   sci_cnrms  = cnrms/mean_o
   sci_rms    = rms/mean_o

   cnrmsS     = cnrms.astype('S5')
   rmsS       = rms.astype('S5')
   bias       = mean_sp-mean_o
   biasS      = bias.astype('S5')
   sci_cnrms  = sci_cnrms.astype('S5')
   sci_rms    = sci_rms.astype('S5')

   text(0.5,9.2, 'CNRMS     ='+cnrmsS)
   text(0.5,8.6,'RMS       ='+rmsS)
   text(0.5,8.0, 'BIAS      ='+biasS)
   text(0.5,7.4,'SCI CNRMS ='+sci_cnrms)
   text(0.5,6.8,  'SCI RMS   ='+sci_rms)

   plt.xlabel('Hs WWM  [m]')
   plt.ylabel('Hs Buoy [m]')

   plt.subplot(323)
   plt.axis([xmin_time,xmax_time,ymin_tm02,ymax_tm02])
   plt.plot(time,tm02_o,'o',color='red')
   plt.plot(time,tm02_sp,'*',color='blue')
   plt.ylabel('Tm02 Buoy [m]')

   plt.subplot(324)
   polycoeffs1 = np.polyfit(tm02_o, tm02_sp, 1)
   yfit1 = sc.polyval(polycoeffs1, tm02_sp)
   def fit(x):
    return (x)
    xfit = array ( [amin(tm02_sp), amax(tm02_sp) ] )
   plt.plot(tm02_o,tm02_sp,'o')
   plt.axis([xmin_tm02,xmax_tm02,ymin_tm02,ymax_tm02])
   plt.xlabel('Tm02 WWM  [s]')
   plt.ylabel('Tm02 Buoy [s]')
   sum_o      = tm02_o.sum()
   sum_sp     = tm02_sp.sum()
   mean_o     = sum_o/tm02_o.size
   mean_sp    = sum_sp/tm02_sp.size
   tmp_crnms  = ((tm02_o-mean_o)-(tm02_sp-mean_sp))**2.
   tmp_rms    = (tm02_o-tm02_sp)**2.
   tmp2_crnms = tmp_crnms.sum()/tm02_o.size
   tmp2_rms   = tmp_rms.sum()/tm02_o.size
   cnrms      = np.sqrt(tmp2_crnms)
   rms        = np.sqrt(tmp2_rms)
   sci_cnrms  = cnrms/mean_o
   sci_rms    = rms/mean_o

   cnrmsS = cnrms.astype('S5')
   rmsS   =   rms.astype('S5')
   bias   = mean_sp-mean_o
   biasS  = bias.astype('S5')
   sci_cnrms = sci_cnrms.astype('S5')
   sci_rms = sci_rms.astype('S5')
   
   text(0.5,9.4, 'CNRMS     ='+cnrmsS)
   text(0.5,8.8,'RMS       ='+rmsS)
   text(0.5,8.2, 'BIAS      ='+biasS)
   text(0.5,7.6,'SCI CNRMS ='+sci_cnrms)
   text(0.5,7,  'SCI RMS   ='+sci_rms)

   plt.subplot(325)
   plt.axis([xmin_time,xmax_time,ymin_tp,ymax_tp])
   plt.plot(time,tp_o,'o',color='red')
   plt.plot(time,tp_sp,'*',color='blue')
   plt.xlabel('TIME  [d]')
   plt.ylabel('Tp Buoy [s]')

   plt.subplot(326)
   polycoeffs1 = np.polyfit(tp_o,tp_sp,1)
   yfit1 = sc.polyval(polycoeffs1, tp_sp)
   def fit(x):
     return (x)
   xfit = array ( [amin(tp_sp), amax(tp_o) ] )
   plt.plot(tp_o,tp_sp,'o')
   plt.axis([xmin_tp,xmax_tp,ymin_tp,ymax_tp])

   sum_o   = tp_o.sum()
   sum_sp  = tp_sp.sum()
   mean_o  = sum_o/tp_o.size
   mean_sp = sum_sp/tp_sp.size
   tmp_crnms  = ((tp_o-mean_o)-(tp_sp-mean_sp))**2.
   tmp_rms    = (tp_o-tp_sp)**2.
   tmp2_crnms = tmp_crnms.sum()/tp_o.size
   tmp2_rms   = tmp_rms.sum()/tp_o.size
   cnrms      = np.sqrt(tmp2_crnms)
   rms        = np.sqrt(tmp2_rms)
   sci_cnrms  = cnrms/mean_o
   sci_rms    = rms/mean_o

   cnrmsS = cnrms.astype('S5')
   rmsS   =   rms.astype('S5')
   bias   = mean_sp-mean_o
   biasS  = bias.astype('S5')
   sci_cnrms = sci_cnrms.astype('S5')
   sci_rms = sci_rms.astype('S5')

   text(0.5,13.7,'CNRMS     ='+cnrmsS)
   text(0.5,12.9,'RMS       ='+rmsS)
   text(0.5,12.1, 'BIAS      ='+biasS)
   text(0.5,11.3,'SCI CNRMS ='+sci_cnrms)
   text(0.5,10.5,  'SCI RMS   ='+sci_rms)

#cb = plt.colorbar()
#cb.set_label('log(counts)')
   plt.xlabel('Hs WWM  [m]')
   plt.ylabel('Hs Buoy [m]')

   legend(loc=2)
   #plt.show()
   fig = plt.gcf()
   fig.set_size_inches(8.,12.)
   #plt.savefig('test2png.png',dpi=100)
   #print FNAME[i]
   savefig(FNAME[i]+'.png',dpi=200)
   plt.clf()
