# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class testIceAttenuation(util.wwmTestTemplate):

  def test(self):
    # launching schismWWM
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=4, endParam=4)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')
    self.fail('ah')

    ###########################
    ### CHECKING THE OUTPUT ###
    ###########################
    ds = nc.Dataset('outputs/schout_10.nc')
    # checking the time
    tmnc = ds.variables['time']
    lasttime = nc.num2date(tmnc[-1], tmnc.units, 'standard')
    self.assertEqual(datetime(2000, 1, 11), lasttime)
    
    # checking that hs is (approximately) uniform in the last time step
    hs = ds.variables['WWM_1'][-1,:]
    self.assertTrue(np.all(hs < 3.9))
    self.assertTrue(np.all(hs > 3.6))
   
    # checking that tm01 is (approximately) uniform in the last time step
    tm = ds.variables['WWM_2'][-1,:]
    self.assertTrue(np.all(tm < 6.95))
    self.assertTrue(np.all(tm > 6.90))

    # closing the file
    ds.close()
    


    


if __name__ == '__main__':
  unittest.main()
