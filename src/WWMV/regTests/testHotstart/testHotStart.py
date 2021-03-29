# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class test2DExplicitPropagation(util.wwmTestTemplate):

  def test(self):
    import pdb; pdb.set_trace()
    # first run: cold start, writing output after 12 h and after 24 h
    os.system('ln -sf param_cold.nml param.nml; ln -sf wwminput_cold.nml wwminput.nml')
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=2, endParam=2)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')
    # combining the hotstart file
    combineCmd = util.getCombineHotstartCommand(iTimeStep=108)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_hotstartXX did not end correctly. Failing')

    # moving the output of the test run
    os.system('cp -r outputs outputs_cold')


    ### CHECKING THE HOTSTART FILE #



    # second run: running from hour 12 with restart.
    os.system('rm fort.* *.bin *.dat *.site wwmcheck.nml')
    os.system('rm outputs/schout_*.nc')
    os.system('ln -sf param_hot.nml param.nml; ln -sf wwminput_hot.nml wwminput.nml')
    os.system('ln -sf outputs_cold/hotstart_it=108.nc hotstart.nc')
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=2, endParam=2)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')
    

    #######################################################
    ## CHECKING THE OUTPUTS of the 2 runs are identical ###
    ##   only Hs should suffice  ##########################
    #######################################################
    ds0 = nc.Dataset('outputs_cold/schout_2.nc')
    hsNoHot = ds0.variables['WWM_1'][-1,:]
    tmNoHot = ds0.variables['time'][:]
    ds0.close()

    ds1 = nc.Dataset('outputs/schout_2.nc')
    hsHot = ds1.variables['WWM_1'][-1,:]
    tmHot = ds1.variables['time'][:]
    ds1.close()

    np.testing.assert_almost_equal(hsNoHot, hsHot)
    np.testing.assert_almost_equal(tmNoHot, tmHot)
    


    


if __name__ == '__main__':
  unittest.main()
