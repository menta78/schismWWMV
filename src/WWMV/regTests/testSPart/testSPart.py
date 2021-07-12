# author: Lorenzo Mentaschi
import os
import unittest
import netCDF4 as nc
import numpy as np
from datetime import datetime

from wwmTestUtil import wwmTestUtil as util


class testSPart(util.wwmTestTemplate):

  def test(self):
    # launching schismWWM
    launchCmd = util.getScWwmLaunchCommand(nProcesses=2)
    exitCode = os.system(launchCmd)
    if exitCode != 0:
      self.fail('schismWWM did not end correctly. Failing')
    # combining the output
    combineCmd = util.getCombineCommand(bgnParam=1, endParam=1)
    exitCode = os.system(combineCmd)
    if exitCode != 0:
      self.fail('combine_outputXX did not end correctly. Failing')

    self.fail('')

    # closing the file
    ds.close()
    


    


if __name__ == '__main__':
  unittest.main()
