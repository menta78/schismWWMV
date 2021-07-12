import numpy as np
from scipy.io import FortranFile
from datetime import datetime

import xarray as xr
from matplotlib.tri import Triangulation, LinearTriInterpolator

from alphaBetaLab import abTriangularMesh as trmsh


dateToLoad = datetime(2000, 1, 2)
timeIndex = -1
outNc = 'outputs/schout_0000_1.nc'
pt = [6., 34.]



def getSpecAtDate(dateToLoad):
  dateToLoadStr = dateToLoad.strftime('%Y%m%d.%H%M%S')

  rkind = 'float64'
  
  fl = FortranFile('P-1.sp2d', 'r')
  numsig = fl.read_ints()[0]
  numdir = fl.read_ints()[0]
  
  spsig = fl.read_reals(dtype=rkind)
  tm = 2*np.pi/spsig
  spdir = fl.read_reals(dtype=rkind)
  
  isum = fl.read_ints()
  
  datestr = ''
  while datestr != dateToLoadStr:
    dateAsByte = fl.read_record(dtype='byte')
    datestr = dateAsByte.tostring().decode('ascii')
    
    dep = fl.read_reals(dtype=rkind)
    curtxy = fl.read_reals(dtype=rkind)
    spec = np.reshape(fl.read_reals(dtype=rkind), [numdir, numsig])
    acloc = np.reshape(fl.read_reals(dtype=rkind), [numdir, numsig])

  fl.close()

  return spsig, spdir, spec



def getPeaksAtPoint(dateToLoad):
  ds = xr.open_dataset(outNc)
  p1hs = ds['WWM_P1HS'][timeIndex,:].values
  p2hs = ds['WWM_P2HS'][timeIndex,:].values
  p1tm = ds['WWM_P1TM01'][timeIndex,:].values
  p2tm = ds['WWM_P2TM01'][timeIndex,:].values
  p1dirm = ds['WWM_P1DM'][timeIndex,:].values
  p2dirm = ds['WWM_P2DM'][timeIndex,:].values
  ds.close()
  
  msh = trmsh.loadFromGr3File('hgrid.gr3')
  ndid = list(msh.nodes.keys())
  ndid.sort()
  xs = np.array([msh.nodes[k][0] for k in ndid])
  ys = np.array([msh.nodes[k][1] for k in ndid])

  triObj = Triangulation(xs,ys)
  intpltr = LinearTriInterpolator(triObj, p1hs)
  p1hsPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p1tm)
  p1tmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p1dirm)
  p1dirmPt = intpltr(pt[0], pt[1])

  intpltr = LinearTriInterpolator(triObj, p2hs)
  p2hsPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p2tm)
  p2tmPt = intpltr(pt[0], pt[1])
  intpltr = LinearTriInterpolator(triObj, p2dirm)
  p2dirmPt = intpltr(pt[0], pt[1])
  
  spsig, spdir, spec = getSpecAtDate(dateToLoad)
  tts = (spsig/2/np.pi)**-1
  drs = (np.pi/2-spdir)/np.pi*180
  
  plt.pcolor(tts, drs, spec)
  plt.scatter([p1tmPt, p2tmPt], [p1dirmPt-180, p2dirmPt-180], s=[p1hsPt, p2hsPt])

  plt.show()


