import numpy as np
from gplot import *
from calcspec import *

def list2mask(filename, wd=4., wl=0.01, merge=True):
   """ convert a line list to a spectrum mask """
   # wd  - [km/s] width applied to line position
   # wl  - [km/s] width of transition region for linear interpolation
   # Example:
   # >>> from masktools import *
   # >>> mask = list2mask('../lib/Th_LP07.dat', merge=False)
   # >>> gplot(mask, 'w lp')
   # >>> mask = list2mask('../lib/Th_LP07.dat')
   # >>> ogplot(mask, 'w lp')
   c = 299792.458
   hw = wd / c / 2. # half width
   wdd = hw + wl / c
   linepos = np.genfromtxt(filename, dtype=None, usecols=(0))
   maskl = linepos[:, None] * np.array([[1-wdd, 1-hw, 1+hw, 1+wdd]])

   if merge:
      # merge overlapping regions
      # note those groups might contain more than one linepos
      #    ____ ____            _________
      #  _/    X    \__   =>  _/         \__
      # might be also realised with scipy.measurement.label (for one dimension)
      # http://dragly.org/2013/03/25/working-with-percolation-clusters-in-python/
      idx, = np.where(maskl[:-1,3] >= maskl[1:,0])
      last = np.ones_like(idx)  # default: next idx element idx + 1 is the last in a group

      for i in range(len(idx)-2,0,-1): # do it reverse
         if idx[i] + 1 == idx[i+1]:
            last[i] += last[i+1]
      maskl[idx,2:] = maskl[idx+last,2:]
      maskl = np.delete(maskl, idx+1, axis=0)

   maskf = np.tile([0., 1, 1, 0], len(maskl))
   mask = np.vstack((maskl.ravel(),maskf)).T
   return mask
