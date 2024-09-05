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

def mask2taperingweights(mask, taper_halfwidth=10, downweighting_factor=1e-2):
    '''
    Creates tapering weights from a binary mask.

    Parameters
    ----------
    mask : ndarray
        Binary mask indicating the region of interest.
    taper_halfwidth : int, optional
        Half-width of the tapering region. The default is 10.
    downweighting_factor : float, optional
        Downweighting factor applied to tapering weights. The default is 1e-2.

    Returns
    -------
    ndarray
        Tapering weights for downweighting points around the masked region.
    '''

    # Create a binary array from the mask
    taper = np.clip(mask, 0, 1)

    # Convolve with a flat hat to extend the size of the mask
    flat_hat = np.ones(taper_halfwidth)
    taper = np.convolve(taper, flat_hat, mode='same')
    taper = np.clip(taper, 0, 1)

    # Convolve with a triangular kernel to taper the weights
    triangle = 1 - np.abs(np.linspace(-1, 1, taper_halfwidth + 1))
    triangle /= np.sum(triangle)
    taper = np.convolve(taper, triangle, mode='same')

    # Set weights to 1 at the edges due to edge effect of convolution
    taper[:taper_halfwidth] = 1
    taper[-taper_halfwidth:] = 1

    # Apply downweighting factor to tapering weights
    taper = 1 + taper * (1/downweighting_factor - 1)

    # Invert to obtain the final weights
    w_taper = 1 / taper

    # Uncomment the following lines for debugging or visualization
    # gplot(np.clip(mask, 0, 1))
    # gplot(w_taper)
    # pause(o)

    return w_taper
