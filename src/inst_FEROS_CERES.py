from __future__ import print_function

from read_spec import *
from astropy.time import Time

# Instrument parameters
name = __name__[5:]

obsloc = dict(lat=-29.2584, lon=-70.7345, elevation=2400.)   # lasilla
#obsname = "lasilla" # for barycorrpy

pat = '*[vs][ep].fits'  # _wave, _sp

iomax = 25
pmax =  4096 - 300

maskfile = 'telluric_mask_carm_short.dat'

# Instrument read functions
def scan(self, s, pfits=True):
   """
   Returns
   -------
   namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
   Example
   -------
   >>> read_carm_vis(filename)

   """
   HIERARCH = 'HIERARCH '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   self.header = hdr = hdulist[0].header
   self.drs = hdr.get('PIPELINE', 'DRS')

   if self.drs == 'CERES':
      self.instname = hdr['INST'] + "_CERES"
      self.drsberv = hdr.get('BARYCENTRIC CORRECTION (KM/S)', np.nan)
      self.drsbjd = hdr.get('MBJD', np.nan) + 2400000.5  # same as MJD!?
      self.dateobs = hdr['HIERARCH SHUTTER START DATE'] + 'T' + hdr['HIERARCH SHUTTER START UT']
      self.mjd = hdr.get('HIERARCH MJD')
      self.drift = np.nan
      self.e_drift = np.nan
      self.fileid = self.dateobs
      self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))
      self.timeid = self.fileid
      self.ccf.rvc = hdr.get('RV', np.nan)
      self.ccf.err_rvc = hdr.get('RV_E', np.nan)

      self.ra = hdr['HIERARCH RA']
      self.de = hdr['HIERARCH DEC']
      self.airmass = hdr.get('HIERARCH TARG AIRMASS START', np.nan)
      self.exptime = hdr['HIERARCH TEXP (S)']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5
      # estimate SNR
      o = 15   # reference order
      fo = hdulist[0].section[1][o]
      eo = 1/np.sqrt(hdulist[0].section[2][o])
      self.sn55 = np.nanmedian(np.abs(fo/eo))
      hdr['OBJECT'] = hdr['HIERARCH TARGET NAME']

def data(self, orders=None, pfits=True):
   hdulist = self.hdulist
   # read order data
   if self.drs == 'CERES':
      #f = hdulist[0].section[1][orders]   # with blaze (can require deg=4)
      f = hdulist[0].section[3][orders]   # deblazed
      w = hdulist[0].section[0][orders]   # / (1.00004) only need for FIES? Or, was it fixed?
      #e = 1 / np.sqrt(hdulist[0].section[2][orders])   # with blaze
      e = 1 / np.sqrt(hdulist[0].section[4][orders])   # deblazed
      w = w / (1 + self.drsberv/299_792.458)   # undo berv correction

   bpmap = np.isnan(f).astype(int)        # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[f == 0] |= flag.nan
      bpmap[e == 0] |= flag.nan

   return w, f, e, bpmap
