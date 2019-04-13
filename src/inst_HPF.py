from read_spec import *

# Instrument parameters
name = 'HPF'
obsname = "hpf" # for barycorrpy
       # wikipedia
obsloc = dict(lat=30.681444,    # 30d40'53" N
              lon=-104.014722,  # 104d00'53" W
              elevation=2026)
pat = '*.fits'
pmax = 2048 - 300
iomax = 28

maskfile = 'telluric_mask_carm_short.dat'


# Instrument read functions
def scan(self, s, orders=None, pfits=True, verb=True):
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
   if 1:
      self.header = hdr = hdulist[0].header
      self.instname = hdr['INSTRUME']
      self.drsberv = hdr.get('BERV', np.nan)
      self.drsbjd = hdr.get('BJD', np.nan) + 2400000
      self.dateobs = hdr['DATE-OBS']
      # dateobs is used by MH, but date-obs seems more reliable from FILENAME
      # CARACAL computes mjd-obs also from FILENAME
      self.mjd = hdr.get('JD_FW18') - 2400000.5
      # for HPF spectra the drift is already included in the wavelength solution
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', hdr.get(HIERARCH+'CARACAL DRIFT RV', np.nan))
      self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT RVERR', np.nan))
      self.sn55 = hdr.get('SNR 36', 50)
      self.fileid = hdr.get('DATE-OBS', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.timeid = self.fileid
      self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))

      self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
      self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)

      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['ITIME']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5

def data(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      f = hdulist['Sci Flux'].section[orders]
      w =  hdulist['Sci Wavl'].section[orders]
      e = (hdulist['Sci Variance'].section[orders])**-0.5
      bpmap = np.isnan(f).astype(int)            # flag 1 for nan

      with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

      return w, f, e, bpmap


