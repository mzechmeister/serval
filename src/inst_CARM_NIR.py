from __future__ import print_function

from read_spec import *
from read_spec import Inst
# Instrument parameters

name = 'CARM_NIR'
obsname = 'ca'
obsloc = dict(lat=37.2236, lon= -2.5463, elevation=2168.)

iomax = 28
iomax *= 2 # reshaping
pmax = 1800

oset = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 28, 29, 31, 46, 48, 50, 52]
coset = sorted(set(oset) | {0, 2, 12, 13, 16, 22, 23, 24, 25, 26, 27, 30, 32, 33, 34, 45, 47, 49, 51, 53, 54, 55})

# these are the old setting pre atm cal
oset_atmspec = sorted(set(range(iomax)) - {36, 37, 38, 39, 40, 41})   # these are the orders that are not used for telluric modeling
coset_atmspec = oset_atmspec   # these are the orders that are not used for telluric modeling

default_fib = 'A'  # default fiber for science spectra

maskfile = 'telluric_mask_nir4.dat'
#atmspec = 'atm_carm_nir.fits'
#atmspec_mask = 'telluric_mask_CARM_NIR_0.25_limit.dat'   # the 0.25 in the filename shows the transmission limit for the telluric lines
                                                         # lines that cannot be corrected; needed for CARM NIR. Other limits are also provided under /lib.
                                                         # In first tests, the 0.25 limit seems to be a good compromise between masking too many lines and
                                                         # not masking enough lines that are not corrected properly and yield the best results for Barnard's Star.
skyfile = 'sky_carm_nir'

pat = '*-nir_%(fib)s.fits *-nir_%(fib)s-tac.fits'   # => nir_A.fits, nir_B.fits

atm_cal_order = [30, 32, 33]   # these are the best orders; they are fitted simultaneously, o30 is mainly for O2 and o32 and o33 are for H2O.

def scan(self, s, pfits=True):
   """
   SYNTAX: read_carm_nir(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   HIERARCH = 'HIERARCH '
   self.hdulist = hdulist = pyfits.open(s, memmap=False) # slow 30 ms, memmap see https://github.com/mzechmeister/serval/issues/32
   hdr = hdulist[0].header
   self.header = hdr
   self.instname = hdr['INSTRUME'][0:4]+'_NIR'
   #data,hdr = fitsio.read(s,header=True)

   self.drsberv = hdr.get('HIERARCH CARACAL BERV', np.nan)
   # BJD for stars, JD for for calibration products
   self.drsbjd = hdr.get('HIERARCH CARACAL BJD', hdr.get('MJD-OBS',0.0))
   self.drsbjd = hdr.get('HIERARCH CARACAL BJD', np.nan) + 2400000
   self.dateobs = hdr['DATE-OBS']
   self.dateobs = hdr['FILENAME'].replace("h",":").replace("m",":")
   self.dateobs = self.dateobs[4:8]+"-"+self.dateobs[8:10]+"-"+self.dateobs[10:21]
   self.mjd = hdr.get('HIERARCH CARACAL MJD-OBS')
   if not self.mjd:
      import warnings
      warnings.warn("Warning: keyword HIERARCH CARACAL MJD-OBS not found! This was implemented in CARACAL v2.00."+
                    "Please use lastest products.")
   if isinstance(self.drsbjd, str): self.drsbjd = 0.0   # workaround for MJD-OBS bug @2016-Jan
   self.fox = HIERARCH+'CARACAL FOX XWD' in hdr
   # check LED for NIR This order can be affect
   r = hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan) /  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan)
   if r > 1.5:
      self.flag |= sflag.led
      print(r, hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan), hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan), "# This spectrum could be affected by LED; or back background")
   self.sn55 = min(hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan),  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan))
   if self.dateobs[:10] in ('2016-01-13', '2016-01-14', '2016-01-22', '2016-01-23'):
      self.sn55 = min(self.sn55, 10.) # bug in NIR fits file
   self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
   #self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan)
   #self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan)
   self.drift = hdr.get(HIERARCH+'CARACAL SERVAL FP RV', hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan))
   self.e_drift = hdr.get(HIERARCH+'CARACAL SERVAL FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan))
   if np.isnan(self.drift):
       self.drift = hdr.get(HIERARCH+'CARACAL GUESS FP RV', np.nan)

   self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
   self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)

   #fileid = str(hdr.ascardlist()['MJD-OBS'])     # read the comment
   self.fileid = hdr.get('FILENAME', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
   self.calmode = hdr.get(HIERARCH+'CARACAL FIB','')
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]
   self.timeid = self.fileid
   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr.get('AIRMASS', 0.0)
   self.exptime = hdr['EXPTIME']
   self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
   if self.exptime: self.tmmean /= self.exptime   # normalise
   if self.tmmean==0: self.tmmean = 0.5


def data(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      #f = hdulist['SPEC'].data
      # "data" atribute seems to open again the fits file. For large data set (GJ273) this lead to "error: too many files open". So use "section"
      # reshape orders to half-orders; bad hack to stick with section (should we avoid data.reshape?)
      hdulist['SPEC']._axes = hdulist['WAVE']._axes = hdulist['SIG']._axes = [2040, 56]
      f = 1.*hdulist['SPEC'].section[orders]
      w = hdulist['WAVE'].section[orders]
      e = hdulist['SIG'].section[orders]
      hdulist['SPEC']._axes = hdulist['WAVE']._axes = hdulist['SIG']._axes = [4080, 28]

      bpmap = np.isnan(f).astype(int)   # flag 1 for nan
      if self.fox:
         # scale spectrum
         e = e * 100000.
         f = f * 100000.
      else:
         # fib B, linear extraction
         e = 0*f + np.sqrt(np.nanmedian(f)) # unweighted maybe for HCL
         bpmap[f>300000] |= flag.sat
         #bpmap[:,1:] |= flag.sat * (f>300000)[:,:-1]
         #bpmap[:,:-1] |= flag.sat * (f>300000)[:,1:]
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      return w, f, e, bpmap

