from read_spec import *
from read_spec import Inst

# https://www.eso.org/sci/facilities/paranal/instruments/espresso/ESPRESSO_User_Manual_P105_v1.0.pdf
# ftp://ftp.eso.org/pub/dfs/pipelines/instruments/espresso/espdr-reflex-tutorial-1.3.2.pdf

# Instrument parameters
pat = "*_0003.fits"

name = inst = __name__[5:]
obsname = 'paranal' # for barycorrpy
obsloc = dict(lat=-24.6268, lon=-70.4045, elevation=2648.) # HIERARCH ESO TEL3 GEO*

iomax = 170
pmax = 8800

def scan(self, s, pfits=True, verb=False):
   drs = self.drs
   if '.tar' in s:
      s = file_from_tar(s, inst=inst, fib=self.fib, pfits=pfits)

   if 1:
      self.HIERARCH = HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + 'ESO '
      HIERQC = HIERINST + 'QC '
      k_tmmean = HIERINST + 'OCS EM OBJ1 TMMEAN'
      self.HIERDRS = HIERDRS = HIERINST + 'DRS '
      k_sn55 = HIERQC + 'ORDER55 SNR'
      k_berv = HIERQC + 'BERV'
      k_bjd = HIERQC + 'BJD'

      self.hdulist = hdulist = pyfits.open(s)
      hdr = self.hdulist[0].header

      self.instname = hdr['INSTRUME']
      if self.instname != 'ESPRESSO':
         pause('\nWARNING: inst should be ESPRESSO, but got: '+self.inst+'\nSee option -inst for available inst.') 

      self.airmass = hdr.get(HIERINST+'TEL3 AIRM START', np.nan)
      self.exptime = hdr['EXPTIME']
      self.mjd = hdr['MJD-OBS']
      self.dateobs = hdr['DATE-OBS']
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

#HIERARCH ESO TEL3 GEOELEV = 2648. / [m] Elevation above sea level               
#HIERARCH ESO TEL3 GEOLAT = -24.6268 / [deg] Tel geo latitute (+=North)          
#HIERARCH ESO TEL3 GEOLON = -70.4045 / [deg] Tel geo longitude (+=East)          
      self.obs.lon = hdr['HIERARCH ESO TEL3 GEOLON']
      self.obs.lat = hdr['HIERARCH ESO TEL3 GEOLAT']
      self.obs.elevation = hdr['HIERARCH ESO TEL3 GEOELEV']

      self.tmmean = hdr.get(k_tmmean, 0.5)

      self.drsbjd = hdr.get(k_bjd)
      self.drsberv = hdr.get(k_berv, np.nan)
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.drift = hdr.get(HIERQC+'DRIFT MEAN', np.nan)

      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

      self.timeid = ffileid = hdr['ARCFILE'][6:29]
      self.calmode = hdr.get(HIERINST+'INS3 CALSEL NAME','NOTFOUND')

      if hdr[HIERINST+'INS MODE'] == 'EGGS':
         self.flag |= sflag.eggs

      hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
      self.header = self.hdr = hdr # self.header will be set to None


def data(self, orders, pfits=True):
   hdr = self.hdr
   if not hasattr(self, 'hdulist'):
      scan(self, self.filename)
   f = self.hdulist['SCIDATA'].section[orders]
   e = self.hdulist['ERRDATA'].section[orders]
   w = self.hdulist['WAVEDATA_A'].section[orders]
   bpmap = 1 * (self.hdulist['QUALDATA'].section[orders] > 0)

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
      bpmap[f > 300000] |= flag.sat    # unchecked estimate for saturation level

      w = airtovac(w)
      return w, f, e, bpmap


