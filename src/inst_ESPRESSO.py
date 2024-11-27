from read_spec import *
from read_spec import Inst
from calcspec import redshift

# https://www.eso.org/sci/facilities/paranal/instruments/espresso/ESPRESSO_User_Manual_P105_v1.0.pdf
# ftp://ftp.eso.org/pub/dfs/pipelines/instruments/espresso/espdr-reflex-tutorial-1.3.2.pdf

name = inst = __name__[5:]
drs = name[9:]   # strip ESPRESSO

# ESPRESSO     2.0 and 2.5. Version 3.1 untested. Cf. ['ESO PRO REC1 PIPE ID']
# ESPRESSO     HR_1x1 (espdr/1.3.0) [170x9111]
# EXPRESSO_HR2 HR_4x2 (binning 2x) [170x4545]
# EXPRESSO_MR  4UT medium resolution; both slices are already merged into one order [85x4545]

# Instrument parameters
pat = "*_0003.fits *.tar *ES_S2DA*.fits"   # splitted by white space

obsname = 'paranal'   # for barycorrpy
obsloc = dict(lat=-24.6268, lon=-70.4045, elevation=2648.)   # HIERARCH ESO TEL3 GEO*

oset = np.s_[6:] #if not join_slice else np.s_[6/2:]

user_mode = 'MULTIMR' if 'MR' in drs else 'SINGLEHR'   # medium (85 slices) or high (170 slices) resolution

iomax = 85 if user_mode == 'MULTIMR' else 170

biny = 2 if 'HR2' in drs else 1
pmin = 1600
pmax = 7500

if biny == 2:
   pmin = 800
   pmax = 3700

join_slice = (user_mode != 'MULTIMR') and False   # experimental and only drs<2.0.0
# The idea is to create one high S/N template per order, instead of two noisy templates.
# Runs slower, because wavelengths needs to be re-sorted.

if join_slice:
   oset = np.s_[6/2:]
   iomax /= 2
   pmin *= 2
   pmax *= 2

def scan(self, s, pfits=True, verb=False):
   if '.tar' in s:
      # ES_S2BA has blaze
      # ES_S2DA is deblazed
      # ES_PFBA is blaze function
      s = file_from_tar(s, inst=inst, fib=self.fib, pat="ES_S2DA", pfits=pfits)
   if 1:
      self.HIERARCH = HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + 'ESO '
      HIERQC = HIERINST + 'QC '
      k_tmmean = HIERINST + 'OCS EM OBJ3 TMMEAN'
      self.HIERDRS = HIERDRS = HIERINST + 'DRS '
      k_sn55 = HIERQC + ('ORDER55 SNR' if user_mode == 'MR' else 'ORDER110 SNR') # @ 573 nm
      k_berv = HIERQC + 'BERV'
      k_bjd = HIERQC + 'BJD'

      self.hdulist = hdulist = pyfits.open(s)
      hdr = self.hdulist[0].header

      self.instname = hdr['INSTRUME'] + name[8:]
      if self.instname[:8] != 'ESPRESSO':
         pause('\nWARNING: inst should be ESPRESSO, but got: '+self.inst+'\nSee option -inst for available inst.') 

      # Get unit telescope used (1, 2, 3, or 4)
      # TELESCOP= 'ESO-VLT-U1'         / ESO <TEL>
      # HIERARCH ESO OCS TEL NO =    1 / Number of active telescopes (1 - 4).
      # HIERARCH ESO OCS TEL1 ST =   T / Availability of the telescope.
      # HIERARCH ESO OCS TEL2 ST =   F / Availability of the telescope.
      # HIERARCH ESO OCS TEL3 ST =   F / Availability of the telescope.
      # HIERARCH ESO OCS TEL4 ST =   F / Availability of the telescope.
      utnum = hdr['TELESCOP'].replace('ESO-VLT-U', '')

      # self.airmass = hdr.get(HIERINST+'TEL3 AIRM START', np.nan)
      self.airmass = hdr.get(HIERINST+'TEL{} AIRM START'.format(utnum), np.nan)
      self.exptime = hdr['EXPTIME']
      self.mjd = hdr['MJD-OBS']
      self.dateobs = hdr['DATE-OBS']
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

#HIERARCH ESO TEL3 GEOELEV = 2648. / [m] Elevation above sea level               
#HIERARCH ESO TEL3 GEOLAT = -24.6268 / [deg] Tel geo latitute (+=North)          
#HIERARCH ESO TEL3 GEOLON = -70.4045 / [deg] Tel geo longitude (+=East)          
      self.obs.lon = hdr['HIERARCH ESO TEL{} GEOLON'.format(utnum)]
      self.obs.lat = hdr['HIERARCH ESO TEL{} GEOLAT'.format(utnum)]
      self.obs.elevation = hdr['HIERARCH ESO TEL{} GEOELEV'.format(utnum)]

      self.tmmean = hdr.get(k_tmmean, 0.5)
      if not (0.25 < self.tmmean < 0.75): print('WARNING:')

      self.drsbjd = hdr.get(k_bjd)
      self.drsberv = hdr.get(k_berv, np.nan)
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.drift = hdr.get(HIERQC+'DRIFT MEAN', np.nan)

      if 1:
         ccf = self.ccf
         hdrccf = hdr # pyfits.open(s.replace("_0003.fits", "_0011.fits"))[0].header
         hdrccf['PIPEFILE']
         ccf.rvc = hdrccf['HIERARCH ESO QC CCF RV']
         ccf.err_rvc = hdrccf['HIERARCH ESO QC CCF RV ERROR']
         ccf.fwhm = hdrccf['HIERARCH ESO QC CCF FWHM']
         ccf.e_fwhm = hdrccf['HIERARCH ESO QC CCF FWHM ERROR']
         ccf.contrast = hdrccf['HIERARCH ESO QC CCF CONTRAST']
         ccf.e_contrast = hdrccf['HIERARCH ESO QC CCF CONTRAST ERROR']
         ccf.mask = hdrccf['HIERARCH ESO PRO REC1 PARAM6 VALUE']   # mask_table_id
         ccf.rvguess = hdrccf['HIERARCH ESO PRO REC1 PARAM3 VALUE']   # rv_center
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

      self.timeid = ffileid = hdr['ARCFILE'][6:29]
      self.calmode = hdr.get(HIERINST+'INS3 CALSEL NAME','NOTFOUND')

      if hdr['ESO INS MODE'] != user_mode:    # SINGLEHR, MULTIMR, UHR
         self.flag |= sflag.config
         print('wrong mode:', user_mode, hdr['ESO INS MODE'])
      if hdr['ESO DET BINY'] != biny:
         self.flag |= sflag.config
         print('wrong binning:', biny, hdr['ESO DET BINY'])


      hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
      self.header = self.hdr = hdr # self.header will be set to None


def data(self, orders, pfits=True):
   hdr = self.hdr
   if not hasattr(self, 'hdulist'):
      scan(self, self.filename)

   waveext = 'WAVEDATA_VAC_BARY' if 'WAVEDATA_VAC_BARY' in self.hdulist else 'WAVEDATA_A'
   if join_slice:
       # join slices
       w = self.hdulist[waveext].data
       dim = w.shape
       dim = dim[0]/2, dim[1]*2
       w = w.reshape(dim)
       ii = np.argsort(w[orders], axis=-1)
       oo = np.arange(np.shape(w)[0])[orders,np.newaxis]
       w = w[oo,ii]
       f = self.hdulist['SCIDATA'].data.reshape(dim)[oo,ii]
       e = self.hdulist['ERRDATA'].data.reshape(dim)[oo,ii]
       bpmap = 1 * (self.hdulist['QUALDATA'].data.reshape(dim)[oo,ii] > 0)
   else:
       f = self.hdulist['SCIDATA'].section[orders]
       e = self.hdulist['ERRDATA'].section[orders]
       w = self.hdulist[waveext].section[orders]
       bpmap = 1 * (self.hdulist['QUALDATA'].section[orders] > 0)

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
      bpmap[f > 300000] |= flag.sat    # unchecked estimate for saturation level

   #w = airtovac(w)
   w = redshift(w, ve=self.drsberv, wlog=False)   # ESPRESSO wavelengths are BERV corrected. Undo the correction here.

   return w, f, e, bpmap


