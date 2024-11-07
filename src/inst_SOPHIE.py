from read_spec import *
from read_spec import Inst
# Instrument parameters

name = inst = __name__[5:]
obsname = 'ohp' # for barycorrpy
obsloc = dict(lat=43.9308, lon=5.7133, elevation=650)
# 43.9308, 5.7133, 650 m
# HIERARCH OHP TEL ALT       =                 80.54 /                            
# HIERARCH OHP TEL GEOELEV   =                650.00 /                            
# HIERARCH OHP TEL GEOLAT    =             +43.92944 /                            
# HIERARCH OHP TEL GEOLON    =               +5.7125 /                            
iomax = 39
pat = '*_e2ds.fits'

#maskfile = 'telluric_mask_carm_short.dat'


def scan_ccf(self):
   ccffiles = glob.glob(self.filename[:-10]+'*_ccf.fits')  # replace e2ds with *_ccf
   ccf = self.ccf
   if ccffiles:
      ccf.filename = ccffiles[0]
      hdr = pyfits.getheader(ccf.filename)
      HIERDRS = 'HIERARCH OHP DRS '
      ccf.rvc = hdr.get(HIERDRS+'CCF RV', np.nan)   # [km/s], not clear whether this value is drift corrected?
      ccf.contrast = hdr.get(HIERDRS+'CCF CONTRAST', np.nan)
      ccf.fwhm = hdr.get(HIERDRS+'CCF FWHM', np.nan)
      ccf.mask = hdr.get(HIERDRS+'CCF MASK', np.nan)
      ccf.err_rvc = hdr.get(HIERDRS+'CCF ERR', np.nan)  # [km/s]
      ccf.bis = hdr.get(HIERDRS+'CCF SPAN', np.nan)   # [km/s]
      self.drsberv = hdr.get(HIERDRS+'BERV', np.nan)  # [km/s]
      self.drsbjd = hdr.get(HIERDRS+'BJD', np.nan)

def scan(self, s, pfits=True, verb=False):
   """
   SYNTAX: read_harps(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   http://atlas.obs-hp.fr/sophie/intro.html
   """
   drs = self.drs
   if isinstance(s, str) and '.gz' in s:
      # if s is isinstance(s,tarfile.ExFileObject) then s.position will change !? resulting in:
      # *** IOError: Empty or corrupt FITS file
      pfits = True

   if 1:
      HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + 'OHP '
      k_tmmean = HIERINST + 'INS PM TMEAN'
      # In old HARPN the keyword is different and the value absolute
      #k_tmmean = {'HARPS': HIERINST + 'INS DET1 TMMEAN', 'HARPN': HIERINST + 'EXP1 TMMEAN'}[inst]
      if drs:
         self.HIERDRS = HIERDRS = HIERINST + 'DRS '
         k_sn55 = HIERDRS + 'CAL EXT SN31'
         k_berv = HIERDRS + 'BERV'
         k_bjd = HIERDRS + 'BJD'

      if pfits is True:  # header with pyfits
         self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
         hdr = self.hdulist[0].header # pyfits.getheader(s)
      elif pfits==2:     # a faster version
         args = ('INSTRUME', 'OBJECT', 'MJD-OBS', 'DATE-OBS', 'OBS TARG NAME', 'EXPTIME',
                 'MJD-OBS', 'FILENAME', 'RA', 'DEC', k_tmmean, HIERINST+'DPR TYPE',
                 HIERINST+'DPR TECH', HIERINST+'INS MODE', HIERINST+'OBS TARG NAME')
         args += (k_bjd, k_berv, k_sn55)
         if drs:
            args += (HIERDRS+'BLAZE FILE', HIERDRS+'DRIFT RV USED',
                     HIERDRS+'CAL TH DEG LL', HIERDRS+'CAL LOC NBO',
                     HIERDRS+'CAL TH COEFF LL')

         hdr = imhead(s, *args)
         self.hdu = getext(s)
      else:
         #hdr = fitsio.read_header(s) no faster?
         self.f, hdrio = fitsio.read(s, header=True)
         hdr = dict((key, val.strip() if type(val) is str else val) for key,val in dict(hdrio).iteritems())
         HIERARCH = ''

      #self.drs = 'DRS CAL LOC NBO' in "".join(hdr.keys())  # check DRS or FOX
      self.instname = hdr['INSTRUME'] #if self.drs else 'HARPS'
      if self.instname != name:
         pause('\nWARNING: inst should be HARPS or HARPN, but got: '+self.inst+'\nSee option -inst for available inst.') 
      self.HIERARCH = HIERARCH

      version = float(hdr[HIERINST+'DRS VERSION'])
      self.airmass = hdr.get('AIRMASS', np.nan)
      #self.exptime = hdr['EXPTIME']
      self.exptime = hdr[HIERINST+'CCD UIT']
      self.mjd = hdr[HIERINST+'OBS MJD']
      if version < 0.50:
          self.mjd -= 2400000.5  # now it is really MJD

      self.hdulist[0].verify('silentfix')
      if 'rounded' in hdr.comments[HIERINST+'OBS DATE START']:
          print('WARNING: proprietary data, dates are rounded.')
          self.dateobs = hdr[HIERINST+'OBS DATE START'] + 'T00:00:00.000'
      else:
          self.dateobs = hdr[HIERINST+'OBS DATE START']
          if self.dateobs[-4] == ":":  # old data (DRS <v0.50, ~<2007) format ms with colon: T00:00:00:000
              self.dateobs = self.dateobs[:-4] + "." + self.dateobs[-3:]

      self.ra = hdr.get(HIERINST+'TEL ALPHA', hdr.get('HHIERARCH OHP TEL ALPHA')) # "HH" typo e.g. '2011-10-03T03:46:09.243' in v0.50
      self.de = hdr.get(HIERINST+'TEL DELTA', hdr.get('HHIERARCH OHP TEL DELTA'))
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

      self.obs.lon = obsloc['lon']
      self.obs.lat = obsloc['lat']

      if k_tmmean not in hdr:
         warnings.warn('Warning: old HARPN data? Setting tmmean to 0.5!')
      self.tmmean = hdr.get(k_tmmean, 0.5)

      self.drsbjd = hdr.get(k_bjd)
      if self.drsbjd is None:
         self.drsbjd = hdr.get('MJD-OBS')
      #if self.drsbjd: # comment out because the sa cannot be calculated with str
         #self.drsbjd = repr(self.drsbjd)
      self.drsberv = hdr.get(k_berv, np.nan)
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.blaze = hdr.get(HIERDRS+'BLAZE FILE', 0)
      self.drift = hdr.get(HIERDRS+'DRIFT RV', np.nan)
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan
      else:
         self.e_drift = hdr.get(HIERDRS+'DRIFT NOISE', np.nan)
      self.drift = 0   # reset, because drifts are already included in the wavemap
      # self.e_drift is kept to account now for wavemap uncertainty

      if self.instname == 'HARPS':
         # read the comment
         #if pfits==2: fileid = hdr['DATE-OBS']; self.timeid=fileid #ok!? not always
         if pfits:
            if hasattr(hdr, 'comments'):   # pfits==2 or pyfits.__version__>'2.' https://github.com/astropy/pyregion/issues/20
               fileid = hdr.comments['MJD-OBS']
            else:
               fileid = str(hdr.ascardlist()['MJD-OBS'])
         else: fileid = hdrio.get_comment('MJD-OBS')
         self.timeid = fileid[fileid.index('(')+1 : fileid.index(')')]
      elif self.instname == 'HARPN':
         self.timeid = fileid = hdr['FILENAME'][6:29]
         hdr['OBJECT'] = hdr[HIERINST+'OBS TARG NAME'] # HARPN has no OBJECT keyword
      elif self.instname == 'SOPHIE':
         self.timeid = fileid = self.dateobs # hdr['DATE']

      #calmode = hdr.get('IMAGETYP',0).split(",")[:2]
      calmode = hdr.get(HIERINST+'DPR TYPE','NOTFOUND').split(',')[:2]
      self.calmode = ','.join(calmode)
      calmodedict = {'STAR,WAVE': 'OBJ,CAL', 'STAR,DARK': 'OBJ,SKY'}
      if self.calmode in calmodedict: self.calmode = calmodedict[self.calmode]

      hdr['OBJECT'] = hdr.get(HIERINST+'TARG NAME', 'FOX')
      self.header = self.hdr = hdr # self.header will be set to None

   scan_ccf(self)


def data(self, orders, pfits=True):
   hdr = self.hdr
   drs = self.drs
   if 1:  # read order data
      if hasattr(self, 'hdu'):   # read directly
         f = self.hdu.getdata(o=orders)
         if not drs:
            e = self.hdu.getdata('SIG', o=orders)
            w = self.hdu.getdata('WAVE', o=orders)
      else:
         if not hasattr(self, 'hdulist'):
            scan(self, self.filename)
         f = self.hdulist[0 if drs else 'SPEC'].section[orders]
         if not drs:
            e = self.hdulist['WAVE'].section[orders]
            w = self.hdulist['SIG'].section[orders]

      if not drs:
         f *= 100000
         e *= 100000

      bpmap = np.isnan(f).astype(int)   # flag 1 for nan
      if not drs: bpmap[e==0] |= flag.nan
      if drs:
         # print " applying wavelength solution ", file
         # omax = self.hdu['SPEC'].NAXIS1
         omax = hdr[self.HIERDRS+'CAL LOC NBO'] # 72 for A and 71 for B
         d = hdr[self.HIERDRS+'CAL TH DEG LL']
         xmax = 4077
         x = np.empty((d+1, xmax), 'int64')
         x[0].fill(1)                               # x[0,*] = x^0 = 1,1,1,1,1,...
         x[1] = np.arange(xmax)                     #        = x^1 = 0,1,2,3,4,...
         for i in range(1,d): x[i+1] = x[i] * x[1]  #        = x^i
         if not hasattr(self, 'A'):
         #A = np.array([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))],dtype='float64').reshape(omax,d+1) #slow 30 ms
            self.A = np.reshape([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))], (omax,d+1)) #slow 30 ms
         w = np.dot(self.A[orders], x)  # wavelength lambda
         #w = self.hdulist['WAVE_A'].section[orders] # better? stored as float32, 10m/s numeric difference
         e = np.sqrt(np.where(bpmap, 0., 5**2 * 6 + np.abs(f, dtype=float)))
         blaze = self.hdulist['BLAZE_A'].section[orders]
         f = f / blaze
         e = e / blaze

      with np.errstate(invalid='ignore'):
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f > 300000] |= flag.sat    # estimate for saturation level:
                                       # HARPS.2004-10-03T01:30:44.506.fits:
                                       # last order: e2ds_B: 346930 (x=2158) raw: 62263 (y=1939)

      w = airtovac(w)
      return w, f, e, bpmap


