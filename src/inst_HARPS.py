from read_spec import *
from read_spec import Inst
# Instrument parameters

name = inst = __name__[5:]
obsname = {'HARPS': 'eso', 'HARPN': 'lapalma'}[name] # for barycorrpy
iomax = {'HARPS': 72, 'HARPN': 69}[name]

#maskfile = 'telluric_mask_carm_short.dat'

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

   """
   drs = self.drs
   if '.tar' in s:
      s = file_from_tar(s, inst=inst, fib=self.fib, pfits=pfits)

   if isinstance(s, str) and '.gz' in s:
      # if s is isinstance(s,tarfile.ExFileObject) then s.position will change !? resulting in:
      # *** IOError: Empty or corrupt FITS file
      pfits = True

   if 1:
      HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + {'HARPS': 'ESO ', 'HARPN': 'TNG '}[inst]
      k_tmmean = {'HARPS': HIERINST + 'INS DET1 TMMEAN', 'HARPN': HIERINST + 'EXP_METER_A EXP CENTROID'}[inst]
      # In old HARPN the keyword is different and the value absolute
      #k_tmmean = {'HARPS': HIERINST + 'INS DET1 TMMEAN', 'HARPN': HIERINST + 'EXP1 TMMEAN'}[inst]
      if drs:
         self.HIERDRS = HIERDRS = HIERINST + 'DRS '
         k_sn55 = HIERDRS + 'SPE EXT SN55'
         k_berv = HIERDRS + 'BERV'
         k_bjd = HIERDRS + 'BJD'
      else:
         k_sn55 = HIERARCH + 'FOX SNR 55'
         k_berv = 'E_BERV'
         k_bjd = 'E_BJD'

      if pfits is True:  # header with pyfits
         self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
         hdr = self.hdulist[0].header # pyfits.getheader(s)
      elif pfits==2:     # a faster version
         args = ('INSTRUME', 'OBJECT', 'MJD-OBS', 'DATE-OBS', 'OBS TARG NAME', 'EXPTIME',
                 'MJD-OBS', 'FILENAME', 'RA', 'DEC', k_tmmean, HIERINST+'DPR TYPE',
                 HIERINST+'DPR TECH', HIERINST+'INS MODE', HIERINST+'OBS TARG NAME')
         args += (k_bjd, k_berv, k_sn55)
         # args += (HIERINST+'OBS PI-COI NAME', HIERINST+'OBS PROG ID')
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
      if self.instname not in ('HARPS', 'HARPN'):
         pause('\nWARNING: inst should be HARPS or HARPN, but got: '+self.inst+'\nSee option -inst for available inst.') 
      self.HIERARCH = HIERARCH

      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.mjd = hdr['MJD-OBS']
      self.dateobs = hdr['DATE-OBS']
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

      self.obs.lon = -70.7345
      self.obs.lat = -29.2584

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
      self.drift = hdr.get(HIERDRS+'DRIFT RV USED', np.nan)
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

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
      #calmode = hdr.get('IMAGETYP',0).split(",")[:2]
      calmode = hdr.get(HIERINST+'DPR TYPE','NOTFOUND').split(',')[:2]
      self.calmode = ','.join(calmode)
      calmodedict = {'STAR,WAVE': 'OBJ,CAL', 'STAR,DARK': 'OBJ,SKY'}
      if self.calmode in calmodedict: self.calmode = calmodedict[self.calmode]

      if hdr[HIERINST+'DPR TECH'] == 'ECHELLE,ABSORPTION-CELL':
         self.flag |= sflag.iod
      if hdr[HIERINST+'INS MODE'] == 'EGGS':
         self.flag |= sflag.eggs

      hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
      self.header = self.hdr = hdr # self.header will be set to None


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
         xmax = 4096
         x = np.empty((d+1, xmax), 'int64')
         x[0].fill(1)                               # x[0,*] = x^0 = 1,1,1,1,1,...
         x[1] = np.arange(xmax)                     #        = x^1 = 0,1,2,3,4,...
         for i in range(1,d): x[i+1] = x[i] * x[1]  #        = x^i
         if not hasattr(self, 'A'):
         #A = np.array([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))],dtype='float64').reshape(omax,d+1) #slow 30 ms
            self.A = np.reshape([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))], (omax,d+1)) #slow 30 ms
         w = np.dot(self.A[orders], x)  # wavelength lambda
         e = np.sqrt(np.where(bpmap, 0., 5**2 * 6 + np.abs(f, dtype=float)))

      with np.errstate(invalid='ignore'):
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f > 300000] |= flag.sat    # estimate for saturation level:
                                       # HARPS.2004-10-03T01:30:44.506.fits:
                                       # last order: e2ds_B: 346930 (x=2158) raw: 62263 (y=1939)

      w = airtovac(w)
      return w, f, e, bpmap


