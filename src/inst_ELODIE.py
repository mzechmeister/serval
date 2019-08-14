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
iomax = 67
pat = '*_s2d.fits'

#maskfile = 'telluric_mask_carm_short.dat'

def scan_ccf(self):
   ccffiles = glob.glob(self.filename[:-10]+'*_ccf.fits')  # replace e2ds with *_ccf
   ccf = self.ccf
   if ccffiles:
      ccf.filename = ccffiles[0]
      hdr = pyfits.getheader(ccf.filename)
      ccf.rvc = hdr.get('CTRFIT', np.nan)   # [km/s], not clear whether this value is drift corrected?
      #ccf.err_rvc = hdr.get('CTRSIG', np.nan)  # [km/s]
      ccf.contrast = hdr.get('AMPFIT', np.nan) # to be normalised with CTEFIT?
      ccf.fwhm = hdr.get('SIGFIT', np.nan)
      ccf.mask = hdr.get('MASNAM', np.nan)
      self.drsberv = hdr.get('BERV', np.nan)  # [km/s]
      #self.drsbjd = hdr.get('BJD', np.nan)


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

   developed for
   VERSION = '1.3     '           / TACOS_VERSION

   http://www.obs-hp.fr/www/guide/elodie/tgmet-paper/NODE4.HTML
   http://www.obs-hp.fr/www/guide/elodie/manuser.html
   """
   drs = self.drs
   if isinstance(s, str) and '.gz' in s:
      # if s is isinstance(s,tarfile.ExFileObject) then s.position will change !? resulting in:
      # *** IOError: Empty or corrupt FITS file
      pfits = True

   if 1:
      k_tmmean = 'INS PM TMEAN'
      if drs:
         k_sn55 = 'SN'
         k_berv = 'BERV'
         k_bjd = 'BJD'

      self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
      hdr = self.hdulist[0].header # pyfits.getheader(s)

      self.instname = hdr['INSTRUME'] #if self.drs else 'HARPS'
      if self.instname != name:
         pause('\nWARNING: inst should be HARPS or HARPN, but got: '+self.inst+'\nSee option -inst for available inst.') 

      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.mjd = hdr['MJD-OBS']   # observ. start
      self.hdulist[0].verify('silentfix')
      self.dateobs = hdr['DATE-OBS']

      self.ra = hdr['ALPHA']
      self.de = hdr['DELTA']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S')

      self.obs.lon = obsloc['lon']
      self.obs.lat = obsloc['lat']

      if k_tmmean not in hdr:
         warnings.warn('Warning: old HARPN data? Setting tmmean to 0.5!')
      self.tmmean = hdr.get(k_tmmean, 0.5)

      self.drsbjd = hdr.get(k_bjd)
      if self.drsbjd is None:
         self.drsbjd = hdr.get('MJD-OBS')

      self.drsberv = hdr.get('BERV', np.nan)  # [km/s]
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.blaze = hdr.get('BLAZE FILE', 0)
      self.drift = hdr.get('DRIFT RV', np.nan)
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

      self.timeid = fileid = hdr['H_AFR000']

      self.calmode = hdr.get('IMATYP', 'UNKNOWN')

      hdr['OBJECT'] = hdr.get('OBJECT')
      self.header = self.hdr = hdr # self.header will be set to None

   scan_ccf(self)

def data(self, orders, pfits=True):
   hdr = self.hdr
   drs = self.drs
   if 1:  # read order data
      if hasattr(self, 'hdu'):   # read directly
         f = self.hdu.getdata(o=orders)
      else:
         if not hasattr(self, 'hdulist'):
            scan(self, self.filename)
         f = self.hdulist[0].section[orders]

      bpmap = np.isnan(f).astype(int)   # flag 1 for nan

      if drs:
         omax = int(hdr['NORDER']) # 67
         deg_o = int(hdr['DEGOLL'])
         deg_x = int(hdr['DEGXLL'])
         xmax = int(hdr['NAXIS1']) # 1024
         x = np.arange(xmax)                     #        = x^1 = 0,1,2,3,4,...
         xn = x / (xmax-1.) * 2 - 1.
         Ti = np.polynomial.chebyshev.chebvander(xn, deg_x)
         if not hasattr(self, 'ATj'):
         #A = np.array([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))],dtype='float64').reshape(omax,d+1) #slow 30 ms
            A = np.zeros((deg_x+1, deg_o+1))
            A[0,0] = 0.25 * hdr['COELL1']
            A[1:,0] = [0.5 * hdr['COELL'+str(1+(deg_o+1)*i)] for i in range(1, deg_x+1)]
            A[0,1:] = [0.5 * hdr['COELL'+str(1+j)] for j in range(1, deg_o+1)]
            for i in range(1, deg_x+1):
                A[i,1:] = [hdr['COELL'+str(1+j+(deg_o+1)*i)] for j in range(1, deg_o+1)]
            o = np.arange(omax)
            m = hdr['LLOFFSET'] + hdr['LLPARS'] * o
            minv = 1 / m
            minvn = (minv - minv.min()) / (minv.max() - minv.min()) * 2 - 1.
            Tj = np.polynomial.chebyshev.chebvander(minvn, deg_o)
            self.ATj = A.dot(Tj.T).T / m[:,np.newaxis]
            #.T.dot(Ti.T) / m[orders,np.newaxis]  #  1 / m sum A Ti Tj

         # wavelength lambda
         # np.einsum('ij,xi,oj->ox', self.A, Ti, Tj[orders])
         #w = A.dot(Tj[orders].T).T.dot(Ti.T) / m[orders,np.newaxis]  #  1 / m sum A Ti Tj
         w = self.ATj[orders].dot(Ti.T)

         blaze = self.hdulist['BLZ'].section[orders]
         e = np.sqrt(np.where(bpmap, 0., 5**2 * 6 + np.abs(f*blaze, dtype=float))) / blaze
#         pause()
      with np.errstate(invalid='ignore'):
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f > 300000] |= flag.sat    # estimate for saturation level:

      w = airtovac(w)
      return w, f, e, bpmap


