from read_spec import *

# Instrument parameters
name = __name__[5:]
obsloc = dict(lat=42.9333, lon= 0.1333, elevation=2869.4)

pat = '*.fits'

iomax = 39
pmax =  7800 - 300
#oset = ':66'

maskfile = 'telluric_mask_carm_short.dat'


# Instrument read functions
def scan(self, s, pfits=True):
   hdulist = self.hdulist = pyfits.open(s)
   self.header = hdr = hdulist[0].header
 
   self.instname = hdr['INSTRUM']
   if self.instname == 'NARVAL':
       self.instname = 'NEONARVAL'   # fix unreliable keyword self.drsberv = hdr.get('BERV', np.nan)
   self.drsberv = hdr.get('BERV', np.nan)
   self.drsbjd = hdr.get('BTLA', np.nan)
   self.dateobs = hdr['DATE-FIT']
   self.mjd = hdr.get('DATE_JUL')
   self.drift = np.nan
   self.e_drift = np.nan
   self.sn55 = hdr.get('SNR', 50)   # always 50
   self.fileid = self.dateobs
   self.calmode = "%s,%s,%s" % (hdr.get('OBSTYPE', ''), hdr.get('STOKNAME', ''), hdr.get('P_NAME2', '')) #
   self.timeid = self.fileid

   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr.get('AIRMASS', np.nan)
   self.exptime = hdr['EXPTIME']
   self.tmmean = 0.5 # hdr.get(HIERARCH+'CARACAL TMEAN', 0.5)

def data(self, orders=None, pfits=True):
   hdulist = self.hdulist
   # read order data
   f = hdulist[1].data['Intensity'].reshape(39,7800)[::-1][orders]
   w_air = hdulist[1].data['Wavelength1'].reshape(39,7800)[::-1][orders]
   e = hdulist[1].data['Error'].reshape(39,7800)[::-1][orders]
   w = airtovac(w_air)

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan

   return w, f, e, bpmap


