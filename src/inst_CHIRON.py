from read_spec import *
from astropy.time import Time

#Instrument parameters
name = __name__[5:]
obsloc = dict(lat=-30.169286, lon=-70.806789, elevation=2200.)   # approximate height - don't know exact elevation

pat = 'achi*.fits'

iomax = 59   # Number of orders
pmax = 3200 - 200   # Snip off the pixels at the ends of the order

# No maskfile for telluric lines

# Instrument read functions
def scan(self, s, pfits=True, verb = False):
   hdulist = self.hdulist = pyfits.open(s)
   self.header = hdr = hdulist[0].header
 
   self.instname = hdr['FPA']
   if self.instname == 'CHIRON  ':
       self.instname = 'CHIRON'
   self.drsberv = hdr.get('BERV', np.nan)
   self.drsbjd = hdr.get('BTLA', np.nan)
   self.dateobs = hdr['DATE-OBS']
   self.mjd = hdr.get('DATE_JUL', Time(self.dateobs, format='isot', scale='utc').mjd)
   self.drift = np.nan
   self.e_drift = np.nan
   self.sn55 = 10 * hdr.get('EMNETINT', np.nan)   # a guess for S/N from exposure meter
   self.fileid = hdr['OBSID']
   self.calmode = "%s,%s,%s" % (hdr.get('OBSTYPE', ''), hdr.get('STOKNAME', ''), hdr.get('P_NAME2', ''))
   self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr['AIRMASS']
   self.exptime = hdr['EXPTIME']
   self.tmmean = 0.5   # Flux weighted mean point
   # Probably, it should be computed from EMMNWB  = '2025-10-07T09:03:46.419' / mean time with bckgrd subtraction
   self.timeid = self.dateobs

def data(self, orders, pfits=True):
    hdulist = self.hdulist

    hdu_data = hdulist[0]

    w = hdu_data.section[orders,:,0].astype(float)   # Uh, wavelengths only in float32!?
    f = hdu_data.section[orders,:,1].astype(float)
    e = np.sqrt(1.**2 + np.abs(f))   # RON ~4.66 and gain ~1.3 should be used

    bpmap = 1 * np.isnan(f)
    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg

    return w, f, e, bpmap
