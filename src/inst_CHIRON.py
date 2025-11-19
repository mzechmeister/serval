from read_spec import *
from astropy.time import Time

#Instrument parameters
name = __name__[5:]
obsloc = dict(lat = -30.169286, lon = -70.806789, elevation = 2200.) #approximate height - don't know exact elevation

pat = 'achi*.fits'

iomax = 59 #Number of orders
pmax = 3200 - 200 #Snip off the pixels at the ends of the order

#What is oset? Do I need to "join slices" for CHIRON if I'm evaluating each order separately?

#No maskfile for telluric lines

# Instrument read functions
def scan(self, s, pfits=True, verb = False):
   hdulist = self.hdulist = pyfits.open(s)
   self.header = hdr = hdulist[0].header
 
   self.instname = hdr['FPA'] #Name of spectrograph
   if self.instname == 'CHIRON  ':
       self.instname = 'CHIRON'
   self.drsberv = hdr.get('BERV', default = np.nan) #Barycentric earth radial velocity
   self.drsbjd = hdr.get('BTLA', default = np.nan) #What does this mean?
   self.dateobs = hdr['DATE-OBS']
   self.mjd = hdr.get('DATE_JUL', default = Time(self.dateobs, format = 'isot', scale = 'utc').mjd)
   self.drift = np.nan
   self.e_drift = np.nan
   self.sn55 = 10 * hdr.get('EMNETINT', default = np.nan)   # a guess for S/N?
   self.fileid = hdr['OBSID']
   self.calmode = "%s,%s,%s" % (hdr.get('OBSTYPE', ''), hdr.get('STOKNAME', ''), hdr.get('P_NAME2', '')) #I have no idea what this means
   self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr['AIRMASS']
   self.exptime = hdr['EXPTIME']
   self.tmmean = 0.5 #Flux weighted mean point (default = 0.5 for most spectrographs)
   self.timeid = self.dateobs

def data(self, orders, pfits=True): #Need to edit this section (What do I write when this spectrograph has multiple orders?)
    hdulist = self.hdulist

    hdu_data = hdulist[0]

    w = hdu_data.section[orders,:,0].astype(float)   # Uh, wavelengths only in float32!?
    f = hdu_data.section[orders,:,1].astype(float)
    e = np.sqrt(1.**2 + np.abs(f))   # RON 4.66 and gain 1.3 should be used

    bpmap = 1 * np.isnan(f)
    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg

    return w, f, e, bpmap
