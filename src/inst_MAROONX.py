from read_spec import *
from read_spec import Inst
from calcspec import redshift
import h5py
from astropy.time import Time

# Instrument parameters
pat = "*_0003.fits"

name = inst = __name__[5:]
obsname = 'geminiN'   # for barycorrpy
obsloc = dict(lat=19.82396, lon=155.46984, elevation=4213.)   # https://en.wikipedia.org/wiki/Gemini_Observatory

oset = np.s_[6:] #if not join_slice else np.s_[6/2:]
iomax =24
pmin = 300
pmax = 3800

join_slice = False   # experimental
# The idea is to create one high S/N template per order, instead of two noisy templates.
# Runs slower, because wavelengths needs to be re-sorted.

oabs_list = map(str, range(122,98,-1))

if join_slice:
   oset = np.s_[6/2:]
   iomax /= 2
   pmin *= 2
   pmax *= 2

def scan(self, s, pfits=True, verb=False):
   drs = self.drs

   self.hdf = h5py.File(s)
   timetag = s[-33:-18]
   self.utc = datetime.datetime.strptime(timetag, '%Y%m%dT%H%M%S')
   self.exptime = float(s[-8:-4])

   self.instname = 'MAROONX'

   self.airmass = 1.
   self.dateobs = self.utc.isoformat()
   self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd
   self.ra = None #hdr['RA']
   self.de = None # hdr['DEC']

   self.obs.lon = obsloc['lon']
   self.obs.lat = obsloc['lat']
   self.obs.elevation = obsloc['elevation']

   self.tmmean = 0.5

   self.drsbjd = self.mjd
   self.drsberv = 0.0
   self.sn55 = 20.0
   self.drift = 0.0

   self.timeid = timetag
   self.calmode = 'NOTFOUND'

   self.header = self.hdr = pyfits.Header(
       {'SIMPLE': True,
        'NAXIS': 2,
        'NAXIS1': 1000,
        'NAXIS2': 2000,
        'OBJECT': 'NONE'})  # self.header will be set to None


def data(self, orders, pfits=True):
   hdr = self.hdr
   if not hasattr(self, 'hdf'):
      scan(self, self.filename)
   hdf = self.hdf

   if join_slice:
       # Not implemented!
       # join slices
       w = hdf['wavelength_solution']['fiber_2'].data.reshape(170/2, 9141*2)
       ii = np.argsort(w[orders], axis=-1)
       oo = np.arange(np.shape(w)[0])[orders,np.newaxis]
       w = w[oo,ii]
       f = self.hdulist['SCIDATA'].data.reshape(170/2, 9141*2)[oo,ii]
       e = self.hdulist['ERRDATA'].data.reshape(170/2, 9141*2)[oo,ii]
       bpmap = 1 * (self.hdulist['QUALDATA'].data.reshape(170/2, 9141*2)[oo,ii] > 0)
   else:
       if orders == np.s_[:]:
          w = np.vstack([hdf['wavelength_solution']['fiber_2'][o] for o in oabs_list])
          f = np.vstack([hdf['optimal_extraction']['fiber_2'][o] for o in oabs_list])
          e = 1/np.sqrt(np.vstack([hdf['optimal_var']['fiber_2'][o] for o in oabs_list]))
       else:
          w = np.array(hdf['wavelength_solution']['fiber_2'][oabs_list[orders]])
          f = np.array(hdf['optimal_extraction']['fiber_2'][oabs_list[orders]])
          e = 1/np.sqrt(hdf['optimal_var']['fiber_2'][oabs_list[orders]])

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
      bpmap[f > 300000] |= flag.sat    # unchecked estimate for saturation level

   return w, f, e, bpmap


