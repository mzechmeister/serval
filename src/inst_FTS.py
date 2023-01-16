from read_spec import *
from opus2py import get_interferogram_metadata, read_spectrum, get_data_set_x_data
# http://www.astro.physik.uni-goettingen.de/~kwessel/opus2py/opus2py.tar.gz
from astropy.time import Time

# Instrument parameters
name = __name__[5:]
obsloc = dict(lat=30.671666694444447, lon=255.97833330555557, elevation=2075)   # elp: McDonald http://www.astropy.org/astropy-data/coordinates/sites.json

pat = '*.0'

iomax = 68       # The FTS will be split in so main chunks
pmax =  24000
oset = '25:67'

maskfile = 'telluric_mask_carm_short.dat'


# Instrument read functions
def scan(self, filename, pfits=True):
    hdulist = meta = get_interferogram_metadata(filename)

    self.header = hdr = pyfits.PrimaryHDU().header
    hdr['OBJECT'] = 'unknown'
    self.instname = 'FTS'

    timetag = meta.date + b"T" + meta.time[0:12]
    self.utc = datetime.datetime.strptime(timetag.decode(), '%d/%m/%YT%H:%M:%S.%f')
    self.dateobs = self.utc.isoformat()[:-3]
    self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd

    self.drsberv = 0.
    self.drsbjd = self.mjd + 2_450_000.5

    self.drift = np.nan
    self.e_drift = np.nan
    self.sn55 = 50   # hdr.get('SNR', 0)
    self.fileid = self.dateobs
    self.calmode = ''
    self.timeid = self.fileid

    self.ra = 0. # hdr['RA']
    self.de = 0. # hdr['DEC']
    self.airmass = hdr.get('AIRMASS', np.nan)
    self.exptime = 0.  # hdr['EXPTIME']
    self.tmmean = 0.5   # Is there a weighted midpoint? Date)

def data(self, orders=None, pfits=True):

    # read order data
    f = np.array(read_spectrum(self.filename))
    w = 1e8 / np.array(get_data_set_x_data(self.filename, 'SAMPLED-SPECTRUM'))

    chunk_size = f.size // iomax

    f = f[:iomax*chunk_size].reshape(iomax, -1)[orders]
    w = w[:iomax*chunk_size].reshape(iomax, -1)[orders]

    e = 0*f + np.median(f, axis=f.ndim-1, keepdims=True) / 50  # arbitrary choice

    bpmap = np.isnan(f).astype(int)   # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
