from read_spec import *

# Instrument parameters
name = __name__[5:]
obsloc = dict(lat=42.9333, lon= 0.1333, elevation=2869.4)

pat = '*.fits.fz'

iomax = 68
pmax =  3800
oset = ':52'

maskfile = 'telluric_mask_carm_short.dat'


# Instrument read functions
def scan(self, filename, pfits=True):
    hdulist = self.hdulist = pyfits.open(filename)
    self.header = hdr = hdulist[0].header

    self.instname = hdr['INSTRUME']
    if self.instname.startswith('fa'):
        self.instname = 'NRES'
    self.drsberv = hdr.get('BARYCORR', np.nan) / 1000   # [km/s]
    self.drsbjd = hdr.get('TCORR', np.nan)
    self.dateobs = hdr['DATE-OBS']
    self.mjd = hdr.get('MJD-OBS')
    self.drift = np.nan
    self.e_drift = np.nan
    self.sn55 = hdr.get('SNR', 0)
    self.fileid = hdr['ORIGNAME']
    self.calmode = hdr.get('OBJECTS', '')
    self.timeid = self.fileid

    self.ra = hdr['RA']
    self.de = hdr['DEC']
    self.airmass = hdr.get('AIRMASS', np.nan)
    self.exptime = hdr['EXPTIME']
    self.tmmean = 0.5   # Is there weighted midpoint?   TCORR = Exposure Mid-Time (Barycentric Julian Date)

    self.ccf.rvc = hdr.get('RV', np.nan)
    self.ccf.err_rvc = hdr.get('RVERR', np.nan)

def data(self, orders=None, pfits=True):
    hdulist = self.hdulist
    # read order data
    fib = np.s_[::-2]   # assuming two fibres, but there could be three?
    f = hdulist['SPECTRUM'].data['flux'][fib][orders]
    w = hdulist['SPECTRUM'].data['wavelength'][fib][orders]
    e = hdulist['SPECTRUM'].data['uncertainty'][fib][orders]
    w[w==0] = np.nan

    bpmap = np.isnan(f).astype(int)   # flag 1 for nan
    bpmap |= np.isnan(w)

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap

