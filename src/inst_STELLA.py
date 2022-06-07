from read_spec import *
#from read_spec import Inst
# Instrument parameters
from readmultispec import *
from astropy.time import Time

name = 'STELLA'
obsname = 'teide'
obsloc = dict(lat=28.3, lon= -16.5097, elevation=2390.)

# slicer = 'all'
# slicer = 'odd'
slicer = 'even'

if slicer == 'all':
    iomax = 164 # NAXIS2
else:
    iomax = 82 # NAXIS2
snmax = 500
oset = ':'

# maskfile = 'telluric_mask_carm_short.dat'

pat = '*_botzfxsEcd.fits'


def scan(self, s, pfits=True):
    """
    Returns
    -------
    namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
    Example
    -------
    >>> read_carm_vis(filename)

    """
    hdulist = self.hdulist = pyfits.open(s)
    self.s = s  #propagate the filename

    self.header = hdr = hdulist[0].header
    self.instname = hdr['INSTRUME']
    self.instname = 'STELLA'
    self.drsberv = hdr.get('VCORRECT', np.nan)
    self.drsbjd = hdr.get('HJD', np.nan)
    self.dateobs = hdr['DATE-OBS']
    self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd

    self.drift = hdr.get('RVDRIFT', np.nan)
    self.e_drift = hdr.get('RVDRERR', np.nan)
    self.sn55 = hdr.get('SNEST', np.nan)

    self.fileid = hdr.get('DATE-OBS', 0)
    self.timeid = self.fileid
    self.calmode = 'Th-Ar'

    self.ccf.rvc = hdr.get('VHELIO', np.nan)
    self.ccf.err_rvc = np.nan

    # self.ra = hdr['RA_HMS']#falta pasar a grados
    # self.de = hdr['DEC_DMS']#falta pasar a grados
    self.ra = hdr['OBJRA']
    self.de = hdr['OBJDEC']

    self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

    self.obs.lon = obsloc['lon']
    self.obs.lat = obsloc['lat']

    self.airmass = hdr.get('AIRMASS', np.nan)
    self.exptime = hdr['EXPTIME']
    self.tmmean = 0.5


def data(self, orders, pfits=True):
    hdums = readmultispec(self.s)  #hdumultispec
    if slicer == 'all':
        m = np.arange(0, hdums['wavelen'][::-1].shape[0], 1)
    if slicer == 'odd':
        m = np.arange(0, hdums['wavelen'][::-1].shape[0], 2) + 1
    if slicer == 'even':
        m = np.arange(0, hdums['wavelen'][::-1].shape[0], 2)

    f = hdums['flux'][0,:,:][::-1][m,:][orders,:]
    w = hdums['wavelen'][::-1][m,:][orders,:]
    e = hdums['flux'][2,:,:][::-1][m,:][orders,:]

    bpmap = np.isnan(f).astype(int)            # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
