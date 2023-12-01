from read_spec import *
#from read_spec import Inst
# Instrument parameters

name = 'SPIRou'
obsname = 'cfht'
obsloc = dict(lat=19.8253, lon= -155.4689, elevation=4213.)

iomax = 49 # NAXIS2
snmax = 500
oset = ':'

maskfile = 'telluric_mask_carm_short.dat'

pat = '*t.fits'   # Corrected from tellurics


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
    self.header = hdr = hdulist[0].header
    hdr1 = hdulist[1].header
    self.instname = hdr['INSTRUME']
    self.drsberv = hdr1.get('BERV', np.nan)
    self.drsbjd = hdr1.get('BJD', np.nan)
    self.dateobs = hdr['DATE-OBS'] + 'T' + hdr['UTC-OBS'][:7]
    self.mjd = hdr['MJD-OBS']
    self.drift = hdr1.get('RV_DRIFT', np.nan)
    self.e_drift = np.nan
    self.sn55 = hdr.get('SPEMSNR', np.nan) # Estimated, where? how?

    self.fileid = hdr.get(self.dateobs, 0)
    self.timeid = self.fileid
    self.calmode = hdr.get('DPRTYPE', '')

    self.ccf.rvc = hdr.get('CCFMNRV', np.nan)
    self.ccf.err_rvc = hdr.get('DVRMS_CC', np.nan)*1e-3

    self.ra = hdr['RA_DEG']
    self.de = hdr['DEC_DEG']

    self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S')

    self.obs.lon = obsloc['lon']
    self.obs.lat = obsloc['lat']

    self.airmass = hdr.get('AIRMASS', np.nan)
    self.exptime = hdr['EXPTIME']
    self.tmmean = 0.5


def data(self, orders, pfits=True):
    hdulist = self.hdulist

    f = hdulist['FluxAB'].section[orders]
    w = hdulist['WaveAB'].section[orders]
    e = np.ones_like(w)
    blaze = hdulist['BlazeAB'].section[orders]
    f = f / blaze
    e = 0*f + np.median(f, axis=f.ndim-1, keepdims=True) / 50  # arbitrary choice

    bpmap = np.isnan(f).astype(np.uint64)            # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
