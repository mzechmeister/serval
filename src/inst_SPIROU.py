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

oset = [0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 40, 42, 43, 44, 45, 46, 47]


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
    self.sn55 = hdr1.get('EXTSN033', 'SPEMSNR')  # Order 33 ~1674nm
    # self.sn55 = hdr.get('SPEMSNR', np.nan) # Estimated, where? how?

    self.fileid = self.dateobs
    self.timeid = self.fileid
    self.calmode = hdr.get('DPRTYPE', '')

    # to extract the RV and the drift value it is necessary to open the ccf file
    try:
        hdr2 = pyfits.open(s[:-6] + 'v.fits')[1].header # if the file has the same name
        self.drift = hdr2.get('RV_DRIFT', np.nan)
        self.ccf.rvc = hdr2.get('RV_OBJ', np.nan)
        self.ccf.err_rvc = hdr2.get('DVRMS_CC', np.nan)*1e-3
    except:
        self.drift = np.nan
        self.ccf.rvc = np.nan
        self.ccf.err_rvc = np.nan
    self.e_drift = np.nan

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
    w = 10*hdulist['WaveAB'].section[orders]   # nm => Angstrom
    e = np.ones_like(w)
    blaze = hdulist['BlazeAB'].section[orders]
    f = f / blaze
    e = 0*f + np.nanmedian(f, axis=f.ndim-1, keepdims=True) / self.sn55  # divided by the S/N

    bpmap = np.isnan(f).astype(int)            # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
