from astropy.io import fits 
from read_spec import *
from astropy.coordinates import EarthLocation
from astropy.time import Time

name = inst = __name__[5:]

pat = '*.fits'  # _wave, _sp

ldt = EarthLocation.of_site("DCT")
obsloc = dict(lat=ldt.lat.value, lon=ldt.lon.value, elevation=ldt.height.value)

iomax = 86
oset = '37:76'
pmax = 7920 - 500

def scan(self, s, pfits=True):
    hdu = hdulist = self.hdulist = fits.open(s)
    self.header = self.hdr = hdr0 = hdu[0].header
    self.instname = hdr0['INSTRMNT']
    self.drsbjd = hdu[1].header.get('BARYMJD', np.nan) + 2_400_000.5
    self.drsberv = 0 # wavelength are aready in barycentric
    self.dateobs = hdr0['DATE-SHT']   # Time shutter opened
    self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd

    self.fileid = hdr0['OBS_ID']
    self.calmode = hdu[1].header['WAVE-CAL']
    self.ra = hdr0['RA']
    self.de = hdr0['DEC']
    self.airmass = hdr0['AIRMASS']
    self.exptime = float(hdr0['AEXPTIME'])   # floats are strings in header?
    self.timeid = self.dateobs.replace(" ", "T")
    self.sn55 = 55

    if 'tellurics' not in hdulist[1].data.names:
        self.flag |= sflag.lowSN

def data(self, orders=None, pfits=True):
    datatbl = self.hdulist[1].data

    w_air = datatbl['bary_excalibur'][orders] # This is the barycentric excalibur corrected wavelengths.
    if np.isnan(datatbl['bary_excalibur'][37,2000]):
        w_air = datatbl['bary_wavelength'][orders]
        #print("bary_excalibur is NaN. Using bary_wavelength.")
    s = datatbl['spectrum'][orders]
    continuum_model = datatbl['continuum'][orders]
    e = uncertainty = datatbl['uncertainty'][orders]
    tellurics = datatbl['tellurics'][orders]

    f = s / continuum_model / tellurics
    #e = e / continuum_model / tellurics  <= right?
#    w_air = datatbl['wavelength_'+beam].reshape(37,-1).astype('float64')[::-1][orders]
    w = airtovac(w_air)
    bpmap = np.isnan(f).astype(int)      # flag 1 for nan
    bpmap |= np.isnan(w)                 # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
