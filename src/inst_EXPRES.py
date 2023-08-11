from astropy.io import fits 
from read_spec import *
from astropy.coordinates import EarthLocation

name = inst = __name__[5:]

pat = '*.fits'  # _wave, _sp

ldt = EarthLocation.of_site("DCT")
obsloc = dict(lat=ldt.lat.value, lon=ldt.lon.value, elevation=ldt.height.value)

iomax = 86
oset = '37:76'
pmax = 7920 - 500

def scan(self, s, pfits=True):
    hdu = hdulist = self.hdulist = fits.open(s)
    self.instname = hdu[0].header['INSTRMNT']
    self.header = self.hdr = hdu[0].header
    self.drsbjd = hdu[1].header.get('BARYMJD', np.nan) + 2_450_000.5  # 0.5 right?
    self.drsberv = 0 # wavelength are aready in barycentric
    self.dateobs = hdu[0].header['DATE-OBS']
    
    self.fileid = hdu[0].header['OBS_ID']
    
    self.calmode = hdu[1].header['WAVE-CAL']
    self.ra = hdu[0].header['RA']
    self.de = hdu[0].header['DEC']
    self.airmass = hdu[0].header['AIRMASS']
    self.exptime = hdu[0].header['AEXPTIME']
    self.timeid = self.dateobs.replace(" ", "T")
    self.sn55 = 55
 
def data(self, orders=None, pfits=True):
    hdulist = self.hdulist
 
    w_air = hdulist[1].data['bary_excalibur'][orders] # This is the barycentric excalibur corrected wavelengths.
    s = hdulist[1].data['spectrum'][orders]
    continuum_model = hdulist[1].data['continuum'][orders]
    e = uncertainty = hdulist[1].data['uncertainty'][orders]
    tellurics = hdulist[1].data['tellurics'][orders]

    f = s / continuum_model / tellurics
    #e = e / continuum_model / tellurics  <= right?
#    w_air = hdulist[1].data['wavelength_'+beam].reshape(37,-1).astype('float64')[::-1][orders]
    w = airtovac(w_air)
    bpmap = np.isnan(f).astype(int)      # flag 1 for nan
    bpmap |= np.isnan(w)                 # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

    return w, f, e, bpmap
