from read_spec import *
from astropy.time import Time

# Instrument parameters
name = __name__[5:]   # FIES (iraf), FIES_CERES

# obsname = (28.75728, -17.88508, 2382) from wiki
obsname = "not" # for barycorrpy
obsname = "lapalma" # for barycorrpy

pat = '*[vs][ep].fits'  # _wave, _sp

iomax = {'FIES': 79, 'FIES_CERES': 76}[name]
pmax =  {'FIES': 2048, 'FIES_CERES': 2102}[name] - 300
if name == 'FIES_CERES':
   oset = ':66'   # CERES background subtraction problem in blue orders

name = 'FIES'

maskfile = 'telluric_mask_carm_short.dat'


# Instrument read functions
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
   HIERARCH = 'HIERARCH '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   self.header = hdr = hdulist[0].header
   self.drs = hdr.get('PIPELINE', 'DRS')

   if self.drs == 'CERES':
      self.instname = hdr['INST']
      self.drsberv = hdr.get('BARYCENTRIC CORRECTION (KM/S)', np.nan)
      self.drsbjd = hdr.get('MBJD', np.nan) + 2400000.5  # same as MJD!?
      self.dateobs = hdr['HIERARCH SHUTTER START DATE'] + 'T' + hdr['HIERARCH SHUTTER START UT']
      self.mjd = hdr.get('HIERARCH MJD')
      self.drift = np.nan
      self.e_drift = np.nan
      self.fileid = self.dateobs
      self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))
      self.timeid = self.fileid
      self.ccf.rvc = hdr.get('RV', np.nan)
      self.ccf.err_rvc = hdr.get('RV_E', np.nan)

      self.ra = hdr['HIERARCH RA']
      self.de = hdr['HIERARCH DEC']
      self.airmass = hdr.get('HIERARCH TARG AIRMASS START', np.nan)
      self.exptime = hdr['HIERARCH TEXP (S)']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5
      # estimate SNR
      f = hdulist[0].section[1]
      self.snr = np.nanmedian(np.abs(f[:,1:]/(f[:,1:]- f[:,:-1])), axis=1)
      self.sn55 = self.snr[55]
      hdr['OBJECT'] = hdr['HIERARCH TARGET NAME']
   else:
      # IRAF header, *_wave.fits
      self.instname = hdr['INSTRUME']
      self.drsberv = hdr.get('BERV', np.nan)
      self.drsbjd = hdr.get('HJD', np.nan)
      self.dateobs = hdr['DATE-OBS']
      self.mjd = hdr.get('JD', Time(self.dateobs, scale='utc').jd) - 2400000.5
      # for HPF spectra the drift is already included in the wavelength solution
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', hdr.get(HIERARCH+'CARACAL DRIFT RV', np.nan))
      self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT RVERR', np.nan))
      self.sn55 = hdr.get('SNR 36', 50)
      self.fileid = hdr.get('DATE-OBS', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.timeid = self.fileid
      self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]

      self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
      self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)

      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5
      # estimate SNR
      f = hdulist[0].section
      if len(f.hdu.shape) == 3:
          # We have cube. Assuming f[0] flux, f[1] similar flux (linear) f[2] error
          f = f[0]
      self.snr = np.median(np.abs(f[:,1:]/(f[:,1:]- f[:,:-1])), axis=1)
      self.sn55 = self.snr[55]

def data(self, orders=None, pfits=True):
   hdulist = self.hdulist
   # read order data
   if self.drs == 'CERES':
      f = hdulist[0].section[1][orders]
      w = hdulist[0].section[0][orders] / (1.00004)
      e = 1/np.sqrt(hdulist[0].section[2][orders])
   else:
      # stupid iraf header
      gg = readmultispec(self.filename, reform=True, quiet=True)
      # "".join(self.header['WAT2_0*'].values()).split("spec")
      w = airtovac(gg['wavelen'][orders])
      f = hdulist[0].section
      if len(f.hdu.shape) == 3:
          e = f[2][orders]
          f = f[0][orders]
      else:
          f = f[orders]
          e = f*0 + 1/self.snr[orders, np.newaxis] *0.5

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan

   return w, f, e, bpmap




# helper scripts to read iraf stored wavelengths
# from
# https://github.com/kgullikson88/General/blob/master/readmultispec.py

"""readmultispec.py
Read IRAF (echelle) spectrum in multispec format from a FITS file.
Can read most multispec formats including linear, log, cubic spline,
Chebyshev or Legendre dispersion spectra.
Usage: retdict = readmultispec(fitsfile, reform=True)
Inputs:
fitfile     Name of the FITS file
reform      If true (the default), a single spectrum dimensioned
            [4,1,NWAVE] is returned as flux[4,NWAVE].  If false,
            it is returned as a 3-D array flux[4,1,NWAVE].
Returns a dictionary with these entries:
flux        Array dimensioned [NCOMPONENTS,NORDERS,NWAVE] with the spectra.
            If NORDERS=1, array is [NCOMPONENTS,NWAVE]; if NCOMPONENTS is also
            unity, array is [NWAVE].  (This can be changed
            using the reform keyword.)  Commonly the first dimension
            is 4 and indexes the spectrum, an alternate version of
            the spectrum, the sky, and the error array.  I have also
            seen examples where NCOMPONENTS=2 (probably spectrum and
            error).  Generally I think you can rely on the first element
            flux[0] to be the extracted spectrum.  I don't know of
            any foolproof way to figure out from the IRAF header what the
            various components are.
wavelen     Array dimensioned [NORDERS,NWAVE] with the wavelengths for
            each order.
header      The full FITS header from pyfits.
wavefields  [NORDERS] List with the analytical wavelength
            description (polynomial coefficients, etc.) extracted from
            the header.  This is probably not very useful but is
            included just in case.
History:
Created by Rick White based on my IDL readechelle.pro, 2012 August 15
Apologies for any IDL-isms that remain!
"""

import numpy as np
import astropy.io.fits as pyfits



def nonlinearwave(nwave, specstr, verbose=False):
    """Compute non-linear wavelengths from multispec string
    
    Returns wavelength array and dispersion fields.
    Raises a ValueError if it can't understand the dispersion string.
    """

    fields = specstr.split()
    if int(fields[2]) != 2:
        raise ValueError('Not nonlinear dispersion: dtype=' + fields[2])
    if len(fields) < 12:
        raise ValueError('Bad spectrum format (only %d fields)' % len(fields))
    wt = float(fields[9])
    w0 = float(fields[10])
    ftype = int(fields[11])
    if ftype == 3:

        # cubic spline

        if len(fields) < 15:
            raise ValueError('Bad spline format (only %d fields)' % len(fields))
        npieces = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if verbose:
            print 'Dispersion is order-%d cubic spline' % npieces
        if len(fields) != 15 + npieces + 3:
            raise ValueError('Bad order-%d spline format (%d fields)' % (npieces, len(fields)))
        coeff = np.asarray(fields[15:], dtype=float)
        # normalized x coordinates
        s = (np.arange(nwave, dtype=float) + 1 - pmin) / (pmax - pmin) * npieces
        j = s.astype(int).clip(0, npieces - 1)
        a = (j + 1) - s
        b = s - j
        x0 = a ** 3
        x1 = 1 + 3 * a * (1 + a * b)
        x2 = 1 + 3 * b * (1 + a * b)
        x3 = b ** 3
        wave = coeff[j] * x0 + coeff[j + 1] * x1 + coeff[j + 2] * x2 + coeff[j + 3] * x3

    elif ftype == 1 or ftype == 2:

        # chebyshev or legendre polynomial
        # legendre not tested yet

        if len(fields) < 15:
            raise ValueError('Bad polynomial format (only %d fields)' % len(fields))
        order = int(fields[12])
        pmin = float(fields[13])
        pmax = float(fields[14])
        if verbose:
            if ftype == 1:
                print 'Dispersion is order-%d Chebyshev polynomial' % order
            else:
                print 'Dispersion is order-%d Legendre polynomial (NEEDS TEST)' % order
        if len(fields) != 15 + order:
            # raise ValueError('Bad order-%d polynomial format (%d fields)' % (order, len(fields)))
            if verbose:
                print 'Bad order-%d polynomial format (%d fields)' % (order, len(fields))
                print "Changing order from %i to %i" % (order, len(fields) - 15)
            order = len(fields) - 15
        coeff = np.asarray(fields[15:], dtype=float)
        # normalized x coordinates
        pmiddle = (pmax + pmin) / 2
        prange = pmax - pmin
        x = (np.arange(nwave, dtype=float) + 1 - pmiddle) / (prange / 2)
        p0 = np.ones(nwave, dtype=float)
        p1 = x
        wave = p0 * coeff[0] + p1 * coeff[1]
        for i in range(2, order):
            if ftype == 1:
                # chebyshev
                p2 = 2 * x * p1 - p0
            else:
                # legendre
                p2 = ((2 * i - 1) * x * p1 - (i - 1) * p0) / i
            wave = wave + p2 * coeff[i]
            p0 = p1
            p1 = p2

    else:
        raise ValueError('Cannot handle dispersion function of type %d' % ftype)

    return wave, fields


def readmultispec(fitsfile, reform=True, quiet=False):
    """Read IRAF echelle spectrum in multispec format from a FITS file
    
    Can read most multispec formats including linear, log, cubic spline,
    Chebyshev or Legendre dispersion spectra
    
    If reform is true, a single spectrum dimensioned 4,1,NWAVE is returned
    as 4,NWAVE (this is the default.)  If reform is false, it is returned as
    a 3-D array.
    """

    fh = pyfits.open(fitsfile)
    try:
        header = fh[0].header
        flux = fh[0].data
    finally:
        fh.close()
    temp = flux.shape
    nwave = temp[-1]
    if len(temp) == 1:
        nspec = 1
    else:
        nspec = temp[-2]

    # first try linear dispersion
    try:
        crval1 = header['crval1']
        crpix1 = header['crpix1']
        cd1_1 = header['cd1_1']
        ctype1 = header['ctype1']
        if ctype1.strip() == 'LINEAR':
            wavelen = np.zeros((nspec, nwave), dtype=float)
            ww = (np.arange(nwave, dtype=float) + 1 - crpix1) * cd1_1 + crval1
            for i in range(nspec):
                wavelen[i, :] = ww
            # handle log spacing too
            dcflag = header.get('dc-flag', 0)
            if dcflag == 1:
                wavelen = 10.0 ** wavelen
                if not quiet:
                    print 'Dispersion is linear in log wavelength'
            elif dcflag == 0:
                if not quiet:
                    print 'Dispersion is linear'
            else:
                raise ValueError('Dispersion not linear or log (DC-FLAG=%s)' % dcflag)

            if nspec == 1 and reform:
                # get rid of unity dimensions
                flux = np.squeeze(flux)
                wavelen.shape = (nwave,)
            return {'flux': flux, 'wavelen': wavelen, 'header': header, 'wavefields': None}
    except KeyError:
        pass

    # get wavelength parameters from multispec keywords
    try:
        wat2 = header['wat2_*']
        count = len(wat2)
    except KeyError:
        raise ValueError('Cannot decipher header, need either WAT2_ or CRVAL keywords')

    # concatenate them all together into one big string
    watstr = []
    for i in range(len(wat2)):
        # hack to fix the fact that older pyfits versions (< 3.1)
        # strip trailing blanks from string values in an apparently
        # irrecoverable way
        # v = wat2[i].value
        v = wat2[i]
        v = v + (" " * (68 - len(v)))  # restore trailing blanks
        watstr.append(v)
    watstr = ''.join(watstr)

    # find all the spec#="..." strings
    specstr = [''] * nspec
    for i in range(nspec):
        sname = 'spec' + str(i + 1)
        p1 = watstr.find(sname)
        p2 = watstr.find('"', p1)
        p3 = watstr.find('"', p2 + 1)
        if p1 < 0 or p1 < 0 or p3 < 0:
            raise ValueError('Cannot find ' + sname + ' in WAT2_* keyword')
        specstr[i] = watstr[p2 + 1:p3]

    wparms = np.zeros((nspec, 9), dtype=float)
    w1 = np.zeros(9, dtype=float)
    for i in range(nspec):
        w1 = np.asarray(specstr[i].split(), dtype=float)
        wparms[i, :] = w1[:9]
        if w1[2] == -1:
            raise ValueError('Spectrum %d has no wavelength calibration (type=%d)' %
                             (i + 1, w1[2]))
            # elif w1[6] != 0:
            #    raise ValueError('Spectrum %d has non-zero redshift (z=%f)' % (i+1,w1[6]))

    wavelen = np.zeros((nspec, nwave), dtype=float)
    wavefields = [None] * nspec
    for i in range(nspec):
        # if i in skipped_orders:
        #    continue
        verbose = (not quiet) and (i == 0)
        if wparms[i, 2] == 0 or wparms[i, 2] == 1:
            # simple linear or log spacing
            wavelen[i, :] = np.arange(nwave, dtype=float) * wparms[i, 4] + wparms[i, 3]
            if wparms[i, 2] == 1:
                wavelen[i, :] = 10.0 ** wavelen[i, :]
                if verbose:
                    print 'Dispersion is linear in log wavelength'
            elif verbose:
                print 'Dispersion is linear'
        else:
            # non-linear wavelengths
            wavelen[i, :], wavefields[i] = nonlinearwave(nwave, specstr[i],
                                                         verbose=verbose)
        wavelen *= 1.0 + wparms[i, 6]
        if verbose:
            print "Correcting for redshift: z=%f" % wparms[i, 6]
    if nspec == 1 and reform:
        # get rid of unity dimensions
        flux = np.squeeze(flux)
        wavelen.shape = (nwave,)
    return {'flux': flux, 'wavelen': wavelen, 'header': header, 'wavefields': wavefields}
 
 
