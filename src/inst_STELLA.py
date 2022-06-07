from read_spec import *
#from read_spec import Inst
# Instrument parameters
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
pmax =  2167 - 300

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
from astropy.io import fits as pyfits


def nonlinearwave(nwave, specstr, printTXT=False, verbose=False):
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
        if printTXT==False:#verbose==True:
            print ('Dispersion is order-%d cubic spline' % npieces)
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
        if printTXT!=False:#verbose:
            if ftype == 1:
                print( 'Dispersion is order-%d Chebyshev polynomial' % order)
            else:
                print ('Dispersion is order-%d Legendre polynomial (NEEDS TEST)' % order)
        if len(fields) != 15 + order:
            # raise ValueError('Bad order-%d polynomial format (%d fields)' % (order, len(fields)))
            if printTXT!=False:#verbose:
                print ('Bad order-%d polynomial format (%d fields)' % (order, len(fields)))
                print ("Changing order from %i to %i" % (order, len(fields) - 15))
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


def readmultispec(fitsfile, reform=True, quiet=False, printTXT=False, verbose=False):
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
                if printTXT!=False:#not quiet:
                    print ('Dispersion is linear in log wavelength')
            elif dcflag == 0:
                if printTXT!=False:#not quiet:
                    print ('Dispersion is linear')
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
                if printTXT!=False:#verbose==True:
                    print ('Dispersion is linear in log wavelength')
            elif printTXT!=False:#verbose==True:
                print ('Dispersion is linear')
        else:
            # non-linear wavelengths
            wavelen[i, :], wavefields[i] = nonlinearwave(nwave, specstr[i],
                                                         verbose=verbose)
        wavelen *= 1.0 + wparms[i, 6]
        if printTXT!=False:#verbose==True:
            print ("Correcting for redshift: z=%f" % wparms[i, 6])
    if nspec == 1 and reform:
        # get rid of unity dimensions
        flux = np.squeeze(flux)
        wavelen.shape = (nwave,)
    return {'flux': flux, 'wavelen': wavelen, 'header': header, 'wavefields': wavefields}

