#! /usr/bin/env python
__author__ = 'Mathias Zechmeister'
__version__ = '2018-01-10'

import datetime
import glob
import gzip
import os
import sys
import tarfile
import time
import warnings
from collections import namedtuple

import pyfits
#import fitsio
import numpy as np

import brv_we14idl
import brv_we14html
from pause import pause, stop
from gplot import *
import sunrise

class nameddict(dict):
   """
   Examples
   --------
   >>> nameddict({'a':1, 'b':2})
   {'a': 1, 'b': 2}
   >>> x = nameddict(a=1, b=2)
   >>> x.a
   1
   >>> x['a']
   1
   >>> x.translate(3)
   ['a', 'b']

   """
   __getattr__ = dict.__getitem__

   def translate(self, x):
      return [name for name,f in self.items() if (f & x) or f==x==0]

# bpmap flagging
flag = nameddict(
   ok=       0, # good pixel
   nan=      1, # nan flux in pixel of spectrum
   neg=      2, # significant negative flux in pixel of spectrum (f < -3*f_err < 0 && )
   sat=      4, # too high flux (saturated)
   atm=      8, # telluric at wavelength of spectrum
   sky=     16, # sky emission at wavelength of spectrum
   out=     32, # region outside the template
   clip=    64, # clipped value
   lowQ=   128, # low Q in stellar spectrum (mask continuum)
   badT=   256, # bad corresponding region in the template spectrum
)

# spectrum flagging
sflag = nameddict(
   iod=      2,
   eggs=     4,
   dist=    16, # coordinates too much off
   daytime= 32, # not within nautical twilight
   lowSN=   64, # too low S/N
   hiSN=   128, # too high S/N
   led=    256, # LED on during observation (CARM_VIS)
   rvnan=  512
)


# bpmap flags
flag_cosm = flag.sat  # @ FEROS for now use same flag as sat
def_wlog = True
brvref = ['DRS', 'MH', 'WEhtml', 'WEidl', 'WE']

class Spectrum:
   """
   Examples
   --------
   >>> from read_spec import *
   >>> x = Spectrum('/media/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/HARPS/red_DRS/gj551/HARPS.2013-05-06T01:54:55.600.tar', fib='A',orders=-1,drs=1)
   >>> x = Spectrum('/home/zechmeister/data/harps/gj588/data/HARPS.2013-05-11T08:46:58.520_e2ds_A.fits')

   """
   brvref = 'MH' # can be overwritten by serval.py depending on user input
   def __init__(self, filename, inst='HARPS', pfits=True, orders=None, wlog=def_wlog, drs=False, fib=None, targ=None, verb=False):
      """
      pfits : fits reader
             True: pyfits (safe, e.g. to get the full header of the start refence spectrum)
             2: my own simple pure python fits reader, faster than pyfits
             It was developed for HARPS e2ds, but the advent ADP.*.tar was difficult.

      """
      self.filename = filename
      self.drsname = 'DRS'
      self.drs = drs
      self.flag = 0
      self.bflag = 0
      self.tmmean = 0
      self.f = None
      self.w = None
      self.e = None
      self.bpmap = None
      self.mod = None
      self.fib = fib
      self.fo = None
      self.wo = None
      self.eo = None
      self.bo = None
      self.drift = np.nan
      self.e_drift = np.nan
      self.utc = None
      self.ra = None
      self.de = None
#      self.obs = type('specdata', (object,), {'lat': None, 'lon': None})
      self.obs = type('specdata', (object,), dict(lat=None, lon=None))
      self.airmass = np.nan
      if '.gz' in filename: pfits=True

      self.ccf = type('ccf',(), dict(rvc=np.nan, err_rvc=np.nan, bis=np.nan, fwhm=np.nan, contrast=np.nan, mask=0, header=0))

      read_spec(self, filename, inst=inst, pfits=pfits, verb=verb)   # scan fits header
      if inst != self.inst:
         pause('WARNING:', filename, 'from', self.inst, ', but mode is', inst)
      self.obj = self.header['OBJECT']

      ### Barycentric correction ###

      if targ and targ.name == 'cal':
         self.bjd, self.berv = self.drsbjd, 0.
      elif targ and targ.ra:  # unique coordinates
         if self.brvref == 'MH':
            # fastest version
            #sys.path.append(os.environ['HOME']+'/programs/BarCor/')
            sys.path.insert(1, sys.path[0]+os.sep+'BarCor')            # bary now in src/BarCor
            import bary
            self.bjd, self.berv = bary.bary(self.dateobs, targ.ra, targ.de, inst, epoch=2000, exptime=self.exptime*2* self.tmmean, pma=targ.pmra, pmd=targ.pmde)
         if self.brvref in ('WEhtml', 'WEidl', 'WE'):
            # Wright & Eastman (2014) via online or idl request
            # cd /home/raid0/zechmeister/idl/exofast/bary
            # export ASTRO_DATA=/home/raid0/zechmeister/
            #idl -e 'print, bjd2utc(2457395.24563d, 4.58559072, 44.02195596), fo="(D)"'
            jd_utc = [self.mjd + 2400000.5 + self.exptime*self.tmmean/24./3600]
            ra = (targ.ra[0] + targ.ra[1]/60. + targ.ra[2]/3600.) * 15  # [deg]
            de = (targ.de[0] + np.copysign(targ.de[1]/60. + targ.de[2]/3600., targ.de[0]))       # [deg]
            obsname = {'CARM_VIS':'ca', 'CARM_NIR':'ca', 'FEROS':'eso', 'HARPS':'eso', 'HARPN':'lapalma'}[inst]
            if self.brvref == 'WE':
               # pure python version
               import brv_we14py
               self.bjd, self.berv = brv_we14py.bjdbrv(jd_utc=jd_utc[0], ra=ra, dec=de, obsname=obsname, pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0])
            elif self.brvref == 'WEhtml':
               self.bjd = brv_we14html.utc2bjd(jd_utc=jd_utc, ra=ra, dec=de)
               self.berv = brv_we14html.bvc(jd_utc=jd_utc, ra="%s+%s+%s"%targ.ra, dec="%s+%s+%s"%targ.de, obsname='ca', pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0], raunits='hours', deunits='degrees')[0]
            else:
               self.bjd, self.berv = brv_we14idl.bjdbrv(jd_utc=jd_utc[0], ra=ra, dec=de, obsname=obsname, pmra=targ.pmra, pmdec=targ.pmde, parallax=0., rv=0., zmeas=[0])

            self.berv /= 1000.   # m/s to km/s
      else:
         self.bjd, self.berv = self.drsbjd, self.drsberv
         self.brvref = self.drsname

      if self.fib == 'B':
         self.berv = np.nan
         self.drsberv = np.nan

      self.header['HIERARCH SERVAL BREF'] = (self.brvref, 'Barycentric code')
      self.header['HIERARCH SERVAL BJD'] = (self.bjd, 'Barycentric Julian Day')
      self.header['HIERARCH SERVAL BERV'] = (self.berv, '[km/s] Barycentric correction')
      #self.header['HIERARCH SERVAL RA'] = (targ.ra[0], 'Barycentric code')
      #self.header['HIERARCH SERVAL DE'] = (targ.ra[0], 'Barycentric code')
      #self.header['HIERARCH SERVAL PMRA'] = (targ.ra[0], 'Barycentric code')
      #self.header['HIERARCH SERVAL PMDE'] = (targ.ra[0], 'Barycentric code')

      if self.utc:
         date = self.utc.year, self.utc.month, self.utc.day
         sunris = sunrise.sun(*date, lon=self.obs.lon, lat=self.obs.lat)
         sunset = sunrise.sun(*date, lon=self.obs.lon, lat=self.obs.lat, rise=False)
         ut = self.utc.hour + self.utc.minute/60. +  self.utc.second/3600.
         utend = ut + self.exptime/3600.
         dark = ((sunset < ut <= utend < sunris) or   # |****S__UT__R***|
                 (sunris < sunset < ut <= utend) or   # |__R*******S__UT|
                 (utend < sunris < sunset < ut)  or   # |T__R********S_U|
                 (ut < utend < sunris < sunset))      # |UT__R********S_|
         #if not dark:
            #self.flag |= sflag.daytime

      if orders is not None: self.read_data(orders=orders, wlog=wlog)

   def __get_item__(self, order):
      # spectrum needs to be dict like
      return self.get_data(order)

   def get_data(self, orders=np.s_[:], wlog=def_wlog, verb=False, **kwargs):
      """Returns only data."""
      o = orders
      if self.w:
         w, f, e, b = self.w[o], self.f[o], self.e[o], self.bpmap[o]
      else:
         w, f, e, b = read_spec(self, self.filename, inst=self.inst, orders=orders, **kwargs)
         w = np.log(w) if wlog else w.astype(np.float)
         f = f.astype(float)
         e = e.astype(float)
         self.bflag |= np.bitwise_or.reduce(b.ravel())
      return type('specdata', (object,),
               dict(w=w, f=f, e=e, bpmap=b, berv=self.berv, o=o))

   def read_data(self, verb=False, **kwargs):
      """Read only data, no header, store it as an attribute."""
      data = self.get_data(**kwargs)
      if isinstance(kwargs.get('orders'), int):
         self.wo, self.fo, self.eo, self.bo = data.w, data.f, data.e, data.bpmap
      else:
         self.w, self.f, self.e, self.bpmap = data.w, data.f, data.e, data.bpmap



def read_spec(self, s, inst='HARPS', plot=False, **kwargs):
   #print s, inst
   if '.tar' in s: s = file_from_tar(s, inst=inst, fib=self.fib, **kwargs)
   if 'HARP' in inst:  sp = read_harps(self, s, inst=inst, **kwargs)
   elif inst == 'CARM_VIS':  sp = read_carm_vis(self, s, **kwargs)
   elif inst == 'CARM_NIR':  sp = read_carm_nir(self, s, **kwargs)
   elif inst == 'FEROS': sp = read_feros(self, s, **kwargs)
   elif inst == 'FTS': sp = read_fts(self, s, **kwargs)
   else:
      return None

   if plot:
      gplot_set("set xlabel 'wavelength'; set ylabel 'intensity'")
      for o in range(len(sp.w)):
         gplot(sp.w[o], sp.f[o], " w lp t 'order %i'"%o)
         pause(o)
   return sp

def write_template(filename, flux, wave, header=None, hdrref=None, clobber=False):
   if not header and hdrref: header = pyfits.getheader(hdrref)
   hdu = pyfits.PrimaryHDU(header=header)
   warnings.resetwarnings() # supress nasty overwrite warning http://pythonhosted.org/pyfits/users_guide/users_misc.html
   warnings.filterwarnings('ignore', category=UserWarning, append=True)
   hdu.writeto(filename, clobber=clobber, output_verify='fix')
   warnings.resetwarnings()
   warnings.filterwarnings('always', category=UserWarning, append=True)

   if isinstance(flux, np.ndarray):
      pyfits.append(filename, flux)
      pyfits.append(filename, wave)
   else:
      # pad arrays with zero to common size
      maxpix = max(arr.size for arr in flux if isinstance(arr, np.ndarray))
      flux_new = np.zeros((len(flux), maxpix))
      wave_new = np.zeros((len(flux), maxpix))
      for o,arr in enumerate(flux):
          if isinstance(arr, np.ndarray): flux_new[o,:len(arr)] = arr
      for o,arr in enumerate(wave):
          if isinstance(arr, np.ndarray): wave_new[o,:len(arr)] = arr
      pyfits.append(filename, flux_new)
      pyfits.append(filename, wave_new)

   pyfits.setval(filename, 'EXTNAME', value='SPEC', ext=1)
   pyfits.setval(filename, 'EXTNAME', value='WAVE', ext=2)
   #fitsio.write(filename, flux)

def write_res(filename, datas, extnames, header='', hdrref=None, clobber=False):
   if not header and hdrref: header = pyfits.getheader(hdrref)
   hdu = pyfits.PrimaryHDU(header=header)
   warnings.resetwarnings() # supress nasty overwrite warning http://pythonhosted.org/pyfits/users_guide/users_misc.html
   warnings.filterwarnings('ignore', category=UserWarning, append=True)
   hdu.writeto(filename, clobber=clobber, output_verify='fix')
   warnings.resetwarnings()
   warnings.filterwarnings('always', category=UserWarning, append=True)

   for i,extname in enumerate(extnames):
     data = datas[extname]
     if isinstance(data, np.ndarray):
        pyfits.append(filename, data)
     else:
        1/0

     pyfits.setval(filename, 'EXTNAME', value=extname, ext=i+1)
   #fitsio.write(filename, flux)

def write_fits(filename, data, header='', hdrref=None, clobber=True):
   if not header and hdrref: header = pyfits.getheader(hdrref)
   warnings.resetwarnings() # supress nasty overwrite warning http://pythonhosted.org/pyfits/users_guide/users_misc.html
   warnings.filterwarnings('ignore', category=UserWarning, append=True)
   pyfits.writeto(filename, data, header, clobber=clobber, output_verify='fix')
   warnings.resetwarnings()
   warnings.filterwarnings('always', category=UserWarning, append=True)

def read_template(filename):
   hdu = pyfits.open(filename)
   #return hdu[1].data, hdu[0].data  # wave, flux
   return hdu[2].data, hdu[1].data, hdu[0].header  # wave, flux

def read_harps_ccf(s):
   ccf = namedtuple('ccf', 'rvc err_rvc bis fwhm contrast mask')
   if ".tar" in s:
      tar = tarfile.open(s)
      extr = None
      for member in tar.getmembers():
         if 'A.fits' in member.name:
            if '_ccf_' in member.name and not extr: extr = member
            if '_bis_' in member.name: extr = member; is_bis = 1  # prefer bis
      if not extr: return ccf(0,0,0,0,0,0)
      s = tar.extractfile(extr)
   else:
      s = glob.glob(s.replace("_e2ds","_ccf_*"))
      if s: s = s[0]
      else: return ccf(*[np.nan]*6)
      #else: return ccf(0,0,0,0,0,0)
   # ccf = namedtuple('ccf', 'rvc err_rvc bis fwhm contrast mask header')
   HIERARCH = 'HIERARCH '

   if 1:
      hdr = imhead(s, HIERARCH+'ESO DRS CCF RVC', HIERARCH+'ESO DRS CCF CONTRAST', HIERARCH+'ESO DRS CCF FWHM', HIERARCH+'ESO DRS CCF MASK', HIERARCH+'ESO DRS DVRMS',HIERARCH+'ESO DRS BIS SPAN')
   elif 0:
      if ".tar" in s:
         s = tar.extractfile(extr)
         hdulist = pyfits.open(s)
         hdr = hdulist[0].header
         tar.close()
         #hdr = pyfits.getheader(s) # doesn't work for file like object?
   else:
      tar.extract(extr, path='tarfits')
      os.system('mv tarfits/* tmp.fits ')
      data,hdr = fitsio.read('tmp.fits',header=1)
      HIERARCH = ''
      tar.close()

   rvc = hdr[HIERARCH+'ESO DRS CCF RVC']   # [km/s]
   contrast = hdr.get(HIERARCH+'ESO DRS CCF CONTRAST', np.nan)
   fwhm = hdr[HIERARCH+'ESO DRS CCF FWHM']
   mask = hdr[HIERARCH+'ESO DRS CCF MASK']
   e_rvc = hdr.get(HIERARCH+'ESO DRS DVRMS', np.nan) / 1000.   # [km/s]
   bis = hdr.get(HIERARCH+'ESO DRS BIS SPAN', np.nan)
   return ccf(rvc, e_rvc, bis, fwhm, contrast, mask)

def read_harps(self, s, inst='HARPS', orders=None, pfits=True, verb=False):
   """
   SYNTAX: read_harps(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   drs = self.drs
   if isinstance(s, str) and '.gz' in s:
      # if s is isinstance(s,tarfile.ExFileObject) then s.position will change !? resulting in:
      # *** IOError: Empty or corrupt FITS file
      pfits = True

   if orders is None or self.header is None or (pfits==True and not hasattr(self, 'hdulist')):
      HIERARCH = 'HIERARCH '
      HIERINST = HIERARCH + {'HARPS': 'ESO ', 'HARPN': 'TNG '}[inst]
      k_tmmean = {'HARPS': HIERINST + 'INS DET1 TMMEAN', 'HARPN': HIERINST +  'EXP_METER_A EXP CENTROID'}[inst]
      if drs:
         self.HIERDRS = HIERDRS = HIERINST + 'DRS '
         k_sn55 = HIERDRS + 'SPE EXT SN55'
         k_berv = HIERDRS + 'BERV'
         k_bjd = HIERDRS + 'BJD'
      else:
         k_sn55 = HIERARCH + 'FOX SNR 55'
         k_berv = 'E_BERV'
         k_bjd = 'E_BJD'

      if pfits is True:  # header with pyfits
         self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
         hdr = self.hdulist[0].header # pyfits.getheader(s)
      elif pfits==2:     # a faster version
         args = ('INSTRUME', 'OBJECT', 'DATE-OBS', 'OBS TARG NAME', 'EXPTIME',
                 'MJD-OBS', 'FILENAME', 'RA', 'DEC', k_tmmean, HIERINST+'DPR TYPE',
                 HIERINST+'DPR TECH', HIERINST+'INS MODE', HIERINST+'OBS TARG NAME')
         args += (k_bjd, k_berv, k_sn55)
         if drs:
            args += (HIERDRS+'BLAZE FILE', HIERDRS+'DRIFT RV USED',
                     HIERDRS+'CAL TH DEG LL', HIERDRS+'CAL LOC NBO',
                     HIERDRS+'CAL TH COEFF LL')

         hdr = imhead(s, *args)
         self.hdu = getext(s)
      else:
         #hdr = fitsio.read_header(s) no faster?
         self.f, hdrio = fitsio.read(s, header=True)
         hdr = dict((key, val.strip() if type(val) is str else val) for key,val in dict(hdrio).iteritems())
         HIERARCH = ''

      #self.drs = 'DRS CAL LOC NBO' in "".join(hdr.keys())  # check DRS or FOX
      self.inst = hdr['INSTRUME'] #if self.drs else 'HARPS'
      if self.inst not in ('HARPS', 'HARPN'):
         pause('\nWARNING: inst should be HARPS or HARPN, but got: '+self.inst+'\nSee option -inst for available inst.') 
      self.HIERARCH = HIERARCH

      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.dateobs = hdr['DATE-OBS']
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S.%f')

      self.obs.lon = -70.7345
      self.obs.lat = -29.2584

      self.tmmean = hdr[k_tmmean]

      self.drsbjd = hdr.get(k_bjd)
      if self.drsbjd is None:
         self.drsbjd = hdr.get('MJD-OBS')
      #if self.drsbjd: # comment out because the sa cannot be calculated with str
         #self.drsbjd = repr(self.drsbjd)
      self.drsberv = hdr.get(k_berv, np.nan)
      self.sn55 = hdr.get(k_sn55, np.nan)
      self.blaze = hdr.get(HIERDRS+'BLAZE FILE', 0)
      self.drift = hdr.get(HIERDRS+'DRIFT RV USED', np.nan)
      if abs(self.drift) > 1000:
         # sometimes there are crazy drift values ~2147491.59911, e.g. 2011-06-15T08:11:13.465
         self.drift = np.nan

      if self.inst == 'HARPS':
         # read the comment
         #if pfits==2: fileid = hdr['DATE-OBS']; self.timeid=fileid #ok!? not always
         if pfits:
            if hasattr(hdr, 'comments'):   # pfits==2 or pyfits.__version__>'2.' https://github.com/astropy/pyregion/issues/20
               fileid = hdr.comments['MJD-OBS']
            else:
               fileid = str(hdr.ascardlist()['MJD-OBS'])
         else: fileid = hdrio.get_comment('MJD-OBS')
         self.timeid = fileid[fileid.index('(')+1 : fileid.index(')')]
      elif self.inst == 'HARPN':
         self.timeid = fileid = hdr['FILENAME'][6:29]
         hdr['OBJECT'] = hdr[HIERINST+'OBS TARG NAME'] # HARPN has no OBJECT keyword
      #calmode = hdr.get('IMAGETYP',0).split(",")[:2]
      calmode = hdr.get(HIERINST+'DPR TYPE','NOTFOUND').split(',')[:2]
      self.calmode = ','.join(calmode)
      calmodedict = {'STAR,WAVE': 'OBJ,CAL', 'STAR,DARK': 'OBJ,SKY'}
      if self.calmode in calmodedict: self.calmode = calmodedict[self.calmode]

      if hdr[HIERINST+'DPR TECH'] == 'ECHELLE,ABSORPTION-CELL':
         self.flag |= sflag.iod
      if hdr[HIERINST+'INS MODE'] == 'EGGS':
         self.flag |= sflag.eggs

      hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
      self.header = hdr

   hdr = self.header

   if verb: print "read_harps:", self.timeid, hdr['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode

   if orders is not None:  # read order data
      if pfits == 2:   # read directly
         f = self.hdu.getdata(o=orders)
         if not drs:
            e = self.hdu.getdata('SIG', o=orders)
            w = self.hdu.getdata('WAVE', o=orders)
      elif pfits:
         f = self.hdulist[0 if drs else 'SPEC'].section[orders]
         if not drs:
            e = self.hdulist['WAVE'].section[orders]
            w = self.hdulist['SIG'].section[orders]

      if not drs:
         f *= 100000
         e *= 100000

      bpmap = np.isnan(f).astype(int)   # flag 1 for nan
      if not drs: bpmap[e==0] |= flag.nan
      if drs:
         # print " applying wavelength solution ", file
         # omax = self.hdu['SPEC'].NAXIS1
         omax = hdr[self.HIERDRS+'CAL LOC NBO'] # 72 for A and 71 for B
         d = hdr[self.HIERDRS+'CAL TH DEG LL']
         xmax = 4096
         x = np.empty((d+1, xmax), 'int64')
         x[0].fill(1)                               # x[0,*] = x^0 = 1,1,1,1,1,...
         x[1] = np.arange(xmax)                     #        = x^1 = 0,1,2,3,4,...
         for i in range(1,d): x[i+1] = x[i] * x[1]  #        = x^i
         if not hasattr(self, 'A'):
         #A = np.array([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))],dtype='float64').reshape(omax,d+1) #slow 30 ms
            self.A = np.reshape([hdr[self.HIERDRS+'CAL TH COEFF LL'+str(i)] for i in range(omax*(d+1))], (omax,d+1)) #slow 30 ms
         w = np.dot(self.A[orders], x)  # wavelength lambda
         e = np.sqrt(np.where(bpmap, 0., 5**2 * 6 + np.abs(f, dtype=float)))

      with np.errstate(invalid='ignore'):
         bpmap[f < -3*e] |= flag.neg      # flag 2 for zero and negative flux
         bpmap[f > 300000] |= flag.sat    # estimate for saturation level:
                                       # HARPS.2004-10-03T01:30:44.506.fits:
                                       # last order: e2ds_B: 346930 (x=2158) raw: 62263 (y=1939)

      w = airtovac(w)
      return w, f, e, bpmap


def read_carm_vis(self, s, orders=None, pfits=True, verb=True):
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
   hdulist = pyfits.open(s) # slow 30 ms
   if orders is None:
      self.header = hdr = hdulist[0].header
      if 'HIERARCH CARACAL DRIFT FP REF' in hdr: del hdr['HIERARCH CARACAL DRIFT FP REF']
      self.inst = hdr['INSTRUME'][0:4] + '_VIS'
      self.drsberv = hdr.get('HIERARCH CARACAL BERV', np.nan)
      # BJD for stars, JD for for calibration products
      self.drsbjd = hdr.get('HIERARCH CARACAL BJD', np.nan) + 2400000
      #if isinstance(self.drsbjd, str):
         #self.drsbjd = 0.0   # workaround for MJD-OBS bugs (empty or missing fractional digits) @2016-Jan
      self.dateobs = hdr['DATE-OBS']
      # dateobs is used by MH, but date-obs seems more reliable from FILENAME
      # CARACAL computes mjd-obs also from FILENAME
      self.dateobs = hdr['FILENAME'].replace("h",":").replace("m",":")
      self.dateobs = self.dateobs[4:8]+"-"+self.dateobs[8:10]+"-"+self.dateobs[10:21]
      self.mjd = hdr.get('HIERARCH CARACAL MJD-OBS')
      if not self.mjd:
         import warnings
         warnings.warn("Warning: keyword HIERARCH CARACAL MJD-OBS not found! This was implemented in CARACAL v2.00."+
                       "Please use lastest products.")
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', hdr.get(HIERARCH+'CARACAL DRIFT RV', np.nan))
      self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT RVERR', np.nan))
      self.fox = HIERARCH+'CARACAL FOX XWD' in hdr
      self.sn55 = hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 36', np.nan)
      sn25 = hdr.get(HIERARCH+'CARACAL FOX SNR 25', np.nan)
      sn30 = hdr.get(HIERARCH+'CARACAL FOX SNR 30', np.nan)
      if sn25 > 70 and sn25 > 10*sn30: # hig
         self.flag |= sflag.led
         #print(sn25, sn30, self.sn55, )

      self.fileid = hdr.get('FILENAME', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.timeid = self.fileid
      self.calmode = hdr.get('SOURCE', '') #.split("_")[3] #fileid[fileid.index('(')+1:fileid.index(')')]
      self.calmode = hdr.get(HIERARCH+'CARACAL FIB', '')
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

   if orders is not None:  # read order data
      f = hdulist['SPEC'].section[orders]
      w =  hdulist['WAVE'].section[orders]
      e = hdulist['SIG'].section[orders]
      bpmap = np.isnan(f).astype(np.uint64)            # flag 1 for nan

      bpmap0 = np.zeros((61,4096), dtype=np.uint64)
      bpmap0[14:38,[2453-3, 2453-2, 2453-1, 2453, 2453+1, 2453+2, 2453+3]] |= 1
      bpmap0[14:38,1643] |= 1   # ghost of hotspot tail
      bpmap0[14:38,2459] |= 1   # spikes of hotspot satellite (bug not correct due to bug in v2.00)
      bpmap0[15:41,3374] |= 1   # displaced column; ignore by marking as nan
      bpmap0[28,3395:3400] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[34,838:850] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis.fits
      bpmap0[34,2035:2044] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[34,3150:3161] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis.fits
      bpmap0[35,403:410] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,754:759] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[35,1083:1093] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,1944:1956] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,2710:2715] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[35,3050:3070] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,3706:3717] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,3706:3717] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[36,303:308] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,312:317] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,1311:1315] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,1325:1329] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[37,1326:1343] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[39,1076:1082] |= flag.sky # car-20170626T02h00m17s-sci-gtoc-vis
      bpmap0[39,1204:1212] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,1236:1243] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[39,1463:1468] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,2196:2203] |= flag.sky # car-20160520T03h10m13s-sci-gtoc-vis.fits
      bpmap0[39,2493:2504] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,3705:3717] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[40,2765:2773] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[40,3146:3153] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[40,3556:3564] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,486:491] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,495:501] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,1305:1315] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,480:490] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,1316:1330] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,2363:2368] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,2375:2382] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[44,3355:3361] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,311:321] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[46,835:845] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,1156:1171] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[46,1895:1905] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,2212:2232] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[47,2127:2133] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2218:2223] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2260:2266] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2313:2319] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3111:3116] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3267:3272] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3316:3321] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3432:3438] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[47,3480:3488] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3658:3665] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[49,1008:1017] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,2532:2544] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,3046:3056] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,3574:3588] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      # interpolate bad columns, they mess up a lot the creation of the template
      # We do this only when we read all order (preRVs), but not for coadding (single orders)
      if orders == np.s_[:]:
         # the hotspot
         #f[14:38,2453-3: 2453+4] = f[14:38,2453-4][:,np.newaxis] +  (f[14:38,2453+4]-f[14:38,2453-4]).reshape(24,1)*np.arange(1,8).reshape(1,7)/8.
         #f[15:41,3374] = np.nan # 0.5*(f[15:41,3374-1]+f[15:41,3374+1])
         pass


      bpmap |= bpmap0[orders]
      bpmap = bpmap.astype(int)

      with np.errstate(invalid='ignore'):
         # arrgh, in newer numpy version comparison with nan raises a warning
        if self.fox:
           e *= 10.
           f *= 10.
        else:
           e = np.sqrt(5.*10 + np.abs(f))
           bpmap[f>300000] |= flag.sat
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

      #e[bpmap > 0] = 1000. * np.max(e[bpmap == 0]) # set to zero for preRVs
      return w, f, e, bpmap

   if verb: print "read_carm_vis:", self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode


def read_carm_nir(self, s, orders=None, pfits=True, verb=True):
   """
   SYNTAX: read_carm_nir(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   HIERARCH = 'HIERARCH '
   hdulist = pyfits.open(s) # slow 30 ms
   if orders is None:
      hdr = hdulist[0].header
      self.header = hdr
      self.inst = hdr['INSTRUME'][0:4]+'_NIR'
      #data,hdr = fitsio.read(s,header=True)

      self.drsberv = hdr.get('HIERARCH CARACAL BERV', np.nan)
      # BJD for stars, JD for for calibration products
      self.drsbjd = hdr.get('HIERARCH CARACAL BJD', hdr.get('MJD-OBS',0.0))
      self.drsbjd = hdr.get('HIERARCH CARACAL BJD', np.nan) + 2400000
      self.dateobs = hdr['DATE-OBS']
      self.dateobs = hdr['FILENAME'].replace("h",":").replace("m",":")
      self.dateobs = self.dateobs[4:8]+"-"+self.dateobs[8:10]+"-"+self.dateobs[10:21]
      self.mjd = hdr.get('HIERARCH CARACAL MJD-OBS')
      if not self.mjd:
         import warnings
         warnings.warn("Warning: keyword HIERARCH CARACAL MJD-OBS not found! This was implemented in CARACAL v2.00."+
                       "Please use lastest products.")
      if isinstance(self.drsbjd, str): self.drsbjd = 0.0   # workaround for MJD-OBS bug @2016-Jan
      self.fox = HIERARCH+'CARACAL FOX XWD' in hdr
      # check LED for NIR This order can be affect
      r = hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan) /  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan)
      if r >1.5: print r, hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan),  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan)
      self.sn55 = min(hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan),  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan))
      if self.dateobs[:10] in ('2016-01-13', '2016-01-14', '2016-01-22', '2016-01-23'):
         self.sn55 = min(self.sn55, 10.) # bug in NIR fits file
      self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan)
      self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan)
      self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
      self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)

   #fileid = str(hdr.ascardlist()['MJD-OBS'])     # read the comment
      self.fileid = hdr.get('FILENAME', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.calmode = hdr.get(HIERARCH+'CARACAL FIB','')
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]
      self.timeid = self.fileid
      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.airmass = hdr.get('AIRMASS', 0.0)
      self.exptime = hdr['EXPTIME']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean==0: self.tmmean = 0.5

   if orders is not None:  # read order data
      bp = np.array([[1,230],[1,313],[1,633],[1,882],[1,1330],[1,1411],[1,1513],[1,1546],[1,1584],[1,1613],[1,1771],[1,1840],
                        [2,230],[2,1015],[2,1024],[2,1070],[2,1089],[2,1146],[2,1158],[2,1257],[2,2256],[2,1674],[2,1695],[2,1707],[2,1771],[2,1799],[2,1842],
                        [3,818],[3,941],[3,1183],[3,1323],[3,1464],[3,1598],[3,1692],[3,1762],
                        [4,231],[4,237],[4,329],[4,331],[4,231],[4,237],[4,329],[4,331],[4,947],[4,1222],[4,1248],[4,1280],[4,1385],[4,1613],[4,1675],[4,1813],
                        [5,319],[5,407],[5,442],[5,452],[5,749],[5,832],[5,972],[5,1140],[5,1514],[5,1525],[5,1547],[5,1624],[5,1650],
                        [6,927],[6,1100],[6,1215],[6,1270],[6,1299],[6,1343],[6,1483],[6,1507],[6,1545],[6,1738],
                        [7,281],[7,660],[7,661],[7,837],[7,969],[7,1080],[7,1143],[7,1290],[7,1426],[7,1570],[7,1741],
                        [8,237],[8,266],[8,528],[8,564],[8,682],[8,730],[8,800],[8,1063],[8,1195],[8,1229],[8,1548],[8,1592],[8,1610],[8,1660],[8,1716],[8,1758],[8,1826],
                        [9,313],[9,323],[9,496],[9,582],[9,618],[9,627],[9,775],[9,887],[9,1359],[9,1387],[9,1693],[9,1799],
                        [10,175],[10,202],[10,339],[10,441],[10,541],[10,575],[10,693],[10,758],[10,862],[10,871],[10,918],[10,1174],[10,1381],[10,1500],[10,1845],
                        [11,296],[11,352],[11,420],[11,745],[11,746],[11,1152],[11,1237],[11,1395],[11,1495],[11,1643],[11,1681],[11,1771],[11,1831],
                        [12,194],[12,245],[12,274],[12,301],[12,309],[12,314],[12,333],[12,335],[12,401],[12,418],[12,497],[12,665],[12,764],[12,815],[12,935],[12,994],[12,1097],[12,1176],[12,1217],[12,1390],[12,1624],[12,1873],
                        [13,201],[13,510],[13,734],[13,777],[13,844],[13,979],[13,1238],[13,1463],[13,1501],[13,1919],
                        [14,250],[14,282],[14,339],[14,415],[14,482],[14,527],[14,766],[14,852],[14,875],[14,926],[14,967],[14,1121],[14,1134],[14,1143],[14,1237],[14,1238],[14,1774],[14,1798],
                        [15,247],[15,284],[15,400],[15,691],[15,902],[15,951],[15,961],[15,1036],[15,1047],[15,1152],[15,1153],[15,1179],[15,1185],[15,1243],[15,1247],[15,1449],[15,1531],[15,1532],
                        [16,228],[16,282],[16,363],[16,324],[16,363],[16,409],[16,520],[16,763],[16,877],[16,968],[16,1309],[16,1467],[16,1516],[16,1855],
                        [17,212],[17,553],[17,691],[17,760],[17,945],[17,1020],[17,1067],[17,1350],[17,1359],[17,1659],[17,1450],[17,1453],[17,1478],[17,1589],[17,1696],[17,1796],[17,1719],[17,1760],[17,1807],[17,1855],[17,2011],
                        [18,338],[18,597],[18,792],[18,818],[18,942],[18,975],[18,1028],[18,1149],[18,1311],[18,1348],[18,1355],[18,1604],[18,1724],[18,1825],
                        [19,198],[19,209],[19,244],[19,273],[19,1242],[19,1372],[19,1382],[19,1418],[19,1462],[19,1550],[19,1580],
                        [20,450],[20,624],[20,655],[20,966],[20,1042],[20,1065],[20,1075],[20,1136],[20,1560],[20,1591],[20,1651],[20,1655],[20,1657],
                        [21,301],[21,437],[21,450],[21,545],[21,581],[21,642],[21,766],[21,775],[21,782],[21,806],[21,851],[21,853],[21,1067],[21,1097],[21,1125],[21,1152],[21,1249],[21,1398],[21,1403],[21,1459],[21,1581],[21,1620],[21,1823],[21,1865],[21,1924],
                        [22,383],[22,391],[22,407],[22,441],[22,644],[22,680],[22,843],[22,1007],[22,1021],[22,1022],[22,1023],[22,1031],[22,1290],[22,1944],[22,1983],[22,1996],
                        [23,304],[23,716],[23,742],[23,761],[23,821],[23,935],[23,945],[23,1029],[23,1266],[23,1560],[23,1586],
                        [24,356],[24,473],[24,502],[24,599],[24,688],[24,704],[24,744],[24,758],[24,767],[24,1019],[24,1268],[24,1270],[24,1576],[24,1639],[24,1751],[24,1860],
                        [25,190],[25,430],[25,827],[25,936],[25,938],[25,1375],[25,1508],[25,1583],[25,1670],[25,1681],[25,1714],[25,1724],[25,1926],[25,1927],[25,1974],
                        [26,802],[26,889],[26,1191],[26,1291],[26,1435],[26,1508],[26,1636],[26,1740],[26,1930],
                        [1,2733],[1,2654],[1,2672],[2,3644], #SCA1
                        [3,2765],[3,3543],[3,3544],[3,3624],
                        [4,2881],[4,2687],
                        [7,3346],
                        [8,2404],[8,2405],
                        [9,2813],[9,2814],
                        [10,2455],[10,3675],
                        [12,3152],
                        [13,2845],[13,3572],
                        [14,2590],[14,2652],[14,3130],
                        [15,2473],[15,3064],[15,3083],[15,3738],[15,3739],
                        [16,3131],[16,3231],[16,3130],[16,3499],[16,3626],
                        [17,2518],[17,2519],[17,2814],[17,3202],[17,3325],[17,3740],
                        [18,2383],[18,2514],[18,2739],[18,3390],[18,3568],[18,3744],
                        [19,3321],[19,3332],[19,3481],[19,3619],
                        [21,2715],[21,3110],[21,3250],[21,3673],
                        [23,3333],[23,3606],[23,3722],
                        [24,2722],[24,3205],[24,3310],[24,3552],
                        [25,2755],[25,3058],[25,3786],
                        [27,2440],[27,3616],[27,3792]
                        ]).T
      if orders == np.s_[:]:
         f = hdulist['SPEC'].data
         dim = f.shape
         if 1: # reshaping
            dim = (dim[0]*2,dim[1]/2)
            f = f.reshape(dim)
            bp[0], bp[1] = bp[0]*2, bp[1]%dim[1]

         w =  hdulist['WAVE'].data.reshape(dim)
         e = hdulist['SIG'].data.reshape(dim)
         bpmap = np.isnan(f).astype(int)            # flag 1 for nan
         # a hand made bad pixel list car-20160218T18h48m52s-sci-gtoc-nir.fits
         bpmap[bp[0],bp[1]] = 1
      # interpolate bad columns, they mess up a lot the creation of the template
      # We do this only when we read all order (preRVs), but not for coadding (single orders)
         if 1: # int
            g =1.0*f
            f[bp[0],bp[1]] = (f[bp[0],bp[1]-1]+f[bp[0],bp[1]+1]) / 2
            #o=17; gplot(g[o],'w l,',f[o],  'w lp')
            #pause()
      else:
         orders, det = divmod(orders,2)
         f = hdulist['SPEC'].section[orders]
         if 1:
            sl = [None, None]
            sl[1-det] = f.shape[0]/2
            sl = slice(*sl)
            f = f[sl]
            bp[0], bp[1] = bp[0]*2, bp[1]%f.shape[0]

         w =  hdulist['WAVE'].section[orders][sl]
         e = hdulist['SIG'].section[orders][sl]
         bpmap = np.isnan(f).astype(int)            # flag 1 for nan
         bp = bp[:,np.where(bp[0]==orders)[0]]
         if len(bp):
           bpmap[bp[1]] = 1
           f[bp[1]] = (f[bp[1]-1]+f[bp[1]+1]) / 2
         #pause()
      if self.fox:
         e *= 100000.
         f *= 100000.
#         e = 1000+0*f
      else:
         # fib B, linear extraction
         e = 0*f + np.sqrt(np.median(f)) # unweighted maybe for HCL
         bpmap[f>300000] |= flag.sat
         #bpmap[:,1:] |= flag.sat * (f>300000)[:,:-1]
         #bpmap[:,:-1] |= flag.sat * (f>300000)[:,1:]
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      return w, f, e, bpmap

   if verb: print "read_carm_nir:", self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode


def read_csfs_vis(self, s, orders=None, pfits=True, verb=True):
   """
   SYNTAX: read_carm_vis(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
   """
   HIERARCH = 'HIERARCH '
   hdulist = pyfits.open(s) # slow 30 ms
   if orders is None:
      hdr = hdulist[0].header
      self.header = hdr
      self.inst = hdr['INSTRUME'][0:4]+'_VIS'
   #data,hdr = fitsio.read(s,header=True)

      self.drsberv = -hdr.get('FIBARV',0)/1000.
      self.drsberv = -hdr.get(HIERARCH+'CAHA SIMU RVA',0)/1000.
      self.drsberv = hdr.get('E_BERV',0)  # @Ansgar simulation
      #self.drsbjd = hdr.get('MJD-OBS',0)
      self.drsberv = 0
      self.drsbjd = hdr.get('OBS_DAYS',0)
      self.drsbjd = hdr.get('E_BJD',0)

      self.sn55 = hdr.get(HIERARCH+'CARACAL FOX SNR 36', np.nan)
      self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
      self.drift = hdr.get(HIERARCH+'CARACAL DRIFT RV', np.nan)
      self.fox = HIERARCH+'CARACAL FOX XWD' in hdr

   #fileid = str(hdr.ascardlist()['MJD-OBS'])     # read the comment
      self.fileid = hdr.get('FILENAME',0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.calmode = hdr.get('SOURCE',0) #.split("_")[3] #fileid[fileid.index('(')+1:fileid.index(')')]
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]
      self.timeid = self.fileid
      self.exptime = hdr['EXPTIME']

   if orders is not None:  # read order data
      f = hdulist['SPEC'].section[orders]
      bpmap = np.isnan(f).astype(int)            # flag 1 for nan
      w =  hdulist['WAVE'].section[orders]
      e = hdulist['SIG'].section[orders]
      if self.fox:
         e *= 100000.
         f *= 100000.
      else:
         e = 1000+0*f
         bpmap[self.f>300000] |= flag.sat
         # flag neighbours
         bpmap[:,1:] |= flag.sat * (f>300000)[:,:-1]
         bpmap[:,:-1] |= flag.sat * (f>300000)[:,1:]
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      return w, f, e, bpmap

   if verb: print "read_carm_vis:", self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode


def read_feros(self, s, orders=None, pfits=True, verb=True):
   """
   SYNTAX: read_feros(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   hdulist = pyfits.open(s) # slow 30 ms
   if orders is None or self.header is None:
      HIERARCH = 'HIERARCH '
      hdr = hdulist[0].header
      self.inst = hdr['INSTRUME']

      #data,hdr = fitsio.read(s,header=True)
      indberv = [i for i, s in enumerate(hdr.values()) if 'BARY_CORR' in str(s)][0] + 1
      self.drsberv = float(hdr[indberv])
      # FEROS drsberv are wrong and already applied to the spectra
      # remove it from the spectrum
      self.drsbjd = hdr.get('MJD-OBS',np.nan) #hdr.get(HIERARCH+'ESO DRS BJD',0)
      if self.drsbjd: self.drsbjd += 2400000.5
      if self.fib=='B': self.drsberv = np.nan

      self.tmmean = 0.5 # hdr['HIERARCH ESO INS DET1 TMMEAN'] not available
      self.exptime = hdr['EXPTIME']
      self.dateobs = hdr['ARCFILE'][6:-5] # hdr['DATE-OBS']
      # FEROS has no S/N estimates; judge with airmass and seeing
      #print hdr['DATASUM'], self.exptime, hdr['HIERARCH ESO TEL AIRM START'], hdr['HIERARCH ESO TEL AMBI FWHM START']
      seefhwm = hdr['HIERARCH ESO TEL AMBI FWHM START']  # if -1 assume 3.0
      self.sn55 = 100 / hdr['HIERARCH ESO TEL AIRM START'] / (seefhwm if seefhwm>0 else 3.0)
      self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
      self.drift = hdr.get(HIERARCH+'ESO DRS DRIFT RV USED', np.nan)

      #fileid = str(hdr.ascardlist()['MJD-OBS'])     # read the comment
      fileid = hdr.get('ARCFILE',0) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.calmode = hdr.get('ORIGFILE',0).split("_")[3] #fileid[fileid.index('(')+1:fileid.index(')')]
      calmodedict = {'objcal': 'OBJ,ThAr', 'objsky': 'OBJ,SKY'}
      if self.calmode in calmodedict: self.calmode = calmodedict[self.calmode]
      self.timeid = self.dateobs #fileid

      if self.calmode != 'OBJ,ThAr' and hdr.get('FILENAME') != 'rebinned1.bdf': self.flag = 1
      #  hdr['OBJECT'] = hdr.get('OBJECT','FOX')
      self.header = hdr
   hdr = self.header
   if verb: print "read_feros:", self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode

   if orders is not None:  # read order data
      if self.filename.endswith('.mt'):
         realname = os.path.realpath(self.filename)
         data = [pyfits.getdata(realname.replace('01.mt','%02i.mt'%(o+1))) for o in range(39)]
         fsize = max([len(fo) for fo in data])
         f = np.zeros((39,fsize)).astype(float)
         for o,fo in enumerate(data): f[o,0:len(fo)] = fo* 100000.
      else:
         f = hdulist[0].data * 100000.
      bpmap = np.isnan(f).astype(int)  # flag 1 for nan
      bpmap[f < 0.0001] |= flag.neg    # flag 2 for zero and negative flux

      # print " applying wavelength solution ", file
      CRPIX1 = -49.
      CDELT1  = 0.03
      w = 0 * f.astype(float)
      his = list(hdr.get_history())
      wind = [i for i,x in enumerate(his) if 'WSTART' in x][0]
      wstart = ' '.join(his[wind+1:wind+14])
      wstart = np.fromstring(wstart, dtype=float, sep=' ')
      for i in range(39):
         #wstart = 3.527250000000000E+03
         w[i] = wstart[i] + CDELT1*np.arange(w[0].size)
      c = 299792.4580  # [km/s]
      w *= (1 - self.drsberv/c)
      w = airtovac(w)
         #from gplot import  *; gplot(lam[i],data[i])
      #pause()
         #lam = hdulist['WAVE'].data
      e = np.ones_like(f) * 100.
      e[bpmap==0] = 1.0*10. #np.sqrt(data[bpmap==0]) # slow 50 ms

      if self.fib!='B':  # filter for cosmics
         look = (None,)
         kap = 4.5
         for o in range(39):
            idx, = np.where(bpmap[o]==0)
            if 1:
               # using a robust 4.5 68 percentile clipping
               #hh = np.argsort(f[o][idx])
               #ii = hh[len(hh)*0.98:]  # reject 2%
               p1sig = (100-68.2) / 2
               p16, p50, p84 = np.percentile(f[o][idx],(p1sig,50,100-p1sig))
               rbsig = (p84-p16) / 2 # robust five sigma
               ii = idx[f[o][idx] > p50+kap*rbsig]
               if o in look:
                  gplot(w[o],f[o],',',w[o][idx],f[o][idx],',',w[o][ii],f[o][ii], ",%f, %f, %f, %f"%(p16, p50, p84, p50+kap*rbsig))
            if 0:
               # moving robust 4.5sigma clipping
               csum = np.cumsum(f[o][idx])
               mvmean = (csum[400:] - csum[:-400]) / 400
               #mvmean = np.lib.pad(mvmean, (400,400), 'edge')
               mvmean = np.concatenate((np.zeros(200) + mvmean[0], mvmean,np.zeros(200) + mvmean[-1]))
               csum = np.cumsum(abs(f[o][idx]-mvmean))
               mvrbsig = (csum[400:] - csum[:-400]) / 400
               mvrbsig = np.concatenate((np.zeros(200) + mvrbsig[0], mvrbsig,np.zeros(200) + mvrbsig[-1]))
               ii = idx[f[o][idx] > mvmean+kap*mvrbsig]
               if o in look:
                  gplot(w[o],f[o],',',w[o][idx],f[o][idx],mvmean,mvmean+kap*mvrbsig,',"" us 1:3, "" us 1:4,',w[o][ii],sf[o][ii])
            bpmap[o][ii] |= flag_cosm
            if o in look:
               pause(o)
      hdulist.close()
      return w, f, e, bpmap
   hdulist.close()


def read_fts(self,s, orders=None, filename=None, pfits=True, verb=True):
   """
   SYNTAX: read_fts(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   HIERARCH = 'HIERARCH '
   if orders is None:
      hdr = {'OBJECT': 'Iod'}
      self.header = hdr
      self.inst = 'FTS'
      self.drsberv = hdr.get('bla', 0)
      self.fileid = os.path.basename(s) #fileid[fileid.index('(')+1:fileid.index(')')]
      self.drsbjd = 0.0 if '.fits' in s else float(self.fileid.split('.')[1].split('_')[0])
      with open('/home/data1/fts/Lemke/2015_I/2015-01-29/DPT/parameter_Q_Si_I2_001.txt') as finfo:
         for line in finfo:
            if self.fileid.replace('_ScSm.txt','\t') in line:
               line = line.split(); date=line[2].split('/'); time=line[3].split(':')
               #from subprocess import call
               #call(['bash','date2jd', date[2], date[1], date[0]] +time)
               from subprocess import Popen, PIPE
               p = Popen(['bash','date2jd', date[2], date[1], date[0]] +time, stdin=PIPE, stdout=PIPE, stderr=PIPE)
               output, err = p.communicate()
               rc = p.returncode
               self.drsbjd = float(output)
      self.sn55 = 10
      self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
      self.drift = hdr.get(HIERARCH+'ESO DRS DRIFT RV USED', np.nan)
      self.calmode = hdr.get('SOURCE', 0) #.split("_")[3] #fileid[fileid.index('(')+1:fileid.index(')')]
      self.timeid = self.fileid
      self.exptime = 0
   if verb: print "read_fts:", self.timeid, self.header['OBJECT'], self.drsbjd, self.sn55, self.drsberv, self.drift, self.flag, self.calmode

   if orders is not None:  # read order data
      nord = 70    # some arbitary shaping to 70x10000
      nw = 700000
      if '.fits' in s:
         hdulist = pyfits.open(s)
         w, f = hdulist[0].data[:,1600000:]
      else:
         #w, f = np.loadtxt(s, skiprows=1600000, unpack=True) ; too slow 45s vs 5.2s
         data = np.fromfile(s, sep=' ', count=2*(1600000+nw))[2*1600000:]
         w = 1e7 / data[::2]
         f = data[1::2]

      w = 10 * w[:nw].reshape(nord,nw/nord)
      f = f[:nw].reshape(nord,nw/nord)*4000
      bpmap = np.isnan(f).astype(int)            # flag 1 for nan

      xmax = w.size
      e = np.ones_like(w)
      return w, f, e, bpmap


tarmode = 5
# 0  - extract physically (works for all, slow)
#      tarfits directory not cleaned
# 1  - treat tarfile as normal file which has offset, seek() and open()
#      but no size limit function, therefore not supported by pyfits
# 2  - as 1, but with try to overload build-in open (not working for both)
# 3  - tar.extractfile (use its build-in functions for read, offset and size, slower read!?) not supported by pyfits
# 5  - use tar.extractfile, works with pyfits, doesn't work with myfits (np.fromfile accepts only real files and files object, not tar)
# 6  - open as normal file, does not work

def file_from_tar(s, inst='HARPS', fib=None, **kwargs):
   """
   Returns a file-like object.
   Reading header and data should use the same fits reader

   >>> s = '/home/astro115/carmenes/data/HARPS/DRS/gj876/ADP.2014-12-07T19:00:30.953.tar'
   >>> tar = tarfile.open(s)
   >>> extr = tar.getmember('data/reduced/2014-12-06/HARPS.2014-12-07T00:35:29.894_e2ds_A.fits')

   """
   pat = {'HARPS': {'A': '_e2ds_A.fits', 'B': '_e2ds_B.fits'}[fib],
          'FEROS': {'A': '.1061.fits', 'B': '.1062.fits'}[fib]} [inst]
   tar = tarfile.open(s)
   for member in tar.getmembers():
       if pat in member.name: extr = member

   if tarmode in (1,2,3) and kwargs.get('pfits') == 2:
      # We could use tar.extractfile(extr) but this requires reopen the tar file
      # and does not support "with open()".
      # Instead we will use offset, seek() and open().
      # but both seems to be slower than physical extraction
      extr.mother = s
      s = extr
   elif tarmode in (0,5) and kwargs.get('pfits') == 2:
      # Open the tar file as a normal file and store the position and size of the fits file
      # as attributes.
      s = type('tarobj', (file,), {})(s, mode='rb')   # since class, because the in the buildin file class we cannot set new attributes
      s.mem = extr.name
      s.offset_data = extr.offset_data
      s.size = extr.size
      tar.close()
   elif tarmode == 5:
      # does not work with np.fromfile
      s = tar.extractfile(extr)
   else:
      print 'extract'
      tar.extract(extr, path='tarfits')   # extract physically
      s = 'tarfits/'+extr.name
   #tar.close()
   return s


class imhead(dict):
   """
   Returns fitsheader as a dict.

   Keyword arguments:
   filename -- filename of the fitsfile

   Returns
   -------
   imhead: dict of values and dict of comments
   EXTHDRSZ: size of the header [bytes]

   Examples
   --------
   >>> x = imhead(filename)
   >>> x['OBJECT']
   >>> x.comments['OBJECT']

   """
   def __init__(self, s, *args, **kwargs):
      hdr = {}
      self.comments = {}
      extpos = kwargs.get('extpos', 0)  # position of the extension with the fitsfile
      count = kwargs.get('count', -1)
      #if not len(args):
         #args = ''   # read all
         #stop( '\n\nimhead: there must be are args; reading all keys not yet supported\n')
      args += ('NAXIS', 'BITPIX')
      NR = 0

      if isinstance(s, tarfile.TarInfo):
         if tarmode==1:
            fi = open(s.mother)
         if tarmode==3:
            fi = tarfile.open(s.mother).extractfile(s)
      elif isinstance(s, (tarfile.ExFileObject, file)):   # tarmode == 5,6
         fi = s
      elif s.endswith('.gz'):   # normal file
         fi = gzip.open(s)
      else:   # normal file
         fi = open(s)

      if not extpos and hasattr(s, 'offset_data'):
         extpos = s.offset_data
      #pause()

      #with open(s) as fi:
      if 1:
         fi.seek(extpos)
         for card in iter(lambda:fi.read(80), ''):   # read in 80 byte blocks
            NR += 1
            if card.startswith('END '): break
            if card.startswith(args):
               #key, val, comment = card.replace("= ","/ ",1).split("/ ")
               key, val, comment = card.replace(' /','= ',1).split('= ')
               hdr[key.strip()] = val.strip("' ") if "'" in val else float(val) if '.' in val else int(val)
               #hdr[key.strip()] = val.strip("' ") if any(i in val for i in "'TF") else float(val) if '.' in val and not 'NaN.' in val else int(val) if val.strip(" -").isdigit() else val
               self.comments[key.strip()] = comment
               count -= 1
               if count==0: args = () # all found; do not check cards anymore; only read to end
         #NR = (fi.tell()-extpos) / 80

      hsz = 2880 * ((NR-1)/36 + 1)
      dsz = 0     # EXTDATSZ
      self.NAXIS = hdr.get('NAXIS', 0)
      if self.NAXIS:
         self.BITPIX = hdr['BITPIX']
         self.NAXIS1 = hdr['NAXIS1']
         dsz = abs(self.BITPIX)   # data size
         dsz *= self.NAXIS1
         if self.NAXIS > 1:
            self.NAXIS2 = hdr['NAXIS2']
            dsz *= self.NAXIS2
         dsz = ((dsz/8-1)/2880+1) * 2880

      self.EXTHDRSZ = hsz   # add hdr size
      self.EXTDATA = extpos + hsz
      self.EXTEND = extpos + hsz + dsz
      super(imhead, self).__init__(hdr)


class getext(dict):
   """
   Fast reading of single orders.

   Uses imhead to scan the extension headers and to get
   the data types and sizes.

   Examples
   --------
   >>> hdu = getext(filename)
   >>> hdu['SPEC']
   >>> hdu.getdata('SPEC', o=13)

   """
   def __init__(self, s):
      """
      Extracts from header for each extension the data position, type, and shape
      and saves the extension information as dictionary.

      """
      self.fileobj = s
      ext = {}
      filepos = 0
      fileend = -1
      i = 0
      if isinstance(s, tarfile.TarInfo) and s.tarmode:
         self.fileobj = s.mother
      if hasattr(s, 'offset_data'):
         filepos = s.offset_data
         fileend = s.offset_data + s.size
      #pause()

      while True:
         exthdr = imhead(self.fileobj, 'EXTNAME', extpos=filepos)
         if exthdr.EXTEND <= filepos: break   # no new records
         filepos = exthdr.EXTEND
         ext[i] = ext[exthdr.get('EXTNAME', i)] = exthdr
         if filepos == fileend: break   # end of filemember in tarfile
         i += 1
      super(getext, self).__init__(ext)

   def getdata(self, extname=None, o=np.s_[:]):
      if extname is None:  # search first extension with data
         extname = 0
         while not self[extname].NAXIS: extname += 1
      ext = self[extname]
      dtype = {-32: '>f4', -64: '>f8'}[ext.BITPIX]
      dsize = {-32: 4, -64: 8}[ext.BITPIX]
      was_open = True
      if isinstance(self.fileobj, (tarfile.ExFileObject, file)):
         funit = self.fileobj
      else:
         funit = open(self.fileobj)
         was_open = False

#      with open(self.fileobj) as funit:
      if 1:
         if o == np.s_[:]: # read all
            funit.seek(ext.EXTDATA)
            data = np.fromfile(funit, dtype=dtype, count=ext.NAXIS1*ext.NAXIS2).reshape((ext.NAXIS2,ext.NAXIS1))
         else:
            funit.seek(ext.EXTDATA+o*ext.NAXIS1*dsize)
            data = np.fromfile(funit, dtype=dtype, count=ext.NAXIS1)
      if not was_open:
         funit.close()
      return data


def bary(obj,bjd, exptime):
   from subprocess import Popen, PIPE, STDOUT
   #bjd = 2450000.5+ obj['MJD-OBS']
   #bjd = obj['HIERARCH ESO OBS START']
   #obj = obj['OBJECT']
   #pause(bjd)
   YY, MM, DD, hh, mm, ss = bjd.replace('T', ' ').replace(':', ' ').replace('-', ' ').split(' ')
   # check alias names
   with open('./src/bary/star_alias.txt','r') as f:
      for line in f:
         if obj in line:   obj=line.split()[0]
   # prepare input file
   s = ["Select observatory (LASILLA / PARANAL):",
        "LASILLA",
        "number of observations:",
        "1",
        "coord. (0) or object(1)?",
        "1",
        "Objectname in list:",
        obj,
        "Epoch:",
        "2000.0",
        "R.A. (hh mm ss)",
        "00 00 00",
        "Dec. (deg mm ss)",
        "00 00 00",
        "Date (dd mm yyyy)",
        DD+" "+MM+" "+YY, #"13 08 2000",
        "Time (hh mm ss)",
        hh+" "+mm+" "+ss, #"07 35 49",
        "Exptime (sec)",
        str(int(exptime))]
   #pause(exptime)
   p = Popen(['./src/bary/eph'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
   eph = p.communicate(input="\n".join(s)+"\n")[0]
   if  'Star not in starlist' in eph: stop(obj, 'Star not in starlist')
   berv = float(eph.split("\n")[-3].lstrip(" DELTAV=").rstrip("m/s"))
   return berv

def airtovac(wave_air):
   """
taken from idl astrolib
;+
; NAME:
;       AIRTOVAC
; PURPOSE:
;       Convert air wavelengths to vacuum wavelengths 
; EXPLANATION:
;       Wavelengths are corrected for the index of refraction of air under 
;       standard conditions.  Wavelength values below 2000 A will not be 
;       altered.  Uses relation of Ciddor (1996).
;
; CALLING SEQUENCE:
;       AIRTOVAC, WAVE_AIR, [ WAVE_VAC]
;
; INPUT/OUTPUT:
;       WAVE_AIR - Wavelength in Angstroms, scalar or vector
;               If this is the only parameter supplied, it will be updated on
;               output to contain double precision vacuum wavelength(s). 
; OPTIONAL OUTPUT:
;        WAVE_VAC - Vacuum wavelength in Angstroms, same number of elements as
;                 WAVE_AIR, double precision
;
; EXAMPLE:
;       If the air wavelength is  W = 6056.125 (a Krypton line), then 
;       AIRTOVAC, W yields an vacuum wavelength of W = 6057.8019
;
; METHOD:
;	Formula from Ciddor 1996, Applied Optics 62, 958
;
; NOTES:
;       Take care within 1 A of 2000 A.   Wavelengths below 2000 A *in air* are
;       not altered.
; REVISION HISTORY
;       Written W. Landsman                November 1991
;       Use Ciddor (1996) formula for better accuracy in the infrared 
;           Added optional output vector, W Landsman Mar 2011
;       Iterate for better precision W.L./D. Schlegel  Mar 2011
;-
   """

   wave_vac = wave_air * 1.0
   g = wave_vac > 2000     #Only modify above 2000 A

   if np.sum(g):
      for iter in [0, 1]:
         if isinstance(g, np.ndarray):
            sigma2 = (1e4/wave_vac[g])**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac[g] = wave_air[g] * fact              #Convert Wavelength
         else: # scalar version
            sigma2 = (1e4/wave_vac)**2.     #Convert to wavenumber squared
            # Compute conversion factor
            fact = 1. + 5.792105e-2 / (238.0185 - sigma2) + \
                               1.67917e-3 / (57.362 - sigma2)
            wave_vac = wave_air * fact              #Convert Wavelength

   return wave_vac


if __name__ == "__main__":
   if not 'debug' in sys.argv:
      x=Spectrum(*sys.argv[1:], inst='HARPS', pfits=2)
      x.read_data()
   else:
      sys.argv.remove('debug')
      try:
         read_spec(*sys.argv[1:])
      except:
         import pdb, sys
         e, m, tb = sys.exc_info()
         pdb.post_mortem(tb)


