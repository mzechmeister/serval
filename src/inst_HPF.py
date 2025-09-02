#! /usr/bin/env python

# imports
from read_spec import *
from astropy.time import Time
import numpy as np
import cspline as spl

# Instrument parameters
name = 'HPF'
obsname = "hpf" # for barycorrpy
       # wikipedia
obsloc = dict(lat=30.681444,    # 30d40'53" N
              lon=-104.014722,  # 104d00'53" W
              elevation=2026)
R = 53000. # resolving power

pat = '*.fits'
pmax = 2048 - 300
iomax = 28
oset = "[4,5,6,14,15,16,17,18,26]"
coset = "[4,5,6,14,15,16,17,18,26]"

maskfile = 'telluric_mask_carm_short.dat'
skyfile = 'sky_mask_5_sigma.dat'
blazefile = 'hpf_blaze_spectra.fits'  # https://github.com/grzeimann/Goldilocks_Documentation/blob/master/hpf_blaze_spectra.fits


# Instrument read functions
def scan(self, s, orders=None, pfits=True, verb=True):
   """
   Returns
   -------
   namedtuple('spectrum', 'w f berv bjd drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           drift - Used RV Drift
           sn55  - S_N order center55
   Example
   -------
   >>> scan(filename)

   """
   HIERARCH = 'HIERARCH '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   self.header = hdr = hdulist[0].header
   self.instname = hdr['INSTRUME']
   self.drsberv = hdr.get('BERV', np.nan)
   self.drsbjd = hdr.get('BJD', np.nan) + 2400000
   self.dateobs = hdr['DATE-OBS']
   self.mjd = Time(self.dateobs, format='isot', scale='utc').mjd
   self.sn55 = hdr.get('SNR', 50)
   self.fileid = hdr.get('DATE-OBS', 0) 
   self.timeid = self.fileid

   self.calmode = "%s,%s,%s" % (hdr.get('SCI-OBJ', ''), hdr.get('CAL-OBJ', ''), hdr.get('SKY-OBJ', ''))

   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr.get('AIRMASS', np.nan)
   self.exptime = hdr['ITIME']
   self.tmmean = 0.5

   if 'Goldilocks' in self.filename:
      self.drift = hdr.get(HIERARCH+'LRVCORR', np.nan)
      # no error given (https://github.com/grzeimann/Goldilocks_Documentation)
      self.e_drift = 0
   
   elif 1:
      # for HPF spectra the drift is already included in the wavelength solution
      self.drift = 0
      self.e_drift = 0

   if self.mjd > 59731:
       # add estimate of drift offset post 59731 downtime
      self.drift = self.drift - 60

def data(self, orders, pfits=True):
   # function handles HPF spectra reduced by the Goldilocks and PSU pipelines
   if 'Goldilocks' in self.filename:
      return data_goldilocks(self, orders, pfits=True)
   elif 1:
      return data_psu(self, orders, pfits=True)
      
def data_psu(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      # science
      w,f,e = hdulist['Sci Wavl'].data, hdulist['Sci Flux'].data, hdulist['Sci Variance'].data
      w,f,e = w[orders],f[orders],e[orders]
      e = e**0.5

      # sky
      wsky,sky,esky = hdulist['Sky Wavl'].data, hdulist['Sky Flux'].data, hdulist['Sky Variance'].data
      wsky,sky,esky = wsky[orders],sky[orders],esky[orders]
      esky = esky**0.5            

      # f, e = deblaze(f,e,orders) # deblaze spectra -> already is deblazed
      #sky, esky = deblaze(sky,esky,orders,channel=0) # deblaze sky -> already is deblazed
      #f, e = skysub_orders([w,f,e,wsky,sky,esky]) # sky subtraction disabled

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan   

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      bpmap[np.isnan(e)] |= flag.nan

      return w, f, e, bpmap

def data_goldilocks(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      # science
      w,f,e = hdulist['Sci Wavl'].data, hdulist['Sci Flux'].data, hdulist['Sci Error'].data
      w,f,e = w[orders],f[orders],e[orders]      

      # sky
      wsky,sky,esky = hdulist['Sky Wavl'].data, hdulist['Sky Flux'].data, hdulist['Sky Error'].data
      wsky,sky,esky = wsky[orders],sky[orders],esky[orders]

      f, e = deblaze(f,e,orders) # deblaze spectra      
      sky, esky = deblaze(sky,esky,orders,channel=0) # deblaze sky      
      #f, e = skysub_orders([w,f,e,wsky,sky,esky]) # sky subtraction disabled

   bpmap = np.isnan(f).astype(int)            # flag 1 for nan

   with np.errstate(invalid='ignore'):
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      bpmap[np.isnan(e)] |= flag.nan
     
      return w, f, e, bpmap

# Sky subtraction function and wrapper
def polyfit(x,y,w,deg):
    '''
    Wrapper for polyfit, so it can deal with nans.
    '''
    # mask bad pixel
    m1 = ~np.any(np.isnan([x,y,w]),axis=0)
    m2 = ~np.any(np.isinf([x,y,w]),axis=0)
    m = np.logical_and(m1,m2)
    m3 = ~np.any(np.array([x,y,w])==0,axis=0)
    m = np.logical_and(m,m3)
    return np.polyfit(x[m],y[m],deg,w=w[m])


def splinefit(x,y,w):
    '''
    Wrapper for splinefit, so it can deal with nans.
    '''
    # mask bad pixel
    m1 = ~np.any(np.isnan([x,y,w]),axis=0)
    m2 = ~np.any(np.isinf([x,y,w]),axis=0)
    m = np.logical_and(m1,m2)
    m3 = ~np.any(np.array([x,y,w])==0,axis=0)
    m = np.logical_and(m,m3)
    return spl.ucbspl_fit(x[m],y[m],w[m],K=len(x),lam=1e-2)


def skysub_wfit(wsci,sci,esci,wsky,sky,esky):
    '''
    Sky subtraction with fit of sky emission
    
    Performs cubic spline interpolation of sky emission.
    Fits scaling of the sky channel to the sci channel.
    Subtracts scaled sky channel from sci channel. 
    
    Returns sky subtracted flux and flux error.
    '''

    # spline interpolation of sky emission to science wavelength
    # ignore overread
    splfit = splinefit(wsky[10:-10],
                       sky[10:-10],
                       esky[10:-10]**-2)
    sky_fct = splfit.to_spl()
    sky_int = sky_fct(wsci)
    # keep uncertainty
    esky_int = esky

    # polyfit sci continuum
    z = polyfit(wsci[10:-10],
                (sci-sky_int)[10:-10],
                (esci**2+esky**2)[10:-10]**-1,
                3)
    p = np.poly1d(z)
    scicont = p(wsci)

    # polyfit scaling (deg=0 -> simple scaling factor)
    z2 = polyfit(wsci[10:-10],
                 ((sci-scicont)/sky_int)[10:-10],
                 ((esci/sky_int)**2+(esky_int*scicont*sky_int**-2)**2)[10:-10]**-1,
                 0)
    scal = np.poly1d(z2)
    
    # subtract scaled sky emission
    flux = sci-scal(wsci)*sky_int
    eflux = (esci**2+(scal(wsci)*esky_int)**2)**0.5
 
    if 0:
       import matplotlib.pylab as plt
       plt.figure(figsize=(12,8))
       plt.plot(wsci[10:-10],(sci-sky_int)[10:-10],c='k',marker='.',lw=0)
       x = np.linspace(np.min(wsci),np.max(wsci),10000)
       plt.plot(x,p(x),c='r',lw=1,alpha=1)
       plt.ylabel('Flux [abtr.]')
       plt.xlabel('Wavelength [angstrom]')
       plt.show()

       plt.figure(figsize=(12,8))
       plt.plot(wsci[10:-10],((sci-scicont)/sky_int)[10:-10],c='k',marker='.',lw=0)
       x = np.linspace(np.min(wsci),np.max(wsci),10000)
       plt.plot(x,scal(x),c='r',lw=1,alpha=1)
       plt.ylabel('Flux [abtr.]')
       plt.xlabel('Wavelength [angstrom]')
       plt.show()


       plt.figure(figsize=(12,8))
       plt.plot(wsci,scicont,c='orange',marker='.',alpha=0.5,label='sci cont')
       plt.plot(wsci,sci,c='k',marker='.',alpha=0.5,label='sciflux')
       plt.plot(wsci,sky_int,c='gray',marker='.',alpha=0.5,label='skyflux (interpolated)')
       plt.plot(wsci,scal(wsci)*sky_int+scicont,c='b',marker='.',alpha=0.2,label='model skyflux+starflux')
       plt.plot(wsci,sci-sky,c='green',marker='.',alpha=0.5,label='sciflux-skyflux')
       plt.plot(wsci,flux,c='r',marker='.',alpha=0.5,label='model starflux')

       plt.legend()

       plt.ylabel('Flux [abtr.]')
       plt.xlabel('Wavelength [angstrom]')
       plt.show()

    return np.array([flux,eflux],dtype=np.float64)

def skysub(wsci,sci,esci,wsky,sky,esky):
    '''
    Sky subtraction
    
    Performs cubic spline interpolation of sky emission.
    Returns sky subtracted flux and flux error.
    '''

    # spline interpolation of sky emission to science wavelength
    # ignore overread
    splfit = splinefit(wsky[10:-10],
                       sky[10:-10],
                       esky[10:-10]**-2)
    sky_fct = splfit.to_spl()
    sky_int = sky_fct(wsci)
    # keep uncertainty
    esky_int = esky
    
    # subtract interpolated sky emission
    flux = sci-sky_int
    eflux = (esci**2+esky_int**2)**0.5
    
    return np.array([flux,eflux],dtype=np.float64)

    
def skysub_orders(dat):
    '''
    Wrapper to perform sky subtraction order wise.
    '''

    # prepare data for sky subtraction
    if len(np.shape(dat)) == 2:
        # special case: single order 

        # sky subtraction
        f,e = skysub_wfit(*tuple(np.array(dat)))
        return f,e
    else:
        dat = np.array(dat).swapaxes(0, 1)

        # sky subtraction: loop over orders
        f,e = np.array([skysub_wfit(*tuple(d)) for d in dat]).swapaxes(0, 1)
        return f,e

# Deblazing function
def deblaze(f,e,orders,channel=1):
    '''
    Perform deblazing using the blazefile.
    '''

    # load blaze function from lib
    path = os.path.dirname(os.path.abspath(__file__)).replace('src','lib') + os.sep
    b = pyfits.getdata(path + blazefile)[int(2-channel)::3][::-1][orders]

    # normalize blaze function
    b = np.array([bb/np.nanmedian(bb) for bb in b]) 

    # devide flux by blaze function
    f = f/b
    e = e/b
    
    return f,e
