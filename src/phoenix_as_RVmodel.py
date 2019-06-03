#! /usr/bin/python

try:
   import astropy.io.fits as pyfits
except:
   import pyfits
import sys
import os

wmax=8000
wmin=3500

def readphoenix(fitsphoe, wmin=3500, wmax=8000):
   #fitsphoe = '/home/astro56/husser/PhoenixModels/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte03900-4.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
   #fitsphoe = '/home/astro115/carmenes/serval/phoenix_tmp/lte03400-5.00-0.0_carm.fits'
   #fitsphoe = '/home/astro115/carmenes/serval/phoenix_tmp/lte03400-5.00-0.0_carm.fits'
   print 'readphoenix', fitsphoe
   hdulist = pyfits.open(fitsphoe)
   flux = hdulist[0].data
   #hdulist = pyfits.open('/home/astro56/husser/PhoenixModels/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
   hdulist = pyfits.open(os.path.dirname(fitsphoe)+os.sep+'/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
   wave = hdulist[0].data
   hdulist.close()
   #if flux.shape!=wave.shape:
      #flux = flux.T
   imap = (wmin<wave) * (wave<wmax)
   w = wave[imap]
   f = flux[imap]
   return w, f

def phoe_vac():
   from read_harps import airtovac
   hdulist = pyfits.open('/home/astro56/husser/PhoenixModels/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
   wavevac = airtovac(hdulist[0].data)
   hdu = pyfits.PrimaryHDU(wavevac)
   hdu.writeto('/home/astro43/zechmeister/carmenes/serval/lib/WAVE_VAC_PHOENIX-ACES-AGSS-COND-2011.fits')
