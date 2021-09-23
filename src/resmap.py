from __future__ import print_function
import argparse
import glob
from astropy.io import fits
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
argopt = parser.add_argument   # function short cut
argopt('tag', help='Tag, output directory and file prefix (e.g. Object name).')
argopt('-o', help='Orders.', action='store_true')
argopt('-P',  help='Sort phase folding with by period in RV', action='store_true')
argopt('-a', help='Sort by airmass', action='store_true')
argopt('-berv', help='Sort by airmass in BERV', action='store_true')
argopt('-frame', nargs='?', help='ccd, earth, bary, star')
argopt('-ext', nargs='?', help='extension', default='res')
argopt('-data', help='Compute res+fmod.', action='store_true')
argopt('-ratio', help='Compute res/fmod.', action='store_true')
argopt('-sort', nargs='?', help='Criteria to sort spectra (default: filename)')


args = parser.parse_args()

tag = args.tag.rstrip('/')
ext = args.ext

print(tag+'/res/*.fits')
filenames = np.array(sorted(glob.glob(tag+'/res/*.fits')))
N = len(filenames)

r = []
w = []
o = 25
airmass = []
berv = []
snr = []
# HIERARCH SERVAL RVC = 270.0132768714515 / [m/s] RV drift corrected              

for f in filenames:
    print(f) #hdu.info())
    hdu = fits.open(f)
    airmass.append(hdu[0].header['airmass'])
    berv.append(hdu[0].header['HIERARCH SERVAL BERV'])
    snr.append(hdu[0].header['HIERARCH CARACAL FOX SNR 55'])
    r += [[]]
    if f == filenames[0]:
        w0 = hdu['wave'].data
    for o in range(25, 43):
#        print(o)
        d = hdu[ext].data[o]
        if args.ratio:
            d = 1 + d / hdu['fmod'].data[o]
        if args.data:
            d = d + hdu['fmod'].data[o]
            d /= np.nanmedian(d)
        if args.frame:
            z = 0
            if args.frame == 'bary':
                 z = berv[-1]/3e5
            w = hdu['wave'].data[o]
            d = np.interp(w0[o], w+z, d)
        # nothing is done for fram ccd 
        r[-1] += [d]


print('resmap.fits')
r = np.array(list(map(np.array, r)))

if args.a:
    ii = np.argsort(airmass)
    r = r[ii]
    print(zip(filenames[ii], np.array(airmass)[ii]))
if args.sort == 'berv' or args.berv:
    ii = np.argsort(berv)
    r = r[ii]
    print(zip(filenames[ii], np.array(airmass)[ii]))
if args.sort == 'snr':
    ii = np.argsort(snr)
    r = r[ii]

fits.writeto('resmap.fits', np.swapaxes(r,0,1), overwrite=True)

from ds9 import *
ds9('resmap.fits', port='resmap')
import time
#time.sleep(3)
from pause import *
pause()
ds9.text(0*np.arange(N), np.arange(N), list(filenames[ii]), port='resmap')


