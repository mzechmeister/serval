#! /usr/bin/python
from __future__ import division, print_function

# To merge echelle templates from serval into a s1d spectrum.
# Requires deblazed spectra.
#
# usage: ./merge_ech.py J20450+444/J20450+444.tpl.fits
#
# The wavelengths are stored with the three header keywords

import numpy as np
from astropy.io import fits
import sys

from pause import *

tplname = sys.argv[1]

if 1:
    hdu = fits.open(tplname)
    f0 = hdu['SPEC'].data
    w0 = hdu['WAVE'].data   # these are logarithmic wavelengths!
    ok = w0 > 0
    w = w0[ok]
    f = f0[ok]
    ii = np.argsort(w)
    w = w[ii]
    f = f[ii]

import cspline as spl
nx = len(w0[0])
x = np.linspace(-1, 1, nx)
# downweight to zero towards the egde with a parabola (to reduce jumps)
we = np.array([1-x**2]*len(w0))[ok][ii]

# Number of knots (reduced for the fit by a factor of 2 to minimise oscillations)
K = int((w[-1]-w[0])/ (w[1]-w[0]) /2)

 
pspllam = None
uspl = spl.ucbspl_fit(w, f, we, K=K, lam=pspllam)
wk, fk = uspl.osamp(2)

# gplot(w[:50000], f[:50000],  'w l ,', wk[:50000],fk[:50000], 'lc 3')

filename = tplname.split('/')[-1].replace(".fits","") + ".s1d.fits"
hdr = hdu[0].header
hdr['CDELT1'] = wk[1] - wk[0]
hdr['CRVAL1'] = wk[0]
hdr['CRPIX1'] = 1.
sp = fits.PrimaryHDU(header=hdr, data=fk)
fits.HDUList([sp]).writeto(filename, overwrite=True)
print("=>", filename)

#pause()

