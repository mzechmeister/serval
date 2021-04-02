#! /usr/bin/env python
from __future__ import division

# merge orders of a template
# only for blaze-normalised spectra
# sampling currently hardcoded

from astropy.io import fits
from gplot import *
import numpy as np
import cspline

def tpl1d(filename, out=None, show=False):
    hdu = fits.open(filename)
    hdr = hdu[0].header

    w, f = hdu['wave'].data.astype(float), hdu['spec'].data.astype(float)
    ok = w != 1.

    # parabolic window to downweight of edges to reduce egde effects
    nx = ok[0].size
    x = np.arange(nx)/nx - 0.5
    y = 1 - x**2/0.25
    ymap = ok * y
    spl = cspline.ucbspl_fit(w[ok], f[ok], w=ymap[ok], K=int((w[ok].max()-w[ok].min())*3e5/0.2), lam=0.1)

    if show:
        #gplot(w[ok], f[ok], 'w l, ', spl.osamp(1), 'w l')
        gplot(np.exp(w[ok]), f[ok], 'w l t "2d", ', np.exp(spl.osamp(1)[0]), spl.osamp(1)[1], 'w l t "1d"')
        from pause import pause; pause()
    tpl = np.rec.array(spl.osamp(1), names='lnwave,flux')

    if not out:
        out = filename.replace('.fits','') + '_1d.fits'

    fits.writeto(out, tpl, hdr[20:], overwrite=True)

    print("-> "+out)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    parser.add_argument('out', nargs='?')
    parser.add_argument('-show', action='store_true')
    args = parser.parse_args()
    tpl1d(**vars(args))

# ./tpl1d.py /home/astro115/carmenes/data/svn/serval/CARM_VIS/J12479+097/template.fits J12479+097_tpl.fits
