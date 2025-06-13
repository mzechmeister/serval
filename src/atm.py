import numpy as np
from astropy.io import fits

airmass_to_c1 = lambda x: x

def load(filename):
    global tpl1, tpl2, tplo1, tplo2, is_echelle, airmass_to_c1

    # velocity offset correction
    v1 = 1.
    if filename.endswith('stdatm.fits'):
        v1 = 1.000001

    # load atm data
    hdu = fits.open(filename)

    if filename.endswith('stdAtmos_vis.fits'):
        # from lib/atmos/stdAtmos_vis.fits
        d = hdu[1].data
        w = d['lambda']
        f1 = d.O2
        f2 = d.H2O

        # setup a model
        ok = np.isfinite(f1)
        tplo1 = tpl1 = v1*w[ok], 1.*f1[ok]

        ok = np.isfinite(f2)
        tplo2 = tpl2 = 1.*w[ok], 1.*f2[ok]
        is_echelle = False
    else:
        w = hdu['WAVE'].data
        f1 = hdu['O2'].data
        f2 = hdu['H2O'].data
        c0 = hdu[0].header['C1_REF']      # linear relation to transform airmass to coefficients
        c1 = hdu[0].header['C1_SCALE']
        airmass_to_c1 = np.poly1d((1.1*c1, c0))
        # setup a model
        tpl1 = v1*w, 1.*f1
        tpl2 = 1.*w, 1.*f2
        is_echelle = True

    #f1 = f1 / np.nanpercentile(f1[44], 90)
    #f2 = f2 / np.nanpercentile(f2[44], 90)


def fit_atm_par(u, f, a1=None, o=None):
    '''guess atm parameters

    a1: value for a fixed relative air mass
    '''
    
    ok = np.isfinite(f)
    w2 = np.exp(u[ok])

    if is_echelle:
        atm1 = np.interp(w2, tpl1[0][o], tpl1[1][o])
        atm2 = np.interp(w2, tpl2[0][o], tpl2[1][o])
    else:
        atm1 = np.interp(w2, *tpl1)
        atm2 = np.interp(w2, *tpl2)

    A = np.c_[atm1*0+1, np.log(atm1), np.log(atm2)]
    y = f[ok]
    if a1 is not None:
        A = A[:, [0,2]]
        c1 = airmass_to_c1(a1)
        y /= atm1 ** c1

    atm_par = np.linalg.lstsq(A, np.log(y), rcond=None)[0]
    
    if a1 is not None:
        atm_par = [atm_par[0], c1, atm_par[1]]

    return atm_par

def _calc_atm(uo, atm_par):
    # the normalisation scaling is ignored
    w2 = np.exp(uo)
    atm1 = np.interp(w2, *tplo1)   # only linear interpolation so far
    atm2 = np.interp(w2, *tplo2)
    yatmo = atm1**atm_par[1] * atm2**atm_par[2]
    return yatmo

def calc_atm(uo, atm_par, order=None):
    # o is the info, from which echelle order the model is to be used
    if is_echelle:
        global tplo1, tplo2
        yatmo = 1 + 0*uo
        for o in (range(len(tpl1[0])) if order is None else [order]):
            if o >= len(tpl1[0]):
                # model not available
                continue
            tplo1 = tpl1[0][o], tpl1[1][o]
            tplo2 = tpl2[0][o], tpl2[1][o]
            if order is None:
                yatmo[o] = _calc_atm(uo[o], atm_par)
            else:
                yatmo = _calc_atm(uo, atm_par)
        return yatmo
    else:
        return _calc_atm(uo, atm_par)
