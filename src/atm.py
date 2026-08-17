import numpy as np
from astropy.io import fits

airmass_to_c1 = lambda x: x

def load(filename):
    ''' Load atmospheric model from file

    Standard files that ship with SERVAL and currently supported are
    'stdatm.fits'
        Description: Standard atmospheric model
    'stdAtmos_vis.fits'
        Description: Standard atmospheric model from MOLECFIT that still needs a convolution
    'atm_carm_nir.fits' or 'atm_carm_vis.fits'
        Description: Data driven atmospheric model for the CARMENES VIS or NIR, derived from a large number of telluric standard spectra. The model is provided as a function of airmass and echelle order.

    Also a user can provide their own atmospheric model in the format of a FITS file with the following structure:
    - The first extension (index 0) should contain the header with the following keywords:
        - C1_REF: Reference value for the linear relation to transform airmass to coefficients
        - C1_SCALE: Scale value for the linear relation to transform airmass to coefficients
    - The second extension (index 1) should contain the data with the following columns:
        - WAVE: Wavelength values
        - O2: Oxygen transmission values
        - H2O: Water vapor transmission values

    Parameters
    ----------
    filename : str
        path to the atmospheric model file

    '''

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


def fit_atm_par(ln_wave, f, a1=None, o=None):
    ''' Determine the atmospheric parameters from a given spectrum and 2 templates. If the airmass is known, it can be provided to reduce the number of free parameters.

    Parameters
    ----------
    ln_wave: np.ndarray
        natural log of the wavelength, i.e. log(wavelength)
    f : np.ndarray
        flux of the observed spectrum
    a1 : float, optional
        airmass of the observation; if known it can be provided to reduce the number of free parameters
    o : int, list, slice or np.ndarray, optional
        echelle order

    Returns
    -------
    atm_par : np.ndarray
        fitted atmospheric parameters
    '''

    ok = np.isfinite(f)
    tpl_ok = np.isfinite(tpl1[1][o]) & np.isfinite(tpl2[1][o])
    w2 = np.exp(ln_wave[ok])

    if is_echelle:
        # quick workaround for the NIR, tpl_ok maps a 2D array to 1D, it breaks if orders overlap, so be careful in the future
        # paste together o30 and o33 in NIR
        atm1 = np.interp(w2, tpl1[0][o][tpl_ok], tpl1[1][o][tpl_ok])
        atm2 = np.interp(w2, tpl2[0][o][tpl_ok], tpl2[1][o][tpl_ok])

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
        # insert the airmass coefficient c1 into the parameter array
        atm_par = np.array([atm_par[0], c1, atm_par[1]])

    return atm_par

def _calc_atm(ln_wave, atm_par):
    ''' Wrapper function to calculate the atmospheric transmission from the templates and the atmospheric parameters

    Parameters
    ----------
    ln_wave : np.ndarray
        natural log of the wavelength, i.e. log(wavelength)
    atm_par : np.ndarray
        atmospheric parameters

    Returns
    -------
    yatmo : np.ndarray
        calculated atmospheric transmission
    '''
    # the normalisation scaling is ignored
    wave = np.exp(ln_wave)
    atm1 = np.interp(wave, *tplo1)   # only linear interpolation so far
    atm2 = np.interp(wave, *tplo2)
    yatmo = atm1**atm_par[1] * atm2**atm_par[2]

    yatmo = np.where(yatmo < 0.02, np.nan, yatmo)
    #yatmo = np.nan_to_num(yatmo, nan=1.0, posinf=1.0, neginf=1.0)
    return yatmo

def calc_atm(ln_wave, atm_par, order=None):
    ''' Calculate the atmospheric transmission from the templates and the atmospheric parameters.
        If the order is provided, only that order is calculated, otherwise all orders are calculated.

    Parameters
    ----------
    ln_wave : np.ndarray
        natural log of the wavelength, i.e. log(wavelength)
    atm_par : np.ndarray
        atmospheric parameters
    order : int, optional
        echelle order, if provided only that order is calculated, otherwise all orders are calculated
    Returns
    -------
    yatmo : np.ndarray
        calculated atmospheric transmission
    '''

    # o is the info, from which echelle order the model is to be used
    if is_echelle:
        global tplo1, tplo2

        # create an array to hold the atmospheric transmission for all orders
        yatmo = np.ones_like(ln_wave)

        for o in (range(len(tpl1[0])) if order is None else [order]):
            if o >= len(tpl1[0]):
                # model not available
                continue
            tpl_ok = np.isfinite(tpl1[0][o]) & np.isfinite(tpl2[0][o])

            tplo1 = tpl1[0][o][tpl_ok], tpl1[1][o][tpl_ok]
            tplo2 = tpl2[0][o][tpl_ok], tpl2[1][o][tpl_ok]

            if order is None:
                yatmo[o] = _calc_atm(ln_wave[o], atm_par)
            else:
                yatmo = _calc_atm(ln_wave, atm_par)

        return yatmo
    else:
        return _calc_atm(ln_wave, atm_par)
