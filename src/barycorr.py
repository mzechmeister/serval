# -*- coding: utf-8 -*-
# Python 2 compatibility
from __future__ import division

# Required packages
#import requests
import urllib2
from numpy import array, ndarray


"""
Python routines that query Jason Eastman's web applets for barycentric
velocity and time correction (barycorr.pro, utc2bjd.pro and bjd2utc.pro)

When using one of the services provided through this module, please cite the
corresponding paper:
    Wright and Eastman (2014), PASP 126, pp. 838–852 [BVC calculation]
    Eastman et al. (2010), PASP 122, pp. 935–946 [BJD calculations]

The Python interface is written by René Tronsgaard (Aarhus University) and may
be used, modified or redistributed without restrictions. However, please do not
remove this message.

More info: http://github.com/tronsgaard/barycorr

See also: 
    http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.html
    http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html
    http://astroutils.astronomy.ohio-state.edu/time/bjd2utc.html
    
MZ: requests replaced by urllib2
MZ: barycorr.py renamed to brv_we14html.py
"""

__version__ = '1.1'
__all__ = ['bvc', 'utc2bjd', 'bjd2utc']

# Speed of light
_c = 299792458.


def bvc(jd_utc, ra, dec, obsname=None, lat=None, lon=None, elevation=None,
        pmra=0.0, pmdec=0.0, parallax=0.0, rv=0.0, zmeas=0.0,
        epoch=2451545.0, tbase=0.0, raunits='degrees', deunits='degrees'):
    """
    Query the web interface for barycorr.pro and compute the barycentric
    velocity correction.
    Keyword obsname refers to observatory.pro in the IDL Astronomy User Library

    See also: http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.html

    :param jd_utc: Julian date (UTC)
    :param ra: RA (J2000) [deg]
    :param dec: Dec (J2000) [deg]
    :param obsname: Observatory name (overrides coordinates if set)
    :param lat: Observatory latitude  [deg]
    :param lon: Observatory longitude (E) [+/-360 deg]
    :param elevation: Observatory elevation [m]
    :param pmra: Proper motion (RA*cos(Dec)) [mas/yr]
    :param pmdec: Proper motion (Dec) [mas/yr]
    :param parallax: Parallax [mas]
    :param rv: Radial velocity (within 100 km/s) [m/s]
    :param zmeas: Measured redshift
    :param epoch: Epoch (default 2448348.56250, J2000)
    :param tbase: Baseline subtracted from times (default 0.0)
    :return: Barycentric correction for zmeas
    """

    # Check if there are multiple values of jd_utc
    if not isinstance(jd_utc, (list, ndarray)):
        jd_utc = [jd_utc]

    # Prepare GET parameters
    params = {
        'JDS': ','.join(map(repr, jd_utc)),
        'RA': ra,
        'DEC': dec,
        'PMRA': pmra,
        'PMDEC': pmdec,
        'PARALLAX': parallax,
        'RV': rv,
        'ZMEAS': '0.0',  # Applied manually below
        'EPOCH': epoch,
        'TBASE': tbase,
        'RAUNITS': raunits,
        'DEUNITS': deunits
    }

    # Set observatory
    if obsname is not None:
        params['OBSNAME'] = obsname
    elif None not in (lat, lon, elevation):
        params['LAT'] = lat
        params['LON'] = lon
        params['ELEVATION'] = elevation
    else:
        raise BarycorrError('Observatory location not set')

    # Query the web server
    result = _query_webserver(
        'http://astroutils.astronomy.ohio-state.edu/exofast/barycorr.php',
        params,
        len(jd_utc)
    )

    # Apply the correction (workaround to allow multiple values for z_meas)
    zmeas = array(zmeas)
    zb = result / _c
    return _c * ((1. + zmeas) * (1. + zb) - 1)  # Corrected velocity


def utc2bjd(jd_utc, ra, dec):
    """
    Query the web interface for utc2bjd.pro and compute the barycentric
    Julian Date for each value in jd_utc.

    See also: http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html

    :param jd_utc: Julian date (UTC)
    :param ra: RA (J2000) [deg]
    :param dec: Dec (J2000) [deg]
    :return: BJD(TDB) at ~20 ms accuracy (observer at geocenter)
    """

    # Check if there are multiple values of jd_utc
    if not isinstance(jd_utc, (list, ndarray)):
        jd_utc = [jd_utc]

    # Prepare GET parameters
    params = {
        'JDS': ','.join(map(repr, jd_utc)),
        'RA': ra,
        'DEC': dec,
        'FUNCTION': 'utc2bjd'
    }

    # Query the web server
    return _query_webserver(
        'http://astroutils.astronomy.ohio-state.edu/time/convert.php',
        params,
        len(jd_utc)
    )


def bjd2utc(bjd_tdb, ra, dec):
    """
    Query the web interface for bjd2utc.pro and compute the Julian Date (UTC)
    for each value in bjd_tdb.

    See also: http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html

    :param bjd_tdb: Barycentric Julian Date (TDB)
    :param ra: RA (J2000) [deg]
    :param dec: Dec (J2000) [deg]
    :return: JD(UTC) at ~20 ms accuracy (observer at geocenter)
    """

    # Check if there are multiple values of jd_utc
    if not isinstance(bjd_tdb, (list, ndarray)):
        bjd_tdb = [bjd_tdb]

    # Prepare GET parameters
    params = {
        'JDS': ','.join(map(repr, bjd_tdb)),
        'RA': ra,
        'DEC': dec,
        'FUNCTION': 'bjd2utc'
    }

    # Query the web server
    return _query_webserver(
        'http://astroutils.astronomy.ohio-state.edu/time/convert.php',
        params,
        len(bjd_tdb)
    )


def _query_webserver(server_url, params, expected_length):
    """
    Query server_url with params and return results if length matches
    expected_length (internal function).
    """

    # Fire the HTTP request
    try:
#        r = requests.get(server_url, params=params)
        g = "&".join([ "%s=%s"%x for x in params.items()])
        #print(server_url+'?'+g)
        r = urllib2.urlopen(server_url+'?'+g)
    except:
        raise BarycorrError(
            'Could not connect to web server ({})'.format(server_url)
        )

    # Convert multiline string output to numpy float array
    try:
#        result = [float(x) for x in r.text.splitlines() if len(x) > 0]
        result = [float(x) for x in r.read().splitlines() if len(x) > 0]
        if len(result) != expected_length:
            raise BarycorrError(
                'Unexpected length of result\n{}'.format(r.url)
            )
        if expected_length == 1:
            return result[0]
        else:
            return array(result)
    except:
        raise BarycorrError(
            'Could not parse output from web server:\n{}'.format(r.url)
        )


class BarycorrError(BaseException):
    pass
