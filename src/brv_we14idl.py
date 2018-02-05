import os
import sys
import subprocess


# barycorr.sav was created by starting IDL plain, running barycorr and bjd2utc once and then storing all compiled routines

# IDL> SAVE, /ROUTINES, FILENAME='barycorr.sav'


# cp ~/idl/exofast/bary/barycorr.sav .
# cp ~/idl/exofast/bary/JPLEPH.405 .


'''
# on command line this works
idl -quiet -IDL_STARTUP ""
bjd2utc=0 & restore, 'barycorr.sav'&!EXCEPT=0  & setenv, "ASTRO_DATA=/home/raid0/zechmeister/"  & e = {pmra:2888.92, pmdec:410.10, px:278.76, obsname:"ca"} & bjd = bjd2utc(2457395.24563d, 4.58559072, 44.02195596) & brv = barycorr(bjd, 4.58559072, 44.02195596,0., _extra=e, epoch=2448349.0625) 
print, bjd, brv
'''

# the following will not work. "-e" executes only one line, but function bjd2utc is not recognised when restoring in the same line.
'''
idl = subprocess.check_output(['idl','-quiet', '-IDL_STARTUP', '""', '-e', 'restore, "barycorr.sav"& !EXCEPT=0  & setenv, "ASTRO_DATA=/home/raid0/zechmeister/"  & e = {pmra:2888.92, pmdec:410.10, px:278.76, obsname:"ca"} & bjd=bjd2utc(2457395.24563d, 4.58559072, 44.02195596) & brv = barycorr(bjd, 4.58559072, 44.02195596,0., _extra=e, epoch=2448349.0625) & print, bjd, brv'])
    , shell=True, stdin=subprocess.PIPE,
                universal_newlines=True, bufsize=0)  # This line is needed for python3! Unbuffered and to 
a=subprocess.check_output(['idl','-quiet', '-IDL_STARTUP', '""', '-e', 'print,5'])

subprocess.check_output(['idl','-quiet', '-IDL_STARTUP', '""', '-e', 'restore, "barycorr.sav"& !EXCEPT=0  & setenv, "ASTRO_DATA=/home/raid0/zechmeister/"  & e = {pmra:2888.92, pmdec:410.10, px:278.76, obsname:"ca"} \n bjd=bjd2utc(2457395.24563d, 4.58559072, 44.02195596) \n brv = barycorr(bjd, 4.58559072, 44.02195596,0., _extra=e, epoch=2448349.0625) & print, bjd, brv'])

# We need to send stdin so we use Popen and communicate
idl = subprocess.Popen(['idl', '-quiet', '-IDL_STARTUP', '""'], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
idl.stdin.write('restore, "barycorr.sav"& !EXCEPT=0  & setenv, "ASTRO_DATA=/home/raid0/zechmeister/"  & e = {pmra:2888.92, pmdec:410.10, px:278.76, obsname:"ca"} \n bjd=bjd2utc(2457395.24563d, 4.58559072, 44.02195596) \n brv = barycorr(bjd, 4.58559072, 44.02195596,0., _extra=e, epoch=2448349.0625) & print, bjd, brv\nexit\n')
stdout, stderr = idl.communicate()
bjd, brv = stdout.split()[-2:]
'''

def bjdbrv(jd_utc, ra, dec, obsname=None, lat=None, lon=None, elevation=None,
        pmra=0.0, pmdec=0.0, parallax=0.0, rv=0.0, zmeas=0.0,
        epoch=2451545.0, tbase=0.0, raunits='degrees', deunits='degrees'):
   """
   Wrapper to idl for barycorr.pro and utc2bjd. Computes the barycentric
   velocity correction and julian date in one call.
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

   Example:
   --------
   >>> from brv_we14idl import bjdbrv
   >>> print bjdbrv(2457395.24563, 4.585590721, 44.02195596, 'ca')
   (2457395.247062143, -23684.543818733157)

   """
   params = {
        'PMRA': pmra,
        'PMDEC': pmdec,
        'PARALLAX': parallax,
        'RV': rv,
        'ZMEAS': 0.0,  # Applied manually below
        'EPOCH': epoch,
        'TBASE': tbase,
        'RAUNITS': "'"+raunits+"'",
        'DEUNITS': "'"+deunits+"'"
   }

   # Set observatory
   if obsname:
      params['OBSNAME'] = "'"+obsname+"'"
   elif None not in (lat, lon, elevation):
      params['LAT'] = lat
      params['LON'] = lon
      params['ELEVATION'] = elevation
   else:
      raise BarycorrError('Observatory location not set')

   cmd = ['e = {' + (", ".join(["%s:%s"% p for p in params.items()]))+'}',
          'bjd = utc2bjd(%sd, %sd, %sd)' % (jd_utc, ra, dec),
          'brv = barycorr(%sd, %sd, %sd, 0., _extra=e, epoch=2448349.0625)' % (jd_utc, ra, dec),
          'print, bjd, brv, FORMAT="(F,F)"']

   stdout = _query_idl(cmd)
   bjd, brv = stdout.split()[-2:]
   return float(bjd), float(brv)

def bvc(jd_utc, ra, dec, obsname=None, lat=None, lon=None, elevation=None,
        pmra=0.0, pmdec=0.0, parallax=0.0, rv=0.0, zmeas=0.0,
        epoch=2451545.0, tbase=0.0, raunits='degrees', deunits='degrees'):
   """
   Wrapper to idl for barycorr.pro and compute the barycentric
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

   Example:
   --------
   >>> from brv_we14idl import bvc
   >>> print bvc(2457395.24563, 4.585590721,  44.02195596, 'ca')
   -23684.5438187

; INPUTS:
;  JDUTC    - The full julian date in UTC, (scalar or vector), e.g.,
;             2450000d0. Must be double precision.
;  RA/DEC   - The ICRS coordinates of the target (scalar degrees), in EPOCH
;             (default J2000). See Fig 10 for required precision.
;  ZMEAS    - The measured redshift (e.g., the result of cross correlation
;             with template spectrum). Scalar or vector.
;
; OPTIONAL INPUTS:
;
;  OBSNAME  - The name of the observatory, as an input to
;             OBSERVATORY.PRO to retreive the latitude, longitude, and
;             altitude of the observing station. Either this,
;             latitude, longitude, and altitude, or R_ITRF must be
;             specified. If set, lat, long, alt, and r_itrf are
;             ignored.
;             See Fig 14 for required precision.
;  LAT      - Latitude (geodetic) of the observatory, in degrees with
;             respect to the WGS84 datum (scalar)
;  LONG     - West longitude of the observatory, in degrees with
;             respect to the WGS84 datum (scalar)
;  ALT      - Altitude of the observatory, in meters (scalar)
;  R_ITRF   - Three element array containing the XYZ coordinates of
;             observatory relative to geocenter, in meters. Only used
;             if obsname, lat, long, and alt are not specified.
;  EPOCH    - The epoch of the coordinates in JD_TDB, default is
;             julday(1,1,2000,12) = 2451545d0 => J2000
;             Overridden by HIP keyword.
;  TBASE    - The baseline that has been subtracted from the input
;             JD_UTCs, for higher precision times.
;
;             **** THESE ARE REQUIRED FOR CM/S PRECISION ****
;  PMRA     - Proper motion in RA, in mas/year (scalar)
;             Must be in units of arc (as is typical for most
;             modern catalogs); i.e. PMRA = d(RA)/dt * cos(DEC)
;  PMDEC    - Proper motion in dec, in mas/year (scalar)
;             See Fig 13 for required precision.
;  PX       - Parallax of target, in mas (scalar) (See Fig 11)
;  RV       - Radial velocit of target, in m/s
;             ** only ~100 km/s precision required for cm/s precision
;             decade timescales for nearby disk star (see Fig 12)**
   """
   params = {
        'PMRA': pmra,
        'PMDEC': pmdec,
        'PARALLAX': parallax,
        'RV': rv,
        'ZMEAS': '0.0',  # Applied manually below
        'EPOCH': epoch,
        'TBASE': tbase,
        'RAUNITS': "'"+raunits+"'",
        'DEUNITS': "'"+deunits+"'"
   }

   # Set observatory
   if obsname:
      params['OBSNAME'] = "'"+obsname+"'"
   elif None not in (lat, lon, elevation):
      params['LAT'] = lat
      params['LON'] = lon
      params['ELEVATION'] = elevation
   else:
      raise BarycorrError('Observatory location not set')

   cmd = ['e = {' + (", ".join(["%s:%s"% p for p in params.items()]))+'}',
          'brv = barycorr(%sd, %sd, %sd, 0., _extra=e)' % (jd_utc, ra, dec),
          'print, brv, FORMAT="(F,F)"']

   stdout = _query_idl(cmd)
   brv = stdout.split()[-1]
   return float(brv)

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

    # Prepare parameters
    params = {
      #  'JDS': ','.join(map(repr, jd_utc)),
        'RA': ra,
        'DEC': dec,
        'FUNCTION': 'utc2bjd'
    }

    cmd = ['print, utc2bjd(%sd, %sd, %sd), FORMAT="(F)"'% (jd_utc, ra, dec)]

    bjd = float(_query_idl(cmd))

    return bjd


def bjd2utc(bjd_tdb, ra, dec):
    """
    Query idl for bjd2utc.pro and compute the Julian Date (UTC)
    for bjd_tdb.

    See also: http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html

    :param bjd_tdb: Barycentric Julian Date (TDB)
    :param ra: RA (J2000) [deg]
    :param dec: Dec (J2000) [deg]
    :return: JD(UTC) at ~20 ms accuracy (observer at geocenter)
    """

    ## Check if there are multiple values of jd_utc
    #if not isinstance(bjd_tdb, (list, ndarray)):
        #bjd_tdb = [bjd_tdb]

    cmd = ['print, bjd2utc(%sd, %sd, %sd), FORMAT="(F)"'% (bjd_tdb, ra, dec)]

    utc = float(_query_idl(cmd))

    return utc


def _query_idl(cmd):
   """
   Query idl with cmd and return results.
   """
   barycorrsav = os.path.join(os.path.dirname(__file__), "barycorr.sav")
   if not os.environ.get('ASTRO_DATA'):
      raise EnvironmentError('environment variable ASTRO_DATA is not set')

   idl = subprocess.Popen(['idl', '-quiet', '-IDL_STARTUP', '""'], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
   cmd = ['restore, "%s"'%barycorrsav
         , '!EXCEPT=0'   # suppress overflow messages
         #, 'setenv, "ASTRO_DATA=/home/raid0/zechmeister/"'
         ] + cmd +['exit\n']

   #print '\n'.join(cmd)
   idl.stdin.write('\n'.join(cmd))
   stdout, stderr = idl.communicate()
   return stdout

#print bjdbrv(2457395.24563, 4.585590721,  44.02195596, 'ca')
