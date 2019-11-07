import barycorrpy
from astropy.time import Time
from astropy import time, coordinates as coord, units as u

# https://github.com/shbhuk/barycorrpy/issues/27
# https://github.com/shbhuk/barycorrpy/wiki/11.-USNO-IERS-servers-are-down
# https://github.com/astropy/astropy/issues/9427
from astropy.utils.iers import conf as iers_conf
iers_conf.iers_auto_url = 'https://astroconda.org/aux/astropy_mirror/iers_a_1/finals2000A.all'
iers_conf.auto_max_age = None

def bjdbrv(jd_utc, ra, dec, obsname=None, lat=0., lon=0., elevation=None,
        pmra=0., pmdec=0., parallax=0., rv=0., zmeas=0.,
        epoch=2451545.0, tbase=0., **kwargs):
   """
   Wrapper to barycorrpy.py and utc2bjd. Computes the barycentric
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
   >>> from brv_we14py import bjdbrv
   >>> print bjdbrv(2457395.24563, 4.585590721,  44.02195596, 'ca')
   (2457395.247062386, -23684.54364462639)

   """
   # translation obsname_idl obsname_py
   if obsname=='ca':
      lat = 37.2236
      lon = -2.5463
      elevation = 2168.
   if obsname=='eso':
      #obsname = 'lasilla'
      lat = -29.2584
      lon = -70.7345
      elevation = 2400.
   if obsname=='lapalma':
       lat = 28.754000
       lon = -17.88905555
       elevation = 2387.2

   # Barycentric Julian Date
   # adapted from http://docs.astropy.org/en/stable/time/#barycentric-and-heliocentric-light-travel-time-corrections
   targ = coord.SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
   loc = coord.EarthLocation.from_geodetic(lon, lat, height=elevation)
   #times = time.Time(jd_utc, format='jd', scale='utc', location=loc)
   #ltt_bary = times.light_travel_time(targ)
   JDUTC = Time(jd_utc, format='jd', scale='utc')
   if JDUTC.isscalar:
      ltt_bary = JDUTC.light_travel_time(targ, location=loc)
      # Does not work vectorised with numpy 1.14
      # *** TypeError: For this input type lists must contain either int or Ellipsis
      # https://github.com/astropy/astropy/issues/7051
      bjd = JDUTC.tdb + ltt_bary
   else:
      bjd = [(jdutc.tdb + jdutc.light_travel_time(targ, location=loc)).value for jdutc in JDUTC]

   brv, warning_and_error, status = barycorrpy.get_BC_vel(JDUTC, ra=ra, dec=dec, epoch=epoch, pmra=pmra,
                   pmdec=pmdec, px=parallax, lat=lat, longi=lon, alt=elevation, **kwargs)

   return (bjd.value, brv[0]) if JDUTC.isscalar else (bjd, brv)


# print bjdbrv(2457395.24563, 4.585590721,  44.02195596, 'ca', leap_update=False)

