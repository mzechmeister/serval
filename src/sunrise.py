def sun(YYYY, MM, DD, lon=-74.3, lat=40.9, zenith='nautical', rise=True):
   '''
   Sunrise/Sunset Algorithm

   Parameters
   ----------
   YYYY, MM, DD : date of sunrise/sunset
   lon, lat : latitude, longitude for location for sunrise/sunset
      NOTE: longitude is positive for East and negative for West
   zenith : Sun's zenith for sunrise/sunset
      'offical'      = 90 deg 50'
      'civil'        = 96 deg
      'nautical'     = 102 deg
      'astronomical' = 108 deg

   Source
   ------
   http://williams.best.vwh.net/sunrise_sunset_algorithm.htm
   Almanac for Computers, 1990
   published by Nautical Almanac Office
   United States Naval Observatory
   Washington, DC 20392

   Example
   -------
   Worked example (from book):
   June 25, 1990:  25, 6, 1990
   Wayne, NJ:      40.9, -74.3
   Office zenith:  90 50' cos(zenith) = -0.01454
   >>> sun(1990, 6, 25, zenith='offical')
   9.441398443612975
   >>> sun(2005, 1, 1, lon=-70.7345, lat=-29.2584)  # HARPS
   8.78873517683543
   >>> sun(2005, 1, 1, lon=-70.7345, lat=-29.2584, rise=False)  # HARPS
   0.7622523167627335
   '''
   from math import floor, sin, cos, tan, asin, acos, atan, radians, degrees

   zenith = radians({'offical': 90 + 50/60.,
                     'civil': 96,
                     'nautical': 102,
                     'astronomical': 108}[zenith])

   # 1. first calculate the day of the year
   N1 = floor(275 * MM / 9)
   N2 = floor((MM + 9) / 12)
   N3 = 1 + floor((YYYY - 4 * floor(YYYY/4) + 2) / 3)
   N = N1 - N2*N3 + DD - 30

   # 2. convert the longitude to hour value and calculate an approximate time
   lngHour = lon / 15.

   if rise: # rising time is desired:
      t = N + (6-lngHour) / 24
   else:    # setting time is desired:
      t = N + (18-lngHour) / 24

   # 3. calculate the Sun's mean anomaly
   M = 0.9856*t - 3.289   # [deg]

   # 4. calculate the Sun's true longitude
   L = M + 1.916 * sin(radians(M)) + 0.020 * sin(2*radians(M)) + 282.634
   L = L % 360   # NOTE: L potentially needs to be adjusted into the range [0,360)

   # 5a. calculate the Sun's right ascension
   RA = degrees(atan(0.91764 * tan(radians(L))))
   RA = RA % 360 # NOTE: RA potentially needs to be adjusted into the range [0,360)

   # 5b. right ascension value needs to be in the same quadrant as L
   Lquadrant  = floor( L/90) * 90
   RAquadrant = floor(RA/90) * 90
   RA = RA + (Lquadrant - RAquadrant)

   # 5c. right ascension value needs to be converted into hours
   RA = RA / 15

   # 6. calculate the Sun's declination
   sinDec = 0.39782 * sin(radians(L))
   cosDec = cos(asin(sinDec))

   # 7a. calculate the Sun's local hour angle
   lat = radians(lat)
   cosH = (cos(zenith) - sinDec*sin(lat)) / (cosDec*cos(lat))

   if cosH > 1:
      print('the sun never rises on this location (on the specified date)')
      return None
   elif cosH < -1:
      print('the sun never sets on this location (on the specified date)')
      return None

   # 7b. finish calculating H and convert into hours
   H = degrees(acos(cosH))
   if rise:
      H = 360 - H
   H = H / 15.

   # 8. calculate local mean time of rising/setting
   T = H + RA - 0.06571*t - 6.622

   # 9. adjust back to UTC
   UT = T - lngHour
   UT = UT % 24 # NOTE: UT potentially needs to be adjusted into the range [0,24)

   return UT
   '''

   Example:
   UT = 4.488 - -4.953
      = 9.441
      = 9h 26m

10. convert UT value to local time zone of latitude/longitude

   localT = UT + localOffset

   Example:
   localT = 9h 26m + -4
          = 5h 26m
          = 5:26 am EDT
   '''

if __name__ == "__main__":
   import doctest
   doctest.testmod()

