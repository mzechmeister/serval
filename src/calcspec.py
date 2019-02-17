import numpy as np
from scipy import interpolate

from pause import pause
import cubicSpline
from read_spec import def_wlog


c = 299792.4580   # [km/s]

def redshift(x, vo=0., ve=0.):
   """
   x: The measured wavelength.
   v: Speed of the observer [km/s].
   ve: Speed of the emitter [km/s].

   Returns:
      The emitted wavelength l'.

   Notes:
      f_m = f_e (Wright & Eastman 2014)

   """
   if np.isnan(vo): vo = 0     # propagate nan as zero (@calibration in fib B)
   a = (1.0+vo/c) / (1.0+ve/c)
   if def_wlog:
      return x + np.log(a)   # logarithmic
      #return x + a          # logarithmic + approximation v << c
   else:
      return x * a
      #return x / (1.0-v/c)

def dopshift(x, v=0.):
   """
   Convenience function for redshift.

   x: The true measured wavelength.
   v: Speed of the emitter [km/s].

   Returns:
      The emitted wavelength l'.

   """
   return redshift(x, ve=v)

def barshift(x, v=0.):
   """
   Convenience function for redshift.

   x: The measured wavelength.
   v: Speed of the observer [km/s].

   Returns:
      The true wavelengths at the barycentre.

   """
   return redshift(x, vo=v)

def calcspec(x2, v, *a, **kwargs):
   """
   Doppler shifts the template and applies the polynom a using the Horner schema.
   p(x,a) * f(x, v)

   v : Radial velocity.
   a : Polynomial coefficients.
   fmod : If present, polynomial is applied on this model, otherwise it is loaded again.
   retpoly : If present the computed polynom is returned (instead of fmod * poly).
       This is needed e.g. to normalise the errors.
   """
   #fmod=cubicSpline.evalSpline(globvar.xt,globvar.yt,globvar.kt,x2*(1.0-v/c))
   #fmod=interpolate.splev(x2*(1.0-v/c),globvar.tck,der=0)
   #fmod=interpolate.splev(dopshift(x2,v),globvar.tck,der=0)
   #fmod = spline_ev(dopshift(x2,v), globvar.tck)
   # static variables calcspec.tck calcspec.wcen must be preset
   if not kwargs.get('retpoly'):
      fmod = kwargs.get('fmod', calcspec.tpl(dopshift(x2,v))) # parsed fmod or new
   x2 = x2 - calcspec.wcen
   # Use a Horner schema to evaluate the polynom
   # poly = a_0 + a_1 x + a_2 x^2 + ... = a_0 + x (a_1 + x (a_2+ ...))
   poly = a[-1]
   #for b in a[-2::-1]: poly = b+x2*poly  # add up the contributions
   for b in a[-2::-1]:  # add up the contributions p = ((a_3*x+a_2)*x+a_1)*x+a_0
      poly *= x2
      poly += b
   return poly if kwargs.get('retpoly') else fmod * poly

def qfacmask(lnw,f,sz=26,thres=3000.,plot=False):
   '''mask the parts of the spectrum which have a a Q-factor higher than thres
      calculate moving Q factor'''
   lnf = np.log(f)
   qi = (np.gradient(lnf)/np.gradient(lnw))**2 # ((lnf[1:]-lnf[:-1])/(lnw[1:]-lnw[:-1]))**2
   Q = np.zeros(len(qi)+sz) # includes padded zeros
   Q[sz/2:-sz/2] = np.cumsum(qi)
   Q[-sz/2:]=Q[-sz/2-1] # pad the last; the first have zero
   num = np.empty_like(qi)
   num[sz/2:-sz/2] = sz
   num[:sz/2] = np.arange(sz/2,sz)
   num[-sz/2:] = np.arange(sz/2,sz)[::-1]
   Qb = np.sqrt((Q[sz:]-Q[:-sz])/num)   # normalise
   idx = Qb>thres  # core index
   reg = True*idx
   for i,msk in enumerate(idx):
      if msk: reg[max(i-sz/2,1):i+sz/2] = True
   if plot:
      from gplot import *
      gplot(lnw,f,Qb," ,'' us 1:3 axis x1y2, "+str(thres)+" axis x1y2 t 'threshold'"); 
      ogplot(lnw[reg], f[reg],",", lnw[idx], f[idx])
      pause()
   return reg

