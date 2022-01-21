#! /usr/bin/env python
from __future__ import print_function

__author__ = 'Mathias Zechmeister'
__version__ = '2021-03-31'

description = '''
SERVAL - SpEctrum Radial Velocity AnaLyser (%s)
     by %s
''' % (__version__, __author__)

import argparse
import copy
import ctypes
from ctypes import c_void_p, c_double, c_int
# from datetime import datetime imported with read_spec!
import glob
import os
import resource
import stat as os_stat
import sys
import time

import importlib

import numpy as np
from numpy import std,arange,zeros,where, polynomial,setdiff1d,polyfit,array, newaxis,average
from scipy import interpolate, optimize
from scipy.optimize import curve_fit

from gplot import *
from pause import pause, stop
from wstat import wstd, wmean, wrms, rms, mlrms, iqr, wsem, nanwsem, nanwstd, naniqr, quantile
from read_spec import *   # flag, sflag, def_wlog
from calcspec import *
from targ import Targ
import cubicSpline
import cspline as spl
import masktools
import phoenix_as_RVmodel
from chi2map import Chi2Map

gplot2 = Gplot() # for a second plot window
gplot.bar(0).colors('classic')
gplot.tmp = None

if 'gplot_set' in locals():
   raise ImportError('Please update new gplot.py.')

if tuple(map(int,np.__version__.split('.'))) > (1,6,1):
   np.seterr(invalid='ignore', divide='ignore') # suppression warnings when comparing with nan.

def lam2wave(l, wlog=def_wlog):
   return np.log(l) if wlog else l

def wave2lam(w, wlog=def_wlog):
   return np.exp(w) if wlog else w

servalsrc = os.path.dirname(os.path.realpath(__file__)) + os.sep
servaldir = os.path.dirname(os.path.realpath(servalsrc)) + os.sep
servallib = servaldir + 'lib' + os.sep
os.environ['ASTRO_DATA'] = servalsrc

ptr = np.ctypeslib.ndpointer
_pKolynomial = np.ctypeslib.load_library('polyregression.so', servalsrc)
_pKolynomial.polyfit.restype = c_double
_pKolynomial.polyfit.argtypes = [
   ptr(dtype=np.float),  # x2
   ptr(dtype=np.float),  # y2
   ptr(dtype=np.float),  # e_y2
   ptr(dtype=np.float),  # fmod
   #ptr(dtype=np.bool),  # ind
   c_double,             # ind
   c_int, c_double, c_int,  # n, wcen, deg
   ptr(dtype=np.float),  # p
   ptr(dtype=np.float),  # lhs
   ptr(dtype=np.float)   # pstat
]
_pKolynomial.interpol1D.argtypes = [
   ptr(dtype=np.float),  # xn
   ptr(dtype=np.float),  # yn
   ptr(dtype=np.float),  # x
   ptr(dtype=np.float),  # y
   c_int, c_int          # nn, n
]

c = 299792.4580   # [km/s] speed of light

def nans(*args, **kwargs):
   return np.nan * np.empty(*args, **kwargs)

# default values
postiter = 3     # number of iterations for post rvs (postclip=3)
sp, fmod = None, None    # @getHalpha
apar = zeros(3)      # parabola parameters
astat = zeros(3*2-1)
alhs = zeros((3, 3))


def Using(point, verb=False):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    if verb: print('%s: usertime=%s systime=%s mem=%s mb' % (point,usage[0],usage[1],
                (usage[2]*resource.getpagesize())/1000000.0))
    return (usage[2]*resource.getpagesize())/1000000.0

class Logger(object):

   def __init__(self):
      self.flush = sys.stdout.flush
      self.terminal = sys.stdout
      self.logfile = None # open(logfilename, "a")
      self.logbuf = ''

#   def flush(self):
#      # dummy for function astropy iers download; progress bar will not be shown
#      pass

   def isatty(self): return False     # python + astropy+_download_file_from_source

   def write(self, message):   # fork the output to stdout and file
      self.terminal.write(message)
      if self.logfile:
         self.logfile.write(message)
      else:
         self.logbuf += message

   def logname(self, logfilename):
       self.logfile = open(logfilename, 'a')
       self.logfile.write(self.logbuf)
       print('logging to', logfilename)

def minsec(t): return '%um%.3fs' % divmod(t, 60)   # format time


class Vectorise(list):
    '''
    Bundles instances and gets and sets attributes on all instances.
    '''
    __getattr__ = lambda self, attr: np.array([getattr(s, attr) for s in self])

    def __setattr__(self, attr, vals):
        if hasattr(vals, '__len__'):
            [setattr(s, attr, val) for (s,val) in zip(self,vals)]
        else:
            [setattr(s, attr, vals) for s in self]


class interp:
   """interpolation similar to interpolate.interp1d but faster
   array must be sorted; 1D arrays only !!!
   """
   def __init__(self, x, y) :
      self.x = 1 * x # we like real arrays
      self.y = 1 * y
   def __call__(self, xx):
      yy = 0 * xx
      _pKolynomial.interpol1D(self.x, self.y, xx, yy, self.x.size, xx.size)
      return yy


class Tpl:
   def __init__(self, wk, fk, initfunc, evalfunc, mask=False, berv=None, vsini=None):
      '''
      wk : barycentric corrected wavelength
      '''
      ii = slice(None)
      if mask:
        ii = np.isfinite(fk)
      self.wk = self.wk0 = wk[ii]
      self.fk = self.fk0 = fk[ii]   # broadened and unbroadened flux
      self.vsini = vsini
      if vsini and self.wk[1]-self.wk[0]:
         # does not handle gaps! Only for serval templates. Phoenix is not log-uniform sampled.
         self.wk, self.fk = rotbroad(self.wk0, self.fk0, vsini)
      self.berv = berv
      self.initfunc = initfunc
      self.funcarg = self.initfunc(self.wk, self.fk)
      self.evalfunc = evalfunc
   def __call__(self, w, der=0):
      return self.evalfunc(w, self.funcarg, der=der)
   def mskatm(self, w, msk):
      # mask regions (atm, stellar) in template
      # need target velocity and velocity range
      if msk and self.berv:
 #         bb[o][tellmask(barshift(ww[o],-spt.berv))>0.01] |= flag.atm   # mask as before berv correction
#         bb[o][skymsk(barshift(ww[o],-spt.berv))>0.01] |= flag.sky   # mask as before berv correction
     # mask as before berv correction
         return msk(barshift(w,-self.berv)) > 0.01
      else:
         return slice(0, 0) # empty
   def rotbroad(self, vsini=0):
      if vsini > 0:
         # broaden template
         self.wk, self.fk = rotbroad(self.wk0, self.fk0, vsini)

         # update funcarg for evaluation
         self.funcarg = self.initfunc(self.wk, self.fk)


def rotbroad(x, f, v):
    '''
    Broaden a spectrum by rotation.

    Parameters
    ----------
    x : Uniform ln(lambda) grid.
    f : Spectrum (sampled uniformly in log(lambda), i.e. velocity).
    v : Rotational velocity vsini [km/s]

    Returns
    -------
    x : ln(lambda) truncated to valid range.
    f : Broaden spectrum.

    '''
    # f(x)= x<-1 ? -pi/4 : x>1? pi/4 : (sqrt(1-x**2)*x+asin(x))/2; pl sqrt(1-x**2), dx=0.8, (f(x+dx/2)-f(x-dx/2))/dx, f(x)

    dv = (x[1]-x[0]) * c   # velocity step [km/s]
    k = int(v / dv)        # kernel sampling

    # rotation broadening kernel
    A = np.sqrt(1 - (np.arange(-k,k+1.)/k)**2)
    A /= sum(A)    # normalise kernel to unity area

    return x[k:-k], np.convolve(f, A, mode='valid')

    frot = 0
    for i in range(-k, k+1):
        frot += A[v+i] * f[v-i:-v-i-1]  # -0 does not work
    return frot


def analyse_rv(obj, postiter=1, fibsuf='', oidx=None, safemode=False, pdf=False):
   """
   """
   print(obj+'/'+obj+'.rvc'+fibsuf+'.dat')
   allrv = np.genfromtxt(obj+'/'+obj+'.rvo'+fibsuf+'.dat')
   allerr = np.genfromtxt(obj+'/'+obj+'.e_rvo'+fibsuf+'.dat')
   sbjd = np.genfromtxt(obj+'/'+obj+'.rvo'+fibsuf+'.dat', dtype=('|S33'), usecols=[0]) # as string
   snr = np.genfromtxt(obj+'/'+obj+'.snr'+fibsuf+'.dat')

   if np.size(allrv) == 1:
      return   # just one line, e.g. drift

   bjd, RVc_old, e_RVc_old, RVd, e_RVd, RV_old, e_RV_old, BRV, RVsa = np.genfromtxt(obj+'/'+obj+'.rvc'+fibsuf+'.dat', unpack=True)

   orders, = np.where(np.sum(allerr[:,5:]>0, 0))   # orders with all zero error values
   if oidx is not None:
      omiss = set(oidx) - set(orders)
      if omiss: pause('WARNING: orders', omiss,'missing')
      else: orders = np.array(sorted(set(orders) & set(oidx)))

   omap = orders[:,newaxis]
   rv, e_rv = allrv[:,5+orders], allerr[:,5+orders]

   RV, e_RV = nanwsem(rv, e=e_rv, axis=1)
   RVc = RV - np.nan_to_num(RVd) - np.nan_to_num(RVsa)
   e_RVc = np.sqrt(e_RV**2 + np.nan_to_num(e_RVd)**2)

   # post RVs (re-weightening)
   ok = e_rv > 0
   ordmean = average(rv*0, axis=0)
   gplot.key('tit "'+obj+'"')
   for i in range(1+postiter): # centering, first loop to init, other to clip
      RVp, e_RVp = nanwsem(rv-ordmean, e=e_rv*ok, axis=1)
      orddisp = rv - ordmean - RVp[:,newaxis]
      ordstd, d_ordmean = wstd(orddisp, e_rv*ok, axis=0)
      ordmean += d_ordmean                # update ordmean
      ok &= np.abs(orddisp-d_ordmean) <= 3*ordstd  # clip and update mask
      if np.isnan(RVp.sum()) or np.isnan(e_RVp.sum()):
         print('WARNING: nan in post RVs. Maybe to few measurements. Re-setting to originals.\n')
         RVp, e_RVp = RV, e_RV
         break
      else:
         gplot(orddisp.T, ' matrix us (%s-1):3 t ""' % "".join(['$1==%s?%s:' % io for io in enumerate(orders)]))
         ogplot(orders, ordstd, ' w lp lt 3 t "", "" us 1:(-$2) w lp t ""')
      if 0: pause(i)
   rvp = rv - ordmean
   pdf=0
   if pdf:
      gplot.term('pdfcairo; set out "%s.pdf"'% obj)
   if 1: # chromatic slope
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"chromatic index [m/s/Np]"')
      gplot('"'+obj+'/'+obj+'.srv'+fibsuf+'.dat" us ($1-2450000):4:5 w e pt 7')
      if not safemode: pause('chromatic index')
      gplot.xlabel('"RV [m/s]"').ylabel('"chromatic index [m/s/Np]"')
      gplot('"" us 2:4:3:5:($1-2450000)  w xyerr pt 7 palette')
      if not safemode: pause('correlation RV - chromatic index')

   snro = snr[:,2+orders]
   print("total SNR:", np.sum(snro**2)**0.5)
   if 1: # plot SNR
      gplot.reset().xlabel('"Order"; set ylabel "SNR"; set ytics nomirr; set y2label "total SNR"; set y2tics; set yrange[0:]; set y2range[0:]')
      gplot(snro, 'matrix us ($2+%i):3' % np.min(orders), flush='')
      ogplot(np.sum(snro**2,axis=1)**0.5,' us (%i):1 axis x1y2 t "total SNR"'% (np.min(orders)+len(orders)))
      if not safemode: pause()
      gplot.reset()

   # store post processed rvs
   RVpc = RVp - np.nan_to_num(RVd) - np.nan_to_num(RVsa)
   e_RVpc = np.sqrt(e_RVp**2 + np.nan_to_num(e_RVd)**2)
   #unit_rvp = [open(obj+'/'+obj+'.post'+fibsuf+'.dat', 'w'), open(obj+'/'+obj+'.post.badrv'+fibsuf+'.dat', 'w')]
   np.savetxt(obj+'/'+obj+'.post'+fibsuf+'.dat', list(zip(sbjd, RVpc, e_RVpc, RVp, e_RVp, RVd, e_RVd, BRV, RVsa)), fmt='%s')

   print('Statistic on dispersion in RV time series for', obj)
   print('        %10s %10s %10s %10s %10s'% ('mlrms [m/s]', 'std [m/s]', 'wstd [m/s]', 'iqr [m/s]', 'wiqr [m/s]'))
   print('RVmed: %10.4f' % std(allrv[:,3]))
   print('RV:   '+' %10.4f'*5 % (mlrms(RV,e_RV)[0], std(RV), wstd(RV,e_RV)[0], iqr(RV,sigma=True), iqr(RV, w=1/e_RV**2, sigma=True)))
   print('RVp:  '+' %10.4f'*5 % (mlrms(RVp,e_RVp)[0], std(RVp), wstd(RVp,e_RVp)[0], iqr(RVp,sigma=True), iqr(RVp, w=1/e_RVp**2, sigma=True)))
   print('RVc:  '+' %10.4f'*5 % (mlrms(RVc,e_RVc)[0], std(RVc), wstd(RVc,e_RVc)[0], iqr(RVc,sigma=True), iqr(RVc, w=1/e_RVc**2, sigma=True)))
   print('RVpc: '+' %10.4f'*5 % (mlrms(RVpc,e_RVpc)[0], nanwstd(RVpc), nanwstd(RVpc,e=e_RVpc), naniqr(RVpc,sigma=True), iqr(RVpc, w=1/e_RVpc**2, sigma=True)))
   print('median internal precision', np.median(e_RV))
   print('Time span [d]: ', bjd.max()-bjd.min())

   gplot.reset().xlabel('"BJD - 2 450 000"; set ylabel "RV [m/s]"')
   gplot('"'+obj+'/'+obj+'.rvc'+fibsuf+'.dat" us ($1-2450000):2:3 w e pt 7 t "rvc %s"'%obj)
   if not safemode: pause('rvc')

   gplot(bjd, allrv[:,3],' us ($1-2450000):2 t "RVmedian" lt 4', flush='')
   ogplot(bjd, RV, e_RV, ' us ($1-2450000):2:3 w e t "RV"', flush='')
   ogplot(bjd, RVp, e_RVp, ' us ($1-2450000):2:3 w e t "RVp"', flush='')
   ogplot(bjd, RVc, e_RVc, ' us ($1-2450000):2:3 w e pt 7 lt 1 t "RVc"', flush='')
   ogplot(bjd, RVpc, e_RVpc,' us ($1-2450000):2:3 w e pt 6 lt 7 t "RVpc"')
   if not safemode: pause('rv')

   if 0:  # error comparison
      gplot_set('set xlabel "Order"')
      gplot(orders, ordstd, ' t "external error"')
      ogplot(orders, np.mean(e_rv, axis=0), ' t "internal error"')
      if not safemode: pause('ord disp')
   if 1:  # Order dispersion
      gplot.reset().xlabel('"Order"')
      #gplot('"'+filename,'" matrix every ::%i::%i us ($1-5):3' % (omin+5,omax+5))
      #ogplot(allrv[:,5:71],' matrix every ::%i us 1:3' %omin)
      # create ord, rv,e_rv, bb
      #ore = [(o+omin,orddisp[n,o],e_rv[n,o],~ok[n,o]) for n,x in enumerate(orddisp[:,0]) for o,x in enumerate(orddisp[0,:])]
      #ore = [(o,)+x for row in zip(orddisp,e_rv,~ok) for o,x in enumerate(zip(*row),omin)]
      ore = [np.tile(orders,orddisp.shape).ravel(), orddisp.ravel(), e_rv.ravel(), ~ok.ravel()]
      gplot(*ore + ['us 1:2:($3/30) w xe'], flush='')
      if not ok.all(): ogplot('"" us 1:2:($3/30/$4) w xe', flush='') # mark 3-sigma outliners
      #gplot(orddisp,' matrix  us ($1+%i):3' % omin, flush='')
      #gplot('"'+filename,'" matrix every ::'+str(omin+5)+' us ($1-5):3')
      #ogplot(ordmean,' us ($0+%i):1 w lp lt 3 pt 3  t "ord mean"' %omin, flush='')
      ogplot(orders, ordstd,' w lp lt 3 pt 3  t "1 sigma"', flush='')
      ogplot('"" us 1:(-$2) w lp lt 3 pt 3  t ""', flush='')
      ogplot('"" us 1:($2*3) w lp lt 4 pt 3  t "3 sigma"', flush='')
      ogplot('"" us 1:(-$2*3) w lp lt 4 pt 3  t ""', flush='')
      ogplot('"" us ($1+0.25):(0):(sprintf("%.2f",$2)) w labels rotate t""', flush='')
      ogplot(*ore+[ 'us 1:2:($3)  w point pal pt 6'])
      if not safemode: pause('ord scatter')
   print('mean ord std:', np.mean(ordstd), ', median ord std:',np.median(ordstd))

   if pdf:
      gplot.out()

   return allrv

def nomask(x):   # dummy function, tellurics are not masked
   return 0. * x


def lineindex(l, r1, r2):
   if np.isnan(r1[0]) or np.isnan(r2[0]):
      return np.nan, np.nan
   s = l[0] / (r1[0]+r2[0]) * 2
   # error propagation
   e = s * np.sqrt((l[1]/l[0])**2 + (r1[1]**2+r2[1]**2)/(r1[0]+r2[0])**2)
   return s, e



lines = {
         'Halpha': (6562.808, -15.5, 15.5),   # Kuerster et al. (2003, A&A, 403, 1077)
         'Halpha': (6562.808, -40., 40.),
         'Haleft': (6562.808, -300., -100.),
         'Harigh': (6562.808, 100, 300),
         'Haleft': (6562.808, -500., -300.),
         'Harigh': (6562.808, 300, 500),
         'CaI':    (6572.795, -15.5, 15.5),   # Kuerster et al. (2003, A&A, 403, 1077)
         'CaH':    (3968.470, -1.09/2./3968.470*c, 1.09/2./3968.470*c), # Lovis et al. (2011, arXiv1107.5325)
         'CaK':    (3933.664, -1.09/2./3933.664*c, 1.09/2./3933.664*c), # Lovis et al.
         'CaIRT1': (8498.02, -15., 15.),      # NIST + my definition
         'CaIRT1a': (8492, -40, 40),          # My definition
         'CaIRT1b': (8504, -40, 40),          # My definition
         'CaIRT2': (8542.09, -15., 15.),      # NIST + my definition
         'CaIRT2a': (8542.09, -300, -200),    # My definition
         'CaIRT2b': (8542.09, 250, 350),      # My definition, +50 due telluric
         'CaIRT3': (8662.14, -15., 15.),      # NIST + my definition
         'CaIRT3a': (8662.14, -300, -200),    # NIST + my definition
         'CaIRT3b': (8662.14, 200, 300),      # NIST + my definition
         'NaD1':   (5889.950943, -15., 15.),  # NIST + my definition
         'NaD2':   (5895.924237, -15., 15.),  # NIST + my definition
         'NaDref1': (5885, -40, 40),          # My definition
         'NaDref2': ((5889.950+5895.924)/2, -40, 40),   # My definition
         'NaDref3': (5905, -40, 40)           # My definition
}

def get_o_of_line(typ, wavemap):
   # find the orders of a spectral line
   wcen, dv1, dv2 = lines.get(typ, (None, None, None))
   wcen = lam2wave(airtovac(wcen))
   if not wcen:
      return None
   o = np.where((np.nanmin(wavemap,axis=1) < wcen) & (wcen < np.nanmax(wavemap, axis=1)))[0]
   if o.size > 0:
      o = o[0]   # Take the first. If another order is preferred, it should be overwritten with the inst file.
      o = np.argmin(np.nanmax(np.abs(wavemap-wcen), axis=1))   # take the one with which best centered
      return o
   else:
      return None

def getHalpha(v, typ='Halpha', line_o=None, rel=False, plot=False):
   """
   v [km/s]
   sp,fmod as global variables !
   deblazed sp should be used !
   """
   o = line_o.get(typ)
   if o is None: return np.nan, np.nan
   wcen, dv1, dv2 = lines[typ]
   wcen = lam2wave(airtovac(wcen))

   ind = o, (wcen+(v-sp.berv+dv1)/c < sp.w[o]) & (sp.w[o] < wcen+(v-sp.berv+dv2)/c)

   if rel: # index relative to template
      # relative index is not a good idea for variable lines
      if np.isnan(fmod[ind]).all(): return np.nan, np.nan
      if np.isnan(fmod[ind]).any(): pause()
      mod = fmod[ind]
      res = sp.f[ind] - mod
      I = np.mean(res/fmod[ind])
      e_I = rms(np.sqrt(sp.f[ind])/fmod[ind]) / np.sqrt(np.sum(ind[1])) # only for photon noise
   else:   # mean absolute flux ("box fitting" .... & ref reg =>  absolute index EW width)
      I = np.mean(sp.f[ind])
      e_I = 1. / ind[1].sum() * np.sqrt(np.sum(sp.e[ind]**2))
      mod = I + 0*sp.f[ind]

   if plot==True or typ in plot:
      print(typ, I, e_I)
      gplot(sp.w[o], fmod[o], sp.f[o], 'us 1:2 w l lt 3 t "template", "" us 1:3 lt 1 t "obs"')
      ogplot(sp.w[ind], sp.f[ind], mod, 'lt 1 pt 7 t "'+typ+'", "" us 1:3 w l t "model", "" us 1:($2-$3) t "residuals"')
      pause()

   return I, e_I


def polyreg(x2, y2, e_y2, v, deg=1, retmod=True):   # polynomial regression
   """Returns polynomial coefficients and goodness of fit."""
   fmod = calcspec(x2, v, 1.)  # get the shifted template
   if 0: # python version
      ind = fmod>0.01     # avoid zero flux, negative flux and overflow
      p,stat = polynomial.polyfit(x2[ind]-calcspec.wcen, y2[ind]/fmod[ind], deg-1, w=fmod[ind]/e_y2[ind], full=True)
      SSR = stat[0][0]
   else: # pure c version
      pstat = np.empty(deg*2-1)
      p = np.empty(deg)
      lhs = np.empty((deg, deg))
      # _pKolynomial.polyfit(x2, y2, e_y2, fmod, ind, ind.size,globvar.wcen, rhs.size, rhs, lhs, pstat, chi)
      ind = 0.0001
      # no check for zero division inside _pKolynomial.polyfit!
      #pause()
      SSR = _pKolynomial.polyfit(x2, y2, e_y2, fmod, ind, x2.size, calcspec.wcen, deg, p, lhs, pstat)
      if SSR < 0:
         ii, = np.where((e_y2<=0) & (fmod>0.01))
         print('WARNING: Matrix is not positive definite.', 'Zero or negative yerr values for ', ii.size, 'at', ii)
         p = 0*p
         #pause(0)

      if 0: #abs(v)>200./1000.:
         gplot(x2,y2,calcspec(x2, v, *p), ' w lp,"" us 1:3 w l, "" us 1:($2-$3) t "res [v=%f,SSR=%f]"'%(v, SSR))
         pause('v, SSR: ',v,SSR)
   if retmod: # return the model
      return p, SSR, calcspec(x2, v, *p, fmod=fmod)
   return p, SSR

def gauss(x, a0, a1, a2, a3):
   z = (x-a0) / a2
   y = a1 * np.exp(-z**2 / 2) + a3 #+ a4 * x + a5 * x**2
   return y

def optidrift(ft, df, f2, e2=None):
   """
   Scale derivative to the residuals.

   Model:
      f(v) = A*f - A*df/dv * v/c

   """
   # pre-normalise
   #A = np.dot(ft, f2) / np.dot(ft, ft)   # unweighted (more robust against bad error estimate)
   A = np.dot(1/e2**2*ft, f2) / np.dot(1/e2**2*ft, ft)
   fmod = A * ft
   v = -c * np.dot(1/e2**2*df, f2-fmod) / np.dot(1/e2**2*df, df) / A #**2
   if 0:
      # show RV contribution of each pixel!
      print('median', np.median(-(f2-fmod)/df*c/A*1000), ' v', v*1000)
      gplot(-(f2-fmod)/df*c/A*1000, e2/df*c/A*1000, 'us 0:1:2 w e, %s' % (v*1000))
      pause()
   e_v = c / np.sqrt(np.dot(1/e2**2*df, df)) / A
   fmod = fmod - A*df*v/c
   #gplot(f2, ',', ft*A, ',', f2-fmod,e2, 'us 0:1:2 w e')
   #gplot(f2, ',', ft*A, ',', f2-A*ft,e2, 'us 0:1:2 w e')
   #gplot(f2, ft*A, A*ft-A*df/c*v, A*ft-A*df/c*0.1,' us 0:1, "" us 0:2, "" us 0:3, "" us 0:4, "" us 0:($1-$3) w lp,  "" us 0:($1-$4) w lp lt 7')
   return  type('par',(),{'params': np.append(v,A), 'perror': np.array([e_v,1.0])}), fmod

def CCF(wt, ft, x2, y2, va, vb, e_y2=None, keep=None, plot=False, ccfmode='trapeze'):
   # CCF is on the same level as the least square routine fit_spec
   if keep is None: keep = np.arange(len(x2))
   if e_y2 is None: e_y2 = np.sqrt(y2)
   vgrid = np.arange(va, vb, v_step)

   # CCF is a data compression/smoothing/binning
   SSR = np.arange(vgrid.size, dtype=np.float64)

   def model_box(x, v):
       idx = np.searchsorted(wt+v/c, x)
       return ft[idx] > 0

   for i in SSR:
      if ccfmode == 'box':
         # real boxes (zero order interpolation)
         #idx = np.searchsorted(wt+vgrid[i]/c, x2[keep])
         #SSR[i] = np.dot(y2[keep], ft[idx]>0) / np.sum(ft[idx]>0)
         y2mod = model_box(x2[keep], vgrid[i])
         SSR[i] = np.dot(y2[keep], y2mod) / np.sum(y2mod)
         # gplot(np.exp(x2), y2, ft[idx]*y2.max(),'w lp, "" us 1:3 w l lt 3')
      if ccfmode == 'trapeze':
         stop()
         # gplot(np.exp(x2), y2,'w lp,', np.exp(x2[keep]),ft[idx]*5010, 'w l lt 3')

   # Peak analysis with gaussian fit
   # The Gaussian is a template for the second level product CCF
   try:
      params, covariance = curve_fit(gauss, vgrid, SSR, p0=(0., SSR.max()-SSR.min(), 2.5, SSR.min()))
      perror = np.diag(covariance) if np.isfinite(covariance).min() else 0*params

      SSRmod = gauss(vgrid, *params)
      SSRmin = gauss(params[0], *params)
   except:
      params = [0,0,0,0]
      perror = [0,0,0,0]
      SSRmod = 0
      SSRmin = 0

   if 0:
      # old CCF_binless version
      # uses just the mask which has a fixed window size
      idx = np.searchsorted(wt, x2[keep]) #v=0
      ind = (ft[idx]>0) & (ft[idx-1]>0)
      iidx = idx[ind] - 1 # previous/lower index
      xv2 = (x2[keep[ind]]-0.5*(wt[iidx+1,0]+wt[iidx])) * c
      yv2 = y2[keep[ind]] / ft[iidx]
      ev2 = e_y2[keep[ind]] / ft[iidx]
   elif ccfmode == 'binless':
      # get the line center from the mask
      # grab the data around the line center with the window size
      # phase fold all windows to velocity space by subtracting the line center
      # assume windows do not overlap
      wsz = 6.0 # [km/s] 3.0 for HARPS, 6.0 for FEROS
      linecen = 0.5 * (wt[1::4]+wt[2::4]) # get the linecenter form the mask
      lineflux = 0.5 * (ft[1::4]+ft[2::4])
      idx = np.searchsorted(linecen, x2[keep]) # find index of nearest largest line
      idx = idx - ((linecen[idx]-x2[keep]) > (x2[keep]-linecen[idx-1])) # index of the nearest (smaller or larger) line
      ind = np.abs(x2[keep]-linecen[idx]) < wsz/c # index for values with the window
      xv2 = (x2[keep[ind]]-linecen[idx[ind]]) * c # phase fold
      yv2 = y2[keep[ind]] / lineflux[idx[ind]]
      ev2 = e_y2[keep[ind]] / lineflux[idx[ind]]

      # normalise window flux with the flux mean in each window
      iidx = idx[ind] - 1
      lidx = iidx - iidx.min()
      nlines = lidx.max() + 1
      linenorm = np.zeros((nlines,2))
      for i,y in zip(lidx, yv2): linenorm[i] += np.array([y,1])
      linenorm = (linenorm[:,0] / linenorm[:,1])[lidx]  # mean = sum / n
      yv2norm = yv2 / linenorm
      #yv2 = yv2norm

   fmod = 0
   if 0 or plot and ccfmode=='binless':
      # try a robust gaussfit via iterative clipping
      iidx = idx[ind] - 1
      fmod = 0 * y2
      try:
         #par2, cov2 = curve_fit(gauss, xv2, yv2, p0=(0,1,2.5,0))
         nclip = 1
         good = np.arange(len(xv2))
         for it in range(nclip+1):
            par2, cov2 = curve_fit(gauss, xv2[good], yv2norm[good], p0=(0,1,2.5,0))
            perr2 = np.diag(cov2)
            SSR2mod = gauss(vgrid, *par2)
            yv2mod = gauss(xv2, *par2)
            normfac = params[1] / par2[1]
            fmod[keep[ind]] = gauss(xv2, *par2)
            # clip
            kap = 3 * np.std(yv2norm[good]-yv2mod[good])
            out = np.abs(yv2norm - yv2mod) > kap
            rml = np.unique(lidx[out])
            good = np.sum(lidx == rml[:,newaxis],0)==0 # ind for lines not reject (good windows)
            if plot:
               gplot(xv2, yv2norm, yv2mod,good, ', "" us 1:3, "" us 1:($2/($4==0)),', vgrid,SSR2mod,SSR2mod-kap,SSR2mod+kap, 'w l lt 7, "" us 1:3 w l lt 1, "" us 1:4 w l lt 1')
            #pause()
      except:
         par2 = [0]
         perr2 = [0,0,0]
         normfac = 0
         SSR2mod = vgrid*0

   if plot&1:
      gplot(vgrid, SSR, SSRmod, 'w lp pt 7 t "CCF", "" us 1:3 w l lt 3 lw 3 t"CCF fit %.2f +/- %.2f m/s"'%(params[0]*1000,perror[0]*1000))
   if plot&1 and ccfmode=='binless':
      # insert a -1 row for broken lines in plot
      gg = []
      for g,diff in zip(np.dstack((xv2, normfac*yv2norm, normfac*ev2, iidx-iidx.min()))[0], np.ediff1d(iidx, to_begin=0)): gg += [g*0-1,g] if diff else [g]
      gg = np.array(gg)

      gplot.palette('model HSV rgbformulae 3,2,2; set bar 0')
      gplot(gg, 'us 1:($2/($4>-1)):4 w lp palette pt 7, "" us 1:($2/($4>-1)):3 w e pt 1 lt 1,', vgrid, normfac*SSR2mod, 'w l lc 3 lw 3 t"%.2f +/- %.2f m/s"'%(par2[0]*1000,perr2[2]*1000))
      #gplot(xv2, normfac*yv2norm, normfac*ev2, iidx-iidx.min(), np.ediff1d(iidx, to_begin=0), 'us 1:2:4 w lp palette pt 7, "" us 1:2:3 w e pt 1 lt 1,', vgrid, normfac*SSR2mod, 'w l lc 0 lw 3 t"%.2f +/- %.2f m/s"'%(par2[0]*1000,perr2[2]*1000))
      ogplot(vgrid, SSR, SSRmod, 'w lp lt 7 pt 7 t "CCF", "" us 1:3 w l lt 7 lw 3 t"CCF fit %.2f +/- %.2f m/s"'%(params[0]*1000,perror[0]*1000))

      pause(params[0]*1000, perror[0]*1000)
      print(params[0]*1000, perror[0]*1000)

   if plot&1: pause()
   if plot&2:
      idx = np.searchsorted(wt, x2)
      idx0 = np.searchsorted(wt-vgrid[0]/c, x2)
      gplot(np.exp(x2), y2, ft[idx]/10, y2*(ft[idx]>0), y2*(ft[idx0]>0), 'w lp t "x2,y2 input spectrum", "" us 1:3 w l lt 3 t "mask (v=0)", "" us 1:4 lt 3 t "x2,y2*mask", "" us 1:5  lt 4 pt 1 t "x2,y2*mask(v=va)')
      pause()
   stat = {'std': 1, 'snr': 0}
   params[0] *= -1
   return type('par',(),{'params': params, 'perror':perror, 'ssr':SSRmin, 'niter':0}), fmod, (vgrid, SSR, SSRmod), stat

def SSRstat(vgrid, SSR, dk=1, plot=2):
   '''analyse peak
   plot: 0 - never
         1 - always
         2 - on exception
   '''
   k = SSR[dk:-dk].argmin() + dk   # best point (exclude borders)
   vpeak = vgrid[k-dk:k+dk+1]
   SSRpeak = SSR[k-dk:k+dk+1] - SSR[k]
   dv = vgrid[dk] - vgrid[0]
   # interpolating parabola a0+a1*x+a2*x**2 (direct solution) through the three pixels in the minimum
   a = np.array([0, (SSR[k+dk]-SSR[k-dk])/(2*dv), (SSR[k+dk]-2*SSR[k]+SSR[k-dk])/(2*dv**2)])  # interpolating parabola for even grid
   v = (SSR[k+dk]-SSR[k-dk]) / (SSR[k+dk]-2*SSR[k]+SSR[k-dk]) * 0.5 * dv

   v = vgrid[k] - a[1]/2./a[2]   # position of parabola minimum
   e_v = np.nan
   if -1 in SSR:
      print('opti warning: bad ccf.')
   elif a[2] <= 0:
      print('opti warning: a[2]=%f<=0.' % a[2])
   elif not vgrid[0] <= v <= vgrid[-1]:
      print('opti warning: v not in [va,vb].')
   else:
      e_v = 1. / a[2]**0.5
   if plot&1 or (plot&2 and np.isnan(e_v)):
      gplot.yrange('[*:%f]' % SSR.max())
      gplot(vgrid, SSR-SSR[k], " w lp, v1="+str(vgrid[k])+", %f+(x-v1)*%f+(x-v1)**2*%f," % tuple(a), [v,v], [0,SSR[1]], 'w l t "%f km/s"'%v)
      ogplot(vpeak, SSRpeak, ' lt 1 pt 6; set yrange [*:*]')
      pause(v)
   return v, e_v, a

def opti(va, vb, x2, y2, e_y2, p=None, vfix=False, plot=2):
   """ vfix to fix v for RV constant stars?
   performs a mini CCF; the grid stepping
   returns best v and errors from parabola curvature
   plot: 0 - never
         1 - always
         2 - on exception
   """
   vgrid = np.arange(va, vb, v_step)
   nk = len(vgrid)

   SSR = np.empty(nk)
   for k in range(nk):
      p, SSR[k] = polyreg(x2, y2, e_y2, vgrid[k], len(p), retmod=False)

   # analyse the CCF peak fitting
   v, e_v, a = SSRstat(vgrid, SSR, plot=plot)

   if np.isnan(e_v):
      v = vgrid[nk//2]   # actually it should be nan, but may the next clipping loop or plot use vcen
      print(" Setting  v=" % v)
   if vfix: v = 0.
   p, SSRmin, fmod = polyreg(x2, y2, e_y2, v, len(p))   # final call with v
   if p[0] < 0:
      e_v = np.nan
      print("Negative scale value. Setting  e_v= %f" % e_v)

   if plot&1 or (plot&2 and np.isnan(e_v)):
      gplot(x2, y2, fmod, ' w lp, "" us 1:3 w lp lt 3')
      pause(v)
   return type('par', (), {'params': np.append(v,p), 'perror': np.array([e_v,1.0]), 'ssr': (vgrid,SSR)}), fmod


def optivsini(va, vb, v_step, vtpl, x2, y2, e_y2, tpl, p=None, plot=2):
   """
   performs a mini CCF; the grid stepping
   returns best vsini and errors from parabola curvature
   """
   warnings.filterwarnings("ignore", category=DeprecationWarning) 
   vgrid = np.arange(va, vb, v_step)
   nk = len(vgrid)

   SSR = np.empty(nk)

   # set up wcen
   calcspec.wcen = np.mean(x2)

   for k in range(nk):

      # set up template
      calcspec.tpl = copy.deepcopy(tpl)

      # broaden template
      calcspec.tpl.rotbroad(vgrid[k])

      # fitting
      p, SSR[k] = polyreg(x2, y2, e_y2, vtpl, len(p), retmod=False)

   # analyse the CCF peak fitting
   v, e_v, a = SSRstat(vgrid, SSR, plot=plot)

   if np.isnan(e_v):
      v = vgrid[nk//2]   # actually it should be nan, but may the next clipping loop or plot use vcen
      print(" Setting  v=" % v)
   if p[0] < 0:
      e_v = np.nan
      print("Negative scale value. Setting  e_v= %f" % e_v)

   if np.argmin(SSR) == 0:
      # if the minimum is at zero, we cannot use parabola
      # take vstep as uncertainty
      v = 0
      e_v = v_step


   # template with fitted vsini
   calcspec.tpl = copy.deepcopy(tpl)

   # broaden template
   calcspec.tpl.rotbroad(v)

   # final call with v
   p, SSRmin, fmod = polyreg(x2, y2, e_y2, vtpl, len(p))

   if plot&1 or (plot&2 and np.isnan(e_v)):
      gplot(x2, y2, fmod, ' w lp, "" us 1:3 w lp lt 3')
      pause(v)
   return type('par', (), {'params': np.append(v,p), 'perror': np.array([e_v,1.0]), 'ssr': (vgrid,SSR)}), fmod


def fitspec(tpl, w, f, e_f=None, v=0, vfix=False, clip=None, nclip=1, keep=None, indmod=np.s_[:], v_step=True, df=None, plot=2, deg=3, chi2map=False):
   """
   Performs the robust least square fit via iterative clipping.

   vfix : boolean, optional
       v is fixed. For RV constant stars or in coadding when only the background polynomial is computed.
   indmod : Index range to be finally calculated-
   clip : Kappa sigma clipping value.
   nclip : Number of clipping iterations (default: 0 if clip else 1).
   df : Derivative for drift measurement.
   v_step : Number of v steps (only background polynomial => v_step = false).

   """
   if keep is None: keep = np.arange(len(w))
   if clip is None: nclip = 0   # number of clip iterations; default 1
   calcspec.wcen = np.mean(w[keep])
   calcspec.tpl = tpl

   p = np.array([v, 1.] + [0]*deg)   # => [v,1,0,0,0]
   fMod = np.nan * w
   #fres = 0.*w     # same size
   for n in range(nclip+1):
      if df is not None:
         '''drift mode: scale derivative to residuals'''
         par, fModkeep = optidrift(tpl[1].take(keep,mode='clip'), df.take(keep,mode='clip'), f.take(keep,mode='clip'),
                                 e_f.take(keep,mode='clip'))
      elif v_step:
         '''least square mode'''
         par, fModkeep = opti(v+v_lo, v+v_hi, w.take(keep,mode='clip'), f.take(keep,mode='clip'),
                              e_f.take(keep,mode='clip'), p[1:], vfix=vfix, plot=plot)
         ssr = par.ssr
      else:
         '''only background polynomial'''
         p, SSR, fModkeep = polyreg(w.take(keep,mode='clip'), f.take(keep,mode='clip'), e_f.take(keep,mode='clip'), v, len(p)-1)
         par = type('par',(),{'params': np.append(v,p), 'ssr': SSR})

      p = par.params
      par.niter = n
      if 1:
         # exclude model regions with negative flux
         # can occur e.g. in background in ThAr, or due to tellurics
         ind = fModkeep > 0
         keep = keep[ind]      # ignore the pixels modelled with negative flux
         fModkeep = fModkeep[ind]
         # all might be negative (for low/zero S/N data)
      #fres[keep] = (f[keep] - fMod[keep]) / fMod[keep]**0.5
      #res_std = rms(fres[keep])     # residual noise / photon noise
      # use error predicted by model
      #fres = (f.take(keep,mode='clip')-fModkeep) / fModkeep**0.5
      # use external errors
      fres = (f.take(keep,mode='clip')-fModkeep) / e_f.take(keep,mode='clip')
      res_std = rms(fres)     # residual noise / photon noise
      if n < nclip:
         ind = np.abs(fres) <= clip*res_std
         nreject = len(keep) - np.sum(ind)
         if nreject: keep = keep[ind]   # prepare next clip loop
         # else: break
      if len(keep)<10: # too much rejected? too many negative values?
         print("too much rejected, skipping")
         break
      if 0 and np.abs(par.params[0]*1000)>20:
         if df:
            fMod = ft * p[1]     # compute also at bad pixels
         else:
            fMod = calcspec(w, *p)     # compute also at bad pixels
         gplot.y2tics().ytics('nomir; set y2range [-5:35];')
         gplot(w,fMod,' w lp pt 7 ps 0.5 t "fmod"')
         gplot+(w[keep],fMod[keep],' w lp pt 7 ps 0.5 t "fmod[keep]"')
         gplot+(w,f,' w lp pt 7 ps 0.5 t "f"')
         #ogplot(w[keep],fres[keep],' w lp pt 7 ps 0.5 lc rgb "red" axis x1y2 t "residuals"',flush='')
         ogplot(w[keep],fres,' w lp pt 7 ps 0.5 lc rgb "black" axis x1y2, 0 w l lt 2 axis x1y2 t"", '+str(res_std)+' w l lt 1 axis x1y2, '+str(-res_std)+ ' w l lt 1 t "" axis x1y2')
         pause('large RV', par.params[0]*1000)

   stat = {"std": res_std, 'ssrmin': fres.sum(), "snr": np.mean(fModkeep)/np.mean(np.abs(f.take(keep,mode='clip')-fModkeep))}
   #pause(stat["snr"], wmean(fModkeep)/wrms(f.take(keep,mode='clip')-fModkeep), np.median(fModkeep)/np.median(np.abs(f.take(keep,mode='clip')-fModkeep)) )
   if df is not None:
      fMod[indmod] = tpl[1][indmod]*p[1] - df[indmod]*p[1]*p[0]/c  # compute also at bad pixels
   else:
      fMod[indmod] = calcspec(w[indmod], *p)   # compute also at bad pixels

   if chi2map:
      return par, fMod, keep, stat, ssr
   else:
      return par, fMod, keep, stat



def serval():

   if not bp: sys.stdout = Logger()

   global obj, targ, oset, coset, last, tpl, sp, fmod, reana, inst, fib, look, looki, lookt, lookp, lookssr, lookvsini, pmin, pmax, debug, pspllam, kapsig, nclip, atmfile, skyfile, atmwgt, omin, omax, ptmin, ptmax, driftref, deg, targrv, tplrv, tplvsini

   outdir = obj + '/'
   fibsuf = '_B' if inst=='FEROS' and fib=='B' else ''

   print(description)

   ### SELECT INSTRUMENTAL FORMAT ###
   # general default values
   pat = getattr(inst, 'pat', '*tar')  # default search pattern
   maskfile = servallib + getattr(inst, 'maskfile', 'telluric_mask_atlas_short.dat')

   # instrument specific default values
   iomax = inst.iomax

   if inst.name in ['CARM_VIS', 'CARM_NIR', 'FEROS']:
      if fib == '': fib = 'A'
   elif 'HARP' in inst.name:
      if fib == '': fib = 'A'
      if fib == 'B': iomax = 71
   elif inst.name == 'FEROS':
      iomax = 38
      if fib == 'B': maskfile = servallib + 'feros_mask_short.dat'
      ptomin = np.array([1800, 2000, 1800, 2000, 2000, 1600, 1500, 1400, 1100, 1000,
                         1000, 1000, 1000,  900,  800,  800,  800,  600,  500,  500,
                          400,  400,  300,  100,  100,  100,  100,  100,  100,  100,
                          100,  100,  100,  100,  100,  100,  100,  100,  100])
      ptomax = np.array([3100, 4000, 4000, 4000, 4000, 4500, 4600, 4600, 4800, 4800,
                         5000, 5200, 5200, 5300, 5600, 5700, 5800, 6000, 6100, 6300,
                         6500, 6600, 6700, 6800, 7100, 7200, 7400, 7700, 7900, 8100,
                         8400, 8600, 8900, 9200, 9500, 9800,10200,10600,11100])
      pomin = ptomin + 300
      pomax = ptomax - 300
      pmin = pomin.min()
      pmax = pomin.max()
   elif inst.name == 'FTS':
      iomax = 70
      pmin = 300
      pmax = 50000/5 - 500

   pat = pat % {'fib': fib}

   ptmin = pmin - 100   # oversize the template
   ptmax = pmax + 100

   orders = np.arange(iomax)[oset]
   corders = np.arange(iomax)[coset]

   orders = sorted(set(orders) - set(o_excl))
   corders = sorted(set(corders) - set(co_excl))
   omin = min(orders)
   omax = max(orders)
   comin = min(corders)
   comax = max(corders)

   if reana:
      x = analyse_rv(obj, postiter=postiter, fibsuf=fibsuf, oidx=orders)
      exit()

   os.system('mkdir -p '+obj)

   t0 = time.time()
   print("start time:", time.strftime("%Y-%m-%d %H:%M:%S %a", time.localtime()))

   if fib == 'B' or (ccf is not None and 'th_mask' in ccf):
      pass


   ### SELECT TARGET ###
   if tplrv=='auto' and (tpl is None or isinstance(tpl, int)):  tplrv = targrv
   if targrv=='auto' and (tpl is None or isinstance(tpl, int)): targrv = tplrv
   try:
      targrv_src, targrv = 'user', float(targrv) # is a float, so user input
   except:
      targrv_src, targrv = targrv, None
   try:
      tplrv_src, tplrv = 'usertpl', float(tplrv)
   except:
      tplrv_src, tplrv = tplrv, None


   targ = type('TARG', (), dict(name=targ, plx=targplx, rv=targrv, sa=float('nan')))
   targ.ra, targ.de = targrade
   targ.pmra, targ.pmde = targpm

   if targ.name == 'cal':
      print('no barycentric correction (calibration)')
   elif targ.ra and targ.de or targ.name:
      targ = Targ(targ.name, targrade, targpm, plx=targplx, rv=targrv, csv=obj+'/'+obj+'.targ.csv')
      print(' using sa=', targ.sa, 'm/s/yr', 'ra=', targ.ra, 'de=', targ.de, 'pmra=', targ.pmra, 'pmde=', targ.pmde)
   else:
      print('using barycentric correction from DRS')

   # choose the interpolation type
   spltype = 3 # 3=> fast version
   spline_cv = {1: interpolate.splrep, 2: cubicSpline.spl_c,  3: cubicSpline.spl_cf }[spltype]
   spline_ev = {1: interpolate.splev,  2: cubicSpline.spl_ev, 3: cubicSpline.spl_evf}[spltype]

   print(dir_or_inputlist)
   print('tpl=%s pmin=%s nset=%s omin=%s omax=%s' % (tpl, pmin, nset, omin, omax))

   ''' SELECT FILES '''
   files = []
   for patk in pat.split():   # search with multiple suffices
       files += glob.glob(dir_or_inputlist + os.sep + patk)
   files = sorted(files)

   isfifo = not files and (os_stat.S_ISFIFO(os.stat(dir_or_inputlist).st_mode))
   if os.path.isfile(dir_or_inputlist) or isfifo:
      if dir_or_inputlist.endswith(('.txt', '.lis')) or isfifo:
         with open(dir_or_inputlist) as f:
            print('getting filenames from file (',dir_or_inputlist,'):')
            for line in f:
               line = line.split()   # remove comments
               if line:
                  if os.path.isfile(line[0]):
                     files += [line[0]]
                     print(len(files), files[-1])
                  else:
                     print('skipping:', line[0])
         if fib == 'B':
           files = [f.replace('_A.fits','_B.fits') for f in files]
           print('renaming', files)
      else:
         # case if dir_or_inputlist is one fits-file
         files = [dir_or_inputlist]

   # handle tar, e2ds, fox
   drs = bool(len(files))
   if 'HARPS' in inst.name and not drs:  # fox
      files = sorted(glob.glob(dir_or_inputlist+'/*[0-9]_'+fib+'.fits'))
   if 'CARM' in inst.name and not files:
      files = sorted(glob.glob(dir_or_inputlist+'/*pho*_'+fib+'.fits'))
   if 'FEROS' in inst.name:
      files += sorted(glob.glob(dir_or_inputlist+'/f*'+('1' if fib=='A' else '2')+'0001.mt'))
   if 'FTS' in inst.name:
      files = sorted(glob.glob(dir_or_inputlist+'/*ap08.*_ScSm.txt'))
      files = [s for s in files if '20_ap08.1_ScSm.txt' not in s and '20_ap08.2_ScSm.txt' not in s and '001_08_ap08.193_ScSm.txt' not in s ]

   files = np.array(files)[nset]
   nspec = len(files)
   if not nspec:
      print("no spectra found in", dir_or_inputlist, 'or using ', pat, inst)
      exit()
   # expand slices to index arrays
   if look: look = np.arange(iomax)[look]
   if lookt: lookt = np.arange(iomax)[lookt]
   if lookp: lookp = np.arange(iomax)[lookp]
   if lookssr: lookssr = np.arange(iomax)[lookssr]
   if lookvsini: lookvsini = np.arange(iomax)[lookvsini]

   if outfmt or outchi: os.system('mkdir -p '+obj+'/res')
   with open(outdir+'lastcmd.txt', 'w') as f:
      print(' '.join(sys.argv), file=f)
   with open('cmdhistory.txt', 'a') as f:
      print(' '.join(sys.argv), file=f)

   badfile = open(outdir + obj + '.flagdrs' + fibsuf + '.dat', 'w')
   infofile = open(outdir + obj + '.info' + fibsuf + '.csv', 'w')
   bervfile = open(outdir + obj + '.brv' + fibsuf + '.dat', 'w')
   prefile = outdir + obj + '.pre' + fibsuf + '.dat'
   rvofile = outdir + obj + '.rvo' + fibsuf + '.dat'
   e_rvofile = outdir + obj + '.e_rvo' + fibsuf + '.dat'
   snrfile = outdir + obj + '.snr' + fibsuf + '.dat'
   chifile = outdir + obj + '.chi' + fibsuf + '.dat'
   halfile = outdir + obj + '.halpha' + fibsuf + '.dat'
   nadfile = outdir + obj + '.nad' + fibsuf + '.dat'
   irtfile = outdir + obj + '.cairt' + fibsuf + '.dat'
   dlwfile = outdir + obj + '.dlw' + fibsuf + '.dat'
   e_dlwfile = outdir + obj + '.e_dlw' + fibsuf + '.dat'

   # (echo 0 0 ; awk '{if($2!=x2){print x; print $0}; x=$0; x2=$2;}' telluric_mask_atlas.dat )> telluric_mask_atlas_short.dat
   #################################
   ### Loading and prepare masks ###
   #################################
   mask = None
   tellmask = nomask
   skymsk = nomask

   if fib == 'B':
      if atmfile == 'auto':
         atmfile = None
      if skyfile == 'auto':
         skyfile = None

   if ccf and 'th_mask' not in ccf:
      atmfile = None

   if atmspec:
      # standard telluric spectrum
      if inst.name != 'CARM_VIS':
         pause('Only implemented for CARM_VIS')
      import astropy.io.fits as pyfits
      #hdu = pyfits.open('/home/astro115/carmenes/tellurics/stdatmos_vis/stdatmos_vis30a090rh0780p000t.fits')
      hdu = pyfits.open(atmspec)
      wt, ft = lam2wave(hdu[1].data.field(0).astype(np.float64)), hdu[1].data.field(1).astype(np.float64)
      atmmod = spl.ucbspl_fit(wt, ft, K=int(ft.size/2))
      tidx = slice(None, -5)

   if atmfile:
      if atmfile != 'auto':
         maskfile = atmfile
      if 'mask_ne' in atmfile:
         maskfile = servallib + atmfile

      print('maskfile', maskfile)
      mask = np.genfromtxt(maskfile)

      if 'telluric_mask_atlas_short.dat' in maskfile:
         lcorr = 0.000009  # Guillems mask needs this shift of 2.7 km/s
         mask[:,0] = airtovac(mask[:,0]) * (1-lcorr)
      if 'th_mask' in maskfile: # well, it is not atmosphere, but ...
         mask[:,1] = mask[:,1] == 0  # invert mask; actually the telluric mask should be inverted (so that 1 means flagged and bad)

   if skyfile:
      if skyfile=='auto' and inst.name=='CARM_NIR':
         skyfile = servallib + 'sky_carm_nir'
         sky = np.genfromtxt(skyfile)
         skymsk = interp(lam2wave(sky[:,0]), sky[:,1])


   msksky = [0] * iomax
   if 1 and inst.name=='CARM_VIS':
      import astropy.io.fits as pyfits
      msksky = flag.sky * pyfits.getdata(servallib + 'carm_vis_tel_sky.fits')

   if msklist: # convert line list to mask
      mask = masktools.list2mask(msklist, wd=mskwd)
      mask[:,1] = mask[:,1] == 0  # invert mask; actually the telluric mask should be inverted (so that 1 means good)

   if mask is None:
      print('using telluric mask: NONE')
   else:
      # make the discrete mask to a continuous mask via linear interpolation
      tellmask = interp(lam2wave(mask[:,0]), mask[:,1])
      print('using telluric mask: ', maskfile)

   if 0:
      mask2 = np.genfromtxt('telluric_add.dat')
      # DO YOU NEED THIS: mask2[:,0] = airtovac(mask2[:,0])  ??
      i0 = 0 #where(mask[:,0]<mask2[0][0])[0][-1]
      for ran in mask2:
         while (mask[i0,0]<ran[0]): i0 += 1
         #if mask[i0-1,0]==0:
         mask = np.insert(mask,i0,[ran[0]-0.0000001,0.0],axis=0) # insert
         mask = np.insert(mask,i0+1,[ran[0],1.0],axis=0) # insert
         #while (mask[i0,0]<ran[1]): np.delete(mask,i0,axis=0)
         #if mask[i0-1,0]==0:
         mask = np.insert(mask,i0+2,[ran[1],1.0],axis=0) # insert
         mask = np.insert(mask,i0+3,[ran[1]+0.0000001,0.0],axis=0) # insert


   ################################
   ### READ FITS FILES ############
   ################################
   if nspec > 500:
      ulimit = resource.getrlimit(resource.RLIMIT_OFILE)
      if ulimit[0] < 4096:
         print("Many files! Adapting ulimit from %s to (4096, 4096)." % (ulimit,))
         resource.setrlimit(resource.RLIMIT_OFILE, (4096,4096))

   splist = Vectorise([])
   spi = None
   SN55best = 0.
   print("    # %*s %*s OBJECT    BJD        SN  DRSBERV  DRSdrift flag calmode" % (-len(inst.name)-6, "inst_mode", -len(os.path.basename(files[0])), "timeid"))

   for n,filename in enumerate(files):   # scanning fitsheader
      print('%3i/%i ' % (n+1, nspec), end='')
      sp = Spectrum(filename, inst=inst, pfits=2 if 'HARPS' in inst.name else True, drs=drs, fib=fib, targ=targ, verb=True)
      splist.append(sp)
      sp.sa = targ.sa / 365.25 * (sp.bjd-splist[0].bjd)
      sp.header = None   # saves memory(?), but needs re-read (?)
      if inst.name in ('HARPS', 'HARPN') and drs: sp.ccf = read_harps_ccf(filename)
      if any(x in filename for x in n_excl):
          sp.flag |= sflag.user
      if sp.sn55 < snmin or np.isnan(sp.sn55): sp.flag |= sflag.lowSN
      if sp.sn55 > snmax or np.isnan(sp.sn55): sp.flag |= sflag.hiSN
      if distmax and sp.ra and sp.de:
         # check distance for mis-pointings
         # yet no proper motion included
         ra = (targ.ra[0] + (targ.ra[1]/60 + targ.ra[2]/3600)* np.sign(targ.ra[0]))*15   # [deg]
         de = targ.de[0] + (targ.de[1]/60 + targ.de[2]/3600)* np.sign(targ.de[0])   # [deg]
         dist = np.sqrt(((sp.ra-ra)*np.cos(np.deg2rad(sp.de)))**2 + (sp.de-de)**2) * 3600
         if dist > distmax: sp.flag |= sflag.dist
      if not sp.flag:
         if SN55best < sp.sn55 < snmax:
            SN55best = sp.sn55
            spi = n
      else:
         print(sp.bjd, sp.ccf.rvc, sp.ccf.err_rvc, sp.timeid, sp.flag, file=badfile)
      print(sp.bjd, sp.berv, sp.drsbjd, sp.drsberv, sp.drift, sp.timeid, sp.tmmean, sp.exptime, sp.berv_start, sp.berv_end, file=bervfile)

   badfile.close()
   bervfile.close()
   bp or sys.stdout.logname(obj+'/log.'+obj)

   t1 = time.time() - t0
   print(nspec, "spectra read (%s)\n" % minsec(t1))

   obsloc = getattr(inst, 'obsloc', {})
   if obsloc and targ.ra and moonsep:
       # The computation have a lot of overhead and are done in vectorised fashion.
       print('Calculating moon separation')
       try:
           from astropy.time import Time
           from astropy.coordinates import SkyCoord, EarthLocation, AltAz
           from astropy.coordinates import solar_system_ephemeris, get_moon, get_sun
           import astropy.units as u
           from astropy.utils.iers import conf as iers_conf
           iers_conf.iers_auto_url = 'https://datacenter.iers.org/data/9/finals2000A.all'

           dateobs = Time(splist.dateobs)
           loc = EarthLocation.from_geodetic(lat=obsloc['lat'], lon=obsloc['lon'], height=obsloc['elevation'])
           aa = AltAz(location=loc, obstime=dateobs)
           sc = SkyCoord("%i:%i:%s" %targ.ra, "%i:%i:%s" %targ.de, unit=(u.hourangle, u.deg), obstime=dateobs)
           with solar_system_ephemeris.set('builtin'):
               sun = get_sun(dateobs)
               moon = get_moon(dateobs, loc)

           splist.sunalt = sun.transform_to(aa).alt.deg
           splist.secz = sc.transform_to(aa).secz    # airmass estimate
           splist.moonsep = moon.separation(sc).deg
           splist.moonphase = moon.separation(sun).deg   # https://astroplan.readthedocs.io/en/latest/_modules/astroplan/moon.html
           #splist.moonillumination = splist.moonphase
           splist.flag |= sflag.daytime * (splist.sunalt > sunalt)
           splist.flag |= sflag.moon * (splist.moonsep < moonsep)
       except Exception as e:
            print(e)
            pause('Warning: Could not compute moon separation. Use "-moonsep 0" to skip.')

   # filter for the good spectra
   check_daytime = True
   spoklist = Vectorise([])
   for sp in splist:
      print(sp.timeid, sp.bjd, sp.berv, sp.sn55, sp.obj, sp.exptime, sp.ccf.mask, sp.flag, sp.airmass, sp.ra, sp.de, sp.sunalt, sp.moonsep, sp.moonphase, sep=';', file=infofile)
      if sp.flag & (sflag.nosci|sflag.config|sflag.iod|sflag.dist|sflag.lowSN|sflag.hiSN|sflag.led|sflag.user):
         print('bad spectrum:', sp.timeid, sp.obj, sp.calmode, 'sn: %s flag: %s %s' % (sp.sn55, sp.flag, sflag.translate(sp.flag)))
      else:
         if sp.flag & (check_daytime*sflag.daytime|sflag.moon):
             # we calcuate RVs, but dLW can be biased upwards (not used in coadding)
             print('poor spectrum:', sp.timeid, sp.obj, sp.calmode, 'sn: %s flag: %s %s' % (sp.sn55, sp.flag, sflag.translate(sp.flag)))
         spoklist += [sp]

   infofile.close()

   nspecok = len(spoklist)
   if not nspecok:
      print("WARNING: no good spectra")
      if not safemode: pause()   # ???

   snrmedian = np.median([sp.sn55 for sp in spoklist])
   with open(outdir+obj+'.drs.dat', 'w') as myunit:
      for sp in spoklist:
          print(sp.bjd, sp.ccf.rvc*1000., sp.ccf.err_rvc*1000., sp.ccf.fwhm, sp.ccf.bis, sp.ccf.contrast, sp.timeid, file=myunit)  #, sp.hdr['HIERARCH ESO OBS PROG ID'], sp.hdr['HIERARCH ESO OBS PI-COI NAME']

   ################################
   ### Template canditates ########
   ################################
   if spi is None:
      print('No highest S/N found; selecting first file as reference')
      spi = 0

   if last:
      tpl =  outdir + obj + '.tpl%s.fits' % fibsuf
   elif tpl is None:
      tpl = spi   # choose highest S/N spectrum

   if isinstance(tpl, int):
      spi = tpl
   spt = files[spi]   # splist[tpl]

   # In any case we want to read the highest S/N spectrum (e.g. @wfix).
   # Use pyfits to get the full header
   spt = Spectrum(spt, inst=inst, pfits=True, orders=np.s_[:], drs=drs, fib=fib, targ=targ)
   # Estimate a Q for each order to identify fast rotators
   # Estimation similar as in Bouchy01
   # yet tellurics are not masked
   Wi = (spt.f[:,2:]-spt.f[:,:-2]) / (spt.w[:,2:]-spt.w[:,:-2]) / spt.e[:,1:-1]   # Eq.(8)
   dv = c / np.sqrt(np.nansum(Wi**2, axis=1))   # theoretical RV precision Eq.(10)
   sn = np.sqrt(np.nansum((spt.f[:,1:-1] / spt.e[:,1:-1])**2, axis=1))   # total SNR over order, corresponds to sqrt(Ne) in Bouchy
   Q = c / dv / sn
   #Q = np.sqrt(np.nansum(Wi**2, axis=1)) / np.sqrt(np.nansum((spt.f[:,1:-1] / spt.e[:,1:-1])**2, axis=1) # a robust variant
   #gplot(Q)

   print('median SN:', snrmedian)
   print('template:', spt.timeid, 'SN55', spt.sn55, '#', spi, ' <e_rv>=%0.2fm/s, <Q>=%s' % (np.median(dv)*1000, np.median(Q)))

   # find the order where to measure the line indices
   line_o = {}
   for l in lines:   # Halpha, NaD, CaIRT
      line_o[l] = get_o_of_line(l, spt.w)

   nord = len(spt.w)

   RV = None
   if vtfix:
      # RVs are zero
      np.savetxt(prefile, np.array([(sp.bjd, 0., 0.) for sp in spoklist]))


   ################################
   ### Select high S_N template ###
   ################################
   print()
   ntpix = ptmax - ptmin
   pixxx = arange((ntpix-1)*4+1) / 4.
   osize = len(pixxx)
   ww = [0] * nord
   ff = [0] * nord


   if inst.name == 'FEROS':
      ntopix = ptomax - ptomin

   is_ech_tpl = True   # echelle or continuous spectrum
   TPL = [None] * nord
   TPLrv = None

   if 1:
      to = time.time()
      if ccf:
         ccfmask = np.loadtxt(servallib + ccf)
      elif driftref:
         print(driftref)
         spt = Spectrum(driftref, inst=inst, pfits=True, orders=np.s_[:], drs=drs, fib=fib, targ=targ)
         ww, ff = spt.w, spt.f
      elif isinstance(tpl, str):
         print("restoring template: ", tpl)
         try:
            if 'phoe' in tpl or 'PHOENIX-ACES-AGSS-COND' in tpl:
               # PHOENIX template must be downloaded into servallib
               ww, ff = phoenix_as_RVmodel.readphoenix(servallib + ('lte03900-5.00-0.0_carm.fits' if 'phoe' in tpl else tpl), wmin=np.exp(np.nanmin(spt.w)), wmax=np.exp(np.nanmax(spt.w)))
               ww = lam2wave(ww)
               is_ech_tpl = False
               TPL = [Tpl(ww, ff, spline_cv, spline_ev, vsini=tplvsini)] * nord
               ww = [lam2wave(ww)] * nord   # @vsiniauto
               ff = [ff] * nord
               TPLrv = 0.
            elif tpl.endswith('.tpl%s.fits'%fibsuf) or os.path.isdir(tpl):
               # last option
               # read a spectrum stored order wise
               print("tplvsini", tplvsini)
               ww, ff, head = read_template(tpl+(os.sep+os.path.basename(tpl.rstrip(os.sep))+'.tpl.fits' if os.path.isdir(tpl) else ''))
               TPL = [Tpl(wo, fo, spline_cv, spline_ev, vsini=tplvsini) for wo,fo in zip(ww,ff)]
               if 'HIERARCH SERVAL COADD NUM' in head:
                  print('HIERARCH SERVAL COADD NUM:', head['HIERARCH SERVAL COADD NUM'])
                  if omin<head['HIERARCH SERVAL COADD COMIN']: pause('omin to small')
                  if omax>head['HIERARCH SERVAL COADD COMAX']: pause('omax to large')
               TPLrv = head['HIERARCH SERVAL TARG RV']
            else:
               # a user specified other observation/star
               spt = Spectrum(tpl, inst=inst, pfits=True, orders=np.s_[:], drs=drs, fib=fib, targ=targ)
               #ww, ff = barshift(spt.w,spt.berv), spt.f
               TPL = [Tpl(wo, fo, spline_cv, spline_ev, mask=True, berv=spt.berv) for wo,fo in zip(barshift(spt.w,spt.berv),spt.f)]
               TPLrv = spt.ccf.rvc
         except:
            print('ERROR: could not read template:', tpl)
            exit()

         if inst.name == 'FEROS':
            www = [0] * len(ww)
            fff = [0] * len(ww)
            for o in range(len(ww)):
               ind = ww[o] > 0 # remove padded zeros
               if ind.any():
                  www[o] = ww[o][ind]
                  fff[o] = ff[o][ind]
            ww = www
            ff = fff
      else:
         '''set up a spline template from spt'''
         TPLrv = spt.ccf.rvc
         for o in sorted(set(orders) | set(corders)):
            if inst.name == 'FEROS':
               ptmin = ptomin[o]
               ptmax = ptomax[o]
               pixxx = arange((ntopix[o]-1)*4+1)/4.
               osize = len(pixxx)

            pixx, = where((np.isfinite(spt.w) & np.isfinite(spt.f) & np.isfinite(spt.e))[o,ptmin:ptmax])
            idx = pixx + ptmin
            #ind = spt.bpmap[o,idx] == 0  # let out zero errors, interpolate over

            # Smoothing with bspline
            # the number of knots is halfed.
            smod = spl.ucbspl_fit(barshift(spt.w[o,idx],spt.berv), spt.f[o,idx], K=int(idx.size*ofac/2), e_yk=True, lam=0.00001)

            # Conversion to cardinal spline and then to fast spline
            smod_spl = smod.to_spl()
            #gplot(smod.osamp(10),',',wwo, ffo, ',', wwo,spline_ev(ww5,kko))
            xk, kko = smod_spl.xk, smod_spl.a
            yk = smod_spl()
            kko = xk, yk, kko[1]/(xk[1:]-xk[:-1]), np.append(2*kko[2]/(xk[1:]-xk[:-1])**2, 0), kko[3]/(xk[1:]-xk[:-1])**3
            TPL[o] = Tpl(xk, yk, spline_cv, spline_ev, mask=True, berv=spt.berv)
            TPL[o].funcarg = kko # replace with original spline

            if 0 or o==-50:
               #gplot(ww[o],ff[o], ',', barshift(spt.w[o,ptmin:ptmax],spt.berv), spt.f[o,ptmin:ptmax])
               gplot(barshift(spt.w[o,ptmin:ptmax],spt.berv), spt.f[o,ptmin:ptmax], spt.e[o,idx],'w e pt 7,', smod.osamp(10), 'w l lc 3')
               pause(o)

            #bb[o][ww[o]<=spt.w[o,idx[0]]] |= flag.nan   # flag egdes where nan flux might occur
            #bb[o][ww[o]>=spt.w[o,idx[-1]]] |= flag.nan

            #if inst .name== 'FEROS':
               # works with list so we do it here
               #bb[o][ff[o]<0] |= flag.neg

         #if inst.name != 'FEROS':
            #bb[ff<0] |= flag.neg
         #ind = (bb&flag.neg)==0
         #gplot(barshift(spt.w[o,ptmin:ptmax],spt.berv),spt.f[o,ptmin:ptmax])
         #ogplot(ww[o],ff[o]); pause()

   if vsiniauto:
      # if use automatic determination of vsini

      # save original template (without broadening)
      TPL0 = copy.deepcopy(TPL)

      # set up array for vsini
      VSINI = np.nan * np.empty([nord,2])


   rvdrs = np.array([sp.ccf.rvc for sp in spoklist])
   targrvs = {'simbad': targ.rv,
              'tpl': TPLrv,
              'drsspt': spt.ccf.rvc,
              'drsmed': np.median(rvdrs[np.isfinite(rvdrs)]),
              'user': targrv,
              'usertpl': tplrv
             }

   if targrv_src=='auto':
      if targ.name == 'cal':
         targrv_src = 'user'
         targrvs['user'] = 0
      elif np.isfinite(targrvs['drsspt']):
         targrv_src = 'drsspt'
      else:
         print('DRS RV is NaN in spt, trying median')
         if np.isfinite(targrvs['drsmed']):
            targrv_src = 'drsmed'
         elif targrvs['simbad'] is not None:
            print('DRS RV is NaN in all spec, simbad RV')
            targrv_src = 'simbad'
         else:
            print('DRS RV is NaN in all spec, simbad RV')
            targrv_src = 'usertpl'  # e.g. phoenix + online mode

   targrv = targrvs.get(targrv_src, 0)
   print('setting targ RV to: %s km/s (%s)' % (targrv, targrv_src))

   if targrv is None or np.isnan(targrv):
       raise ValueError(
          '\n  Absolute RV is not available from header or simbad or user.' +
          '\n  It is required to measure line indices.' +
          '\n  You could specify e.g. "-targrv 0", but re-check your estimate on some prominent lines e.g. with "-looki".'
       )

   if tplrv_src=='auto':
      if np.isfinite(TPLrv):
          # for external templates take value from fits header
          tplrv_src = 'tpl'
      else:  # e.g. simbad
          tplrv_src = targrv_src
   # no action for: tplrv_src=='usertpl'

   tplrv = targrvs.get(tplrv_src, 0)
   print('setting tpl RV to:  %s km/s (%s)' % (tplrv, tplrv_src))


   if skippre or vtfix:
      # restore the pre RVs
      if os.path.isfile(prefile):
         bjd, RV, e_RV = np.genfromtxt(prefile, unpack=True)
      else:
         pause('pre RV file', prefile, 'does not exist')



   for iterate in range(1, niter+1):

      print('\nIteration %s / %s (%s)' % (iterate, niter, obj))

      if RV is not None:
         ################################
         ### create high S_N template ###
         ################################
         tpl = outdir + obj + '.tpl' + fibsuf + '.fits'

         ntset = len(spoklist[tset])
         spt.header['HIERARCH SERVAL OFAC'] = (ofac, 'oversampling factor per raw pixel')
         spt.header['HIERARCH SERVAL PSPLLAM'] = (pspllam, 'smoothing value of the psline')
         spt.header['HIERARCH SERVAL UTC'] = (datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"), 'time of coadding')

         wk = [[np.nan]] * nord
         fk = [[np.nan]] * nord
         ek = [[np.nan]] * nord
         bk = [[0]] * nord

         for o in corders:
            print("coadding o %02i: " % o, end='')     # continued below in iteration loop
            npix = len(spt.w[o])
            osize = len(spt.w[o][ptmin:ptmax]) - 1   # -1 for compatibility with previous version
            nk = int(osize * ofac)

            wmod = nans((ntset,npix))
            mod = zeros((ntset,npix))
            emod = zeros((ntset,npix))
            bmod = zeros((ntset,npix), dtype=int)
            for n,sp in enumerate(spoklist[tset]):
             '''get the polynomials'''
             if not sp.flag:  # every flagged spectrum is excluded; more stringent than in RV loop
               if cache:
                   sp.read_data()
               sp = sp.get_data(pfits=2, orders=o)
               if atmspec:
                  ft = atmmod(sp.w)
                  sp.f = sp.f / ft
                  sp.e = sp.e / ft

               if inst.name == 'FEROS':
                  thisnpix = len(sp.bpmap)
                  # spectra can have different size
                  if thisnpix < npix:
                     sp.bpmap = np.append(sp.bpmap, np.ones(npix-thisnpix))
                     sp.e = np.append(sp.e, np.ones(npix-thisnpix))
                     sp.w = np.append(sp.w, sp.w[-1]*np.ones(npix-thisnpix))
                     sp.f = np.append(sp.f, np.ones(npix-thisnpix))
                  if thisnpix > npix:
                     sp.bpmap = sp.bpmap[:npix-thisnpix]
                     sp.e = sp.e[:npix-thisnpix]
                     sp.w = sp.w[:npix-thisnpix]
                     sp.f = sp.f[:npix-thisnpix]

               bmod[n] = sp.bpmap | msksky[o]
               bmod[n][tellmask(sp.w)>0.01] |= flag.atm
               bmod[n][skymsk(sp.w)>0.01] |= flag.sky
               # see https://github.com/mzechmeister/serval/issues/19#issuecomment-452661455
               # note in this step the RVs have reverted signs.
               bmod[n][tellmask(dopshift(redshift(sp.w, vo=sp.berv, ve=RV[n]/1000.), spt.berv))>0.01] |= flag.badT
               bmod[n][skymsk(dopshift(redshift(sp.w, vo=sp.berv, ve=RV[n]/1000.), spt.berv))>0.01] |= flag.badT

               w2 = redshift(sp.w, vo=sp.berv, ve=RV[n]/1000.)   # correct also for stellar rv
               i0 = np.searchsorted(w2, TPL[o].wk[0]) - 1   # w2 must be oversized
               if i0<0:
                  i0 = 0
               ie = np.searchsorted(w2, TPL[o].wk[-1])
               pind, = where(bmod[n][i0:ie] == 0)
               bmod[n][:i0] |= flag.out
               bmod[n][ie:] |= flag.out
               if np.sum(sp.f[pind]<0) > 0.4*pind.size:
                  print('too many negative data points in n=%s, o=%s, RV=%s; skipping order' % (n, o, RV[n]))
                  wmod[n] = np.nan
                  mod[n] = np.nan
                  emod[n] = np.nan
                  continue

               if inst.name == 'FEROS':
                   uind = pind*1
                   hh = np.argsort(sp.f[i0:ie])
                   ii = hh[0:len(hh)*0.98]
                   pind = np.intersect1d(pind, ii)
                   #gplot(sp.w[i0:ie], sp.f[i0:ie],',',sp.w[i0:ie][uind],sp.f[i0:ie][uind],',',sp.w[i0:ie][pind], sp.f[i0:ie][pind])
               if not len(pind):
                  print('no valid points in n=%s, o=%s, RV=%s; skipping order' % (n, o, RV[n]))
                  if not safemode: pause()
                  break

               # get poly from fit with mean RV
               par, fmod, keep, stat = fitspec(TPL[o],
                  w2[i0:ie], sp.f[i0:ie], sp.e[i0:ie], v=0, vfix=True, keep=pind, v_step=False, clip=kapsig, nclip=nclip, deg=deg, plot=0)   # RV  (in dopshift instead of v=RV; easier masking?)
               poly = calcspec(w2, *par.params, retpoly=True)
               wmod[n] = w2
               mod[n] = sp.f / poly   # be careful if  poly<0
               emod[n] = sp.e / poly
               if 0:# o in lookt: #o==-29:
                  gplot-(w2, sp.f, poly, ' t "sp.f", "" us 1:3 w l t "poly"')
                  gplot<(w2[i0:ie][keep], sp.f[i0:ie][keep], sp.e[i0:ie][keep], 'w e')
                  gplot+(w2[i0:ie], (sp.f / poly)[i0:ie], ' w l t "spf.f/poly",', TPL[o].wk, TPL[o].fk, ' w l t "TPL"')
                  pause(n)
              #(fmod<0) * flag.neg
            bmod &= ~flag.badT   # take out badT flag from observation, needed only in normalisation
            ind = (bmod&(flag.nan+flag.neg+flag.out)) == 0 # not valid
            tellind = (bmod&(flag.atm+flag.sky)) > 0                  # valid but down weighted
            #emod[tellind] *= 1000
            ind *= emod > 0.0
            we = 0*mod
            we[ind] = 1. / emod[ind]**2
            if atmfile and ('UNe' in atmfile or 'UAr' in atmfile or 'ThNe' in atmfile or 'ThAr' in atmfile): # old downweight scheme
               we[ind] += 0.000000001
               we[tellind] /= 10       # downweight
            elif atmwgt: # down weight with a constant factor
               # for low SN spectra or variing SN between observation, e.g. Trappist-1
               we[tellind] *= atmwgt   # downweight
            elif 0: # down weight with line depth
               # for high SN spectra, deep absorption line
               #we[tellind] *= (mod[tellind]/np.percentile(mod[ind+~tellind],95)).clip(0.02,1)**4 / 10
               fcont = np.abs(np.percentile(mod[ind&~tellind],95)*1.1)
               #fcont = quantile(mod[ind&~tellind], 0.95, w=1/emod[ind&~tellind])*1.1
               print(fcont)
               #we[tellind] = 1/(5*emod[tellind]**2 + (mod[tellind]-fcont)**2)
               #we[tellind] = 0.1/emod[tellind]**2   # for low S/N # keeps the relative S/N properties of the data
               #we[tellind] = 1/(emod[tellind]**2 + (fcont*np.log(abs(mod[tellind]/fcont).clip(1e-6)))**2) # for high S/N
               we[tellind] = 1/(emod[tellind]**2 + (fcont*np.log((np.sqrt(mod[tellind]**2+emod[tellind]**2)/fcont).clip(1e-6)))**2) # for high S/N
               #we[tellind] *= (mod[tellind]/fcont).clip(0.02,1)**4 / 10
               # old error vs new
               #gplot(wmod[ind],mod[ind], 1/np.sqrt(we[ind]), emod[ind], 'us 1:2:3 w e, "" us 1:2:4 w e')
            else:
               we[tellind] = 0.1 / ntset / (emod[tellind]**2 + np.median(emod[ind])**2)
            ind0 = ind*1

            n_iter = 2
            if inst.name == 'FEROS': n_iter = 3
            for it in range(n_iter+1):   # clip 5 sigma outliers
               #ind2 = ww[o] < wmod[ind].max()  # but then the complement for ind2 in ff is from spt!!!
                                               # maybe extrapolation better

               # B-spline fit for co-adding
               #if o == 46:stop()
               mu, e_mu = None, None
               if pmu and pe_mu:
                  # use mean as an estimate for continuum of absoption spectrum
                  mu = wmean(mod[ind], w=we[ind]) if pmu is True else pmu
                  # deviation of mu should be as large or larger than mu
                  e_mu = pe_mu * mu  # 5 *mu

               smod, ymod = spl.ucbspl_fit(wmod[ind], mod[ind], we[ind], K=nk, lam=pspllam, mu=mu, e_mu=e_mu, e_yk=True, retfit=True)

               #yfit = ww[o]* 0 # np.nan
               #ind2 &= (ww[o]> smod.xmin) & (ww[o]< smod.xmax)
               #yfit[ind2] = smod(ww[o][ind2])
               ww[o], yfit = smod.osamp(1.000000001*(4*osize+1)/smod.K)  # store with (+1) for compatibility
               wko = smod.xk     # the knot positions
               fko = smod()      # the knot values
               eko = smod.e_yk   # the error estimates for knot values
               dko = smod.dk()   # ~second derivative at knots

               #gplot(ww[o][ind2], yfit, ',', wko, fko,',',ww[o][ind2], gg)

               # normalised residuals
               # including telluric and original values gives conservative limits
               res = (mod[ind]-ymod) / emod[ind]
               # original values gives conservative limits
               #res = (mod[ind]-ymod) * np.sqrt(we[ind])
               #sig = std(res)
               sig = std(res[~tellind[ind]])
               # iqr(res, sig=True) untested, will be slower (in normal mode) than std but fewer iterations
               #gplot(wmod[ind], res,', %s lt 3, %s lt 3, %s lt 2, %s lt 2' %(-sig,sig,-ckappa[0]*sig, ckappa[1]*sig))

               if np.isnan(sig):
                  msg ='nan err_values in coadding. This may happen when data have gaps e.g. due masking or bad pixel flaging. Try the -pspline option'
                  if safemode:
                     print(msg)
                     exit()
                  pause(msg)
                  gplot(wmod[ind], mod[ind], we[ind])
               if 1:
                  # flexible sig
                  # sig = np.sqrt(spl._ucbspl_fit(wmod[ind], res**2, K=nk/5)[0])  # not good, van be negative
                  # get fraction of the data to each knot for weighting
                  G, kkk = spl._cbspline_Bk(wmod[ind], nk//5)
                  chik = np.zeros(nk//5+2)   # chi2 per knot
                  normk = np.zeros(nk//5+2)  # normalising factor to compute local chi2_red
                  for k in range(4):
                     normk += np.bincount(kkk+k, G[k], nk//5+2)
                     chik += np.bincount(kkk+k, res**2 * G[k], nk//5+2)

                  vara = spl.ucbspl(chik/normk, wmod[ind].min(), wmod[ind].max())
                  sig = np.sqrt(vara(wmod[ind]))

                  if 0 and o==33:
                     # show adaptive clipping threshold
                     gplot(wmod[ind], res, -sig, sig,  -sig*ckappa[0], sig*ckappa[1], ', "" us 1:3 w l lt 3, "" us 1:4 w l lt 3, ""  us 1:5 w l lt 7, ""  us 1:6 w l lt 7')
                     pause('look ',o)
               # pause('look ',o)

               okmap = np.array([True] * len(res))
               if ckappa[0]: okmap *= res > -ckappa[0]*sig
               if ckappa[1]: okmap *= res < ckappa[1]*sig
               okmap[tellind[ind]] = True # Oh my god. Do not reject the tellurics based on emod. That likely gives gaps and then nans.

               print("%.5f (%d) " % (np.median(sig), np.sum(~okmap)), end='')
               #gplot(wmod[ind], res,',', wmod[ind][tellind[ind]], res[tellind[ind]])
               #pause()
               if it < n_iter: ind[ind] *=  okmap

            if ofacauto:
               # BIC to get optimal knot spacing (smoothing)
               chired = []
               BIC = []
               K = np.logspace(np.log10(10), np.log10(ntpix), dtype=int)
               for Ki in K:
                  smod, ymod = spl.ucbspl_fit(wmod[ind], mod[ind], we[ind], K=Ki, lam=pspllam, mu=mu, e_mu=e_mu, e_yk=True, retfit=True)
                  chi = ((mod[ind] - ymod)**2*we[ind]).sum()
                  chired += [ chi / (we[ind].size-Ki)]
                  BIC += [ chi + np.log(we[ind].size)*Ki]

               Ko = K[np.argmin(BIC)]
               smod, ymod = spl.ucbspl_fit(wmod[ind], mod[ind], we[ind], K=Ko, lam=pspllam, mu=mu, e_mu=e_mu, e_yk=True, retfit=True)
               print("K=%d " % Ko, end='')

               if 0:
                  gplot2(K, BIC, 'w lp,', Ko, min(BIC), 'lc 3 pt 7')
                  gplot(wmod[ind], mod[ind], 'w d,', smod.osamp(10), 'w lp ps 0.3,', smod.xk, smod(), 'w p')
                  #pause()

            if vsiniauto:
               # vsini steps
               vs_hi = 150
               vs_step = 1

               # set up data for fitting
               ### cut 200 km/s at the edges, buffer for rotbroadening
               okmap = np.where((bmod==0) & (wmod > TPL0[o].wk[0]+200/c) & (wmod < TPL0[o].wk[-1]-200/c))
               ### sort by wavelength (for spl_evf in calcspec)
               sind = np.argsort(wmod[okmap])

               x = wmod[okmap][sind]
               y = mod[okmap][sind]
               yerr = emod[okmap][sind]

               # fit the template
               par, fModkeep = optivsini(0, vs_hi, vs_step, 0, #tplrv,
                                       x, y, yerr,
                                       TPL0[o],
                                       [1,0,0,0], plot=False)

               # chi2, vsini
               ssr = par.ssr
               vsini = par.params[0]
               e_vsini = par.perror[0]
               VSINI[o] = [vsini, e_vsini]

               # return best fitting vsini
               if not np.isnan(e_vsini):
                  print("\nvsini = %0.5f +/- %0.5f km/s (order %i)\n" % (vsini,e_vsini,o), end='')

               # plotting gplot
               if o in lookvsini:
                  gplot.xlabel('"ln(wavelength)"').ylabel('"flux"')
                  gplot(TPL0[o].wk, TPL0[o].fk, 'w l lc 2 t "tpl",',
                        x, y, 'lc 1 pt 1 ps 0.3 t "%s [o=%s]",'%(obj,o),
                        rotbroad(TPL0[o].wk, TPL0[o].fk, vsini), 'w l t "tpl (vsini=%.2f +/- %.2f km/s)"'% (vsini,e_vsini))
                  gplot2.xlabel('"vsini [km/s]"').ylabel('"SSR"')
                  gplot2(ssr, 'w lp t "",', vsini, min(ssr[1]), 'lc 3 pt 7 t "vsini = %.2f +/- %.2f km/s"' % (vsini,e_vsini))

                  pause()

            # estimate the number of valid points for each knot
            edges = 0.5 * (wko[1:]+wko[:-1])
            edges = np.hstack((edges[0]+2*(wko[0]-edges[0]), edges, edges[-1]+2*(wko[-1]-edges[-1])))
            bko,_ = np.histogram(wmod[ind], bins=edges, weights=(bmod[ind]==0)*1.0)

            '''estimate S/N for each spectrum and then combine all S/N'''
            sn = []
            yymod = mod * 0
            yymod[ind] = ymod
            for n,sp in enumerate(spoklist[tset]):
               if sp.sn55 < 400 and ind[n].any():
                  spt.header['HIERARCH COADD FILE %03i' % (n+1)] = (sp.timeid, 'rv = %0.5f km/s' % (-RV[n]/1000.))
                  iind = (n, ind[n])   # a short-cut for the indexing
                  signal = wmean(mod[iind], 1/emod[iind]**2)  # the signal
                  noise = wrms(mod[iind]-yymod[iind], emod[iind])   # the noise
                  sn.append(signal/noise)
            sn = np.sum(np.array(sn)**2)**0.5

            print(' S/N: %.5f' % sn)
            if not np.isnan(sn):
               spt.header['HIERARCH SERVAL COADD SN%03i' % o] = (float("%.3f" % sn), 'signal-to-noise estimate')
            if ofacauto:
               spt.header['HIERARCH SERVAL COADD K%03i' % o] = (Ko, 'optimal knot number')

            # plot the model and spt
            if o in lookt:
               gplot.bar(0).key("tit '%s order %s'"% (obj,o))
               gplot.mxtics().mytics().xlabel("'ln {/Symbol l}'")
               # estimate a good x2tics spacing
               i10, f10 = divmod(np.log10((np.exp(max(TPL[o].wk))-np.exp(min(TPL[o].wk)))/5), 1)
               dx2 = [1,2,5][abs(10**f10-[1,2,5]).argmin()] * 10**i10
               gplot.bind('f "dx2=(GPVAL_X2_MAX-GPVAL_X2_MIN)/5; i10=10**floor(log10(dx2)); dx2=dx2/i10; set x2tics i10*(dx2<1.5?1:dx2<4?2:5);        replot"')
               gplot.x2label("'{/Symbol l} [A]'").x2tics(dx2).mx2tics().link('x via exp(x) inverse log(x)').xtics("nomirr")
               #gplot(wmod[ind],mod[ind], 1/np.sqrt(we[ind]), emod[ind], 'us 1:2:3 w e lt 2, "" us 1:2:4  w e pt 7 ps 0.5 lt 1')
               #hasflag = lambda array,flag: (array&flag)==flag
               #hasflags = lambda array,flags: [hasflag(array,flag) for flag in flags]
               #bmod[ind<ind0] |= flag.clip
               #linecolor = np.array([0,1,6,4,4,5])[np.argmax([0*bmod, (bmod&~flag.badT)==0] + hasflags(bmod, [flag.atm, flag.out, flag.nan, flag.clip]), axis=0)]
               #gplot2.bar(0)(wmod.ravel(), mod.ravel(), emod.ravel(), linecolor.ravel(), ' us 1:2:3:4  w e pt 7 ps 0.5 lc var t "data"')
               gplot-(wmod[ind], mod[ind], emod[ind], ' w e pt 7 ps 0.5 t "data"')
               gplot<(TPL[o].wk, TPL[o].fk, 'us 1:2 w lp lt 2 ps 0.5 t "spt"')
               gplot<(ww[o], yfit,' us 1:2 w l lt 3 t "template"')
               if (~ind).any():
                  gplot<(wmod[~ind], mod[~ind], emod[~ind].clip(0,mod[ind].max()/20),'us 1:2:3 w e lt 4 pt 7 ps 0.5 t "flagged"')
               if (ind<ind0).any():
                  gplot<(wmod[ind<ind0], mod[ind<ind0], emod[ind<ind0].clip(0,1000),' w e lt 5 pt 7 ps 0.5 t "clipped"')
               if tellind.any():
                  gplot<(wmod[tellind],mod[tellind],' us 1:2 lt 6 pt 7 ps 0.5 t "atm"')
               gplot+(wko, fko, eko.clip(0,fko.max()), 'w e lt 3 pt 7 t "template knots"')
               if 0: # overplot normalised residuals
                  gplot_set('set y2tics; set ytics nomir; set y2range [-5*%f:35*%f]; set bar 0.5'%(sig,sig))
                  ogplot(wmod[ind], res,' w p pt 7 ps 0.5 lc rgb "black" axis x1y2')
               pause('lookt ',o)
               gplot.reset()
               # plot relative residuals
               if 0:
                  gplot(wmod[ind], res, -sig*ckappa[1], -sig, sig, sig*ckappa[1], 'us 1:2, "" us 1:3 w l lt 2, "" us 1:4 w l lt 3, "" us 1:5 w l lt 3, "" us 1:6 w l lt 2')
                  pause('lookt ',o)

            # apply the fit
            #ff[o][ind2] = yfit
            if ofacauto:
               # replace template.fits with optimal knot spacing (smoothing) for RVs
               yfit = ww[o] * 0 # np.nan
               ind2 = (ww[o]> smod.xmin) & (ww[o]< smod.xmax)
               yfit[ind2] = smod(ww[o][ind2])
               # pause()

            if not vsiniauto:
                ff[o] = yfit
                TPL[o] = Tpl(ww[o], ff[o], spline_cv, spline_ev)
                wk[o] = wko
                fk[o] = fko
                ek[o] = eko
                bk[o] = bko

         if vsiniauto:
            # apply fitted rotational broadening
            for o in orders:
               # update template with median vsini
               TPL[o].rotbroad(np.nanmedian(VSINI[:,0]))
            # do not store broadened phoenix template (is_ech_tpl=false)
         else:
            if isinstance(ff, np.ndarray) and np.isnan(ff.sum()): stop('nan in template')
            spt.header['HIERARCH SERVAL COADD OMIN'] = (omin, 'minimum order for RV')
            spt.header['HIERARCH SERVAL COADD OMAX'] = (omax, 'maximum order for RV')
            spt.header['HIERARCH SERVAL COADD COMIN'] = (comin, 'minimum coadded order')
            spt.header['HIERARCH SERVAL COADD COMAX'] = (comax, 'maximum coadded order')
            spt.header['HIERARCH SERVAL COADD NUM'] = (nspecok, 'number of spectra used for coadd')
            if targ.rv is not None and np.isfinite(targ.rv):
               spt.header['HIERARCH SERVAL TARG RV CSV'] = (targ.rv, '[km/s] RV from targ.csv')
            if np.isfinite(targrvs['drsspt']):
               spt.header['HIERARCH SERVAL TARG RV DRSSPT'] = (targrvs['drsspt'], '[km/s] DRS RV of spt')
            if np.isfinite(targrvs['drsmed']):
               spt.header['HIERARCH SERVAL TARG RV DRSMED'] = (targrvs['drsmed'], '[km/s] median DRS RV')
            spt.header['HIERARCH SERVAL TARG RV'] = (targrv, '[km/s] RV used')
            spt.header['HIERARCH SERVAL TARG RV SRC'] = (targrv_src, 'Origin of TARG RV')

            # Oversampled template
            write_template(tpl, ff, ww, spt.header, hdrref='', clobber=1)
            # Knot sampled template
            write_res(outdir+obj+'.fits', {'SPEC':fk, 'SIG':ek, 'WAVE':wk, 'NMAP':bk}, tfmt, spt.header, hdrref='', clobber=1)
            os.system("ln -sf " + os.path.basename(tpl) + " " + outdir + "template.fits")
            print('\ntemplate written to ', tpl)
            if 0: os.system("ds9 -mode pan '"+tpl+"[1]' -zoom to 0.08 8 "+tpl+"  -single &")

         # end of coadding
      print("time: %s\n" % minsec(time.time()-to))

      lstarmask = nomask # only 1D arrays
      do_reg = False
      if do_reg: reg = np.ones(ff.shape,dtype=bool)

      if ccf:
         pass
      elif is_ech_tpl:
         for o in orders:
            ii = np.isfinite(ff[o])
            if do_reg:
               print('q factor masking order ', o)
               reg[o] = qfacmask(ww[o],ff[o]) #, plot=True)

      if do_reg:
         idx = reg!=-1
         aa = ww[idx]
         bb = 1.0 * reg[idx]
         aaorg = 1.0 * aa.flatten()
         bborg = 1.0 * bb.flatten()
         idx = np.argsort(aa.flatten())
         aa = aa.flatten()[idx]
         bb = bb.flatten()[idx]

         #pos1, = where(bb)
         #for i, posi in enumerate(pos1[:-1]):
         #if (aa[pos1[i+1]]-aa[posi]) <(3./100000.): bb[posi:pos1[i+1]] = 1.0
         #idx1g = aa[idx1[1:]]-aa[idx1[:-1]]<3./(100000.)

         #compress mask
         def compress(aa,bb):
            '''keep only points with pre- or post-gradient'''
            idx = np.empty_like(bb,dtype=bool)
            idx[1:] = bb[1:]-bb[:-1] != 0      # pre-gradient
            idx[:-1] += bb[:-1]-bb[1:] != 0    # or post-gradient
            idx[[0,-1]] = True                 # keep always the endpoints
            return aa[idx],bb[idx]

         aa, bb = compress(aa, bb)

         # remove too small masked continuum regions (likely in order overlap)
         idx, = where(bb==1)
         dw = aa[idx[1:]]-aa[idx[:-1]]
         iidx = dw < (3./100000.)          # 3 resolution elements
         iidx = np.concatenate((idx[iidx]+1,idx[1:][iidx]-1))
         bb[iidx] = 1
         aa, bb = compress(aa, bb)
         #gplot(aaorg, bborg, 'w lp'); ogplot(aa,bb, 'w lp')
         #gplot(aa,bb); ogplot(aa[idx],bb[idx], 'w lp')
         #lstarmask = interp(ww.flatten()[idx], reg.flatten()[idx]*1.0) # only 1D arrays
         aa[0]=1 # extend the range
         aa[-1]=12
         lstarmask = interp(aa,bb) # only 1D arrays
         #lstarmask = interpolate.interp1d(aa, bb)


      ### Least square fitting
      results = dict((sp.timeid,['']*nord) for sp in splist)
      table = nans((7,nspec))
      bjd, RV, e_RV, rvm, rvmerr, RVc, e_RVc = table   # get pointer to the columns
      CRX, e_CRX = nans((2,nspec))   # get pointer to the columns
      mlRV, e_mlRV = nans((2,nspec))
      mlRVc, e_mlRVc = nans((2,nspec))
      mlCRX, e_mlCRX = nans((2,nspec))
      tCRX = np.rec.fromarrays(nans((5,nspec)), names='CRX,e_CRX,a,e_a,l_v' )   # get pointer to the columns
      xo = nans((nspec,nord))

      snr = nans((nspec,nord))
      rchi = nans((nspec,nord))
      Nok = nans((nspec,nord))


      meas_index = targrv_src and 'B' not in fib #and not 'th_mask' in ccf
      meas_CaIRT = meas_index and line_o['CaIRT1']
      meas_NaD = meas_index and line_o['NaD1']

      if meas_index:
         halpha = []
         haleft = []
         harigh = []
         cai = []
         cak = []
         cah = []
         irt1 = []
         irt1a = []
         irt1b = []
         irt2 = []
         irt2a = []
         irt2b = []
         irt3 = []
         irt3a = []
         irt3b = []
         nad1 = []
         nad2 = []
         nadr1 = []   # NaD Ref 1
         nadr2 = []
         nadr3 = []

      rvccf, e_rvccf = zeros((nspec,nord)), zeros((nspec,nord))
      diff_rv = bool(driftref)
      chi2map = [None] * nord
      #chi2map = nans((nord, int(np.ceil((v_hi-v_lo)/ v_step))))
      chi2map = nans((nord, len(np.arange(targrv-tplrv+v_lo, targrv-tplrv+v_hi, v_step))))
      diff_width = not (ccf or diff_rv)
      RV, e_RV = nans((2, nspec))
      rv, e_rv = nans((2, nspec, nord))
      dLW, e_dLW = nans((2, nspec)) # differential width change
      dlw, e_dlw = nans((2, nspec, nord)) # differential width change

      print('Iteration %s / %s (%s)' % (iterate, niter, obj))
      print("RV method: ", 'CCF' if ccf else 'DRIFT' if diff_rv else 'LEAST SQUARE')

      for n,sp in enumerate(spoklist):
         #if sp.flag:
            #continue
            # introduced for drift measurement? but then Halpha is not appended and writing halpha.dat will fail
         #sp = copy.deepcopy(sp)  # to prevent attaching the data to spoklist
         if not cache:
             sp = copy.copy(sp)  # to prevent attaching the data to spoklist
         if sp.filename.endswith('.gz') and sp.header:
            # for gz and if deepcopy and if file already open (after coadding header still present) this will result in "AttributeError: 'GzipFile' object has no attribute 'offset'"
            # deepcopy probably does not copy everything properly
            sp.header = None
         sp.read_data()
         if atmspec:
            # Divide out a standard atmosphere
            ft = atmmod(sp.w[tidx])
            sp.f[tidx] /= ft
            sp.e[tidx] /= ft
         bjd[n] = sp.bjd

         if wfix: sp.w = spt.w
         fmod = sp.w * np.nan
         for o in orders:
            w2 = sp.w[o]
            x2 = np.arange(w2.size)
            f2 = sp.f[o]
            e2 = sp.e[o]
            b2 = sp.bpmap[o] | msksky[o]

            if inst.name == 'FEROS':
               pmin = pomin[o]
               pmax = pomax[o]

            #if inst.name=='FEROS' and fib!='B':  # filter
               #hh = np.argsort(sp.f[o]); ii=hh[0:len(hh)*0.98]; pind=np.intersect1d(pind, ii)
            b2[:pmin] |= flag.out
            b2[pmax:] |= flag.out
            b2[tellmask(w2)>0.01] |= flag.atm
            b2[skymsk(w2)>0.01] |= flag.sky
            b2[(tellmask(barshift(w2, -spt.berv+sp.berv+(tplrv-targrv)))>0.01)!=0] |= flag.badT   # flag for bad template
            b2[(skymsk(barshift(w2, -spt.berv+sp.berv+(tplrv-targrv)))>0.01)!=0] |= flag.badT
            #if inst.name == 'HARPS':
               #b2[lstarmask(barshift(w2,sp.berv))>0.01] |= flag.lowQ
               #pause()
            pind = x2[b2==0]
            if not pind.size: continue

            wmod = barshift(w2, np.nan_to_num(sp.berv))   # berv can be NaN, e.g. calibration FP, ...
            if debug:   # check the input
               gplot(dopshift(ww[o],tplrv), ff[o], ',', dopshift(wmod,targrv), f2)
               pause(o)
            rchio = 1

            if ccf:
               '''METHOD CCF'''
               f2 *= b2==0
               par, f2mod, vCCF, stat = CCF(lam2wave(ccfmask[:,0]), ccfmask[:,1], wmod, f2, targrv+v_lo, targrv+v_hi, e_y2=e2, keep=pind, plot=(o in look)+2*(o in lookssr), ccfmode=ccfmode)

               rvccf[n,o] = par.params[0] * 1000
               e_rvccf[n,o] = par.perror[0] * 1000
               keep = pind
               rchio = 1
               #pause(rvccf[n,o], e_rvccf[n,o])
            elif diff_rv:
               '''METHOD DRIFT MEASUREMENT'''
               # Correlate the residuals with the first derivative.
               # spt.w[o] and wmod are ignored, i.e. correction for sp.berv
               # We estimate the gradients using finite differences. NaN will propagate to neighbours!?
               #stop()
               #toosharp = spt.bpmap[o] * 0
               # flag template with too sharp pixels and likely cosmics or overspill
               spt.bpmap[o][1:] |= (0.15*spt.f[o][1:] > np.abs(spt.f[o][:-1])) * flag.sat
               spt.bpmap[o][:-1] |= (0.15*spt.f[o][:-1] > np.abs(spt.f[o][1:])) * flag.sat

               b2[np.where(spt.bpmap[o][1:])[0]] |= flag.badT    # flag pixels with bad neighbours resulting in bad gradients
               b2[np.where(spt.bpmap[o][:-1])[0]+1] |= flag.badT
               # flag also the overnext, espcially if using derivative from spline ( ringing, CARM_VIS overspill)
               b2[np.where(spt.bpmap[o][2:] & flag.sat)[0]] |= flag.badT
               b2[np.where(spt.bpmap[o][:-2] & flag.sat)[0]+2] |= flag.badT
               pind = x2[b2==0]
               if 1:
                  # robust pre-filtering
                  rr = sp.f[o]/spt.f[o]   # spectrum ratio
                  qq = quantile(rr[pind], [0.25,0.5,0.75])
                  clip_lo = qq[1] - 5*(qq[1]-qq[0])
                  clip_hi = qq[1] + 5*(qq[2]-qq[1])
                  #gplot(pind, rr[pind], ', %s,%s,%s,%s,%s' % (tuple(qq)+(clip_lo, clip_hi)))
                  b2[(b2==0) & ((rr < clip_lo) | (rr > clip_hi))] |= flag.clip
                  # print 'pre-clipped:', np.sum((b2 & flag.clip) >0)
                  sp.f[o][pind]/spt.f[o][pind]
                  pind = x2[b2==0]
                  #pause('pre-clip', o)

               if 0:
                  # derivative from gradient (one side numerial derivative)
                  dy = np.gradient(spt.f[o], np.gradient(spt.w[o]))
                  ddy = np.gradient(dy, np.gradient(spt.w[o]))
               elif 1:
                  # symmetric numerical derivative
                  # may underestimate gradients
                  #dy = 0 * spt.f[o]
                  #dy[1:-1] = (spt.f[o][2:]-spt.f[o][:-2]) / (spt.w[o][2:]-spt.w[o][:-2])
                  if tuple(map(int,np.__version__.split("."))) < (1,13):
                     dy = np.gradient(spt.f[o], np.gradient(spt.w[o]))
                  else:
                     dy = np.gradient(spt.f[o], spt.w[o])
                  ddy = 0 * spt.f[o]
                  # ddy[1:-1] = (spt.f[o][2:]-2*spt.f[o][1:-1]+spt.f[o][:-2]) / ((spt.w[o][2:]-spt.w[o][:-2])**2 / 4) # assumes dw_i ~ dw_(i+1)
                  # ddy = ((f(x+h2)-f(x))/h2 - (f(x)-f(x-h1))/h1) / ((h1+h2)/2)
                  ddy[1:-1] = ((spt.f[o][2:]-spt.f[o][1:-1])/(spt.w[o][2:]-spt.w[o][1:-1]) - (spt.f[o][1:-1]-spt.f[o][:-2]) / (spt.w[o][1:-1]-spt.w[o][:-2])) / ((spt.w[o][2:]-spt.w[o][:-2])/2)
               else:
                  # first and second derivative from spline interpolation
                  # may overestimate gradients
                  # ringing may influence gradients
                  kkk = interpolate.splrep(spt.w[o], spt.f[o])
                  dy = interpolate.splev(spt.w[o], kkk, der=1)
                  ddy = interpolate.splev(spt.w[o], kkk, der=2)
               if 0:#o==29:
                  gplot(f2, 'w p,', spt.f[o], 'w lp,', pind, f2[pind])
                  pause(o)

               par, f2mod, keep, stat = fitspec((spt.w[o],spt.f[o]), wmod,f2,e2, v=targrv/1000, clip=kapsig, nclip=nclip,keep=pind, df=dy, plot=(o in lookssr)|(2*(not safemode)))

               e_vi = np.abs(e2/dy)*c*1000.   # velocity error per pixel
               e_vi_min = 1/ np.sqrt(np.sum(1/e_vi[keep]**2)) # total velocity error (Butler et al., 1996)
               #print(np.abs(e_vi[keep]).min(), e_vi_min)
               rchio = 1
               # compare both gradients
               #gplot(keep,e_vi[keep])
               #gplot(spt.f[o][keep], dy[keep]/((sp.f[o][2:]-sp.f[o][:-2]) / (sp.w[o][2:]-sp.w[o][:-2]))[keep-1])
               #gplot(keep, dy[keep]/((sp.f[o][2:]-sp.f[o][:-2]) / (sp.w[o][2:]-sp.w[o][:-2]))[keep-1]*50, ',', keep,sp.f[o][keep],spt.f[o][keep], ', "" us 1:3')
               #pause()

               #v = -c * np.dot(1/e2[keep]**2*dy[keep], (f2-f2mod)[keep]) / np.dot(1/e2[keep]**2*dy[keep], dy[keep]) / A**2
               #dlwo = c**2 *np.dot(1/e2[keep]**2*ddy[keep], (f2-f2mod)[keep]) / np.dot(1/e2[keep]**2*ddy[keep], ddy[keep]) / A**2
               #e_dlwo = c**2 * np.sqrt(1 / np.dot(1/e2[keep]**2*ddy[keep], ddy[keep]) / A**2)
               #rchi = rms(((f2-f2mod)-dlwo/c**2*ddy)[keep]/e2[keep])
               #e_dlwo *= 1000*rchi
               #print(par.params[1],par.params[0], v, dlwo*1000, e_dlwo)

            else:
               '''DEFAULT METHOD: LEAST SQUARE'''
               if 0:
                  gplot(ww[o], ff[o], 'w l,', wmod,f2/np.mean(f2)*np.mean(ff[o]),'w lp')
                  pause(o)

               # pause()
               if o==41: pind=pind[:-9]   # @CARM_NIR?
               par, f2mod, keep, stat, chi2mapo = fitspec(TPL[o], wmod, f2, e2, v=targrv-tplrv, clip=kapsig, nclip=nclip, keep=pind, indmod=np.s_[pmin:pmax], plot=(o in lookssr)|(2*(not safemode)), deg=deg, chi2map=True)

               if diff_width:
                  '''we need the model at the observation and oversampled since we need the second derivative including the polynomial'''
                  #f2mod = calcspec(wmod, *par.params) #calcspec does not work when w < wtmin
                  #ftmod_tmp = calcspec(ww[o], *par.params) #calcspec does not work when w < wtmin
                  #wo = TPL[o]
                  #i0 = np.searchsorted(dopshift(ww[o],par.params[0]), ww[o][0])
                  #i1 = np.searchsorted(dopshift(ww[o],par.params[0]), ww[o][-1]) - 1
                  #ftmod_tmp = 0*ww[o]
                  #ftmod_tmp[i0:i1] = calcspec(ww[o][i0:i1], *par.params) #calcspec does not work when w < wtmin
                  #ftmod_tmp[0:i0] = ftmod_tmp[i0]
                  #ftmod_tmp[i1:] = ftmod_tmp[i1-1]
                  poly = calcspec(wmod, *par.params, retpoly=True)


                  '''estimate differential changes in line width ("FWHM")
                     correlate the residuals with the second derivative
                  We want to scale the residuals (without the poly) to the second derivative 
                  f = p * t
                  t = f/p
                  a*ddt ~ y/p - t = (y-pt)/p = r/p
                  r = p*a*ddt
                  '''
                  dy = poly * TPL[o](dopshift(wmod, par.params[0]), der=1)
                  ddy = poly * TPL[o](dopshift(wmod, par.params[0]), der=2)

                  if not def_wlog:
                     dy *= wmod
                     ddy *= wmod**2
                  v = -c * np.dot(1/e2[keep]**2*dy[keep], (f2-f2mod)[keep]) / np.dot(1/e2[keep]**2*dy[keep], dy[keep])

                  moon = 0
                  if moon:
                     # simple moon contamination
                     a1, a0 = np.polyfit(f2mod[keep], f2[keep], 1)
                     print(a1, a0)
                     if 0:
                        gplot(w2[keep], f2mod[keep],'w lp,',w2[keep], f2[keep], ' w lp')
                        gplot(f2mod[keep], f2[keep], w2[keep],' palette , x, %s+%s*x' %(a0,a1) )
                        gplot(w2[keep], f2mod[keep],'w lp,',w2[keep], f2[keep], ' w lp,',w2[keep], (f2[keep]-a0)/a1, ' w lp')
                        # residuals
                        gplot(w2[keep], f2[keep]-f2mod[keep], ' w lp, ', w2[keep], f2[keep]-a1*f2mod[keep]-a0, 'w lp lc 3')
                     #pause(o)

                     f2 = f2c = (f2-a0)/a1 # correct observation

                     dlwo = c**2 * np.dot(1/e2[keep]**2*ddy[keep], (f2c-f2mod)[keep]) / np.dot(1/e2[keep]**2*ddy[keep], ddy[keep])
                     e_dlwo = c**2 * np.sqrt(1 / np.dot(1/e2[keep]**2, ddy[keep]**2))
                  else:
                     dlwo = c**2 * np.dot(1/e2[keep]**2*ddy[keep], (f2-f2mod)[keep]) / np.dot(1/e2[keep]**2*ddy[keep], ddy[keep])
                     e_dlwo = c**2 * np.sqrt(1 / np.dot(1/e2[keep]**2, ddy[keep]**2))
                  drchi = rms(((f2-f2mod) - dlwo/c**2*ddy)[keep] / e2[keep])
                  #print(par.params[1],par.params[0], v, dlwo*1000, e_dlwo)
                  if np.isnan(dlwo) and not safemode: pause()
                  dlw[n,o] = dlwo * 1000       # convert from (km/s) to m/s km/s
                  e_dlw[n,o] = e_dlwo * 1000 * drchi
                  if 0:
                     gplot(wmod,(f2-f2mod),e2, 'us 1:2:3 w e,', wmod[keep],(f2-f2mod)[keep], 'us 1:2')
                     ogplot(wmod,dlwo/c**2*ddy); #ogplot(wmod,f2, 'axis x1y2')
                     pause(o, 'dLW', dlw[n,o])

            fmod[o] = f2mod
            if par.perror is None: par.perror = [0.,0.,0.,0.]
            results[sp.timeid][o] = par
            rv[n,o] = rvo = par.params[0] * 1000. # to [km/s] - sp.drift
            snr[n,o] = stat['snr']
            rchi[n,o] = stat['std']
            Nok[n,o] = len(keep)

            if not diff_rv:
               vgrid = chi2mapo[0]
               chi2map[o] = chi2mapo[1] # chi2mapo[1].min() - (a[0]+a[1]*v+a[2]*v**2)
               # pause(o, chi2map[o].min())
               # for
               ##chi2map[o] = chi2mapo[1].min() + ((chi2mapo[0]-par.params[0])/ (par.perror[0]))**2
               ##chi2map[o] =  2*np.log(np.pi*par.perror[0]) - ((chi2mapo[0]-par.params[0])/ (par.perror[0]))**2
               #chi2map[o] = ((chi2mapo[0]-par.params[0])/ (par.perror[0]))**2
               #pause(o)

            e_rv[n,o] = par.perror[0] * stat['std'] * 1000
            if verb: print("%s-%02u  %s  %7.2f +/- %5.2f m/s %5.2f %5.1f it=%s %s" % (n+1, o, sp.timeid, rvo, par.perror[0]*1000., stat['std'], stat['snr'], par.niter, np.size(keep)))

            clipped = np.sort(list(set(pind).difference(set(keep))))
            if len(clipped):
               b2[clipped] = flag.clip
            if o in lookp or (o in look and iterate==niter) or (not safemode and (abs(rvo/1000-targrv+tplrv)>rvwarn and not sp.flag) or debug>1):
               if def_wlog: w2 = np.exp(w2)
               res = np.nan * f2
               res[pmin:pmax] = (f2[pmin:pmax]-f2mod[pmin:pmax]) / e2[pmin:pmax]  # normalised residuals
               b = str(stat['std'])
               gplot.key('left Left rev samplen 2 tit "%s (o=%s, v=%.2f+/-%.2f m/s)"'%(obj,o,rvo, e_rv[n,o]))
               gplot.ytics('nomirr; set y2tics; set y2range [-5*%f:35*%f]; set bar 0.5'%(rchio, rchio))
               gplot.put('i=1; bind "$" "i = i%2+1; xlab=i==1?\\"pixel\\":\\"wavelength\\"; set xlabel xlab; set xra [*:*]; print i; repl"')
               gplot-('[][][][-5:35]', x2, w2, f2, e2.clip(0.,f2.max()), 'us (column(i)):3:4 w errorli t "'+sp.timeid+' all"')
               gplot<(x2,w2, f2, ((b2==0)|(b2==flag.clip))*0.5, 1+4*(b2==flag.clip), 'us (column(i)):3:4:5 w p pt 7 lc var ps var t "'+sp.timeid+' telluric free"')
               gplot<(x2,w2, f2mod,(b2==0)*0.5, 'us (column(i)):3:4 w lp lt 3 pt 7 ps var t "Fmod"')
               gplot<(x2,w2, res, b2, "us (column(i)):3:4 w lp pt 7 ps 0.5 lc var axis x1y2 t 'residuals'")
               # legend with translation of bpmap, plot dummy using NaN function
               gplot<(", ".join(["NaN w p pt 7 ps 0.5 lc "+str(f) +" t '"+str(f)+" "+",".join(flag.translate(f))+"'" for f in np.unique(b2)]))

               # zero line and 1 sigma level
               gplot<("0 axis x1y2 lt 3 t'',"+b+" axis x1y2 lt 1,-"+b+" axis x1y2 lt 1 t ''")

               if atmspec:
                  gplot<(x2, w2, atmmod(np.log(w2))*40-5, 'us (column(i)):3 w l lt rgb 0x999999  axis x1y2')
                  gplot<(x2,w2, atmmod(np.log(w2)) * ((b2&flag.atm)==flag.atm)*40-5, 'us (column(i)):3 w filledcurve x1 fs transparent solid 0.5 noborder lc 9 axis x1y2 t "tellurics"')
               else:
                  gplot<(x2,w2, ((b2&flag.atm)!=flag.atm)*40-5, 'us (column(i)):3 w filledcurve x2 fs transparent solid 0.5 noborder lc 9 axis x1y2 t "tellurics"')
               gplot+(x2,w2, ((b2&flag.sky)!=flag.sky)*40-5, 'us (column(i)):3 w filledcurve x2 fs transparent solid 0.5 noborder lc 6 axis x1y2 t "sky"')
               pause('large RV ' if abs(rvo/1000-targrv+tplrv)>rvwarn else 'look ', o, ' rv = %.3f +/- %.3f m/s   rchi = %.2f' %(rvo, e_rv[n,o], rchi[n,o]))

         # end loop over orders

         # ind = setdiff1d(where(e_rv[n]>0.)[0],[71]) # do not use the failed and last order
         ind, = where(np.isfinite(e_rv[n])) # do not use the failed and last order
         rvm[n], rvmerr[n] = np.median(rv[n,ind]), std(rv[n,ind])
         if len(ind) > 1: rvmerr[n] /= (len(ind)-1)**0.5

         # Mean RV
         RV[n], e_RV[n] = wsem(rv[n,ind], e=e_rv[n,ind])
         RVc[n] = RV[n] - np.nan_to_num(sp.drift) - np.nan_to_num(sp.sa)
         e_RVc[n] = np.sqrt(e_RV[n]**2 + np.nan_to_num(sp.e_drift)**2)
         print(n+1, '/', nspec, sp.timeid, sp.bjd, RV[n], e_RV[n])

         # Chromatic trend
         if len(ind)>1:
            # scipy version
            #np.polynomial.polyval(x,[a,b])
            def func(x, a, b): return a + b*x #  np.polynomial.polyval(x,a)
            #
            # x = np.mean(np.exp(spt.w) if def_wlog else spt.w, axis=1)    # lambda
            # x = 1/np.mean(np.exp(spt.w) if def_wlog else spt.w, axis=1)  # 1/lambda
            x = np.nanmean(spt.w if def_wlog else np.log(spt.w), axis=1)  # ln(lambda)
            xc = np.nanmean(x[ind])   # only to center the trend fit
            # fit trend with curve_fit to get parameter error
            pval, cov = curve_fit(func, x[ind]-xc, rv[n][ind], [0.0, 0.0], e_rv[n][ind])
            perr = np.sqrt(np.diag(cov))
            #print(cov, pval)
            #pval, cov = np.polyfit(x[ind]-xc, rv[n][ind], 1, w=1/e_rv[n][ind], cov=True)
            #print(cov, pval)

            l_v = np.exp(-(pval[0]-RV[n])/pval[1]+xc)
            CRX[n], e_CRX[n], xo[n] = pval[1], perr[1], x
            tCRX[n] = CRX[n], e_CRX[n], pval[0], perr[0], l_v

            # manual regression
            lhs = np.array([1/e_rv[n][ind], (x[ind]-xc)/e_rv[n][ind]]).T
            cov = np.linalg.inv(np.dot(lhs.T, lhs))
            #e_CRX[n] = np.sqrt(cov[1,1])

            # with scaling
            scale = np.sqrt((lhs*lhs).sum(axis=0))
            cov = np.linalg.inv(np.dot((lhs/scale).T, lhs/scale)) / np.outer(scale, scale)

            #pause(e_CRX[n])

            #coli,stat = polynomial.polyfit(arange(len(rv[n]))[ind],rv[n][ind], 1, w=1./e_rv[n][ind], full=True)
            if 0:   # show trend in each order
               gplot.log('x; set autoscale xfix; set xtic add (0'+(",%i"*10)%tuple((np.arange(10)+1)*1000)+')')
               gplot(np.exp(x[ind]), rv[n][ind], e_rv[n][ind], ' us 1:2:3 w e pt 7, %f+%f*log(x/%f), %f' % (RV[n], pval[1],l_v,RV[n]))
               pause()

         if len(ind)>1 and not diff_rv:
            # ML version of chromatic trend
            oo = ~np.isnan(chi2map[:,0]) & ~np.isnan(rchi[n])

            gg = Chi2Map(chi2map, (v_lo, v_step), RV[n]/1000, e_RV[n]/1000, rv[n,oo]/1000, e_rv[n,oo]/1000, orders=oo, keytitle=obj+' ('+inst.name+')\\n'+sp.timeid, rchi=rchi[n], No=Nok[n], name='')
            mlRV[n], e_mlRV[n] = gg.mlRV, gg.e_mlRV

            mlRVc[n] = mlRV[n] - np.nan_to_num(sp.drift) - np.nan_to_num(sp.sa)
            e_mlRVc[n] = np.sqrt(e_mlRV[n]**2 + np.nan_to_num(sp.e_drift)**2)

            if lookmlRV:
               gg.plot()
               pause(n, mlRV[n], e_mlRV[n])

            try:
               mlCRX[n], e_mlCRX[n] = gg.mlcrx(x, xc, oo)
            except:
               print("Could not compute mlCRX for n=", n)

            if lookmlCRX:
               gg.plot_fit()
               pause(n, CRX[n], mlCRX[n])


         # Line Indices
         vabs = tplrv + RV[n]/1000.
         kwargs = {'line_o':line_o, 'plot':looki}
         if meas_index:
            halpha += [getHalpha(vabs, 'Halpha', **kwargs)]
            haleft += [getHalpha(vabs, 'Haleft', **kwargs)]
            harigh += [getHalpha(vabs, 'Harigh', **kwargs)]
            cai += [getHalpha(vabs, 'CaI', **kwargs)]
            cak += [getHalpha(vabs, 'CaK', **kwargs)] if inst.name=='HARPS' else [(np.nan,np.nan)]
            cah += [getHalpha(vabs, 'CaH', **kwargs)] if inst.name=='HARPS' else [(np.nan,np.nan)]
         if meas_CaIRT:
            irt1 +=  [getHalpha(vabs, 'CaIRT1', **kwargs)]
            irt1a += [getHalpha(vabs, 'CaIRT1a', **kwargs)]
            irt1b += [getHalpha(vabs, 'CaIRT1b', **kwargs)]
            irt2 +=  [getHalpha(vabs, 'CaIRT2', **kwargs)]
            irt2a += [getHalpha(vabs, 'CaIRT2a', **kwargs)]
            irt2b += [getHalpha(vabs, 'CaIRT2b', **kwargs)]
            irt3 +=  [getHalpha(vabs, 'CaIRT3', **kwargs)]
            irt3a += [getHalpha(vabs, 'CaIRT3a', **kwargs)]
            irt3b += [getHalpha(vabs, 'CaIRT3b', **kwargs)]
         if meas_NaD:
            nad1 += [getHalpha(vabs, 'NaD1', **kwargs)]
            nad2 += [getHalpha(vabs, 'NaD2', **kwargs)]
            nadr1 += [getHalpha(vabs, 'NaDref1', **kwargs)]
            nadr2 += [getHalpha(vabs, 'NaDref2', **kwargs)]
            nadr3 += [getHalpha(vabs, 'NaDref3', **kwargs)]

         if diff_width:
            ind, = where(np.isfinite(e_dlw[n]))
            dLW[n], e_dLW[n] = wsem(dlw[n,ind], e=e_dlw[n,ind])

         if 0: # plot RVs of all orders
            gplot.key('title "rv %i:  %s"' %(n+1,sp.timeid))
            gplot(orders,rv[n,orders],e_rv[n,orders],rvccf[n,orders],e_rvccf[n,orders],'us 1:2:3 w e, "" us 1:4:5 w e t "ccf", {0}, {0} - {1}, {0}+{1} lt 2'.format(RV[n],e_RV[n]))
            pause('rvo')
         if 0: # plot dLW of all orders
            gplot.key('title "dLW %i:  %s"' %(n+1,sp.timeid))
            gplot(orders, dlw[n][orders], e_dlw[n][orders], ' w e, %f, %f, %f lt 2'%(dLW[n]+e_dLW[n],dLW[n],dLW[n]-e_dLW[n]))
            pause('dlw')

         if outfmt and not np.isnan(RV[n]):   # write residuals
            data = {'fmod': fmod, 'wave': sp.w, 'spec': sp.f,
                    'err': sp.e, 'bpmap': sp.bpmap, 'waverest': redshift(sp.w, vo=sp.berv, ve=RV[n]/1000.)}
            outfile = os.path.basename(sp.filename)
            outfile = os.path.splitext(outfile)[0] + outsuf
            if 'res' in outfmt: data['res'] = sp.f - fmod
            if 'ratio' in outfmt: data['ratio'] = sp.f / fmod

            sph = Spectrum(sp.filename, inst=inst, pfits=True, drs=drs, fib=fib, targ=targ).header
            sph['HIERARCH SERVAL RV'] = (RV[n], '[m/s] Radial velocity')
            sph['HIERARCH SERVAL E_RV'] = (e_RV[n], '[m/s] RV error estimate')
            sph['HIERARCH SERVAL RVC'] = (RVc[n], '[m/s] RV drift corrected')
            sph['HIERARCH SERVAL E_RVC'] = (e_RVc[n], '[m/s] RVC error estimate')
            write_res(outdir+'res/'+outfile, data, outfmt, sph, clobber=1)

         if outchi and not np.isnan(RV[n]):   # write residuals
            gplot.palette('defined (0 "blue", 1 "green", 2 "red")')
            gplot.xlabel('"v [m/s]"; set ylabel "chi^2/max(chi^2)"; set cblabel "order"')
            #gplot(chi2map, ' matrix us ($1*%s+%s):3:2 w l palette'%(v_step, v_lo))
            if 0:
               gplot(chi2map.T /chi2map.max(axis=1), ' matrix us ($1*%s+%s):3:2 w l palette'%(v_step, v_lo))
            outfile = os.path.basename(sp.filename)
            outfile = os.path.splitext(outfile)[0] + '_chi2map.fits'
            hdr = spt.header[0:10]
            hdr.insert('COMMENT', ('CDELT1', v_step))
            hdr['CTYPE1'] = 'linear'
            hdr['CUNIT1'] = 'km/s'
            hdr['CRVAL1'] = v_lo
            hdr['CRPIX1'] = 1
            hdr['CDELT1'] = v_step
            hdr['CTYPE2'] = 'linear'
            hdr['CRVAL2'] = 1
            hdr['CRPIX2'] = 1
            hdr['CDELT2'] = 1
            write_fits(outdir+'res/'+outfile, chi2map, hdr+spt.header[10:])

         if n>0 and not safemode:
            # plot time series
            gplot(bjd-2450000, RV, e_RV, ' us 1:2:3 w e pt 7') # explicitly specify columns to deal with NaNs
         #pause()

      # end of loop over spectra

      if iterate < niter:
         np.savetxt(prefile, list(zip(bjd[:nspecok], RV[:nspecok], e_RV[:nspecok])), fmt="%s")
         continue
      # write final results
      rvfile = outdir+obj+fibsuf+'.dat'
      rvcfile = outdir+obj+'.rvc'+fibsuf+'.dat'
      crxfile = outdir+obj+'.crx'+fibsuf+'.dat'
      mlcfile = outdir+obj+'.mlc'+fibsuf+'.dat' # maximum likehood estimated RVCs and CRX
      srvfile = outdir+obj+'.srv'+fibsuf+'.dat' # serval top-level file
      vsinifile = outdir+obj+'.vsini'+fibsuf+'.dat'
      rvunit = [open(rvfile, 'w'), open(outdir+obj+'.badrv'+fibsuf+'.dat', 'w')]
      rvounit = [open(rvofile, 'w'), open(rvofile+'bad', 'w')]
      rvcunit = [open(rvcfile, 'w'), open(rvcfile+'bad', 'w')]
      crxunit = [open(crxfile, 'w'), open(crxfile+'bad', 'w')]
      mlcunit = [open(mlcfile, 'w'), open(mlcfile+'bad', 'w')]
      srvunit = [open(srvfile, 'w'), open(srvfile+'bad', 'w')]
      mypfile = [open(e_rvofile, 'w'), open(e_rvofile+'bad', 'w')]
      snrunit = [open(snrfile, 'w'), open(snrfile+'bad', 'w')]
      chiunit = [open(chifile, 'w'), open(chifile+'bad', 'w')]
      dlwunit = [open(dlwfile, 'w'), open(dlwfile+'bad', 'w')]
      e_dlwunit = [open(e_dlwfile, 'w'), open(e_dlwfile+'bad', 'w')]
      halunit = irtunit = nadunit = []
      if meas_index:
         halunit = [open(halfile, 'w'), open(halfile+'bad', 'w')]
      if meas_CaIRT:
         irtunit = [open(irtfile, 'w'), open(irtfile+'bad', 'w')]
      if meas_NaD:
         nadunit = [open(nadfile, 'w'), open(nadfile+'bad', 'w')]

      if vsiniauto:
         with open(vsinifile, 'w') as vsiniunit:
            print('#Median vsini [km/s]:', file=vsiniunit)
            print(np.nanmedian(VSINI[:,0]),np.sqrt(2/np.pi)*np.nanstd(VSINI[:,0]), file=vsiniunit)
            print('\n#order', 'vsini[km/s]','error[km/s]', file=vsiniunit)
            for o in range(nord):
               print(o, VSINI[o,0], VSINI[o,1], file=vsiniunit)

      for n,sp in enumerate(spoklist):
         if np.isnan(rvm[n]): sp.flag |= sflag.rvnan
         rvflag = int((sp.flag&(sflag.config+sflag.iod+sflag.rvnan)) > 0)
         if rvflag: 'nan RV for file: '+sp.filename
         print(sp.bjd, RVc[n], e_RVc[n], file=rvunit[int(rvflag or np.isnan(sp.drift))])
         print(sp.bjd, RV[n], e_RV[n], rvm[n], rvmerr[n], *rv[n], file=rvounit[rvflag])
         print(sp.bjd, RV[n], e_RV[n], rvm[n], rvmerr[n], *e_rv[n], file=mypfile[rvflag])
         print(sp.bjd, RVc[n], e_RVc[n], sp.drift, sp.e_drift, RV[n], e_RV[n], sp.berv, sp.sa, file=rvcunit[rvflag])
         print(sp.bjd, *(list(tCRX[n])+list(xo[n])), file=crxunit[rvflag])
         print(sp.bjd, RVc[n], e_RVc[n], CRX[n], e_CRX[n], dLW[n], e_dLW[n], file=srvunit[rvflag])
         print(sp.bjd, mlRVc[n], e_mlRVc[n], mlCRX[n], e_mlCRX[n], dLW[n], e_dLW[n], file=mlcunit[rvflag])
         print(sp.bjd, dLW[n], e_dLW[n], *dlw[n], file=dlwunit[rvflag])
         print(sp.bjd, dLW[n], e_dLW[n], *e_dlw[n], file=e_dlwunit[rvflag])
         print(sp.bjd, np.nansum(snr[n]**2)**0.5, *snr[n], file=snrunit[rvflag])
         print(sp.bjd, *rchi[n], file=chiunit[rvflag])
         if meas_index:
            print(sp.bjd, *(lineindex(halpha[n],harigh[n],haleft[n]) + halpha[n] + haleft[n] + harigh[n] + lineindex(cai[n],harigh[n],haleft[n])), file=halunit[rvflag])  #,cah[n][0],cah[n][1]
         if meas_CaIRT:
            print(sp.bjd, *(lineindex(irt1[n], irt1a[n], irt1b[n]) + lineindex(irt2[n], irt2a[n], irt2b[n]) + lineindex(irt3[n], irt3a[n], irt3b[n])), file=irtunit[rvflag])
         if meas_NaD:
            print(sp.bjd, *(lineindex(nad1[n],nadr1[n],nadr2[n]) + lineindex(nad2[n],nadr2[n],nadr3[n])), file=nadunit[rvflag])

      for ifile in rvunit + rvounit + rvcunit + snrunit + chiunit + mypfile + crxunit \
                 + srvunit + mlcunit + dlwunit +e_dlwunit + halunit + irtunit + nadunit:
         ifile.close()

      t2 = time.time() - t0
      print()
      print(nspec, 'spectra processed', rvfile+"  (total %s, compu %s)\n" %(minsec(t2), minsec(t2-t1)))

   # end of iterate loop

   if not driftref and nspec>1:
      x = analyse_rv(obj, postiter=postiter, fibsuf=fibsuf, safemode=safemode)
      if safemode<2: pause('TheEnd')

   bp or sys.stdout.logfile.close()


def arg2slice(arg):
   """Convert string argument to a slice."""
   # We want four cases for indexing: None, int, list of ints, slices.
   # Use [] as default, so 'in' can be used.
   if isinstance(arg, str):
      arg = eval('np.s_['+arg+']')
   return [arg] if isinstance(arg, int) else arg


if __name__ == "__main__":
   insts = [os.path.basename(i)[5:-3] for i in glob.glob(servalsrc+'inst_*.py')]

   # check first the instrument with preparsing
   preparser = argparse.ArgumentParser(add_help=False)
   preparser.add_argument('args', nargs='*')
   preparser.add_argument('-inst', help='instrument', default='HARPS', choices=insts)
   preargs, _ =  preparser.parse_known_args()

   inst = preargs.inst
   inst = importlib.import_module('inst_'+inst)

   # instrument specific default (see inst_*.py)
   pmin = getattr(inst, 'pmin', 300)
   pmax = getattr(inst, 'pmax', 3800)
   oset = getattr(inst, 'oset', {'FEROS':'10:'}.get(inst.name,':'))
   coset = getattr(inst, 'coset', None)
   ofac = getattr(inst, 'ofac', 1.)
   snmax = getattr(inst, 'snmax', 400.)

   default = " (default: %(default)s)."
   epilog = """\

   Default parameters are for """+inst.name+""".

   usage example:
   %(prog)s tag dir_or_filelist -targ gj699 -snmin 10 -oset 40:
   """
   parser = argparse.ArgumentParser(description=description, epilog=epilog, add_help=False, formatter_class=argparse.RawTextHelpFormatter)
   argopt = parser.add_argument   # function short cut
   argopt('obj', help='Tag, output directory and file prefix (e.g. Object name).')
   argopt('dir_or_inputlist', help='Directory name with reduced data fits/tar or a file listing the spectra (only suffixes .txt or .lis accepted).', nargs='?')
   argopt('-targ', help='Target name requested in simbad for coordinates, proper motion, parallax and absolute RV.')
   argopt('-targrade', help='Target coordinates: [ra|hh:mm:ss.sss de|de:mm:ss.sss].', nargs=2, default=[None,None])
   argopt('-targpm', help='Target proper motion: pmra [mas/yr] pmde [mas/yr].', nargs=2, type=float, default=[0.0,0.0])
   argopt('-targplx', help='Target parallax', type=float, default='nan')
   argopt('-targrv', help='[km/s] Target RV guess (for index measures) [float, "drsspt", "drsmed", "targ", None, "auto"]. None => no measure; targ => from simbad, hdr; auto => first from headers, second from simbad))', default={'CARM_NIR':None, 'else':'auto'})
   argopt('-atmmask', help='Telluric line mask ('' for no masking)'+default, default='auto', dest='atmfile')
   argopt('-atmwgt', help='Downweighting factor for coadding in telluric regions'+default, type=float, default=None)
   argopt('-atmspec', help='Telluric spectrum  (in fits format, e.g. lib/stdatmos_vis30a090rh0780p000t.fits) to correct spectra by simple division.'+default, type=str, default=None)
   argopt('-brvref', help='Barycentric RV code reference', choices=brvrefs, type=str, default='WE')
   argopt('-msklist', help='Ascii table with vacuum wavelengths to mask.', default='') # [flux and width]
   argopt('-mskwd', help='[km/s] Broadening width for msklist.', type=float, default=4.)
   argopt('-mskspec', help='Ascii 0-1 spectrum.'+default, default='')
   argopt('-cache',  help='store the spectra in memory instead rereading (intended for few spectra in ascii format)', action='store_true')
   argopt('-ccf',  help='mode ccf [with files]', nargs='?', const='th_mask_1kms.dat', type=str)
   argopt('-ccfmode', help='type for ccf template', nargs='?', default='box',
                      choices=['box', 'binless', 'gauss', 'trapeze'])
   argopt('-coset', help='index for order in coadding (default: oset)'+default, default=coset, type=arg2slice)
   argopt('-co_excl', help='orders to exclude in coadding (default: o_excl)', type=arg2slice)
   argopt('-ckappa', help='kappa sigma (or lower and upper) clip value in coadding. Zero values for no clipping'+default, nargs='+', type=float, default=(4.,4.))
   argopt('-deg',  help='degree for background polynomial', type=int, default=3)
   argopt('-distmax', help='[arcsec] Max distance telescope position from target coordinates.', nargs='?', type=float, const=30.)
   argopt('-driftref', help='reference file for drift mode', type=str)
   argopt('-fib',  help='fibre', choices=['', 'A', 'B', 'AB'], default='')
   argopt('-inst', help='instrument '+default, default='HARPS', choices=insts)
   argopt('-nset', '-iset', help='slice for file subset (e.g. 1:10, ::5)', default=':', type=arg2slice)
   argopt('-n_excl',  help='pattern for files to exclude', nargs='+', default=[])
   argopt('-kapsig', help='kappa sigma clip value'+default, type=float, default=3.0)
   argopt('-last', help='use last template (-tpl <obj>/template.fits)', action='store_true')
   argopt('-look', help='slice of orders to view the fit [:]', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-looki', help='list of indices to watch', nargs='*', choices=sorted(lines.keys()), default=[]) #, const=['Halpha'])
   argopt('-lookt', help='slice of orders to view the coadd fit [:]', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-lookp', help='slice of orders to view the preRV fit [:]', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-lookssr', help='slice of orders to view the ssr function [:]', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-lookvsini', help='slice of orders to view vsini fit [:]', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-lookmlRV', help='chi2map and master', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-lookmlCRX', help='chi2map and CRX fit ', nargs='?', default=[], const=':', type=arg2slice)
   argopt('-moonsep', help='Flag threshold for min. moon separation'+default, type=float, default=7.)
   argopt('-nclip', help='max. number of clipping iterations'+default, type=int, default=2)
   argopt('-niter', help='number of RV iterations'+default, type=int, default=2)
   argopt('-oset', help='index for order subset (e.g. 1:10, ::5)'+default, default=oset, type=arg2slice)
   argopt('-o_excl', help='Orders to exclude (e.g. 1,10,3)'+default, default=[], type=arg2slice)
   #argopt('-outmod', help='output the modelling results for each spectrum into a fits file',  choices=['ratio', 'HARPN', 'CARM_VIS', 'CARM_NIR', 'FEROS', 'FTS'])
   argopt('-ofac', help='oversampling factor in coadding'+default, default=ofac, type=float)
   argopt('-ofacauto', help='automatic knot spacing with BIC.', action='store_true')
   argopt('-outchi', help='output of the chi2 map', nargs='?', const='_chi2map.fits')
   argopt('-outfmt', help='output format of the fits file (default: None; const: fmod err res wave)', nargs='*', choices=['wave', 'waverest', 'err', 'fmod', 'res', 'spec', 'bpmap', 'ratio'], default=None)
   argopt('-outsuf', help='output suffix', default='_mod.fits')
   argopt('-pmin', help='Minimum pixel'+default, default=pmin, type=int)
   argopt('-pmax', help='Maximum pixel'+default, default=pmax, type=int)
   argopt('-pspline', help='pspline as coadd filter [smooth value]', nargs='?', const=0.0000001, dest='pspllam', type=float)
   argopt('-pmu', help='analog to GP mean. Default no GP penalty. Without the mean in each order. Otherwise this value.', nargs='?', const=True, type=float)
   argopt('-pe_mu', help='analog to GP mean deviation', default=5., type=float)
   argopt('-reana', help='flag reanalyse only', action='store_true')
   argopt('-rvwarn', help='[km/s] warning threshold in debug'+default, default=2., type=float)
   argopt('-safemode', help='does not pause or stop, optional level 1  2 (reana)', nargs='?', type=int, const=1, default=False)
   argopt('-skippre', help='Skip pre-RVs.', action='store_true')
   argopt('-skymsk', help='Sky emission line mask ('' for no masking)'+default, default='auto', dest='skyfile')
   argopt('-snmin', help='minimum S/N (considered as not bad and used in template building)'+default, default=10, type=float)
   argopt('-snmax', help='maximum S/N (considered as not bad and used in template building)'+default, default=snmax, type=float)
   argopt('-sunalt', help='Flag threshold for max. altitude of Sun above the horizon'+default, type=float, default=-12.)
   argopt('-tfmt', help='output format of the template. NMAP is a an estimate for the number of good data points for each knot. DDSPEC is the second derivative for cubic spline reconstruction. (default: SPEC SIG WAVE)', nargs='*', choices=['SPEC', 'SIG', 'WAVE', 'NMAP', 'DDSPEC'], default=['SPEC', 'SIG', 'WAVE'])
   argopt('-tpl',  help="template filename or directory, if None or integer a template is created by coadding, where highest S/N spectrum or the filenr is used as start tpl for the pre-RVs", nargs='?')
   argopt('-tplrv', help='[km/s] template RV. By default taken from the template header and set to 0 km/s for phoe tpl. [float, "tpl", "drsspt", "drsmed", "targ", None, "auto"]', default='auto')
   argopt('-tplvsini', help='[km/s] Rotational velocity to broaden template.', type=float)
   argopt('-tset',  help="slice for file subset in template creation", default=':', type=arg2slice)
   argopt('-verb', help='verbose', action='store_true')
   v_lo, v_hi, v_step = -5.5, 5.6, 0.1
   argopt('-vsiniauto', help='automatic vsini computation. Requires external serval template (-tpl)', action='store_true')
   argopt('-vrange', help='[km/s] velocity grid around targrv (v_lo, v_hi, v_step)'+default, nargs='*', default=(v_lo, v_hi, v_step), type=float)
   argopt('-vtfix', help='fix RV in template creation', action='store_true')
   argopt('-wfix', help='fix wavelength solution', action='store_true')
   argopt('-debug', help='debug flag', nargs='?', default=0, const=1)
   argopt('-bp',   help='break points, e.g.: read_spec.py:125', nargs='*', type=str)
   argopt('-pdb',  help='debug post_mortem', action='store_true')
   argopt('-cprofile', help='profiling', action='store_true')
   # use add_help=false and re-add with more arguments
   argopt('-?', '-h', '-help', '--help',  help='show this help message and exit', action='help')
   #parser.__dict__['_option_string_actions']['-h'].__dict__['option_strings'] += ['-?', '-help']


   for i, arg in enumerate(sys.argv):   # allow to parse negative floats
      if len(arg) and arg[0]=='-' and arg[1].isdigit(): sys.argv[i] = ' ' + arg
   print(sys.argv)

   args = parser.parse_args()
   globals().update(vars(args))

   inst = importlib.import_module('inst_'+inst)
   Spectrum.brvref = brvref

   if tpl and tpl.isdigit(): tpl = int(tpl)
   oset = arg2slice(oset)
   if isinstance(o_excl, dict): o_excl = arg2slice(o_excl[inst.name]) if inst.name in o_excl else []
   if isinstance(targrv, dict): targrv = targrv[inst.name] if inst.name in targrv else targrv['else']
   if coset is None: coset = oset
   if co_excl is None: co_excl = o_excl

   if skippre or vtfix or last or isinstance(tpl, str) or driftref:
      niter -= 1

   if dir_or_inputlist is None:
      ## execute last command
      #with open(obj+'/lastcmd.txt') as f:
         #lastcmd = f.read()
      #os.system(lastcmd)
      x = analyse_rv(obj, postiter=postiter)
      exit()

   #if targ is None: targ = obj   # sometimes convenient, but not always (targsa request is nasty)
   if len(vrange) == 1:
      v_lo, v_hi = -vrange[0], vrange[0]
   elif len(vrange) == 2:
      v_lo, v_hi = vrange
   elif len(vrange) == 3:
      v_lo, v_hi, v_step = vrange
   elif len(vrange) > 3:
      pause('too many args for -vrange')

   if len(ckappa) == 1:
      ckappa = ckappa * 2 # list with two entries ;)

   if outfmt == []:
      outfmt = ['fmod', 'err', 'res', 'wave']

   if cprofile:
      sys.argv.remove('-cprofile')
      os.system('python -m cProfile -s time -o speed.txt $SERVAL/src/serval.py '+" ".join(sys.argv[1:]))
      os.system('~/python/zechmeister/gprof2dot.py -f pstats speed.txt|  dot -Tsvg -o callgraph.svg')
      print("speed.txt")
      exit()

   if bp:
      # bp seems to slow down astropy imports (moon)
      try:
         from StringIO import StringIO ## for Python 2
      except ImportError:
         from io import StringIO ## for Python 3
      import pdb

      # initial continue must be given via stdin, emulate this here
      bp = ('break %s\n' * len(bp)) % tuple(bp)
      sys.stdin = StringIO(bp+'cont')
      pdb.set_trace()
      sys.stdin = sys.__stdin__
      print('logging turned off')
      # Reset the Logger function, forking the output + set_trace yields "^[[A" in command history

   # sys.modules['keyring'] = sys.modules['os'] # astroquery> keyrings slows down

   if vsiniauto:
      if niter < 3:
         niter = 3

   try:
      serval()
   except:
      if args.pdb:
         print('ex')
         import pdb
         e, m, tb = sys.exc_info()
         sys.stdout = sys.__stdout__  # resets the Logger
         pdb.post_mortem(tb)
      else:
         raise


