from __future__ import division

__author__ = 'Mathias Zechmeister'
__version__ = '2018-05-04'

import numpy as np
from scipy.linalg import solve_banded, solveh_banded
import os.path
try:
   # for debugging
   from gplot import *
   from pause import pause
   from ds9 import ds9
except:
   pass

from ctypes import c_int, c_long, c_double

ptr_double = np.ctypeslib.ndpointer(dtype=np.float)
ptr_int = np.ctypeslib.ndpointer(dtype=np.int)
 
_cbspline = np.ctypeslib.load_library(os.path.join(os.path.dirname(__file__), 'cbspline'), '.')

_cbspline.cbspl_Bk.argtypes = [ptr_double,    # x
                                 ptr_double,  # G
                                 ptr_int,     # kk
                                 c_int,       # n
                                 c_long]      # fix
_cbspline.lhsbnd_fill.argtypes = [ptr_double, # BTBbnd
                                 ptr_double,  # G
                                 ptr_double,  # w
                                 ptr_int,     # kk
                                 c_int,       # n
                                 c_int]       # nk
_cbspline.add_DTD.argtypes = [ptr_double,     # BTBbnd
                                 c_int,       # nk
                                 c_double,    # lam
                                 c_int]       # pord
_cbspline.rhs_fill.argtypes = [ptr_double,    # BTy
                                 ptr_double,  # G
                                 ptr_double,  # wy
                                 ptr_int,     # kk
                                 c_int]       # n
_cbspline.cholsol.argtypes = [ptr_double,     # A
                                 ptr_double,  # b
                                 c_int]       # n
_cbspline.cholbnd.argtypes = [ptr_double,     # Abnd
                                 ptr_double,  # b
                                 c_int,       # n
                                 c_int]       # bw

_cbspline.cholbnd_upper.argtypes = _cbspline.cholbnd.argtypes 
_cbspline.bandsol.argtypes = _cbspline.cholbnd.argtypes

def cbspline_Bk(x, K, xmin=None, xmax=None, fix=True):
   '''
   Uniform cubic B-spline with direct and vectorized computation and compressed storage.

   (K knots, K-1 intervals each covered by D+1 Bsplines, K+2D construction knots for Bsplines
   K+D-1 regression bsplines

   Parameters
   ----------
   x : float or array_like
       Data point in units of knots.
   K : integer
       Number of uniform knots.
   xmin : float
       Position of the first knot.
   xmax : float
       Position of the last knot.
   fix : boolean
       Fixes the end to avoid dummy, recommended for regression.

   Returns
   -------
   knot, coeff
       The four non-zero B-spline coefficients and the knot indices of the first coefficient (compressed format).

   Examples
   --------
   >>> x = np.arange(0, 9.01, 0.1)
   >>> B, k = cbspline_Bk(x, 10)
   >>> gplot(x, B, ', "" us 1:3, "" us 1:4, "" us 1:5')

   Depending on round off the edge point is attributed to the last or previous last knot.

   >>> cbspline_Bk([0,22], 100)
   (array([[0.16666667, 0.        ],
          [0.66666667, 0.16666667],
          [0.16666667, 0.66666667],
          [0.        , 0.16666667]]), array([ 0, 98]))
   >>> cbspline_Bk([22], 100, xmin=0, fix=False)
   (array([[0.16666667],
          [0.66666667],
          [0.16666667],
          [0.        ]]), array([99]))
   >>> cbspline_Bk([23], 100, xmin=0, fix=False)
   (array([[4.78309876e-43],
          [1.66666667e-01],
          [6.66666667e-01],
          [1.66666667e-01]]), array([98]))

   '''
   x = np.asarray(x)
   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   x = (K-1)/(xmax-xmin) * (x-xmin)

   kk, pp = divmod(x, 1)
   G = np.empty((4, x.size), order='F')

   if fix:
      idx, = np.where(kk==K-1)
      if idx.size:
         kk[idx] -= 1
         pp[idx] = 1   # it should be zero before

   if 0:
      ### direct computation with p
      G[0] = (1-pp)**3/6     # = ((1-pp)-pp*(1-pp))/2
      G[1] = (3*pp**3 - 6*pp**2 + 4)/6
      G[2] = (-3*pp**3 + 3*pp**2 + 3*pp + 1)/6
      G[3] = pp**3/6         # = ((pp+1)*pp-pp)/2
   elif 1:
      ### direct and symmetric computation with p and q
      qq = 1 - pp
      G[0] = qq**3/6
      G[3] = pp**3/6
      G[1] = 3*G[3] - pp**2 + 4/6   # -3*pp**2 (1-pp)  -3*pp**2 + 4
      G[2] = 3*G[0] - qq**2 + 4/6   # -3*pp**3 + 3*pp**2 + 3*pp + 1
   else: #faster? untested
      ### direct and symmetric computation with p and q
      G[0] = 1 - pp
      G[2] = -G[0]*G[0];  G[1] = -pp*pp    # -qq^2; -pp^2
      G[0] *= G[2];       G[3] = G[1]*pp   # -qq^3; -pp^3
      G[1:3] += 4/6                      # -pp^2+4/6
      G[2] -= 0.5*G[0]                     # -(-qq^3/2)
      G[1] -= 0.5*G[3]

      G[0] *= -1/6
      G[3] *= -1/6

   return G, kk.astype(int)

def _cbspline_Bk(x, K, xmin=None, xmax=None, fix=True):
   '''C implementation of cbspline_Bk.'''
   x = np.asarray(x)
   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   x = (K-1)/(xmax-xmin) * (x-xmin)

   G = np.empty((4, x.size), order='F')   # "Fortran" because in C the outer loop is over x
   kk = np.empty(x.size, dtype=int)
   _cbspline.cbspl_Bk(x, G, kk, x.size, fix*(K-2))
   return G, kk

def bk2bknat(G, kk, K):
   '''
   Convert normal B-splines to natural B-splines (bk to bk_nat).

   >>> x = np.r_[-2:12:0.1]
   >>> B, k = cbspline_Bk(x, 10, 0, 10)
   >>> bk2bknat(B, k, 10)
   >>> gplot(x, B, ', "" us 1:3, "" us 1:4, "" us 1:5')

   '''
   idx, = np.where(kk==0)
   if idx.size:
      G[2,idx] -= G[0,idx]
      G[1,idx] += 2 * G[0,idx]
      G[0,idx] = 0
   idx, = np.where(kk==K-2)
   if idx.size:
      G[1,idx] -= G[3,idx]
      G[2,idx] += 2 * G[3,idx]
      G[3,idx] = 0
   idx, = np.where(kk==K-1) # for points at the edge knot
   if idx.size:
      G[0,idx] -= G[2,idx]
      G[1,idx] += 2 * G[2,idx]
      G[2,idx] = 0


class spl:
   '''
   Cardinal cubic spline.

   Examples
   --------
   >>> a = np.zeros(10+2); a[5] = 6
   >>> s = ucbspl(a).to_spl()
   >>> x = np.r_[-2:10:0.1]

   Compare spline and bspline

   >>> gplot(x, s(x), ',', x, s.to_cbspl()(x), ',', s.xk, s())

   '''
   def __init__(self, a, b, c, d, xmin=0., xmax=None):
      self.a = a, b, c, d
      self.K = K = a.size + 1   # there is one more knot than intervals
      self.xmin = xmin
      self.xmax = K if xmax is None else xmax
      self.xk = np.linspace(xmin, self.xmax, num=K)

   def __call__(self, x=None, der=0, border='extrapolate'):
      if x is None:
         # for last knot append end of last interval a+b*1+c*1^2+d*1^3
         return np.r_[self.a[0], np.sum(zip(*self.a)[-1])]
      else:
         x = (self.K-1)/(self.xmax-self.xmin) * (x-self.xmin)
         k, p = divmod(x, 1)
         k = k.astype(np.int)
         # handling of locations outside
         if border == 'extrapolate':
            # further options could be natural
            ii, = np.where(k < 0)
            p[ii], k[ii] = x[ii], 0
            ii, = np.where(k > self.K-2)
            p[ii], k[ii] = x[ii]-(self.K-2), self.K-2
         # use Horner schema y = a + x (b+ x (c + xd)) in a slightly different way
         y = self.a[-1][k]   # init with d (not 0.)
         for ci in self.a[-2::-1]:
            y *= p; y += ci[k]
         return y

   def osamp(self, ofac=1):
      '''
      Compute an oversampled spline curve.

      Parameters
      ----------
      ofac : oversampling factor.

      Returns
      -------
      x : ndarray
         Oversampled postions.
      y : ndarray
         Spline values.
      '''
      x = np.linspace(self.xmin, self.xmax, num=self.K*ofac)
      return x, self(x)

   def to_cbspl(self):
      '''
      Coefficient transformation from normal spline to Bspline.

      The (K-1)x4 spline coefficients cover K-1 intervals and K knots.
      Here the K+2 B-spline coefficients are computed.
      '''
      a, b, c, d = self.a
      a = np.r_[a[0] - b[0] + 2/3*c[0],
                a - 1/3*c,
                a[-2:] + 2*b[-2:] + 11/3*c[-2:] + 6*d[-2:]]
      return ucbspl(a, self.xmin, self.xmax)


class ucbspl:
   '''
   Uniform cubic B-spline evaluations.

   Parameters
   ----------
   a : array_like
      B-spline coefficients.

   Examples
   --------

   >>> a = np.zeros(10+2); a[5] = 6
   >>> cs = ucbspl(a)
   >>> cs()
   array([0., 0., 0., 1., 4., 1., 0., 0., 0., 0.])
   >>> cs(cs.xk)
   array([0., 0., 0., 1., 4., 1., 0., 0., 0., 0.])
   >>> cs(5), cs([4.3,5.5])
   (array(1.), array([3.541, 0.125]))
   >>> x = np.r_[:10:0.1]
   >>> gplot(x, cs(x), 'w lp,', cs.xk, cs(), 'pt 7 lt 3')

   Convert to a normal spline and get the coefficients.

   >>> #cs.to_spl()
   (array([0., 0., 0., 1., 4., 1., 0., 0., 0.]), array([ 0.,  0.,  0.,  3.,  0., -3.,  0.,  0.,  0.]), array([ 0.,  0.,  0.,  3., -6.,  3.,  0.,  0.,  0.]), array([ 0.,  0.,  1., -3.,  3., -1.,  0.,  0.,  0.]))

   Interpolant

   >>> y = np.zeros(10); y[5] = 1.
   >>> cs = ucbspl_fit(y)
   >>> xx = np.r_[-2:10:0.01]
   >>> gplot(xx, cs(xx), 'w l,', cs.xk, cs(), 'pt 7 lt 1')

   '''
   def __init__(self, a, xmin=0., xmax=None):
      self.a = np.array(a)
      self.K = K = self.a.size - 2
      self.xmin = xmin
      self.xmax = K-1 if xmax is None else xmax
      # knot positions (uniform knot grid)
      self.xk = np.linspace(xmin, self.xmax, num=K)

   def __call__(self, x=None, der=0):
      a = self.a
      if x is None:
         # simplified evaluation for knots only
         return 1./6 * (a[:-2] + 4*a[1:-1] + a[2:])
      else:
         B, kk = cbspline_Bk(x, self.K, self.xmin, self.xmax)
         # spline evalution  y_i = sum_k a_k*B_k(x_i)
         y = 0.
         for k in [0, 1, 2, 3]: y += a[kk+k] * B[k]
         # y = np.array([np.dot(G[i],a[kki:kki+4]) for i,kki in enumerate(kk)])
         # y = np.sum(G * a.k[kk[:,np.newaxis]+np.arange(4).T], axis=1)
         # y = np.einsum('ij,ij->i', G, a.k[kk[:,np.newaxis]+np.arange(4).T])
         return y.reshape(kk.shape)

   def osamp(self, ofac=1):
      '''
      Compute an oversampled spline curve.

      Parameters
      ----------
      ofac : oversampling factor.

      Returns
      -------
      x : ndarray
         Oversampled postions.
      y : ndarray
         Spline values.
      '''
      x = np.linspace(self.xmin, self.xmax, num=self.K*ofac)
      return x, self(x)

   def to_spl(self):
      '''
      Convert to cardinal cubic spline.

      The K+2 cubic B-spline coefficients are converted to spline cofficients
      in the K-1 intervals.
      '''
      a = self.a 
      c0 = 1./6 * (a[:-3] + 4*a[1:-2] + a[2:-1])   # = yk
      c1 = 0.5 * (-a[:-3] + a[2:-1])
      c2 = 0.5 * (a[:-3] - 2*a[1:-2] + a[2:-1])
      c3 = 1./6 * (-a[:-3] + 3*a[1:-2] - 3*a[2:-1] + a[3:])
      return spl(c0, c1, c2, c3, self.xmin, self.xmax)

   def dk(self):
      a = self.a
      return 6*a[:-2] - 12*a[1:-1]+ 6*a[2:]


class v_cspl:
   '''
   Variance prediction with the covariance matrix.
   '''
   def __init__(self, cov, xmin=0., xmax=None):
      self.cov = cov
      self.xmin = xmin
      self.xmax = xmax
      self.K = self.cov.shape[0] - 2
   def __call__(self, x=None):
      B, kk = cbspline_Bk(x, self.K, self.xmin, self.xmax)
      v_f = 0.
      for k in range(4):
         for j in range(4): 
            v_f += B[k] * self.cov[kk+k,kk+j] * B[j]
      return v_f.reshape(kk.shape)


def ucbspl_fit(x, y=None, w=None, K=10, xmin=None, xmax=None, lam=0., pord=2, mu=None, e_mu=None, nat=True, retfit=False, var=False, e_yk=False, cov=False, plot=False, c=True):
   '''
   Fit a uniform cubic spline to data.

   Parameters
   ----------
   x : array_like
      Data position.
   y : array_like
      Data values.
   w : array_like
      Data weights (1/e_y^2).
   K : integer
      Number of uniform knots.
   xmin : float
      Position of the first knot (default minimum of x).
   xmax : float
      Position of the last knot (default maximum of x).
   lam : float
      Penality, smoothing value (default 0, i.e. spline regression).
   pord : int
      Penality order (default 2, i.e. second derivative/curvature).
   mu : array_like
      Analog to mean value of a Gaussian process.
   e_mu : array_like
      Deviation of mu.
   nat : boolean
      Natural spline. (Unclear how to treat this for pord!=2.)
   e_yk : boolean
      Error estimates for the knot values.
      In unweighted case, it estimates the np.sqrt((N-1)/(K-1)).
   retfit: boolean
      If true the returned tuple contains the predicted values y(x) at the data position x.
   var : boolean
      Variance of the model prediction.
   cov : boolean
      If false band matrices are used for efficient solution of equation system with band solver.
      If true covariances are estimated using matrix inversion.

   Returns
   -------
   ucbspl
      The spline model.
   ucbspl, yfit, varmod, v_cspl
      The content of tuple depends on keywords retfit, var, cov.

   Examples
   --------
   >>> x = np.r_[0:100., 400:500, 600:700]
   >>> xx = np.r_[0:x[-1]:0.1]
   >>> y = (np.sin(0.1*x)+1)*1000

   Comparison of simple and proper variance estimation

   >>> spl, v_spl, c = ucbspl_fit(x, y, w=1./1000**2, K=50, lam=0.00000001, pord=1, nat=0, var=True, cov=True)
   >>> gplot(xx, spl(xx), np.sqrt(v_spl(xx)), np.sqrt(c(xx)), ' us 1:($2-$3):($2+$3) with filledcurves fill transparent solid 0.5, "" us 1:($2-$4):($2+$4) with filledcurves fill transparent solid 0.5 lt 3, ""  w l lt 7,', spl.xk, spl(), ' lt 7 ps 0.5 pt 7,', x, y, ' pt 7 lt 1')

   >>> spl, v_spl = ucbspl_fit(x, y, w=0.00001, K=50, mu=800, e_mu=800, nat=0, var=True)
   >>> gplot(x, y, ' pt 7')
   >>> ogplot(xx, spl(xx), np.sqrt(v_spl(xx)), '  w l, "" us 1:($2-$3):($2+$3) with filledcurves fill transparent solid 0.5')

   Comparison of penality order

   >>> gplot(x, y, ' pt 7')
   >>> for e_mu in [5, 1, 0.1]:
   ...    spl = ucbspl_fit(x, y, K=50, mu=500, e_mu=e_mu)
   ...    ogplot(xx, spl(xx), 'w l t "mu=0.5 +/- %s"'%e_mu)
   >>> for lam in [0.1, 1, 10, 100, 10000]:
   ...    spl = ucbspl_fit(x, y, K=50, lam=lam, pord=1, nat=0)
   ...    ogplot(xx,spl(xx), 'w l t "lam=%s",' % lam, spl.xk, spl(), ' ps 0.3 pt 6 lt 3 t ""')
   >>> for lam in [0.0000001, 0.01, 0.1, 1, 10, 100000]:
   ...    spl = ucbspl_fit(x, y, K=50, lam=lam, pord=2, nat=0)
   ...    ogplot(xx,spl(xx), 'w l t "lam=%s",' % lam, spl.xk, spl(), ' ps 0.3 pt 6 lt 3 t ""')
   >>> for lam in [1/0.1**2, 1, 1/5.**2]:
   ...    spl = ucbspl_fit(x, y, K=50, lam=lam, pord=0, nat=0)
   ...    ogplot(xx, spl(xx), 'w l t "lam=%s"' % lam)

   '''
   if y is None:    # uniform data
      y = x
      x = np.linspace(xmin or 0, y.size-1 if xmax is None else xmax, num=y.size)
   if w is None:
      w = 1
      wy = y
   else:
      wy = w * y
   wy = np.ascontiguousarray(wy)

   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   G, kk = [cbspline_Bk, _cbspline_Bk][c](x, K, xmin, xmax)

   if nat:
      Gorg = 1 * G
      bk2bknat(G, kk, K)

   nk = K + 2
   BTy = np.zeros(nk)
   BTBbnd = np.zeros((4, nk))

   if 0:
      # the slow way without band storage
      print 'bspline2'
      B = bspline2((x-x.min())/(x.max()-x.min())*K, K, D=3)
      n = B.shape[1]
      D = np.diff(np.eye(n),n=pord).T
      print 'BTy'
      BTy = np.dot(B.T, y)
      print 'BTB'
      BTB = np.dot(B.T,B) + lam*np.dot(D.T,D)
      #BTB = B.T.dot(B)+lam*D.T.dot(D)
      #a = np.linalg.solve(BTB, BTy) # too slow
      a = solve_banded((3,3), bandmat(BTB,bw=7), BTy)
      print a
      yfit = np.dot(B,a)

   if c:
      # C version for band matrix
      _cbspline.rhs_fill(BTy, G, wy, kk, kk.size)
      # w.size=1 not broadcasted in lhsbnd_fill 
      _cbspline.lhsbnd_fill(BTBbnd, G, w if np.size(w)>1 else np.full_like(y, w), kk, kk.size, nk)
      if nat:
         BTy = 1*BTy[1:K+1]
         BTBbnd = 1*BTBbnd[:,1:K+1]
      if lam:
         _cbspline.add_DTD(BTBbnd, BTy.size, lam, pord)
   else:
      # Python version for band matrix
      # compute sum_jk B_j(x_i) B_k(x_i) w_i
      for k in range(4):
         BTy += np.bincount(kk+k, wy*G[k], nk)
         #for j in range(k+1): 
            #BTBbnd[3-k+j] += np.bincount(kk+k, w*G[k]*G[j], nk) # fill upper
         #for j in range(k+1): 
            #BTBbnd[k-j] += np.bincount(kk+j, w*G[k]*G[j], nk)  # fill lower
         for j in range(k,4):
            BTBbnd[j-k] += np.bincount(kk+k, w*G[k]*G[j], nk)  # fill lower

      if nat:
         BTy = BTy[1:K+1] * 1
         BTBbnd = BTBbnd[:,1:K+1] *1

      if lam:
         # Add penalty lam*DTD
         print lam
         if pord == 0:
            # diagonal
            BTBbnd[0] += lam
         elif pord == 1:
            #  1  2 ...  2  2  1
            # -1 -1 ... -1 -1
            # diagonal
            BTBbnd[0,[0,-1]] += lam
            BTBbnd[0,1:-1] += 2*lam
            # subdiagonals
            BTBbnd[1,:-1] -= lam
         elif pord == 2:
            #   1  5  6 ...  6  5  1
            #  -2 -4 -4 ... -4 -2
            #   1  1  1 ...  1
            # diagonal
            BTBbnd[0,[0,-1]] += lam
            BTBbnd[0,[1,-2]] += 5 * lam
            BTBbnd[0,2:-2] += 6 * lam
            # first subdiagonal
            BTBbnd[1,[0,-2]] -= 2*lam
            BTBbnd[1,1:-2] -= 4 * lam
            # second subdiagonal
            BTBbnd[2,:-2] += lam
         else:
            # generic version
            D = np.diff(np.eye(nk), n=pord)
            DTD = lam * np.dot(D,D.T)
            for k in range(4):
               BTBbnd[k,:-k] += np.diag(DTD, k)

   #ds9(BTBbnd)
   #pause()

   if mu is not None and e_mu:
      # GP like penality with mu and variance 
      BTy += mu / e_mu**2
      BTBbnd[0] += 1. / e_mu**2

   if cov:
      # invert matrix
      BTB = np.diag(BTBbnd[0])
      for i in range(1,4): BTB += np.diag(BTBbnd[i,:-i], i) + np.diag(BTBbnd[i,:-i], -i)
      covmat = np.linalg.inv(BTB)
      BTy = covmat.dot(BTy)
   else:
      if c: # python version
         #BTB = np.diag(BTBbnd[0])
         #tt = BTy*1
         #for i in range(1,4): BTB += np.diag(BTBbnd[i,:-i], i) + np.diag(BTBbnd[i,:-i], -i)
         #ds9(BTB)
         #print _cbspline.cholsol(BTB, tt, BTy.size)
         _cbspline.cholbnd(BTBbnd, BTy, BTy.size, 3)
         #_cbspline.cholbnd_upper(gg, uu,  BTy.size, 3)
      elif 1: # python version
         solveh_banded(BTBbnd, BTy, lower=True, overwrite_ab=True, overwrite_b=True)
      else:
         # old c version with Gauss elimination
         # symmetric part
         BTBbnd[4,:-1] = BTBbnd[2,1:]
         BTBbnd[5,:-2] = BTBbnd[1,2:]
         BTBbnd[6,:-3] = BTBbnd[0,3:]
         #pause()
         _cbspline.bandsol(BTBbnd, BTy, BTy.size, 7)
         # only works with "1*" (copy needed: np.ctypeslib.as_ctypes(BTBbnd): TypeError: strided arrays not supported)

   a = BTy
   if nat:
      G = Gorg
      a = np.r_[2*a[0]-a[1], a, 2*a[K-1]-a[K-2]]

   mod = ucbspl(a, xmin, xmax)
   out = mod
   if retfit or var or cov:
      out = out,

   if retfit:
      # same as yfit = mod(x), but re-using G[k]
      yfit = 0.
      for k in [0,1,2,3]: yfit += a[kk+k] * G[k]
      out += yfit,

   # error estimation for knot coefficients ak
   #    sig(ak)^-2 = sum_i Bk(x_i) * w_i
   if var or e_yk:
      wa = 0.  # = sum_i Bk(x_i) * w_i =  1/var(ak)
      for k,Gk in enumerate(G): # short loop, but bincount could also have overhead
         wa += np.bincount(kk+k, Gk*w, nk)
      if mu is not None and e_mu:
         wa +=  1. / e_mu**2
      with np.errstate(divide='ignore'): # yet lam is not handled
         vara = 1. / wa
      varmod = ucbspl(vara, xmin, xmax)
      if var:
         out += varmod,

   if cov:
      # compute   sum_i,j B_j(x) * cov_jk B_k(x)
      out += v_cspl(covmat, xmin, xmax),
      #gplot(x, v_f(x), varmod(x), ', "" us 1:3,', varmod.xk, varmod(), v_f(varmod.xk),', "" us 1:3')

   # error estimation for knot values
   # var(yk) = 1/6 var(ak-1) + 4/6 var(ak) + 1/6 var(ak+1)
   if e_yk:
      mod.e_yk = np.sqrt(varmod())

   if plot:
      gplot(mod.osamp(20), ' w l lt 1,',  # oversample the knots
            mod.xk, mod(), ' lt 1 pt 7,',
            x, y, mod(x), ' lt 3, "" us 1:3 w l lt 2')

   return out


### The following functions are more for demonstration, rather than efficient use.
def SolveBanded(A, D, bw=3):
    # Find the diagonals
    ud = np.insert(np.diag(A,1), 0, 0) # upper diagonal
    d = np.diag(A) # main diagonal
    ld = np.insert(np.diag(A,-1), len(d)-1, 0) # lower diagonal
    # simplified matrix
    ab = np.matrix([
        ud,
        d,
        ld,
    ])
    #ds9(ab)
    return solve_banded((1, 1), ab, D )

def bandmat(A, bw=3):
   '''Convert a matrix to a band matrix.'''
   Abnd = np.zeros((bw, A.shape[0]))
   Abnd[bw/2,:] = np.diag(A)
   for i in range(1, bw/2+1):
       Abnd[bw/2-i,i:] = np.diag(A, i) # upper diagonal
       Abnd[bw/2+i,:-i] = np.diag(A, -i) # lower diagonal
   return Abnd


def bspline(x, k, d=3):
   '''
   Uniform B-spline with Cox-de Boor recursion formula.
   Python implementation

   Parameters
   ----------
   x : float, interpolation point
   k : integer, knot
   d : integer, degree

   Examples
   --------
   >>> bspline(7, 5)
   0.6666666666666666
   >>> bspline(1.5, 1, d=0)
   1.0
   >>> bspline(1.5, 1)
   0.020833333333333332
   '''
   if d>0:
      return (x-k)/d*bspline(x,k,d-1) + (k+d+1-x)/d*bspline(x,k+1,d-1)
   return float(k<=x<k+1)
   #return 1.0 if k<=x<k+1 else 0.0

def bspline1(x, K, D=3):
   '''
   Uniform B-spline with for-loop for a scalar x.
   Numpy implementation.

   Parameters
   ----------
   x : float, scalar (x <= K)
   K : integer, number of knots
   D : integer, Bspline degree

   '''
   #dx = (xr-xl)/nk
   k = np.arange(K+D-1)
   p = x - k
   B = ((0<=p) * (p<1)).astype(float)  # =Bk0
   for d in range(1, D+1):
#      B = p[:-d]*B[:-1]/d + (d+1-p[:-d])/d*B[1:]         # Bkd
      B = (p[:-d]*B[:-1] - (p[d:]-1)*B[1:])/d         # Bkd
      print B
   print zip(k,B)
   return B

def bspline2(x,K,D=3):
   '''
   Uniform B-spline with for-loop and vectorized.

   Parameters
   ----------
   x - vector
   K - number of knots
   D - Bspline degree
   (K knots, K-1 intervals each covered by D+1 Bsplines, K+2D construction knots for Bsplines
   K+D-1 regression bsplines

   Returns
   -------
      array(K+D,x.size)

   Example
   -------
   >>> x = np.arange(40)/40.0+0.1
   >>> gplot(x, bspline2(x, 10, D=3))
   '''
   k = np.arange(-D, K+D)  # knot number
   p = x-k[:,np.newaxis]
   B = ((0<=p) * (p<1)).astype(float)  # =Bk0
   for d in range(1,D+1):  # slow B contains many zero
      #pause(B[:,-1].size,d)
      B = (p[:-d]*B[:-1] + (1-p[d:])*B[1:])/d         # Bkd
   #print zip(k,B)
   #ds9(B)
   #pause(d)
   return B.T


def cbspline(x, k):
   '''
   Uniform cubic B-spline with direct computation for scalar.
   python implementation

   Parameters
   ----------
   x : interpolation point
   k : knot

   '''
   if x<k: return 0
   x = x-k
   if x<1: return x**3/6                       # 1/6*(x-k)^3                                     if k   <= s < k+1
   if x<2: x = 2-x; return ((x-2)*3*x*x+4)/6   # 1/6*[-3(x-k-1)^3 + 3(x-k-1)^2 + 3(x-k-1) + 1]   if k+1 <= s < k+2
   if x<3: x -= 2; return ((x-2)*3*x*x+4)/6    # 1/6*[3(x-k-2)^3 - 6(x-k-2)^2 + 4]               if k+2 <= s < k+3
   if x<4: return (4-x)**3/6                   # 1/6*[1-(x-k-3)]^3                               if k+3 <= s < k+4
   return 0                                    # otherwise

def cbspline_v00(x,k):
   '''
   Uniform cubic B-spline direct
   python implementation
   x - interpolation point
   k - knot
   d - degree
   '''
   if x<k: return 0
   x = x-k
   if x<1: return x**3/6                   # 1/6*(x-k)^3                                     if k   <= s < k+1
   x -= 1
   if x<1: return (((-x+1)*x+1)*3*x+1)/6   # 1/6*[-3(x-k-1)^3 + 3(s-k-1)^2 + 3(s-k-1) + 1]   if k+1 <= s < k+2
   x -= 1
   if x<1: return ((x-2)*3*x**2+4)/6       # 1/6*[3(x-k-2)^3 - 6(s-k-2)^2 + 4]               if k+2 <= s < k+3
   x -= 1
   if x<1: return (1-x)**3/6               # 1/6*[1-(x-k-3)]^3                               if k+3 <= s < k+4
   return 0                                # otherwise

def Bspline(x, k, d=3):
   '''uniform B-spline with Cox-de Boor recursion formula
   x - interpolation vector
   k - knot
   d - degree
   python implementation

   Examples
   --------
   gplot(t, Bspline(t,0,d=1), ' w l,', t, Bspline(t,1,d=1), ' w l')
   ogplot(t, Bspline(t,0,d=2), ' w l,', t, Bspline(t,1,d=2), ' w l')
   ogplot(t, Bspline(t,0,d=3), ' w l,', t, Bspline(t,1,d=3), ' w l')
   ogplot(t, Bspline(t,0,3)+Bspline(t,1,3), ' w l')

   '''
   return np.array([bspline(xi,k,d) for xi in x])

def cBspline(x,k):
   '''uniform B-spline with Cox-de Boor recursion formula
   x - interpolation vector
   k - knot
   d - degree
   python implementation
   '''
   return np.array([cbspline(xi,k) for xi in x])



def example():
   # a noisy time serie
   x = np.arange(10490)/40.0 + 0.1
   yorg = np.sin(x*1.65) + 0.1*np.sin(x*1.65*6.32)
   y = yorg + 0.021*np.random.randn(x.size)

   mod, yfit = ucbspl_fit(x, y, w=1+0*y, K=1800, lam=0.11, retfit=True)

   gplot(x, y, yorg, ' lt 2,"" us 1:3 w l  lt 1 t "org"')
   ogplot(x, yfit, yfit-yorg,'w l,"" us 1:3')


   A = np.array([
    [20, -5,  0,  0],
    [-5, 15, -5,  0],
    [ 0, -5, 15, -5],
    [ 0,  0, -5, 10]
   ])
   D = np.array([1100, 100, 100, 100])
   print "SciPy - Solve Banded"
   print "A:", A
   print "D:", D
   print "Result:", SolveBanded(A, D)
   pause()

if __name__ == "__main__":
   #pause()
   import doctest
   doctest.testmod()
   #exec(doctest.script_from_examples(ucbspl_fit.__doc__))
   #exec(doctest.script_from_examples(ucbspl.__doc__))
   #exec(doctest.script_from_examples(cbspline_Bk.__doc__))
   #exec(doctest.script_from_examples(bk2bknat.__doc__))
   #exec(doctest.script_from_examples(bspline.__doc__))

   #import pdb;  pdb.set_trace() #run("pass", globals(), locals()); #pdb.set_trace()
   t = np.arange(400)/40.0
   Bspline(t,0,1)

   pause()
   example()


'''
      x    x    x    x    x    x    x    x
     /|\  / \  / \  / \  / \  / \  / \  /|\
    / | \/   \/   \/   \/   \/   \/   \/ | \
   /  | /\   /\   /\   /\   /\   /\   /\ |  \
  /   |/  \ /  \ /  \ /  \ /  \ /  \ /  \|   \
 +----+----+----+----+----+----+----+----+----+
-1    0    1    2    3    4    5    6    7    8

'''

