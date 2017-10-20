__author__ = 'Mathias Zechmeister'
__version__ = '2016-08-16'

import numpy as np
from scipy.linalg import solve_banded
from gplot import *
import os.path
from pause import pause

from ctypes import c_int, c_double
import ctypes

ptr = np.ctypeslib.ndpointer

_cbspline = np.ctypeslib.load_library(os.path.join(os.path.dirname(__file__), 'cbspline'), '.')

#_cbspline.polyfit.restype = ctypes.c_int
_cbspline.cbsplcoeff.argtypes = [ ptr(dtype=np.float), # x
                                 ptr(dtype=np.float),  # G
                                 ptr(dtype=np.int),    # kk
                                 c_int]    #n
_cbspline.lhsbnd_fill.argtypes = [ptr(dtype=np.float), # lhsbnd
                                 ptr(dtype=np.float),  # G
                                 ptr(dtype=np.float),  # w
                                 ptr(dtype=np.int),    # kk
                                 c_int,    # n
                                 c_int,    # bw
                                 c_int,    # nk
                                 c_double] #lam
_cbspline.rhs_fill.argtypes = [ptr(dtype=np.float),    # rhs
                                 ptr(dtype=np.float),  # G
                                 ptr(dtype=np.float),  # w
                                 ptr(dtype=np.int),    # kk
                                 c_int]                # n
_cbspline.bandsol.argtypes = [ptr(dtype=np.float),     # rhs
                                 ptr(dtype=np.float),  # G
                                 #ctypes.POINTER(c_int), # n
                                 #ctypes.POINTER(c_int)] # bw
                                 c_int, # n
                                 c_int] # bw


def cbspline_Bk(x, K, xmin=None, xmax=None):
   '''
   Uniform cubic B-spline with direct and vectorized computation and compressed storage.

   Parameters
   ----------
   x : float or array_like
       Data point in units of knots.
   xmin : float
       Position of the first knot.
   xmax : float
       Position of the last knot.

   (K knots, K-1 intervals each covered by D+1 Bsplines, K+2D construction knots for Bsplines
   K+D-1 regression bsplines

   Returns
   -------
   tuple : tuple
       The four non-zero B-spline coefficients and the knots of the first coefficient  (compressed format).

   Examples
   --------
   >>> x = np.arange(100.)/10
   >>> B, k = cbspline_Bk(x)
   >>> gplot(x, B[:,0])

   '''
   x = np.asarray(x)
   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   x = (K-1.)/(xmax-xmin) * (x-xmin)

   kk, pp = divmod(x, 1)
   G = np.empty((4, x.size))

   if 0:
      ### direct computation with p
      G[0] = (1-pp)**3/6     # = ((1-pp)-pp*(1-pp))/2
      G[1] = (3*pp**3 - 6*pp**2 + 4)/6
      G[2] = (-3*pp**3 + 3*pp**2 + 3*pp + 1)/6
      G[3] = pp**3/6         # = ((pp+1)*pp-pp)/2
   elif 0:
      ### direct and symmetric computation with p and q
      qq = 1 - pp
      G[0]= qq**3
      G[1]= 3*pp**3 - 6*pp**2 + 4 # -3*pp**2 (1-pp)  -3*pp**2 + 4
      G[2]= 3*qq**3 - 6*qq**2 + 4 # -3*pp**3 + 3*pp**2 + 3*pp + 1
      G[3]= pp**3
      G /= 6
   else: #faster? untested
      ### direct and symmetric computation with p and q
      G[0] = 1 - pp
      G[2] = -G[0]*G[0];  G[1] = -pp*pp    # -qq^2; -pp^2
      G[0] *= G[2];       G[3] = G[1]*pp   # -qq^3; -pp^3
      G[1:3] += 4.0/6                      # -pp^2+4/6
      G[2] -= 0.5*G[0]                     # -(-qq^3/2)
      G[1] -= 0.5*G[3]

      G[0] *= -1./6
      G[3] *= -1./6

   return np.array(G.T, order='C'), kk.astype(int)

def _cbspline_Bk(x, K, xmin=None, xmax=None):
   '''as cbspline_Bk (c implementation) '''
   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   x = (K-1.)/(xmax-xmin) * (x-xmin)

   G = np.empty((x.size, 4))
   kk = np.empty(x.size, dtype=int)
   _cbspline.cbsplcoeff(x, G, kk, x.size)
   return G, kk

def bk2bknat(G, kk, K):
   '''Convert normal B-splines to natural B-splines (bk to bk_nat)'''
   idx, = np.where(kk==0)
   if idx.size:
      G[idx,2] -= G[idx,0]
      G[idx,1] += 2 * G[idx,0]
      G[idx,0] = 0
   idx, = np.where(kk==K-2)
   if idx.size:
      G[idx,1] -= G[idx,3]
      G[idx,2] += 2 * G[idx,3]
      G[idx,3] = 0
   idx, = np.where(kk==K-1)
   if idx.size:
      G[idx,0] -= G[idx,2]
      G[idx,1] += 2 * G[idx,2]
      G[idx,2] = 0


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


def pspline(x, y, K=10, lam=0, pord=2, w=None, reta=False,retmod=True):
   '''fit penalised uniform cubic Bsplines
   x,y - data to be fitted
   w  - weights of data (1/y_err^2)
   K - number of uniform knots
   lam - penality (default 0, i.e. spline regression)
   pord - penality order (default 2, i.e. second derivative/curvature)
   Returns: the model

   '''
   #B = bspline2(x/(x.max()-x.min())*(K-1),K,D=3)
   #gplot(x/(x.max()-x.min())*(K-1), x*0 )
   if w is None:
      wy = y
   else:
      ww = w.size/np.sum(w)*w
      wy = w*y

   if 0: # the slow way without band storage
      print 'bspline2'
      B = bspline2((x-x.min())/(x.max()-x.min())*K,K,D=3)
      n = B.shape[1]
      D = np.diff(np.eye(n),n=pord).T
      print 'rhs'
      rhs = np.dot(B.T,y)
      print 'lhs'
      lhs = np.dot(B.T,B)+lam*np.dot(D.T,D)
      #lhs = B.T.dot(B)+lam*D.T.dot(D)
      #a = np.linalg.solve(lhs, rhs) # too slow
      a = solve_banded((3,3),bandmat(lhs,bw=7), rhs)
      print a
      ymod = np.dot(B,a)

   if 1: # with band storage lhsbnd
      print 'cbspline_Bk'
      #G0, kk0 = cbspline_Bk((x-x.min())/(x.max()-x.min()).astype(float)*K)
      G, kk = _cbspline_Bk((x-x.min())/(x.max()-x.min()).astype(float)*K)

      nk = K+3
      print 'rhs ',
      rhs = np.zeros(nk+1)
      if 1:
         _cbspline.rhs_fill(rhs,G,wy,kk,kk.size,7)
      else:
         for i,kki in enumerate(kk):
            rhs[kki:kki+4] += G[i]*wy[i] # rhs[-1] = 0 !!!!

      print 'lhs ',
      if 0:
         lhs = np.zeros((nk+1,nk+1))
         if w is None:
            for i,kki in enumerate(kk):
               lhs[kki:kki+4,kki:kki+4] += G[i]*G[i,:,np.newaxis]
         else:
            for i,kki in enumerate(kk):
               lhs[kki:kki+4,kki:kki+4] += G[i]*G[i,:,np.newaxis]*ww[i]
         lhsbnd = bandmat(lhs[:-1,:-1],bw=7).T
      else:
         lhsbnd = np.zeros((7,nk+1))
#         bklind = 3+np.arange(4)+(6*np.arange(4))[:,np.newaxis]
         rr=np.arange(4)
         bklind = (nk+1)*(rr-rr[:,np.newaxis]) + rr[:,np.newaxis]

         if 1 or w is None:
            for i,kki in enumerate(kk): lhsbnd.flat[3*(nk+1)+kki+bklind] += G[i]*G[i,:,np.newaxis]
         else:
            if 0:
               _cbspline.lhsbnd_fill(lhsbnd,G,w,kk,kk.size,7,rhs.size,0.0)
            else:
               for i,kki in enumerate(kk): lhsbnd.flat[3*(nk+1)+kki+bklind] += w[i]*G[i]*G[i,:,np.newaxis]
         if lam!=0:      # penalty with bklind
            hh = np.diff(np.eye(5),n=pord)
            hh = lam*np.dot(hh,hh.T)
            lhsbnd.flat[3*(nk+1)+bklind[:2,:2]] += hh[:2,:2]
            lhsbnd.flat[3*(nk+1)+(nk+1-3)+bklind[:2,:2]] += hh[-2:,-2:]
            lhsbnd[1,2:-1] += 1*lam
            lhsbnd[2,2:-2] -= 4*lam
            lhsbnd[3,2:-3] += 6*lam
            lhsbnd[-3,1:-3] -= 4*lam
            lhsbnd[-2,:-3] += 1*lam
         lhsbnd = lhsbnd[:,:-1]
         #ds9(lhsbnd)
      print 'solving'
      a = solve_banded((3,3), lhsbnd, rhs[:-1])
      a = np.append(a,0)
      print a
      ymod = 0*y
      for i,kki in enumerate(kk):
          ymod[i] = np.dot(G[i],a[kki:kki+4])
   #gplot(np.linspace(x.min(),x.max(),K+1),a[1:-1], 'w l'); ogplot(x,y,ymod, ', "" us 1:3 w l')

   if reta: return ymod,np.linspace(x.min(),x.max(),K+1),a[1:-1]
   return ymod


class spl:
   '''cardinal cubic spline'''
   def __init__(self, a, b, c, d, xmin=0., xmax=None):
      self.a = a, b, c, d
      self.K = a.size
      self.xmin = xmin
      self.xmax = self.K if xmax is None else xmax
      self.dx = float(self.xmax-self.xmin) / (a.size-4)

   def __call__(self, x=None, der=0):
      if x is None:
         return self.a[0]
      else:
         x = (self.K-1.)/(self.xmax-self.xmin) * (x-self.xmin)
         k, p = divmod(x, 1)
         k = k.astype(np.int)
         # use Horner schema y = a + x (b+ x (c + xd)) in a slightly different way
         y = self.a[-1][k]   # init with d (not 0.)
         for ci in self.a[-2::-1]:
            y *= p; y += ci[k]
         return y

   def to_cbspl(self):
      '''coefficient transformation from normal spline to Bspline'''
      a, b, c, d = self.a
      a1 = 1/3.*(3*a[0] - 3*b[0] + 2*c[0])
      a2 = 1/3.*(3*a + 6*b + 11*c + 18*d)
      a = np.append(a1, a - 1/3.*c)
      a = np.append(a, a2[-2:])
      return a

class ucbspl:
   '''
   Uniform cubic B-spline evaluations.

   Parameters
   ----------
   a : array_like
       B-spline coefficients

   Examples
   --------

   >>> a = np.zeros(9+4); a[5] = 6
   >>> cs = ucbspl(a)
   >>> cs()
   array([ 0.,  0.,  0.,  1.,  4.,  1.,  0.,  0.,  0.,  0.])
   >>> cs(cs.xk())
   array([ 0.,  0.,  0.,  1.,  4.,  1.,  0.,  0.,  0.])
   >>> cs(5), cs([4.3,5.5])
   (array(2.875), array([ 3.905191,  1.157125]))
   >>> x = np.arange(100.)/10
   >>> gplot(x, cs(x))

   convert to a normal spline and get the coefficients

   >>> cs.to_spl().a

   '''
   def __init__(self, a, xmin=0., xmax=None):
      self.a =  np.array(a, ndmin=1, copy=0)
      self.K = a.size - 3
      self.xmin = xmin
      self.xmax = self.K if xmax is None else xmax
      self.dx = float(self.xmax-self.xmin) / (a.size-4)

   def __call__(self, x=None, der=0):
      a = self.a
      if x is None:
         # simplified evaluation for knots only
         return 1./6 * (a[:-3] + 4*a[1:-2] + a[2:-1])
      else:
         Bx, kk = cbspline_Bk(x, self.K, self.xmin, self.xmax)
         # spline evalution  y_i = sum_k a_k*B_k(x_i)
         y = 0.
         for m in [0,1,2,3]: y += a[kk+m] * Bx[:,m]
         # ymod = np.array([np.dot(G0[i],a.k[kki:kki+4]) for i,kki in enumerate(kk)])
         # ymod = np.sum(G * a.k[kk[:,np.newaxis]+np.arange(4).T], axis=1)
         # ymod = np.einsum('ij,ij->i', G, a.k[kk[:,np.newaxis]+np.arange(4).T])
         return y.reshape(kk.shape)

   def xk(self):
      # return the knot positions (uniform knot grid)
      return self.xmin + self.dx * np.arange(self.K)

   def to_spl(self):
      '''convert to cardinal cubic spline'''
      # one dummy index, two
      a = self.a
      c0 = 1./6 * (a[:-3] + 4*a[1:-2] + a[2:-1])   # = yk
      c1 = 0.5 * (-a[:-3] + a[2:-1])
      c2 = 0.5 * (a[:-3] - 2*a[1:-2] + a[2:-1])
      c3 = 1./6 * (-a[:-3] + 3*a[1:-2] - 3*a[2:-1] + a[3:])
      return spl(c0, c1, c2, c3, self.xmin, self.xmax)

   def dk(self):
      a = self.a
      return 6*a[:-2] - 12*a[1:-1]+ 6*a[2:]


def _ucbspl_fit(x, y=None, w=None, K=10, xmin=None, xmax=None, lam=0., pord=2, nat=True, retmod=True, reta=False, vara=False, sigk=False, plot=False):
   '''
   Fit a uniform cubic spline to data.

   ymod = _cbspline_fit ( [x, ] y, w=w, K=K, natural=kwnatural, pord=pord, a=a, vara=vara, xout=xout, , chol=chol)

   Parameters
   ----------
   x : array_like
   y : array_like
   w : array_like
      weights
   K : integer
       Number of knots.
   lam : float, smoothing value
   nat : natural spline
   sigk : Error estimates for the knot values.
          In unweighted case, it estimates the np.sqrt((N-1)/(K-1)).

   Examples
   --------
   >>> import pspline
   >>> from pspline import _cbspline_fit
   >>> import numpy as np
   >>> reload(pspline); from pspline import _cbspline_fit; y = np.sin(np.arange(10)*0.1); ymod = _cbspline_fit(np.arange(10.), y, plot=True, nat=True)
   >>> reload(pspline); from pspline import _cbspline_fit; y = np.sin(np.arange(100)*0.1); ymod = _cbspline_fit(np.arange(100.), y, plot=True, nat=True)

   '''
   if y is None:    # regular data
      y = x
      x = np.arange(y.size)
   if w is None:   # not working with lhsbnd_fill w.size=1
      w = np.ones_like(y)
      wy = y
   else:
      ww = w.size/np.sum(w)*w
      wy = w * y
   if lam is None: lam = 0

   if xmin is None: xmin = x.min()
   if xmax is None: xmax = x.max()

   G, kk = _cbspline_Bk(x, K, xmin, xmax)
#   G, kk = _cbspline_Bk(x, xmin, dx=(xmax-xmin)/float(K-1))

   if nat:
      Gorg = 1 * G
      bk2bknat(G, kk, K)

   nk = K + 3
   rhs = np.zeros(nk)
   lhsbnd = np.zeros((7,nk))

   _cbspline.rhs_fill(rhs, G, wy, kk, kk.size)
   _cbspline.lhsbnd_fill(lhsbnd, G, w, kk, kk.size, 7, nk, lam)

   if nat:
      rhs = rhs[1:K+1]
      lhsbnd = lhsbnd[:,1:K+1]
      G = Gorg
   else:
      rhs = rhs[:nk-1]
      lhsbnd = lhsbnd[:,:nk-1]

   if 0: # python version
      solve_banded((3,3), lhsbnd, rhs, overwrite_ab=True, overwrite_b=True)
   else: # c version
      _cbspline.bandsol(1*lhsbnd, rhs, rhs.size, 7)
      #_cbspline.bandsol(1*lhsbnd, rhs, c_int(rhs.size), c_int(7)) # only works with "1*" (copy needed: np.ctypeslib.as_ctypes(lhsbnd): TypeError: strided arrays not supported)

   if nat:
      a = np.hstack((2*rhs[0]-rhs[1], rhs, 2*rhs[K-1]-rhs[K-2], 0.))
   else:
      a = np.hstack((rhs, 0.))

   mod = ucbspl(a, xmin, xmax)
   out = []
   if retmod:
      ymod = np.array([np.dot(G[i],a[kki:kki+4]) for i,kki in enumerate(kk)])
      out = out + [ymod]

   if reta:
      out = out + [mod]

   # error estimation for knot coefficients ak
   #    sig(ak)^-2 = sum_i Bk(x_i) * w_i
   if vara or sigk:
      wB = np.zeros(nk)  # = sum_i Bk(x_i) * w_i =  1/var(ak)
      for i,kki in enumerate(kk):
         wB[kki:kki+4] += w[i] * G[i]
      idx = np.where(wB > 0.)
      vara = ucbspl(np.zeros(nk), xmin, xmax)
      vara.a[idx] = 1. / wB[idx]

   # error estimation for knot values
   # var(yk) = 1/6 var(ak-1) + 4/6 var(ak) + 1/6 var(ak+1)
   if sigk:
      xk = mod.xk() # uniform knot grid
      yk = mod()
      mod.sigk = np.sqrt(vara())

   if plot:
      xk = mod.xk()   # uniform knot grid
      yk = mod()
      xxk = xmin + (xmax-xmin)/(K*20-1) * np.arange(K*20)   # oversample the knots
      yyk = mod(xxk)
      gplot(xxk, yyk, ' w l lt 1,', xk, yk, ' lt 1 pt 7,', x, y, ymod, ' lt 3, "" us 1:3 w l lt 2')

   return out


### The following functions are more for demonstration, rather than efficient use.

def bspline(x, k, d=3):
   '''
   Uniform B-spline with Cox-de Boor recursion formula.
   Python implementation

   Parameters
   ----------
   x : float, interpolation point
   k : integer, knot
   d : integer, degree
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
   >>> bspline2(x, 10, D=3)

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
   print 'done'
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

   ymod, = _ucbspl_fit(x,y,w=1+0*y,K=1800,lam=0.11)

   gplot(x, y, yorg, ' lt 2,"" us 1:3 w l  lt 1 t "org"')
   ogplot(x, ymod, ymod-yorg,'w l,"" us 1:3')


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
   #with open('.pdbrc', 'w') as bp:
        #print >>bp,'break 320'

   #import pdb;  pdb.set_trace() #run("pass", globals(), locals()); #pdb.set_trace()
   bspline(1.5, 1,0)
   bspline(1.5,1,3)
   t = np.arange(400)/40.0
   Bspline(t,0,1)

   from gplot import *
   gplot(t ,Bspline(t,1,1), ' w l,', t, Bspline(t,0,2), ' w l')
   ogplot(t, Bspline(t,1,2), ' w l,', t, Bspline(t,1,2), ' w l')
   ogplot(t, Bspline(t,0,3), ' w l,', t, Bspline(t,1,3), ' w l')
   ogplot(t, Bspline(t,0,3)+Bspline(t,1,3), ' w l')
   #pause()
   example()


'''
     x    x    x    x    x    x    x
    / \  / \  / \  / \  / \  / \  / \
   /   \/   \/   \/   \/   \/   \/   \
  /    /\   /\   /\   /\   /\   /\    \
 /    /  \ /  \ /  \ /  \ /  \ /  \    \
+----+----+----+----+----+----+----+----+
0    1    2    3    4    5    6    7    8

'''

