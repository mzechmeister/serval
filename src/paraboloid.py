#from __future__ import print_function

import numpy as np
from gplot import *
from pause import *
#def pause():
#   pass

xzip = zip; zip = lambda *x: list(xzip(*x))

class paraboloid():
   '''
   Create and evaluate the paraboloid.

# init: p = paraboloid(W); paraboloid(W, Xc, zc)
# call: z = p(*X)
# p.W, p.matrix, p.c paraboloid coeffs
# p.center (root)
# p(p.center), p.max, p.min, p.extremum

   Parameters
   ----------
   W : Symmetric matrix describing the paraboloid (offsets, curvature).
   X : List of parameter vectors to evalute the paraboloid.

   Methods
   -------
   __call__ : Returns list of paraboloid values.

   Attributes
   ----------
   dim : dimension

   Notes
   -----
   The symmetric matrix W for the paraboloid is
      z = W00 + W01*x1 + W02*x2+ W11*x1**2 + W22*x2**2 + 2*W12*x1*x2
        = Wmn xmi * xni
        = Xi W Xi.T
   So the first row/column contains the constant and linear terms and the rest
   the quadratic and cross-terms.

   Examples
   --------
   Create a parabola (1D paraboloid) with minimum value 10
   and curvature 3 centered at 5:
   z(x) = 10 + 3 (x-5)^2
        = 10 + 3*5^2 - 2*3*5x + 3x^2

   >>> z = paraboloid([[3]], xc=[5], zc=10)
   >>> z.dim
   1
   >>> z.W
   array([[ 85, -15],
          [-15,   3]])

   Check function value at minimum

   >>> z(5)
   array([10])

   Create a symmetric 3x3 matrix and a paraboloid of dimension two

   >>> W = np.diag([0, 0.05**-2, 0.3**-2])
   >>> W[0,1] = W[1,0] = 100
   >>> z = paraboloid(W)
   >>> z.dim
   2

   Create sampling points and display.

   >>> x = np.random.rand(10000)*2-1
   >>> y = np.random.rand(10000)*20-10
   >>> gplot(x,y, z(x,y), 'us  1:2:3 palette')

   '''
   def __init__(self, W, xc=None, zc=0):
      '''
      Parameters
      ----------
      W : full paraboloid matrix or curvature matrix if xc is specified
      xc : vector of paraboloid extremum
      zc : paraboloid value at xc
      '''
      W = np.array(W)

      if not np.array_equal(W, W.T):
         W = (W.T + W) / 2
         # raise Exception('Matrix must be symmetric.')

      if xc is not None:
         # the reversion operation, shift to xc
         W0 = -W.dot(xc)
         z0 = zc + np.diag(W).dot(np.square(xc))
         W = np.vstack((np.hstack((z0, W0)),
                        np.vstack((W0, W )).T))
         # W = np.block([[z0, W0],[W0[:,np.newaxis], W]])
      self.W = W
      self.dim = len(W) - 1

   def __call__(self, *X):
      if self.dim == len(X):
         # for de-centered paraboloid prepend dimension
         X = np.vstack((X[0]*0+1, X))

      z = np.einsum('ik,ij,jk->k', X, self.W, X)
      return z

   def center(self, xc=None):
      '''
      Extracts the curvature matrix and the center of the paraboloid

      >>> W = np.diag([0, 0.05**-2, 0.3**-2])
      >>> W[0,1] = W[1,0] = 100
      >>> W
      array([[  0.        , 100.        ,   0.        ],
             [100.        , 400.        ,   0.        ],
             [  0.        ,   0.        ,  11.11111111]])
      >>> Wc, xc = paraboloid(W).center()
      >>> xc
      array([-0.25,  0.  ])
      >>> paraboloid(Wc, xc=xc).W
      array([[ 25.        , 100.        ,  -0.        ],
             [100.        , 400.        ,   0.        ],
             [ -0.        ,   0.        ,  11.11111111]])

      '''
      W = self.W
      if xc is None:
         # center to minimum
         xc = np.linalg.solve(W[1:,1:], -W[0,1:])
         return W[1:,1:]*1, xc
      else:
         # the reversion operation, shift to xc
         W0 = -W.dot(xc)
         Wraw = np.hstack((np.r_[0,W0][:,np.newaxis], np.vstack((W0, W))))
         return  Wraw


def fit_paraboloid(X, z, offset=True):
   '''
   Fit a paraboloid z(X).

   This can be used approximate the extrema of chi2 and likelihood samples.

   Parameters
   ----------
   X : position vectors
   z : values
   offset : if true, the offset (Xc,zc) is fitted

   Returns
   -------
   Best fitting fit_paraboloid z(X) = sum_mn Wmn * xm * xn using linear least square.

   Example
   -------
   Create sampling points

   >>> x = np.random.rand(10000)*2-1
   >>> y = np.random.rand(10000)*20-10
   >>> X = x, y

   Create paraboloid matrix

   >>> W = np.array([[1/0.05**2, 0.5*1/0.05/0.3],  [0.5*1/0.05/0.3,  1/0.3**2]])
   >>> #V = np.array([[0.05**2, 0.5*0.05*0.3],  [0.5*0.05*0.3,  1*0.3**2]])
   >>> #W = np.linalg.inv(V)
   >>> #z_i = W11 *x1_i**2 + W22 * x2_i**2 + 2* W12 * x1_i * x2_i
       #    = Wmn xmi * xni = Xi W Xi.T
   >>> z = paraboloid(W, xc=[-0.5, 1])
   >>> gplot(x,y,z(*X), 'us 1:2:3 palette')
   >>> W
   array([[400.        ,  33.33333333],
          [ 33.33333333,  11.11111111]])
   >>> zmod = fit_paraboloid([x-0.5,y+2], z(x,y), offset=False)
   >>> #gplot+(x, y, zmod([x-0.5,y+1]), ' us 1:2:3 palette pt 6')
   >>> zmod = fit_paraboloid(X, z(x,y))
   >>> gplot+(X, zmod(*X), ' us 1:2:3 palette pt 6')
   >>> W; zmod.center()
   array([[400.        ,  33.33333333],
          [ 33.33333333,  11.11111111]])
   (array([[400.        ,  33.33333333],
          [ 33.33333333,  11.11111111]]), array([-0.5,  1. ]))

   Model:
     z = X W X.T
     data z_i, X
     min sum (z_i - z)^2 = sum (z_i - X_i W X_i)
     Cmn : 0 = sum (z_i - X_inm Wnm X_inm) * X_i
           z_i X_i =  X_i W X_i X_i

   '''
   X = np.array(X)
   if offset:
      # prepend a column with ones
      X = np.vstack((X[0]*0+1, X))
   k = np.argmin(z)

   dX = X #- X[...,k][:,np.newaxis]
   dz = z #- z[k]

   n = len(dX) 
   idx_triu = np.triu_indices(n)   # indices for the upper-triangle of an (n, m) array

   # create basis x^2 xy y^2 ...
   GG = (dX * dX[:,np.newaxis,:])[idx_triu]
   # RHS
   xLx = GG.dot(dz)
   # LHS
   xxxx = GG.dot(GG.T)
   Hu = xLx.dot(np.linalg.inv(xxxx))  # or np.linalg.solve(xxxx, xLx)
   H = np.zeros((n, n))
   H[idx_triu] = Hu
   H = (H + H.T) / 2

   return paraboloid(H)


def covmat_fit(X, z, **kwargs):
   '''
   Estimate covariance matrix.

   Parameters
   ----------
   X : samples
   z : likelihood values

   Example
   -------
   >>> x = np.random.rand(1000)*2-1
   >>> y = np.random.rand(1000)*20-10
   >>> X = x, y
   >>> V = np.array([[0.05**2, 0.5*0.05*0.3],  [0.5*0.05*0.3, 0.3**2]])
   >>> W = np.linalg.inv(V)
   >>> z = paraboloid(W, xc=[0,0])
   >>> cov = covmat_fit(X, z(*X))
   >>> cov.Va
   array([[0.0025, 0.0075],
          [0.0075, 0.09  ]])
   >>> cov.e_a
   array([0.05, 0.3 ])
   >>> cov.C
   array([[1. , 0.5],
          [0.5, 1. ]])
   >>> cov.Xc
   array([2.72004641e-17, 3.21409566e-16])

   '''
   H = fit_paraboloid(X, z).W
   #print(np.array2string(Va,precision=2))
   return covmat(H, **kwargs) # Va, e_a


class covmat():
   '''
   Class create cov matrix from full paraboloid matrix.

   Methods
   init: V = covmat(W)
   diag: V.diag V.trace
   sigma: C = V.sigma
   norm: C = V.norm
   corner:

   Attributes
   ----------
   dim :
   N : Number data contributed to W
   inv:  W = V.inv
   p.project, p.marginalise
   contor:
   showC: gnuplot image
   show:

   Notes
   -----
   Covariance matrix:
           s_1^2 s_12   s_1n       s_1^2 g_12*s_1*s_2 g_1n*s_1*s_n
   Cov a = s_21  s_2^2  s_2n   =         s_2^2
           s_n1  s_n2   s_n^2

   Pearson coefficient:
   g_xy = cov(x,y) / sqrt(Var(x)*Var(y))
        = s_xy / (s_x*s_y)

   '''
   def __init__(self, W, N=None):
      self.W = W
      self.p = paraboloid(W)
      self.N = N
      # extract the center position (optimum) and the curvatures matrix
      Vinv, self.Xc = self.p.center()
      self.min = self.p(*zip([1]+ list(self.Xc)))
      self.Va = Va = np.linalg.inv(Vinv)
      if N:
         DOF = N - len(Va)
         chi2red = self.min / DOF        # as in curve_fit
         #chi2red = self.min / (DOF-2.)   # as in np.polyfit
         Va *= chi2red
      self.e_a = np.sqrt(np.diag(Va))
      self.C = (1/self.e_a) * Va * (1/self.e_a[np.newaxis].T)

   def contor(self):
      pass

   def corner(self, histogram=True, labels=True):
      '''
      Pseudo corner plots to visualise the covariance matrix.

      Yep, frequentist can also do those correlation plots,
      albeit only linear correlation which may look boring.

      Parameters
      ----------
      labels : list for axis labels
      '''
      W = self.W
      N = len(W) - 1
      if labels is True:
         labels = ["a_%s"%(i+1) for i in range(N)]

      Va = self.Va
      #gplot.stdout=1
      gplot.term("wxt size 640,640")
      gplot.multiplot()
      gplot.put('''
         $line << EOD
         -1
         1
EOD
         N = %s.
         set tics front
         set lmar 0
         set rmar 0
         set tmar 0
         set bmar 0
         gap = 0.02
         x0 = 0.10
         y0 = 0.08
         set size (1-x0)/N-gap,(1-y0)/N-gap  # should be square
         unset colorbox
      ''' % N)
      gplot.mxtics().mytics()
      gplot.key("noautotitle")
      #gplot.unset("key")

      for i in range(1,N+1):
         # the first column ist plotted first
         j=i
         gplot.put('m=%s; n=%s' %(i,j))
         gplot.origin("x0+(1-x0)*(m-1)/N, y0+(1-y0)*(1-n/N)")
         if i==1:
            gplot.ylabel("'%s'" % labels[0])
         if i>1:
            gplot.ytics('format ""')
            gplot.unset('ylabel')

         if histogram:
            # histogram
            gplot.xtics('format "" rotate')
            if i==N:
               gplot.xtics('format "%g"').xlabel("'%s'" % labels[i-1])

            gplot.var(a=self.Xc[i-1])
            gplot.var(e_a=self.e_a[i-1])
            gplot("[a-2*e_a:a+2*e_a] exp(-(x-a)**2/e_a**2/2) t '%g +/- %g'" % (self.Xc[i-1], self.e_a[i-1]), flush='')
            #gplot<("[a-2*e_a:a+2*e_a] exp(-(x-a)**2/e_a**2/2) t 'asdfasdf'")
            gplot<("$line us (a+e_a):($1>0) w l lc 'grey' dt 3")
            gplot<("$line us (a-e_a):($1>0) w l lc 'grey' dt 3")
            gplot+("$line us (a):($1>0) w l lc 'black' dt 2")
         for j in range(i+1,N+1):
           gplot.put('m=%s; n=%s' %(i,j))
           gplot.origin("x0+(1-x0)*(m-1)/N,y0+(1-y0)*(1-n/N)")
           if i==1:
              gplot.ylabel("'%s'" % labels[j-1])
           gplot.xlabel("'%s'" % labels[j-1] if i==N else '')
           xc, yc = self.Xc[[i-1,j-1]] #W[0,[i,j]]
           rho = self.C[i-1,j-1] #W[i,j] / np.sqrt(W[i,i]*W[j,j])
           #print(xc, yc, rho, self.e_a)
           if j==N:
              gplot.xtics('format "%g"')
              gplot.xlabel("'%s'" % labels[i-1])

           gplot.var(rho=rho, e_b=self.e_a[j-1], b=self.Xc[j-1])
           gplot.put("x(t,g) = e_a*(sqrt(1+g)*cos($1)-sqrt(1-g)*sin($1))/sqrt(2); y(t,g)=e_b*(sqrt(1+g)*cos($1)+sqrt(1-g)*sin($1))/sqrt(2)")
           if 1:
              # map
              gplot.samples(30).isosample(30)
              # marginalise likelihoods
              Sij = np.array([[Va[i-1,i-1], Va[i-1,j-1]],
                              [Va[j-1,i-1], Va[j-1,j-1]]])
              invSij = np.linalg.inv(Sij)
              gplot.put("C11=%s; C22=%s; C12=%s" %  (invSij[0,0], invSij[1,1], invSij[0,1]))
              #gplot.put("C11=%s; C22=%s; C12=%s" %  (W[i,i], W[j,j], W[i,j]))
              #gplot("[a-2*e_a:a+2*e_a][b-2*e_b:b+2*e_b] '++' us 1:2:(C11*($1-a)**2+C22*($2-b)**2+2*C12*($1-a)*($2-b)) w image", flush='')
              gplot.urange("[a-2*e_a:a+2*e_a]").vrange("[b-2*e_b:b+2*e_b]")
              gplot("[a-2*e_a:a+2*e_a][b-2*e_b:b+2*e_b] '++' us 1:2:(C11*($1-a)**2+C22*($2-b)**2+2*C12*($1-a)*($2-b)) w image", flush='')
           if 1:
              # contors
              gplot<("u=1, [0:2*pi] '+' us (a + x($1,rho)):(b+y($1,rho)) w l lc 'white', [0:2*pi] '+' us (a + 2*x($1,rho)):(b+2*y($1,rho)) w l lc 'white'")
           if 1:
              # best fit
              gplot<(xc,yc, " pt 5 lc 'white'")
           if 1:
              # one sigma limits
              gplot<("$line us (a+2*e_a*$1):(b+e_b) w l lc 'white' dt 3")
              gplot<("$line us (a+2*e_a*$1):(b-e_b) w l lc 'white' dt 3")
              gplot<("$line us (a+2*e_a*$1):(b) w l lc 'white' dt 2")
              gplot<("$line us (a+e_a):(b+2*e_b*$1) w l lc 'white' dt 3")
              gplot<("$line us (a-e_a):(b+2*e_b*$1) w l lc 'white' dt 3")
              gplot+("$line us (a):(b+2*e_b*$1) w l lc 'white' dt 2")
           #pause(j,i)
      gplot.unset('multiplot')



if __name__=='__main__':
   import sys
   try:
      case = float(sys.argv[1])
   except:
      case = None
   if not case:
      # test all examples 
      # ok in python 2, but failing in python3 but due to formaatting
      import doctest
      doctest.testmod()
   elif case == 1:
      # chi2map for a linear fit and comparison with curve_fit
      from scipy.optimize import curve_fit

      # model
      trend = lambda x, a, b: a + b*x   # = np.polynomial.polyval(x, a)
      a0, a1 = 2, 2

      # data
      x = np.arange(20)
      e = 2 + np.random.rand(20)
      y = trend(x, a0, a1) + 2*np.random.normal(0, e)   # scale e by 2
      gplot(x, y, e, "w e pt 7")

      # curve_fit
      a_cf, cov_cf = curve_fit(trend, x, y, [0.0, 0.0], e)
      print(a_cf); cov_cf

      # paraboloid 
      # parameter and chi2 samples
      P = [[1,0], [0,0], [-1,0], [0,1], [0,-1], [1,-1]]
      # chi2map
      z = [np.dot((y-trend(x,*p))**2, 1/e**2) for p in P]

      cov = covmat_fit(zip(*P), z, N=len(x))
      #cov = covmat_fit(zip(*P), z)
      gplot+("%s+%s*x w l lc 3"% tuple(cov.Xc))

      print(a0, a1)
      p = zip(a_cf, np.sqrt(np.diag(cov_cf)))
      print("curve_fit:  %s +/- %s;" % p[0], "%s +/- %s" % p[1])
      p = zip(cov.Xc, cov.e_a)
      print("paraboloid: %s +/- %s;" % p[0], "%s +/- %s" % p[1])
      print(cov_cf)
      print(cov.Va)
      f = trend(x, *cov.Xc)
      X = [0*x+1, x]
      # uncertainty in prediction
      varf = np.einsum('ik,ij,jk->k', X, cov.Va, X)
      gplot+(x, f, np.sqrt(varf), 'us 1:($2+$3) w l lc 5, "" us 1:($2-$3) w l lc 5')

      pause()
   elif case == 2:
      # comparison linear fit: polyfit vs linalg vs paraboloid
      a = 20
      b = 200
      x = 17 + np.random.rand(100)*2
      e_y = 10.*np.ones_like(x)
      x -= np.mean(x)
      y = a + b*x + np.random.normal(scale=e_y, size=100)
      (b_pf,a_pf), cov = np.polyfit(x,y, 1, w=1/e_y, cov=True)
      gplot(x, y, e_y, a_pf+x*b_pf, ', "" us 1:4 w l')

      chi2map = [(ai, bi, np.sum((y-(ai+bi*x))**2/e_y**2)) for ai in np.arange(0,40) for bi in np.arange(180,220)]
      ai, bi, chi2map = zip(*chi2map)

      cm = covmat_fit([ai,bi], chi2map)

      # best fit parameters
      print('polyfit    a, b', a_pf, b_pf)
      print('paraboloid a, b', cm.Xc)

      lhs = np.array([1/e_y, x/e_y]).T
      cov_ls = np.linalg.inv(np.dot(lhs.T, lhs))

      print('paraboloid cov\n', cm.Va)
      print('linreg cov\n', cov_ls)

      # with scaling
      scale = np.sqrt((lhs*lhs).sum(axis=0))
      cov_lss = np.linalg.inv(np.dot((lhs/scale).T, lhs/scale)) / np.outer(scale, scale)
      print('linreg cov + scale\n', cov_lss)

      # polyfit with fac = chi2min/DOF where DOF = N-order-2 (!), order = deg + 1
      print('linreg cov + scale + DOF\n', cm.min/(len(x)-2-2.) * cov_lss)
      print('polyfit cov\n', np.flipud(np.fliplr(cov)))

      # gplot.splot(ai, bi, chi2map)
      pause()
   else:
      # Demo of corner
      # gnuplot-python bug: some plot/elements disappear when setting keytitles
      x1 = np.random.rand(1000)*2-1
      x2 = np.random.rand(1000)*20-10
      x3 = np.random.rand(1000)*20-10
      X = x1, x2, x3

      # parameter variances
      s = np.array([0.05, 0.3, 99.])

      # normalised covariance matrix
      C = np.eye(len(s))
      C[0,1] = C[1,0] = 0.08
      C[0,2] = C[2,0] = -0.8
      C[1,2] = C[2,1] = 0.5

      # covariance matrix
      V = C * s * s[np.newaxis].T  # s * C * s[np.newaxis].T not exactly symmetric
      print("V\n", V)

      # weight matrix / inverse covariance matrix
      W = np.linalg.inv(V)
      W = (W + W.T) / 2 # ensure extract symmetry
      z = paraboloid(W, xc=[50, 2, 150])
      print("W\n", W)

      gg = covmat(fit_paraboloid(X, z(*X)).W)
      gg.corner()
      pause()
      #gg.p(*gg.Xc) # not tricky enough
      #gg.p(*zip([1]+ list(gg.Xc)))

