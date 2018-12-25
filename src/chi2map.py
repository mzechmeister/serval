import numpy as np
try:
   import pyfits
except:
   import astropy.io.fits as pyfits

import cspline as spl
import paraboloid
from gplot import *
from pause import *
from wstat import wsem, wmean


def SSRstat(vgrid, SSR, dk=1, plot='maybe'):
   '''taken from serval.py and modified with SSRv (chi2 minimum)'''
   v_step = vgrid[1] - vgrid[0]
   # analyse peak
   k = SSR[dk:-dk].argmin() + dk   # best point (exclude borders)
   vpeak = vgrid[k-dk:k+dk+1]
   SSRpeak = SSR[k-dk:k+dk+1] - SSR[k]
   # interpolating parabola (direct solution) through the three pixels in the minimum
   a = np.array([0, (SSR[k+dk]-SSR[k-dk])/(2*v_step), (SSR[k+dk]-2*SSR[k]+SSR[k-dk])/(2*v_step**2)])  # interpolating parabola for even grid
   v = (SSR[k+dk]-SSR[k-dk]) / (SSR[k+dk]-2*SSR[k]+SSR[k-dk]) * 0.5 * v_step

   SSR_v = SSR[k] - a[1]**2/a[2]/4. # a1*v+a2*v**2   # position of parabola minimum = a1*v+a2*v**2 = - a[1]^2/2./a[2]+a[1]^2/4./a[2] = -a[1]^2/a[2]/4.
   v = vgrid[k] - a[1]/2./a[2]   # parabola minimum
   e_v = np.nan
   if -1 in SSR:
      print 'opti warning: bad ccf.'
   elif a[2] <= 0:
      print 'opti warning: a[2]=%f<=0.' % a[2]
   elif not vgrid[0] <= v <= vgrid[-1]:
      print 'opti warning: v not in [va,vb].'
   else:
      e_v = 1. / a[2]**0.5
   if (plot==1 and np.isnan(e_v)) or plot==2:
      gplot.yrange('[*:%f]'%SSR.max())
      gplot(vgrid, SSR-SSR[k], " w lp, v1="+str(vgrid[k])+", %f+(x-v1)*%f+(x-v1)**2*%f," % tuple(a), [v,v], [0,SSR[1]], 'w l t "%f km/s"'%v)
      ogplot(vpeak, SSRpeak, ' lt 1 pt 6; set yrange [*:*]')
      pause(v)
   return v, e_v, a, SSR_v


class Chi2Map:
   '''
   Plot, mlRV, and mlCRX from chi2-maps.

   name - file time stamp
   '''
   def __init__(self, chi2map, vrange, RV, e_RV, rv, e_rv, orders, keytitle, rchi=None, No=None, name=''):
      self.args = RV, e_RV, rv, e_rv, orders, keytitle, rchi
      self.No = No # Number of orders
      self.chi2map = chi2map
      self.name = name
      self.info = name.split('/')[-1]
      self.vrange = v_lo, v_step = vrange
      self.vgrid = v_lo + v_step * np.arange(chi2map.shape[1])

      # normalised master map
      self.mCCF = chi2map[orders].sum(axis=0)

      # master CCF from rescaled order-wise CCF
      if rchi is not None:
         #self.rCCF
         self.mCCF = (chi2map[orders]/rchi[orders][:,np.newaxis]**2).sum(axis=0)

      self.nmCCF = self.mCCF / self.mCCF.max()

      self.mlRV, self.e_mlRV, _, _ = SSRstat(self.vgrid, self.mCCF, dk=1, plot=0)
      self.mlRV, self.e_mlRV = self.mlRV*1000, self.e_mlRV*1000

      # Derive again the RV from the CCF minima
      v, e_v, _, self.SSRv = zip(*[SSRstat(self.vgrid, cmo, dk=1, plot=0) for cmo in chi2map[orders]])
      V, e_V = wsem(np.array(v),e=np.array(e_v), rescale=0)
      V, e_V = wsem(np.array(v),e=np.array(e_v))
      if 0:
         print V, e_V
         print RV, e_RV
         gplot(rv, e_rv, 'us 0:1:2 w e,', v, e_v,  'us 0:1:2 w e lc 3, %s,%s' % (RV, V))
         gplot(rv, e_rv, 'us 0:1:2 w e,', v, e_v*rchi[1:][orders],  'us 0:1:2 w e lc 3, %s,%s' % (RV, V))
         pause()

   def plot(self):
      RV, e_RV, rv, e_rv, orders, keytitle, rchi = self.args
      chi2map = self.chi2map[orders]

      # Plot all maps normalised
      gplot.key('Left center rev top title "%s"' % keytitle)
      gplot.palette('defined (0 "blue", 1 "green", 2 "red")')
      gplot.xlabel('"v [km/s]"; set ylabel "chi^2 / max(chi^2)"; set cblabel "order"')
      #gplot(chi2map, ' matrix us ($1*%s+%s):3:2 w l palette'%self.vrange)
      gplot((chi2map/chi2map.max(axis=1)[:,np.newaxis]).T, ' matrix us (%s+$1*%s):3:2 w l palette t "%s"'% (self.vrange +(self.info,)), flush='')

      #ogplot(self.rCCF/self.rCCF.max(), ' us (%s+$0*%s):1 w l lt 4 lw 2 t "total"'%self.vrange, flush='')
      # the master CCF
      ogplot(self.nmCCF, ' us (%s+$0*%s):1 w l lt 7 lw 2 t "lnL"'%self.vrange, flush='')
      # the parabola  for the master
      ogplot('((x-%s)**2/%s**2+%s)/%s lc 7 dt 2 t "lnL (parabola %.5g +/- %.5g m/s)"'% (self.mlRV/1000, self.e_mlRV/1000., self.mCCF.min(), self.mCCF.max(), self.mlRV, self.e_mlRV), flush='')
      ogplot("((x-%s)**2/%s**2+%s)/%s lc 9 dt 3 t 'RV,e\_RV (parabola)'"%(RV, e_RV, self.mCCF.min(), self.mCCF.max()), flush='')

      ogplot([RV]*2, [0,1], [e_RV]*2, 'us 1:2:3 w xerrorlines pt 2 lt 9 t "%.5g +/- %.5g m/s"'% (RV*1000, e_RV*1000), flush='')
      ogplot(rv, (chi2map /chi2map.max(axis=1)[:,np.newaxis]).min(axis=1), 'us 1:2 lt 7 pt 1 t "chi^2_{min}"')
      #pause()

   def mlcrx(self, x, xc, ind):
      '''Maximum likelihood'''
      self.xarg = x, xc, ind

      chi2map = self.chi2map #- self.chi2map.min(axis=1)[:,np.newaxis]
      v_lo, v_step =self.vrange
      vgrid = v_lo + v_step * np.arange(chi2map.shape[1])

      if 0:
         # test data
         oo = ~np.isnan(chi2map[:,0])
         rv[i] = (x-8.6)*0.6*1000 + np.sin(x*100)*100+50
         e_rv[i] = rv[i]*0+10
         rv[i][~oo] =np.nan
         gplot(x,rv[i],e_rv[i], " w e")
         #pval, cov = curve_fit(func, x[ind]-xc, rv[i][ind], np.array([0.0, 0.0]), e_rv[i][ind])
         #crx[i] = pval[1]
         #chi2map = (1./2*(vgrid[:,np.newaxis]-rv[i]/1000)**2/(e_rv[i]/1000)**2).T

      def plot_input(*args):
         RV, e_RV, rv, e_rv, orders, keytitle, rchi = args
         gplot.palette('defined (0 "blue", 1 "green", 2 "red")')
         mm = np.nanmean(chi2map, axis=0)

         gplot.multiplot("layout 1,2")
         # Plot order,RV,chi2map and overplot RV and master map
         # easy to plot as image
         gplot.xlabel('"order"').ylabel('"RV [km/s]"').cblabel('"{/Symbol c^2}"')
         gplot(chi2map, 'matrix us 1:($2*%s+%s):3  w image, '%(v_step, v_lo),
               rv, e_rv, 'us 0:1:2 w e pt 6 lc 7, ',
               mm,mm, 'us (-$1-5):($2*%s+%s):3 matrix w image'%(v_step, v_lo),',',
               [-5.5],[RV],[e_RV], 'us 1:2:3 w e pt 7 lc 7')
         gplot.xlabel('"RV [km/s]"').ylabel('""').cblabel('"order"')
         gplot(chi2map.T, 'matrix us ($1*%s+%s):3:2  w l palette , '%(v_step, v_lo),
               chi2map.min(axis=1), rv, e_rv, 'us 2:1:3 w xe pt 7, ',
               mm, 'us ($0*%s+%s):1 w l'%(v_step, v_lo),' lc 7 lw 2,',
               [mm.min()], [RV], [e_RV], 'us 2:1:3 w xe pt 7 lc 7 ')
         gplot.unset("multiplot")

      if 0:
         plot_input(*self.args)
         pause()

      # Spline interpolate each chi2map
      self.ll = ll = x[ind]
      rchi = self.args[-1]
      self.cs = cs = [spl.ucbspl_fit(vgrid, chi_o, K=len(chi_o)) for chi_o in (chi2map[ind]-np.array(self.SSRv)[:,np.newaxis])/rchi[ind][:,np.newaxis]**2]

      self.vv = vv = np.arange(vgrid[30], vgrid[-30], 0.01)

      if 0:
      # discrete brute force minisation
         lnLvb = []
         bb = np.arange(-0.8, 0.8, 0.005)
         for bbi in bb:
            lnLv = []
            for x_o,cs_o in zip(ll, cs):
               lnLv += [cs_o(vv+(x_o-xc)*bbi)]
            lnLvb += [lnLv]
         lnLvb = np.array(lnLvb)
         lnLvbint = lnLvb.sum(axis=1).T
         ii = np.unravel_index(lnLvbint.argmin(), lnLvbint.shape)
         crx_i = bb[ii[1]] * 1000   # km to m/s
         zz = vv[ii[0]] * 1000

      if 1:
      # paraboloid
         #imin,jmin = np.unravel_index(np.argmin(lnLvbint), lnLvbint.shape)
         #sl = np.s_[imin-10:imin+11, jmin-10:jmin+11]
         #X = np.tile(vv, bb.shape+(1,)).T[sl].ravel()
         #Y = np.tile(bb, vv.shape+(1,))[sl].ravel()
         #z = lnLvbint[sl].ravel()  # ln L = chi^2/2

         #cov = paraboloid.covmat_fit([X,Y], z, N=self.No[ind].sum())
         #cov =bbi paraboloid.covmat_fit([X,Y], z)

         # setup a grid around the start guess mlRv and mlCRX=0
         par = map(np.ravel, np.meshgrid(np.arange(self.mlRV/1000.-1, self.mlRV/1000.+1,0.2), np.arange(-500,500,50)/1000.))
         # par = map(np.ravel, np.meshgrid(np.arange(self.mlRV/1000.-0.1, self.mlRV/1000.+.1,0.02), np.arange(-500,500,50)/1000.))

         z = [np.sum([cs_o(vvi+(x_o-xc)*bbi) for x_o,cs_o in zip(ll, cs)]) for vvi,bbi in zip(*par)]

         # gplot.splot(par, z, ' palette')
         #pause()
         cov = paraboloid.covmat_fit(zip(*par), z, N=len(cs))

         zmod = cov.p(zip(*par))
         #pause("\n", cov.Va*1000**2, cov.Xc*1000)
         #cov = paraboloid.covmat_fit(par, z)
         #cov = paraboloid.covmat_fit(par, z, N=self.No[ind].sum())

         v, b = cov.Xc
         e_v, e_b = cov.e_a
         #gplot.splot(par, z, ' palette,', v,b, cov.min); pause("v", v, "b", b)
         #gplot.splot(par, z, zmod, ' palette, "" us 1:2:4 palette pt 6 , "" us 1:2:(($3-$4)*10),', v,b, cov.min); pause("v", v, "b", b)

         if -1<v<1 and -0.5<b<0.5 and np.isfinite(cov.e_a).all():
            # refit closer to minimum
            par = map(np.ravel, np.meshgrid(np.arange(v-3*e_v, v+3*e_v, e_v/2), np.arange(b-3*e_b, b+3*e_b, e_b/2)))
            z = [np.sum([cs_o(vvi+(x_o-xc)*bbi) for x_o,cs_o in zip(ll, cs)]) for vvi,bbi in zip(*par)]
            zmod = cov.p(zip(*par))
            cov = paraboloid.covmat_fit(zip(*par), z, N=len(cs))

         self.crx = cov.Xc[1] * 1000
         self.e_crx = cov.e_a[1] * 1000
         #print self.crx, self.e_crx

      if 0:
      # scipy optimisation
         def lnL(a):
            return np.sum([cs_o(a[0]+a[1]*(x_o-xc)) for x_o,cs_o in zip(ll, cs)])

         from scipy.optimize import fmin
         pp = fmin(lnL, [0,0], disp=False)
         #print pp * 1000
         self.zz = pp[0] * 1000
         self.crx = pp[1] * 1000
         self.e_crx = None # e_crx[i]


      #SSRstat(self.vgrid, self.mCCF, dk=1, plot=0)[0:2]
      #ii = np.argmin(self.mCCF); ss = np.s_[ii-3:ii+3]
      #uu = paraboloid.covmat_fit([self.vgrid[ss]], self.mCCF[ss])
      #uu.Xc, uu.e_a

      return self.crx, self.e_crx
      #print crx_i, crx[i], lnL([vv[ii[0]],bb[ii[1]]]), lnLvbint.min()
      #gplot(x,rv[i],e_rv[i], " w e,", x,((x-xc)*crx_i/1000.)*1000+zz, ',', x,(x-xc)* pval[1]+pval[0])

   # plot best fits:
   def plot_fit(self):
         RV, e_RV, rv, e_rv, orders, keytitle, rchi = self.args
         x, xc, ind = self.xarg
         vv = self.vv
         lnL0 = np.array([cs_o(vv) for cs_o in self.cs])
         gplot.palette('defined (0 "blue", 1 "green", 2 "red")')
         gplot.multiplot("layout 1,2")
         gplot.xlabel('"v [m/s]"; set ylabel "chi^2 / max(chi^2)"; set cblabel "order"')
         nlnL0 = lnL0.T /lnL0.max(axis=1)
         mm = nlnL0.sum(axis=1)

         #pause()
         # left panel: the lnL parabolas in RV[km/s]
         gplot(nlnL0, ' matrix us ($1*%s+%s):3:2 w l palette,'%(vv[1]-vv[0], vv[0]),
               nlnL0.min(axis=0), rv, e_rv, 'us 2:1:3 w xe pt 6 lc 7, ',
               mm/mm.max(), 'us ($0*%s+%s):1 w l'%(vv[1]-vv[0], vv[0]),' lc 7 lw 2,',
               #[(mm/mm.max()).min()], [RV[i]/1000], [e_RVc[i]/1000], 'us 2:1:3 w xe pt 7 lc 7 '
               )

         # right panel: the map ln(lam) vs RV [m/s]
         gplot.xlabel('"ln(wavelength)"; set ylabel "velocity v [m/s]"; set cblabel "chi^2 / max(chi^2)"')
         gplot.yrange('[%s:%s]'% ((rv*1000-e_rv*1000).min(), 1000*(rv+e_rv).max()))
         # lnL as maps ,width  shoul be addjust
         # often some maps are not drawn
         gplot("a='0 "+" ".join(map(str,self.ll))+"'"+', for [i=2:%s+1]'%len(lnL0), vv*1000, lnL0/lnL0.max(axis=1)[:,np.newaxis], ' us (word(a,i)+0):1:(0.002):(%s):i with boxxy palette fs solid t "",'%((vv[1]-vv[0])*1000/2),
               self.ll, rv*1000, e_rv*1000, " w e pt 7 lc 7 t 'RV_o',", #x,(x-xc)* self.crx+pval[0], 'w l lc 1,',
               x, ((x-xc)*self.crx/1000.)*1000+self.zz, 'w l lc 3 t "%s"' % self.crx)
         gplot.unset("multiplot")
         gplot.yrange('[*:*]')

         #pause()
         if 0:
            gplot.xlabel('"v [m/s]"; set ylabel "chi^2/max(chi^2)"; set cblabel "order"')
            #gplot(chi2map, ' matrix us ($1*%s+%s):3:2 w l palette'%(v_step, v_lo))
            gplot(chi2map.T /chi2map.max(axis=1), ' matrix us ($1*%s+%s):3:2 w l palette'%(v_step, v_lo))

            gplot(lnLvbint, 'matrix us 1:3:2  w l palette,', uu.sum(axis=0), "us ($0-30):1")
            #gplot(chi2map, 'matrix us 1:2:3 w l palette')
            ij=19; gplot(lnLvb[ij], 'matrix us 1:2:3 w image,', np.argmin(lnLvb[ij],axis=1))
            gplot.splot(np.log(np.array(lnLvbint).T-lnLvbint.min()), 'matrix us 1:2:3 w l palette')
      #plot_fit()


def fromfile(name, *args):
   chi2map =  np.array(pyfits.getdata(name), np.float)
   hdr = pyfits.getheader(name)
   vrange = hdr['CRVAL1'], hdr['CDELT1'] # -15, 0.1
   return Chi2Map(chi2map, vrange, *args, name=name)


if __name__ == '__main__':
   pass
