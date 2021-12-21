#! /usr/bin/python
from __future__ import division, print_function

import sys
import argparse

import numpy as np
from gplot import *
from pause import *
from wstat import nanwsem, wmean, mlrms, wstd
try:
   import gls
except:
   print('Cannot import gls')

try:
    import astropy.io.fits as pyfits
except:
    print('Cannot import astropy.io.fits')

import chi2map

__author__ = 'Mathias Zechmeister'
__version__ = '2021-03-31'

description = '''
SERVAL - SpEctrum Radial Velocity AnaLyser (%s)
     by %s
''' % (__version__, __author__)


# fix of keepdim capability of genfromtxt
genfromtxt2d = lambda *x,**y: np.atleast_2d(np.genfromtxt(*x,**y))


class srv:
   '''
   Analysis of SERVAL products.
   
   '''

   def __init__(self, obj, fibsuf='', oidx=None, safemode=False, pdf=False, plotrvo=True, cen=False):
      '''
      Load all SERVAL products into an object.

      obj : Folder with the SERVAL products.
      '''
      self.dir = obj
      self.tag = tag = obj.rstrip('/').split('/')[-1]
      self.pre = pre = obj+'/'+tag
      self.rms = np.nan

      print(pre+'.rvc'+fibsuf+'.dat')
      self.allrv = genfromtxt2d(pre+'.rvo'+fibsuf+'.dat')
      sbjd = np.genfromtxt(pre+'.rvo'+fibsuf+'.dat', dtype=('|S33'), usecols=[0])   # as string
      self.snr = genfromtxt2d(pre+'.snr'+fibsuf+'.dat')
      self.dlw = genfromtxt2d(pre+'.dlw'+fibsuf+'.dat')

      try:
         self.allerr = genfromtxt2d(pre+'.e_rvo'+fibsuf+'.dat')
         self.e_dlw = genfromtxt2d(pre+'.e_dlw'+fibsuf+'.dat')
      except:
         self.allerr = genfromtxt2d(pre+'.rvo'+fibsuf+'.daterr')
         print('No e_dlw. Consider updating and/or rerunning serval')
         self.e_dlw = 0 * self.dlw
      self.rchi = genfromtxt2d(pre+'.chi'+fibsuf+'.dat')
      try:
         self.halpha = genfromtxt2d(pre+'.halpha.dat').T
      except:
         pass
      try:
         self.tpre = genfromtxt2d(pre+'.pre'+fibsuf+'.dat')
      except:
         print('warning: %s not found' % pre+'.pre'+fibsuf+'.dat')
      self.dLW, self.e_dLW = self.dlw.T[[1,2]]

      self.brv = genfromtxt2d(pre+'.brv'+fibsuf+'.dat')
      # info includes also flagged files; exclude them based on unpairable bjd
      # (due to different formatting use bjd from brv.dat not info.csv.)
      bjd = np.atleast_1d(np.genfromtxt(pre+'.brv.dat', usecols=[0]))

      nn = [n for (n,t) in enumerate(bjd) if t in self.allrv[:,0]]
      try:
          self.info = np.genfromtxt(pre+'.info.csv', dtype=None, delimiter=';', encoding=None)[nn]
          self.info.dtype.names = ('timeid', 'bjd', 'berv', 'sn55', 'obj', 'exptime', 'ccf.mask', 'flag', 'airmass', 'ra', 'de', 'sunalt', 'moonsep', 'moonphase')
          self.sunalt = self.info['sunalt']
          self.moonsep = self.info['moonsep']
          self.moonphase = self.info['moonphase']
          self.flag = self.info['flag']
          self.info = self.info['timeid']
      except:
          self.info = np.atleast_1d(np.genfromtxt(pre+'.info.cvs', dtype=str, usecols=[0], delimiter=';'))[nn]
          print('suffix cvs was renamed to csv. Update by rerunning serval.')

      if not self.info.ndim:
         self.info = self.info[np.newaxis]
      info = " ".join(self.info)
      self.inst = ''
      if '-vis.fits' in info: self.inst = 'CARM_VIS'
      if '-nir.fits' in info: self.inst = 'CARM_NIR'
      if '_e2ds' in info: self.inst = 'HARPS'

      gplot.key_tit('"" noenh')
      self.keytitle = self.tag
      if self.inst:
         self.keytitle += ' (' + self.inst.replace('_', ' ') + ')'

      self.tcrx = genfromtxt2d(pre+'.crx.dat').T

      self.N = len(self.allrv)
      #if self.N == 1:
         #return   # just one line, e.g. drift

      self.tsrv = genfromtxt2d(pre+'.srv.dat').T
      self.trvc = self.bjd, RVc_old, e_RVc_old, RVd, e_RVd, RV_old, e_RV_old, BRV, RVsa \
                = genfromtxt2d(pre+'.rvc'+fibsuf+'.dat').T
      self.drs = genfromtxt2d(pre+'.drs.dat')
      try:
         self.tmlc = genfromtxt2d(pre+'.mlc'+fibsuf+'.dat')
      except:
         pass

      with np.errstate(invalid='ignore'):
         orders, = np.where(np.sum(self.allerr[:,5:]>0, 0))   # orders with all zero error values
      self.orders = orders
      if oidx is not None:
         omiss = set(oidx) - set(orders)
         if omiss: pause('WARNING: orders', omiss,'missing')
         else: orders = np.array(sorted(set(orders) & set(oidx)))

      omap = orders[:,np.newaxis]
      self.rv, self.e_rv = self.allrv[:,5+orders], self.allerr[:,5+orders]
      self.RV, self.e_RV = nanwsem(self.rv, e=self.e_rv, axis=1)
      self.has_d = ~np.isnan(RVd) * 1   # has a drift value, HARPS has no e_RVd

      # drift corrected
      self.rvc = self.rv - np.nan_to_num(RVd[:,np.newaxis]) - np.nan_to_num(RVsa[:,np.newaxis])
      self.RVc = self.RV - np.nan_to_num(RVd) - np.nan_to_num(RVsa)
      if cen:
          off = np.nanmedian(self.RVc)
          self.RVc -= off
          self.rvc -= off
      self.e_RVc = np.sqrt(self.e_RV**2 + np.nan_to_num(e_RVd)**2)

   def plot_dlw(self, *args, **kwargs):
      '''Show RVs over order for each observation.'''
      bjd, dLW, e_dLW = self.bjd, self.dLW, self.e_dLW
      arg = ''
      if not self.has_d.all():
         arg += 'us 1:2:3 w e pt 6 lt 7 t "dLW no drift"'
      if self.has_d.any():
         if arg: arg += ', "" '
         arg += 'us 1:2:($3/$4):8 w e pt 7 palette t "dLW"'
      hypertext = ', "" us 1:2:(sprintf("No: %d\\nID: %s\\nBJD: %f\\ndLW: %f +/- %f\\nflag: %d\\nsunalt: %.2f deg\\nmoonsep: %.2f deg\\nmoonphase: %.2f deg",$0+1, stringcolumn(5),$1, $2, $3, $6, $7, $8, $9)):8 w labels hypertext point pt 7 palette t ""'
      arg += hypertext
      arg += ', "" us 1:2:7 palette pt 6 ps 1.2'

      gplot.key('right Right top tit "%s" right'%(self.keytitle))
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"dLW [1000 (m/s)^2]"').cblabel('"     Sun altitude [deg]    Moon separation [deg]"')
      gplot.palette('defined (-30 "black", -5 "yellow", -4.9 "white", 4.9 "white", 5 "yellow", 30 "black")')\
            .cbrange('[-30:30]')

      args += (",", self.brv.T, ' us ($1-2450000):(-($10-$9)) lc 2 pt 12 t "-dBERV [m/s]"')
      gplot(bjd-2450000, dLW, e_dLW, self.has_d, self.info, self.flag, np.nan_to_num(self.sunalt), np.nan_to_num(self.moonsep), np.nan_to_num(self.moonphase), arg, *args, **kwargs)
      if len(args)==3: pause('dLW ', self.tag)

   def plot_halpha(self):
      '''Show Halpha time series.'''
      bjd, halpha, e_halpha = self.halpha[0:3]
      arg = 'us 1:2:3 w e pt 6 lt 7 t "no drift"'
      hypertext = ', "" us 1:2:(sprintf("No: %d\\nID: %s\\nBJD: %f\\nHalpha: %f +/- %f",$0+1, stringcolumn(5),$1, $2, $3)) w labels hypertext point pt 0  lt 1 t "",'
      arg += hypertext

      gplot.key('tit "%s"'%(self.keytitle))
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"Halpha index"')

      gplot(bjd-2450000, halpha, e_halpha, e_halpha, self.info, arg)
      pause('halpha ', self.tag)

   def drsrv(self):
      '''Show RVs over order for each observation.'''
      drsbjd, drsRVc, drse_RVc = self.drs.T[0:3]
      drsRVc = drsRVc - np.nan_to_num(self.trvc[8])
      bjd, RVc, e_RVc = self.bjd, self.RVc, self.e_RVc
      if 1:
         gplot.key('tit "%s"'%self.keytitle)
         gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]"')
         #gplot(';set ytics nomirr; set y2tics;', flush='')
         gplot(drsbjd-2450000, drsRVc, drse_RVc, 'us 1:2:3 w e pt 7 lt 1 t "DRS (rms = %1.3g m/s)"'% mlrms(drsRVc, e=drse_RVc)[0])
#         ogplot(bjd-2450000, RVc-np.median(RVc)+np.median(drsRVc), e_RVc, ' us 1:2:3 w e pt 7 lt 1 t "RVc-med(RVc)+med(DRS)"')
         ogplot(bjd-2450000, RVc-np.median(RVc)+np.nanmedian(drsRVc), e_RVc, ' us 1:2:3 w e pt 5 lt 3 t "\\nRVc (rms = %1.3g m/s)\\nrms(diff) = %.2f m/s"'%(self.mlrms[0], mlrms(drsRVc - RVc, e_RVc)[0]))
         pause('drsrv ', self.tag) # , ""  us 1:2:($3/$4) w e pt 7 lt 1 

   def mlc(self):
      '''Show RVs over order for each observation.'''
      _, mlRVc, e_mlRVc = self.tmlc.T[0:3]
      bjd, RVc, e_RVc = self.bjd, self.RVc, self.e_RVc
      if 1:
         gplot.key('tit "%s"'%self.keytitle)
         gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]"')
         gplot(bjd-2450000, RVc, e_RVc, ' us 1:2:3 w e pt 5 lc 1 t "RVc (rms = %1.3g m/s)'%self.mlrms[0])
         ogplot(bjd-2450000,self.tmlc.T[1:], 'us 1:2:3 w e pt 7 lc 3 t "\\nmlRVc (rms = %1.3g m/s)\\nrms(diff)=%.2f m/s"'% (mlrms(mlRVc, e=e_mlRVc)[0], mlrms(mlRVc - RVc, e_RVc)[0]))
         pause('rv ', self.tag) # , ""  us 1:2:($3/$4) w e pt 7 lt 1 

   def plot_chi2map(self, maxnorm=True, suf='_A_chi2map.fits'):
      n = 0
      #for n in range(self.N):
      while 0 <= n < self.N:
         name = self.dir +'/res/' + self.info[n][:-5] + suf

         chimap = chi2map.fromfile(name, self.RV[n]/1000., self.e_RV[n]/1000., self.rv[n]/1000., self.e_rv[n]/1000., self.orders, self.keytitle, self.rchi[n])

         chimap.plot()
         print('mlRV', chimap.mlRV, chimap.e_mlRV)
         print('RV  ', self.RV[n], self.e_RV[n])

         nn = pause('%i/%i %s %s'% (n+1, self.N, self.bjd[n], self.info[n]))
         try:
            n = int(nn)
         except:
            if nn in ('-', "D", "B"):
               n -= 1
            elif nn in ('^'):
               n = 0
            elif nn in ('$'):
               n = self.N - 1
            else:
               n += 1

   def plot_mlcrx(self, maxnorm=True, suf='_A_chi2map.fits'):
      n = 0
      #for n in range(self.N):
      while 0 <= n < self.N:
         name = self.dir +'/res/' + self.info[n][:-5] + suf
         chimap = chi2map.fromfile(name, self.RV[n]/1000., self.e_RV[n]/1000., self.rv[n]/1000., self.e_rv[n]/1000., self.orders, self.keytitle, self.rchi[n])

         x = self.tcrx[6:].T[n]
         ind = self.orders
         xc = np.mean(x[ind])
         crxml, e_crxml = chimap.mlcrx(x, xc, ind)
         print('crxml', crxml)
         print('crx  ', self.tcrx[1,n], self.tcrx[2,n])
         chimap.plot_fit()

         nn = pause('%i/%i %s %s'% (n+1, self.N, self.bjd[n], self.info[n]))
         try:
            n = int(nn)
         except:
            if nn in ('-', "D", "B"):
               n -= 1
            elif nn in ('^'):
               n = 0
            elif nn in ('$'):
               n = self.N - 1
            else:
               n += 1
      #pause()

   def gls(self, fend=1.):
      '''
      Period analysis of RV data.
      '''
      gg = gls.Gls((self.bjd -2450000, self.RVc, self.e_RVc), plot=1, fend=fend)
      gg.info()
      gg.plot(block=True)
      #pause()

   def kcita(self):
      import csv
      f = open('/home/astro115/carmenes/data/carmencita.081.csv', 'r')
      reader = csv.DictReader(f)
      x = [row for i, row in enumerate(reader) if row["#Karmn\t\t"]== self.tag+'\t']
      if x:
         x = x[0]
         print('#Karm', x["#Karmn\t\t"], ' | GJ'+x[" GJ\t\t"], ' | '+x[" Name\t\t\t"])
         print('vsini [km/s]:', x[" vsini_kms-1\t"], "(%s)"%x[" Ref28\t"])
         print('Prot [km/s]:', x[" P_d \t\t"], x[" eP_d\t\t"], "(%s)"%x[" Ref29\t"])
         print('SpT :', x[" SpT\t\t"])

   def targ(self):
     try:
      try:
          with open(self.pre+'.targ.csv') as f:
              self.line = f.read()
      except:
          with open(self.pre+'.targ.cvs') as f:
              self.line = f.read()
      line = self.line.split(';')        # ['gj699', "NAME Barnard's star", ' 17 57 ...]
      print(line[0],"; ", line[1])
      line = " ".join(line[2:]).split()
      self.ra = tuple(map(float,line[0:3]))  # rammss = (14.,29.,42.94)
      self.de = tuple(map(float,line[3:6]))  # demmss = (-62.,40.,46.16)
      self.pmra = float(line[6].replace("~","0."))             # pma = -3775.75
      self.pmde = float(line[7].replace("~","0."))             # pmd = 765.54
      self.plx = float(line[11].replace("~","nan"))
      self.rvabs = float(line[16].replace("~","nan"))
     except:
        pass


   def stat(self):
      '''
      Periodanalysis of RV data.
      '''
      self.mlrms = mlrms(self.RVc, e=self.e_RVc)
      print('N:     ', self.N)
      print('T [d]:     ', self.bjd.max() - self.bjd.min())
      print('wrms_RVc [m/s]:   %.2f\njitter [m/s]: %.2f' % self.mlrms)
      #pause()

   def plotrvno(self, rvobkg=False):
      '''Show RVs over order for each observation.'''
      bjd, RVc, e_RVc = self.bjd, self.RVc, self.e_RVc
      crx, e_crx = self.tcrx[1:3]
      lam_o = np.exp(self.tcrx[6:].T)
      lnlv = np.log(self.tcrx[5])
      print("use '$' to toggle between orders and wavelength")
      n = 0
      while 0 <= n < self.N:
         gplot.key('tit "%s\\n%s %s %s"'%(self.tag, n+1, bjd[n], self.info[n]))
         if 1:
            gplot.multiplot('layout 2,1')
            # top panel
            gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]"')
            gplot.unset('log')
            arg = ''
            if not self.has_d.all():
               arg += 'us 1:2:3 w e pt 6 lt 7 t "RV no drift",'
            if self.has_d.any():
               if arg: arg += ' "" '
               arg += 'us 1:($2/$4):3 w e pt 7 lt 7 t "RVc",'
            hypertext = '"" us 1:2:(sprintf("No: %d\\nID: %s\\nBJD: %f\\nRV: %f +/- %f",$0+1, stringcolumn(5),$1, $2, $3)) w labels hypertext point pt 0  lt 1 t "",'

            gplot(bjd-2450000, RVc, e_RVc, self.has_d, self.info, arg,
                  hypertext,
                  bjd[n]-2450000, RVc[n], e_RVc[n], 'us 1:2:3 w e pt 7 t "n=%s"'%(n+1))

         # bottom panel
         gplot.xlabel('"wavelength"').ylabel('"RV [m/s]"')
         gplot.autoscale('noextend')
         gplot.put('i=2; bind "$" "i = i%2+1; xlab=i==1?\\"order\\":\\"wavelength\\"; set xlabel xlab; set xra [*:*]; print i; if (i==1) {unset log} else {set log x} ; repl"')
         RVmod = crx[n]*(np.log(lam_o[n,self.orders])-lnlv[n])+RVc[n]
         RVlow = (crx[n]-e_crx[n])*(np.log(lam_o[n,self.orders])-lnlv[n])+RVc[n]
         RVupp = (crx[n]+e_crx[n])*(np.log(lam_o[n,self.orders])-lnlv[n])+RVc[n]
         hypertext = ', "" us i:3:(sprintf("No: %d\\nID: %s\\nBJD: %f\\nRV: %f +/- %f\\n' % (n+1,self.info[n],bjd[n], RVc[n], e_RVc[n])+ 'o: %d\\nrv[o]: %f +/- %f", $1, $3, $4)) w labels hypertext point pt 0 lt 1 t ""'

         plot = gplot
         if rvobkg:
             # plot all rvos as a background
             o = np.tile(self.orders, (self.rvc.shape[0],1))
             gplot-(o.ravel(), lam_o[:,self.orders].ravel(), self.rvc.ravel(), ' us i:3 pt 7 lc rgb "#99cccccc" t "all"')
             plot = gplot.oplot

         plot(self.orders, lam_o[n,self.orders], self.rvc[n], self.e_rv[n], RVmod, RVlow, RVupp, ' us i:3:4 w e lc 1 pt 7 t "n = %i"'%(n+1)+hypertext+', "" us i:5 w l lt 2 t "CRX = %g +/- %g m/s/Np", "" us i:6:7 w  filledcurves lt 2 fs transparent solid 0.20 t "",  %s lt 3 t "%g +/- %g m/s", "+" us 1:(%s):(%s) w filledcurves lt 3 fs transparent solid 0.20 t ""' %(crx[n], e_crx[n], RVc[n], RVc[n], e_RVc[n], RVc[n]-e_RVc[n], RVc[n]+e_RVc[n]))

         gplot.unset('multiplot')
         nn = pause('%i/%i %s %s'% (n+1, self.N, bjd[n], self.info[n]))
         try:
            n += int(nn)
         except:
            if nn in ('-', "D", "B"):
               n -= 1
            elif nn in ('^'):
               n = 0
            elif nn in ('$'):
               n = self.N - 1
            else:
               n += 1

   def plotrvo(self):
      bjd, RVc, e_RVc, RVd, e_RVd, RV, e_RV, BRV, RVsa = self.trvc
      allrv = self.allrv
      bjd, RV, e_RV, rv, e_rv = allrv[:,0], allrv[:,1], allrv[:,2], allrv[:,5:], self.allerr[:,5:]
      bjdmap = np.tile(bjd[:,np.newaxis], rv.shape[1])
      omap = np.tile(np.arange(len(rv.T)), rv.shape[0]).T

      # Drift and sa yet no applied to rvo
      # gplot(bjd, RV, e_RV, 'us 1:2:3 w e pt 7',)
      # gplot(bjdmap.ravel(), rv.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 palette, "'+obj+'/'+obj+'.rvo.dat'+'" us 1:2:3 w e pt 7 lt 7')
      rvc = rv - (np.nan_to_num(RVd) + np.nan_to_num(RVsa))[:,np.newaxis]
      gplot.palette('define (0 "blue", 1 "green", 2 "red")').key('tit "%s"'% self.keytitle)
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]')
      # not drift corrected
      #gplot(bjdmap.ravel(), rv.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 palette,', bjd, RV, e_RV, 'us 1:2:3 w e pt 7 lt 7')
      # drift corrected
      arg = ''
      hypertext = ''
      if not self.has_d.all():
         arg += 'us 1:2:3 w e pt 6 lt 7 t "RV no drift"'
      if self.has_d.any():
         if arg: arg += ', "" '
         arg += 'us 1:2:($3/$4) w e pt 7 lt 7 t "RVc"'
      if 1:
         hypertext = ' "" us 1:2:(sprintf("o: %d\\nBJD: %f\\nRV: %f", $4, $1, $2)) w labels hypertext point pt 0 lt 1 t "",'
         arg += ', "" us 1:2:(sprintf("No: %d\\nID: %s\\nBJD: %f\\nRV: %f+/-%f",$0+1, stringcolumn(5),$1, $2, $3)) w labels hypertext point pt 0  lt 1 t ""'

      omax = rv.shape[1]
      o = omap.ravel()[np.isfinite(rvc.ravel())].min()
      while 0 <= o < omax:
         gplot-(bjdmap.ravel()-2450000, rvc.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 ps 0.5 palette t "RV_o",'+hypertext, bjd-2450000, RVc, e_RVc, self.has_d, self.info, arg) # 'us 1:2:3 w e pt 6 lt 7, ""  us 1:2:($3/$4) w e pt 7 lt 7')
         gplot+(bjdmap[:,o]-2450000, rvc[:,o], e_rv[:,o], 'us 1:2:3 w e pt 12 lc 9 t "%s"'%o)
         oo = pause('%i/%i'% (o, omax))
         try:
            o += int(oo)
         except:
            if oo in ('-', "D", "B"):
               o -= 1
            elif oo in ('^',):
               o = 0
            elif oo in ('$',):
               o = self.N - 1
            else:
               o += 1

      pause('rvo ', self.tag)

   def plot_dlwo(self):
      bjd, RVc, e_RVc, RVd, e_RVd, RV, e_RV, BRV, RVsa = self.trvc
      allrv = self.allrv
      bjd, RV, e_RV, rv, e_rv = allrv[:,0], allrv[:,1], allrv[:,2], allrv[:,5:], self.allerr[:,5:]
      dlw = self
      (bjd, dLW, e_dLW), dlw, e_rv = self.dlw.T[[0,1,2]], self.dlw[:,3:], self.allerr[:,5:]
      bjdmap = np.tile(bjd[:,np.newaxis], rv.shape[1])
      omap = np.tile(np.arange(len(self.dlw.T)), dlw.shape[0]).T

      rvc = rv - (np.nan_to_num(RVd) + np.nan_to_num(RVsa))[:,np.newaxis]
      gplot.palette('define (0 "blue", 1 "green", 2 "red")').key('tit "%s"'% self.keytitle)
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"dLW [1000 (m/s)^2]')
      # not drift corrected
      #gplot(bjdmap.ravel(), rv.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 palette,', bjd, RV, e_RV, 'us 1:2:3 w e pt 7 lt 7')
      # drift corrected
      arg = ''
      if not self.has_d.all():
         arg += 'us 1:2:3 w e pt 6 lt 7 t "dLWo (no drift)"'
      if self.has_d.any():
         if arg: arg += ', "" '
         arg += 'us 1:2:($3/$4) w e pt 7 lt 7 t "dLWo"'
      gplot(bjdmap.ravel()-2450000, dlw.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 ps 0.5 palette t "RV_o",', bjd-2450000, dLW, e_dLW, self.has_d, arg) # 'us 1:2:3 w e pt 6 lt 7, ""  us 1:2:($3/$4) w e pt 7 lt 7')

      pause('dLWo ', self.tag)

   def plot_dlwno(self):
      '''Show RVs over order for each observation.'''
      (bjd, dLW, e_dLW), dlw, e_dlw = self.dlw.T[[0,1,2]], self.dlw[:,3:], self.e_dlw[:,3:]
      lam_o = np.exp(self.tcrx[6:].T)
      n = 0
      while 0 <= n < self.N:
         gplot.key('tit "%s\\n%s %s %s"'%(self.tag, n+1, bjd[n], self.info[n]))
         if 1:
            gplot.multiplot('layout 2,1')
            # bottom panel
            self.plot_dlw(",", [bjd[n]-2450000], [dLW[n]], 'us 1:2 lc "#ff0000" pt 4 ps 1. t "n=%s (%s)"'%(n+1, self.info[n]), flush=' \n')

         gplot.xlabel('"wavelength {/Symbol l} [A]"')\
              .ylabel('"dLW [1000(m/s)^2]"')
         gplot.put('i=2; bind "$" "i = i%2+1; xlab=i==1?\\"order\\":\\"wavelength\\"; set xlabel xlab; set xra [*:*]; print i; if (i==1) {unset log} else {set log x} ; repl"')
         hypertext = ', "" us i:3:(sprintf("No: %d\\nID: %s\\nBJD: %f\\ndLW: %f +/- %f\\n' % (n+1,self.info[n],bjd[n], dLW[n], e_dLW[n])+ 'o: %d\\ndlw[o]: %f +/- %f", $1, $3, $4)) w labels hypertext point pt 0 lt 1 t ""'
         gplot(self.orders, lam_o[n,self.orders], dlw[n,self.orders], e_dlw[n,self.orders], ' us i:3:4 w e pt 7 t "dLW_{%s,o} (%s)"' % (n+1, self.info[n]) + hypertext,
              ', %s t "dLW_{%s} = %5g +/- %5g", "+" us 1:(%s):(%s) w filledcurves lt 3 fs transparent solid 0.20 t ""' % (dLW[n], n+1, dLW[n], e_dLW[n], dLW[n]- e_dLW[n], dLW[n]+ e_dLW[n]))
         gplot.unset('multiplot')
         nn = pause('%i/%i %s %s'% (n+1,self.N, bjd[n], self.info[n]))
         try:
            n += int(nn)
         except:
            if nn in ('-', "D", "B"):
               n -= 1
            elif nn in ('^'):
               n = 0
            elif nn in ('z'):
               break
            elif nn in ('$'):
               n = self.N - 1
            else:
               n += 1

   def plotrv(self):
      '''Show RVs over order for each observation.'''
      bjd, RVc, e_RVc = self.bjd, self.RVc, self.e_RVc
      arg = ''
      if not self.has_d.all():
         arg += 'us 1:2:3 w e pt 6 lt 7 t "RV no drift"'
      if self.has_d.any():
         if arg: arg += ', "" '
         arg += 'us 1:($2/$4):3 w e pt 7 lt 7 t "RVc"'
      hypertext = ', "" us 1:2:(sprintf("No: %d\\nID: %s\\nBJD: %f\\nRV: %f +/- %f",$0+1, stringcolumn(5),$1, $2, $3)) w labels hypertext point pt 0  lt 1 t ""'
      arg += hypertext
      gplot.key('tit "%s (rms = %.3g m/s)"' % (self.keytitle, self.mlrms[0]))
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]"')
      gplot(bjd-2450000, RVc, e_RVc, self.has_d, self.info, arg)
      pause('rv ', self.tag)

   def postrv(self, postiter=1, fibsuf='', oidx=None, safemode=False, pdf=False, plotrvo=True):
      """
      """
      #pause()

      rv, e_rv = self.allrv[:,5+self.orders], self.allerr[:,5+self.orders]
      RV, e_RV = nanwsem(rv, e=e_rv, axis=1)

      allrv = self.allrv
      allerr = self.allerr

      bjd, RVc, e_RVc, RVd, e_RVd, RV, e_RV, BRV, RVsa = self.trvc

      bjd, RV, e_RV, rv, e_rv = allrv[:,0], allrv[:,1], allrv[:,2], allrv[:,5:], allerr[:,5:]
      bjdmap = np.tile(bjd[:,np.newaxis], rv.shape[1])
      omap = np.tile(np.arange(len(rv.T)), rv.shape[0]).T


      #gplot(bjd, RV, e_RV, 'us 1:2:3 w e pt 7',)
      #gplot(bjdmap.ravel(), rv.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 palette, "'+obj+'/'+obj+'.rvo.dat'+'" us 1:2:3 w e pt 7 lt 7')
      rvc = rv - (np.nan_to_num(RVd) + np.nan_to_num(RVsa))[:,np.newaxis]

      gplot.palette('define (0 "blue", 1 "green", 2 "red")').key('tit "%s"'% self.keytitle)
      gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]')
      # not drift corrected
      gplot(bjdmap.ravel()-2450000, rv.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 palette,', bjd-2450000, RV, e_RV, 'us 1:2:3 w e pt 7 lt 7')
      # drift corrected
      #gplot(bjdmap.ravel(), rvc.ravel(), e_rv.ravel(), omap.ravel(), 'us 1:2:4 w p pt 7 ps 0.5 palette,', bjd, RVc, e_RVc, 'us 1:2:3 w e pt 7 lt 7')

      pause(self.tag)

      # minimise (rv_io - rv_i - m_o)^2 / (e_rv_io^2 + e_rv_i^2 + e_rv_o^2)
      e_n = 0
      e_o = 0
      m_o = 0
      for e_n in [0,0.1,1,2,5,10,20,30,50,100,200]:
         e_r = np.sqrt(e_rv**2 + e_n**2 + e_o**2)
         RV, e_RV = nanwsem(rv-m_o, e=e_r, axis=1)
         RV = RV[:,np.newaxis]
         #stop()
         m_o = (wmean(rv-RV, 1/e_r**2, axis=0)*0) [np.newaxis,:]  # order offset
         r = rv - RV - m_o # residuals
         stop()
         e_n = np.sqrt((np.nansum((r**2-e_rv**2-e_o**2)/e_r**2, axis=1)/np.nansum(1/e_r**2, axis=1)).clip(min=0))[:,np.newaxis] # keepdims not in v1.6.1
         e_o = np.sqrt((np.nansum((r**2-e_rv**2-e_n**2)/e_r**2, axis=0)/np.nansum(1/e_r**2, axis=0)).clip(min=0))[np.newaxis,:]
         L = - 0.5* np.nansum(r**2/e_r**2) - 0.5* np.nansum(np.log(2.*np.pi*e_r**2))
         print(L, e_n.ravel())
         gplot.unset(' multiplot')
         gplot.datafile('missing "nan"')
         gplot.multiplot('layout 1,2;')
         gplot.xlabel('"order"').ylabel('"rvo - RV [m/s]"')
         gplot(rv - RV, ' matrix,', m_o.ravel(), 'us 0:1')
         gplot.xlabel('"obs No."').ylabel('"rvo - RV [m/s]"')
         gplot((rv - RV).T, ' matrix')
         gplot.unset('multiplot')
         pause()

      # ind, = where(np.isfinite(e_rv[i])) # do not use the failed and last order
      RV, e_RV = nanwsem(rv, e=e_rv, axis=1)
      RVc = RV - np.nan_to_num(RVd) - np.nan_to_num(RVsa)
      e_RVc = np.sqrt(e_RV**2 + np.nan_to_num(e_RVd)**2)
      pause()

   def prerv(self):
      '''Show RVs over order for each observation.'''
      prebjd, preRV, pree_RVc = self.tpre.T[0:3]
      preRVc = preRV - np.nan_to_num(self.trvc[8])
      bjd, RVc, e_RVc = self.bjd, self.RVc, self.e_RVc
      if 1:
         gplot.key('tit "%s"'%self.keytitle)
         gplot.xlabel('"BJD - 2 450 000"').ylabel('"RV [m/s]"')
         #gplot(';set ytics nomirr; set y2tics;', flush='')
         gplot(prebjd-2450000, preRVc, pree_RVc, 'us 1:2:3 w e pt 5 lt 3 t "preRVc (rms = %1.3g m/s)"'% mlrms(preRVc, e=pree_RVc)[0])
#         ogplot(bjd-2450000, RVc-np.median(RVc)+np.median(drsRVc), e_RVc, ' us 1:2:3 w e pt 7 lt 1 t "RVc-med(RVc)+med(DRS)"')
         ogplot(bjd-2450000, RVc-np.median(RVc)+np.median(preRVc), e_RVc, ' us 1:2:3 w e pt 7 lt 1 t "RVc (rms = %1.3g m/s)"'%self.mlrms[0])
         pause('rv ', self.tag) # , ""  us 1:2:($3/$4) w e pt 7 lt 1 

   def spaghetti(self):
         o = np.tile(self.orders, (self.rvc.shape[0],1))
         bjd = np.tile(self.bjd,  (self.rvc.shape[1],1)).T

         # color bjd with lc variable => same color for same floor(bjd)
         #gplot(ovis, vis.rvc.ravel(), bjdvis.ravel(), ' lc variable pt 7 ps 0.5,', onir+60, nir.rvc.ravel(), bjdnir.ravel(), ' lc variable pt 7 ps 0.5')
         gplot.key('tit "%s"' % self.keytitle)
         gplot.palette('defined')
         gplot.xlabel('"Order"').ylabel('"RV [m/s]"').cblabel('"BJD - 2 450 000"')
         gplot(o.ravel(), self.rvc.ravel(), bjd.ravel(), 'us 1:2:($3-2450000) palette')
         pause('spaghetti %s' % self.keytitle)

   def plot_vsini(self):
         gplot.ylabel('"rotation velocity vsini [km/s]"')
         gplot.xlabel('"order index o"')
         gplot('"%s" us (0.5):(med=$1):(e_vsini=$2) every :::0::0' % (self.pre+'.vsini.dat') )
         gplot('med+e_vsini w filledcurves y=med-e_vsini fc rgb "#e3f3ff" t "", med lc 3 t sprintf("median = %.2f +/- %.2f km/s", med, e_vsini)', ', "" us 1:2:3 every :::1::1  w e lc 1 pt 6 t "%s"' % self.tag)
         pause('vsini %s' % self.keytitle)

   def disp(self):
      orddisp = self.rv - self.RV[:,np.newaxis]
      #ok &= np.abs(orddisp-d_ordmean) <= 3*ordstd  # clip and update mask
      ordstd, d_ordmean = wstd(orddisp, self.e_rv, axis=0)
      #ordmean += d_ordmean                # update ordmean
      #orddisp
      gplot.reset()
      gplot.tit('"Order dispersion %s";'%self.keytitle)
      gplot.xlabel('"Order"').ylabel('"RV_{n,o} - RV_n [m/s]"')
      gplot(orddisp.T, ' matrix us (%s-1):3 t ""' % "".join(['$1==%s?%s:' % io for io in enumerate(self.orders)]))
      ogplot(self.orders, ordstd, ' w lp lt 3 t "", "" us 1:(-$2) w lp t ""')
      pause('disp')
      ##gplot('"'+filename,'" matrix every ::%i::%i us ($1-5):3' % (omin+5,omax+5))
      ##ogplot(allrv[:,5:71],' matrix every ::%i us 1:3' %omin)
      ## create ord, rv,e_rv, bb
      ##ore = [(o+omin,orddisp[n,o],e_rv[n,o],~ok[n,o]) for n,x in enumerate(orddisp[:,0]) for o,x in enumerate(orddisp[0,:])]
      ##ore = [(o,)+x for row in zip(orddisp,e_rv,~ok) for o,x in enumerate(zip(*row),omin)]
      #ore = [np.tile(orders,orddisp.shape).ravel(), orddisp.ravel(), e_rv.ravel(), ~ok.ravel()]
      #gplot(*ore + ['us 1:2:($3/30) w xe'], flush='')
      #if not ok.all(): ogplot('"" us 1:2:($3/30/$4) w xe', flush='') # mark 3-sigma outliners
      ##gplot(orddisp,' matrix  us ($1+%i):3' % omin, flush='')
      ##gplot('"'+filename,'" matrix every ::'+str(omin+5)+' us ($1-5):3')
      ##ogplot(ordmean,' us ($0+%i):1 w lp lt 3 pt 3  t "ord mean"' %omin, flush='')
      #ogplot(orders, ordstd,' w lp lt 3 pt 3  t "1 sigma"', flush='')
      #ogplot('"" us 1:(-$2) w lp lt 3 pt 3  t ""', flush='')
      #ogplot('"" us 1:($2*3) w lp lt 4 pt 3  t "3 sigma"', flush='')
      #ogplot('"" us 1:(-$2*3) w lp lt 4 pt 3  t ""', flush='')
      #ogplot('"" us ($1+0.25):(0):(sprintf("%.2f",$2)) w labels rotate t""', flush='')
      #ogplot(*ore+[ 'us 1:2:($3)  w point pal pt 6'])
      #if not safemode: pause('ord scatter')

   def ls(self, suf='_A_mod.fits'):
      '''Show last square fit form fits file.'''
      #prebjd, preRV, pree_RVc = self.tpre.T[0:3]
      #for n in range(self.N):
      gplot.key("Left left rev bottom title '%s'" % self.keytitle)
      '''
             if def_wlog: w2 = np.exp(w2)
            res = np.nan * f2
            res[pmin:pmax] = (f2[pmin:pmax]-f2mod[pmin:pmax]) / e2[pmin:pmax]  # normalised residuals
            b = str(stat['std'])
            #gplot_set('set key left Left rev samplen 2 tit "%s (o=%s, v=%.2fm/s)"'%(obj,o,rvo))
            #gplot_set('set ytics nomirr; set y2tics; set y2range [-5*%f:35*%f]; set bar 0.5'%(rchi,rchi))
            #gplot_set('i=1; bind "$" "i = i%2+1; xlab=i==1?\\"pixel\\":\\"wavelength\\"; set xlabel xlab; set xra [*:*]; print i; repl"')
            #gplot('[][][][-5:35]',x2, w2, f2, e2.clip(0.,f2.max()), 'us (column(i)):3:4 w errorli t "'+sp.timeid+' all"', flush='')
            ogplot(x2,w2, f2, ((b2==0)|(b2==flag.clip))*0.5, 1+4*(b2==flag.clip), 'us (column(i)):3:4:5 w p pt 7 lc var ps var t "'+sp.timeid+' telluric free"', flush='')
            ogplot(x2,w2, f2mod,(b2==0)*0.5, 'us (column(i)):3:4 w lp lt 3 pt 7 ps var t "Fmod"', flush='')
            ogplot(x2,w2, res, b2, "us (column(i)):3:4 w lp pt 7 ps 0.5 lc var axis x1y2 t 'residuals'", flush='')
            # legend with translation of bpmap, plot dummy using NaN function
            ogplot(", ".join(["NaN w p pt 7 ps 0.5 lc "+str(f) +" t '"+str(f)+" "+",".join(flag.translate(f))+"'" for f in np.unique(b2)]), flush='')

            ogplot("0 axis x1y2 lt 3 t'',"+b+" axis x1y2 lt 1,-"+b+" axis x1y2 lt 1 t ''", flush='')

            ogplot(x2,w2, ((b2&flag.atm)!=flag.atm)*40-5, 'us (column(i)):3 w filledcurve x2 fs transparent solid 0.5 noborder lc 9 axis x1y2 t "tellurics"', flush='')
            ogplot(x2,w2, ((b2&flag.sky)!=flag.sky)*40-5, 'us (column(i)):3 w filledcurve x2 fs transparent solid 0.5 noborder lc 6 axis x1y2 t "sky"')
            pause('large RV ' if abs(rvo/1000-rvguess+vref)>rvwarn else 'look ', o, ' rv = %.3f +/- %.3f m/s   rchi = %.2f' %(rvo, e_rv[i,o],rchi2[i,o]))
      # end loop over orders
      '''
      n = 0
      o = 40
      while True:
         name = self.dir +'/res/' + self.info[n][:-5] + suf
         hdu = pyfits.open(name)
         resmap = hdu['res'].data
         fmod = hdu['fmod'].data
         e_f = hdu['err'].data
         wave = hdu['wave'].data
         flux = fmod + resmap
         bpmap = 0 *wave
         rvo = 1
         if 1:
            gplot_set('set xlabel "wavelength [A]"; set ylabel "flux"')
            gplot_set('set key default left Left rev samplen 2 tit "%s (o=%s, v=%.2fm/s)"'%(self.tag,o,rvo))
            gplot_set('set ytics nomirr; set y2tics; set y2range [-5*%f:35*%f]; set bar 0.5'%(1,1))
            gplot_set('i=1; bind "$" "i = i%2+1; xlab=i==1?\\"pixel\\":\\"wavelength\\"; set xlabel xlab; set xra [*:*]; print i; repl"')
            x2 = wave[o]
            w2 = wave[o]
            f2 = flux[o]
            f2mod = fmod[o]
            e2 = e_f[o]
            b2= bpmap[o]
            res = resmap[o]/ f2
            gplot('[][][][-5:35]',x2, w2, f2, e2.clip(0.,f2.max()), 'us (column(i)):3:4 w errorli t "'+self.info[n]+' all"', flush='')
            #ogplot(x2,w2, f2, ((b2==0)|(b2==flag.clip))*0.5, 1+4*(b2==flag.clip), 'us (column(i)):3:4:5 w p pt 7 lc var ps var t "'+sp.timeid+' telluric free"', flush='')
            ogplot(x2,w2, f2mod,(b2==0)*0.5, 'us (column(i)):3:4 w lp lt 3 pt 7 ps var t "Fmod"', flush='')
            ogplot(x2,w2, res/e2, b2, "us (column(i)):3:4 w lp pt 7 ps 0.5 lc var axis x1y2 t 'residuals'")
            #pause()
            # legend with translation of bpmap, plot dummy using NaN function
            #ogplot(", ".join(["NaN w p pt 7 ps 0.5 lc "+str(f) +" t '"+str(f)+" "+",".join(flag.translate(f))+"'" for f in np.unique(b2)]), flush='')
            no = pause('%i/%i %s %s %s'% (n+1, self.N, o, self.bjd[n], self.info[n]))
            try:
               o = int(no)
            except:
               if no in ('-', "D", "B"):
                  n -= 1
               elif no in ('e'):
                  break
               elif no in ('^'):
                  n = 0
               elif no in ('$'):
                  n = self.N - 1
               elif no in (')'):
                  o += 1
               elif no in ('('):
                  o -= 1
               else:
                  n += 1
            n = np.clip(n, 0, self.N-1)
            o = np.clip(o, 0, resmap.shape[0] - 1)


   def xcorr(self, block=True):
      '''Correlation plots.'''
      import matplotlib
      if (matplotlib.get_backend() != "TkAgg"):
         matplotlib.use("TkAgg")

      import matplotlib.pylab as mlp

      # Drift and sa yet no applied to rvo
      #rvc = rv - (RVd + RVsa)[:,np.newaxis]
      bjd, RVc, e_RVc = self.trvc[0:3]
      bjd, RV, e_RV, crx, e_crx, dwid, e_dwid = self.tsrv

      col = mlp.cm.rainbow(mlp.Normalize()(bjd))
      datstyle = dict(fmt='r.', capsize=0, elinewidth=1., markeredgewidth=0.5, markeredgecolor='k')
      def plot_ecol(plt, x, y, e_y=None, e_x=None):
           # script for scatter plot with errorbars and time color-coded
           datstyle = dict(color=col, marker='.', edgecolor='k', linewidth=0.5, zorder=2)
           if e_y is not None:
              errstyle = dict(yerr=e_y, marker='', ls='', elinewidth=0.5)
              if matplotlib.__version__ < '2.' :
                 errstyle['capsize'] = 0.
                 datstyle['s'] = 8**2   # requires square size !?
              else:
                 errstyle['ecolor'] = col
              #_, _, (c,) = plt.errorbar(x, y, xerr=e_x, **errstyle)
              plt.errorbar(x, y, xerr=e_x, **errstyle)
              if matplotlib.__version__ < '2.':
                 c.set_color(col)
           plt.scatter(x, y, **datstyle)

      fig = mlp.figure()
      fig.subplots_adjust(hspace=0.10, wspace=0.08, right=0.98, top=0.98)

      # BJD-RV
      ax1 = fig.add_subplot(3, 2, 1)
      #ax.set_title("Normalized periodogram")
      mlp.setp(ax1.get_xticklabels(), visible=False)
      ax1.set_ylabel("RV [m/s]")
      plot_ecol(ax1, bjd-2450000, RVc, e_RVc)
#      ax1.errorbar(bjd-2450000, RVc, e_RVc, ecolor=col, color=col, **datstyle)

      # BJD-COL
      ax3 = fig.add_subplot(3, 2, 3, sharex=ax1)
      mlp.setp(ax3.get_xticklabels(), visible=False)
      ax3.set_ylabel("chromatic index")
      #ax3.errorbar(bjd-2450000, crx, e_crx, **datstyle)
      plot_ecol(ax3, bjd-2450000, crx, e_crx)

      # BJD-DLW
      ax5 = fig.add_subplot(3, 2, 5, sharex=ax1)
      ax5.set_xlabel("BJD - 2 450 000")
      ax5.set_ylabel("dLW")
      #ax5.errorbar(bjd-2450000, dwid, e_dwid, **datstyle)
      plot_ecol(ax5, bjd-2450000, dwid, e_dwid)

      # RV-COL
      ax4 = fig.add_subplot(3, 2, 4, sharey=ax3)
      mlp.setp(ax4.get_xticklabels(), visible=False)
      mlp.setp(ax4.get_yticklabels(), visible=False)

      #np.fit(RVc, crx)
      a, cov_a = np.polyfit(RVc, crx, 1, cov=True)
      e_a = np.sqrt(cov_a[0,0])
      liney = [RVc.min(), RVc.max()]
      linex = np.polyval(a, liney)
      ax4.plot(liney, linex, label='kappa: %.4g+/-%.4g'% (a[0], e_a), zorder=3, color='k')
      #ax4.errorbar(RVc, crx, e_crx, xerr=e_RVc, **datstyle)
      plot_ecol(ax4, RVc, crx, e_crx, e_x=e_RVc)
      ax4.legend(loc='upper right', frameon=False, framealpha=1, fontsize='small')

      # RV-DLW
      ax6 = fig.add_subplot(3, 2, 6, sharex=ax4, sharey=ax5)
      mlp.setp(ax6.get_yticklabels(), visible=False)
      ax6.set_xlabel("RV [m/s]")
      #ax6.errorbar(RVc, dwid, e_dwid, xerr=e_RVc, **datstyle)
      plot_ecol(ax6, RVc, dwid, e_dwid, e_x=e_RVc)

      for x in fig.axes:
         x.tick_params(direction='in', which='both', top=True, right=True)
         x.minorticks_on()

      if hasattr(mlp.get_current_fig_manager(), 'toolbar'):
         # check seems not needed when "TkAgg" is set
         mlp.get_current_fig_manager().toolbar.pan()
      if block:
         print("Close the plot to continue.")
      else:
         mlp.ion()

      mlp.show()

      pause(obj)

if __name__ == "__main__":
   '''
   Example:
   '''
   default = " (default: %(default)s)."
   epilog = """\
   usage example:
   %(prog)s tag dir_or_filelist -targ gj699 -snmin 10 -oset 40:
   """
   parser = argparse.ArgumentParser(description=description, epilog=epilog, add_help=False)
   argopt = parser.add_argument   # function short cut
   argopt('tags', nargs='*', help='Tag, output directory and file prefix')
   argopt('-cen', help='center RVs to zero median', action='store_true')
   argopt('-chi2map', help='plot the chi2map', action='store_true')
   argopt('-disp', help='plot order dispersion', action='store_true')
   argopt('-dlw', help='plot dLW', action='store_true')
   argopt('-dlwo', help='plot dLW_o colorcoded', action='store_true')
   argopt('-dlwno', help='plot dLW and the dLW_o for spectrum n in a lower panel', action='store_true')
   argopt('-drs', help='plot DRS RV vs RVc', action='store_true')
   argopt('-gls', help='GLS periodogram', action='store_true')
   argopt('-i', help='interactive task selection', action='store_true')
   argopt('-halpha', help='plot Halpha index', action='store_true')
   argopt('-ls', help='ls from fits file', action='store_true')
   argopt('-mlc', help='plot mlRVc vs RVc', action='store_true')
   argopt('-mlcrx', help='plot fitting of mlcrx', action='store_true')
   argopt('-pre', help='plot preRV vs RVc', action='store_true')
   argopt('-postrv', help='kappa sigma clip value', action='store_true')
   argopt('-rv', help='plot rv', action='store_true')
   argopt('-rvno', help='plot rv and the rvo for spectrum n in a lower panel', action='store_true')
   argopt('-rvnno', help='plot rv and the rvo for spectrum n in a lower panel with all rvo as background', action='store_true')
   argopt('-rvo', help='plot rvo colorcoded', action='store_true')
   argopt('-spaghetti', help='plot o-rvno colorcoded', action='store_true')
   argopt('-vsini', help='plot measured vsini for each order (see serval option -vsiniauto)', action='store_true')
   argopt('-x', help='cross plot'+default, action='store_true')
   argopt('-?', '-h', '-help', '--help',  help='show this help message and exit', action='help')

   args = parser.parse_args()

   for tag in args.tags:
      obj = srv(tag, plotrvo='plotrvo' in sys.argv, cen=args.cen)
      obj.targ()
      obj.stat()
      if not True in args.__dict__.values():
         obj.plotrv()
      elif args.i:
         while True:
            g = pause('next 0, q: quit')
            if g=='r': obj.plotrv()
            if g=='1': obj.drsrv()
            if g=='x': obj.xcorr()
            if g=='g': obj.gls()
            if g=='4': obj.plotrvno(rvobkg=args.rvnno)
            if g=='5': obj.plotrvo()
            if g=='6': obj.postrv()
      else:
         if args.rv:
            obj.plotrv()
         if args.dlw:
            obj.plot_dlw()
         if args.disp:
            obj.disp()
         if args.drs:
            obj.drsrv()
         if args.pre:
            obj.prerv()
         if args.spaghetti:
            obj.spaghetti()
         if args.ls:
            obj.ls()
         if args.halpha:
            obj.plot_halpha()
         if args.x:
            obj.xcorr()
         if args.gls:
            obj.gls()
         if args.mlc:
            obj.mlc()
         if args.rvno or args.rvnno:
            obj.plotrvno(rvobkg=args.rvnno)
         if args.chi2map:
            obj.plot_chi2map()
         if args.mlcrx:
            obj.plot_mlcrx()
         if args.rvo:
            obj.plotrvo()
         if args.vsini:
            obj.plot_vsini()
         if args.dlwo:
            obj.plot_dlwo()
         if args.dlwno:
            obj.plot_dlwno()
         if args.postrv:
            obj.postrv()


