#!/usr/bin/python
import tempfile
import os
import sys
import numpy as np
from subprocess import Popen, PIPE, STDOUT
import argparse

try:
   import pyfits
except:
   import astropy.io.fits as pyfits

path = os.path.dirname(__file__)
if path=='': path='.'

def bary_harps(file, ra=None, dec=None, epoch=2000, pma=0.0, pmd=0.0):
   ''' compute barycentric correction for HARPS spctrum using coordinates from
    input or fits header
    ra - rammss e.g. -142942.94 (for 14:29:42.94)
    dec - demmss e.g. -624046.16 (-62:40:46.16)
    pma - [mas/yr] proper motion in alpha (-3775.75)
    pmd - [mas/yr] proper motion in delta (765.54)
    epoch - epoch of coordinates
   '''
   head = pyfits.getheader(file)
   if ra is None:
      ra = head['HIERARCH ESO TEL TARG ALPHA']
      dec = head['HIERARCH ESO TEL TARG DELTA']
      epoch = head['HIERARCH ESO TEL TARG EPOCH']
      pma = head['HIERARCH ESO TEL TARG PMA']*1000
      pmd = head['HIERARCH ESO TEL TARG PMD']*1000

   exptime = head['EXPTIME'] * 2*head['HIERARCH ESO INS DET1 TMMEAN']
   dateobs = head['DATE-OBS']
   #bjd = head['HIERARCH ESO DRS BERV']
   #berv = head['HIERARCH ESO DRS BJD']

   x = str(ra).split('.')
   rammss = float(x[0][:-4]), float(x[0][-4:-2]), float(x[0][-2:]+'.'+x[1])
   # rammss = divmod(ra,10000); rammss = (rammss[0],) + divmod(rammss[1],100)
   x = str(dec).split('.') 
   demmss = float(x[0][:-4]), float(x[0][-4:-2]), float(x[0][-2:]+'.'+x[1])
   #print ra, dec, pma, pmd
   return bary(dateobs, rammss, demmss, 14, epoch, exptime, pma, pmd) #, bjd, berv

def bary(dateobs, rammss, demmss, inst, epoch=2000, exptime=0.0, pma=0.0, pmd=0.0, obsloc=None):
   ''' baricentric correction with BarCor
   rammss - (ra,mm,ss) or 'ra:mm:ss' or ra
   demmss - (de,mm,ss)

   Examples
   --------
   >>> bary.bary('2016-02-04T05:55:21', (17.0, 57.0, 48.4997994034), (4.0, 41.0, 36.111354228), 'CARM_VIS', epoch=2000, exptime=450.1, pma=-802.803, pmd=10362.542)
   array([2.45742275e+06, 1.90159600e+01])

   >>> bary.bary('2016-02-04T05:55:21', (17.0, 57.0, 48.4997994034), (4.0, 41.0, 36.111354228), 'PASS', epoch=2000, exptime=450.1, pma=-802.803, pmd=10362.542, obsloc={'lat':-2.5463, 'lon':37.2236, 'elevation':2168.})

   '''
   instnum = {'HARPS':14, 'HARPN':17,  'HPF':18, 'CARM_VIS':16, 'CARM_NIR':16, 'FIES':16}.get(inst, 0)

   if type(rammss) is str:
      if ':' in rammss:
         rammss = tuple(map(float,rammss.split(':')))
   if type(demmss) is str:
      if ':' in demmss:
         demmss = tuple(map(float,demmss.split(':')))

   par = tempfile.NamedTemporaryFile()
   dat = tempfile.NamedTemporaryFile()
   res = tempfile.NamedTemporaryFile()
   resfile = res.name
   res.close()
   #print "%10.4f%10.4f%10.4f"%rammss + "%10.4f%10.4f%10.4f"%demmss + "%5i%5i%10.4f%10.4f\n"%(epoch, 14, pma, pmd)
   par.write("this\n")
   par.write("%10.4f%10.4f%10.4f"%rammss + "%10.4f%10.4f%10.4f"%demmss + "%5i%5i%10.4f%10.4f\n"%(epoch, instnum, pma, pmd))
   if not instnum:
       # pass via obs location via file (instead fortran hardcoding)
       par.write("%s\n%s %s %s\n" % (inst, obsloc['lon'], obsloc['lat'], obsloc['elevation']))
   par.seek(0)
   par.name

   dat.write('1 ' + dateobs.replace(':',' ').replace('-',' ').replace('T',' ') + ' ' +str(exptime))
   dat.seek(0)

   p = Popen([path+'/bary.e'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
   #p = Popen(['/home/raid0/zechmeister/programs/BarCor/bary.e'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
   eph = p.communicate(input="\n".join([par.name,dat.name,res.name])+"\n"); eph
   result = np.loadtxt(res.name,skiprows=8)
   os.remove(res.name)
   result[-2] = result[-2] + 2400000.
   return result[-2:]

def main(files, **kwargs):
   for file in files:
      print file, fmt % tuple(bary_harps(file,**kwargs))

if __name__ == "__main__":
   parser = argparse.ArgumentParser(
      description='baricentric radial velocity correction',
      epilog="example:\n ./bary.py /media/5748244b-c6d4-49bb-97d4-1ae0a8ba6aed/data/HARPS/red_obj/gj551/HARPS.2013-05-08T02:39:50.479_A.fits\n")
   argopt = parser.add_argument # function short cut
   argopt('files', nargs='+',  help='filenames')
   argopt('-c',   nargs=2, help='coord rammss.sss [-]demmss.ss', default=[None, None], type=float, action='store')
   argopt('-pm',  nargs=2, help='[mas/yr] proper motions pmra pmde',default=[0.0, 0.0],  type=float, action='store')
   argopt('-fmt', help='output format', default="%.8f %f", type=str)

   args = parser.parse_args()
   globals().update(vars(args))

   sys.exit(main(files, ra=c[0], dec=c[1], pma=pm[0], pmd=pm[1]))
   #sys.exit(main(sys.argv[1:]))

