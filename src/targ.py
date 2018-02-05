import urllib, urllib2
import os
import sys

# python version for the shell query
# star=`sed 's/ /_/g; s/+/%2B/g' <<<$@`
# echo %Info: encode name: $star >&2
# wget -q -O - 'http://simbad.u-strasbg.fr/simbad/sim-script?submit=submit+script&script=output+script%3Doff%0D%0Aoutput+console%3Doff%0D%0Aformat+object+form1+%22%25OBJECT+%3A+%25IDLIST%281%29+%3A+%25COO%28A+D%29+%25PM%28A+D+[E]%29+%25PLX\n%22%0D%0Aquery+'$star | tee | sed 's/:[^:]*//; q' | awk '{$2=": nan"; if ($14!="~") $2=": "22.98*(($9/1000)^2+($10/1000)^2)/$14; print}'

# See also:
#    http://simbad.u-strasbg.fr/simbad/sim-fscript
#
#    output script=off
#    output console=off
#    format object form1 "%OBJECT : %IDLIST(1) : %COO(A D) %PM(A D [E]) %PLX\n"
#    query gj699

def simbad_query(targ):
   '''
   Make a request to simbad to get RADE, PM and PLX.
   
   targ : str, identifier name for simbad request.

   Example:
   --------
   >>> simbad_query('GJ699')
   "GJ699 : NAME Barnard's star : 17 57 48.49803 +04 41 36.2072 -798.58  10328.12  [1.72 1.22 0] 548.31 [1.51] A 2007A&A...474..653V\n\n"
   '''

   query = urllib.quote("\n".join([
      'output script=off',
      'output console=off',
      'format object form1 "%OBJECT;%IDLIST(1);%COO(A D);%PM(A D [E]);%PLX;%RV"',
      'query '+targ]))
   # urllib2.urlopen('http://simbad.u-strasbg.fr/simbad/sim-script', 'submit=submit+script&script=output+script%3Doff%0D%0Aoutput+console%3Doff%0D%0Aformat+object+form1+%22%25OBJECT+%3A+%25IDLIST%281%29+%3A+%25COO%28A+D%29+%25PM%28A+D+[E]%29+%25PLX\n%22%0D%0Aquery+gj699')

   result = urllib2.urlopen('http://simbad.u-strasbg.fr/simbad/sim-script', 'submit=submit+script&script='+query).read()
   return result

class Targ:
   '''Properties of the target.

   Store coordinate information, etc. as attributes.

   Example:
   --------
   >>> targ = Targ('GJ699')
   targ.py: Requesting simbad for 'gj699'
   >>> targ.ra
   (17.0, 57.0, 48.49803)

   '''
   def __init__(self, name, rade=(None, None), pm=(None, None), plx=None, rv=None, sa=float('nan'), cvs=None):
      self.name = name
      self.sa = sa
      self.ra, self.de = rade
      self.pmra, self.pmde = pm
      if self.ra and self.de:
         # apply user input from command line
         self.ra = tuple(map(float,self.ra.split(':')))
         self.de = tuple(map(float,self.de.split(':')))
      else:
         # look for name, try to read from file. If not found or different object then make a request to simbad
         if not self.fromfile(cvs) or not self.line.startswith(self.name+";"):
            self.query()
            self.tofile(cvs)
         self.assignAttr(self.line)
      if self.pmra and self.plx:
         self.sa = 22.98 * ((self.pmra/1000)**2+(self.pmde/1000)**2) / self.plx

   def fromfile(self, filename):
      '''Restore info from a file.'''
      self.line = None
      if os.path.exists(filename):
         print "targ.py: restoring '%s' from %s" % (self.name, filename)
         with open(filename) as f:
            self.line = f.read()
      return self.line

   def query(self):
      print "targ.py: requesting simbad for '%s'" % self.name
      self.line = simbad_query(self.name)

   def assignAttr(self, line):
      # parse the request
      line = self.line.split(';')[2:]        # ['gj699', "NAME Barnard's star", ' 17 57 ...]
      line = " ".join(line).split()
      self.ra = tuple(map(float,line[0:3]))  # rammss = (14.,29.,42.94)
      self.de = tuple(map(float,line[3:6]))  # demmss = (-62.,40.,46.16)
      self.pmra = float(line[6].replace("~","0."))             # pma = -3775.75
      self.pmde = float(line[7].replace("~","0."))             # pmd = 765.54
      self.plx = float(line[11].replace("~","nan"))
      self.rv = float(line[16].replace("~","nan"))

   def tofile(self, filename=None):
      if filename:
         with (open(filename, 'w') if filename else sys.stdout) as f:
            print >>f, self.line
         print 'storing in', filename
      else:
         print self.line


if __name__ == "__main__":
   name = 'gj699'
   if len(sys.argv): name = sys.argv[1]
   targ = Targ(name)
   #targ = Targ('gj699', fromfilename='bla')
   print targ.sa, targ.pmra, targ.pmde, targ.plx


