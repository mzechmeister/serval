from read_spec import *
from arrays import Arrays

# Instrument parameters
name = __name__[5:]
obsloc = dict(lat=42.9333, lon= 0.1333, elevation=2869.4)
obsname = None
pat = '*.fits'

iomax = 40
#pmax =  7800 - 300
pmax =  -300
#oset = ':66'
snmax = 1900

maskfile = 'telluric_mask_carm_short.dat'

# Instrument read functions
def scan(self, s, pfits=True):
   with open(s[:-1]+'meta') as f:
       hdrcards = f.read()
   hdr = dict(x.replace("'", "").split(" = ") for x in hdrcards.split("\n"))

   self.instname = hdr['Instrument']
   self.drsberv = float(hdr.get('Helvel', np.nan)[:-5])  # strip km/s
   self.drsbjd = hdr.get('BTLA', np.nan)
   self.dateobs = hdr['Observation date']
   self.mjd = float(hdr.get('Julian Date'))
   self.drift = np.nan
   self.e_drift = np.nan
   self.sn55 = float(hdr.get('SnrMax', 50))   
   self.fileid = self.dateobs
   self.calmode = "%s,%s,%s" % (hdr.get('Instrumental Mode', ''), hdr.get('Stokes', ''), hdr.get('P_NAME2', ''))
   self.timeid = self.fileid

   self.ra = float(hdr['RA'])
   self.de = float(hdr['DEC'])
   self.airmass = float(hdr.get('Airmass', np.nan))
   self.exptime = float(hdr['Texposure'].strip(' sec '))
   self.tmmean = 0.5

   self.header = pyfits.Header({'HIERARCH '+k:v for k,v in hdr.items()})
   self.header['OBJECT'] = hdr['Object name']

def data(self, orders=None, pfits=True):
   # w_air, f, e = np.loadtxt(self.filename, skiprows=2, unpack=True) # very slow
   with open(self.filename) as f:
       f.readline()
       f.readline()
       _data = np.fromfile(f, sep=" ")

   w_air, f, e = _data.reshape(-1,3).T
   w = airtovac(w_air*10*(1.000-self.drsberv/3e5 ))  # or /(1+berv)?
 
   # find the jumps and split into orders
   idx = list(np.where(abs(np.diff(w)) > 0.1)[0]+1)
   oidx = [slice(*sl) for sl in zip([0]+idx, idx+[None])]

   idx = [0] + idx + [len(f)]
   if 0:
       # padding orders with nans
       nx = np.diff(idx).max()
       no = len(oidx)
       f_ = f
       w_ = w
       e_ = e
   
       self.f = f = np.zeros((no,nx)) * np.nan
       self.w = w = np.zeros((no,nx)) * np.nan
       self.e = e = np.zeros((no,nx)) * np.nan
   
       for o in np.arange(no):
           w[o,:idx[o+1]-idx[o]] = w_[oidx[o]]
           f[o,:idx[o+1]-idx[o]] = f_[oidx[o]]
           e[o,:idx[o+1]-idx[o]] = e_[oidx[o]]

   else:
       self.w = w = Arrays([w[oi] for oi in oidx])
       self.f = f = Arrays([f[oi] for oi in oidx])
       self.e = e = Arrays([e[oi] for oi in oidx])

   self.bpmap = bpmap = np.isnan(f).astype(int)            # flag 1 for nan
   with np.errstate(invalid='ignore'):
       bpmap[f < -3*e] |= flag.neg
       bpmap[e==0] |= flag.nan

   return w[orders], f[orders], e[orders], bpmap[orders]


