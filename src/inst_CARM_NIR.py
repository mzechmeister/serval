from __future__ import print_function

from read_spec import *
from read_spec import Inst
# Instrument parameters

name = 'CARM_NIR'
obsname = 'ca'
obsloc = dict(lat=37.2236, lon= -2.5463, elevation=2168.)

iomax = 28
iomax *= 2 # reshaping
pmax = 1800

oset = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 28, 29, 31, 46, 48, 50, 52]
coset = sorted(set(range(iomax)) - {17,18,19,20,21,35,36,37,38,39,40,41,42,43,44})

maskfile = 'telluric_mask_carm_short.dat'
maskfile = 'telluric_mask_nir4.dat'


def scan(self, s, pfits=True):
   """
   SYNTAX: read_carm_nir(filename)
   OUTPUT: namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55

   """
   HIERARCH = 'HIERARCH '
   self.hdulist = hdulist = pyfits.open(s) # slow 30 ms
   hdr = hdulist[0].header
   self.header = hdr
   self.instname = hdr['INSTRUME'][0:4]+'_NIR'
   #data,hdr = fitsio.read(s,header=True)

   self.drsberv = hdr.get('HIERARCH CARACAL BERV', np.nan)
   # BJD for stars, JD for for calibration products
   self.drsbjd = hdr.get('HIERARCH CARACAL BJD', hdr.get('MJD-OBS',0.0))
   self.drsbjd = hdr.get('HIERARCH CARACAL BJD', np.nan) + 2400000
   self.dateobs = hdr['DATE-OBS']
   self.dateobs = hdr['FILENAME'].replace("h",":").replace("m",":")
   self.dateobs = self.dateobs[4:8]+"-"+self.dateobs[8:10]+"-"+self.dateobs[10:21]
   self.mjd = hdr.get('HIERARCH CARACAL MJD-OBS')
   if not self.mjd:
      import warnings
      warnings.warn("Warning: keyword HIERARCH CARACAL MJD-OBS not found! This was implemented in CARACAL v2.00."+
                    "Please use lastest products.")
   if isinstance(self.drsbjd, str): self.drsbjd = 0.0   # workaround for MJD-OBS bug @2016-Jan
   self.fox = HIERARCH+'CARACAL FOX XWD' in hdr
   # check LED for NIR This order can be affect
   r = hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan) /  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan)
   if r > 1.5:
      self.flag |= sflag.led
      print(r, hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan), hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan), "# This spectrum could be affected by LED; or back background")
   self.sn55 = min(hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 16', np.nan),  hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 17', np.nan))
   if self.dateobs[:10] in ('2016-01-13', '2016-01-14', '2016-01-22', '2016-01-23'):
      self.sn55 = min(self.sn55, 10.) # bug in NIR fits file
   self.blaze = '' #hdr[HIERARCH+'ESO DRS BLAZE FILE']
   #self.drift = hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan)
   #self.e_drift = hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan)
   self.drift = hdr.get(HIERARCH+'CARACAL SERVAL FP RV', hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan))
   self.e_drift = hdr.get(HIERARCH+'CARACAL SERVAL FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan))
   self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
   self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)

   #fileid = str(hdr.ascardlist()['MJD-OBS'])     # read the comment
   self.fileid = hdr.get('FILENAME', 0) #fileid[fileid.index('(')+1:fileid.index(')')]
   self.calmode = hdr.get(HIERARCH+'CARACAL FIB','')
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]
   self.timeid = self.fileid
   self.ra = hdr['RA']
   self.de = hdr['DEC']
   self.airmass = hdr.get('AIRMASS', 0.0)
   self.exptime = hdr['EXPTIME']
   self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
   if self.exptime: self.tmmean /= self.exptime   # normalise
   if self.tmmean==0: self.tmmean = 0.5


def data(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      bp = np.array([[1,230],[1,313],[1,633],[1,882],[1,1330],[1,1411],[1,1513],[1,1546],[1,1584],[1,1613],[1,1771],[1,1840],
                        [2,230],[2,1015],[2,1024],[2,1070],[2,1089],[2,1146],[2,1158],[2,1257],[2,2256],[2,1674],[2,1695],[2,1707],[2,1771],[2,1799],[2,1842],
                        [3,818],[3,941],[3,1183],[3,1323],[3,1464],[3,1598],[3,1692],[3,1762],
                        [4,231],[4,237],[4,329],[4,331],[4,231],[4,237],[4,329],[4,331],[4,947],[4,1222],[4,1248],[4,1280],[4,1385],[4,1613],[4,1675],[4,1813],
                        [5,319],[5,407],[5,442],[5,452],[5,749],[5,832],[5,972],[5,1140],[5,1514],[5,1525],[5,1547],[5,1624],[5,1650],
                        [6,927],[6,1100],[6,1215],[6,1270],[6,1299],[6,1343],[6,1483],[6,1507],[6,1545],[6,1738],
                        [7,281],[7,660],[7,661],[7,837],[7,969],[7,1080],[7,1143],[7,1290],[7,1426],[7,1570],[7,1741],
                        [8,237],[8,266],[8,528],[8,564],[8,682],[8,730],[8,800],[8,1063],[8,1195],[8,1229],[8,1548],[8,1592],[8,1610],[8,1660],[8,1716],[8,1758],[8,1826],
                        [9,313],[9,323],[9,496],[9,582],[9,618],[9,627],[9,775],[9,887],[9,1359],[9,1387],[9,1693],[9,1799],
                        [10,175],[10,202],[10,339],[10,441],[10,541],[10,575],[10,693],[10,758],[10,862],[10,871],[10,918],[10,1174],[10,1381],[10,1500],[10,1845],
                        [11,296],[11,352],[11,420],[11,745],[11,746],[11,1152],[11,1237],[11,1395],[11,1495],[11,1643],[11,1681],[11,1771],[11,1831],
                        [12,194],[12,245],[12,274],[12,301],[12,309],[12,314],[12,333],[12,335],[12,401],[12,418],[12,497],[12,665],[12,764],[12,815],[12,935],[12,994],[12,1097],[12,1176],[12,1217],[12,1390],[12,1624],[12,1873],
                        [13,201],[13,510],[13,734],[13,777],[13,844],[13,979],[13,1238],[13,1463],[13,1501],[13,1919],
                        [14,250],[14,282],[14,339],[14,415],[14,482],[14,527],[14,766],[14,852],[14,875],[14,926],[14,967],[14,1121],[14,1134],[14,1143],[14,1237],[14,1238],[14,1774],[14,1798],
                        [15,247],[15,284],[15,400],[15,691],[15,902],[15,951],[15,961],[15,1036],[15,1047],[15,1152],[15,1153],[15,1179],[15,1185],[15,1243],[15,1247],[15,1449],[15,1531],[15,1532],
                        [16,228],[16,282],[16,363],[16,324],[16,363],[16,409],[16,520],[16,763],[16,877],[16,968],[16,1309],[16,1467],[16,1516],[16,1855],
                        [17,212],[17,553],[17,691],[17,760],[17,945],[17,1020],[17,1067],[17,1350],[17,1359],[17,1659],[17,1450],[17,1453],[17,1478],[17,1589],[17,1696],[17,1796],[17,1719],[17,1760],[17,1807],[17,1855],[17,2011],
                        [18,338],[18,597],[18,792],[18,818],[18,942],[18,975],[18,1028],[18,1149],[18,1311],[18,1348],[18,1355],[18,1604],[18,1724],[18,1825],
                        [19,198],[19,209],[19,244],[19,273],[19,1242],[19,1372],[19,1382],[19,1418],[19,1462],[19,1550],[19,1580],
                        [20,450],[20,624],[20,655],[20,966],[20,1042],[20,1065],[20,1075],[20,1136],[20,1560],[20,1591],[20,1651],[20,1655],[20,1657],
                        [21,301],[21,437],[21,450],[21,545],[21,581],[21,642],[21,766],[21,775],[21,782],[21,806],[21,851],[21,853],[21,1067],[21,1097],[21,1125],[21,1152],[21,1249],[21,1398],[21,1403],[21,1459],[21,1581],[21,1620],[21,1823],[21,1865],[21,1924],
                        [22,383],[22,391],[22,407],[22,441],[22,644],[22,680],[22,843],[22,1007],[22,1021],[22,1022],[22,1023],[22,1031],[22,1290],[22,1944],[22,1983],[22,1996],
                        [23,304],[23,716],[23,742],[23,761],[23,821],[23,935],[23,945],[23,1029],[23,1266],[23,1560],[23,1586],
                        [24,356],[24,473],[24,502],[24,599],[24,688],[24,704],[24,744],[24,758],[24,767],[24,1019],[24,1268],[24,1270],[24,1576],[24,1639],[24,1751],[24,1860],
                        [25,190],[25,430],[25,827],[25,936],[25,938],[25,1375],[25,1508],[25,1583],[25,1670],[25,1681],[25,1714],[25,1724],[25,1926],[25,1927],[25,1974],
                        [26,802],[26,889],[26,1191],[26,1291],[26,1435],[26,1508],[26,1636],[26,1740],[26,1930],
                        [1,2733],[1,2654],[1,2672],[2,3644], #SCA1
                        [3,2765],[3,3543],[3,3544],[3,3624],
                        [4,2881],[4,2687],
                        [7,3346],
                        [8,2404],[8,2405],
                        [9,2813],[9,2814],
                        [10,2455],[10,3675],
                        [12,3152],
                        [13,2845],[13,3572],
                        [14,2590],[14,2652],[14,3130],
                        [15,2473],[15,3064],[15,3083],[15,3738],[15,3739],
                        [16,3131],[16,3231],[16,3130],[16,3499],[16,3626],
                        [17,2518],[17,2519],[17,2814],[17,3202],[17,3325],[17,3740],
                        [18,2383],[18,2514],[18,2739],[18,3390],[18,3568],[18,3744],
                        [19,3321],[19,3332],[19,3481],[19,3619],
                        [21,2715],[21,3110],[21,3250],[21,3673],
                        [23,3333],[23,3606],[23,3722],
                        [24,2722],[24,3205],[24,3310],[24,3552],
                        [25,2755],[25,3058],[25,3786],
                        [27,2440],[27,3616],[27,3792]
                        ]).T
      if orders == np.s_[:]:
         #f = hdulist['SPEC'].data
         # "data" atribute seems to open again the fits file. For large data set (GJ273) this lead to "error: too many files open". So use "section"
         f = 1.*hdulist['SPEC'].section[orders]
         dim = f.shape
         if 1: # reshaping
            dim = (dim[0]*2,dim[1]//2)
            f = f.reshape(dim)
            bp[0], bp[1] = bp[0]*2, bp[1]%dim[1]

         #w = hdulist['WAVE'].data.reshape(dim)
         #e = hdulist['SIG'].data.reshape(dim)
         w = hdulist['WAVE'].section[orders].reshape(dim)
         e = hdulist['SIG'].section[orders].reshape(dim)
         bpmap = np.isnan(f).astype(int)            # flag 1 for nan
         # a hand made bad pixel list car-20160218T18h48m52s-sci-gtoc-nir.fits
         bpmap[bp[0],bp[1]] = 1
      # interpolate bad columns, they mess up a lot the creation of the template
      # We do this only when we read all order (preRVs), but not for coadding (single orders)
         if 1: # int
            #g = 1.0 * f
            f[bp[0],bp[1]] = (f[bp[0],bp[1]-1]+f[bp[0],bp[1]+1]) / 2
            #o=17; gplot(g[o],'w l,',f[o],  'w lp')
            #pause()
      else:
         orders, det = divmod(orders,2)
         f = 1.*hdulist['SPEC'].section[orders]
         if 1:
            # split and renumber the orders and indidexs
            sl = [None, None]
            sl[1-det] = f.shape[0]/2
            sl = slice(*sl)
            f = f[sl]
            bp[0], bp[1] = bp[0]*2, bp[1]%f.shape[0]

         w = hdulist['WAVE'].section[orders][sl]
         e = hdulist['SIG'].section[orders][sl]

         bpmap = np.isnan(f).astype(int)            # flag 1 for nan
         bp = bp[:,np.where(bp[0]==orders)[0]]
         if len(bp):
           bpmap[bp[1]] = 1
           f[bp[1]] = (f[bp[1]-1]+f[bp[1]+1]) / 2
         #pause()
      if self.fox:
         # scale spectrum
         e = e * 100000.
         f = f * 100000.
      else:
         # fib B, linear extraction
         e = 0*f + np.sqrt(np.nanmedian(f)) # unweighted maybe for HCL
         bpmap[f>300000] |= flag.sat
         #bpmap[:,1:] |= flag.sat * (f>300000)[:,:-1]
         #bpmap[:,:-1] |= flag.sat * (f>300000)[:,1:]
      bpmap[f < -3*e] |= flag.neg
      bpmap[e==0] |= flag.nan
      return w, f, e, bpmap

