from read_spec import *
#from read_spec import Inst
# Instrument parameters

name = 'CARM_VIS'
obsname = 'ca'
obsloc = dict(lat=37.2236, lon= -2.5463, elevation=2168.)

iomax = 61 # NAXIS2
snmax = 500
oset = '10:52'

maskfile = 'telluric_mask_carm_short.dat'

def scan(self, s, pfits=True):
   """
   Returns
   -------
   namedtuple('spectrum', 'w f berv bjd blaze drift timeid sn55 ')
           w    - wavelength
           f    - flux
           berv - Barycentric Earth Radial Velocity
           bjd  - Barycentric Julian Day
           blaze - Blaze filename
           drift - Used RV Drift
           sn55  - S_N order center55
   Example
   -------
   >>> read_carm_vis(filename)

   """
   HIERARCH = 'HIERARCH '
   hdulist = self.hdulist = pyfits.open(s) # slow 30 ms
   if 1:
      self.header = hdr = hdulist[0].header
      if 'HIERARCH CARACAL DRIFT FP REF' in hdr: del hdr['HIERARCH CARACAL DRIFT FP REF']
      self.instname = hdr['INSTRUME'][0:4] + '_VIS'
      self.drsberv = hdr.get('HIERARCH CARACAL BERV', np.nan)
      # BJD for stars, JD for for calibration products
      self.drsbjd = hdr.get('HIERARCH CARACAL BJD', np.nan) + 2400000
      #if isinstance(self.drsbjd, str):
         #self.drsbjd = 0.0   # workaround for MJD-OBS bugs (empty or missing fractional digits) @2016-Jan
      self.dateobs = hdr['DATE-OBS']
      # dateobs is used by MH, but date-obs seems more reliable from FILENAME
      # CARACAL computes mjd-obs also from FILENAME
      self.dateobs = hdr['FILENAME'].replace("h",":").replace("m",":")
      self.dateobs = self.dateobs[4:8]+"-"+self.dateobs[8:10]+"-"+self.dateobs[10:21]
      self.mjd = hdr.get('HIERARCH CARACAL MJD-OBS')
      if not self.mjd:
         import warnings
         warnings.warn("Warning: keyword HIERARCH CARACAL MJD-OBS not found! This was implemented in CARACAL v2.00."+
                       "Please use lastest products.")
      self.drift = hdr.get(HIERARCH+'CARACAL SERVAL FP RV', hdr.get(HIERARCH+'CARACAL DRIFT FP RV', np.nan))
      self.e_drift = hdr.get(HIERARCH+'CARACAL SERVAL FP E_RV', hdr.get(HIERARCH+'CARACAL DRIFT FP E_RV', np.nan))
      self.fox = HIERARCH+'CARACAL FOX XWD' in hdr
      self.sn55 = hdr.get(HIERARCH+'CARACAL '+('FOX' if self.fox else 'LXT')+' SNR 36', np.nan)   # @ 746nm
      sn25 = hdr.get(HIERARCH+'CARACAL FOX SNR 25', np.nan)
      sn30 = hdr.get(HIERARCH+'CARACAL FOX SNR 30', np.nan)
      if sn25 > 70 and sn25 > 10*sn30: # hig
         self.flag |= sflag.led
         #print(sn25, sn30, self.sn55, )

      self.fileid = hdr.get('FILENAME', 0)
      self.timeid = self.fileid
      self.calmode = hdr.get('SOURCE', '') #.split("_")[3] #fileid[fileid.index('(')+1:fileid.index(')')]
      self.calmode = hdr.get(HIERARCH+'CARACAL FIB', '')
   #calmodedict = {'objcal':'OBJ,CAL','objsky':'OBJ,SKY'}
   #if calmode in calmodedict: calmode = calmodedict[calmode]

      self.ccf.rvc = hdr.get(HIERARCH+'CARACAL SERVAL RV', np.nan)
      self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL SERVAL E_RV', np.nan)
      #self.ccf.rvc = hdr.get(HIERARCH+'CARACAL CCF RV', np.nan)
      #self.ccf.err_rvc = hdr.get(HIERARCH+'CARACAL CCF E_RV', np.nan)

      self.ra = hdr['RA']
      self.de = hdr['DEC']
      self.airmass = hdr.get('AIRMASS', np.nan)
      self.exptime = hdr['EXPTIME']
      self.tmmean = hdr.get(HIERARCH+'CARACAL TMEAN', 0.0)
      if self.exptime: self.tmmean /= self.exptime   # normalise
      if self.tmmean == 0: self.tmmean = 0.5


def data(self, orders, pfits=True):
   hdulist = self.hdulist
   if 1:  # read order data
      f = hdulist['SPEC'].section[orders]
      w = hdulist['WAVE'].section[orders]
      e = hdulist['SIG'].section[orders]
      bpmap = np.isnan(f).astype(np.uint64)            # flag 1 for nan

      bpmap0 = np.zeros((61,4096), dtype=np.uint64)
      bpmap0[14:38,[2453-3, 2453-2, 2453-1, 2453, 2453+1, 2453+2, 2453+3]] |= 1
      bpmap0[14:38,1643] |= 1   # ghost of hotspot tail
      bpmap0[14:38,2459] |= 1   # spikes of hotspot satellite (bug not correct due to bug in v2.00)
      bpmap0[15:41,3374] |= 1   # displaced column; ignore by marking as nan
      bpmap0[28,3395:3400] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[34,838:850] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis.fits
      bpmap0[34,2035:2044] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[34,3150:3161] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis.fits
      bpmap0[35,403:410] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,754:759] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[35,1083:1093] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,1944:1956] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,2710:2715] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[35,3050:3070] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,3706:3717] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[35,3706:3717] |= flag.sky # car-20160803T22h46m41s-sci-gtoc-vis
      bpmap0[36,303:308] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,312:317] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,1311:1315] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[36,1325:1329] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[37,1326:1343] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[39,1076:1082] |= flag.sky # car-20170626T02h00m17s-sci-gtoc-vis
      bpmap0[39,1204:1212] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,1236:1243] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[39,1463:1468] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,2196:2203] |= flag.sky # car-20160520T03h10m13s-sci-gtoc-vis.fits
      bpmap0[39,2493:2504] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[39,3705:3717] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[40,2765:2773] |= flag.sky # car-20170419T03h27m48s-sci-gtoc-vis
      bpmap0[40,3146:3153] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[40,3556:3564] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,486:491] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,495:501] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[41,1305:1315] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,480:490] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,1316:1330] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,2363:2368] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[42,2375:2382] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[44,3355:3361] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,311:321] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[46,835:845] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,1156:1171] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[46,1895:1905] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[46,2212:2232] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[47,2127:2133] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2218:2223] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2260:2266] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,2313:2319] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3111:3116] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3267:3272] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3316:3321] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3432:3438] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[47,3480:3488] |= flag.sky # car-20160714T00h18m29s-sci-gtoc-vis
      bpmap0[47,3658:3665] |= flag.sky # car-20170509T03h05m21s-sci-gtoc-vis
      bpmap0[49,1008:1017] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,2532:2544] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,3046:3056] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      bpmap0[49,3574:3588] |= flag.sky # car-20160701T00h49m36s-sci-gtoc-vis.fits
      # interpolate bad columns, they mess up a lot the creation of the template
      # We do this only when we read all order (preRVs), but not for coadding (single orders)
      if orders == np.s_[:]:
         # the hotspot
         #f[14:38,2453-3: 2453+4] = f[14:38,2453-4][:,np.newaxis] +  (f[14:38,2453+4]-f[14:38,2453-4]).reshape(24,1)*np.arange(1,8).reshape(1,7)/8.
         #f[15:41,3374] = np.nan # 0.5*(f[15:41,3374-1]+f[15:41,3374+1])
         pass


      bpmap |= bpmap0[orders]
      bpmap = bpmap.astype(int)

      with np.errstate(invalid='ignore'):
        # arrgh, in newer numpy version comparison with nan raises a warning
        if self.fox:
           e = e * 10.
           f = f * 10.
        else:
           e = np.sqrt(5.*10 + np.abs(f))
           bpmap[f>300000] |= flag.sat
        bpmap[f < -3*e] |= flag.neg
        bpmap[e==0] |= flag.nan

      return w, f, e, bpmap

