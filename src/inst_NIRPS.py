from read_spec import *
from calcspec import redshift

inst    = 'NIRPS'
name    = 'NIRPS'
obsname = 'lasilla'
obsloc  = dict(lat=-29.2567, lon=-70.7337, elevation=2400)

iomax   = 71
snmax   = 500
oset    = ':'

maskfile = 'telluric_mask_carm_short.dat'

pat = '*.tar *_S2D_SKYSUB_TELLSUB_A.fits'



def scan(self, s, pfits=True, verb=False):

    drs = self.drs
    if '.tar' in s:
        s = file_from_tar(s, inst=inst, fib=self.fib, pfits=pfits)
      
    # Support *.tar containers via the helper in read_spec
    if isinstance(s, str) and s.endswith('.tar'):
        fib = getattr(self, 'fib', 'A') or 'A'
        s = file_from_tar(s, inst=name, fib=fib, pat='_S2D_%s' % fib, pfits=pfits)


    HIERARCH = 'HIERARCH '
    HIERINST = HIERARCH + 'ESO '
    k_tmmean = HIERINST + 'DET EM1 TMMEAN'
    
    self.HIERQC = HIERQC = HIERINST + 'QC '

    if pfits is True:  # header with pyfits
        self.hdulist = hdulist = pyfits.open(s)
        hdr = self.hdulist[0].header

    self.instname = hdr['INSTRUME']
    self.HIERARCH = HIERARCH

    self.airmass = np.nanmean([hdr.get('HIERARCH ESO TEL AIRM END', np.nan), hdr.get('HIERARCH ESO TEL AIRM START', np.nan)])
    self.exptime = hdr['EXPTIME']
    self.mjd = hdr['MJD-OBS']
    self.dateobs = hdr['DATE-OBS']
    self.ra = hdr['RA']
    self.de = hdr['DEC']
    self.utc = datetime.datetime.strptime(self.dateobs, '%Y-%m-%dT%H:%M:%S')

    self.obs.lon = -70.7345
    self.obs.lat = -29.2584

    if k_tmmean not in hdr:
        warnings.warn('Warning: Setting tmmean to 0.5!')
    self.tmmean = hdr.get(k_tmmean, 0.5)

    self.DRS = hdr.get('HIERARCH ESO PRO REC1 PIPE ID')

    self.drsbjd = hdr.get(HIERQC + 'BJD', np.nan)
    self.drsberv = hdr.get(HIERQC + 'BERV', np.nan)

    self.sn55 = hdr.get(HIERQC+'ORDER59 SNR', np.nan)
    self.blaze = hdr.get(HIERINST+'PRO REC1 CAL23 NAME', 0)
    self.drift = hdr.get(HIERINST+'PRO REC1 PARAM16 VALUE', np.nan)

    self.fileid = self.dateobs
    self.timeid = self.fileid
    self.calmode = hdr.get('DPRTYPE', 'NOT FOUND')
    
    self.drift = np.nan
    self.ccf.rvc = hdr.get('HIERARCH ESO QC CCF RV', np.nan)
    self.ccf.err_rvc = hdr.get('HIERARCH ESO QC CCF RV ERROR', np.nan)

    hdr['OBJECT'] = hdr.get('OBJECT', 'FOX')
    self.header = self.hdr = hdr # self.header will be set to None
    
    
def data(self, orders, pfits=True):
    hdulist = self.hdulist
    
    f = hdulist['SCIDATA'].section[orders]
    e = hdulist['ERRDATA'].section[orders]
    w = redshift(hdulist['WAVEDATA_VAC_BARY'].section[orders], ve=self.drsberv, wlog=False)   # undo drsberv correction

    bpmap = np.isnan(f).astype(int)   # flag 1 for nan

    with np.errstate(invalid='ignore'):
        bpmap[f < -3*e] |= flag.neg
        bpmap[f > 300000] |= flag.sat
        bpmap[e==0] |= flag.nan
    
    return w, f, e, bpmap
