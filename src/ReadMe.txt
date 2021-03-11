Format Description

default format (*.dat): columns separated by blanks time series
*.cvs:                  columns separated by semicolon

File Summary:
--------------------------------------------------------------------------------
FileName       Explanations
--------------------------------------------------------------------------------
obj.rvc.dat    Radial velocity (drift and secular acceleration corrected, nan-drifts treated as zero)
obj.brv.dat    Barycentric earth radial velocity
obj.cairt.dat  CAII IRT line index (requires absolute RVs)
obj.chi.dat    sqrt(reduced chi^2) (orderwise)
obj.crx.dat    Chromatic RV index (wavelength depedence of the RV)
obj.dat        Radial velocity (drift and secular acceleration corrected, nan-drift excluded)
obj.dlw.dat    Differential line width (orderwise)
obj.drs.dat    Online radial velocity from fits header
obj.e_dlw.dat  Error for differential line width (orderwise)
obj.fits       Template
obj.halpha.dat Halpha line index (requires absolute RVs)
obj.info.cvs   Infomation from fits file
obj.mlc.dat    RV and CRX averages via maximum likelihood maps (experimental)
obj.nad.dat    NaD line index (requires absolute RVs)
obj.post.dat   post processed RVs with re-weighting of orders (experimental)
obj.pre.dat    RVs measured against highest S/N spectrum, before coadding
obj.rvo.dat    Radial velocity (orderwise, not drift corrected)
obj.snr.dat    Signal-to-noise (orderwise)
obj.srv.dat    Serval products time series compilation
obj.targ.cvs   Target properties from catalog request
lastcmd.txt    Used SERVAL options
log.obj        Plain text logfile of SERVAL stdout
--------------------------------------------------------------------------------
NOTES:
[1]   (BJD,BERV) depending on serval input option, propagated from CARACAL
      or recomputed by SERVAL (for -brvref WE: BJD_TDB, MH: BJD_UTC)
[2]   Reference region for CAIRT1 are 8492 (-40, 40) km/s and 8504 (-40, 40) km/s
[3]   Reference region for CAIRT2 are 8542.09 (-300, -200) km/s and 8542.09 (250, 350) km/s
[4]   Reference region for CAIRT3 are 8662.14 (-300, -200) km/s and 8662.14 (200, 300) km/s
[5]   Reference region for NAD1 are 5885 (-40, 40) km/s and 5892.94 (200, 300) km/s
[6]   Reference region for NAD2 are 5892.94 (200, 300) km/s and 5905 (-40, 40) km/s

Description of file: obj.rvc.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVC       Radial velocity (drift and sa corrected)
     3 D      m/s     E_RVC       Radial velocity error
     4 D      m/s       DRIFT     CARACAL drift measure
     5 D      m/s     E_DRIFT     CARACAL drift measure error
     6 D      m/s       RV        Radial velocity
     7 D      m/s     E_RV        Radial velocity error
     8 D      km/s      BERV      Barycentric earth radial velocity [1]
     9 D      m/s       SADRIFT   Drift due to secular acceleration
--------------------------------------------------------------------------------


Description of file: obj.brv.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      km/s      BERV      Barycentric earth radial velocity [1]
     3 D      ---       DRSBJD    Barycentric Julian date (CARACAL)
     4 D      km/s      DRSBERV   Barycentric earth radial velocity (CARACAL)
     5 D      m/s       DRSDRIFT  CARACAL drift measure in fib B
     6 A      ---       TIMEID    Identifier for file (time stamp)
     7 D      ---       TMMEAN    Flux weighted mean point
     8 D      s         EXPTIME   Exposure time
--------------------------------------------------------------------------------


Description of file: obj.cairt.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       CAIRT1    CaII-IRT1 index (8498.02) (-15.,15.) km/s [2]  # NIST + my definition
     3 D      ---     E_CAIRT1    CAIRT1 error
     4 D      ---       CAIRT2    CaII-IRT2 index (8542.09) (-15.,15.) km/s [3]   # NIST + my definition
     5 D      ---     E_CAIRT2    CAIRT2 error
     6 D      ---       CAIRT3    CaII-IRT3 index (8662.14) (-15.,15.) km/s [4]   # NIST + my definition
     7 D      ---     E_CAIRT3    CAIRT3 error
--------------------------------------------------------------------------------


Description of file: obj.chi.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       RCHI      Overall reduced chi^2
     3 D      ---       RCHIO_00  Sqrt(reduced chi^2) in order 0
     4 D      ---       RCHIO_01  Sqrt(reduced chi^2) in order 1
     5 D      ---       RCHIO_02  Sqrt(reduced chi^2) in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.crx.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       CRX       Chromatic index (Slope over logarithmic wavelength)
     3 D      m/s     E_CRX       Slope error
     4 D      m/s       CRX_OFF   Offset parameter alpha
     5 D      m/s     E_CRX_OFF   Error for offset parameter alpha
     6 D      A         L_V       Wavelength at which the slope intersects RV
     7 D                LNL_00    Central logarithmic wavelength of aperture 1
     8 D                LNL_01    Central logarithmic wavelength of aperture 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVC       Radial velocity (drift and sa corrected)
     3 D      m/s     E_RVC       Radial velocity error
--------------------------------------------------------------------------------


Description of file: obj.dlw.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       dLW       differential line width
     3 D      ---     E_dLW       error for dLW
     4 D      ---       dLWO_00   differential line width in order 0
     5 D      ---       dLWO_01   differential line width in order 1
     6 D      ---       dLWO_02   differential line width in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.drs.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVDRS     Radial velocity (from CARACAL or HARPS DRS)
     3 D      m/s     E_RVDRS     Radial velocity error (from CARACAL or HARPS DRS)
     4 D      km/s      FWHM      Full half maximum (HARPS DRS only)
     5 D      km/s      BIS       Bisector span (HARPS DRS only)
     6 D      %         CONTRAST  Contrast (HARPS DRS only)
     7 A      ---       TIMEID    FILENAME from fits header
--------------------------------------------------------------------------------


Description of file: obj.e_dlw.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       dLW       Differential line width
     3 D      ---     E_dLW       Error for dLW
     4 D      ---     E_dLWO_00   Error for differential line width in order 0
     5 D      ---     E_dLWO_01   Error for differential line width in order 1
     6 D      ---     E_dLWO_02   Error for differential line width in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.halpha.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       HALPHA    Halpha index (6562.808) (-40,40) km/s
     3 D      ---     E_HALPHA    HALPHA error
     4 D      ---       HACEN     Halpha mean flux
     5 D      ---     E_HACEN     HACEN error
     6 D      ---       HALEFT    Halpha reference region (-300,-100) km/s
     7 D      ---     E_HALEFT    HALEFT error
     8 D      ---       HARIGH    Halpha reference region (100,300) km/s
     9 D      ---     E_HARIGH    HARIGH error
    10 D      ---       CAI       Calcium I index
    11 D      ---     E_CAI       CAI error
--------------------------------------------------------------------------------


Description of file: obj.info.cvs
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 A      ---       TIMEID    Identifier for file (time stamp)
     2 D      ---       BJD       Barycentric Julian date [1]
     3 D      ---       SNREF     Signal-to-ratio in reference order (CARM_VIS: 36, CARM_NIR: 16, HARPS: 55)
     4 A      ---       OBJ       Fits keyword OBJECT
     5 D      s         EXPTIME   Exposure time
     6 A      ---       SPT       Spectral type SpT ccf.mask
     7 I256   ---       FLAG      Flag bad spectra/instrument mode (0: ok)
     8 D      ---       AIRMASS   Air mass from fits header
     9 D      deg       RA        Telescope Ascension (degrees)
    10 D      deg       DE        Telescope declination (degrees)
--------------------------------------------------------------------------------


Description of file: obj.mlc.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       MLRVC     ML Radial velocity (drift and sa corrected)
     3 D      m/s     E_MLRVC     ML Radial velocity error
     4 D      m/s       MLCRX     ML Chromatic index (Slope over logarithmic wavelength)
     5 D      m/s     E_MLCRX     error for MLCRX (slope error)
     6 D      m^2/s^2   DLW       Differential Line Width
     7 D      m^2/s^2 E_DLW       Error in DLW
--------------------------------------------------------------------------------


Description of file: obj.nad.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       NAD1      NaD1 index (5889.950943) (-15,15) km/s [5]
     3 D      ---     E_NAD1      NAD1 error
     4 D      ---       NAD2      NaD2 index (5895.924237) (-15,15) km/s [6]
     5 D      ---     E_NAD2      NAD2 error
--------------------------------------------------------------------------------


Description of file: obj.post.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVRC      Radial velocity (reweighted, drift corrected)
     3 D      m/s     E_RVRC      Radial velocity error
     4 D      m/s       RVR       Radial velocity (reweighted)
     5 D      m/s     E_RVR       Radial velocity error
     6 D      m/s       DRIFT     CARACAL drift measure
     7 D      m/s     E_DRIFT     CARACAL drift measure error
     8 D      m/s       BERV      Barycentric earth radial velocity [1]
     9 D      m/s       SADRIFT   Drift due to secular acceleration
--------------------------------------------------------------------------------


Description of file: obj.pre.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVPRE     Radial velocity against highest S/N spectrum
     3 D      m/s     E_RVPRE     Radial velocity error
--------------------------------------------------------------------------------


Description of file: obj.rvo.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RV        Radial velocity (mean, not drift and not sa corrected)
     3 D      m/s     E_RV        Radial velocity error
     4 D      m/s       RVMED     Radial velocity (median)
     5 D      m/s     E_RVMED    Radial velocity (median) error
     6 D      m/s       RVO_00    Radial velocity in order 0
     7 D      m/s       RVO_01    Radial velocity in order 1
     8 D      m/s       RVO_02    Radial velocity in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.rvo.daterr
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RV        Radial velocity (mean, not drift and not sa corrected)
     3 D      m/s     E_RV        Radial velocity error
     4 D      m/s       RVMED     Radial velocity (median)
     5 D      m/s     E_RVMED     Radial velocity (median) error
     6 D      m/s     E_RVO_00    Radial velocity error in order 0
     7 D      m/s     E_RVO_01    Radial velocity error in order 1
     8 D      m/s     E_RVO_02    Radial velocity error in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.snr.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      ---       SNR       Overall signal to noise
     3 D      ---       SNRO_00   Signal to noise in order 0
     4 D      ---       SNRO_01   Signal to noise in order 1
     5 D      ---       SNRO_02   Signal to noise in order 2
etc...
--------------------------------------------------------------------------------


Description of file: obj.srv.dat
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 D      ---       BJD       Barycentric Julian date [1]
     2 D      m/s       RVC       Radial velocity (drift and sa corrected)
     3 D      m/s     E_RVC       Radial velocity error
     4 D      m/s       CRX       Chromatic index (Slope over logarithmic wavelength)
     5 D      m/s     E_CRX       error for CRX (slope error)
     6 D      m^2/s^2   DLW       Differential Line Width
     7 D      m^2/s^2 E_DLW       Error in DLW
--------------------------------------------------------------------------------

Description of file: obj.targ.cvs
--------------------------------------------------------------------------------
Column Format Units     Label     Explanations
--------------------------------------------------------------------------------
     1 A      ---       OBJECT    Name requested in simbad
     2 A      ---       ID        Simbad main identifier
     3 A      ---       COO       Coordinates RA DE
     4 A      mas/yr    PM        Proper motion pmRA pmDE [error ellipse]
     6 A      mas       PLX       Parallax [error] quality bibcode
     7 A      km/s      RV        Absolute RV (wavelength) quality [error] bibcode
--------------------------------------------------------------------------------


Description of file: obj.fits
--------------------------------------------------------------------------------
Extension Format Units     Label     Explanations
--------------------------------------------------------------------------------
      [0] A       ---                primary header (keywords for coadd method, used
                                     spectra and RVs, template SNR for each order)
      [1] D       ---      SPEC      flux
      [2] D       A        WAVE      wavelength


template_post3.fits	template from coadd method post3	same as for template.fits

