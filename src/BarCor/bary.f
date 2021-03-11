      PROGRAM BARCOR
C===========================================================================
C     
C     Last release: August 25, 2009
C
C     BarCor computes barycentric corrections of radial velocities and time
C     for input UT date and time. The program is written in Fortran77.      
C 
C     Input data files:
C
C     a) parameter file (fixed format):
C        1. record: star name (max. 26 characters)
C        2. record: equatorial coordinates of the star (RA - hours, minutes, 
C                seconds; DEC - degrees, ',"), their equinox, numerical code
C                of the observatory (format 6F10.4,2I5). 
C
C     Geographical coordinates of the following observatories are built in
C     the program and can be chosen via the appropriate numerical code:
C
C     1... Ondrejov       2... DAO Victoria    3... Lick        4... Okayama
C     5... David Dunlap   6... Observatoire Haute Provence      7... Crimea
C     8... Zelenchuk 6-m  9... Kitt Peak National Observatory  10... CFHT
C     11.. Tautenburg     12.. WHT                             13... NOT
C     14.. LaSilla ESO 3.6m  70.7345W 29.2584S 2400m  @ MZ from ESO web
C     15.. Paranal ESO UVES  70°24′15″W  24°37′38″S 2635m  @ MZ from Wiki
C                          = -70.404167 -24.627222
C     16.. Calar Alto 3.5m 
C     17.. HARPN GEOELEV     17 53 20.6 W   28 45 14.4 N   2387.2m  @MZ from header
C                          =  -17.88905555  28.754000   2387.2m
C     18.. HPF GEOELEV      104 00 53 W     30 40 53 N     2026.m  @MZ from wiki
C                          = -104.014722  30.681444   2026.m
C
C     b) file with dates of exposure (free format):
C        records containing:
C          number of the spectrum, year, month, day, hours, minutes and seconds
C          for the beginning of the exposure and the exposure time in seconds.
C          Date and time is in UT. 
C
C============================================================================
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION OBS(18),DEL(18),SIR(18),VYS(18),VORBP(3),VROTP(3)
      DOUBLE PRECISION VALS(400),RRD(6),RR(3),SS(3),CVAL(400),ET2(2)
      DOUBLE PRECISION JD,JDBAR,DPREC(3,3),DP(3),DPP(3),DELTAT(538)
      DOUBLE PRECISION PVSUN(6),LMST
      CHARACTER*120 IN,IN2,OUT
      CHARACTER*26 STAR
      CHARACTER*6 NAMS(400),TTL(14,3),CNAM(400)
      INTEGER IPT(39),DENUM
      LOGICAL KM,BARY
      DOUBLE PRECISION MUALP,MUDEL,D_TIME,DS ! @MZ
C
      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/STCOMX/KM,BARY,PVSUN
      COMMON/CHRHDR/CNAM,TTL
C
      PRINT *,'Input file with star name and equatorial coordinates?'
      READ(5,'(A)') IN2
      PRINT *,'Input data file with dates of exposures?'
      READ(5,'(A)') IN
      PRINT *,'Output file of results?'
      READ(5,'(A)') OUT
      OPEN(1,FILE=IN,STATUS='OLD',BLANK='ZERO')
      OPEN(2,FILE=IN2,STATUS='OLD',BLANK='ZERO')
      OPEN(9,FILE=OUT,STATUS='NEW',BLANK='ZERO')
C
C     Reading input parameters
C
      READ(2,101) STAR
      READ(2,102) AL1,AL2,AL3,D1,D2,D3,LEQ,K,MUALP,MUDEL !@MZ Proper motions [mas/yr]
  101 FORMAT(A26)
  102 FORMAT(6F10.4,2I5,2F10.4)
      DEQ=LEQ
C
C
      WRITE(9,200)
  200 FORMAT('  BarCor                                          ',
     *       'Release  25 August 2009') 
C
C     The observatory choice
C
      DATA OBS/8HONDREJOV,8HVICTORIA,8H  LICK  ,8H OKAYAMA,8H DUNLAP ,
     1         8H   OHP  ,8H  KRYM  ,8HZELENCUK,8H  KPNO  ,8H  CFHT  ,
     2         8H  TAUTEN,8H  WHT   ,8H     NOT,
     3         8H ESO3p6m,8H ESOUVES,8HCAHA3p5m,8H HARPN  ,8H HPF    /
      DATA DEL/4.106481481D-2,-3.428240741D-1,-3.379050926D-1,         ! DEL = LONGITUDE/360
     1         3.711018519D-1,-2.206168981D-1,1.5972222D-2,9.4444444D-2,
     2         1.151041600D-1,-3.097222222D-1,-4.3186728395D-1,
     3         3.2531019D-2,-4.9672222D-2,-4.968079D-2,
     4         -1.9648472222D-1,-1.95567130555556D-1,-0.0707305555D-1, !@MZ -70.7345/360, -70.404167/360,  -2.5463/360
     5         -4.9691820972222D-2, -0.288929783333333/                                    ! -17.88905555/360
      DATA SIR/8.712053368D-1,8.468531456D-1,6.517107908D-1,           ! SIR = LATIDUTE/360
     1         6.034281963D-1,7.65549891D-1,7.66781317D-1,7.80743951D-1,
     2         7.618943961D-1,5.57632695D-1,3.460309168D-1,8.8977079D-1,
     3         5.0171337D-1,5.0190918D-1,
     4         -5.10655414D-1,-4.29826109519583D-1,6.4967437944536D-1, !@MZ [rad] -29.2584/180*pi, -24.627222/180*pi, 37.2236/180*pi
     5         5.0185197311845D-1, 5.35492217066259D-1/                                     !          28.754000/180*pi
      DATA VYS/527.0D0,229.0D0,1290.0D0,372.0D0,245.0D0,684.0D0,650.0D0,
     1         2070.0D0,2120.0D0,4215.0D0,341.0D0,2332.0D0,2382.0D0,
     2         2400.0D0, 2635.0D0, 2168.0D0, 2387.2D0, 2026.0D0/ !@MZ
      DATA CD,CS/2.6179938779914943333D-1,1.745329251994329556D-2/ ! @MZ CD = pi/12=2pi/24;  CS = pi/180
C
C     DATA TT = Terrestrial - UT1 Time
C
C       1620-1754
      DATA DELTAT/ 124,119,115,110,106,102,98,95,91,88,85,82,79,77,74,
     *    72,70,67,65,63,62,60,58,57,55,54,53,51,50,49,48,47,46,45,44,
     *    43,42,41,40,38,37,36,35,34,33,32,31,30,28,27,26,25,24,23,22,
     *    21,20,19,18,17,16,15,14,14,13,12,12,11,11,10,10,9,9,9,9,
     *    9,9,9,9,9,10,9,9,9,9,9,9,9,10,10,10,10,10,10,10,
     *    10,10,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,12,
     *    12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,14,14,14,14,
C       1755-1899
     *    14,14,14,15,15,15,15,15,15,15,16,16,16,16,16,16,16,16,16,16,
     *    17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,16,16,16,
     *    16,15,15,14,14,13.7,13.4,13.1,12.9,12.7,12.6,12.5,12.5,12.5,
     *    12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.5,12.4,12.3,12.2,12.0,
     *    11.7,11.4,11.1,10.6,10.2,9.6,9.1,8.6,8.0,7.5,7.0,6.6,6.3,6.0,
     *    5.8,5.7,5.6,5.6,5.6,5.7,5.8,5.9,6.1,6.2,6.3,6.5,6.6,6.8,6.9,
     *    7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.7,7.8,7.8,7.88,7.82,7.54,6.97,
     *    6.40,6.02,5.41,4.10,2.92,1.81,1.61,0.10,-1.02,-1.28,-2.69,
     *    -3.24,-3.64,-4.54,-4.71,-5.11,-5.40,-5.42,-5.20,-5.46,-5.46,
     *    -5.79,-5.63,-5.64,-5.80,-5.66,-5.87,-6.01,-6.19,-6.64,-6.44,
     *    -6.47,-6.09,-5.76,-4.66,-3.74,
C       1900-1961
     *    -2.72,-1.54,-0.02,1.24,2.64,3.86,5.37,6.14,7.75,9.13,10.46,
     *    11.53,13.36,14.65,16.01,17.20,18.24,19.06,20.25,20.95,21.16,
     *    22.25,22.41,23.03,23.49,23.62,23.86,24.49,24.34,24.08,24.02,
     *    24.00,23.87,23.95,23.86,23.93,23.73,23.92,23.96,24.02,24.33,
     *    24.83,25.30,25.70,26.24,26.77,27.28,27.78,28.25,28.71,29.15,
     *    29.57,29.97,30.36,30.72,31.07,31.35,31.68,32.18,32.68,33.15,
     *    33.59,
C       1962-1971
     *          1.813,1.936,2.058,2.142,2.289,2.407,2.551,2.659,
     *          2.847,3.043,3.217,3.350,3.558,3.761,3.964,4.139,
     *          4.360,4.586,4.817,5.006,5.248,5.476,5.700,5.871,
     *          6.111,6.344,6.573,6.775,7.021,7.274,7.521,7.734,
     *          7.997,8.271,8.526,8.722,8.985,9.237,9.503,9.741,
C       1972-1981
     *     10.045,10.347,10.638,10.888,11.189,11.489,11.771,12.016,
     *     12.301,12.553,12.814,13.023,13.292,13.554,13.799,13.999,
     *     14.274,14.545,14.812,15.054,15.336,15.595,15.850,16.062,
     *     16.351,16.651,16.917,17.123,17.402,17.670,17.918,18.112,
     *     18.355,18.582,18.792,18.970,19.196,19.414,19.629,19.777,
C       1982-1991
     * 19.983,20.184,20.391,20.550,20.773,21.036,21.250,21.400,
     * 21.6042,21.7603,21.9018,22.0074,22.1587,22.3058,22.4516,22.5333,
     * 22.6872,22.8157,22.9292,23.0058,23.1382,23.2788,23.3972,23.4816,
     * 23.6356,23.7823,23.9100,23.9771,24.1159,24.2443,24.3857,24.4899,
     * 24.6713,24.8631,25.0386,25.1803,25.3813,25.5871,25.7735,25.9204,
C       1992-1999
     * 26.1252,26.3562,26.5569,26.7146,26.9378,27.1734,27.4010,27.5748,
     * 27.8004,28.0202,28.2171,28.3738,28.6014,28.8437,29.0614,29.2196,
     * 29.4447,29.6292,29.8129,29.9362,30.1110,30.2914,30.4731,30.6087,
     * 30.7818,30.9622,31.1004,31.1582,31.2833,31.3840,31.4801,31.4801,
C       2000-2010(January-December)
     * 31.636,31.636,31.796,31.796,31.906,31.906,32.016,32.016,
     * 32.116,32.116,32.226,32.226,32.286,32.286,32.366,32.366,
     * 32.386,32.386,32.416,32.416,32.436,32.600,32.600,32.600,
     * 32.700,32.800,32.900,33.000,33.000,33.000,33.000,33.000,
     * 33.000,33.000,33.000,33.000,34.000,34.000,34.000,34.000,
     * 34.000,34.000,34.000,34.000/
C
      IF(K.EQ.0) THEN
         !MZ: read obs params from file instead hardcoding
         READ(2,*) IN
         READ(2,*) SIRKA, DELKA, ELEV
         DELKA = DELKA / 360
         SIRKA = SIRKA * (3.1415926535897932D0/180)
         HVEZD = TRANSFER (IN, HVEZD)
      ELSE
         HVEZD=OBS(K)
         DELKA=DEL(K)
         SIRKA=SIR(K)
         ELEV=VYS(K)
      ENDIF
C
  202 FORMAT(10X)
  233 WRITE(9,203) STAR
  203 FORMAT(2X,'List of observations of star ',A26)
      L1=AL1
      L2=AL2
      JD1=D1
      JD2=D2
      JD3=D3
      WRITE(9,204) HVEZD,LEQ,L1,L2,AL3,LEQ,JD1,JD2,JD3
  204 FORMAT(2X,A8,22X,3HRA(,I4,1H),2I3,F5.1,3X,3HDA(,I4,1H),I4,2I3)
      WRITE(9,202)
      WRITE(9,205)
  205 FORMAT(72(1H=))
      WRITE(9,206)
  206 FORMAT(1X,'N.    Date & UT start    exp[s]         ',
     *      '   BJD                RVcorr')
      WRITE(9,205)
      WRITE(9,202)
C
   11 READ(1,*,END=10) NUM,LP,M,D,A1,A2,A3,A4
      J=LP-1900
      D2=SIGN(D2,D1)
      D3=SIGN(D3,D1)
      EXP=A4
      DEN=0 
      DEN=D+((A1+(A2/60)+(A3/3600))/24)
      D=DEN+A4/172800
C
C     Geocentric Julian Date 
C
      JD=GEO(D,M,LP)
C
C     Change of RA and DEC to radians
C
      RA=CD*(AL1+AL2/60.0+AL3/3600.0)
      DEC=CS*(D1+D2/60.0+D3/3600.0)
C
      L=LP-1619
      IF (L .GT. 342) THEN
         L=(L-343)*4+343
         IF (M .GE. 1) TT=DELTAT(L)
         IF (M .GE. 4) TT=DELTAT(L+1)
         IF (M .GE. 7) TT=DELTAT(L+2)
         IF (M .GE. 10) TT=DELTAT(L+3)
      ELSE 
         TT=DELTAT(L)
      END IF  
C
C     Computation of the Earth coordinates and velocities (AU, AU/day)
C
  666 NTARG=3
      NCTR=12
      ET=JD+((TT+32.184D0)/86400.0D0)
      ET2(1)=ET
      ET2(2)=0.0D0
      CALL  Const (NAMS, VALS, SS, NVS)
      CALL  PLEPH (ET,NTARG,NCTR,RRD)
C     Add proper motion to coordinates @MZ
      D_TIME = (JD-GEO(1.0D0,1,int(DEQ))) / 365.25D0   ! @MZ d_time=o_time-1991.25     !yeardiff between obs and HIP
      DS = CS * MUALP * 1.D-3 / 3600.D0  ! @MZ ds=mualp*1.d-3/3600.d0
      DS = DS / DCOS(DEC)                ! @MZ ds=ds/dcos(dein)
      DS = DS * D_TIME                   ! @MZ ds=ds*d_time
      RA = RA + DS                       ! @MZ rightasc=rightasc+ds
      DS = CS * MUDEL * 1.D-3 / 3600.D0  ! @MZ ds=mudel*1.d-3/3600.d0
      DS = DS * D_TIME                   ! @MZ ds=ds*d_time
      DEC = DEC + DS                     ! @MZ decln=decln+ds
      WRITE(6,*) D_TIME,LP,M,RA,DEC-DS, DEC
C
      PI=3.1415926535897932D0
      DP(1)=DCOS(RA)*DCOS(DEC)
      DP(2)=DSIN(RA)*DCOS(DEC)
      DP(3)=DSIN(DEC)
      U=149597870.691D0/86400.0D0
C
C     Local mean sidereal time     
C
      DOD=((A1+(A2/60.0D0)+(A3/3600.0D0))/24.0D0)+A4/172800.0D0
      DLONG=1.0D0-DELKA 
      TU=(JD-2451545.0D0)/36525.0D0
      T=(ET-2451545.0D0)/36525.0D0
      GMST=DOD+(24110.5493771D0+8640184.79447825D0*TU+
     *     307.4771013D0*(T-TU)+0.092772110D0*T**2-
     *     0.0000002926D0*T**3-0.00000199708D0*T**4-
     *     0.000000002454D0*T**5)/3600.0D0/24.0D0
      LMST=GMST-DLONG
      LMST=(LMST-DFLOAT(IDINT(LMST)))*24.0D0
      IF ( LMST .LT. 0.0D0 )  LMST=LMST+24.0D0
      IF ( LMST .GT. 24.0D0 ) LMST=LMST-24.0D0   
      S=LMST*PI/12.0D0
C
C     Earth rotational velocity
C
      TU=(JD-2451545.0D0-DOD)/36525.0D0
      SDD=1.002737909350795D0+5.9006D-11*TU-5.9D-15*TU**2
      FLAT=1.0D0/298.257223563D0
      EFLAT=2*FLAT-FLAT**2
      VRO=ELEV+(6378137.0D0/DSQRT(1.0D0-EFLAT*DSIN(SIRKA)**2))
      VROT=2.0D0*PI*SDD*DCOS(SIRKA)*VRO/24.0D0/3600.0D0/1000.0D0
      VROTP(1)=-VROT*DSIN(S)
      VROTP(2)=VROT*DCOS(S)
      VROTP(3)=0.0D0
C
      IF (DEQ .EQ. 2000) THEN
         DEKV=2451545.0D0
      ELSE IF (DEQ .EQ. 1950) THEN
         DEKV=2433282.42345905D0
      ELSE IF (DEQ .EQ. 1900) THEN 
         DEKV=2415020.31352D0
      END IF
C
C     Barycentric correction of time
C
      C1=499.004782D0/3600.0D0/24.0D0
      DT=(2451545.0D0-DEKV)/36525.0D0
      CALL PREC(DT,DPREC)
      DO 6 N=1,3
         DPP(N)=DP(1)*DPREC(N,1)+DP(2)*DPREC(N,2)+DP(3)*DPREC(N,3)
    6 CONTINUE
      DBTIM=C1*(RRD(1)*DPP(1)+RRD(2)*DPP(2)+RRD(3)*DPP(3))
      JDBAR=JD-2400000D0+DBTIM
C
C     Barycentric correction of radial velocity    
C
      DT=(ET-DEKV)/36525.0D0 
      CALL PREC(DT,DPREC)
      DO 4 N=1,3
         DPP(N)=DP(1)*DPREC(N,1)+DP(2)*DPREC(N,2)+DP(3)*DPREC(N,3)
    4 CONTINUE
C
      DT=(ET-2451545.0D0)/36525.0D0 
      CALL PREC(DT,DPREC)
      DO 5 N=1,3
         RR(N)=(RRD(4)*DPREC(N,1)+RRD(5)*DPREC(N,2)+
     1         RRD(6)*DPREC(N,3))*U
         VORBP(N)=RR(N)+VROTP(N)
    5 CONTINUE
      VRBCOR=DPP(1)*VORBP(1)+DPP(2)*VORBP(2)+DPP(3)*VORBP(3)
C
C     Writing output file
C
      LD=DEN
      LH=A1
      LM=A2
      LS=A3
      LE=A4
      WRITE(9,201) NUM,LP,M,LD,LH,LM,LS,LE,JDBAR,VRBCOR
  201 FORMAT(I5,I5,5I3,I5,10X,F16.7,5X,F11.6)
      GO TO 11
   10 STOP 17
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION GEO(D,M,LP)
C
C     To compute geocentric JD, written by Vondrak (2001)
C
      REAL*8 D
      i=LP
      j=m
      if(m.le.2) then
      i=i-1
      j=j+12
      endif
      ia=i/100
      ib=2-ia+ia/4
      if(LP.ge.1583) goto 1
      if(LP.le.1581) then
      ib=0
      goto 1
      endif
      if(m.lt.10.or.(m.eq.10.and.d.le.5.)) ib=0
    1 GEO=d+ib+int(365.25D0*(i+4716))+int(30.6001D0*(j+1))-1524.5D0
      RETURN
      END
C
C=====================================================================
C
      SUBROUTINE PREC(DT,DPREC)
C
C     To compute the general precession after Capitaine (2003). 
C
      IMPLICIT DOUBLE PRECISION(D)
      DOUBLE PRECISION DPREC(3,3)
      DATA DRAD/206264.8062470963D0/
C
      DZETA=(2.650545D0+2306.083227D0*DT+0.2988499D0*DT**2+
     *        0.01801828D0*DT**3-0.000005971D0*DT**4-
     *        0.0000003173D0*DT**5)/DRAD
      DZ=(-2.650545D0+2306.077181D0*DT+1.0927348D0*DT**2+
     *        0.01826837D0*DT**3-0.000028596D0*DT**4-
     *        0.0000002904D0*DT**5)/DRAD
      DTHET=(2004.191903D0*DT-0.4294934D0*DT**2-0.04182264D0*DT**3-
     *        0.000007089D0*DT**4-0.0000001274D0*DT**5)/DRAD
C
      DSZETA=DSIN(DZETA)
      DCZETA=DCOS(DZETA)
      DSZ=DSIN(DZ)
      DCZ=DCOS(DZ)
      DSTHET=DSIN(DTHET)
      DCTHET=DCOS(DTHET)
C
      DPREC(1,1)=-DSZETA*DSZ+DCZETA*DCZ*DCTHET
      DPREC(1,2)=-DCZETA*DSZ-DSZETA*DCZ*DCTHET
      DPREC(1,3)=-DCZ*DSTHET
      DPREC(2,1)=DSZETA*DCZ+DCZETA*DSZ*DCTHET
      DPREC(2,2)=DCZETA*DCZ-DSZETA*DSZ*DCTHET
      DPREC(2,3)=-DSZ*DSTHET
      DPREC(3,1)=DCZETA*DSTHET
      DPREC(3,2)=-DSZETA*DSTHET
      DPREC(3,3)=DCTHET
      RETURN   
      END      
C
C======================================================================
C     JPL PLANETARY AND LUNAR EPHEMERIS
C======================================================================
c
      SUBROUTINE FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)
C
C  To set the values of NRECL, KSIZE, NRFILE and NAMFIL.
C
      SAVE
      CHARACTER*80 NAMFIL
C
C  NRECL=1 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN S.P. WORDS
C  NRECL=4 IF "RECL" IN THE OPEN STATEMENT IS THE RECORD LENGTH IN BYTES
C
      NRECL=4
      NRFILE=12
      NAMFIL="JPLEPH"
      KSIZE = 2036
      RETURN
      END
C======================================================================
C
      SUBROUTINE PLEPH(ET,NTARG,NCENT,RRD)
C
C     To read the JPL Planetary Ephemeris and give the position and 
C     velocity of the point 'NTARG' with respect to 'NCENT'.
C
C     ET = D.P. JULIAN EPHEMERIS DATE AT WHICH INTERPOLATION IS WANTED.
C     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
C     NCENT = INTEGER NUMBER OF CENTER POINT.
C
C            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
C
C                1 = MERCURY           8 = NEPTUNE
C                2 = VENUS             9 = PLUTO
C                3 = EARTH            10 = MOON
C                4 = MARS             11 = SUN
C                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
C                6 = SATURN           13 = EARTH-MOON BARYCENTER
C                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
C                            15 = LIBRATIONS, IF ON EPH FILE
C
C             (IF NUTATIONS ARE WANTED, SET NTARG = 14. FOR LIBRATIONS,
C              SET NTARG = 15. SET NCENT=0.)
C
C      RRD = OUTPUT 6-WORD D.P. ARRAY CONTAINING POSITION AND VELOCITY
C            OF POINT 'NTARG' RELATIVE TO 'NCENT'. THE UNITS ARE AU AND
C            AU/DAY. FOR LIBRATIONS THE UNITS ARE RADIANS AND RADIANS
C            PER DAY. IN THE CASE OF NUTATIONS THE FIRST FOUR WORDS OF
C            RRD WILL BE SET TO NUTATIONS AND RATES, HAVING UNITS OF
C            RADIANS AND RADIANS/DAY.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RRD(6),ET2Z(2),ET2(2),PV(6,13)
      DIMENSION SS(3),CVAL(400),PVSUN(6),zips(2)
      data zips/2*0.d0/
      LOGICAL BSAVE,KM,BARY
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      INTEGER LIST(12),IPT(39),DENUM
C
      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/STCOMX/KM,BARY,PVSUN
C
C     INITIALIZE ET2 FOR 'STATE' AND SET UP COMPONENT COUNT
C
      ET2(1)=ET
      ET2(2)=0.D0
      GO TO 11
C
C     ENTRY POINT 'DPLEPH' FOR DOUBLY-DIMENSIONED TIME ARGUMENT
C          (SEE THE DISCUSSION IN THE SUBROUTINE STATE)

      ENTRY DPLEPH(ET2Z,NTARG,NCENT,RRD)
      ET2(1)=ET2Z(1)
      ET2(2)=ET2Z(2)
  11  IF(FIRST) CALL STATE(zips,list,pv,pnut)
      FIRST=.FALSE.
  96  IF(NTARG .EQ. NCENT) RETURN
      DO I=1,12
      LIST(I)=0
      ENDDO
C
C     CHECK FOR NUTATION CALL
C
      IF(NTARG.NE.14) GO TO 97
        IF(IPT(35).GT.0) THEN
          LIST(11)=2
          CALL STATE(ET2,LIST,PV,RRD)
          RETURN
        ELSE
          do i=1,4
          rrd(i)=0.d0
          enddo
          WRITE(6,297)
  297     FORMAT(' *****  NO NUTATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF
C
C     CHECK FOR LIBRATIONS
C
  97  do i=1,6
      rrd(i)=0.d0
      enddo
      IF(NTARG.NE.15) GO TO 98
        IF(IPT(38).GT.0) THEN
          LIST(12)=2
          CALL STATE(ET2,LIST,PV,RRD)
          DO I=1,6
          RRD(I)=PV(I,11)
          ENDDO
          RETURN
        ELSE
          WRITE(6,298)
  298     FORMAT(' *****  NO LIBRATIONS ON THE EPHEMERIS FILE  *****')
          STOP
        ENDIF
C
C       FORCE BARYCENTRIC OUTPUT BY 'STATE'
C
  98  BSAVE=BARY
      BARY=.TRUE.
C
C       SET UP PROPER ENTRIES IN 'LIST' ARRAY FOR STATE CALL
C
      DO I=1,2
      K=NTARG
      IF(I .EQ. 2) K=NCENT
      IF(K .LE. 10) LIST(K)=2
      IF(K .EQ. 10) LIST(3)=2
      IF(K .EQ. 3) LIST(10)=2
      IF(K .EQ. 13) LIST(3)=2
      ENDDO
C
C       MAKE CALL TO STATE
C
      CALL STATE(ET2,LIST,PV,RRD)
      IF(NTARG .EQ. 11 .OR. NCENT .EQ. 11) THEN
      DO I=1,6
      PV(I,11)=PVSUN(I)
      ENDDO
      ENDIF
C
      IF(NTARG .EQ. 12 .OR. NCENT .EQ. 12) THEN
      DO I=1,6
      PV(I,12)=0.D0
      ENDDO
      ENDIF
C
      IF(NTARG .EQ. 13 .OR. NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,13)=PV(I,3)
      ENDDO
      ENDIF
C
      IF(NTARG*NCENT .EQ. 30 .AND. NTARG+NCENT .EQ. 13) THEN
      DO I=1,6
      PV(I,3)=0.D0
      ENDDO
      GO TO 99
      ENDIF
C
      IF(LIST(3) .EQ. 2) THEN
      DO I=1,6
      PV(I,3)=PV(I,3)-PV(I,10)/(1.D0+EMRAT)
      ENDDO
      ENDIF
C
      IF(LIST(10) .EQ. 2) THEN
      DO I=1,6
      PV(I,10)=PV(I,3)+PV(I,10)
      ENDDO
      ENDIF
C
  99  DO I=1,6
      RRD(I)=PV(I,NTARG)-PV(I,NCENT)
      ENDDO
      BARY=BSAVE
      RETURN
      END
C======================================================================
C
      SUBROUTINE INTERP(BUF,T,NCF,NCM,NA,IFL,PV)
C
C     To differentiate and interpolate a set of Chebyshev coefficients to give 
C     position and velocity 
C
C     INPUT:
C
C         BUF   1ST LOCATION OF ARRAY OF D.P. CHEBYSHEV COEFFICIENTS OF POSITION
C           T   T(1) IS DP FRACTIONAL TIME IN INTERVAL COVERED BY
C               COEFFICIENTS AT WHICH INTERPOLATION IS WANTED
C               (0 .LE. T(1) .LE. 1).  T(2) IS DP LENGTH OF WHOLE
C               INTERVAL IN INPUT TIME UNITS.
C         NCF   # OF COEFFICIENTS PER COMPONENT
C         NCM   # OF COMPONENTS PER SET OF COEFFICIENTS
C          NA   # OF SETS OF COEFFICIENTS IN FULL ARRAY
C               (I.E., # OF SUB-INTERVALS IN FULL INTERVAL)
C          IFL  INTEGER FLAG: =1 FOR POSITIONS ONLY
C                             =2 FOR POS AND VEL
C
C     OUTPUT:
C
C         PV   INTERPOLATED QUANTITIES REQUESTED.  DIMENSION
C               EXPECTED IS PV(NCM,IFL), DP.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DOUBLE PRECISION BUF(NCF,NCM,*),T(2),PV(NCM,*),PC(18),VC(18)
C
      DATA NP/2/
      DATA NV/3/
      DATA TWOT/0.D0/
      DATA PC(1),PC(2)/1.D0,0.D0/
      DATA VC(2)/1.D0/
C
C       ENTRY POINT. GET CORRECT SUB-INTERVAL NUMBER FOR THIS SET
C       OF COEFFICIENTS AND THEN GET NORMALIZED CHEBYSHEV TIME
C       WITHIN THAT SUBINTERVAL.
C
      DNA=DBLE(NA)
      DT1=DINT(T(1))
      TEMP=DNA*T(1)
      L=IDINT(TEMP-DT1)+1
C
C         TC IS THE NORMALIZED CHEBYSHEV TIME (-1 .LE. TC .LE. 1)
C
      TC=2.D0*(DMOD(TEMP,1.D0)+DT1)-1.D0
C
C       CHECK TO SEE WHETHER CHEBYSHEV TIME HAS CHANGED,
C       AND COMPUTE NEW POLYNOMIAL VALUES IF IT HAS.
C       (THE ELEMENT PC(2) IS THE VALUE OF T1(TC) AND HENCE
C       CONTAINS THE VALUE OF TC ON THE PREVIOUS CALL.)
C
      IF(TC.NE.PC(2)) THEN
        NP=2
        NV=3
        PC(2)=TC
        TWOT=TC+TC
      ENDIF
C
C       BE SURE THAT AT LEAST 'NCF' POLYNOMIALS HAVE BEEN EVALUATED
C       AND ARE STORED IN THE ARRAY 'PC'.
C
      IF(NP.LT.NCF) THEN
        DO 1 I=NP+1,NCF
        PC(I)=TWOT*PC(I-1)-PC(I-2)
    1   CONTINUE
        NP=NCF
      ENDIF
C
C       INTERPOLATE TO GET POSITION FOR EACH COMPONENT
C
      DO 2 I=1,NCM
      PV(I,1)=0.D0
      DO 3 J=NCF,1,-1
      PV(I,1)=PV(I,1)+PC(J)*BUF(J,I,L)
    3 CONTINUE
    2 CONTINUE
      IF(IFL.LE.1) RETURN
C
C       IF VELOCITY INTERPOLATION IS WANTED, BE SURE ENOUGH
C       DERIVATIVE POLYNOMIALS HAVE BEEN GENERATED AND STORED.
C
      VFAC=(DNA+DNA)/T(2)
      VC(3)=TWOT+TWOT
      IF(NV.LT.NCF) THEN
        DO 4 I=NV+1,NCF
        VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)
    4   CONTINUE
        NV=NCF
      ENDIF
C
C       INTERPOLATE TO GET VELOCITY FOR EACH COMPONENT
C
      DO 5 I=1,NCM
      PV(I,2)=0.D0
      DO 6 J=NCF,2,-1
      PV(I,2)=PV(I,2)+VC(J)*BUF(J,I,L)
    6 CONTINUE
      PV(I,2)=PV(I,2)*VFAC
    5 CONTINUE
      RETURN
      END
C======================================================================
C
      SUBROUTINE SPLIT(TT,FR)
C
C     To break a D.P. number into a D.P. INTEGER and a D.P. fractional part.     
C
C       TT = D.P. INPUT NUMBER
C       FR = D.P. 2-WORD OUTPUT ARRAY.
C            FR(1) CONTAINS INTEGER PART
C            FR(2) CONTAINS FRACTIONAL PART
C
C            FOR NEGATIVE INPUT NUMBERS, FR(1) CONTAINS THE NEXT
C            MORE NEGATIVE INTEGER; FR(2) CONTAINS A POSITIVE FRACTION.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION FR(2)
C
C       MAIN ENTRY -- GET INTEGER AND FRACTIONAL PARTS
C
      FR(1)=DINT(TT)
      FR(2)=TT-FR(1)
      IF(TT.GE.0.D0 .OR. FR(2).EQ.0.D0) RETURN
C
C       MAKE ADJUSTMENTS FOR NEGATIVE INPUT NUMBER
C
      FR(1)=FR(1)-1.D0
      FR(2)=FR(2)+1.D0
      RETURN
      END
C======================================================================
C
      SUBROUTINE STATE(ET2,LIST,PV,PNUT)
C
C     To read and interpolate the JPL Planetary Ephemeris file. 
C
C     INPUT:
C
C         ET2   DP 2-WORD JULIAN EPHEMERIS EPOCH AT WHICH INTERPOLATION
C               IS WANTED.  ANY COMBINATION OF ET2(1)+ET2(2) WHICH FALLS
C               WITHIN THE TIME SPAN ON THE FILE IS A PERMISSIBLE EPOCH.
C
C                A. FOR EASE IN PROGRAMMING, THE USER MAY PUT THE
C                   ENTIRE EPOCH IN ET2(1) AND SET ET2(2)=0.
C
C                B. FOR MAXIMUM INTERPOLATION ACCURACY, SET ET2(1) =
C                   THE MOST RECENT MIDNIGHT AT OR BEFORE INTERPOLATION
C                   EPOCH AND SET ET2(2) = FRACTIONAL PART OF A DAY
C                   ELAPSED BETWEEN ET2(1) AND EPOCH.
C
C                C. AS AN ALTERNATIVE, IT MAY PROVE CONVENIENT TO SET
C                   ET2(1) = SOME FIXED EPOCH, SUCH AS START OF INTEGRATION,
C                   AND ET2(2) = ELAPSED INTERVAL BETWEEN THEN AND EPOCH.
C
C        LIST   12-WORD INTEGER ARRAY SPECIFYING WHAT INTERPOLATION
C               IS WANTED FOR EACH OF THE BODIES ON THE FILE.
C
C                         LIST(I)=0, NO INTERPOLATION FOR BODY I
C                                =1, POSITION ONLY
C                                =2, POSITION AND VELOCITY
C
C               THE DESIGNATION OF THE ASTRONOMICAL BODIES BY I IS:
C
C                         I = 1: MERCURY
C                           = 2: VENUS
C                           = 3: EARTH-MOON BARYCENTER
C                           = 4: MARS
C                           = 5: JUPITER
C                           = 6: SATURN
C                           = 7: URANUS
C                           = 8: NEPTUNE
C                           = 9: PLUTO
C                           =10: GEOCENTRIC MOON
C                           =11: NUTATIONS IN LONGITUDE AND OBLIQUITY
C                           =12: LUNAR LIBRATIONS (IF ON FILE)
C
C     OUTPUT:
C
C          PV   DP 6 X 11 ARRAY THAT WILL CONTAIN REQUESTED INTERPOLATED
C               QUANTITIES.  THE BODY SPECIFIED BY LIST(I) WILL HAVE ITS
C               STATE IN THE ARRAY STARTING AT PV(1,I).  (ON ANY GIVEN
C               CALL, ONLY THOSE WORDS IN 'PV' WHICH ARE AFFECTED BY THE
C               FIRST 10 'LIST' ENTRIES (AND BY LIST(12) IF LIBRATIONS ARE
C               ON THE FILE) ARE SET.  THE REST OF THE 'PV' ARRAY
C               IS UNTOUCHED.)  THE ORDER OF COMPONENTS STARTING IN
C               PV(1,I) IS: X,Y,Z,DX,DY,DZ.
C
C               ALL OUTPUT VECTORS ARE REFERENCED TO THE EARTH MEAN
C               EQUATOR AND EQUINOX OF J2000 IF THE DE NUMBER IS 200 OR
C               GREATER; OF B1950 IF THE DE NUMBER IS LESS THAN 200.
C
C               THE MOON STATE IS ALWAYS GEOCENTRIC; THE OTHER NINE STATES
C               ARE EITHER HELIOCENTRIC OR SOLAR-SYSTEM BARYCENTRIC,
C               DEPENDING ON THE SETTING OF COMMON FLAGS (SEE BELOW).
C
C               LUNAR LIBRATIONS, IF ON FILE, ARE PUT INTO PV(K,11) IF
C               LIST(12) IS 1 OR 2.
C
C         NUT   DP 4-WORD ARRAY THAT WILL CONTAIN NUTATIONS AND RATES,
C               DEPENDING ON THE SETTING OF LIST(11).  THE ORDER OF
C               QUANTITIES IN NUT IS:
C
C                        D PSI  (NUTATION IN LONGITUDE)
C                        D EPSILON (NUTATION IN OBLIQUITY)
C                        D PSI DOT
C                        D EPSILON DOT
C
C     COMMON AREA STCOMX:
C
C          KM   LOGICAL FLAG DEFINING PHYSICAL UNITS OF THE OUTPUT
C               STATES. KM = .TRUE., KM AND KM/SEC
C                          = .FALSE., AU AND AU/DAY
C               DEFAULT VALUE = .FALSE.  (KM DETERMINES TIME UNIT
C               FOR NUTATIONS AND LIBRATIONS.  ANGLE UNIT IS ALWAYS RADIANS.)
C
C        BARY   LOGICAL FLAG DEFINING OUTPUT CENTER.
C               ONLY THE 9 PLANETS ARE AFFECTED.
C                        BARY = .TRUE. =\ CENTER IS SOLAR-SYSTEM BARYCENTER
C                             = .FALSE. =\ CENTER IS SUN
C               DEFAULT VALUE = .FALSE.
C
C       PVSUN   DP 6-WORD ARRAY CONTAINING THE BARYCENTRIC POSITION AND
C               VELOCITY OF THE SUN.
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      DIMENSION ET2(2),PV(6,12),PNUT(4),T(2),PJD(4),BUF(1500),
     . SS(3),CVAL(400),PVSUN(6)
      INTEGER LIST(12),IPT(3,13)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      CHARACTER*6 TTL(14,3),CNAM(400)
      CHARACTER*80 NAMFIL
      CHARACTER*1000 PATH ! @MZ
      LOGICAL KM,BARY
C
      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,NUMDE,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
      COMMON/STCOMX/KM,BARY,PVSUN
C
C       ENTRY POINT - 1ST TIME IN, GET POINTER DATA, ETC., FROM EPH FILE
C
      IF(FIRST) THEN
        FIRST=.FALSE.
      CALL FSIZER3(NRECL,KSIZE,NRFILE,NAMFIL)
      IF(NRECL .EQ. 0) WRITE(*,*)'  ***** FSIZER IS NOT WORKING *****'
      IRECSZ=NRECL*KSIZE
      NCOEFFS=KSIZE/2
C
      CALL get_command_argument(0, path) ! @MZ
        OPEN(NRFILE,
     *       FILE=path(1:INDEX(path, '/', .TRUE.)) // NAMFIL,
     *       ACCESS='DIRECT',
     *       FORM='UNFORMATTED',
     *       RECL=IRECSZ,
     *       STATUS='OLD')
      READ(NRFILE,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,
     . ((IPT(I,J),I=1,3),J=1,12),NUMDE,(IPT(I,13),I=1,3)
      READ(NRFILE,REC=2)CVAL
      NRL=0
      ENDIF
C
C     MAIN ENTRY POINT
C
      IF(ET2(1) .EQ. 0.D0) RETURN
      S=ET2(1)-.5D0
      CALL SPLIT(S,PJD(1))
      CALL SPLIT(ET2(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)+.5D0
      PJD(2)=PJD(2)+PJD(4)
      CALL SPLIT(PJD(2),PJD(3))
      PJD(1)=PJD(1)+PJD(3)
C
C       ERROR RETURN FOR EPOCH OUT OF RANGE
C
      IF(PJD(1)+PJD(4).LT.SS(1) .OR. PJD(1)+PJD(4).GT.SS(2)) GO TO 98
C
C       CALCULATE RECORD # AND RELATIVE TIME IN INTERVAL
C
      NR=IDINT((PJD(1)-SS(1))/SS(3))+3
      IF(PJD(1).EQ.SS(2)) NR=NR-1
      T(1)=((PJD(1)-(DBLE(NR-3)*SS(3)+SS(1)))+PJD(4))/SS(3)
C
C       READ CORRECT RECORD IF NOT IN CORE
C
      IF(NR.NE.NRL) THEN
        NRL=NR
        READ(NRFILE,REC=NR,ERR=99)(BUF(K),K=1,NCOEFFS)
      ENDIF
      IF(KM) THEN
      T(2)=SS(3)*86400.D0
      AUFAC=1.D0
      ELSE
      T(2)=SS(3)
      AUFAC=1.D0/AU
      ENDIF
C
C   INTERPOLATE SSBARY SUN
C
      CALL INTERP(BUF(IPT(1,11)),T,IPT(2,11),3,IPT(3,11),2,PVSUN)
      DO I=1,6
      PVSUN(I)=PVSUN(I)*AUFAC
      ENDDO
C
C   CHECK AND INTERPOLATE WHICHEVER BODIES ARE REQUESTED
C
      DO 4 I=1,10
      IF(LIST(I).EQ.0) GO TO 4
      CALL INTERP(BUF(IPT(1,I)),T,IPT(2,I),3,IPT(3,I),
     & LIST(I),PV(1,I))
      DO J=1,6
       IF(I.LE.9 .AND. .NOT.BARY) THEN
       PV(J,I)=PV(J,I)*AUFAC-PVSUN(J)
       ELSE
       PV(J,I)=PV(J,I)*AUFAC
       ENDIF
      ENDDO
   4  CONTINUE
C
C       DO NUTATIONS IF REQUESTED (AND IF ON FILE)
C
      IF(LIST(11).GT.0 .AND. IPT(2,12).GT.0)
     * CALL INTERP(BUF(IPT(1,12)),T,IPT(2,12),2,IPT(3,12),
     * LIST(11),PNUT)
C
C       GET LIBRATIONS IF REQUESTED (AND IF ON FILE)
C
      IF(LIST(12).GT.0 .AND. IPT(2,13).GT.0)
     * CALL INTERP(BUF(IPT(1,13)),T,IPT(2,13),3,IPT(3,13),
     * LIST(12),PV(1,11))
      RETURN
C
  98  WRITE(*,198)ET2(1)+ET2(2),SS(1),SS(2)
 198  format(' ***  Requested JED,',f12.2,
     * ' not within ephemeris limits,',2f12.2,'  ***')
      stop
   99 WRITE(*,'(2F12.2,A80)')ET2,'ERROR RETURN IN STATE'
      STOP
      END
C======================================================================
C
      SUBROUTINE CONST(NAM,VAL,SSS,N)
C
C     To obtain the constants from the ephemeris file. 
C
C     OUTPUT:
C
C       NAM = CHARACTER*6 ARRAY OF CONSTANT NAMES
C       VAL = D.P. ARRAY OF VALUES OF CONSTANTS
C       SSS = D.P. JD START, JD STOP, STEP OF EPHEMERIS
C         N = INTEGER NUMBER OF ENTRIES IN 'NAM' AND 'VAL' ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      CHARACTER*6 NAM(*),TTL(14,3),CNAM(400)
      DOUBLE PRECISION VAL(*),SSS(3),SS(3),CVAL(400),zips(2)
      DOUBLE PRECISION xx(99)
      data zips/2*0.d0/
      INTEGER IPT(3,13),DENUM,list(11)
      logical first
      data first/.true./
C
      COMMON/EPHHDR/CVAL,SS,AU,EMRAT,DENUM,NCON,IPT
      COMMON/CHRHDR/CNAM,TTL
C
C  CALL STATE TO INITIALIZE THE EPHEMERIS AND READ IN THE CONSTANTS
C
      IF(FIRST) CALL STATE(zips,list,xx,xx)
      first=.false.
      N=NCON
      DO I=1,3
      SSS(I)=SS(I)
      ENDDO
      DO I=1,N
      NAM(I)=CNAM(I)
      VAL(I)=CVAL(I)
      ENDDO
      RETURN
      END
C======================================================================
