!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                          ~!
!~ Name:                                                                    ~!
!~    Buoyancy                                                              ~!
!~                                                                          ~!
!~ Description:                                                             ~!
!~    This function computes Convective Available Potential Energy (CAPE)   ~!
!~    with inhibition as a result of water loading given the data required  ~!
!~    to run up a sounding.                                                 ~!
!~                                                                          ~!
!~    Additionally, since we are running up a sounding anyways, this        ~!
!~    function returns the height of the Level of Free Convection (LFC) and ~!
!~    the pressure at the LFC.  That-a-ways, we don't have to run up a      ~!
!~    sounding later, saving a relatively computationally expensive         ~!
!~    routine.                                                              ~!
!~                                                                          ~!
!~ Usage:                                                                   ~!
!~    ostat = Buoyancy ( tK, rh, p, hgt, sfc, CAPE, ZLFC, PLFC, parcel )    ~!
!~                                                                          ~!
!~ Where:                                                                   ~!
!~                                                                          ~!
!~    IN                                                                    ~!
!~    --                                                                    ~!
!~    tK   = Temperature ( K )                                              ~!
!~    rh   = Relative Humidity ( % )                                        ~!
!~    p    = Pressure ( Pa )                                                ~!
!~    hgt  = Geopotential heights ( m )                                     ~!
!~    sfc  = integer rank within submitted arrays that represents the       ~!
!~           surface                                                        ~!
!~                                                                          ~!
!~    OUT                                                                   ~!
!~    ---                                                                   ~!
!~    ostat         INTEGER return status. Nonzero is bad.                  ~!
!~    CAPE ( J/kg ) Convective Available Potential Energy                   ~!
!~    ZLFC ( gpm )  Height at the LFC                                       ~!
!~    PLFC ( Pa )   Pressure at the LFC                                     ~!
!~                                                                          ~!
!~    tK, rh, p, and hgt are all REAL arrays, arranged from lower levels    ~!
!~    to higher levels.                                                     ~!
!~                                                                          ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
  FUNCTION Buoyancy ( nz, tk, rh, p, hgt, sfc, cape, cin, zlfc, plfc, lidx,  &
                      parcel ) result (ostat)

    IMPLICIT NONE

    INTEGER, INTENT ( IN )  :: nz          !~ Number of vertical levels
    INTEGER, INTENT ( IN )  :: sfc         !~ Surface level in the profile
    REAL,    INTENT ( IN )  :: tk   ( nz ) !~ Temperature profile ( K )
    REAL,    INTENT ( IN )  :: rh   ( nz ) !~ Relative Humidity profile ( % )
    REAL,    INTENT ( IN )  :: p    ( nz ) !~ Pressure profile ( Pa )
    REAL,    INTENT ( IN )  :: hgt  ( nz ) !~ Height profile ( gpm )
    REAL,    INTENT ( OUT ) :: cape        !~ CAPE ( J kg^-1 )
    REAL,    INTENT ( OUT ) :: cin         !~ CIN  ( J kg^-1 )
    REAL,    INTENT ( OUT ) :: zlfc        !~ LFC Height ( gpm )
    REAL,    INTENT ( OUT ) :: plfc        !~ LFC Pressure ( Pa )
    REAL,    INTENT ( OUT ) :: lidx        !~ Lifted index
    INTEGER                 :: ostat       !~ Function return status
                                           !~ Nonzero is bad.

    INTEGER, INTENT ( IN  ) :: parcel      !~ Most Unstable = 1 (default)
                                           !~ Mean layer    = 2
                                           !~ Surface based = 3

    !~ Derived profile variables
    !  -------------------------
    REAL                    :: ws   ( nz ) !~ Saturation mixing ratio
    REAL                    :: w    ( nz ) !~ Mixing ratio
    REAL                    :: etheta( nz )!~ Equivalent potential temperature. Modified by Zhixiao.
    REAL                    :: dTvK ( nz ) !~ Parcel / ambient Tv difference
    REAL                    :: buoy ( nz ) !~ Buoyancy
    REAL                    :: tlclK       !~ LCL temperature ( K )
    REAL                    :: plcl        !~ LCL pressure ( Pa )
    REAL                    :: pel         !~ Equilibrium pressure ( Pa ). Modified by Zhixiao.
    REAL                    :: nbuoy       !~ Negative buoyancy
    REAL                    :: pbuoy       !~ Positive buoyancy

    !~ Source parcel information
    !  -------------------------
    REAL                    :: srctK       !~ Source parcel temperature ( K )
    REAL                    :: srcrh       !~ Source parcel rh ( % )
    REAL                    :: srcws       !~ Source parcel sat. mixing ratio
    REAL                    :: srcw        !~ Source parcel mixing ratio
    REAL                    :: srcp        !~ Source parcel pressure ( Pa )
    REAL                    :: srctheta    !~ Source parcel theta ( K )
    REAL                    :: srcthetaeK  !~ Source parcel theta-e ( K )
    INTEGER                 :: srclev      !~ Level of the source parcel
    REAL                    :: spdiff      !~ Pressure difference
    REAL                    :: srce        !~ Equivalent potential temperature ( K ). Modified by Zhixiao.
 
    !~ Parcel variables
    !  ----------------
    REAL                    :: ptK        !~ Parcel temperature ( K )
    REAL                    :: ptvK       !~ Parcel virtual temperature ( K )
    REAL                    :: tvK        !~ Ambient virtual temperature ( K )
    REAL                    :: pw         !~ Parcel mixing ratio

    !~ Other utility variables
    !  -----------------------
    INTEGER                 :: i, j, k    !~ Dummy iterator
    INTEGER                 :: lfclev     !~ Level of LFC
    INTEGER                 :: ellev      !~ Level of EL. Modified by Zhixiao.
    INTEGER                 :: prcl       !~ Internal parcel type indicator
    INTEGER                 :: mlev       !~ Level for ML calculation
    INTEGER                 :: lyrcnt     !~ Number of layers in mean layer
    LOGICAL                 :: flag       !~ Dummy flag
    LOGICAL                 :: wflag      !~ Saturation flag
    REAL                    :: freeze     !~ Water loading multiplier
    REAL                    :: pdiff      !~ Pressure difference between levs 
    REAL                    :: pm, pu, pd !~ Middle, upper, lower pressures
    REAL                    :: lidxu      !~ Lifted index at upper level
    REAL                    :: lidxd      !~ Lifted index at lower level

    !~ Thermo / dynamical constants
    !  ----------------------------
    REAL                    :: Rd         !~ Dry gas constant
       PARAMETER ( Rd = 287.058 )         !~ J deg^-1 kg^-1
    REAL                    :: Cp         !~ Specific heat constant pressure
       PARAMETER ( Cp = 1004.67 )         !~ J deg^-1 kg^-1
    REAL                    :: g          !~ Acceleration due to gravity
       PARAMETER ( g  = 9.80665 )         !~ m s^-2
    REAL                    :: RUNDEF
       PARAMETER ( RUNDEF = -9.999E30 )
    !~ Initialize variables
    !  --------------------
    ostat  = 0
    CAPE   = REAL ( 0 )
    CIN    = RUNDEF !Change CIN filling values from 0 to default filling. CIN should not initially be filled by 0, because 0 means no inhibition energy. Modified by Zhixiao
    ZLFC   = RUNDEF
    PLFC   = RUNDEF

    !~ Look for submitted parcel definition
    !~ 1 = Most unstable
    !~ 2 = Mean layer
    !~ 3 = Surface based
    !  -------------------------------------
    IF ( parcel > 3 .or. parcel < 1 ) THEN
       prcl = 1
    ELSE
       prcl =  parcel
    END IF

    !~ Initalize our parcel to be (sort of) surface based.  Because of
    !~ issues we've been observing in the WRF model, specifically with
    !~ excessive surface moisture values at the surface, using a true
    !~ surface based parcel is resulting a more unstable environment
    !~ than is actually occuring.  To address this, our surface parcel
    !~ is now going to be defined as the parcel between 25-50 hPa
    !~ above the surface. UPDATE - now that this routine is in WRF,
    !~ going to trust surface info. GAC 20140415
    !  ----------------------------------------------------------------

    !~ Compute mixing ratio values for the layer
    !  -----------------------------------------
    DO k = sfc, nz
      ws  ( k )   = SaturationMixingRatio ( tK(k), p(k) )
      w   ( k )   = ( rh(k)/100.0 ) * ws ( k )
      !Removed by Zhixiao. Water vapor mixing ratio (w) is not conserved during parcel lifting processes. We should not use w to define MU layer.
      !thetav(k)   = Theta ( VirtualTemperature ( tK (k), w (k) ), p(k)/100.0 )
      !Added by Zhixiao. Critical modification: We use the model level with maximum equivalent potential temperature (etheta) below 500mb to define the MU layer
      !Because equivalent potential temperature is conserved in dry and moist adiabatic processes.
      etheta(k)   = Thetae( tK(k), p(k)/100.0, rh(k), w(k) )
    END DO

    srclev      = sfc
    srctK       = tK    ( sfc )
    srcrh       = rh    ( sfc )
    srcp        = p     ( sfc )
    srcws       = ws    ( sfc )
    srcw        = w     ( sfc )
    srctheta    = Theta ( tK(sfc), p(sfc)/100.0 )
    srce        = etheta (sfc) !Modified by Zhixiao
 
    !~ Compute the profile mixing ratio.  If the parcel is the MU parcel,
    !~ define our parcel to be the most unstable parcel within the lowest
    !~ 180 mb.
    !  -------------------------------------------------------------------
    mlev = sfc + 1
    !Change initial searching level from the second to first model level. Because we did not compute pdiff, and p(k-1) properties is unnecessary.
    !Modified by Zhixiao.
    DO k = sfc, nz
 
       !~ Identify the last layer within 100 hPa of the surface
       !  -----------------------------------------------------
       pdiff = ( p (sfc) - p (k) ) / REAL ( 100 )
       IF ( pdiff <= REAL (100) ) mlev = k

       !~ If we've made it past the lowest 500 hPa, exit the loop. MU layer is assumed below 500 hPa. Modified by Zhixiao.
       !  -------------------------------------------------------
       IF ( p(k) <= REAL (50000) ) EXIT

       IF ( prcl == 1 ) THEN
          ! Removed by Zhixiao, w can not used for defining MU layer
          !IF ( (p(k) > 70000.0) .and. (w(k) > srcw) ) THEN
          ! Modified by Zhixiao, MU layer is featured by the max etheta
          IF (etheta(k) > srce)  THEN !Modified by Zhixiao.
             srctheta = Theta ( tK(k), p(k)/100.0 )
             srcw = w ( k )
             srclev  = k
             srctK   = tK ( k )
             srcrh   = rh ( k )
             srcp    = p  ( k )
             srce = etheta(k) !Modified by Zhixiao
          END IF
       END IF
 
    END DO
 
    !~ If we want the mean layer parcel, compute the mean values in the
    !~ lowest 100 hPa.
    !  ----------------------------------------------------------------
    lyrcnt =  mlev - sfc + 1
    IF ( prcl == 2 ) THEN
 
       srclev   = sfc
       srctK    = SUM ( tK (sfc:mlev) ) / REAL ( lyrcnt )
       srcw     = SUM ( w  (sfc:mlev) ) / REAL ( lyrcnt )
       srcrh    = SUM ( rh (sfc:mlev) ) / REAL ( lyrcnt )
       srcp     = SUM ( p  (sfc:mlev) ) / REAL ( lyrcnt )
       srctheta = Theta ( srctK, srcp/100. )
 
    END IF
 
    srcthetaeK = Thetae ( srctK, srcp/100.0, srcrh, srcw )
    !~ Calculate temperature and pressure of the LCL
    !  ---------------------------------------------
    tlclK = TLCL ( tK(srclev), rh(srclev) )
    plcl  = p(srclev) * ( (tlclK/tK(srclev))**(Cp/Rd) )
    
    !~ Now lift the parcel
    !  -------------------
 
    buoy  = REAL ( 0 )
    pw    = srcw
    wflag = .false.
    DO k  = srclev, nz
       IF ( p (k) <= plcl ) THEN
 
          !~ The first level after we pass the LCL, we're still going to
          !~ lift the parcel dry adiabatically, as we haven't added the
          !~ the required code to switch between the dry adiabatic and moist
          !~ adiabatic cooling.  Since the dry version results in a greater
          !~ temperature loss, doing that for the first step so we don't over
          !~ guesstimate the instability.
          !  ----------------------------------------------------------------
 
          IF ( wflag ) THEN
             flag  = .false.
 
             !~ Above the LCL, our parcel is now undergoing moist adiabatic
             !~ cooling.  Because of the latent heating being undergone as
             !~ the parcel rises above the LFC, must iterative solve for the
             !~ parcel temperature using equivalant potential temperature,
             !~ which is conserved during both dry adiabatic and
             !~ pseudoadiabatic displacements.
             !  --------------------------------------------------------------
             ptK   = The2T ( srcthetaeK, p(k), flag )
 
             !~ Calculate the parcel mixing ratio, which is now changing
             !~ as we condense moisture out of the parcel, and is equivalent
             !~ to the saturation mixing ratio, since we are, in theory, at
             !~ saturation.
             !  ------------------------------------------------------------
             pw = SaturationMixingRatio ( ptK, p(k) )
 
             !~ Now we can calculate the virtual temperature of the parcel
             !~ and the surrounding environment to assess the buoyancy.
             !  ----------------------------------------------------------
             ptvK  = VirtualTemperature ( ptK, pw )
             tvK   = VirtualTemperature ( tK (k), w (k) )
 
             !~ Modification to account for water loading
             !  -----------------------------------------
             freeze = 0.033 * ( 263.15 - pTvK )
             IF ( freeze > 1.0 ) freeze = 1.0
             IF ( freeze < 0.0 ) freeze = 0.0
 
             !~ Approximate how much of the water vapor has condensed out
             !~ of the parcel at this level
             !  ---------------------------------------------------------
             freeze = freeze * 333700.0 * ( srcw - pw ) / 1005.7
 
             pTvK = pTvK - pTvK * ( srcw - pw ) + freeze
             dTvK ( k ) = ptvK - tvK
             buoy ( k ) = g * ( dTvK ( k ) / tvK )
 
          ELSE
 
             !~ Since the theta remains constant whilst undergoing dry
             !~ adiabatic processes, can back out the parcel temperature
             !~ from potential temperature below the LCL
             !  --------------------------------------------------------
             ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
 
             !~ Grab the parcel virtual temperture, can use the source
             !~ mixing ratio since we are undergoing dry adiabatic cooling
             !  ----------------------------------------------------------
             ptvK  = VirtualTemperature ( ptK, srcw )
 
             !~ Virtual temperature of the environment
             !  --------------------------------------
             tvK   = VirtualTemperature ( tK (k), w (k) )
 
             !~ Buoyancy at this level
             !  ----------------------
             dTvK ( k ) = ptvK - tvK
             buoy ( k ) = g * ( dtvK ( k ) / tvK )
 
             wflag = .true.
 
          END IF
 
       ELSE
 
          !~ Since the theta remains constant whilst undergoing dry
          !~ adiabatic processes, can back out the parcel temperature
          !~ from potential temperature below the LCL
          !  --------------------------------------------------------
          ptK   = srctheta / ( 100000.0/p(k) )**(Rd/Cp)
 
          !~ Grab the parcel virtual temperture, can use the source
          !~ mixing ratio since we are undergoing dry adiabatic cooling
          !  ----------------------------------------------------------
          ptvK  = VirtualTemperature ( ptK, srcw )
 
          !~ Virtual temperature of the environment
          !  --------------------------------------
          tvK   = VirtualTemperature ( tK (k), w (k) )
 
          !~ Buoyancy at this level
          !  ---------------------
          dTvK ( k ) = ptvK - tvK
          buoy ( k ) = g * ( dtvK ( k ) / tvK )
 
       END IF

       !~ Chirp
       !  -----
!          WRITE ( *,'(I15,6F15.3)' )k,p(k)/100.,ptK,pw*1000.,ptvK,tvK,buoy(k)
 
    END DO
 
    !~ Add up the buoyancies, find the LFC
    !  -----------------------------------
    flag   = .false.
    lfclev = -1
    ellev = -1 !Modified by Zhixiao
    DO k = sfc, nz !Modified by Zhixiao
       !~ LFC is defiend as the highest level when negative buyancy turns postive.
       !  -----------------------------------
       IF ( .not. flag .and. buoy (k) > REAL (0) .and. p (k) <= plcl ) THEN !Modified by Zhixiao
          flag = .true.
          lfclev = k
       END IF
       !~ Take the Highest EL as final result. Modified by Zhixiao
       !  ----------------------------------------------------------------
       IF (k >= 2) THEN !Modified by Zhixiao
          IF ( flag .and. buoy (k) < REAL (0) .and. buoy (k-1) >= REAL (0)) THEN !Modified by Zhixiao
             ellev = k !Modified by Zhixiao
          END IF
       END IF
       ! When buoy turns negative again, reset LFC flag and keep the highest LFC as the effective output
       IF (buoy (k) < REAL (0) .and. flag) THEN 
          flag = .false.
       END IF
    END DO
    IF ((ellev >= 0) .and. (lfclev >= 0)) THEN !Modified by Zhixiao
       plfc = p (lfclev)
       pel = p (ellev)
       CIN = REAL ( 0 )
       DO k = sfc+1, nz
          ! Make CAPE and CIN consistent with AMS definition
          ! https://glossary.ametsoc.org/wiki/Convective_available_potential_energy
          ! https://glossary.ametsoc.org/wiki/Convective_inhibition
          IF ( p (k) <= plcl .and. p (k) > plfc) THEN !Modified by Zhixiao
             ! CIN is the vertically integrated negative buoyant energy between LCL and LFC
             CIN  = CIN  + MIN ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
          END IF
          IF ( p (k) <= plfc .and. p (k) > pel) THEN !Modified by Zhixiao
             ! CAPE is the vertically integrated positive buoyant energy between LFC and EL
             CAPE = CAPE + MAX ( buoy (k), 0.0 ) * ( hgt (k) - hgt (k-1) )
          END IF !Modified by Zhixiao         
       END DO !Modified by Zhixiao
    END IF !Modified by Zhixiao
    !~ Calculate lifted index by interpolating difference between
    !~ parcel and ambient Tv to 500mb.
    !  ----------------------------------------------------------
    DO k = sfc + 1, nz

       pm = 50000.
       pu = p ( k )
       pd = p ( k - 1 )

       !~ If we're already above 500mb just set lifted index to 0.
       !~ --------------------------------------------------------
       IF ( pd .le. pm ) THEN
          lidx = 0.
          EXIT

       ELSEIF ( pu .le. pm .and. pd .gt. pm) THEN

          !~ Found trapping pressure: up, middle, down.
          !~ We are doing first order interpolation.  
          !  ------------------------------------------
          lidxu = -dTvK ( k ) * ( pu / 100000. ) ** (Rd/Cp)
          lidxd = -dTvK ( k-1 ) * ( pd / 100000. ) ** (Rd/Cp)
          lidx = ( lidxu * (pm-pd) + lidxd * (pu-pm) ) / (pu-pd)
          EXIT

       ENDIF

    END DO
    !~ Assuming the the LFC is at a pressure level for now
    !  ---------------------------------------------------
    IF ( lfclev > 0 ) THEN
       PLFC = p   ( lfclev )
       ZLFC = hgt ( lfclev )
    END IF

    IF ( PLFC /= PLFC .OR. PLFC < REAL (0) ) THEN
       PLFC = REAL ( -1 )
       ZLFC = REAL ( -1 )
    END IF
 
    IF ( CAPE /= CAPE ) cape = REAL ( 0 )

    IF ( CIN  /= CIN  ) cin  = RUNDEF

    !~ Chirp
    !  -----
!       WRITE ( *,* ) ' CAPE: ', cape, ' CIN:  ', cin
!       WRITE ( *,* ) ' LFC:  ', ZLFC, ' PLFC: ', PLFC
!       WRITE ( *,* ) ''
!       WRITE ( *,* ) ' Exiting buoyancy.'
!       WRITE ( *,* ) ' ==================================== '
!       WRITE ( *,* ) ''
 
END FUNCTION
