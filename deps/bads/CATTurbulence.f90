!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                          ~!
!~ Name:                                                                    ~!
!~    CATTurbulence                                                         ~!
!~                                                                          ~!
!~ Description:                                                             ~!
!~    This function calculates the turbulence index ( TI ) for one layer    ~!
!~    in the atmosphere given the horizontal wind components and the geo-   ~!
!~    potential height of two pressure levels.  The index is computed for   ~!
!~    the layer between the levels using the deformation and convergence    ~!
!~    of the wind field at the top and bottom of the layer and the vertical ~!
!~    wind shear is calculated within the layer.  The equation used for     ~!
!~    calculating TI is given by:                                           ~!
!~                                                                          ~!
!~                                                                          ~!
!~       TI = VWS * ( DEF + CONV )                                          ~!
!~                                                                          ~!
!~    Where:                                                                ~!
!~       VWS  = Vertical wind shear                                         ~!
!~       DEF  = Deformation                                                 ~!
!~       CONV = Convergence                                                 ~!
!~                                                                          ~!
!~ Notes:                                                                   ~!
!~                                                                          ~!
!~ References:                                                              ~!
!~    Ellrod, G.P. and D.J. Knapp, An objective clear-air turbulence        ~!
!~      forecasting technique: verification and operational use, Weather    ~!
!~      and Forecasting, 7, March 1992, pp. 150-165.                        ~!
!~                                                                          ~!
!~ Written:                                                                 ~!
!~    Scott Rentschler, Software Engineering Services                       ~!
!~    Fine Scale Models Team                                                ~!
!~    Air Force Weather Agency, 16WS/WXN                                    ~!
!~    DSN: 271-3331 Comm: (402) 294-3331                                    ~!
!~    scott.rentschler@offutt.af.mil                                        ~!
!~                                                                          ~!
!~ History:                                                                 ~!
!~    1 February 2008 ................... Scott Rentschler, (SES), 2WG/WEA  ~!
!~    INITIAL VERSION                                                       ~!
!~                                                                          ~!
!~    8 July 2009 ....................... Scott Rentschler, (SES), 2WG/WEA  ~!
!~    Adapted for new driver.                                               ~!
!~                                                                          ~!
!~    1 November 2012 ......................... Scott Rentschler, 16WS/WXN  ~!
!~    Modified to accept layer argument, which adds the flexibility to make ~!
!~    the computation for whichever flight level is desired.  Cleaned up    ~!
!~    some of the code and added a couple comments.                         ~!
!~                                                                          ~!
!~    28 August 2013 .................... Braedi Wickard, SEMS/NG/16WS/WXN  ~!
!~    Adapted for use within the Unified Post Processor. UPP can not handle ~!
!~    the layer argument for flight levels, so reverted to hardcoded levels ~!
!~                                                                          ~!
!~    25 April 2014 ............................. Glenn Creighton, 16WS/WXN ~!
!~    Adapted for use within WRF. WRF already computes many of these terms. ~!
!~    Stripped everything down to its bare bones to remove need to compute  ~!
!~    horizontal terms, now using deformation variables already within WRF. ~!
!~                                                                          ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION CATTurbulence ( ugrdbot, ugrdtop, vgrdbot, vgrdtop &
                         ,defor11bot, defor11top, defor12bot, defor12top &
                         ,defor22bot, defor22top, zbot, ztop ) result ( ti )

  IMPLICIT NONE

  !~ Variable declarations
  !  ---------------------
  REAL,    INTENT ( IN )  :: ugrdbot       !~ U-wind bottom of layer
  REAL,    INTENT ( IN )  :: ugrdtop       !~ U-wind top of layer
  REAL,    INTENT ( IN )  :: vgrdbot       !~ V-wind bottom of layer
  REAL,    INTENT ( IN )  :: vgrdtop       !~ V-wind top of layer
  REAL,    INTENT ( IN )  :: defor11bot    !~ 2*du/dx at bottom of layer
  REAL,    INTENT ( IN )  :: defor11top    !~ 2*du/dx at top of layer
  REAL,    INTENT ( IN )  :: defor12bot    !~ du/dy+dv/dx at bottom of layer
  REAL,    INTENT ( IN )  :: defor12top    !~ du/dy+dv/dx at top of layer
  REAL,    INTENT ( IN )  :: defor22bot    !~ 2*dv/dy at bottom of layer
  REAL,    INTENT ( IN )  :: defor22top    !~ 2*dv/dy at top of layer
  REAL,    INTENT ( IN )  :: zbot          !~ Height grid bottom
  REAL,    INTENT ( IN )  :: ztop          !~ Height grid top
  REAL                    :: ti            !~ Turbulence index

  !~ Function utility variables
  !  --------------------------
  REAL    :: dudx, dudx1, dudx2 !~ Wind differentials
  REAL    :: dvdy, dvdy1, dvdy2
  REAL    :: dudz, dvdz

  REAL    :: depth, vws, conv    !~ Depth, vertical wind shear, convergence
  REAL    :: def, shear, stretch !~ Deformation, shear, stretching terms

  !~ Initialize variables.
  !  ----------------------
  ti = REAL ( 0 )

  !~ Compute vertical wind shear
  !  ---------------------------
  depth = ABS ( ztop - zbot )
  dudz  = ( ugrdbot - ugrdtop ) / depth
  dvdz  = ( vgrdbot - vgrdtop ) / depth
  vws   = SQRT ( dudz**2 + dvdz**2  )

  dudx1 = defor11top / 2.
  dudx2 = defor11bot / 2.
  dudx  = ( dudx1 + dudx2 ) / REAL ( 2 )

  dvdy1 = defor22top / 2.
  dvdy2 = defor22bot / 2.
  dvdy  = ( dvdy1 + dvdy2 ) / REAL ( 2 )

  !~ Compute the deformation
  !  -----------------------
  stretch = dudx - dvdy
  shear   = ( defor12top + defor12bot ) / REAL ( 2 )
  def     = SQRT ( stretch**2 + shear**2 )

  !~ Compute the convergence
  !  -----------------------
  conv    = - ( dudx + dvdy )

  !~ Compute the turbulence index
  !  ----------------------------
  ti = vws * ( def + conv ) * 1.0E+07

  IF ( ti /= ti ) ti = REAL ( 0 )
  IF ( ti < 0   ) ti = REAL ( 0 )

END FUNCTION
