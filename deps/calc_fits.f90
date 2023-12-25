!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                          ~!
!~ Name:                                                                    ~!
!~    calc_fits                                                             ~!
!~                                                                          ~!
!~ Description:                                                             ~!
!~    This function computes Fighter Index Thermal Stress values given      ~!
!~    dry bulb temperature, relative humidity, and pressure.                ~!
!~                                                                          ~!
!~ Usage:                                                                   ~!
!~    fitsval = calc_fits ( p, tK, rh )                                     ~!
!~                                                                          ~!
!~ Where:                                                                   ~!
!~    p   = Pressure ( Pa )                                                 ~!
!~    tK  = Temperature ( K )                                               ~!
!~    rh  = Relative Humidity ( % )                                         ~!
!~                                                                          ~!
!~ Reference:                                                               ~!
!~    Stribley, R.F., S. Nunneley, 1978: Fighter Index of Thermal Stress:   ~!
!~    Development of interim guidance for hot-weather USAF operations.      ~!
!~    SAM-TR-78-6. Eqn. 9                                                   ~!
!~                                                                          ~!
!~    Formula:                                                              ~!
!~       FITS = 0.8281*Twb + 0.3549*Tdb + 5.08 (degrees Celsius)            ~!
!~                                                                          ~!
!~    Where:                                                                ~!
!~       Twb = Wet Bulb Temperature                                         ~!
!~       Tdb = Dry Bulb Temperature                                         ~!
!~                                                                          ~!
!~ Written:                                                                 ~!
!~    Scott Rentschler, Software Engineering Services                       ~!
!~    Fine Scale Models Team                                                ~!
!~    Air Force Weather Agency, 16WS/WXN                                    ~!
!~    DSN: 271-3331 Comm: (402) 294-3331                                    ~!
!~    scott.rentschler@offutt.af.mil                                        ~!
!~                                                                          ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_fits ( p, tK, rh ) RESULT ( fits )
 
  implicit none

  !~ Variable declaration
  !  --------------------
  real, intent ( in ) :: p               !~ Pressure ( Pa )
  real, intent ( in ) :: tK              !~ Temperature ( K )
  real, intent ( in ) :: rh              !~ Rel Humidity ( % )
  real                :: fits            !~ FITS index value
 
  !~ Utility variables
  !  --------------------------
  real                :: twb             !~ Wet bulb temperature ( C )
  real                :: wbt
 
  ! ---------------------------------------------------------------------- !
  ! ---------------------------------------------------------------------- !

  !~ Initialize variables
  !  --------------------
  fits = REAL ( 0 )

  !~ Get the wet bulb temperature in degrees Celsius
  !  -----------------------------------------------
  twb =  WetBulbTemp ( p, tK, rh ) - 273.15 

  !~ Compute the FITS index
  !  ----------------------
  fits = 0.8281*twb + 0.3549*( tK - 273.15 ) + 5.08
 
  !~ Convert the index to Kelvin
  !  ---------------------------
  fits = fits + 273.15

END FUNCTION
