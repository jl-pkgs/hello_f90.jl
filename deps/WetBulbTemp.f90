!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                          ~!
!~ Name:                                                                    ~!
!~    WetBulbTemp                                                           ~!
!~                                                                          ~!
!~ Description:                                                             ~!
!~    This function approximates the Wet Bulb Temperature (K) provided      ~!
!~    dry bulb temperature (K), relative humidity (%), and pressure (Pa).   ~!
!~                                                                          ~!
!~ Usage:                                                                   ~!
!~    wbt = WetBulbTemperature ( p, tK, rh )                                ~!
!~                                                                          ~!
!~ Where:                                                                   ~!
!~    p  = Pressure ( Pa )                                                  ~!
!~    tK = Temperature ( K )                                                ~!
!~    rh = Relative Humidity ( % )                                          ~!
!~                                                                          ~!
!~ Reference:                                                               ~!
!~    American Society of Civil Engineers                                   ~!
!~    Evapotraspiration and Irrigation Water Requirements                   ~!
!~    Jensen et al (1990) ASCE Manual No. 70, pp 176-177                    ~!
!~                                                                          ~!
!~ Written:                                                                 ~!
!~    Scott Rentschler, Software Engineering Services                       ~!
!~    Fine Scale Models Team                                                ~!
!~    Air Force Weather Agency                                              ~!
!~    DSM: 271-3331 Comm: (402) 294-3331                                    ~!
!~    scott.rentschler@offutt.af.mil                                        ~!
!~                                                                          ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION WetBulbTemp ( p, tK, rh) result( wbt )

  implicit none

  !~ Variable delclaration
  !  ---------------------
  real, intent ( in ) :: p        !~ Pressure ( Pa )
  real, intent ( in ) :: tK       !~ Temperature ( K )
  real, intent ( in ) :: rh       !~ Relative Humidity ( % )
  real                :: wbt      !~ Wet Bulb Temperature ( K )
 
  !~ Utility variables
  !  -----------------
  real                :: tdK      !~ Dewpoint temperature ( K )
  real                :: tC       !~ Temperature ( C )
  real                :: tdC      !~ Dewpoint temperature ( K )
  real                :: svapr    !~ Saturation vapor pressure ( Pa )
  real                :: vapr     !~ Ambient vapor pressure ( Pa )
  real                :: gamma    !~ Dummy term
  real                :: delta    !~ Dummy term

  ! ---------------------------------------------------------------------- !
  ! ---------------------------------------------------------------------- !          
  !~ Initialize variables
  !  --------------------
  wbt = REAL ( 0 )
  tC  = tK - 273.15

  !~ Compute saturation vapor pressure ( Pa )
  !  ----------------------------------------
  svapr = calc_es ( tK ) * REAL ( 100 )

  !~ Compute vapor pressure
  !  ----------------------
  vapr  = svapr * ( rh / REAL (100) )

  !~ Grab the dewpoint
  !  -----------------
  tdC = calc_Dewpoint ( tC, rh )
  tdK = tdC + 273.15

  !~ Compute dummy terms
  !  -------------------
  gamma = 0.00066 * ( p / REAL (1000) )
  delta = REAL ( 4098 ) * ( vapr / REAL(1000) )  / ( (tC+237.3)**2 )

  !~ Compute the wet bulb temperature
  !  --------------------------------
  wbt = ( ((gamma * tC) + (delta * tdC)) / (gamma + delta) ) + 273.15

END FUNCTION
