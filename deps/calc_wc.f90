!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    calc_wc   
!~
!~ Description:
!~    This function calculates wind chill given temperature ( K ) and
!~    wind speed ( m/s )
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_wc ( tK, wspd ) RESULT ( wc )

  implicit none

  !~ Variable Declarations
  !  ---------------------
  real, intent ( in  ) :: tK
  real, intent ( in  ) :: wspd

  real                 :: tF, wc, wspd_mph

  wspd_mph = wspd * 2.23693629 ! convert to mph
  tF  = (( tK - 273.15 ) * ( REAL (9) / REAL (5) ) ) + REAL ( 32 )

  wc =    35.74                           &
     +  (  0.6215 * tF                  ) &
     -  ( 35.75   *      ( wspd_mph**0.16 ) ) &
     +  (  0.4275 * tF * ( wspd_mph**0.16 ) )

  wc = (( wc - REAL (32) ) * ( REAL (5) / REAL (9) ) ) + 273.15

END FUNCTION
