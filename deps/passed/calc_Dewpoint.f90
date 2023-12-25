!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    calc_Dewpoint
!~
!~ Description:
!~    This function approximates dewpoint given temperature and rh.
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_Dewpoint ( tC, rh) result( Dewpoint )

  implicit none

  !~ Variable Declaration
  !  --------------------
  real, intent ( in ) :: tC
  real, intent ( in ) :: rh
  real                :: Dewpoint
 
  real :: term, es, e1, e, logs, expon

  expon    = ( 7.5*tC ) / ( 237.7+tC )
  es       = 6.112 * ( 10**expon )     ! Saturated vapor pressure
  e        = es * ( rh/100.0 )         ! Vapor pressure
  logs     = LOG10 ( e/6.112 )
  Dewpoint = ( 237.7*logs ) / ( 7.5-logs )

END FUNCTION
