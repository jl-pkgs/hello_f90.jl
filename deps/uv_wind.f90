!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~ 
!~ Name:
!~    uv_wind
!~
!~ Description:
!~    This function calculates the wind speed given U and V wind
!~    components.
!~ 
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION uv_wind ( u, v ) result ( wind_speed )
 
  IMPLICIT NONE
 
  REAL, INTENT(IN) :: u, v
  REAL :: wind_speed

  wind_speed = sqrt( u*u + v*v )

END FUNCTION
