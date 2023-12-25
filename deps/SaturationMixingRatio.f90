!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    SaturationMixingRatio
!~
!~ Description:
!~    This function calculates saturation mixing ratio given the
!~    temperature ( K ) and the ambient pressure ( Pa ).  Uses 
!~    approximation of saturation vapor pressure.
!~
!~ References:
!~    Bolton (1980), Monthly Weather Review, pg. 1047, Eq. 10
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION SaturationMixingRatio ( tK, p ) result ( ws )

  IMPLICIT NONE

  REAL, INTENT ( IN ) :: tK
  REAL, INTENT ( IN ) :: p
  REAL                :: ws

  REAL :: es

  es = 6.122 * exp ( (17.67*(tK-273.15))/ (tK-29.66) )
  ws = ( 0.622*es ) / ( (p/100.0)-es )

END FUNCTION
