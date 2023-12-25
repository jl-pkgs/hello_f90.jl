!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    calc_es 
!~
!~ Description:
!~    This function returns the saturation vapor pressure over water ( hPa )
!~    given temperature ( K ).
!~
!~ References:
!~    The algorithm is due to Nordquist, W.S., 1973: "Numerical approximations
!~    of selected meteorological parameters for cloud physics problems,"
!~    ecom-5475, Atmospheric Sciences Laboratory, U.S. Army Electronics
!~    Command, White Sands Missile Range, New Mexico, 88002
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_es ( tK ) result ( es )

  implicit none

  !~ Variable Declaration
  !  --------------------
  real, intent ( in ) :: tK
  real                :: es
 
  real                :: p1, p2, c1

  p1 = 11.344    - 0.0303998 * tK
  p2 = 3.49149   - 1302.8844 / tK
  c1 = 23.832241 - 5.02808   * ALOG10 ( tK )
  es = 10.**(c1-1.3816E-7*10.**p1+8.1328E-3*10.**p2-2949.076/tK)

END FUNCTION
