!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~ 
!~ Name:
!~    calc_rh
!~
!~ Description:
!~    This function calculates relative humidity given pressure, 
!~    temperature, and water vapor mixing ratio.
!~ 
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_rh ( p, t, qv ) result ( rh )
  
  IMPLICIT NONE
 
  REAL, INTENT(IN) :: p, t, qv
  REAL :: rh

  ! Local
  ! -----
  REAL, PARAMETER :: pq0=379.90516
  REAL, PARAMETER :: a2=17.2693882
  REAL, PARAMETER :: a3=273.16
  REAL, PARAMETER :: a4=35.86
  REAL, PARAMETER :: rhmin=1.
  REAL :: q, qs
  INTEGER :: i,j,k

  ! Following algorithms adapted from WRFPOST
  ! May want to substitute with another later
  ! -----------------------------------------
  ! qs is g/kg, 
  ! t is K
  ! p is Pa
    q=qv/(1.0+qv)
    qs=pq0/p*exp(a2*(t-a3)/(t-a4)) 

    rh=100.*q/qs

    IF (rh .gt. 100.) THEN
      rh=100.
    ELSE IF (rh .lt. rhmin) THEN
      rh=rhmin
    ENDIF

END FUNCTION
