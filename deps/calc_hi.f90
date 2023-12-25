!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    calc_hi   
!~
!~ Description:
!~    This subroutine calculates the heat index.  Requires temperature ( K )
!~    and relative humidity ( % ).
!~ 
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION calc_hi ( Tk, RH ) result ( HI )

  implicit none

  !~ Variable declarations
  !  ---------------------
  real, intent ( in  ) :: Tk
  real, intent ( in  ) :: RH

  real :: tF, tF2, rh2, HI

  !~ If temperature > 70F then calculate heat index, else set it equal
  !~ to dry temperature
  !  -----------------------------------------------------------------
  IF ( Tk > 294.26111 ) THEN

    tF   = ( (Tk - 273.15) * (REAL (9)/REAL (5))  ) + REAL ( 32 )
    tF2  = tF ** 2
    rh2  = RH ** 2

    HI =  -42.379 &
       +  (  2.04901523   * tF              ) &
       +  ( 10.14333127   * RH              ) &
       -  (  0.22475541   * tF  * RH        ) &
       -  (  6.83783E-03  * tF2             ) &
       -  (  5.481717E-02 * rh2             ) &
       +  (  1.22874E-03  * tF2 * RH        ) &
       +  (  8.5282E-04   * tF  * rh2       ) &
       -  (  1.99E-06     * tF2 * rh2       )

    HI = ((HI - REAL (32)) * (REAL (5)/REAL (9))) + 273.15
  ELSE
    HI = Tk
  END IF

END FUNCTION
