!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    Theta
!~
!~ Description:
!~    This function calculates potential temperature as defined by
!~    Poisson's equation, given temperature and pressure ( hPa ).
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION Theta ( t, p )
IMPLICIT NONE

   !~ Variable declaration
   !  --------------------
   REAL, INTENT ( IN ) :: t
   REAL, INTENT ( IN ) :: p
   REAL                :: theta

   REAL :: Rd ! Dry gas constant
   REAL :: Cp ! Specific heat of dry air at constant pressure
   REAL :: p0 ! Standard pressure ( 1000 hPa )

   Rd =  287.04
   Cp = 1004.67
   p0 = 1000.00

   !~ Poisson's equation
   !  ------------------
   theta = t * ( (p0/p)**(Rd/Cp) )

END FUNCTION
