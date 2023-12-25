!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    PWat
!~
!~ Description:
!~    This function calculates precipitable water by summing up the 
!~    water vapor in a column from the first eta layer to model top
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION Pwat  ( nz, qv, qc, dz8w, rho )

  IMPLICIT NONE

   !~ Variable declaration
   !  --------------------
   INTEGER, INTENT ( IN ) :: nz          !~ Number of vertical levels
   REAL, INTENT ( IN )    :: qv   ( nz ) !~ Specific humidity in layer (kg/kg)
   REAL, INTENT ( IN )    :: qc   ( nz ) !~ Cloud water in layer (kg/kg)
   REAL, INTENT ( IN )    :: dz8w ( nz ) !~ Dist between vertical levels (m)
   REAL, INTENT ( IN )    :: rho  ( nz ) !~ Air density (kg/m^3)
   REAL                   :: Pwat        !~ Precipitable water (kg/m^2)
   INTEGER                :: k           !~ Vertical index

   !~ Precipitable water (kg/m^2)
   !  ---------------------------
   Pwat=0
   DO k = 1, nz
     !Based on AMS PWAT defination (https://glossary.ametsoc.org/wiki/Precipitable_water)
     !PWAT is corrected as the column accumulated water vapor rather than water vapor + cloud water.
     !Modified by Zhixiao
     Pwat = Pwat + qv(k) * dz8w(k) * rho(k)
   ENDDO
           
END FUNCTION
