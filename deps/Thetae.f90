!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    Thetae
!~
!~ Description:
!~    This function returns equivalent potential temperature using the 
!~    method described in Bolton 1980, Monthly Weather Review, equation 43.
!~
!~ Input:
!~ - tK    : Temperature ( K )
!~ - p     : Pressure ( hPa )
!~ - rh    : Relative humidity (%)
!~ - mixr  : Mixing Ratio ( g kg^-1)
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION Thetae ( tK, p, rh, mixr )
IMPLICIT NONE

   !~ Variable Declarations
   !  ---------------------
   REAL :: tK        ! Temperature ( K )
   REAL :: p         ! Pressure ( hPa )
   REAL :: rh        ! Relative humidity
   REAL :: mixr      ! Mixing Ratio ( kg kg^-1)
   REAL :: te        ! Equivalent temperature ( K )
   REAL :: thetae    ! Equivalent potential temperature

   REAL, PARAMETER :: R  = 287.04         ! Universal gas constant (J/deg kg)
   REAL, PARAMETER :: P0 = 1000.0         ! Standard pressure at surface (hPa)
   REAL, PARAMETER :: lv = 2.54*(10**6)   ! Latent heat of vaporization
                                          ! (J kg^-1)
   REAL, PARAMETER :: cp = 1004.67        ! Specific heat of dry air constant
                                          ! at pressure (J/deg kg)
   REAL :: tlc                            ! LCL temperature

   !~ Calculate the temperature of the LCL
   !  ------------------------------------
   tlc = TLCL ( tK, rh )

   !~ Calculate theta-e
   !  -----------------
   thetae = (tK * (p0/p)**( (R/Cp)*(1.- ( (.28E-3)*mixr*1000.) ) ) )* &
               exp( (((3.376/tlc)-.00254))*&
                  (mixr*1000.*(1.+(.81E-3)*mixr*1000.)) )
END FUNCTION
