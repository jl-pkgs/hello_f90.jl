!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    The2T.f90
!~
!~ Description:
!~    This function returns the temperature at any pressure level along a
!~    saturation adiabat by iteratively solving for it from the parcel
!~    thetae.
!~
!~ Dependencies:
!~    function thetae.f90
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION The2T ( thetaeK, pres, flag ) result ( tparcel )
IMPLICIT NONE

   !~ Variable Declaration
   !  --------------------
   REAL,    INTENT     ( IN ) :: thetaeK
   REAL,    INTENT     ( IN ) :: pres
   LOGICAL, INTENT ( INOUT )  :: flag
   REAL                       :: tparcel

   REAL :: thetaK
   REAL :: tovtheta
   REAL :: tcheck
   REAL :: svpr, svpr2
   REAL :: smixr, smixr2
   REAL :: thetae_check, thetae_check2
   REAL :: tguess_2, correction

   LOGICAL :: found
   INTEGER :: iter

   REAL :: R     ! Dry gas constant
   REAL :: Cp    ! Specific heat for dry air
   REAL :: kappa ! Rd / Cp
   REAL :: Lv    ! Latent heat of vaporization at 0 deg. C

   R     = 287.04
   Cp    = 1004.67
   Kappa = R/Cp
   Lv    = 2.500E+6

   !~ Make initial guess for temperature of the parcel
   !  ------------------------------------------------
   tovtheta = (pres/100000.0)**(r/cp)
   tparcel  = thetaeK/exp(lv*.012/(cp*295.))*tovtheta

   iter = 1
   found = .false.
   flag = .false.

   DO
      IF ( iter > 105 ) EXIT

      tguess_2 = tparcel + REAL ( 1 )

      svpr   = 6.122 * exp ( (17.67*(tparcel-273.15)) / (tparcel-29.66) )
      smixr  = ( 0.622*svpr ) / ( (pres/100.0)-svpr )
      svpr2  = 6.122 * exp ( (17.67*(tguess_2-273.15)) / (tguess_2-29.66) )
      smixr2 = ( 0.622*svpr2 ) / ( (pres/100.0)-svpr2 )

      !  ------------------------------------------------------------------ ~!
      !~ When this function was orinially written, the final parcel         ~!
      !~ temperature check was based off of the parcel temperature and      ~!
      !~ not the theta-e it produced.  As there are multiple temperature-   ~!
      !~ mixing ratio combinations that can produce a single theta-e value, ~!
      !~ we change the check to be based off of the resultant theta-e       ~!
      !~ value.  This seems to be the most accurate way of backing out      ~!
      !~ temperature from theta-e.                                          ~!
      !~                                                                    ~!
      !~ Rentschler, April 2010                                             ~!
      !  ------------------------------------------------------------------  !

      !~ Old way...
      !thetaK = thetaeK / EXP (lv * smixr  /(cp*tparcel) )
      !tcheck = thetaK * tovtheta

      !~ New way
      thetae_check  = Thetae ( tparcel,  pres/100., 100., smixr  )
      thetae_check2 = Thetae ( tguess_2, pres/100., 100., smixr2 )

      !~ Whew doggies - that there is some accuracy...
      !IF ( ABS (tparcel-tcheck) < .05) THEN
      IF ( ABS (thetaeK-thetae_check) < .001) THEN
         found = .true.
         flag  = .true.
         EXIT
      END IF

      !~ Old
      !tparcel = tparcel + (tcheck - tparcel)*.3

      !~ New
      correction = ( thetaeK-thetae_check ) / ( thetae_check2-thetae_check )
      tparcel = tparcel + correction

      iter = iter + 1
   END DO

   !IF ( .not. found ) THEN
   !   print*, "Warning! Thetae to temperature calculation did not converge!"
   !   print*, "Thetae ", thetaeK, "Pressure ", pres
   !END IF

END FUNCTION
