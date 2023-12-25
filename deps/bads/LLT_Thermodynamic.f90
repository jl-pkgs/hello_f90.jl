!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                           ~!
!~ Name:                                                                     ~!
!~    LLT_Thermodynamic                                                      ~!
!~                                                                           ~!
!~ Description:                                                              ~!
!~    This function computes the thermodynamic term for the low-level        ~!
!~    turbulence algorithm.                                                  ~!
!~                                                                           ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION LLT_Thermodynamic ( nlayer, tK, hgt ) RESULT ( thermodynamic )
IMPLICIT NONE
 
  !~ Variable Declaration
  !  --------------------
  INTEGER, INTENT ( IN )         :: nlayer
  REAL, INTENT ( IN )            :: tK     ( nlayer ) !~ Temperature (K)
  REAL, INTENT ( IN )            :: hgt    ( nlayer ) !~ Heights ( m )
  REAL                           :: thermodynamic

  !~ Internal function variables
  !  ---------------------------
  INTEGER :: i
  REAL    :: lapse
  REAL    :: PI
     PARAMETER ( PI = 3.14159265359 )

  !  --------------------------------------------------------------------  !
  !  --------------------------------------------------------------------  !

  !~ Initialize variables
  !  --------------------
  thermodynamic = REAL ( 0 )

  !~ Compute the lapse rate over the layer.  The sign gets goofy here,
  !~ but works as coded below.
  !  -----------------------------------------------------------------
  lapse = ( tk(1) - tk(nlayer) ) * REAL ( 1000 )
  lapse = lapse / ( hgt(nlayer) - hgt(1) )

  !~ Compute the thermodynamic component
  !  -----------------------------------
  thermodynamic = lapse / REAL ( 10 )
  thermodynamic = ( thermodynamic + REAL (3) ) / REAL ( 2 )
  thermodynamic = SIN ( thermodynamic * PI )
  thermodynamic = ( thermodynamic + REAL (1) ) / REAL ( 2 )

END FUNCTION
