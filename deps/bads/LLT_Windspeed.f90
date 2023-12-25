!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                           ~!
!~ Name:                                                                     ~!
!~    LLT_Windspeed                                                          ~!
!~                                                                           ~!
!~ Description:                                                              ~!
!~    This function computes the dynamic term for the low-level turbulence   ~!
!~    algorithm.                                                             ~!
!~                                                                           ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION LLT_Windspeed ( nlayer, u, v ) RESULT ( dynamic )
  IMPLICIT NONE
 
  !~ Variable Declaration
  !  --------------------
  INTEGER, INTENT ( IN )         :: nlayer
  REAL, INTENT ( IN )            :: u     ( nlayer )
  REAL, INTENT ( IN )            :: v     ( nlayer )
  REAL                           :: dynamic
 
  !~ Internal function variables
  !  ---------------------------
  INTEGER           :: i
  REAL              :: this_windspeed ( nlayer )
  REAL              :: PI
     PARAMETER ( PI = 3.14159265359 )

  !  --------------------------------------------------------------------  !
  !  --------------------------------------------------------------------  !           
  !~ Initialize variables
  !  --------------------
  dynamic = REAL ( 0 )

  !~ Compute the windspeed
  !  ---------------------
  DO i = 1, nlayer
     this_windspeed ( i ) = SQRT ( u(i)**2 + v(i)**2 )
  END DO

  !~ Compute the dynamic term
  !  -------------------------
  dynamic = ( this_windspeed(1)+this_windspeed(nlayer) ) / REAL (20)
  IF ( dynamic > REAL (2) ) dynamic = REAL ( 2 )
  dynamic = ( dynamic + REAL (3) ) / REAL ( 2 )
  dynamic = SIN ( dynamic*PI )
  dynamic = ( dynamic + REAL (1) ) / REAL ( 2 )


END FUNCTION
