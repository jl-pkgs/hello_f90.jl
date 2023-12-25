!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                           ~!
!~ Name:                                                                     ~!
!~    LLT_TrappedWave                                                        ~!
!~                                                                           ~!
!~ Description:                                                              ~!
!~    This function computes the trapped wave term for the low-level         ~!
!~    turbulence algorithm.                                                  ~!
!~                                                                           ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION LLT_TrappedWave ( nlayer, u, v, p ) RESULT ( trapped )
   IMPLICIT NONE

   !~ Variable Declaration
   !  --------------------
   INTEGER, INTENT ( IN )         :: nlayer
   REAL, INTENT ( IN )            :: u     ( nlayer )
   REAL, INTENT ( IN )            :: v     ( nlayer )
   REAL, INTENT ( IN )            :: p     ( nlayer )
   REAL                           :: trapped

   !~ Internal function variables
   !  ---------------------------
   INTEGER           :: i
   REAL              :: du, dv
   REAL              :: scale_fact, this_p
   REAL              :: dudv, this_dudv
   REAL              :: PI
      PARAMETER ( PI = 3.14159265359 )

   !~ Scale parameters
   !  ----------------
   REAL, PARAMETER :: scale_950 = 0.050000  !~ 1/20
   REAL, PARAMETER :: scale_925 = 0.040000  !~ 1/25
   REAL, PARAMETER :: scale_900 = 0.025000  !~ 1/40
   REAL, PARAMETER :: scale_850 = 0.010000  !~ 1/100
   REAL, PARAMETER :: scale_800 = 0.005000  !~ 1/200
   REAL, PARAMETER :: scale_750 = 0.002941  !~ 1/340
   REAL, PARAMETER :: scale_700 = 0.001923  !~ 1/520
   REAL, PARAMETER :: scale_650 = 0.001351  !~ 1/740
   REAL, PARAMETER :: scale_600 = 0.001000  !~ 1/1000
   REAL, PARAMETER :: scale_550 = 0.000800  !~ 1/1250

   !  --------------------------------------------------------------------  !
   !  --------------------------------------------------------------------  !

   !~ Initialize variables
   !  --------------------
   trapped = REAL ( 0 )

   !~ Compute the trapped wave term
   !  ------------------
   dudv = REAL ( 0 )
   DO i = 2, nlayer

     !~ Compute dudv first
     !  ------------------
     du         = u ( i-1 ) - u ( i )
     dv         = v ( i-1 ) - v ( i )

     !~ Scale based on pressure level
     !  -----------------------------
     this_p = p ( i ) / REAL ( 100 )
     IF ( this_p > REAL (950) ) THEN
       scale_fact = scale_950
     ELSE IF ( this_p <= REAL (950) .AND. this_p > REAL (925) ) THEN
       scale_fact = scale_925
     ELSE IF ( this_p <= REAL (925) .AND. this_p > REAL (900) ) THEN
       scale_fact = scale_900
     ELSE IF ( this_p <= REAL (900) .AND. this_p > REAL (850) ) THEN
       scale_fact = scale_850
     ELSE IF ( this_p <= REAL (850) .AND. this_p > REAL (800) ) THEN
       scale_fact = scale_800
     ELSE IF ( this_p <= REAL (800) .AND. this_p > REAL (750) ) THEN
       scale_fact = scale_750
     ELSE IF ( this_p <= REAL (750) .AND. this_p > REAL (700) ) THEN
       scale_fact = scale_700
     ELSE IF ( this_p <= REAL (700) .AND. this_p > REAL (650) ) THEN
       scale_fact = scale_650
     ELSE IF ( this_p <= REAL (650) .AND. this_p > REAL (600) ) THEN
       scale_fact = scale_600
     ELSE IF ( this_p <= REAL (600) ) THEN
       scale_fact = scale_550
     END IF

     this_dudv = ( (du**2)*(dv**2) ) * scale_fact
     IF ( this_dudv > dudv ) dudv = this_dudv

  END DO

  trapped = dudv
  IF ( trapped > REAL ( 1 ) ) trapped = REAL ( 1 )
  trapped = trapped / REAL ( 4 )

END FUNCTION
