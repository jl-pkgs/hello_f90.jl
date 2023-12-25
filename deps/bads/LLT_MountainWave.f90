!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                           ~!
!~ Name:                                                                     ~!
!~    LLT_MountainWave                                                       ~!
!~                                                                           ~!
!~ Description:                                                              ~!
!~    This function computes the mountain wave term for the low-level        ~!
!~    turbulence algorithm.                                                  ~!
!~                                                                           ~!
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION LLT_MountainWave ( nlayer, tdx, tdy, u, v, tK, hgt) &
                                                       RESULT ( MountainWave )
  IMPLICIT NONE

  !~ Variable Declaration
  !  --------------------
  INTEGER, INTENT ( IN )           :: nlayer
  REAL, INTENT ( IN )              :: tdx               !~ Terrain dx
  REAL, INTENT ( IN )              :: tdy               !~ Terrain dy
  REAL, INTENT ( IN )              :: u   ( nlayer )    !~ U components f
  REAL, INTENT ( IN )              :: v   ( nlayer )    !~ V components
  REAL, INTENT ( IN )              :: tK  ( nlayer )    !~ Temperatures (K)
  REAL, INTENT ( IN )              :: hgt ( nlayer )    !~ Heights ( m )
  REAL                             :: MountainWave      !~ Mountain wave term
 
  !~ Internal function variables
  !  ---------------------------
  REAL    :: u_term
  REAL    :: v_term
  REAL    :: uv_term
  REAL    :: lapse
  REAL    :: total_mw, this_total_mw
  REAL    :: this_uv_term
  REAL    :: min_uv_term, cross_terrain, max_total_mw
  INTEGER :: i, j, k

  REAL    :: PI
     PARAMETER ( PI = 3.14159265359 )

  !  --------------------------------------------------------------------  !
  !  --------------------------------------------------------------------  !

  !~ Initialize variables
  !  --------------------
  MountainWave = REAL ( 0 )

  !~ Loop through the layer
  !  ----------------------
  DO i = 2, nlayer

    !~ Wind terrain term
    !  -----------------
    u_term       = ( (u(i-1) + u(i) ) / REAL(2) ) * tdx
    v_term       = ( (v(i-1) + v(i) ) / REAL(2) ) * tdy
    this_uv_term = ( u_term + v_term ) * REAL ( -1 )
    !IF ( uv_term < REAL (0) ) uv_term = REAL ( 0 )
    IF ( min_uv_term < REAL (0) ) min_uv_term = REAL ( 0 )
    IF ( i == 2 ) THEN
      !uv_term = this_uv_term
      min_uv_term = this_uv_term
    ELSE
      !IF ( this_uv_term < uv_term ) uv_term = this_uv_term
      IF ( this_uv_term < min_uv_term ) min_uv_term = this_uv_term
    END IF

    !~ Lapse rate
    !  ----------
    lapse = ( tK (i-1) - tK (i) ) * REAL ( 1000 )
    lapse  = lapse / ABS ( hgt(i)-hgt(i-1) )
    IF ( lapse > REAL (0) ) lapse = REAL ( 0 )
    lapse = lapse * REAL ( -1 )

    this_total_mw = this_uv_term * lapse * REAL ( 40000 )
    IF ( i == 2 ) THEN
      total_mw = this_total_mw
    ELSE
      IF ( this_total_mw > total_mw ) total_mw = this_total_mw
    END IF
 
  END DO

  !min_uv_term = uv_term
  cross_terrain = min_uv_term * REAL ( 500 )

  IF ( min_uv_term < 0.03 ) THEN
    cross_terrain = REAL ( 0 )
  END IF

  IF ( cross_terrain > REAL (50) ) cross_terrain = REAL ( 50 )

  !~ Multiply the lapse (inversion) array and the mountain wave array
  !  ----------------------------------------------------------------
  IF ( total_mw > REAL (50) ) total_mw = REAL ( 50 )

  !~ Add the cross terrain flow and inversion term
  !  ---------------------------------------------
  MountainWave = ( total_mw*(cross_terrain/50.) ) + cross_terrain
  MountainWave = MountainWave / REAL ( 100 )

END FUNCTION
