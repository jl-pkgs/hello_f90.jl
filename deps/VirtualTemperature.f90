!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~
!~ Name:
!~    VirtualTemperature
!~
!~ Description:
!~    This function returns virtual temperature given temperature ( K )
!~    and mixing ratio.
!~
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION VirtualTemperature ( tK, w ) result ( Tv )
IMPLICIT NONE

   !~ Variable declaration
   real, intent ( in ) :: tK !~ Temperature
   real, intent ( in ) :: w  !~ Mixing ratio ( kg kg^-1 )
   real                :: Tv !~ Virtual temperature

   Tv = tK * ( 1.0 + (w/0.622) ) / ( 1.0 + w )

END FUNCTION
