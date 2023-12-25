!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
!~                                                                     
!~ Name:                                                                
!~    tlcl                                                               
!~                                                                        
!~ Description:                                                            
!~    This function calculates the temperature of a parcel of air would have
!~    if lifed dry adiabatically to it's lifting condensation level (lcl).  
!~                                                                          
!~ References:                                                              
!~    Bolton (1980), Monthly Weather Review, pg. 1048, Eq. 22
!~                                                                          
!!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!!
FUNCTION TLCL ( tk, rh )
  
  IMPLICIT NONE
 
  REAL, INTENT ( IN ) :: tK   !~ Temperature ( K )
  REAL, INTENT ( IN ) :: rh   !~ Relative Humidity ( % )
  REAL                :: tlcl
  
  REAL :: denom, term1, term2

  term1 = 1.0 / ( tK - 55.0 )
  IF ( rh > REAL (0) ) THEN
    term2 = ( LOG (rh/100.0)  / 2840.0 )
  ELSE
    term2 = ( LOG (0.001/1.0) / 2840.0 )
  END IF
  denom = term1 - term2
  tlcl = ( 1.0 / denom ) + REAL ( 55 ) 

END FUNCTION
