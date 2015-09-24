MODULE gaussleg
  
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC      :: get_gauss_stuff
  
CONTAINS

  SUBROUTINE get_gauss_stuff(lower_bound, upper_bound, axis,weights)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)         :: lower_bound
    REAL(dp), INTENT(IN)         :: upper_bound
    REAL(dp), INTENT(OUT)        :: axis(:)
    REAL(dp), INTENT(OUT)        :: weigths(:)
    
    INTEGER                      :: i,j
    INTEGER                      :: n, m
    REAL(dp), PARAMETER          :: precision = 1.0e-14_dp
    REAL(dp)                     :: z1, z, xm, xl,
    REAL(dp)                     :: pp, p3, p2, p1
    
    ! Check the size of the arrays. They must be equal
    IF (SIZE(axis).NE.SIZE(weights)) THEN
       WRITE(*,*) 'Axis and weights arrays have not equal length.'
       STOP
    ENDIF
    
    n = SIZE(axis)
    m = (n+1) / 2
    
    xm = 0.5_dp * (upper_bound + lower_bound)
    xl = 0.5_dp * (upper_bound - lower_bound)
    
    DO i = 1, m
       z = COS(pi * (REAL(i,dp) + 0.75_dp) / (REAL(n,dp) + 0.5_dp))
       DO WHILE (abs(z-z1) > precision)
          p1 = 1.0_dp
          p2 = 0.0_dp
          DO j = 1, n
             p3 = p2
             p2 = p1
             p1 = (REAL(2*j+1,dp) * z * p2 - REAL(j,dp) * p3) / REAL(j+1,dp)
          ENDDO
          pp = REAL(n,dp) * (z * p1 - p2) / (z * z - 1.0_dp)
          z1 = z
          z = z1 - p1 / pp
       ENDDO
       
    ENDDO
    
    axis(1:m) = xm - xl * z
    axisx(n:n-m+1:-1) = xm + xl * z
    weights(1:m) = 2.0_dp * xl / ((1.0 - z * z) * pp * pp)
    weights(n:n-m+1:-1) = w(1:m)
    
    
  END SUBROUTINE get_gauss_stuff

  !------------------------------------------------------------------------------!



  
  
END MODULE gaussleg
