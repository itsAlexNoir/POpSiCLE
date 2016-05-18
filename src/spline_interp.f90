!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: spline_interp
!> \author Alejandro de la Calle, Queen's University Belfast,
!> \author Daniel Dundas, Queen's University Belfast,
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Interpolating routines,
!> \details This module contains routines for a unidimensional
!> spline interpolation.
!
!------------------------------------------------------------------------------
MODULE spline_interp
  
  USE constants
  USE tools
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   :: get_scnd_derivatives_spline
  PUBLIC   :: spline_interpolation
  PUBLIC   :: tridag

  ! Interface
  INTERFACE tridag
     MODULE PROCEDURE tridag_par
  END INTERFACE tridag
  
CONTAINS
  
  !**********************************************************************!
  
  SUBROUTINE get_scnd_derivatives_spline(x,y,yp1,ypn,y2)

    IMPLICIT NONE
    REAL(dp), INTENT(IN)           :: x(:),y(:)
    REAL(dp), INTENT(IN)           :: yp1,ypn
    REAL(dp), INTENT(OUT)          :: y2(:)
    
    INTEGER                        :: n
    REAL(dp), DIMENSION(SIZE(x))   :: a,b,c,r
    !---------------------------------------------------

    IF((SIZE(x).NE.SIZE(y)) .OR. (SIZE(x).NE.SIZE(y2))) THEN
       WRITE(*,*) 'Input arrays do not have the same size in get second derivatives.'
       STOP
    ENDIF
    
    n = SIZE(x)
    ! Set up the tridiagonal equations
    c(1:n-1) = x(2:n) - x(1:n-1)
    r(1:n-1) = 6.0_dp * ( (y(2:n) - y(1:n-1)) / c(1:n-1) )
    r(2:n-1) = r(2:n-1) - r(1:n-2)
    a(2:n-1) = c(1:n-2)
    b(2:n-1) = 2.0_dp * (c(2:n-1) + a(2:n-1))
    b(1)     = 1.0_dp
    b(n)     = 1.0_dp
    
    IF (yp1 > 0.99e30_dp) THEN
       r(1) = 0.0_dp
       c(1) = 0.0_dp
    ELSE
       r(1)=(3.0_dp / (x(2)-x(1))) * ((y(2)-y(1)) / (x(2)-x(1)) - yp1)
       c(1)=0.5_dp
    ENDIF
    
    IF (ypn > 0.99e30_dp) THEN
       r(n)=0.0_dp
       a(n)=0.0_dp
    ELSE
       
       r(n) = ( -3.0_dp / (x(n)-x(n-1))) * &
            ((y(n)-y(n-1)) / (x(n)-x(n-1)) - ypn)
       a(n) = 0.5_dp
    ENDIF

    CALL tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    
  END SUBROUTINE get_scnd_derivatives_spline
  
  !*************************************************************!
  
  SUBROUTINE spline_interpolation(xa,ya,y2a,x,y)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)   :: xa(:), ya(:), y2a(:)
    REAL(dp), INTENT(IN)   :: x
    REAL(dp), INTENT(OUT)  :: y
    
    INTEGER                :: khi, klo, n
    INTEGER                :: left, mflag
    REAL(dp)               :: a,b,h
    !----------------------------------------------------!

    IF((SIZE(xa).NE.SIZE(ya)) .OR. (SIZE(xa).NE.SIZE(y2a))) THEN
       WRITE(*,*) 'Input arrays do not have the same size .'
       STOP
    ENDIF

    n = SIZE(xa)
    CALL interv(xa,SIZE(xa),x,left,mflag)
    klo = MAX(MIN(left,n-1),1)
    
    khi = klo + 1
    h   = xa(khi) - xa(klo)
    
    IF (h == 0.0_dp) THEN
       WRITE(*,*) 'Bad xa input in spline interpolation.'
       STOP
    ENDIF
    
    a = (xa(khi) - x) / h
    b = (x - xa(klo)) / h
    
    y = a * ya(klo) + b * ya(khi) + &
         ( (a**3 - a) * y2a(klo) +(b**3 - b) * y2a(khi) ) * &
         (h**2) / 6.0_dp
    
  END SUBROUTINE spline_interpolation

  !*************************************************!
  
  SUBROUTINE tridag_ser(a,b,c,r,u)
    
    IMPLICIT NONE
    REAL(dp), INTENT(IN)         :: a(:), b(:), c(:), r(:)
    REAL(dp), INTENT(OUT)        :: u(:)
    
    REAL(dp), DIMENSION(SIZE(b)) :: gam
    INTEGER                      :: n,j
    REAL(dp)                     :: bet
    !-------------------------------------------------!

    IF((SIZE(a)+1.NE.SIZE(b)) .OR. (SIZE(a)+1.NE.SIZE(c)+1) &
         .OR. (SIZE(a)+1.NE.SIZE(r)) .OR. (SIZE(a)+1.NE.SIZE(u))) THEN
       WRITE(*,*) 'Input arrays do not have the same size in tridag.'
       STOP
    ENDIF
    n = SIZE(b)
    bet=b(1)
    
    IF (bet == 0.0_dp) THEN
       WRITE(*,*) 'Error in tridag_ser. (I)'
       STOP
    ENDIF
    
    u(1) = r(1)/bet
    DO j = 2, n
       gam(j) = c(j-1) / bet
       bet = b(j) - a(j-1) * gam(j)
       IF (bet == 0.0_dp) THEN
          WRITE(*,*) 'Error in tridag. (II)'
          STOP
       ENDIF
       u(j) = (r(j) - a(j-1) * u(j-1)) / bet
    ENDDO
    
    DO  j = n-1, 1, -1
       u(j) = u(j) - gam(j+1) * u(j+1)
    ENDDO
    
  END SUBROUTINE tridag_ser

  !*************************************************!

  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    
    IMPLICIT NONE
    REAL(dp), INTENT(IN)             :: a(:), b(:), c(:), r(:)
    REAL(dp), INTENT(OUT)            :: u(:)
    
    INTEGER, PARAMETER               :: NPAR_TRIDAG=4
    REAL(dp), DIMENSION(size(b)/2)   :: y,q,piva
    REAL(dp), DIMENSION(size(b)/2-1) :: x,z
    REAL(dp), DIMENSION(size(a)/2)   :: pivc
    INTEGER                          :: n,n2,nm,nx 

    !---------------------------------------------------!

    IF((SIZE(a)+1.NE.SIZE(b)) .OR. (SIZE(a)+1.NE.SIZE(c)+1) &
         .OR. (SIZE(a)+1.NE.SIZE(r)) .OR. (SIZE(a)+1.NE.SIZE(u))) THEN
       WRITE(*,*) 'Input arrays do not have the same size in tridag.'
       STOP
    ENDIF

    n = SIZE(b)
    IF (n < NPAR_TRIDAG) THEN
       CALL tridag_ser(a,b,c,r,u)
    ELSE
       IF (MAXVAL(ABS(b(1:n))) == 0.0_dp) THEN
          WRITE(*,*) 'Possible singular matrix in tridag.'
          STOP
       ENDIF
       
       n2 = SIZE(y)
       nm = SIZE(pivc)
       nx = SIZE(x)
       piva = a(1:n-1:2)/b(1:n-1:2)
       pivc = c(2:n-1:2)/b(3:n:2)
       y(1:nm) = b(2:n-1:2) - piva(1:nm) * c(1:n-2:2) - &
            pivc * a(2:n-1:2)
       q(1:nm) = r(2:n-1:2) - piva(1:nm) * r(1:n-2:2) - &
            pivc * r(3:n:2)
       
       IF (nm < n2) THEN
          y(n2) = b(n) - piva(n2) * c(n-1)
          q(n2) = r(n) - piva(n2) * r(n-1)
       ENDIF
       x = -piva(2:n2) * a(2:n-2:2)
       z = -pivc(1:nx) * c(3:n-1:2)
       
       CALL tridag_par(x,y,z,q,u(2:n:2))
       
       u(1) = (r(1)-c(1)*u(2))/b(1)
       u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
            -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
       IF (nm == n2) &
            u(n) = (r(n) - a(n-1) * u(n-1)) / b(n)
    ENDIF
    
  END SUBROUTINE tridag_par
  
  !*************************************************!
  !*************************************************!  
  
END MODULE spline_interp
