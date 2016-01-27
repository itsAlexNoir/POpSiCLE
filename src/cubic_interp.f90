!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: cubic_interp
!> \author Alejandro de la Calle, Queen's University Belfast,
!> \author Daniel Dundas, Queen's University Belfast,
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Interpolating routines,
!> \details This module contains all the subprograms needed for
!> an interpolation of the wavefunction, including routines for
!> multidimensions grid, Many of the routines in this module have
!> been extract form the numerical recipes books,
!
!------------------------------------------------------------------------------
MODULE cubic_interp
  
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   :: initialize_bicubic_matrix
  PUBLIC   :: initialize_tricubic_matrix
  PUBLIC   :: bicubic_interpolate
  !PUBLIC   :: tricubic_interpolate
  PUBLIC   :: delete_bicubic_matrix
  PUBLIC   :: delete_tricubic_matrix
  
  REAL(dp), ALLOCATABLE, PUBLIC  :: wbi(:, :)
  REAL(dp), ALLOCATABLE, PUBLIC  :: wtri(:, :)
 
  
CONTAINS

  !**********************************************************************!
  
  SUBROUTINE initialize_bicubic_matrix( )
    IMPLICIT NONE
    
    ALLOCATE(wbi(16,16))

    OPEN(UNIT=2,FILE='bicoeffs.dat',FORM='formatted')
    READ(2,*) wbi
    CLOSE(2)
    
!!$    DATA wbi /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
!!$         8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
!!$         2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
!!$         2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
!!$         -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
!!$         -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
!!$         -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
    
  END SUBROUTINE initialize_bicubic_matrix
  
  !*************************************************!
  
  SUBROUTINE initialize_tricubic_matrix( )
    IMPLICIT NONE
    
    ALLOCATE(wtri(64,64))
    OPEN(UNIT=2,FILE='tricoeffs.dat',FORM='formatted')
    READ(2,*) wtri
    CLOSE(2)
    
!!$    DATA wtri /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,16*0,-3,0,9,-6,4*0,9,&
!!$         0,-27,18,-6,0,18,-12,2,0,-6,4,4*0,-6,0,18,-12,4,0,-12,8,&
!!$         32*0,3,0,-9,6,4*0,-9,0,27,-18,6,0,-18,12,-2,0,6,-4,4*0,6,0,&
!!$         -18,12,-4,0,12,-8,&
!!$         40*0,9,0,-27,18,-6,0,18,-12,8*0,-6,0,18,-12,4,0,-12,8,&
!!$         8*0,3,0,-9,6,-2,0,6,-4,24*0,-9,0,27,-18,6,0,-18,12,8*0,6,0,&
!!$         -18,12,-4,0,12,-8,&
!!$         2*0,3,-2,6*0,-9,6,2*0,6,-4,18*0,-9,6,6*0,27,-18,2*0,-18,12,&
!!$         2*0,6,-4,6*0,-18,12,2*0,12,-8,&
!!$         34*0,9,-6,6*0,-27,18,2*0,18,-12,2*0,-6,4,6*0,18,-12,2*0,-12,8,&
!!$         42*0,27,-18,2*0,-18,12,10*0,-18,12,2*0,12,-8,&
!!$         10*0,9,-6,2*0,-6,4,26*0,-27,18,2*0,18,-12,10*0,18,-12,2*0,-12,8,&
!!$         16*0,1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,-2,0,6,-4,4*0,6,0,-18,12,-4,&
!!$         0,12,-8,1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
!!$         32*0,-1,0,3,-2,4*0,3,0,-9,6,-2,0,6,-4,1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
!!$         40*0,-3,0,9,-6,2,0,-6,4,8*0,3,0,-9,6,-2,0,6,-4,&
!!$         24*0,3,0,-9,6,-2,0,6,-4,8*0,-6,0,18,-12,4,0,-12,8,8*0,3,0,-9,6,-2,0,6,-4,&
!!$         18*0,3,-2,6*0,-9,6,2*0,6,-4,2*0,-6,4,6*0,18,-12,2*0,-12,8,2*0,3,&
!!$         -2,6*0,-9,6,2*0,6,-4,&
!!$         34*0,-3,2,6*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,2*0,6,-4,&
!!$         42*0,-9,6,2*0,6,-4,10*0,9,-6,2*0,-6,4,&
!!$         26*0,9,-6,2*0,-6,4,10*0,-18,12,2*0,12,-8,10*0,9,-6,2*0,-6,4,&
!!$         4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,20*0,-3,0,9,-6,6,0,-18,12,-3,0,9,&
!!$         -6,4*0,2,0,-6,4,-4,0,12,-8,2,0,-6,4,&
!!$         36*0,3,0,-9,6,-6,0,18,-12,3,0,-9,6,4*0,-2,0,6,-4,4,0,-12,8,-2,0,6,-4,&
!!$         40*0,-3,0,9,-6,3,0,-9,6,8*0,2,0,-6,4,-2,0,6,-4,&
!!$         8*0,-1,0,3,-2,1,0,-3,2,24*0,3,0,-9,6,-3,0,9,-6,8*0,-2,0,6,-4,2,0,-6,4,&
!!$         6*0,3,-2,2*0,-6,4,2*0,3,-2,22*0,-9,6,2*0,18,-12,2*0,-9,6,6*0,6,-4,&
!!$         2*0,-12,8,2*0,6,-4,&
!!$         38*0,9,-6,2*0,-18,12,2*0,9,-6,6*0,-6,4,2*0,12,-8,2*0,-6,4,&
!!$         42*0,-9,6,2*0,9,-6,10*0,6,-4,2*0,-6,4,&
!!$         10*0,-3,2,2*0,3,-2,26*0,9,-6,2*0,-9,6,10*0,-6,4,2*0,6,-4,&
!!$         0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,17*0,-3,6,-3,5*0,9,-18,9,0,-6,12,-6,&
!!$         0,2,-4,2,5*0,-6,12,-6,0,4,-8,4,&
!!$         33*0,3,-6,3,5*0,-9,18,-9,0,6,-12,6,0,-2,4,-2,5*0,6,-12,6,0,-4,8,-4,&
!!$         41*0,9,-18,9,0,-6,12,-6,9*0,-6,12,-6,0,4,-8,4,&
!!$         9*0,3,-6,3,0,-2,4,-2,25*0,-9,18,-9,0,6,-12,6,9*0,6,-12,6,0,-4,8,-4,&
!!$         2*0,-1,1,6*0,3,-3,2*0,-2,2,18*0,3,-3,6*0,-9,9,2*0,6,-6,2*0,-2,2,6*0,&
!!$         6,-6,2*0,-4,4,&
!!$         34*0,-3,3,6*0,9,-9,2*0,-6,6,2*0,2,-2,6*0,-6,6,2*0,4,-4,&
!!$         42*0,-9,9,2*0,6,-6,10*0,6,-6,2*0,-4,4,&
!!$         10*0,-3,3,2*0,2,-2,26*0,9,-9,2*0,-6,6,10*0,-6,6,2*0,4,-4,&
!!$         20*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,4*0,-2,0,6,-4,4,0,-12,8,-2,0,6,&
!!$         -4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,&
!!$         36*0,-1,0,3,-2,2,0,-6,4,-1,0,3,-2,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,&
!!$         40*0,1,0,-3,2,-1,0,3,-2,8*0,-1,0,3,-2,1,0,-3,2,&
!!$         24*0,-1,0,3,-2,1,0,-3,2,8*0,2,0,-6,4,-2,0,6,-4,8*0,-1,0,3,-2,1,0,-3,2,&
!!$         22*0,3,-2,2*0,-6,4,2*0,3,-2,6*0,-6,4,2*0,12,-8,2*0,-6,4,6*0,3,-2,2*0,&
!!$         -6,4,2*0,3,-2,&
!!$         38*0,-3,2,2*0,6,-4,2*0,-3,2,6*0,3,-2,2*0,-6,4,2*0,3,-2,&
!!$         42*0,3,-2,2*0,-3,2,10*0,-3,2,2*0,3,-2,&
!!$         26*0,-3,2,2*0,3,-2,10*0,6,-4,2*0,-6,4,10*0,-3,2,2*0,3,-2,&
!!$         17*0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,0,-2,4,-2,5*0,6,-12,6,0,-4,8,-4,0,1,-2,&
!!$         1,5*0,-3,6,-3,0,2,-4,2,&
!!$         33*0,-1,2,-1,5*0,3,-6,3,0,-2,4,-2,0,1,-2,1,5*0,-3,6,-3,0,2,-4,2,&
!!$         41*0,-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,&
!!$         25*0,3,-6,3,0,-2,4,-2,9*0,-6,12,-6,0,4,-8,4,9*0,3,-6,3,0,-2,4,-2,&
!!$         18*0,-1,1,6*0,3,-3,2*0,-2,2,2*0,2,-2,6*0,-6,6,2*0,4,-4,2*0,-1,1,6*0,3,&
!!$         -3,2*0,-2,2,&
!!$         34*0,1,-1,6*0,-3,3,2*0,2,-2,2*0,-1,1,6*0,3,-3,2*0,-2,2,&
!!$         42*0,3,-3,2*0,-2,2,10*0,-3,3,2*0,2,-2,&
!!$         26*0,-3,3,2*0,2,-2,10*0,6,-6,2*0,-4,4,10*0,-3,3,2*0,2,-2,&
!!$         5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,21*0,-3,6,-3,0,6,-12,6,0,-3,6,-3,5*0,2,&
!!$         -4,2,0,-4,8,-4,0,2,-4,2,&
!!$         37*0,3,-6,3,0,-6,12,-6,0,3,-6,3,5*0,-2,4,-2,0,4,-8,4,0,-2,4,-2,&
!!$         41*0,-3,6,-3,0,3,-6,3,9*0,2,-4,2,0,-2,4,-2,&
!!$         9*0,-1,2,-1,0,1,-2,1,25*0,3,-6,3,0,-3,6,-3,9*0,-2,4,-2,0,2,-4,2,&
!!$         6*0,-1,1,2*0,2,-2,2*0,-1,1,22*0,3,-3,2*0,-6,6,2*0,3,-3,6*0,-2,2,2*0,4,&
!!$         -4,2*0,-2,2,&
!!$         38*0,-3,3,2*0,6,-6,2*0,-3,3,6*0,2,-2,2*0,-4,4,2*0,2,-2,&
!!$         42*0,3,-3,2*0,-3,3,10*0,-2,2,2*0,2,-2,&
!!$         10*0,1,-1,2*0,-1,1,26*0,-3,3,2*0,3,-3,10*0,2,-2,2*0,-2,2,&
!!$         21*0,1,-2,1,0,-2,4,-2,0,1,-2,1,5*0,-2,4,-2,0,4,-8,4,0,-2,4,-2,5*0,1,-2,&
!!$         1,0,-2,4,-2,0,1,-2,1,&
!!$         37*0,-1,2,-1,0,2,-4,2,0,-1,2,-1,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,&
!!$         41*0,1,-2,1,0,-1,2,-1,9*0,-1,2,-1,0,1,-2,1,&
!!$         25*0,-1,2,-1,0,1,-2,1,9*0,2,-4,2,0,-2,4,-2,9*0,-1,2,-1,0,1,-2,1,&
!!$         22*0,-1,1,2*0,2,-2,2*0,-1,1,6*0,2,-2,2*0,-4,4,2*0,2,-2,6*0,-1,1,2*0,2,&
!!$         -2,2*0,-1,1,&
!!$         38*0,1,-1,2*0,-2,2,2*0,1,-1,6*0,-1,1,2*0,2,-2,2*0,-1,1,&
!!$         42*0,-1,1,2*0,1,-1,10*0,1,-1,2*0,-1,1,&
!!$         26*0,1,-1,2*0,-1,1,10*0,-2,2,2*0,2,-2,10*0,1,-1,2*0,-1,1/
    
  END SUBROUTINE initialize_tricubic_matrix
  
  !*************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE get_bicubic_coeffs
  !
  !> \brief Calculate coefficients used in bicubic interpolation,
  !> \details Given arrays y, y1, y2, and y12, each of length 4,
  !> containing the function, gradients, and cross derivative at the
  !> four grid points of a rectangular grid cell (numbered counterclockwise
  !> from the lower left), and given d1 and d2, the length of the grid cell
  !> in the 1- and 2- directions, this routine returns the 4 × 4 table c
  !> that is used by routine bcuint for bicubic interpolation,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE get_bicubic_coeffs(y,y1,y2,y12,d1,d2,wt,c)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)  :: d1,d2
    REAL(dp), INTENT(IN)  :: y(4)
    REAL(dp), INTENT(IN)  :: y1(4),y2(4)
    REAL(dp), INTENT(IN)  :: y12(4)
    REAL(dp), INTENT(IN)  :: wt(16,16)
    REAL(dp), INTENT(OUT) :: c(4,4)
    
    REAL(dp)  :: x(16)
    
    x(1:4)=y
    x(5:8)=y1*d1
    x(9:12)=y2*d2
    x(13:16)=y12*d1*d2
    x=matmul(wt,x)
    c=reshape(x,(/4,4/),order=(/2,1/))
    
  END SUBROUTINE get_bicubic_coeffs
  
  !********************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE get_tricubic_coeffs
  !
  !> \brief Calculate coefficients used in bicubic interpolation,
  !> \details Given arrays y, y1, y2, y3 and y12, y13, y23 and y123
  !> each of length 8, containing the function, gradients, and cross derivative
  !> at the eight grid points of a cubic grid cell (numbered counterclockwise
  !> from the lower left), and given d1, d2 and d3, the length of the grid cell
  !> in the 1-, 2- and 3-directions, this routine returns the 8 × 8 table c
  !> that is used by routine bcuint for bicubic interpolation,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y3  Value of the derivative of the function in the third
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] y13 Value of the cross-derivative at known points
  !> \param[in] y23 Value of the cross-derivative at known points
  !> \param[in] y123 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE get_tricubic_coeffs(y,y1,y2,y3,y12,y13,y23,y123,d1,d2,d3,wt,c)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)  :: d1,d2,d3
    REAL(dp), INTENT(IN)  :: y(8)
    REAL(dp), INTENT(IN)  :: y1(8),y2(8),y3
    REAL(dp), INTENT(IN)  :: y12(8),y13(8),y23(8)
    REAL(dp), INTENT(IN)  :: y123(8)
    REAL(dp), INTENT(IN)  :: wt(64,64)
    REAL(dp), INTENT(OUT) :: c(8,8)
    
    REAL(dp)  :: x(64)
 
    !----------------------------------------------!
    
    x(1:8) = y
    x(9:16) = y1 * d1
    x(17:24) = y2 * d2
    x(25:32) = y3 * d3
    x(33:40) = y12 * d1 * d2
    x(41:48) = y13 * d1 * d3
    x(49:56) = y23 * d2 * d3
    x(57:64) = y123 * d1 * d2 * d3
    
    x=matmul(wt,x)
    c=reshape(x,(/8,8/),order=(/2,1/))
    
  END SUBROUTINE get_tricubic_coeffs
  
  !********************************************************!
  
  !---------------------------------------------------------------------------
  !
  !  SUBROUTINE bicubic_interpolate
  !
  !> \brief Get interpolated values from bicubic interpolation
  !> \details  Bicubic interpolation within a grid square,
  !> Input quantities are y,y1,y2,y12, x1l and x1u, the lower and upper
  !> coordinates of the grid square in the 1- direction; x2l and x2u
  !> likewise for the 2-direction; and x1,x2, the coordinates of the
  !> desired point for the interpolation, The interpolated function value
  !> is returned as ansy, and the interpolated gradient values as
  !> ansy1 and ansy2,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE bicubic_interpolate(y,y1,y2,y12,x1l,x1u,x2l,x2u,&
       x1,x2,ansy,ansy1,ansy2)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: y(4)
    REAL(dp), INTENT(IN) :: y1(4),y2(4),y12(4)
    REAL(dp), INTENT(IN) :: x1l,x1u,x2l,x2u
    REAL(dp), INTENT(IN) :: x1,x2
    REAL(dp), INTENT(OUT):: ansy,ansy1,ansy2
    
    REAL(dp) :: t,u
    REAL(dp) :: d1, d2
    REAL(dp) :: c(4,4)
    INTEGER  :: i
    
    !---------------------------------------------------------------!
    
    d1 = x1u - x1l
    d2 = x2u - x2l
    
    IF ( (d1.EQ.0.0_dp).OR.(d2.EQ.0.0_dp)) THEN
       WRITE(*,*) 'Problem with input values in bicubic interpolate,'
       WRITE(*,*) 'Are the boundary pair equal?'
       STOP
    ENDIF
    
    CALL get_bicubic_coeffs(y,y1,y2,y12,d1,d2,wbi,c)
    
    t = (x1-x1l) / d1
    u = (x2-x2l) / d2
    
    ansy = 0.0_dp
    ansy2 = 0.0_dp
    ansy1 = 0.0_dp
    
    DO i = 4, 1, -1
       ansy = t * ansy + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
       ansy2 = t * ansy2 + (3.0_sp * c(i,4) * u + 2.0_dp * c(i,3)) * u + c(i,2)
       ansy1 = u * ansy1 + (3.0_sp * c(4,i) * t + 2.0_dp * c(3,i)) * t + c(2,i)
    ENDDO
    
    ansy1 = ansy1 / d1
    ansy2 = ansy2 / d2
    
  END SUBROUTINE bicubic_interpolate
  
  !**********************************************************************!
  
  !---------------------------------------------------------------------------
  !
  !  SUBROUTINE tricubic_interpolate
  !
  !> \brief Get interpolated values from bicubic interpolation
  !> \details  Bicubic interpolation within a grid cube,
  !> Input quantities are y,y1,y2,y3,y12,y13,y23,y123,
  !> x1l and x1u, the lower and upper coordinates of the grid square
  !> in the 1- direction; x2l, x2u and x3l and x3u likewise for the
  !> 2 and 3-direction; and x1,x2,x3 the coordinates of the
  !> desired point for the interpolation, The interpolated function value
  !> is returned as ansy, and the interpolated gradient values as
  !> ansy1 and ansy2,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
!!$  SUBROUTINE tricubic_interpolation(y,y1,y2,y3,y12,y13,y23,y123,&
!!$       x1l,x1u,x2l,x2u,x3l,x3u,&
!!$       x1,x2,x3,ansy,ansy1,ansy2,ansy3)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    REAL(dp), INTENT(IN) :: y(8)
!!$    REAL(dp), INTENT(IN) :: y1(8),y2(8),y3(8)
!!$    REAL(dp), INTENT(IN) :: y12(8), y13(8), y23(8)
!!$    REAL(dp), INTENT(IN) :: y123(8)
!!$    REAL(dp), INTENT(IN) :: x1l,x1u,x2l,x2u,x3l,x3u
!!$    REAL(dp), INTENT(IN) :: x1,x2,x3
!!$    REAL(dp), INTENT(OUT):: ansy,ansy1,ansy2,ansy3
!!$    
!!$    REAL(dp) :: t, u, v
!!$    REAL(dp) :: d1, d2, d3
!!$    REAL(dp) :: c(8,8)
!!$    INTEGER  :: i
!!$    
!!$    !---------------------------------------------------------------!
!!$    
!!$    d1 = x1u - x1l
!!$    d2 = x2u - x2l
!!$    d3 = x3u - x3l
!!$    
!!$    IF ( (d1.EQ.0.0_dp).OR.(d2.EQ.0.0_dp).OR.(d3.EQ.0.0_dp)) THEN
!!$       WRITE(*,*) 'Problem with input values in bicubic interpolate,'
!!$       WRITE(*,*) 'Are the boundary pair equal?'
!!$       STOP
!!$    ENDIF
!!$    
!!$    CALL get_tricubic_coeffs(y,y1,y2,y3,y12,y13,y23,y123,&
!!$         d1,d2,d3,wtri,c)
!!$    
!!$    t = (x1-x1l) / d1
!!$    u = (x2-x2l) / d2
!!$    v = (x3-x3l) / d3
!!$    
!!$    ansy = 0.0_dp
!!$    ansy2 = 0.0_dp
!!$    ansy1 = 0.0_dp
!!$    
!!$    DO i = 4, 1, -1
!!$       ansy = t * ansy + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
!!$       ansy2 = t * ansy2 + (3.0_sp * c(i,4) * u + 2.0_dp * c(i,3)) * u + c(i,2)
!!$       ansy1 = u * ansy1 + (3.0_sp * c(4,i) * t + 2.0_dp * c(3,i)) * t + c(2,i)
!!$    ENDDO
!!$    
!!$    ansy1 = ansy1 / d1
!!$    ansy2 = ansy2 / d2
!!$    
!!$  END SUBROUTINE tricubic_interpolation
  
  !**********************************************************************!

  SUBROUTINE delete_bicubic_matrix()
    
    IMPLICIT NONE
    
    DEALLOCATE(wbi)
    
  END SUBROUTINE delete_bicubic_matrix
  
  !********************************************************************!
  
  SUBROUTINE delete_tricubic_matrix
    
    IMPLICIT NONE
    
    DEALLOCATE(wtri)
    
  END SUBROUTINE delete_tricubic_matrix
  
  !********************************************************************!
  !*********************************************!
END MODULE cubic_interp

