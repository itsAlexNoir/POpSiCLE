!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! PROGRAM: interp1d_test
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Example of interpolating routines in 1d.
!> \details This program performs a couple of interpolations
!> comparing differents methods for interpolation> In particular
!> polynomial against spline interpolation.
!>
!------------------------------------------------------------------------------

PROGRAM interp1d_test
  
  USE constants
  
  IMPLICIT NONE
  
  ! Original array axis
  REAL(dp), ALLOCATABLE        :: x(:)
  ! Axis for interpolated function
  REAL(dp), ALLOCATABLE        :: y(:)
  ! Information array. It contains info of 
  ! the interpolating polynomial.
  REAL(dp), ALLOCATABLE        :: c(:)
  ! Original tabulated values array for function f
  REAL(dp), ALLOCATABLE        :: fx(:), dfx2(:)
  ! Intepolated values array
  REAL(dp), ALLOCATABLE        :: fy1(:), fy2(:)
  ! Array to store the derivates at the 
  !interpolated point
  REAL(dp), ALLOCATABLE        :: yp(:)
  ! Extent for x axis
  REAL(dp)                     :: xmin,xmax
  ! Extent for y axis
  REAL(dp)                     :: ymin,ymax
  ! Gridspacing for x
  REAL(dp)                     :: dx
  ! Number of array points
  INTEGER                      :: numpts
  ! Order of the derivatives
  INTEGER                      :: nder
  ! Error estimation
  REAL(dp)                     :: errory1
  ! Auxiliary integers
  INTEGER                      :: i, ierror
  ! Auxiliary arrays
  REAL(dp), ALLOCATABLE        :: work(:)
  ! Reals for timing
  REAL(dp)                     :: start_time, end_time
  REAL(dp)                     :: poly_time, spline_time
  
  !--------------------------------------------------------!
  
  ! Set number of points
  numpts = 500
  ! Set order of the derivatives to calculate
  nder = 0

  ! Set extent of the original axis
  xmin = 0.0_dp
  xmax = 50.0_dp
  ymin = 0.0_dp
  ymax = 50.0_dp

  ALLOCATE(x(numpts))
  ALLOCATE(y(numpts))
  ALLOCATE(c(numpts))
  ALLOCATE(fx(numpts))
  ALLOCATE(dfx2(numpts))
  ALLOCATE(fy1(numpts))
  ALLOCATE(fy2(numpts))
  
  IF (nder .EQ. 0) THEN
     ALLOCATE(work(1))
     ALLOCATE(yp(1))
  ELSE
     ALLOCATE(work(2*numpts))
     ALLOCATE(yp(nder))
  ENDIF
  
  ! We will make the original axis equidistant
  dx = (xmax - xmin) / numpts

  ! Build original axis, x
  DO i = 1, numpts
     x(i) = (i-1) * dx
  ENDDO
  
  ! Tabulated values for the original function
  fx = SIN(x)
  
  ! Build the axis for interp func.
  ! We will use the intrinsic function random_number
  ! for this purpose.
  !CALL init_random_seed()
  
  CALL random_number(y)
  y = ymin + y * (ymax - ymin)
  
  ! Perform polynomial interpolation
  CALL cpu_time(start_time)
  
  CALL dplint(numpts,x,fx,c)
  
  DO i = 1, numpts
     CALL dpolvl(nder,y(i),fy1(i),yp,numpts,x,c,work,ierror)
  ENDDO
  
  CALL cpu_time(end_time)
  
  write(*,*) 'val y:',y(1), y(2)
  write(*,*) 'Value: ',fy1(1), sin(y(1))
  write(*,*) 'val2: ',fy1(2),sin(y(2))
  
  poly_time = end_time - start_time
  
  WRITE(*,*) 'Sine function: '
  WRITE(*,*) 'Polynomial interpolation time (seconds): ', poly_time
  !WRITE(*,*) 'Error in polynomial interpolation: ',errory1
  WRITE(*,*)
  
  
  ! Save to a file
  OPEN(UNIT=101,FORM='formatted',FILE='results/interp1d_test_sin_tab.dat')
  OPEN(UNIT=103,FORM='formatted',FILE='results/interp1d_test_sin_poly.dat')
  
  
  DO i = 1, numpts
     WRITE(101,*) x(i), fx(i)
     WRITE(103,*) y(i), fy1(i)
  ENDDO
  
  CLOSE(101)
  CLOSE(103)
  
!!$  !---------------------------------------!
!!$  !  Second part: Gaussian function
!!$  !---------------------------------------!
!!$  
!!$  fx = EXP(-(x - (xmax-xmin)/2.0 )**2)
!!$  
!!$  ! Perform polynomial interpolation
!!$  CALL cpu_time(start_time)
!!$  
!!$  DO i = 1, numpts
!!$     call polint(x,fx,y(i),fy1(i),errory1)
!!$  ENDDO
!!$  
!!$  CALL cpu_time(end_time)
!!$  
!!$  poly_time = end_time-start_time
!!$  
!!$  ! Now, the cubic spline interpolation
!!$  ! Set second derivatives
!!$  call set_dy2(x,fx,1.0E99_dp,1.0E99_dp,dfx2)
!!$  
!!$  CALL cpu_time(start_time)
!!$  DO i = 1, numpts
!!$     fy2(i) = splint(x,fx,dfx2,y(i))
!!$  ENDDO
!!$  CALL cpu_time(end_time)
!!$  
!!$  spline_time = end_time-start_time
!!$  
!!$  WRITE(*,*) 'Gaussian function: '
!!$  WRITE(*,*) 'Polynomial interpolation time (seconds): ',poly_time
!!$  WRITE(*,*) 'Cubic spline interpolation time (seconds): ', spline_time
!!$  WRITE(*,*) 'Error in polynomial interpolation: ',errory1
!!$  WRITE(*,*)
!!$  
!!$  ! Save to a file
!!$  OPEN(UNIT=101,FORM='formatted',FILE='results/interp1d_test_exp_tab.dat')
!!$  OPEN(UNIT=103,FORM='formatted',FILE='results/interp1d_test_exp_poly.dat')
!!$  OPEN(UNIT=105,FORM='formatted',FILE='results/interp1d_test_exp_spline.dat')
!!$  
!!$  
!!$  DO i = 1, numpts
!!$     WRITE(101,*) x(i), fx(i)
!!$     WRITE(103,*) y(i), fy1(i)
!!$     WRITE(105,*) y(i), fy2(i)
!!$  ENDDO
!!$  
!!$  CLOSE(101)
!!$  CLOSE(103)
!!$  CLOSE(105)

  
END PROGRAM interp1d_test
