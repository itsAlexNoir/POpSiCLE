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
  USE sheppack

  IMPLICIT NONE
  
  ! Original array axis
  REAL(dp), ALLOCATABLE        :: x(:, :), x1(:), x2(:)
  ! Axis for interpolated function
  REAL(dp), ALLOCATABLE        :: y(:), y1(:), y2(:)
  ! Original tabulated values array for function f
  REAL(dp), ALLOCATABLE        :: fx(:), fyy(:, :)
  ! Intepolated values array
  REAL(dp), ALLOCATABLE        :: fy1(:, :), fy2(:, :)
  REAL(dp), ALLOCATABLE        :: fy3(:, :), fy4(:, :)
  ! Polar axes
  REAL(dp), ALLOCATABLE        :: r_non(:, :), theta_non(:, :)
  REAL(dp), ALLOCATABLE        :: r(:), theta(:)
  ! Array to store the derivates at the 
  !interpolated point
  REAL(dp), ALLOCATABLE        :: yp(:)
  ! Extent for x axis
  REAL(dp)                     :: x1min,x1max,x2min,x2max
  ! Extent for y axis
  REAL(dp)                     :: y1min,y1max,y2min,y2max
  REAL(dp)                     :: y1half, y2half
  ! Extent for polar coords
  REAL(dp)                     :: radmax, radmin
  ! Gridspacing for y
  REAL(dp)                     :: dy1, dy2
  ! Gridspacing for polar coords
  REAL(dp)                     :: dr, dtheta
  ! Number of array points
  INTEGER                      :: numxpts, numx1pts, numx2pts
  INTEGER                      :: numypts, numy1pts, numy2pts
  INTEGER                      :: numrpts, numthetapts
  ! QSHEP2D control parameters
  INTEGER                      :: NQ, NR, NW
  REAL(dp)                     :: DX, DY, RMAX, XMIN, YMIN
  REAL(dp)                     :: Q, QX, QY, EPS, EQ, RQ
  INTEGER, ALLOCATABLE         :: LNEXT(:)
  INTEGER, ALLOCATABLE         :: LCELL(:, :)
  REAL(dp), ALLOCATABLE        :: RSQ(:), RW(:)
  REAL(dp), ALLOCATABLE        :: AQ(:, :), AL(:, :)
  ! Error estimation
  REAL(dp)                     :: errory1
  ! Auxiliary integers
  INTEGER                      :: i, j, ierror
  ! Reals for timing
  REAL(dp)                     :: start_time, end_time
  REAL(dp)                     :: qshep2d_time, lshep2d_time
  REAL(dp)                     :: lrshep2d_time, ripple2d_time
  LOGICAL                      :: trial_function
  
  !--------------------------------------------------------!
  
  WRITE(*,*) '!--------------------------------!'
  WRITE(*,*)
  WRITE(*,*) ' interp_scattered2d_test program'
  WRITE(*,*)
  WRITE(*,*) '!--------------------------------!'
  WRITE(*,*)
  WRITE(*,*)  'Starting....'

  !-------------------------------------------------------!
  
  trial_function = .FALSE.

  ! Set number of points
  !numx1pts = 10
  !numx2pts = 10
  !numxpts = numx1pts * numx2pts
  numxpts = 18000
  numy1pts = 80
  numy2pts = 200
  numypts = numy1pts * numy2pts
  
  WRITE(*,*)
  !WRITE(*,*) 'Number of points in x1: ',numx1pts
  !WRITE(*,*) 'Number of points in x2: ',numx2pts
  WRITE(*,*) 'Number of total points in x: ',numxpts
  WRITE(*,*) 'Number of points in y1: ',numy1pts
  WRITE(*,*) 'Number of points in y2: ',numy2pts
  WRITE(*,*) 'Number of total points in y: ',numypts
    
  ! Set extent of the original axis
  x1min = 0.0_dp
  x1max = 50.0_dp
  x2min = 0.0_dp
  x2max = 50.0_dp
  
  y1min = 0.0_dp
  y1max = 50.0_dp
  y2min = 0.0_dp
  y2max = 50.0_dp
  
  ! Assign control parameters for interpolation
  NR = 3
  NQ = 13
  NW = 19
  
  WRITE(*,*)
  WRITE(*,*) 'Allocating arrays...'
  
  
  ALLOCATE(y(2))
  ALLOCATE(y1(numy1pts))
  ALLOCATE(y2(numy2pts))
  ALLOCATE(fyy(numy1pts,numy2pts))
 
  
  ! Build the axis for interp func.
  ! This will be a regular grid.
  
  ! We will make the points equidistant
  ! at the interpolated grid
  dy1 = (y1max - y1min) / numy1pts
  dy2 = (y2max - y2min) / numy2pts
  
  DO i = 1, numy1pts
     y1(i) = (i-1) * dy1      
  ENDDO
  
  DO i = 1, numy2pts
     y2(i) = (i-1) * dy2
  ENDDO
  
 
  IF (trial_function) THEN
     
     ALLOCATE(x(2,numxpts))
     ALLOCATE(x1(numxpts))
     ALLOCATE(x2(numxpts))
     ALLOCATE(fx(numxpts))
     ALLOCATE(fy1(numy1pts,numy2pts)) 
     ALLOCATE(fy2(numy1pts,numy2pts))
     ALLOCATE(fy3(numy1pts,numy2pts))
     ALLOCATE(fy4(numy1pts,numy2pts))
     
     ALLOCATE(LNEXT(numxpts))
     ALLOCATE(LCELL(NR, NR))
     ALLOCATE(RSQ(numxpts))
     ALLOCATE(RW(numxpts))
     ALLOCATE(AQ(5,numxpts))
     ALLOCATE(AL(2,numxpts))
     
     WRITE(*,*) 'Building axes...'
     ! Build original axis, x
     ! We will use the intrinsic function random_number
     ! for this purpose.
     CALL random_seed()
     
     CALL random_number(x1)
     CALL random_number(x2)
     x1 = x1min + x1 * (x1max - x1min)
     x2 = x2min + x2 * (x2max - x2min)
     
     ! Tabulated values for the original function
     fx = SIN(x1) * SIN(x2)
     
     DO j = 1, numy2pts
        DO i = 1, numy1pts
           fyy(i, j) = SIN(y1(i)) * SIN(y2(j))
        ENDDO
     ENDDO
     
     ! Compute the machine precision EPS.
     EPS = EPSILON (EPS)
     EQ = 0.0_dp
     RQ = 0.0_dp
     
     
     !--------------------------------!
     ! Quadratic sheppard algorithm   !
     !--------------------------------!
     
     WRITE(*,*) 'Calculating quadratic sheppack algorithm in 2D...'
     ! Perform Quadratic sheppard algorithm
     CALL cpu_time(start_time)
     
     WRITE(*,*) 'Constructing the interpolant...'
     CALL QSHEP2 (numxpts, x1, x2, fx, NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
          DX, DY, RMAX, RSQ, AQ, ierror)
     IF (ierror /= 0) THEN
        WRITE(*,*) 'Error calculating the interpolant at QSHEP2.'
        STOP
     ENDIF
     
     WRITE(*,*) 'Interpolating the points...'
     DO i = 1, numy1pts
        DO j = 1, numy2pts
           CALL QS2GRD (y1(i), y2(j), numxpts, x1, x2, fx, NR, LCELL, LNEXT, XMIN, &
                YMIN, DX, DY, RMAX, RSQ, AQ, Q, QX, QY, ierror)
           IF (ierror /= 0) THEN
              WRITE (*, *) 'QSHEP2 - ERROR!'
              WRITE (*, *) 'Error in QS2GRD, IER = ', ierror
              STOP
           END IF
           fy1(i,j) = Q
           EQ = MAX (EQ, ABS (fyy(i,j) - Q))
        ENDDO
     ENDDO
     
     CALL cpu_time(end_time)
     
     WRITE(*,*) 'End of QSHEP2D calculation.'
        
     qshep2d_time = end_time - start_time
     
     ! Calculate the ratio EQ / EPS
     RQ = EQ / EPS
        
     !WRITE(*,*) 'Sine function: '
     WRITE(*,*) 'Quadratic sheppack 2D time (seconds): ', qshep2d_time
     WRITE(*,*) 'Max error and Max error/EPS: ',EQ, RQ
     WRITE(*,*)
     
     
     !-----------------------------!
     ! Linear sheppard algorithm   !
     !-----------------------------!
     
     WRITE(*,*) 'Calculating linear sheppack algorithm in 2D...'
     ! Join random points in one matrix
     x(1,:) = x1
     x(2,:) = x2
     
     ! Set error back to zero
     EQ = 0.0
     RQ = 0.0
     
     ! Perform linear sheppard algorithm
     CALL cpu_time(start_time)
     
     WRITE(*,*) 'Constructing the interpolant...'
     CALL LSHEP ( 2, numxpts, x, fx, AL, RW, ierror )
     IF (ierror /= 0) THEN
        WRITE(*,*) 'Error calculating the interpolant at LSHEP.'
        STOP
     ENDIF
     
     DO j = 1, numy2pts
        DO i = 1, numy1pts
           fy2(i,j) = LSHEPVAL( (/ y1(i), y2(j)/), 2, numxpts, x, fx, AL, RW, ierror )
           !WRITE (*, *) 'LSHEPVAL status output = ', ierror
           EQ = EQ + (fy2(i,j) - fyy(i,j))**2
           RQ = RQ + ABS(fy2(i,j) - fyy(i,j))
        ENDDO
     ENDDO
     
     CALL cpu_time(end_time)
     
     WRITE(*,*) 'End of LSHEP calculation.'
     
     lshep2d_time = end_time - start_time
     
     WRITE(*,*) 'Linear sheppack 2D time (seconds): ', lshep2d_time
     WRITE(*,*) 'Error: Root mean: ',SQRT(EQ / numypts)
     WRITE(*,*) 'Error: Square root: ', RQ / numypts
     WRITE(*,*)
     
     
!!$  !------------------------------------!
!!$  ! Linear robust sheppard algorithm   !
!!$  !------------------------------------!
!!$  
!!$  ! Set error back to zero
!!$  EQ = 0.0
!!$  RQ = 0.0
!!$  AL = 0.0
!!$  RW = 0.0
!!$  
!!$  ! Perform linear robust sheppard algorithm
!!$  CALL cpu_time(start_time)
!!$  
!!$  WRITE(*,*) 'Constructing the interpolant...'
!!$  CALL LSHEP ( 2, numxpts, x, fx, AL, RW, ierror, .TRUE. )
!!$  IF (ierror /= 0) THEN
!!$     WRITE(*,*) 'Error calculating the interpolant at LSHEP.',ierror
!!$     STOP
!!$  ENDIF
!!$  
!!$  DO j = 1, numy2pts
!!$     DO i = 1, numy1pts
!!$        fy3(i,j) = LSHEPVAL( (/ y1(i), y2(j)/), 2, numxpts, x, fx, AL, RW, ierror )
!!$        !WRITE (*, *) 'LSHEPVAL status output = ', ierror
!!$        EQ = EQ + (fy3(i,j) - fyy(i,j))**2
!!$        RQ = RQ + ABS(fy3(i,j) - fyy(i,j))
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  CALL cpu_time(end_time)
!!$  
!!$  WRITE(*,*) 'End of LSHEP calculation.'
!!$  
!!$  lrshep2d_time = end_time - start_time
!!$  
!!$  WRITE(*,*) 'Linear robust sheppack 2D time (seconds): ', lrshep2d_time
!!$  WRITE(*,*) 'Error: Root mean: ',SQRT(EQ / numypts)
!!$  WRITE(*,*) 'Error: Square root: ', RQ / numypts
!!$  WRITE(*,*)
        
     !------------------------------------!
     ! Linear sheppard (ripple) algorithm !
     !------------------------------------!
     
     ! Set error back to zero
     EQ = 0.0
     RQ = 0.0
     AL = 0.0
     RW = 0.0
     
     ! Perform linear robust sheppard algorithm
     CALL cpu_time(start_time)
     
     WRITE(*,*) 'Constructing the interpolant...'
     CALL RIPPLE( 2, numxpts, x, fx, AL, RW, ierror )
     IF (ierror /= 0) THEN
        WRITE(*,*) 'Error calculating the interpolant at RIPPLE.',ierror
        STOP
     ENDIF
     
     DO j = 1, numy2pts
        DO i = 1, numy1pts
           fy4(i,j) = LSHEPVAL( (/ y1(i), y2(j)/), 2, numxpts, x, fx, AL, RW, ierror )
           !WRITE (*, *) 'RIPPLE status output = ', ierror
           EQ = EQ + (fy4(i,j) - fyy(i,j))**2
           RQ = RQ + ABS(fy4(i,j) - fyy(i,j))
        ENDDO
     ENDDO
     
     CALL cpu_time(end_time)
     
     WRITE(*,*) 'End of RIPPLE calculation.'
     
     ripple2d_time = end_time - start_time
     
     WRITE(*,*) 'Linear sheppack (ripple) 2D time (seconds): ', ripple2d_time
     WRITE(*,*) 'Error: Root mean: ',SQRT(EQ / numypts)
     WRITE(*,*) 'Error: Square root: ', RQ / numypts
     WRITE(*,*)
     
     !-------------------!
     ! Save to a file
     !-------------------!
     
     OPEN(UNIT=101,FORM='formatted',FILE='results/y1.dat')
     OPEN(UNIT=103,FORM='formatted',FILE='results/y2.dat')
     OPEN(UNIT=105,FORM='formatted',FILE='results/scattered_points.dat')
     
     DO i = 1, numy1pts
        WRITE(101,*) y1(i)
     ENDDO
     
     DO i = 1, numy2pts
        WRITE(103,*) y2(i)
     ENDDO
     
     DO i = 1, numxpts
        WRITE(105,*) x1(i), x2(i)
     ENDDO
     
     CLOSE(101)
     CLOSE(103)
     CLOSE(105)
     
     OPEN(UNIT=101,FORM='formatted',FILE='results/sin2d.dat')
     
     DO j = 1, numy2pts
        DO i = 1, numy1pts
           WRITE(101,*) fyy(i,j), fy1(i,j), fy2(i,j), &
                fy3(i,j), fy4(i,j)
        ENDDO
     ENDDO
     
     CLOSE(101)
     
  ENDIF
  !********************************************************!
     
  !---------------------------------------!
  !  Second part: From cartesian to
  !   polar coordinates
  !---------------------------------------!
  
  WRITE(*,*) 'From cartesian to polar...'

  numrpts = 200
  numthetapts = 300
  
  radmax = 50.0_dp
  radmin = 30.0_dp
  
  ALLOCATE(r_non(numy1pts, numy2pts))
  ALLOCATE(theta_non(numy1pts, numy2pts))
  
  ALLOCATE(r(numrpts))
  ALLOCATE(theta(numthetapts))
  ALLOCATE(fy1(numrpts,numthetapts))
  ALLOCATE(LNEXT(numypts))
  ALLOCATE(LCELL(NR, NR))
  ALLOCATE(RSQ(numypts))
  ALLOCATE(RW(numypts))
  ALLOCATE(AQ(5,numypts))
  ALLOCATE(AL(2,numypts))
  
  DO j = 1, numy2pts
     DO i = 1, numy1pts
        r_non(i,j) = SQRT(y1(i)**2 + y2(j)**2)
        theta_non(i,j) = ATAN2(y2(j), y1(i))
     ENDDO
  ENDDO
  
  dr = (radmax - radmin) / numrpts
  dtheta = 3.14159_dp / 2.0_dp / numthetapts
  
  DO i = 1, numrpts
     r(i) = radmin + i * dr
  ENDDO
  
  DO i = 1, numthetapts
     theta(i) = (i+1) * dtheta
  ENDDO
  
  y1half = 0.0 ! (y1max - y1min) / 2.0 
  y2half = 0.0! (y2max - y2min) / 2.0
  
  DO j = 1, numy2pts
     DO i = 1, numy1pts
        fyy(i,j) = EXP(- ((y1(i) - y1half)/10.)**2 - ((y2(j) - y2half)/10.)**2 )
     ENDDO
  ENDDO
  
  ! Compute the machine precision EPS.
  EPS = EPSILON (EPS)
  EQ = 0.0_dp
  RQ = 0.0_dp
  
  !--------------------------------!
  ! Quadratic sheppard algorithm   !
  !--------------------------------!
  
  WRITE(*,*) 'Calculating quadratic sheppack algorithm in 2D...'
  ! Perform Quadratic sheppard algorithm
  CALL cpu_time(start_time)
  
  WRITE(*,*) 'Constructing the interpolant...'
  CALL QSHEP2 (numypts,RESHAPE(r_non,(/ numypts/)), RESHAPE(theta_non,(/numypts/)), &
       RESHAPE(fyy,(/numypts/)),NQ, NW, NR, LCELL, LNEXT, XMIN, YMIN, &
       DX, DY, RMAX, RSQ, AQ, ierror)
  IF (ierror /= 0) THEN
     WRITE(*,*) 'Error calculating the interpolant at QSHEP2.'
     STOP
  ENDIF
  
  WRITE(*,*) 'Interpolating the points...'
  DO j = 1, numthetapts
     DO i = 1, numrpts
        CALL QS2GRD (r(i), theta(j), numypts, &
             RESHAPE(r_non,(/numypts/)), RESHAPE(theta_non,(/numypts/)), &
             RESHAPE(fyy,(/numypts/)), &
             NR, LCELL, LNEXT, XMIN, &
             YMIN, DX, DY, RMAX, RSQ, AQ, Q, QX, QY, ierror)
        IF (ierror /= 0) THEN
           WRITE (*, *) 'QSHEP2 - ERROR!'
           WRITE (*, *) 'Error in QS2GRD, IER = ', ierror
           !STOP
        END IF
        fy1(i,j) = Q
        EQ = MAX (EQ, ABS (fyy(i,j) - Q))
     ENDDO
  ENDDO
  
  CALL cpu_time(end_time)
  
  WRITE(*,*) 'End of QSHEP2D calculation.'
  
  qshep2d_time = end_time - start_time
  
  ! Calculate the ratio EQ / EPS
  RQ = EQ / EPS
  
  !WRITE(*,*) 'Sine function: '
  WRITE(*,*) 'Quadratic sheppack 2D time (seconds): ', qshep2d_time
  WRITE(*,*) 'Max error and Max error/EPS: ',EQ, RQ
  WRITE(*,*)
  

  !-------------------!
  ! Save to a file
  !-------------------!
  
  IF (trial_function .EQ. .FALSE.) THEN
     OPEN(UNIT=101,FORM='formatted',FILE='results/y1.dat')
     OPEN(UNIT=103,FORM='formatted',FILE='results/y2.dat')
     
     
     DO i = 1, numy1pts
        WRITE(101,*) y1(i)
     ENDDO
     
     DO i = 1, numy2pts
        WRITE(103,*) y2(i)
     ENDDO
     
  ENDIF
  
  OPEN(UNIT=105,FORM='formatted',FILE='results/radius.dat')
  OPEN(UNIT=107,FORM='formatted',FILE='results/theta.dat')
  
  DO i = 1, numrpts
     WRITE(105,*) r(i)
  ENDDO
  
  DO i = 1, numthetapts
     WRITE(107,*) theta(i)
  ENDDO
  
  
  OPEN(UNIT=109,FORM='formatted',FILE='results/pol_coords.dat')
  DO j = 1, numy2pts
     DO i = 1, numy1pts
        WRITE(109,*) r_non(i,j), theta_non(i,j),fyy(i,j)
     ENDDO
  ENDDO
  
  CLOSE(101)
  CLOSE(103)
  CLOSE(105)
  CLOSE(107)
  CLOSE(109)
  
  
  
  OPEN(UNIT=101,FORM='formatted',FILE='results/cart2pol.dat')
  
  DO j = 1, numthetapts
     DO i = 1, numrpts
        WRITE(101,*) fy1(i,j)
     ENDDO
  ENDDO
  
  CLOSE(101)
  
  
  
END PROGRAM interp1d_test
