PROGRAM cart2sph_ex1

  USE popsicle
  
  IMPLICIT NONE

  INTEGER                    :: numxpts
  INTEGER                    :: numypts  
  INTEGER                    :: numzpts
  INTEGER                    :: numrpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts
  INTEGER                    :: numpts
  REAL(dp), ALLOCATABLE      :: x_ax(:)
  REAL(dp), ALLOCATABLE      :: y_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  INTEGER                    :: dims(3)
  
  COMPLEX(dp), ALLOCATABLE   :: Y_lm(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: R_nl(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: cartfunc(:, :, :)

  REAL(dp)                   :: dx, dy, dz
  COMPLEX(dp), ALLOCATABLE   :: sphfunc(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dth(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dphi(:, :, :)
  COMPLEX(dp)                :: ref_value

  REAL(dp)                   :: dr, dtheta, deltaphi
  INTEGER                    :: lmax
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: rpt, thetapt, phipt
  REAL(dp)                   :: efield(3), afield(3)  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: minerror, maxerror
  INTEGER                    :: ix,iy,iz
  INTEGER                    :: ir, itheta, iphi
  INTEGER                    :: iter
  CHARACTER(LEN = 7)         :: rbstr
  CHARACTER(LEN = 5)         :: lmaxstr
  CHARACTER(LEN = 100)       :: filename
  
  !------------------------------------------------------!
  
  WRITE(*,*) '****************************'
  WRITE(*,*) '     cart2sph_ex1           '
  WRITE(*,*) '****************************'
  WRITE(*,*)
  WRITE(*,*)

  ! Set number of points
  numxpts = 251 !191
  numypts = 251 !191
  numzpts = 251 !191
  
  ! Set grid spacing
  dx   = 0.1_dp
  dy   = 0.1_dp
  dz   = 0.1_dp
  
  ! Create axes
  ALLOCATE(x_ax(1:numxpts))
  ALLOCATE(y_ax(1:numypts))
  ALLOCATE(z_ax(1:numzpts))
  
  DO ix = 1, numxpts
     x_ax(ix) = REAL(-numxpts / 2 + ix,dp) * dx
  ENDDO

  DO iy = 1, numxpts
     y_ax(iy) = REAL(-numypts / 2 + iy,dp) * dy
  ENDDO
  
  DO iz = 1, numzpts
     z_ax(iz) = REAL(-numzpts / 2 + iz,dp) * dz  
  ENDDO
  
  
  ! We are going to set the 3d hydrogen function
  ! n=3, l=2, m=0
  ALLOCATE(Y_lm(1:numxpts,1:numypts,1:numzpts))
  ALLOCATE(R_nl(1:numxpts,1:numypts,1:numzpts))
  ALLOCATE(cartfunc(1:numxpts,1:numypts,1:numzpts))
  
  DO iz = 1, numzpts
     DO iy = 1, numypts  
        DO ix = 1, numxpts
           
           rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
           
           IF (rpt.EQ.0) THEN
              thetapt = pi / 2.0_dp
           ELSE
              thetapt = ACOS( z_ax(iz) / rpt )
           ENDIF
           
           phipt = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
           
           ! l=2, m=0
           !Y_lm(ix,iy,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
           !     ( 3.0_dp * COS(thetapt)**2 - 1.0_dp)
           
           ! l=2, m=1
           Y_lm(ix,iy,iz) = -0.5_dp * SQRT(15.0_dp / twopi) * &
                SIN(thetapt) * COS(thetapt) * EXP(ZIMAGONE*phipt)
           
           R_nl(ix,iy,iz) = 1.0_dp
           
           cartfunc(ix,iy,iz) = EXP(-rpt / 2.0_dp) * R_nl(ix,iy,iz) * &
                Y_lm(ix,iy,iz)
           
        ENDDO
     ENDDO
  ENDDO
  
  
  ! Initialize the boundary
  Rboundary = 9.0_dp !7.0_dp
  tolerance = 0.15_dp
  dr = 0.1_dp
  lmax = 10
  dtheta = 0.1_dp
  deltaphi = 0.1_dp
  dims      = (/numxpts, numypts, numzpts/)
  
  WRITE(rbstr,'(I3.3,F0.3)') INT(Rboundary),Rboundary-INT(Rboundary)
  WRITE(lmaxstr,'(I3.3)') lmax
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  
  WRITE(*,*) 'Number of points in x: ',numxpts
  WRITE(*,*) 'Number of points in y: ',numypts
  WRITE(*,*) 'Number of points in z: ',numzpts
  
  WRITE(*,*) 'Grid spacing in x: ',dx
  WRITE(*,*) 'Grid spacing in y: ',dy
  WRITE(*,*) 'Grid spacing in z: ',dz
  
  WRITE(*,*) 'Grid extent in x: ',x_ax(numxpts)
  WRITE(*,*) 'Grid extent in y: ',y_ax(numypts)
  WRITE(*,*) 'Grid extent in z: ',z_ax(numzpts)
  
  WRITE(*,*) 'Boundary at radius: ',Rboundary
  WRITE(*,*) 'Radius tolerance: ',tolerance
  WRITE(*,*)
  WRITE(*,*) '--------------------------'
  
  ! Build interpolant
  WRITE(*,*) 'Creating interpolant...'
  CALL cpu_time(start_time)
  
  !! For scattered interpolation uncomment the subroutine below
  !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
  !     Rboundary, tolerance, 2, dr, lmax, &
  !     numpts, numrpts, numthetapts, numphipts )
  
  CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
       Rboundary, 2, dr, lmax, numrpts, numthetapts, numphipts)

  !CALL initialize_cartesian_surface(x_ax, y_ax, z_ax, dims, &
  !     Rboundary, tolerance, 2, dr, lmax, .TRUE., filename)
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'Interpolant time (seconds): ', interp_time
  WRITE(*,*) 
  
  WRITE(*,*) 'Total number of points to be interpolated: ',numpts
  
  WRITE(*,*) 'Number of radial boundary points: ',numrpts
  WRITE(*,*) 'Number of polar boundary points: ',numthetapts
  WRITE(*,*) 'Number of azimuthal boundary points: ',numphipts
  
  WRITE(*,*) 'Grid spacing in r: ',dr
  WRITE(*,*) 'Maximum angular momenta: ',lmax
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '--------------------------'
  
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dphi(1:numrpts,1:numthetapts,1:numphipts))

  efield = 0.0_dp
  afield = 0.0_dp
  
  ! Interpolate!!
  WRITE(*,*) 'Interpolating boundary...'
  DO iter = 1, 1!0
     CALL cpu_time(start_time)
     
     !! For scattered interpolation uncomment the subroutine below
     !CALL get_cartesian_boundary(cartfunc, sphfunc, sphfunc_dr, &
     !     sphfunc_dth, sphfunc_dphi, 'quadratic')
     
     CALL get_cartesian_boundary(x_ax, y_ax, z_ax, dims, cartfunc, &
          2, sphfunc, sphfunc_dr, sphfunc_dth, sphfunc_dphi)
     
     !CALL get_cartesian_surface(filename, cartfunc, x_ax, y_ax, z_ax, dims, &
     !     2, 0.0_dp , efield, afield, lmax, .TRUE. )
     
     CALL cpu_time(end_time)
     
     interp_time = end_time - start_time
     
     WRITE(*,*) 'Interpolation time (seconds): ', interp_time
     
  ENDDO
  
  DO iphi = 1, numphipts
     DO itheta = 1, numthetapts
        DO ir = 1, numrpts
           ! For l=2, m=0
!!$           ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  &
!!$                0.25_dp * SQRT(5.0_dp / pi) * &
!!$                ( 3.0_dp * COS(theta_boundary(itheta))**2 - 1.0_dp)

           ! For l=2, m=1
           ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  &
                -0.5_dp * SQRT(15.0_dp / twopi) * &
                SIN(theta_boundary(itheta)) * &
                COS(theta_boundary(itheta)) * &
                EXP(ZIMAGONE * phi_boundary(iphi))
           
           maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta, iphi) - &
                ABS(ref_value)))
           minerror = MIN(minerror, ABS(sphfunc(ir,itheta, iphi) -&
                ABS(ref_value)))
        ENDDO
     ENDDO
  ENDDO
  
  WRITE(*,*) 
  WRITE(*,*) '********RESULTS!!!******'
  WRITE(*,*)
  WRITE(*,*) 'Minimum difference in value: '
  WRITE(*,*) minerror
  WRITE(*,*) 'Maximum difference in value: '
  WRITE(*,*) maxerror
  
!!$  ! Save axes
!!$  WRITE(*,*) 
!!$  WRITE(*,*) 'Write data to disk...'
!!$  
!!$  ! X
!!$  OPEN(unit=33,form='formatted',file='./results/x_ax.dat')
!!$  DO ix = 1, numxpts
!!$     WRITE(33,*) x_ax(ix) 
!!$  ENDDO
!!$  CLOSE(33)
!!$  
!!$  ! Y
!!$  OPEN(unit=33,form='formatted',file='./results/y_ax.dat')
!!$  DO iy = 1, numypts
!!$     WRITE(33,*) y_ax(iy) 
!!$  ENDDO
!!$  CLOSE(33)
!!$
!!$  ! Z
!!$  OPEN(unit=33,form='formatted',file='./results/z_ax.dat')
!!$  DO iz = 1, numzpts
!!$     WRITE(33,*) z_ax(iz)
!!$  ENDDO
!!$  CLOSE(33)
!!$  
!!$  ! R
!!$  OPEN(unit=33,form='formatted',file='./results/r_ax.dat')
!!$  DO ir = 1, numrpts
!!$     WRITE(33,*) rpts_boundary(ir)
!!$  ENDDO
!!$  CLOSE(33)
!!$  
!!$  ! Theta
!!$  OPEN(unit=33,form='formatted',file='./results/theta_ax.dat')
!!$  DO itheta = 1, numthetapts
!!$     WRITE(33,*) theta_boundary(itheta)
!!$  ENDDO
!!$  CLOSE(33)
!!$
!!$  ! Theta
!!$  OPEN(unit=33,form='formatted',file='./results/phi_ax.dat')
!!$  DO iphi = 1, numphipts
!!$     WRITE(33,*) phi_boundary(iphi)
!!$  ENDDO
!!$  CLOSE(33)
!!$
!!$  ! Save original wavefunction
!!$  OPEN(unit=33,form='formatted',file='./results/cartfunction.dat')
!!$  
!!$  DO iz = 1, numzpts
!!$     DO iy = 1, numypts
!!$        DO ix = 1, numxpts
!!$           WRITE(33,*) ABS(cartfunc(ix,iy,iz))**2 
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  CLOSE(33)
!!$  
!!$  ! Save spherical wavefunction
!!$  OPEN(unit=33,form='formatted',file='./results/sphfunction.dat')
!!$
!!$  DO iphi = 1, numphipts
!!$     DO itheta = 1, numthetapts
!!$        DO ir = 1, numrpts
!!$           WRITE(33,*) ABS(sphfunc(ir,itheta,iphi))**2 
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  CLOSE(33)
  
  ! Free memory
  DEALLOCATE(x_ax,y_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cartfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth,sphfunc_dphi)

  WRITE(*,*) '¡Se acabó!'
  WRITE(*,*)
  WRITE(*,*)
  
END PROGRAM cart2sph_ex1
