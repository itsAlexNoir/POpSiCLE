PROGRAM cart2sph_ex2

  USE popsicle
  
  IMPLICIT NONE

  INTEGER                    :: numxpts
  INTEGER                    :: numypts  
  INTEGER                    :: numzpts
  INTEGER                    :: numrpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts
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
  REAL(dp)                   :: dr, dtheta, dphi
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: rpt, thetapt, phipt
  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: maxerror
  INTEGER                    :: ix,iy, iz
  INTEGER                    :: ir, itheta, iphi
  
  ! Set number of points
  numxpts = 200
  numypts = 200
  numzpts = 200
  
  ! Set grid spacing
  dx   = 0.01
  dy   = 0.01
  dz   = 0.01
  
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
           thetapt = ACOS( z_ax(iz) / rpt )
           phipt = ATAN( y_ax(iy) / x_ax(ix) )
           
           Y_lm(ix,iy,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
                ( 3.0_dp * (COS(thetapt))**2 - 1.0_dp)
           
           R_nl(ix,iy,iz) = 1.0_dp
           
           cartfunc(ix,iy,iz) = EXP(-rpt / 2.0_dp) * R_nl(ix,iy,iz) * &
                Y_lm(ix,iy,iz)
           
        ENDDO
     ENDDO
  ENDDO
  
  ! Initialize the boundary
  Rboundary = 1.0_dp
  tolerance = 0.05_dp
  dr = 0.01_dp
  dtheta = 0.05_dp
  dims      = (/numxpts, numypts, numzpts/)
  
  ! Build interpolant
  WRITE(*,*) 'Creating interpolant...'
  CALL cpu_time(start_time)
  
  CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
       Rboundary, tolerance, 2, dr, dtheta, dphi, &
       numrpts, numthetapts, numphipts )
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'Interpolant time (seconds): ', interp_time
  WRITE(*,*) 'Number of radial boundary points: ',numrpts
  WRITE(*,*) 'Number of polar boundary points: ',numthetapts
  WRITE(*,*) 'Number of azimuthal boundary points: ',numphipts
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dphi(1:numrpts,1:numthetapts,1:numphipts))
  
  ! Interpolate!!
  WRITE(*,*) 'Interpolating boundary...'
  CALL cpu_time(start_time)
  
  CALL get_cartesian_boundary(cartfunc, sphfunc, sphfunc_dr, &
       sphfunc_dth, sphfunc_dphi, 'quadratic')
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'Interpolation time (seconds): ', interp_time
  

  DO iphi = 1, numphipts
     DO itheta = 1, numthetapts
        DO ir = 1, numrpts
           ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  0.25_dp * SQRT(5.0_dp / pi) * &
                ( 3.0_dp * (COS(theta_boundary(itheta)))**2 - 1.0_dp)
!!$        WRITE(*,*) sphfunc(ir,itheta)
!!$        WRITE(*,*) ref_value 
!!$        WRITE(*,*) 'Diff: ', ABS(sphfunc(ir,itheta) - ref_value)
           maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta, iphi) - ref_value))
        ENDDO
     ENDDO
  ENDDO
  
  WRITE(*,*) 
  WRITE(*,*) '********RESULTS!!!******'
  WRITE(*,*)
  WRITE(*,*) 'Maximum difference in value: '
  WRITE(*,*) maxerror
  
  ! Save axes
  
  ! X
  OPEN(unit=33,form='formatted',file='./results/x_ax.dat')
  DO ix = 1, numxpts
     WRITE(33,*) x_ax(ix) 
  ENDDO
  CLOSE(33)
  
  ! Y
  OPEN(unit=33,form='formatted',file='./results/y_ax.dat')
  DO iy = 1, numypts
     WRITE(33,*) y_ax(iy) 
  ENDDO
  CLOSE(33)

  ! Z
  OPEN(unit=33,form='formatted',file='./results/z_ax.dat')
  DO iz = 1, numzpts
     WRITE(33,*) z_ax(iz)
  ENDDO
  CLOSE(33)
  
  ! R
  OPEN(unit=33,form='formatted',file='./results/r_ax.dat')
  DO ir = 1, numrpts
     WRITE(33,*) rpts_boundary(ir)
  ENDDO
  CLOSE(33)
  
  ! Theta
  OPEN(unit=33,form='formatted',file='./results/theta_ax.dat')
  DO itheta = 1, numthetapts
     WRITE(33,*) theta_boundary(itheta)
  ENDDO
  CLOSE(33)

  ! Theta
  OPEN(unit=33,form='formatted',file='./results/phi_ax.dat')
  DO iphi = 1, numphipts
     WRITE(33,*) phi_boundary(iphi)
  ENDDO
  CLOSE(33)

  ! Save original wavefunction
  OPEN(unit=33,form='formatted',file='./results/cartfunction.dat')
  
  DO iz = 1, numzpts
     DO iy = 1, numypts
        DO ix = 1, numxpts
           WRITE(33,*) ABS(cartfunc(ix,iy,iz))**2 
        ENDDO
     ENDDO
  ENDDO
  
  CLOSE(33)
  
  ! Save spherical wavefunction
  OPEN(unit=33,form='formatted',file='./results/sphfunction.dat')

  DO iphi = 1, numphipts
     DO itheta = 1, numthetapts
        DO ir = 1, numrpts
           WRITE(33,*) ABS(sphfunc(ir,itheta,iphi))**2 
        ENDDO
     ENDDO
  ENDDO
  
  CLOSE(33)
  
  ! Free memory
  DEALLOCATE(x_ax,y_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cartfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth,sphfunc_dphi)
  
END PROGRAM cart2sph_ex2
