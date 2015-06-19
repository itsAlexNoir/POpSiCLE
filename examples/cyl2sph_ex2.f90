PROGRAM cyl2sph_ex2

  USE popsicle
  
  IMPLICIT NONE

  INTEGER                    :: numrhopts
  INTEGER                    :: numzpts
  INTEGER                    :: numrpts
  INTEGER                    :: numthetapts
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  INTEGER                    :: dims(2)
  
  COMPLEX(dp), ALLOCATABLE   :: Y_lm(:, :)
  COMPLEX(dp), ALLOCATABLE   :: R_nl(:, :)
  COMPLEX(dp), ALLOCATABLE   :: cylfunc(:, :)
  REAL(dp)                   :: drho, dz
  COMPLEX(dp), ALLOCATABLE   :: sphfunc(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dr(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dth(:, :)
  COMPLEX(dp)                :: ref_value
  REAL(dp)                   :: dr, dtheta
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: thetapt, rpt
  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: maxerror
  INTEGER                    :: irho,iz
  INTEGER                    :: ir, itheta
  
  ! Set number of points
  numrhopts = 200
  numzpts   = 500
  
  ! Set grid spacing
  drho = 0.01
  dz   = 0.02
  
  ! Create axes
  ALLOCATE(rho_ax(1:numrhopts))
  ALLOCATE(z_ax(1:numzpts))
  
  DO irho =1, numrhopts
     rho_ax(irho) = REAL(irho,dp) * drho
  ENDDO
  
  DO iz = 1, numzpts
     z_ax(iz) = REAL(-numzpts / 2 + iz,dp) * dz  
  ENDDO
  
  
  ! We are going to set the 3d hydrogen function
  ! n=3, l=2, m=0
  ALLOCATE(Y_lm(1:numrhopts,1:numzpts))
  ALLOCATE(R_nl(1:numrhopts,1:numzpts))
  ALLOCATE(cylfunc(1:numrhopts,1:numzpts))
  
  
  DO iz = 1, numzpts
     DO irho = 1, numrhopts
        
        rpt = rho_ax(irho)**2 + z_ax(iz)**2
        thetapt = atan2(z_ax(iz),rho_ax(irho))
        
        Y_lm(irho,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
             ( 3.0_dp * (COS(thetapt))**2 - 1.0_dp)
        
        R_nl(irho,iz) = 1.0_dp
        
        cylfunc(irho,iz) = EXP(-rpt / 2.0_dp) * R_nl(irho,iz) * &
             Y_lm(irho,iz)
        
     ENDDO
  ENDDO
  
  ! Initialize the boundary
  Rboundary = 1.0_dp
  tolerance = 0.05_dp
  dr = 0.01_dp
  dtheta = 0.05_dp
  dims      = (/numrhopts, numzpts/)
  
  ! Build interpolant
  WRITE(*,*) 'Creating interpolant...'
  CALL cpu_time(start_time)
  
  CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
       Rboundary, tolerance, 2, dr, dtheta, numrpts, numthetapts )
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'Interpolant time (seconds): ', interp_time
  WRITE(*,*) 'Number of radial boundary points: ',numrpts
  WRITE(*,*) 'Number of azimuthal boundary points: ',numthetapts
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts))
  
  ! Interpolate!!
  WRITE(*,*) 'Interpolating boundary...'
  CALL cpu_time(start_time)
  
  CALL get_cylindrical_boundary(cylfunc, sphfunc, sphfunc_dr, &
       sphfunc_dth, 'quadratic')
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'Interpolation time (seconds): ', interp_time
  
  
  DO itheta = 1, numthetapts
     DO ir = 1, numrpts
        ref_value =  EXP(-rpts(ir) / 2.0_dp) *  0.25_dp * SQRT(5.0_dp / pi) * &
             ( 3.0_dp * (COS(theta(itheta)))**3 - 1.0_dp)
!!$        WRITE(*,*) sphfunc(ir,itheta)
!!$        WRITE(*,*) ref_value 
!!$        WRITE(*,*) 'Diff: ', ABS(sphfunc(ir,itheta) - ref_value)
        maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta) - ref_value))
     ENDDO
  ENDDO
  
  WRITE(*,*) 
  WRITE(*,*) '********RESULTS!!!******'
  WRITE(*,*)
  WRITE(*,*) 'Maximum difference in value: '
  WRITE(*,*) maxerror
  
  ! Save original wavefunction
  OPEN(unit=33,form='formatted',file='cylfunction.dat')
  
  DO iz = 1, numzpts
     DO irho = 1, numrhopts
        WRITE(33,*) ABS(cylfunc(irho,iz))**2 
     ENDDO
  ENDDO
  
  CLOSE(33)
  
  ! Save spherical wavefunction
  OPEN(unit=33,form='formatted',file='sphfunction.dat')
  
  DO itheta = 1, numthetapts
     DO ir = 1, numrpts
        WRITE(33,*) ABS(sphfunc(ir,itheta))**2 
     ENDDO
  ENDDO
  
  CLOSE(33)


  ! Free memory
  DEALLOCATE(rho_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cylfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth)

END PROGRAM
