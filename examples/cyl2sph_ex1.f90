PROGRAM cyl2sph_ex1
  
  USE popsicle_aux
  
  IMPLICIT NONE
  
  INTEGER                    :: maxrhopts, maxzpts
  INTEGER                    :: maxrho, minrho
  INTEGER                    :: maxz, minz
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp)                   :: rhoalpha, rhobeta
  REAL(dp)                   :: sqrtrho, sqrt1rho
  REAL(dp)                   :: rhopt, zpt
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
  INTEGER                    :: lmax, fdrule
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: rpt, thetapt
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: minerror, maxerror
  INTEGER                    :: irho,iz
  INTEGER                    :: ir, itheta
  INTEGER                    :: iter
  REAL(dp)                   :: efield(3), afield(3)
  CHARACTER(LEN=100)         :: filename

  
  !------------------------------------------------------------
  
  WRITE(*,*)
  WRITE(*,*) '****************************'
  WRITE(*,*) '       Cyl2sph1 example.'
  WRITE(*,*) '****************************'
  WRITE(*,*) 
  WRITE(*,*) 'Opening...'
  
  ! Set number of points
  
  maxrhopts = 500 !80
  maxzpts = 501 !80
  
  minrho = 1
  maxrho = maxrhopts + 1
  
  minz = 1
  maxz = maxzpts + 1
  
  ! Set grid spacing
  drho   = 0.1_dp
  dz   = 0.1_dp
  
  ! Create axes
  ALLOCATE(rho_ax(minrho:maxrho))
  ALLOCATE(z_ax(minz:maxz))

  rhoalpha = 1.0_dp
  rhobeta = 0.0_dp
  
  DO irho = minrho, maxrho
!!$     rhopt          = REAL(irho,dp) * drho
!!$     sqrtrho	    = SQRT(rhopt)
!!$     sqrt1rho	    = SQRT(rhoalpha + rhobeta * rhopt)
!!$     rho_ax(irho)   = rhopt * sqrtrho / sqrt1rho
     rho_ax(irho)   = REAL(irho-1,dp) * drho
  ENDDO
  
  DO iz = minz, maxz
     z_ax(iz) = REAL(-maxzpts / 2 + iz,dp) * dz
  ENDDO
   
  ! We are going to set the 3d hydrogen function
  ! n=3, l=2
  ALLOCATE(Y_lm(minrho:maxrho,minz:maxz))
  ALLOCATE(R_nl(minrho:maxrho,minz:maxz))
  ALLOCATE(cylfunc(minrho:maxrho,minz:maxz))
  
  DO iz = minz, maxz
     DO irho = minrho, maxrho
        
        rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
        
        IF (rpt.EQ.0) THEN
           thetapt = pi / 2.0_dp
        ELSE
           thetapt = ACOS( z_ax(iz) / rpt )
        ENDIF
        
        ! n=2 l=1 m=0
        !Y_lm(irho,iz) = 0.5_dp * SQRT(3.0_dp / pi ) * COS(thetapt)
        
        ! n=3 l=2 m=0
        Y_lm(irho,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
             ( 3.0_dp * COS(thetapt)**2 - 1.0_dp)
        
        R_nl(irho,iz) = 1.0_dp
        
        cylfunc(irho,iz) = EXP( - rpt / 2.0_dp ) * R_nl(irho,iz) * &
             Y_lm(irho,iz)
     ENDDO
  ENDDO
  
  ! Initialize the boundary
  Rboundary = 20.0_dp
  tolerance = 0.15_dp !0.15_dp
  dr = 0.1_dp
  lmax = 100 !8
  fdrule = 6
  dtheta = 0.1_dp
  dims      = (/maxrhopts, maxzpts/)
  filename    = 'results/cyl2sph'
  
  WRITE(*,*) 'Number of points in rho: ',maxrhopts
  WRITE(*,*) 'Number of points in z: ',maxzpts
  
  WRITE(*,*) 'Grid spacing in rho: ',drho
  WRITE(*,*) 'Grid spacing in z: ',dz

  WRITE(*,*) 'Extent in rho: ',rho_ax(maxrhopts)
  WRITE(*,*) 'Extent in z: ',z_ax(maxzpts)
  
  WRITE(*,*) 'Boundary at radius: ',Rboundary
  WRITE(*,*) 'Radius tolerance: ',tolerance
  WRITE(*,*)
  WRITE(*,*) '--------------------------'
  
  ! Build interpolant
  WRITE(*,*) 'Creating interpolant...'
  
  CALL cpu_time(start_time)
  
  !! For scattered interpolation uncomment the subroutine below
  !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
  !     Rboundary, tolerance, 2, dr, lmax )
  
  !! For scattered interpolation uncomment the subroutine below
  CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
       Rboundary, 2, dr, lmax )
  
!!$  CALL initialize_cylindrical_surface(rho_ax(1:maxrhopts+1), &
!!$       z_ax(1:maxzpts+1),(/maxrhopts+1,maxzpts+1/), Rboundary, &
!!$       tolerance, fdrule, dr, lmax, filename )
    
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  WRITE(*,*)
  WRITE(*,*)           '*********************************************'
  WRITE(*,'(A,F11.8)') 'Interpolant time (seconds):       ', interp_time
  WRITE(*,*)           '*********************************************'  
  WRITE(*,*) 
  
!!$  !WRITE(*,*) 'Total number of points to be interpolated: ',numpts
!!$  
!!$  WRITE(*,'(A,I3)')   'Number of radial boundary points: ',numrpts
!!$  WRITE(*,'(A,I3)')   'Number of polar boundary points:  ',numthetapts
!!$  
!!$  WRITE(*,'(A,F9.3)') 'Grid spacing in r:                ',dr
!!$  WRITE(*,'(A,I3)')   'Maximum angular momenta:          ',lmax
!!$  WRITE(*,*)
!!$  WRITE(*,*) '--------------------------'
     
     
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts))
  
  ! Interpolate!!
  WRITE(*,*) 'Interpolating boundary...'

  efield = (/0.0_dp, 0.0_dp, 0.0_dp /)
  afield = (/0.0_dp, 0.0_dp, 0.0_dp /)
  
  DO iter = 1, 10
     
     CALL cpu_time(start_time)
     
     !! For scattered interpolation uncomment the subroutine below
     CALL get_cylindrical_boundary(cylfunc, sphfunc, sphfunc_dr, &
          sphfunc_dth, 'quadratic')
     
     !CALL get_cylindrical_boundary(rho_ax, z_ax, dims, cylfunc, &
     !     6, sphfunc, sphfunc_dr, sphfunc_dth)
     
!!$     CALL get_cylindrical_surface(filename, &
!!$          cylfunc(1:maxrhopts+1,1:maxzpts+1), &
!!$          rho_ax(1:maxrhopts+1),z_ax(1:maxzpts+1), &
!!$          (/maxrhopts+1,maxzpts+1/), fdrule, 0.0_dp , &
!!$          efield, afield, lmax )
     
     CALL cpu_time(end_time)
     
     interp_time = end_time - start_time
     
     WRITE(*,*)
     WRITE(*,*)           '*********************************************'
     WRITE(*,'(A,F11.8)') ' Interpolation time (seconds): ', interp_time
     !WRITE(*,'(A,E20.12)') ' Interpolation time (seconds): ', interp_time
     WRITE(*,*)           '*********************************************'
     WRITE(*,*)
     
  ENDDO

  minerror = 1.0E10_dp
  maxerror = 0.0_dp
  
  DO itheta = 1, numthetapts
     DO ir = 1, numrpts
        
        ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  &
             0.25_dp * SQRT( 5.0_dp / pi) * &
             ( 3.0_dp * COS(theta_boundary(itheta))**2 - 1.0_dp)
        
        maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta) -&
             ref_value))
        minerror = MIN(minerror, ABS(sphfunc(ir,itheta) -&
             ref_value))
        
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
!!$  IF(rank.EQ.0) THEN
!!$     WRITE(*,*) 
!!$     WRITE(*,*) 'Write data to disk...'
!!$  ENDIF
!!$     
!!$  ! RHO
!!$  OPEN(unit=33,form='formatted',file='./results/rho_ax.dat')
!!$  DO irho = minrho, maxrho
!!$     WRITE(33,*) irho, rho_ax(irho) 
!!$  ENDDO
!!$  CLOSE(33)
!!$  
!!$  ! Z
!!$  OPEN(unit=33,form='formatted',file='./results/z_ax.dat')
!!$  DO iz = minz, maxz
!!$     WRITE(33,*) iz,z_ax(iz)
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
!!$  ! Save original wavefunction
!!$  OPEN(unit=33,form='formatted',file='./results/cartfunction.dat')
!!$  
!!$  DO iz = 1, numzpts
!!$     DO irho = 1, numrhopts
!!$        WRITE(33,*) ABS(cylfunc(irho,iz))**2 
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  CLOSE(33)
!!$  
!!$  ! Save spherical wavefunction
!!$  OPEN(unit=33,form='formatted',file='./results/sphfunction.dat')
!!$
!!$     DO itheta = 1, numthetapts
!!$        DO ir = 1, numrpts
!!$           WRITE(33,*) ABS(sphfunc(ir,itheta))**2 
!!$        ENDDO
!!$     ENDDO
!!$  
!!$  CLOSE(33)
  
  ! Free memory
  DEALLOCATE(rho_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cylfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth)

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '¡Se acabó!'
  WRITE(*,*)
  WRITE(*,*)
  
END PROGRAM cyl2sph_ex1
