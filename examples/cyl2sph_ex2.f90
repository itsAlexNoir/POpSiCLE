PROGRAM cyl2sph_ex2
  
  USE MPI
  USE popsicle
  
  IMPLICIT NONE
  
  INTEGER                    :: maxrhopts, maxzpts
  INTEGER                    :: numrhopts, numzpts
  INTEGER                    :: maxrho, minrho
  INTEGER                    :: maxz, minz
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp)                   :: rhoalpha, rhobeta
  REAL(dp)                   :: sqrtrho, sqrt1rho
  REAL(dp)                   :: rhopt, zpt
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  INTEGER                    :: local_dims(2),global_dims(2)
  
  COMPLEX(dp), ALLOCATABLE   :: Y_lm(:, :)
  COMPLEX(dp), ALLOCATABLE   :: R_nl(:, :)
  COMPLEX(dp), ALLOCATABLE   :: cylfunc(:, :)
  REAL(dp)                   :: drho, dz
  COMPLEX(dp), ALLOCATABLE   :: sphfunc(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dr(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dth(:, :)
  COMPLEX(dp)                :: ref_value
  REAL(dp)                   :: dr
  INTEGER                    :: lmax, fdrule
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: rpt, thetapt
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: minerror, maxerror
  INTEGER                    :: irho,iz
  INTEGER                    :: ir, itheta
  REAL(dp)                   :: efield(3), afield(3)
  
  INTEGER                    :: size, rank, comm
  INTEGER                    :: maxprocessor, ierror
  INTEGER                    :: newcomm
  CHARACTER(LEN=100)         :: filename
  
  !------------------------------------------------------------
  
  ! Start up MPI
  CALL MPI_init( ierror )
  
  ! Find out number of processors.
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierror )
  
  maxprocessor = size - 1
  
  ! Find out number of the processor we are working on.
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierror )
  
  IF(rank.EQ.0) THEN
     WRITE(*,*)
     WRITE(*,*) '****************************'
     WRITE(*,*) '       Cyl2sph2 example.'
     WRITE(*,*) ' We go parallel in this one!'
     WRITE(*,*) '****************************'
     WRITE(*,*) 
     WRITE(*,*) 'Opening...'
  ENDIF
  
  ! Set number of points
  numrhopts = 300
  numzpts = 601
  
  maxrhopts = numrhopts
  maxzpts = numzpts / size
  
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
!!$     rhopt = REAL(irho,dp) * drho
!!$     sqrtrho	    = SQRT(rhopt)
!!$     sqrt1rho	    = SQRT(rhoalpha + rhobeta * rhopt)   
!!$     rho_ax(irho)   = rhopt * sqrtrho / sqrt1rho
     rho_ax(irho)  = REAL(irho,dp) * drho
  ENDDO
  
  DO iz = minz, maxz
     z_ax(iz) = REAL(-numzpts / 2 + rank * maxzpts + iz,dp) * dz
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
        
!!$           Y_lm(ix,iy,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
!!$                ( 3.0_dp * COS(thetapt)**2 - 1.0_dp) !* COS(phipt) !EXP(ZIMAGONE*phipt)
        
        Y_lm(irho,iz) = 0.5_dp * SQRT(3.0_dp / pi ) * COS(thetapt)
        R_nl(irho,iz) = 1.0_dp
        
!!$        cylfunc(irho,iz) = EXP(-rpt / 2.0_dp) * R_nl(irho, iz) * &
!!$             Y_lm(irho,iz)
        
        cylfunc(irho,iz) = EXP( ZIMAGONE * 2.0_dp * rpt) * Y_lm(irho,iz)
        
        cylfunc_deriv(irho,iz) = ZIMAGONE * 2.0_dp * &
             EXP( ZIMAGONE * 2.0_dp * rpt) * Y_lm(irho,iz)
                
     ENDDO
  ENDDO
  
  
  
  ! Initialize the boundary
  Rboundary = 25.0_dp
  tolerance = 0.3_dp !0.15_dp
  
  dr = 0.1_dp
  lmax = 10
  fdrule = 2
  external_derivative = .TRUE.
  local_dims  = (/maxrhopts, maxzpts/)
  global_dims = (/numrhopts, numzpts/)
  filename    = 'results/cyl2sph'
  
  IF(rank.EQ.0) THEN
     WRITE(*,*) 'Number of proccessors: ',size
     
     WRITE(*,*) 'Number of points in rho per proc: ',maxrhopts
     WRITE(*,*) 'Number of points in z per proc: ',maxzpts
     
     WRITE(*,*) 'Number of points in rho: ',numrhopts
     WRITE(*,*) 'Number of points in z: ',numzpts
     
     WRITE(*,*) 'Grid spacing in rho: ',drho
     WRITE(*,*) 'Grid spacing in z: ',dz
     
     WRITE(*,*) 'Boundary at radius: ',Rboundary
     WRITE(*,*) 'Radius tolerance: ',tolerance
     WRITE(*,*)
     WRITE(*,*) '--------------------------'
  ENDIF
  
  ! Build interpolant
  IF(rank.EQ.0) &
       WRITE(*,*) 'Creating interpolant...'
  
  CALL cpu_time(start_time)
  
  !! For scattered interpolation uncomment the subroutine below
  !CALL initialize_cylindrical_boundary(rho_ax, z_ax, local_dims, &
  !     Rboundary, tolerance, 2, dr, lmax, rank, size, MPI_COMM_WORLD )
  
!!$  CALL initialize_cylindrical_boundary(rho_ax, z_ax, local_dims, &
!!$       Rboundary, 2, dr, lmax, rank, size, MPI_COMM_WORLD )
  
  CALL initialize_cylindrical_surface(rho_ax(1:maxrhopts+1), &
       z_ax(1:maxzpts+1),(/maxrhopts+1,maxzpts+1/), Rboundary, &
       tolerance, fdrule, dr, lmax, filename, &
       mpi_rank=rank, mpi_size=size, comm=MPI_COMM_WORLD)
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'From rank ',rank,' Interpolant time (seconds): ', interp_time
  WRITE(*,*) 
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts))
  
  ! Interpolate!!
  IF(rank.EQ.0) &
       WRITE(*,*) 'Interpolating boundary...'
  
  efield = (/0.0_dp, 0.0_dp, 0.0_dp /)
  afield = (/0.0_dp, 0.0_dp, 0.0_dp /)
  
  IF(i_am_surface_local(rank)) THEN
     CALL cpu_time(start_time)
     
     !! For scattered interpolation uncomment the subroutine below
     !CALL get_cylindrical_boundary(cylfunc, sphfunc, sphfunc_dr, &
     !     sphfunc_dth, 'quadratic')
     
!!$     CALL get_cylindrical_boundary(rho_ax, z_ax, local_dims, cylfunc, &
!!$          2, sphfunc, sphfunc_dr, sphfunc_dth)
     
     CALL get_cylindrical_surface(filename, &
          cylfunc(1:maxrhopts+1,1:maxzpts+1), &
          rho_ax(1:maxrhopts+1),z_ax(1:maxzpts+1), &
          (/maxrhopts+1,maxzpts+1/), fdrule, 0.0_dp , &
          efield, afield, lmax, mpi_rank=rank, mpi_size=size )
     
     CALL cpu_time(end_time)
     
     interp_time = end_time - start_time
     
     WRITE(*,*) 'From rank ',rank,' Interpolation time (seconds): ', interp_time
     
     DO itheta = 1, numthetapts
        DO ir = 1, numrpts
           IF(i_am_in_local2D(ir,itheta)) THEN
              ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  &
                   !0.25_dp * SQRT(5.0_dp / pi) * &
                   !( 3.0_dp * COS(theta_boundary(itheta))**2 - 1.0_dp)
                   0.5_dp * SQRT(3.0_dp / pi) * costheta_boundary(itheta)
              
              maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta) -&
                   ref_value))
              minerror = MIN(minerror, ABS(sphfunc(ir,itheta) -&
                   ref_value))
              
           ENDIF
        ENDDO
     ENDDO
     
     
     WRITE(*,*) 
     WRITE(*,*) '********RESULTS!!!******'
     WRITE(*,*)
     WRITE(*,*) 'Minimum difference in value: '
     WRITE(*,*) minerror
     WRITE(*,*) 'Maximum difference in value: '
     WRITE(*,*) maxerror
     
  ENDIF
  
  
!!$  ! Save axes
!!$  IF(rank.EQ.0) THEN
!!$     WRITE(*,*) 
!!$     WRITE(*,*) 'Write data to disk...'
!!$  ENDIF
!!$     
!!$  ! RHO
!!$  OPEN(unit=33,form='formatted',file='./results/rho_ax.dat')
!!$  DO irho = 1, numrhopts
!!$     WRITE(33,*) rho_ax(irho) 
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
!!$
  
  ! Free memory
  DEALLOCATE(rho_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cylfunc)
  DEALLOCATE(sphfunc)
  DEALLOCATE(sphfunc_dr, sphfunc_dth)
  
  CALL MPI_finalize( ierror )
  
  IF(rank.EQ.0) THEN
     WRITE(*,*)
     WRITE(*,*)
     WRITE(*,*) '¡Se acabó!'
     WRITE(*,*)
     WRITE(*,*)
  ENDIF
  
END PROGRAM cyl2sph_ex2
