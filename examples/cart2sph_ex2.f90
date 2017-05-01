PROGRAM cart2sph_ex2

  USE MPI
  USE popsicle_aux
  
  IMPLICIT NONE

  INTEGER                    :: maxxpts
  INTEGER                    :: maxypts  
  INTEGER                    :: maxzpts
  INTEGER                    :: numxpts
  INTEGER                    :: numypts  
  INTEGER                    :: numzpts
  INTEGER                    :: maxx, minx
  INTEGER                    :: maxy, miny
  INTEGER                    :: maxz, minz
  REAL(dp), ALLOCATABLE      :: x_ax(:)
  REAL(dp), ALLOCATABLE      :: y_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  INTEGER                    :: local_dims(3),global_dims(3)
  
  COMPLEX(dp), ALLOCATABLE   :: Y_lm(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: R_nl(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: cartfunc(:, :, :)
  REAL(dp)                   :: dx, dy, dz
  COMPLEX(dp), ALLOCATABLE   :: sphfunc(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dth(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dphi(:, :, :)
  COMPLEX(dp)                :: ref_value
  REAL(dp)                   :: dr
  INTEGER                    :: lmax
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: rpt, thetapt, phipt
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: minerror, maxerror
  INTEGER                    :: ix,iy, iz
  INTEGER                    :: ir, itheta, iphi
  
  INTEGER                    :: size, rank, comm
  INTEGER                    :: maxprocessor, ierror
  INTEGER                    :: newcomm

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
     WRITE(*,*) '       Cart2sph2 example.'
     WRITE(*,*) ' We go parallel in this one!'
     WRITE(*,*) '****************************'
     WRITE(*,*) 
     WRITE(*,*) 'Opening...'
  ENDIF
  
  ! Set number of points
  numxpts = 80
  numypts = 80
  numzpts = 80
  
  maxxpts = numxpts / size
  maxypts = numypts
  maxzpts = numzpts

  minx = 1
  maxx = maxxpts + 1
  
  miny = 1
  maxy = maxypts + 1

  minz = 1
  maxz = maxzpts + 1

  ! Set grid spacing
  dx   = 0.1_dp
  dy   = 0.1_dp
  dz   = 0.1_dp
  
  ! Create axes
  ALLOCATE(x_ax(minx:maxx))
  ALLOCATE(y_ax(miny:maxy))
  ALLOCATE(z_ax(minz:maxz))
  
  DO ix = minx, maxx
     x_ax(ix) = REAL(-numxpts / 2 + rank * maxxpts + ix,dp) * dx
  ENDDO

  DO iy = miny, maxy
     y_ax(iy) = REAL(-numypts / 2 + iy,dp) * dy
  ENDDO
  
  DO iz = minz, maxz
     z_ax(iz) = REAL(-numzpts / 2 + iz,dp) * dz  
  ENDDO
  
!!$  write(*,*) rank,x_ax
  
  ! We are going to set the 3d hydrogen function
  ! n=3, l=2, m=0
  ALLOCATE(Y_lm(minx:maxx,miny:maxy,minz:maxz))
  ALLOCATE(R_nl(minx:maxx,miny:maxy,minz:maxz))
  ALLOCATE(cartfunc(minx:maxx,miny:maxy,minz:maxz))
  
  DO iz = minz, maxz
     DO iy = miny, maxy  
        DO ix = minx, maxx
           
           rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
           
           IF (rpt.EQ.0) THEN
              thetapt = pi / 2.0_dp
           ELSE
              thetapt = ACOS( z_ax(iz) / rpt )
           ENDIF
           
           phipt = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
           
!!$           Y_lm(ix,iy,iz) = 0.25_dp * SQRT(5.0_dp / pi) * &
!!$                ( 3.0_dp * COS(thetapt)**2 - 1.0_dp) !* COS(phipt) !EXP(ZIMAGONE*phipt)
           
           Y_lm(ix,iy,iz) = 0.5_dp * SQRT(3.0_dp / pi ) * COS(thetapt)
           R_nl(ix,iy,iz) = 1.0_dp
           
           cartfunc(ix,iy,iz) = EXP(-rpt / 2.0_dp) * R_nl(ix,iy,iz) * &
                Y_lm(ix,iy,iz)
           
        ENDDO
     ENDDO
  ENDDO
  
  
  ! Initialize the boundary
  Rboundary = 1.0_dp
  tolerance = 0.15_dp
  dr = 0.1_dp
  lmax = 8
  local_dims      = (/maxxpts, maxypts, maxzpts/)
  global_dims      = (/numxpts, numypts, numzpts/)
  
  IF(rank.EQ.0) THEN
     WRITE(*,*) 'Number of proccessors: ',size
     
     WRITE(*,*) 'Number of points in x per proc: ',maxxpts
     WRITE(*,*) 'Number of points in y per proc: ',maxypts
     WRITE(*,*) 'Number of points in z per proc: ',maxzpts
     
     WRITE(*,*) 'Number of points in x: ',numxpts
     WRITE(*,*) 'Number of points in y: ',numypts
     WRITE(*,*) 'Number of points in z: ',numzpts
     
     WRITE(*,*) 'Grid spacing in x: ',dx
     WRITE(*,*) 'Grid spacing in y: ',dy
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
  
  CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, local_dims, &
       Rboundary, 2, dr, lmax, rank, size, MPI_COMM_WORLD )
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,*) 'From rank ',rank,' Interpolation time (seconds): ', interp_time
  WRITE(*,*) 'Interpolant time (seconds): ', interp_time
  WRITE(*,*) 
  
  !WRITE(*,*) 'Total number of points to be interpolated: ',numpts
  
  WRITE(*,*) 'Number of radial boundary points: ',numrpts
  WRITE(*,*) 'Number of polar boundary points: ',numthetapts
  WRITE(*,*) 'Number of azimuthal boundary points: ',numphipts
  WRITE(*,*) 'Number of proccessors involved: ',numsurfaceprocs
  
  WRITE(*,*) 'Grid spacing in r: ',dr
  WRITE(*,*) 'Maximum angular momenta: ',lmax
  !WRITE(*,*) 'Grid spacing in theta: ',dtheta
  !WRITE(*,*) 'Grid spacing in phi: ',dphi
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '--------------------------'
  
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(sphfunc_dphi(1:numrpts,1:numthetapts,1:numphipts))
  
  ! Interpolate!!
  IF(rank.EQ.0) &
       WRITE(*,*) 'Interpolating boundary...'
  
  IF(i_am_surface_local(rank)) THEN
     CALL cpu_time(start_time)
     
     CALL get_cartesian_boundary(x_ax, y_ax, z_ax, local_dims, &
          cartfunc, 2, sphfunc, sphfunc_dr, &
          sphfunc_dth, sphfunc_dphi)
     
     CALL cpu_time(end_time)
     
     interp_time = end_time - start_time
     
     WRITE(*,*) 'From rank ',rank,' Interpolation time (seconds): ', interp_time
     
     DO iphi = 1, numphipts
        DO itheta = 1, numthetapts
           DO ir = 1, numrpts
              IF(i_am_in_local3D(ir,itheta,iphi)) THEN
                 ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  &
                      !0.25_dp * SQRT(5.0_dp / pi) * &
                      !( 3.0_dp * COS(theta_boundary(itheta))**2 - 1.0_dp)
                      0.5_dp * SQRT(3.0_dp / pi) * costheta_boundary(itheta)
                 
                 maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta, iphi) -&
                      ref_value))
                 minerror = MIN(minerror, ABS(sphfunc(ir,itheta, iphi) -&
                      ref_value))
                 
              ENDIF
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
     
  ENDIF
  
  
!!$  ! Save axes
!!$  IF(rank.EQ.0) THEN
!!$     WRITE(*,*) 
!!$     WRITE(*,*) 'Write data to disk...'
!!$  ENDIF
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
!!$
  
  ! Free memory
  DEALLOCATE(x_ax,y_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cartfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth,sphfunc_dphi)
  
  CALL MPI_finalize( ierror )
  
  IF(rank.EQ.0) THEN
     WRITE(*,*) '¡Se acabó!'
     WRITE(*,*)
     WRITE(*,*)
  ENDIF
  
END PROGRAM cart2sph_ex2
