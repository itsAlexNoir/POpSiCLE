PROGRAM cy2sph_ex3
  
  USE popsicle_aux
  
  IMPLICIT NONE

  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: maxrhopts, maxzpts
  INTEGER                    :: numrhopts, numzpts
  INTEGER                    :: dims(2)
  INTEGER, ALLOCATABLE       :: halopts(:,:)

  INTEGER                    :: numproc1drho       
  INTEGER                    :: numproc1dz                  
  INTEGER                    :: maxproc1drho       
  INTEGER                    :: maxproc1dz         
  
  REAL(dp)                   :: deltarho, deltaz
  REAL(dp)                   :: Rboundary
  REAL(dp)                   :: deltar
  INTEGER                    :: lmax
  REAL(dp)                   :: tolerance
  INTEGER                    :: fdrule
  
  INTEGER                    :: irho, iz
  INTEGER                    :: irg, irhog, izg
  INTEGER                    :: iprocrho, iprocz
  INTEGER                    :: ipro
  
  COMPLEX(dp), ALLOCATABLE   :: psipro(:, :)
  COMPLEX(dp), ALLOCATABLE   :: psi(:, :)
  
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  REAL(dp), ALLOCATABLE      :: gpts(:), gp(:)
  REAL(dp), ALLOCATABLE      :: hpts(:), hp(:)
  REAL(dp)                   :: rhoalpha, rhobeta
  REAL(dp)                   :: sqrtrho, sqrt1rho
  REAL(dp)                   :: rhopt, zpt

  REAL(dp)                   :: efield(3), afield(3)
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: comp_time
  
  CHARACTER(LEN = 150)       :: filename
  CHARACTER(LEN = 6)         :: ctime, rbstr
  CHARACTER(LEN = 4)         :: cprocessor, lmaxstr
  CHARACTER(LEN=100)         :: data_directory

  !-----------------------------------------------!
  
  WRITE(*,*)
  WRITE(*,*) '*******************'
  WRITE(*,*) '  Cy2sph3 example.'
  WRITE(*,*) '*******************'
  WRITE(*,*) 
  
  WRITE(*,*) 'Opening...'
  
  ! Set the parameters of the grid
  maxrhopts        = 90
  maxzpts	       = 201

  deltarho         = 0.1_dp
  deltaz           = 0.2_dp
  
  numproc1drho     = 4
  numproc1dz       = 10
  
  maxproc1drho = numproc1drho - 1
  maxproc1dz   = numproc1dz - 1 
  
  numrhopts        = numproc1drho * maxrhopts
  numzpts	   = numproc1dz   * maxzpts
  dims = (/ numrhopts, numzpts /)
  
  ! The radius of the boundary
  Rboundary    = 50.0_dp
  tolerance = 1.0_dp !0.75_dp
  deltar = 0.1_dp
  fdrule = 2
  lmax = 10
  
  data_directory = './data/h2p/cylindrical/'
  WRITE(rbstr,'(F6.3)') Rboundary
  WRITE(lmaxstr,'(I3.3)') lmax

  !--------------!
  ! Create grids
  !--------------!
  WRITE(*,*) 'Building meshes ...'

  ALLOCATE(rho_ax(1:numrhopts))
  ALLOCATE(z_ax(1:numzpts))
  ALLOCATE(gpts(1:numrhopts))
  ALLOCATE(gp(1:numrhopts))
  ALLOCATE(hpts(1:numzpts))
  ALLOCATE(hp(1:numzpts))
  
  rhoalpha = 1.0_dp
  rhobeta = 0.0_dp
  
  DO irho = 1, numrhopts
     rhopt = REAL(irho,dp) * deltarho
     rho_ax(irho) = rhopt
     sqrtrho	    = SQRT(rhopt)
     sqrt1rho	    = SQRT(rhoalpha + rhobeta * rhopt)
     
     gpts(irho)	    = rhopt * sqrtrho / sqrt1rho
     
     gp(irho)	    = 0.5_dp * sqrtrho *			       &
          (2.0_dp * rhobeta * rhopt +                 &
          3.0_dp * rhoalpha)
     gp(irho)	    = gp(irho) / (sqrt1rho ** 3)
  ENDDO
  
  DO iz = 1, numzpts
     zpt = (-0.5_dp * (numzpts-1) + REAL(iz-1,dp)) * deltaz
     z_ax(iz) = zpt
     hpts(iz) = zpt
     hp(iz) = 1.0_dp
  ENDDO

!!$  ALLOCATE(halopts(2,2))
!!$  halopts(1,1) = -1
!!$  halopts(2,1) = -1
!!$  halopts(1,2) = numrhopts
!!$  halopts(2,2) = numzpts
  
  !----------------------------------------!
  ! Initialize cylindrical surface stuff
  !----------------------------------------!
  
!!$  CALL initializate_cylindrical_grid(rho_ax,z_ax,(/numrhopts, numzpts/),halopts,&
!!$       Rboundary, 1.0_dp, 2, deltar, lmax)
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  CALL initialize_cylindrical_surface(gpts, hpts, dims, &
       Rboundary, tolerance, fdrule, deltar, lmax, filename, 0)
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
  WRITE(*,'(A43,F9.6)') 'Time spent on initializing the surface(s): ',comp_time
  
  !-----------------------------!
  ! Read wavefunction from disk
  !-----------------------------!
  
  WRITE(*,*) 'Allocating wavefunctions...'
  ALLOCATE(psipro(1:maxrhopts, 1:maxzpts))
  ALLOCATE(psi(1:numrhopts, 1:numzpts))

  WRITE(*,*) 'Read data...'
  
  ipro = -1
  
  DO iprocz = 0, maxproc1dz     
     DO iprocrho = 0, maxproc1drho    
        
        ipro = ipro + 1
        
        WRITE(cprocessor, '(I4.4)') ipro
        
        filename = TRIM(data_directory) // 'psi/' // cprocessor //    &
             '/psi.' // cprocessor // '.dat'
        
        OPEN(UNIT = 10, FORM = 'unformatted', FILE = filename)
        
        READ(10) psipro
        
        CLOSE(UNIT = 10)
        
        DO iz = 1, maxzpts
           DO irho = 1, maxrhopts
              
              irhog = iprocrho * maxrhopts + irho
              izg = iprocz * maxzpts + iz
              
              psi(irhog, izg) = psi(irhog, izg) + psipro(irho, iz)
              
           ENDDO
        ENDDO
        
     ENDDO
  ENDDO

  !----------------------------------------------!
  ! Unscale the original wavefunction
  ! (This has to do with how the original TDSE
  ! where this wavefunction came from works)
  !----------------------------------------------!
  
  DO iz = 1, numzpts
     DO irho = 1, numrhopts
        psi(irho,iz) = psi(irho,iz) &
             / SQRT(gpts(irho) * gp(irho) * hp(iz) )
     ENDDO
  ENDDO
  
  !-------------------------------------------------!
  ! Get cylindrical surface, and write it to a file
  !-------------------------------------------------!
  efield = 0.0_dp
  afield = 0.0_dp
  
  WRITE(*,*) 'Get the surface. Write it to a file...'

  CALL cpu_time(start_time)
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL get_cylindrical_surface(filename, 0, psi, gpts, hpts, dims, &
       fdrule, 0.0_dp , efield, afield, lmax)
  
  CALL cpu_time(end_time)

  comp_time = end_time - start_time

  WRITE(*,'(A38,F8.5)') 'Time spent on getting the surface (s): ',comp_time
  
  !------------------------------!
  ! Now, load it again from file
  !------------------------------!
  WRITE(*,*) 'Deleting surface...'
  CALL delete_surface2D( )
  
  ! Free arrays!
  WRITE(*,*) 'Free arrays!'
  DEALLOCATE(psi,psipro)
  DEALLOCATE(rho_ax,z_ax,gpts,hpts,gp,hp)
  
  WRITE(*,*) '¡Se acabó!'
  WRITE(*,*)
  WRITE(*,*)
  
END PROGRAM cy2sph_ex3
