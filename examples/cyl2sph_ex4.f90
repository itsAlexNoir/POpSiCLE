PROGRAM cyl2sph_ex4
  
  USE MPI
  USE popsicle
  
  IMPLICIT NONE
  
  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: oldmaxrhopts, oldmaxzpts
  INTEGER                    :: oldnumrhopts, oldnumzpts
  INTEGER                    :: newmaxrhopts, newmaxzpts
  INTEGER                    :: newnumrhopts, newnumzpts
  INTEGER                    :: dims_local(2), dims_global(2)
  
  INTEGER                    :: oldnumproc1drho, oldnumproc1dz
  INTEGER                    :: oldmaxproc1drho, oldmaxproc1dz
  INTEGER                    :: newnumproc1drho, newnumproc1dz
  INTEGER                    :: newmaxproc1drho, newmaxproc1dz
  INTEGER                    :: offsetprocrho, offsetprocz
  INTEGER                    :: ratioprocrho, ratioprocz
  
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
  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: comp_time
  REAL(dp)                   :: efield(3), afield(3)
  
  INTEGER                    :: size, rank, comm
  INTEGER                    :: maxprocessor, ierror
  
  CHARACTER(LEN = 150)       :: filename
  CHARACTER(LEN = 6)         :: ctime, rbstr
  CHARACTER(LEN = 4)         :: cprocessor, lmaxstr
  CHARACTER(LEN=100)         :: data_directory

  !-----------------------------------------------!
  
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
     WRITE(*,*) '       Cyl2sph4 example.'
     WRITE(*,*) ' We go parallel in this one!'
     WRITE(*,*) '****************************'
     WRITE(*,*) 
     WRITE(*,*) 'Opening...'
  ENDIF
  
  ! Set the parameters of the grid
  oldmaxrhopts     = 90
  oldmaxzpts	   = 201
  
  deltarho         = 0.1_dp
  deltaz           = 0.2_dp
  
  oldnumproc1drho  = 4
  oldnumproc1dz    = 10
  
  oldmaxproc1drho = oldnumproc1drho - 1
  oldmaxproc1dz   = oldnumproc1dz - 1 

  oldnumrhopts        = oldnumproc1drho * oldmaxrhopts
  oldnumzpts	      = oldnumproc1dz   * oldmaxzpts

  newnumproc1drho     = 1
  newnumproc1dz       = 5
  
  newmaxproc1drho = newnumproc1drho - 1
  newmaxproc1dz   = newnumproc1dz - 1 
  
  ratioprocrho = oldnumproc1drho / newnumproc1drho
  ratioprocz = oldnumproc1dz / newnumproc1dz
  
  offsetprocrho = 0 * ratioprocrho
  offsetprocz = rank * ratioprocz
  
  newmaxrhopts = oldmaxrhopts * ratioprocrho
  newmaxzpts = oldmaxzpts * ratioprocz
  
  newnumrhopts = newmaxrhopts * newnumproc1drho
  newnumzpts = newmaxzpts * newnumproc1dz
    
  dims_local = (/ newmaxrhopts, newmaxzpts /)
  dims_global = (/ newnumrhopts, newnumzpts /) 
  
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
  IF(rank.EQ.0) &
       WRITE(*,*) 'Building meshes ...'
  
  ALLOCATE(rho_ax(1:newmaxrhopts))
  ALLOCATE(z_ax(1:newmaxzpts))
  ALLOCATE(gpts(1:newmaxrhopts))
  ALLOCATE(gp(1:newmaxrhopts))
  ALLOCATE(hpts(1:newmaxzpts))
  ALLOCATE(hp(1:newmaxzpts))
  
  rhoalpha = 1.0_dp
  rhobeta = 0.0_dp
  
  DO irho = 1, newmaxrhopts
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
  
  DO iz = 1, newmaxzpts
     zpt = (-0.5_dp * REAL(newnumzpts-1,dp) + &
          REAL(newmaxzpts * rank + iz-1,dp)) * deltaz
     z_ax(iz) = zpt
     hpts(iz) = zpt
     hp(iz) = 1.0_dp
  ENDDO
  
  !----------------------------------------!
  ! Initialize cylindrical surface stuff
  !----------------------------------------!
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  CALL initialize_cylindrical_surface(gpts, hpts, dims_local, &
       Rboundary, tolerance, fdrule, deltar, lmax, filename, &
       rank, size, MPI_COMM_WORLD)
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
  WRITE(*,'(A50,I2.2,A5,F9.6)') 'Time spent on initializing the surface &
       on proc ',rank,' (s): ',comp_time
  
  !-----------------------------!
  ! Read wavefunction from disk
  !-----------------------------!
  
  IF(rank .EQ. 0) &
       WRITE(*,*) 'Allocating wavefunctions...'
  ALLOCATE(psipro(1:oldmaxrhopts, 1:oldmaxzpts))
  ALLOCATE(psi(1:newmaxrhopts, 1:newmaxzpts))
  
  IF(rank.EQ.0) &
       WRITE(*,*) 'Read data...'
  
  DO iprocz = 0, ratioprocz-1     
     DO iprocrho = 0, ratioprocrho-1    
        
        ipro = ((offsetprocz + iprocz) * oldnumproc1drho) + &
             offsetprocrho + iprocrho
        
        WRITE(cprocessor, '(I4.4)') ipro
        
        filename = TRIM(data_directory) // 'psi/' // cprocessor //    &
             '/psi.' // cprocessor // '.dat'
        
        OPEN(UNIT = 10, FORM = 'unformatted', FILE = filename)
        
        READ(10) psipro
        
        CLOSE(UNIT = 10)
        
        DO iz = 1, oldmaxzpts
           DO irho = 1, oldmaxrhopts
              
              irhog = iprocrho * oldmaxrhopts + irho
              izg = iprocz * oldmaxzpts + iz
              
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
  
  DO iz = 1, newmaxzpts
     DO irho = 1, newmaxrhopts
        psi(irho,iz) = psi(irho,iz) &
             / SQRT(gpts(irho) * gp(irho) * hp(iz) )
     ENDDO
  ENDDO
  
  !-------------------------------------------------!
  ! Get cylindrical surface, and write it to a file
  !-------------------------------------------------!
  IF(i_am_surface(rank).EQ.1) &
       WRITE(*,*) 'Get the surface. Write it to a file...'

  efield = 0.0_dp
  afield = 0.0_dp
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  
  CALL get_cylindrical_surface(filename, psi, gpts, hpts, dims_local, &
       fdrule, 0.0_dp , &
       efield, afield, lmax, rank, size )
  
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
   IF(i_am_surface(rank).EQ.1) &  
  WRITE(*,'(A50,I2.2,A5,F8.5)') 'Time spent on getting the surface&
       & on proc ',rank,' (s): ',comp_time
  
  !------------------------------!
  ! Now, load it again from file
  !------------------------------!
   IF(i_am_surface(rank).EQ.1) &
       WRITE(*,*) 'Deleting surface...'
  CALL delete_surface2D( rank )
  
  ! Free arrays!
  IF(rank.EQ.0) &
       WRITE(*,*) 'Free arrays!'
  DEALLOCATE(psi,psipro)
  DEALLOCATE(rho_ax,z_ax,gpts,hpts,gp,hp)
  
  CALL MPI_finalize( ierror )
  
  IF(i_am_surface(rank).EQ.1) THEN
     WRITE(*,*) '¡Se acabó!'
     WRITE(*,*)
     WRITE(*,*)
  ENDIF
  
END PROGRAM cyl2sph_ex4
