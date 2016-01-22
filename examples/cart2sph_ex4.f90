PROGRAM cart2sph_ex4
  
  USE MPI
  USE popsicle
  
  IMPLICIT NONE
  
  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: oldmaxxpts, oldmaxypts, oldmaxzpts
  INTEGER                    :: oldnumxpts, oldnumypts, oldnumzpts
  INTEGER                    :: newmaxxpts, newmaxypts, newmaxzpts
  INTEGER                    :: newnumxpts, newnumypts, newnumzpts
  INTEGER                    :: minx, miny, minz
  INTEGER                    :: maxx, maxy, maxz
  INTEGER                    :: dims_local(3), dims_global(3)
  
  INTEGER                    :: oldnumproc1dx, oldnumproc1dy, oldnumproc1dz
  INTEGER                    :: oldmaxproc1dx, oldmaxproc1dy, oldmaxproc1dz
  INTEGER                    :: newnumproc1dx, newnumproc1dy, newnumproc1dz
  INTEGER                    :: newmaxproc1dx, newmaxproc1dy, newmaxproc1dz
  INTEGER                    :: offsetprocx, offsetprocy, offsetprocz
  INTEGER                    :: ratioprocx, ratioprocy, ratioprocz
  
  REAL(dp)                   :: deltax, deltay, deltaz
  REAL(dp)                   :: Rboundary
  REAL(dp)                   :: deltar
  INTEGER                    :: lmax
  REAL(dp)                   :: tolerance
  INTEGER                    :: fdrule
  
  INTEGER                    :: ix, iy, iz
  INTEGER                    :: ixg, iyg, izg
  INTEGER                    :: iprocx, iprocy, iprocz
  INTEGER                    :: ipro
  
  COMPLEX(dp), ALLOCATABLE   :: psipro(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: psi(:, :, :)
  
  REAL(dp), ALLOCATABLE      :: x_ax(:)
  REAL(dp), ALLOCATABLE      :: y_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  REAL(dp), ALLOCATABLE      :: xpts(:), xp(:)
  REAL(dp), ALLOCATABLE      :: ypts(:), yp(:)
  REAL(dp), ALLOCATABLE      :: zpts(:), zp(:)
  REAL(dp)                   :: xpt, ypt, zpt
  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: comp_time
  
  INTEGER                    :: size, rank, comm
  INTEGER                    :: maxprocessor, ierror
  
  CHARACTER(LEN = 150)       :: filename
  CHARACTER(LEN = 6)         :: ctime
  CHARACTER(LEN = 7)         :: rbstr
  CHARACTER(LEN = 5)         :: cprocessor
  CHARACTER(LEN = 4)         :: lmaxstr
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
     WRITE(*,*) '       Cart2sph3 example.'
     WRITE(*,*) ' We go parallel in this one!'
     WRITE(*,*) '****************************'
     WRITE(*,*) 
     WRITE(*,*) 'Opening...'
  ENDIF
  
  ! Set the parameters of the grid
  oldmaxxpts       = 30
  oldmaxypts       = 30
  oldmaxzpts	   = 60
  
  deltax           = 0.2_dp
  deltay           = 0.2_dp
  deltaz           = 0.2_dp
  
  oldnumproc1dx    = 4
  oldnumproc1dy    = 4
  oldnumproc1dz    = 2
  
  oldmaxproc1dx   = oldnumproc1dx - 1
  oldmaxproc1dy   = oldnumproc1dy - 1
  oldmaxproc1dz   = oldnumproc1dz - 1 

  oldnumxpts          = oldnumproc1dx * oldmaxxpts
  oldnumypts          = oldnumproc1dy * oldmaxypts 
  oldnumzpts	      = oldnumproc1dz * oldmaxzpts

  newnumproc1dx       = 4
  newnumproc1dy       = 1
  newnumproc1dz       = 1
  
  newmaxproc1dx = newnumproc1dx - 1
  newmaxproc1dy = newnumproc1dy - 1
  newmaxproc1dz = newnumproc1dz - 1 
  
  ratioprocx = oldnumproc1dx / newnumproc1dx
  ratioprocy = oldnumproc1dy / newnumproc1dy
  ratioprocz = oldnumproc1dz / newnumproc1dz
  
  offsetprocx = rank * ratioprocx
  offsetprocy = 0 * ratioprocy  
  offsetprocz = 0 * ratioprocz
  
  newmaxxpts = oldmaxxpts * ratioprocx
  newmaxypts = oldmaxypts * ratioprocy 
  newmaxzpts = oldmaxzpts * ratioprocz

  minx = 1
  maxx = newmaxxpts + 1

  miny = 1
  maxy = newmaxypts + 1
  
  minz = 1
  maxz = newmaxzpts + 1
  
  newnumxpts = newmaxxpts * newnumproc1dx
  newnumypts = newmaxypts * newnumproc1dy
  newnumzpts = newmaxzpts * newnumproc1dz
  
  dims_local = (/ newmaxxpts, newmaxypts, newmaxzpts /)
  dims_global = (/ newnumxpts, newnumypts, newnumzpts /)
  
  ! The radius of the boundary
  Rboundary    = 8.0_dp
  tolerance = 0.25_dp
  deltar = 0.1_dp
  fdrule = 2
  lmax = 10
  
  data_directory = './data/h2p/cartesian/'
  WRITE(rbstr,'(I3.3,F0.3)') INT(Rboundary),Rboundary-INT(Rboundary)
  WRITE(lmaxstr,'(I3.3)') lmax
  
  !--------------!
  ! Create grids
  !--------------!
  IF(rank.EQ.0) &
       WRITE(*,*) 'Building meshes ...'
  
  ALLOCATE(x_ax(minx:maxx))
  ALLOCATE(y_ax(miny:maxy))
  ALLOCATE(z_ax(minz:maxz))
  ALLOCATE(xpts(minx:maxx))
  ALLOCATE(xp(minx:maxx))
  ALLOCATE(ypts(miny:maxy))
  ALLOCATE(yp(miny:maxy))
  ALLOCATE(zpts(minz:maxz))
  ALLOCATE(zp(minz:maxz))
  
  DO ix = minx, maxx
     xpt = (-0.5_dp * REAL(newnumxpts-1,dp) + &
          REAL(newmaxxpts * rank + ix-1,dp)) * deltax
     x_ax(ix) = xpt
     xpts(ix) = xpt
     xp(ix) = 1.0_dp
  ENDDO
  
  DO iy = miny, maxy
     ypt = (-0.5_dp * REAL(newnumypts-1,dp) + &
          REAL(iy-1,dp)) * deltay
     y_ax(iy) = ypt
     ypts(iy) = ypt
     yp(iy) = 1.0_dp
  ENDDO
  
  DO iz = minz, maxz
     zpt = (-0.5_dp * REAL(newnumzpts-1,dp) + &
          REAL(iz-1,dp)) * deltaz
     z_ax(iz) = zpt
     zpts(iz) = zpt
     zp(iz) = 1.0_dp
  ENDDO
    
  !----------------------------------------!
  ! Initialize cylindrical surface stuff
  !----------------------------------------!
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  CALL initialize_cartesian_surface(xpts, ypts, zpts, dims_local, &
       Rboundary, tolerance, fdrule, deltar, lmax, &
       rank, size,MPI_COMM_WORLD, .TRUE., filename)
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
  WRITE(*,'(A50,I2.2,A5,F9.6)') 'Time spent on initializing the surface &
       on proc ',rank,' (s): ',comp_time
  
  !-----------------------------!
  ! Read wavefunction from disk
  !-----------------------------!
  
  IF(rank .EQ. 0) &
       WRITE(*,*) 'Allocating wavefunctions...'
  
  ALLOCATE(psipro(1:oldmaxxpts, 1:oldmaxypts, 1:oldmaxzpts))
  !ALLOCATE(psi(1:newmaxxpts, 1:newmaxypts, 1:newmaxzpts))
  ALLOCATE(psi(minx:maxx, miny:maxy, minz:maxz))
  
  IF(rank.EQ.0) &
       WRITE(*,*) 'Read data...'
  
  DO iprocz = 0, ratioprocz-1     
     DO iprocy = 0, ratioprocy-1    
        DO iprocx = 0, ratioprocx-1    
           
           ipro = ((offsetprocz + iprocz) * oldnumproc1dx * oldnumproc1dy ) + &
                ((offsetprocy + iprocy) * oldnumproc1dx) + &
                offsetprocx + iprocx
           
           WRITE(cprocessor, '(I5.5)') ipro
           
           filename = TRIM(data_directory) // 'psi/' // cprocessor //    &
                '/psi.' // cprocessor // '.dat'
           
           OPEN(UNIT = 10, FORM = 'unformatted', FILE = filename)
           
           READ(10) psipro
           
           CLOSE(UNIT = 10)
           
           DO iz = 1, oldmaxzpts
              DO iy = 1, oldmaxypts
                 DO ix = 1, oldmaxxpts
                    
                    ixg = iprocx * oldmaxxpts + ix
                    iyg = iprocy * oldmaxypts + iy 
                    izg = iprocz * oldmaxzpts + iz
                    
                    psi(ixg, iyg, izg) = psi(ixg, iyg, izg) + psipro(ix, iy, iz)
                    
                 ENDDO
              ENDDO
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
     DO iy = 1, newmaxypts
        DO ix = 1, newmaxxpts
           
           psi(ix,iy,iz) = psi(ix,iy,iz) &
                / SQRT(xp(ix) * yp(iy) * zp(iz) )
        ENDDO
     ENDDO
  ENDDO
  
  !-------------------------------------------------!
  ! Get cylindrical surface, and write it to a file
  !-------------------------------------------------!
  IF(i_am_surface(rank).EQ.1) &
       WRITE(*,*) 'Get the surface. Write it to a file...'
  
 filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr  
  CALL cpu_time(start_time)
  
  CALL get_cartesian_surface(filename, psi, fdrule, 0.0_dp , &
       0.0_dp, 0.0_dp, lmax, rank, .TRUE. )
  
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
  CALL delete_surface3D( rank )
  
  ! Free arrays!
  IF(rank.EQ.0) &
       WRITE(*,*) 'Free arrays!'
  DEALLOCATE(psi,psipro)
  DEALLOCATE(x_ax,y_ax,z_ax)
  DEALLOCATE(xpts,ypts,zpts,xp,yp,zp)
  
  CALL MPI_finalize( ierror )
  
  IF(i_am_surface(rank).EQ.1) THEN
     WRITE(*,*) '¡Se acabó!'
     WRITE(*,*)
     WRITE(*,*)
  ENDIF
  
END PROGRAM cart2sph_ex4
