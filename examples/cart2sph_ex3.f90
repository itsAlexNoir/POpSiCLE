PROGRAM cart2sph_ex2
  
  USE popsicle
  
  IMPLICIT NONE

  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: maxxpts, maxypts, maxzpts
  INTEGER                    :: numxpts, numypts, numzpts
  INTEGER                    :: dims(3)
 
  INTEGER                    :: numproc1dx
  INTEGER                    :: numproc1dy 
  INTEGER                    :: numproc1dz                  
  INTEGER                    :: maxproc1dx
  INTEGER                    :: maxproc1dy 
  INTEGER                    :: maxproc1dz         
  
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
  REAL(dp)                   :: efield(3), afield(3)
  
  CHARACTER(LEN = 150)       :: filename
  CHARACTER(LEN = 6)         :: ctime
  CHARACTER(LEN = 7)         :: rbstr
  CHARACTER(LEN = 5)         :: cprocessor, lmaxstr
  CHARACTER(LEN=100)         :: data_directory

  !-----------------------------------------------!
  
  WRITE(*,*)
  WRITE(*,*) '*******************'
  WRITE(*,*) '  Cart2sph example.'
  WRITE(*,*) '*******************'
  WRITE(*,*) 
  
  WRITE(*,*) 'Opening...'
  
  ! Set the parameters of the grid
  maxxpts          = 30
  maxypts          = 30  
  maxzpts	   = 60

  deltax           = 0.2_dp
  deltay           = 0.2_dp 
  deltaz           = 0.2_dp
  
  numproc1dx       = 4
  numproc1dy       = 4 
  numproc1dz       = 2
  
  maxproc1dx   = numproc1dx - 1
  maxproc1dy   = numproc1dy - 1 
  maxproc1dz   = numproc1dz - 1 
  
  numxpts          = numproc1dx * maxxpts
  numypts          = numproc1dy * maxypts
  numzpts	   = numproc1dz * maxzpts

  dims = (/ numxpts, numypts, numzpts /)
  
  ! The radius of the boundary
  Rboundary    = 8.0_dp
  tolerance = 0.25_dp
  deltar = 0.1_dp
  fdrule = 2
  lmax = 10
  
  data_directory = './data/h2p/cartesian/'
  !WRITE(rbstr,'(F0.3)') Rboundary
  write(rbstr,'(I3.3,F0.3)') INT(Rboundary),Rboundary-INT(Rboundary)
  WRITE(lmaxstr,'(I3.3)') lmax

  !--------------!
  ! Create grids
  !--------------!
  WRITE(*,*) 'Building meshes ...'

  ALLOCATE(x_ax(1:numxpts))
  ALLOCATE(y_ax(1:numypts)) 
  ALLOCATE(z_ax(1:numzpts))
  ALLOCATE(xpts(1:numxpts))
  ALLOCATE(ypts(1:numypts))
  ALLOCATE(zpts(1:numzpts))
  ALLOCATE(xp(1:numxpts))
  ALLOCATE(yp(1:numypts))
  ALLOCATE(zp(1:numzpts))
  
  DO ix = 1, numxpts
     xpt = (-0.5_dp * REAL(numxpts-1,dp) + REAL(ix-1,dp)) * deltax
     x_ax(ix) = xpt
     xpts(ix) = xpt
     xp(ix) = 1.0_dp
  ENDDO
  
  DO iy = 1, numypts
     ypt = (-0.5_dp * REAL(numypts-1,dp) + REAL(iy-1,dp)) * deltay
     y_ax(iy) = ypt
     ypts(iy) = ypt
     yp(iy) = 1.0_dp
  ENDDO
  
  DO iz = 1, numzpts
     zpt = (-0.5_dp * REAL(numzpts-1,dp) + REAL(iz-1,dp)) * deltaz
     z_ax(iz) = zpt
     zpts(iz) = zpt
     zp(iz) = 1.0_dp
  ENDDO
  
  !----------------------------------------!
  ! Initialize cylindrical surface stuff
  !----------------------------------------!
  
  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  CALL initialize_cartesian_surface(xpts, ypts, zpts, dims, &
       Rboundary, tolerance, fdrule, deltar, lmax, .TRUE., filename)
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
  WRITE(*,'(A43,F9.6)') 'Time spent on initializing the surface(s): ',comp_time
  
  !-----------------------------!
  ! Read wavefunction from disk
  !-----------------------------!
  
  WRITE(*,*) 'Allocating wavefunctions...'
  ALLOCATE(psipro(1:maxxpts, 1:maxypts, 1:maxzpts))
  ALLOCATE(psi(1:numxpts, 1:numypts, 1:numzpts))
  
  WRITE(*,*) 'Read data...'
  
  ipro = -1
  
  DO iprocz = 0, maxproc1dz     
     DO iprocy = 0, maxproc1dy
        DO iprocx = 0, maxproc1dx
           
           ipro = ipro + 1
           
           WRITE(cprocessor, '(I5.5)') ipro
           
           filename = TRIM(data_directory) // 'psi/' // cprocessor //    &
                '/psi.' // cprocessor // '.dat'
           
           OPEN(UNIT = 10, FORM = 'unformatted', FILE = filename)
           
           READ(10) psipro
           
           CLOSE(UNIT = 10)
           
           DO iz = 1, maxzpts
              DO iy = 1, maxypts
                 DO ix = 1, maxxpts
                    
                    ixg = iprocx * maxxpts + ix
                    iyg = iprocy * maxypts + iy 
                    izg = iprocz * maxzpts + iz
                    
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
  
  DO iz = 1, numzpts
     DO iy = 1, numypts
        DO ix = 1, numxpts 
           psi(ix, iy, iz) = psi(ix, iy, iz) &
                / SQRT(xp(ix) * yp(iy) * zp(iz) )
        ENDDO
     ENDDO
  ENDDO

  !-------------------------------------------------!
  ! Get cylindrical surface, and write it to a file
  !-------------------------------------------------!
  efield = 0.0_dp
  afield = 0.0_dp
 
  WRITE(*,*) 'Get the surface. Write it to a file...'

  filename = './results/sphfunc.rb' // rbstr // '.lmax' // lmaxstr
  CALL cpu_time(start_time)
  
  CALL get_cartesian_surface(filename, psi, x_ax, y_ax, z_ax, dims, &
       fdrule, 0.0_dp , efield, afield, lmax, .TRUE. )
  
  CALL cpu_time(end_time)
  
  comp_time = end_time - start_time
  
  WRITE(*,'(A38,F9.5)') 'Time spent on getting the surface (s): ',comp_time
  
  !------------------------------!
  ! Now, load it again from file
  !------------------------------!
  WRITE(*,*) 'Deleting surface...'
  CALL delete_surface3D( )
  
  ! Free arrays!
  WRITE(*,*) 'Free arrays!'
  DEALLOCATE(psi,psipro)
  DEALLOCATE(x_ax,y_ax,z_ax)
  DEALLOCATE(xpts,ypts,zpts,xp,yp,zp)
  
  WRITE(*,*) '¡Se acabó!'
  WRITE(*,*)
  WRITE(*,*)
  
END PROGRAM cart2sph_ex2
