PROGRAM cy2sph_ex
  
  USE popsicle
  
  IMPLICIT NONE

  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: maxrhopts, maxzpts
  INTEGER                    :: numrhopts, numzpts
  INTEGER, ALLOCATABLE       :: dims(:)
  INTEGER, ALLOCATABLE       :: halopts(:,:)

  INTEGER                    :: numproc1drho       
  INTEGER                    :: numproc1dz                  
  INTEGER                    :: maxproc1drho       
  INTEGER                    :: maxproc1dz         
  
  REAL(dp)                   :: deltarho, deltaz
  REAL(dp)                   :: rboundary
  REAL(dp)                   :: deltar, deltatheta
  INTEGER                    :: lmax
  
  INTEGER                    :: irho, iz
  INTEGER                    :: irg, irhog, izg
  INTEGER                    :: iprocrho, iprocz
  INTEGER                    :: ipro
  
  
  COMPLEX(dp), ALLOCATABLE   :: psipro(:, :)
  COMPLEX(dp), ALLOCATABLE   :: psi(:, :)
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  
  CHARACTER(LEN = 200)       :: filename
  CHARACTER(LEN = 6)         :: ctime
  CHARACTER(LEN = 4)         :: cprocessor
  CHARACTER(LEN=100)         :: data_directory
  
  !-----------------------------------------------!
  
  WRITE(*,*)
  WRITE(*,*) '*******************'
  WRITE(*,*) '  Cy2sph example.'
  WRITE(*,*) '*******************'
  WRITE(*,*) 
  
  WRITE(*,*) 'Opening...'
  
  ! Set the parameters of the grid
  maxrhopts        = 90
  maxzpts	   = 201

  deltarho         = 0.1_dp
  deltaz           = 0.2_dp
  
  numproc1drho     = 4
  numproc1dz       = 10
  
  maxproc1drho = numproc1drho - 1
  maxproc1dz   = numproc1dz - 1 
  
  numrhopts        = numproc1drho * maxrhopts
  numzpts	   = numproc1dz   * maxzpts
  
  ! The radius of the boundary
  rboundary    = 50.0_dp
  
  deltar = 0.1_dp
  deltatheta = 0.1_dp
  lmax             = 10
  
  data_directory = './data/h2p/cylindrical/'
  
  ! Allocate arrays
  WRITE(*,*) 'Allocating arrays...'
  ALLOCATE(psipro(1:maxrhopts, 1:maxzpts))
  ALLOCATE(psi(1:numrhopts, 1:numzpts))
  ALLOCATE(rho_ax(1:numrhopts))
  ALLOCATE(z_ax(1:numzpts))
  
  ! Initialize and make meshes
  WRITE(*,*) 'Building meshes ...'
  
  DO irho = 1, numrhopts
     rho_ax(irho) = REAL(irho,dp) * deltarho
  ENDDO
  
  DO iz = 1, numzpts
     z_ax(iz) = REAL(iz,dp) * deltaz
  ENDDO
  
  ALLOCATE(halopts(2,2))
  halopts(1,1) = -1
  halopts(2,1) = -1
  halopts(1,2) = numrhopts
  halopts(2,2) = numzpts
  
  !----------------------------------------!
  ! Initialize cylindrical grid for interp
  !----------------------------------------!
  CALL initializate_cylindrical_grid(rho_ax,z_ax,(/numrhopts, numzpts/),halopts,&
       Rboundary, 1.0_dp, 2, deltar, lmax)
  
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
  
  !CALL cylindrical2spherical()
  
  
  
  ! Free arrays!
  DEALLOCATE(psi,psipro,rho_ax,z_ax)
  DEALLOCATE(dims, halopts)

END PROGRAM cy2sph_ex
