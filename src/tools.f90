MODULE tools
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC        :: indx
  PUBLIC        :: get_simulation_parameters
  
CONTAINS
  
  FUNCTION indx(ix, iy, iz, ir, Nx, Ny, Nz, Nr)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)      :: ix, iy, iz, ir
    INTEGER, INTENT(IN)      :: Nx, Ny, Nz, Nr
    INTEGER                  :: indx
    
    indx = Nx*Ny*Nz*(ir-1) + Nx*Ny*(iz-1) &
         + Nx*(iy-1) + ix
    
  END FUNCTION indx
  
  !------------------------------------------!
  
  SUBROUTINE get_simulation_parameters( rank, &
       dims_local, dims_global, grid_rank, &
       Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
       ipx, ipy, ipz, ipr, &
       numprocx, numprocy, numprocz, numprocr )
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: rank
    INTEGER, INTENT(IN)     :: dims_local(:)
    INTEGER, INTENT(IN)     :: dims_global(:)
    INTEGER, INTENT(IN)     :: grid_rank(:)
    INTEGER, INTENT(OUT)    :: Nx, Ny, Nz, Nr
    INTEGER, INTENT(OUT)    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER, INTENT(OUT)    :: ipx, ipy, ipz, ipr
    INTEGER, INTENT(OUT)    :: numprocx, numprocy
    INTEGER, INTENT(OUT)    :: numprocz, numprocr
    
    !-----------------------------------------------!
    
    IF (rank .EQ. 1) THEN
       Nx = dims_local(1)
       Ny = 1
       Nz = 1
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = 1
       Nzgl = 1
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = 0
       ipz = 0
       ipr = 0
       
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = 1
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = 1
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = 0
       ipr = 0
       
    ELSEIF (rank .EQ. 3) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = grid_rank(3)
       ipr = 0
       
    ELSEIF (rank .EQ. 4) THEN
       
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = dims_local(4)
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = dims_global(4)
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = grid_rank(3)
       ipr = grid_rank(4)
       
    ELSE
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF
    
    ! Work out the number of processors per dimension
    numprocx = Nxgl / Nx
    numprocy = Nygl / Ny
    numprocz = Nzgl / Nz
    numprocr = Nrgl / Nr
    
  END SUBROUTINE get_simulation_parameters
  
  
END MODULE tools
