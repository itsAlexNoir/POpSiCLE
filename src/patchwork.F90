MODULE patchwork
  
  USE constants
  USE tools
#if _COM_MPI
 USE MPI
#endif

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC             :: get_subset_coordinates
#if _COM_MPI  
  PUBLIC             :: sum_grid_pieces_xy 
  PUBLIC             :: sum_grid_pieces_xz
  PUBLIC             :: sum_grid_pieces_yz
  PUBLIC             :: sum_grid_pieces_rz
  PUBLIC             :: sum_grid_pieces_xyz
  ! PUBLIC             :: get_subset
#endif

  INTERFACE get_subset_coordinates
     MODULE PROCEDURE get_subset_coordinates_4Dserial  
     MODULE PROCEDURE get_subset_coordinates_3Dserial
     MODULE PROCEDURE get_subset_coordinates_2Dserial
!#if _COM_MPI
!     MODULE PROCEDURE get_subset_coordinates_3Dparallel
!#endif
  END INTERFACE get_subset_coordinates

CONTAINS

  !----------------------------------------------------!
  
  
  SUBROUTINE get_subset_coordinates_4Dserial(xpts,ypts,zpts,rpts,&
       domain,offset,dims)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: xpts(:)
    REAL(dp), INTENT(IN)        :: ypts(:)
    REAL(dp), INTENT(IN)        :: zpts(:)
    REAL(dp), INTENT(IN)        :: rpts(:)
    REAL(dp), INTENT(IN)        :: domain(:, :)
    INTEGER, INTENT(OUT)        :: offset(:)
    INTEGER, INTENT(OUT)        :: dims(:)
    
    ! First, find out the offset of the subset
    offset(1) = MAXLOC(xpts,1,xpts.LE.domain(1,1))
    offset(2) = MAXLOC(ypts,1,ypts.LE.domain(1,2))
    offset(3) = MAXLOC(zpts,1,zpts.LE.domain(1,3))
    offset(4) = MAXLOC(rpts,1,rpts.LE.domain(1,4))
    
    ! Now, the dimensions of the domain required
    dims(1) = MINLOC(xpts,1,xpts.GE.domain(2,1)) &
         - MAXLOC(xpts,1,xpts.LE.domain(1,1)) - 1
    dims(2) = MINLOC(ypts,1,ypts.GE.domain(2,2)) &
         - MAXLOC(ypts,1,ypts.LE.domain(1,2)) - 1
    dims(3) = MINLOC(zpts,1,zpts.GE.domain(2,3)) &
         - MAXLOC(zpts,1,zpts.LE.domain(1,3)) - 1
    dims(4) = MINLOC(rpts,1,rpts.GE.domain(2,4)) &
         - MAXLOC(rpts,1,rpts.LE.domain(1,4)) - 1
    
  END SUBROUTINE get_subset_coordinates_4Dserial
  
  !------------------------------------------!
  
  SUBROUTINE get_subset_coordinates_3Dserial(xpts,ypts,zpts,domain,&
       offset,dims)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: xpts(:)
    REAL(dp), INTENT(IN)        :: ypts(:)
    REAL(dp), INTENT(IN)        :: zpts(:)
    REAL(dp), INTENT(IN)        :: domain(:, :)
    INTEGER, INTENT(OUT)        :: offset(:)
    INTEGER, INTENT(OUT)        :: dims(:)
    
    ! First, find out the offset of the subset
    offset(1) = MAXLOC(xpts,1,xpts.LE.domain(1,1))
    offset(2) = MAXLOC(ypts,1,ypts.LE.domain(1,2))
    offset(3) = MAXLOC(zpts,1,zpts.LE.domain(1,3))
    
    ! Now, the dimensions of the domain required
    dims(1) = MINLOC(xpts,1,xpts.GE.domain(2,1)) &
         - MAXLOC(xpts,1,xpts.LE.domain(1,1)) - 1
    dims(2) = MINLOC(ypts,1,ypts.GE.domain(2,2)) &
         - MAXLOC(ypts,1,ypts.LE.domain(1,2)) - 1
    dims(3) = MINLOC(zpts,1,zpts.GE.domain(2,3)) &
         - MAXLOC(zpts,1,zpts.LE.domain(1,3)) - 1
    
  END SUBROUTINE get_subset_coordinates_3Dserial
  
  !------------------------------------------!

  SUBROUTINE get_subset_coordinates_2Dserial(xpts,ypts,domain,&
       offset,dims)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: xpts(:)
    REAL(dp), INTENT(IN)        :: ypts(:)
    REAL(dp), INTENT(IN)        :: domain(:, :)
    INTEGER, INTENT(OUT)        :: offset(:)
    INTEGER, INTENT(OUT)        :: dims(:)
    
    ! First, find out the offset of the subset
    offset(1) = MAXLOC(xpts,1,xpts.LE.domain(1,1))
    offset(2) = MAXLOC(ypts,1,ypts.LE.domain(1,2))
    
    ! Now, the dimensions of the domain required
    dims(1) = MINLOC(xpts,1,xpts.GE.domain(2,1)) &
         - MAXLOC(xpts,1,xpts.LE.domain(1,1)) - 1
    dims(2) = MINLOC(ypts,1,ypts.GE.domain(2,2)) &
         - MAXLOC(ypts,1,ypts.LE.domain(1,2)) - 1
    
  END SUBROUTINE get_subset_coordinates_2Dserial
  
  !----------------------------------------------------!
  
#if _COM_MPI  
  !------------------------------------------------------------!
  
  SUBROUTINE sum_grid_pieces_xy(densityproc2D, densitytotal2D, rank, dims_local, &
       dims_global, globalcomm, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: densityproc2D(:, :)
    REAL(dp), INTENT(OUT)    :: densitytotal2D(:, :)
    INTEGER, INTENT(IN)      :: rank
    INTEGER, INTENT(IN)      :: dims_local(:)
    INTEGER, INTENT(IN)      :: dims_global(:)
    INTEGER, INTENT(IN)      :: globalcomm
    INTEGER, INTENT(IN)      :: grid_rank(:)
    INTEGER, INTENT(IN)      :: grid_addresses(:)
    
    !----------------------------------!
    
    INTEGER, ALLOCATABLE       :: members(:)
    INTEGER                    :: simgroup, iogroup, ipro
    INTEGER                    :: iogroup_size, indcomm
    INTEGER, ALLOCATABLE       :: sumcomm(:)
    INTEGER                    :: numprocx, numprocy, numprocz, numprocr
    INTEGER                    :: Nx, Ny, Nz, Nr
    INTEGER                    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                    :: ipx, ipy, ipz, ipr
    INTEGER                    :: iprocx, iprocy, iprocz, iprocr
    INTEGER                    :: ierror
    
    !----------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    iogroup_size = numprocz * numprocr
    ALLOCATE(sumcomm(0:numprocx * numprocy-1))
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    indcomm = -1
    DO iprocy = 0, numprocy - 1 
       DO iprocx = 0, numprocx - 1
          
          indcomm = indcomm + 1
          ipro = 0
          members = 0
          DO iprocr = 0, numprocr - 1
             DO iprocz = 0, numprocz - 1
                ipro = ipro + 1
                members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                     iprocr+1,numprocx,numprocy,numprocz,numprocr))
             ENDDO
          ENDDO
          
          CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
          
          CALL MPI_COMM_CREATE( globalcomm, iogroup, sumcomm(indcomm), ierror)
       ENDDO
    ENDDO
    
    ! Sum all the points in this piece of the xz plane
    indcomm = ipx + numprocx * ipy
    CALL MPI_ALLREDUCE(densityproc2D, densitytotal2D, Nx * Ny, &
         MPI_DOUBLE_PRECISION, MPI_SUM, sumcomm(indcomm), ierror)
    
    DEALLOCATE(members,sumcomm)
    
  END SUBROUTINE sum_grid_pieces_xy
  
  !----------------------------------------------------!
  
  SUBROUTINE sum_grid_pieces_xz(densityproc2D, densitytotal2D, rank, dims_local, &
       dims_global, globalcomm, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: densityproc2D(:, :)
    REAL(dp), INTENT(OUT)    :: densitytotal2D(:, :)
    INTEGER, INTENT(IN)      :: rank
    INTEGER, INTENT(IN)      :: dims_local(:)
    INTEGER, INTENT(IN)      :: dims_global(:)
    INTEGER, INTENT(IN)      :: globalcomm
    INTEGER, INTENT(IN)      :: grid_rank(:)
    INTEGER, INTENT(IN)      :: grid_addresses(:)
    
    !----------------------------------!
    
    INTEGER, ALLOCATABLE       :: members(:)
    INTEGER                    :: simgroup, iogroup, ipro
    INTEGER                    :: iogroup_size, indcomm
    INTEGER, ALLOCATABLE       :: sumcomm(:)
    INTEGER                    :: numprocx, numprocy, numprocz, numprocr
    INTEGER                    :: Nx, Ny, Nz, Nr
    INTEGER                    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                    :: ipx, ipy, ipz, ipr
    INTEGER                    :: iprocx, iprocy, iprocz, iprocr
    INTEGER                    :: ierror
    
    !----------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters(rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    iogroup_size = numprocy * numprocr
    ALLOCATE(sumcomm(0:numprocx * numprocz-1))
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    indcomm = -1
    DO iprocz = 0, numprocz - 1 
       DO iprocx = 0, numprocx - 1
          
          indcomm = indcomm + 1
          ipro = 0
          members = 0
          DO iprocr = 0, numprocr - 1
             DO iprocy = 0, numprocy - 1
                ipro = ipro + 1
                members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                     iprocr+1,numprocx,numprocy,numprocz,numprocr))
             ENDDO
          ENDDO
          
          CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
          
          CALL MPI_COMM_CREATE( globalcomm, iogroup, sumcomm(indcomm), ierror)
       ENDDO
    ENDDO
    
    
    ! Sum all the points in this piece of the xz plane
    indcomm = ipx + numprocx * ipz
    CALL MPI_ALLREDUCE(densityproc2D, densitytotal2D, Nx * Nz, &
         MPI_DOUBLE_PRECISION, MPI_SUM, sumcomm(indcomm), ierror)
    
    DEALLOCATE(members,sumcomm)
    
  END SUBROUTINE sum_grid_pieces_xz
  
  !----------------------------------------------------!

    SUBROUTINE sum_grid_pieces_yz(densityproc2D, densitytotal2D, rank, dims_local, &
       dims_global, globalcomm, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: densityproc2D(:, :)
    REAL(dp), INTENT(OUT)    :: densitytotal2D(:, :)
    INTEGER, INTENT(IN)      :: rank
    INTEGER, INTENT(IN)      :: dims_local(:)
    INTEGER, INTENT(IN)      :: dims_global(:)
    INTEGER, INTENT(IN)      :: globalcomm
    INTEGER, INTENT(IN)      :: grid_rank(:)
    INTEGER, INTENT(IN)      :: grid_addresses(:)
    
    !----------------------------------!
    
    INTEGER, ALLOCATABLE       :: members(:)
    INTEGER                    :: simgroup, iogroup, ipro
    INTEGER                    :: iogroup_size, indcomm
    INTEGER, ALLOCATABLE       :: sumcomm(:)
    INTEGER                    :: numprocx, numprocy, numprocz, numprocr
    INTEGER                    :: Nx, Ny, Nz, Nr
    INTEGER                    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                    :: ipx, ipy, ipz, ipr
    INTEGER                    :: iprocx, iprocy, iprocz, iprocr
    INTEGER                    :: ierror
    
    !----------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters(rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    iogroup_size = numprocx * numprocr
    ALLOCATE(sumcomm(0:numprocy * numprocz-1))
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    indcomm = -1
    DO iprocz = 0, numprocz - 1 
       DO iprocy = 0, numprocy - 1
          
          indcomm = indcomm + 1
          ipro = 0
          members = 0
          DO iprocr = 0, numprocr - 1
             DO iprocx = 0, numprocx - 1
                ipro = ipro + 1
                members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                     iprocr+1,numprocx,numprocy,numprocz,numprocr))
             ENDDO
          ENDDO
          
          CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
          
          CALL MPI_COMM_CREATE( globalcomm, iogroup, sumcomm(indcomm), ierror)
       ENDDO
    ENDDO
    
    ! Sum all the points in this piece of the xz plane
    indcomm = ipy + numprocy * ipz
    CALL MPI_ALLREDUCE(densityproc2D, densitytotal2D, Ny * Nz, &
         MPI_DOUBLE_PRECISION, MPI_SUM, sumcomm(indcomm), ierror)
    
    DEALLOCATE(members,sumcomm)
    
  END SUBROUTINE sum_grid_pieces_yz
  
  !----------------------------------------------------!

    SUBROUTINE sum_grid_pieces_rz(densityproc2D, densitytotal2D, rank, dims_local, &
       dims_global, globalcomm, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: densityproc2D(:, :)
    REAL(dp), INTENT(OUT)    :: densitytotal2D(:, :)
    INTEGER, INTENT(IN)      :: rank
    INTEGER, INTENT(IN)      :: dims_local(:)
    INTEGER, INTENT(IN)      :: dims_global(:)
    INTEGER, INTENT(IN)      :: globalcomm
    INTEGER, INTENT(IN)      :: grid_rank(:)
    INTEGER, INTENT(IN)      :: grid_addresses(:)
    
    !----------------------------------!
    
    INTEGER, ALLOCATABLE       :: members(:)
    INTEGER                    :: simgroup, iogroup, ipro
    INTEGER                    :: iogroup_size, indcomm
    INTEGER, ALLOCATABLE       :: sumcomm(:)
    INTEGER                    :: numprocx, numprocy, numprocz, numprocr
    INTEGER                    :: Nx, Ny, Nz, Nr
    INTEGER                    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                    :: ipx, ipy, ipz, ipr
    INTEGER                    :: iprocx, iprocy, iprocz, iprocr
    INTEGER                    :: ierror
    
    !----------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    iogroup_size = numprocx * numprocy
    ALLOCATE(sumcomm(0:numprocr * numprocz-1))
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    indcomm = -1
    DO iprocr = 0, numprocr - 1 
       DO iprocz = 0, numprocz - 1
          
          indcomm = indcomm + 1
          ipro = 0
          members = 0
          DO iprocy = 0, numprocy - 1
             DO iprocx = 0, numprocx - 1
                ipro = ipro + 1
                members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                     iprocr+1,numprocx,numprocy,numprocz,numprocr))
             ENDDO
          ENDDO
          
          CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
          
          CALL MPI_COMM_CREATE( globalcomm, iogroup, sumcomm(indcomm), ierror)
       ENDDO
    ENDDO
    
    ! Sum all the points in this piece of the xz plane
    indcomm = ipz + numprocz * ipr
    CALL MPI_ALLREDUCE(densityproc2D, densitytotal2D, Nr * Nz, &
         MPI_DOUBLE_PRECISION, MPI_SUM, sumcomm(indcomm), ierror)
    
    DEALLOCATE(members,sumcomm)
    
  END SUBROUTINE sum_grid_pieces_rz
  
  !----------------------------------------------------!

    SUBROUTINE sum_grid_pieces_xyz(densityproc3D, densitytotal3D, rank, dims_local, &
       dims_global, globalcomm, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: densityproc3D(:, :, :)
    REAL(dp), INTENT(OUT)    :: densitytotal3D(:, :, :)
    INTEGER, INTENT(IN)      :: rank
    INTEGER, INTENT(IN)      :: dims_local(:)
    INTEGER, INTENT(IN)      :: dims_global(:)
    INTEGER, INTENT(IN)      :: globalcomm
    INTEGER, INTENT(IN)      :: grid_rank(:)
    INTEGER, INTENT(IN)      :: grid_addresses(:)
    
    !---------------------------------------------------------!
    
    INTEGER, ALLOCATABLE       :: members(:)
    INTEGER                    :: simgroup, iogroup, ipro
    INTEGER                    :: iogroup_size, indcomm
    INTEGER, ALLOCATABLE       :: sumcomm(:)
    INTEGER                    :: numprocx, numprocy, numprocz, numprocr
    INTEGER                    :: Nx, Ny, Nz, Nr
    INTEGER                    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                    :: ipx, ipy, ipz, ipr
    INTEGER                    :: iprocx, iprocy, iprocz, iprocr
    INTEGER                    :: ierror
    
    !---------------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters(rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ! Get the members list of processors to be in the sum
    iogroup_size = numprocr
    ALLOCATE(sumcomm(0:numprocx * numprocy * numprocz-1))
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    indcomm = -1
    DO iprocz = 0, numprocz - 1
       DO iprocy = 0, numprocy - 1
          DO iprocx = 0, numprocx - 1
             
             indcomm = indcomm + 1
             ipro = 0
             members = 0
             DO iprocr = 0, numprocr - 1
                ipro = ipro + 1
                members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                     iprocr+1,numprocx,numprocy,numprocz,numprocr))
             ENDDO
             
             ! create group and communicator
             CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
             
             CALL MPI_COMM_CREATE( globalcomm, iogroup, sumcomm(indcomm), ierror)
          ENDDO
       ENDDO
    ENDDO
    
    ! Sum all the points in this piece of the xz plane
    indcomm = ipx + numprocx * ipy + numprocx * numprocy * ipz
    CALL MPI_ALLREDUCE(densityproc3D, densitytotal3D, Nx * Ny * Nz, &
         MPI_DOUBLE_PRECISION, MPI_SUM, sumcomm(indcomm), ierror)
    
    DEALLOCATE(members,sumcomm)
    
  END SUBROUTINE sum_grid_pieces_xyz
  
  !----------------------------------------------------!
  
!!$  SUBROUTINE get_subset(rank, dims_local, dims_global, &
!!$       grid_rank, grid_ranks_x, grid_ranks_y, grid_ranks_z, &
!!$       grid_ranks_r, new_limits, grid_spacings, &
!!$       globalcomm, newcomm, new_dims_global, new_ipgrid, new_grid_addresses)
!!$    
!!$    IMPLICIT NONE
!!$
!!$    INTEGER, INTENT(IN)      :: rank
!!$    INTEGER, INTENT(IN)      :: dims_local(:)
!!$    INTEGER, INTENT(IN)      :: dims_global(:)    
!!$    INTEGER, INTENT(IN)      :: grid_rank(:)
!!$    INTEGER, INTENT(IN)      :: grid_ranks_x(:)
!!$    INTEGER, INTENT(IN)      :: grid_ranks_y(:)
!!$    INTEGER, INTENT(IN)      :: grid_ranks_z(:)
!!$    INTEGER, INTENT(IN)      :: grid_ranks_r(:)
!!$    REAL(dp), INTENT(IN)     :: grid_spacings(:)
!!$    REAL(dp), INTENT(IN)     :: new_limits(:, :)
!!$    INTEGER, INTENT(IN)      :: globalcomm
!!$    INTEGER, INTENT(OUT)     :: newcomm
!!$    INTEGER, INTENT(OUT)     :: new_dims_global(:)
!!$    INTEGER, INTENT(OUT)     :: new_ipgrid(:)
!!$    INTEGER, ALLOCATABLE, INTENT(OUT) :: new_grid_addresses(:)
!!$    
!!$    !----------------------------------!
!!$    
!!$    INTEGER, ALLOCATABLE     :: members(:)
!!$    INTEGER                  :: simgroup, iogroup, ipro
!!$    INTEGER                  :: np, ip, newrank
!!$    INTEGER                  :: numprocx, numprocy, numprocz, numprocr
!!$    INTEGER                  :: newnumprocx, newnumprocy
!!$    INTEGER                  :: newnumprocz, newnumprocr
!!$    INTEGER                  :: Nx, Ny, Nz, Nr
!!$    INTEGER                  :: Nxgl, Nygl, Nzgl, Nrgl
!!$    INTEGER                  :: ipx, ipy, ipz, ipr
!!$    INTEGER                  :: iprocx, iprocy, iprocz, iprocr
!!$    INTEGER                  :: minprocx, maxprocx
!!$    INTEGER                  :: minprocy, maxprocy
!!$    INTEGER                  :: minprocz, maxprocz
!!$    INTEGER                  :: minprocr, maxprocr
!!$    REAL(dp)                 :: dx, dy, dz, dr
!!$    REAL(dp)                 :: xmin, xmax, x0min, x0max
!!$    REAL(dp)                 :: ymin, ymax, y0min, y0max
!!$    REAL(dp)                 :: zmin, zmax, z0min, z0max
!!$    REAL(dp)                 :: rmin, rmax, r0min, r0max
!!$    INTEGER                  :: ierror
!!$    
!!$    !----------------------------------!
!!$    
!!$    
!!$    ! Assign dimensions
!!$    CALL get_simulation_parameters( rank, &
!!$         dims_local, dims_global, grid_rank, &
!!$         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
!!$         ipx, ipy, ipz, ipr, &
!!$         numprocx, numprocy, numprocz, numprocr )
!!$    
!!$    IF(rank .EQ. 1) THEN
!!$       
!!$       x0min = new_limits(1,1)
!!$       x0max = new_limits(2,1)
!!$
!!$       y0min = 0.0_dp
!!$       y0max = 1.0_dp
!!$       
!!$       z0min = 0.0_dp
!!$       z0max = 1.0_dp
!!$       
!!$       r0min = 0.0_dp
!!$       r0max = 1.0_dp
!!$       
!!$       dx = grid_spacings(1)
!!$       dy = 1.0_dp
!!$       dz = 1.0_dp
!!$       dr = 1.0_dp
!!$       
!!$    ELSEIF(rank .EQ. 2) THEN
!!$       
!!$       x0min = new_limits(1,1)
!!$       x0max = new_limits(2,1)
!!$       
!!$       y0min = new_limits(1,2)
!!$       y0max = new_limits(2,2)
!!$       
!!$       z0min = 0.0_dp
!!$       z0max = 1.0_dp
!!$       
!!$       r0min = 0.0_dp
!!$       r0max = 1.0_dp
!!$       
!!$       dx = grid_spacings(1)
!!$       dy = grid_spacings(2)
!!$       dz = 1.0_dp
!!$       dr = 1.0_dp
!!$       
!!$    ELSEIF(rank .EQ. 3) THEN
!!$       
!!$       x0min = new_limits(1,1)
!!$       x0max = new_limits(2,1)
!!$       
!!$       y0min = new_limits(1,2)
!!$       y0max = new_limits(2,2)
!!$       
!!$       z0min = new_limits(1,3)
!!$       z0max = new_limits(2,3)
!!$       
!!$       r0min = 0.0_dp
!!$       r0max = 1.0_dp
!!$       
!!$       dx = grid_spacings(1)
!!$       dy = grid_spacings(2)
!!$       dz = grid_spacings(3)
!!$       dr = 1.0_dp
!!$       
!!$    ELSEIF(rank .EQ. 4) THEN
!!$       
!!$       x0min = new_limits(1,1)
!!$       x0max = new_limits(2,1)
!!$       
!!$       y0min = new_limits(1,2)
!!$       y0max = new_limits(2,2)
!!$       
!!$       z0min = new_limits(1,3)
!!$       z0max = new_limits(2,3)
!!$       
!!$       r0min = new_limits(1,4)
!!$       r0max = new_limits(2,4)
!!$       
!!$       dx = grid_spacings(1)
!!$       dy = grid_spacings(2)
!!$       dz = grid_spacings(3)
!!$       dr = grid_spacings(4)
!!$       
!!$    ENDIF
!!$    
!!$    ipro = -1
!!$    
!!$    DO iprocr = 1, numprocr - 1
!!$       DO iprocz = 1, numprocz - 1
!!$          DO iprocy = 1, numprocy - 1
!!$             DO iprocx = 1, numprocx - 1
!!$                
!!$                ipro = ipro + 1
!!$                
!!$                xmin = REAL(grid_ranks_x(ipro) * Nx,dp) * dx
!!$                xmax = REAL((grid_ranks_x(ipro) + 1) * Nx,dp) * dx
!!$                
!!$                ymin = REAL(grid_ranks_y(ipro) * Ny,dp) * dy 
!!$                ymax = REAL((grid_ranks_y(ipro) + 1) * Ny,dp ) * dy 
!!$                
!!$                zmin = REAL(grid_ranks_z(ipro) * Nz, dp) * dz
!!$                zmax = REAL((grid_ranks_z(ipro) + 1) * Nz, dp) * dz 
!!$                
!!$                rmin = REAL(grid_ranks_r(ipro) * Nr, dp ) * dr
!!$                rmax = REAL((grid_ranks_r(ipro) + 1) * Nr, dp) * dr 
!!$                
!!$                
!!$                IF( (x0min.GE.xmin .AND. x0min.LE.xmax .OR. & 
!!$                     x0min.LE.xmin .AND. x0max.GE.xmax .OR. &
!!$                     x0max .GE. xmin .AND. x0max .LE. xmax) .AND. &
!!$                     (y0min.GE.ymin .AND. y0min.LE.ymax .OR. & 
!!$                     y0min.LE.ymin .AND. y0max.GE.ymax .OR. &
!!$                     y0max .GE. ymin .AND. y0max .LE. ymax) .AND. &
!!$                     (z0min.GE.zmin .AND. z0min.LE.zmax .OR. & 
!!$                     z0min.LE.zmin .AND. z0max.GE.zmax .OR. &
!!$                     z0max .GE. zmin .AND. z0max .LE. zmax) .AND. &
!!$                     (r0min.GE.rmin .AND. r0min.LE.rmax .OR. & 
!!$                     r0min.LE.rmin .AND. r0max.GE.rmax .OR. &
!!$                     r0max .GE. rmin .AND. r0max .LE. rmax) ) THEN
!!$                   
!!$                   np = np + 1
!!$                   
!!$                   minprocx = MIN(grid_ranks_x(ipro),minprocx)
!!$                   maxprocx = MIN(grid_ranks_x(ipro),maxprocx)
!!$                   
!!$                   minprocy = MIN(grid_ranks_y(ipro),minprocy)
!!$                   maxprocy = MIN(grid_ranks_y(ipro),maxprocy)
!!$                   
!!$                   minprocz = MIN(grid_ranks_z(ipro),minprocz)
!!$                   maxprocz = MIN(grid_ranks_z(ipro),maxprocz)
!!$                   
!!$                   minprocr = MIN(grid_ranks_r(ipro),minprocr)
!!$                   maxprocr = MIN(grid_ranks_r(ipro),maxprocr)
!!$                   
!!$                ENDIF
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    ALLOCATE(members(1:np))
!!$    ALLOCATE(new_grid_addresses(1:np))
!!$    
!!$    newnumprocx =  ((maxprocx + 1) - minprocx)
!!$    newnumprocy =  ((maxprocy + 1) - minprocy)
!!$    newnumprocz =  ((maxprocz + 1) - minprocz)
!!$    newnumprocr =  ((maxprocr + 1) - minprocr)
!!$    
!!$    new_dims_global = 1
!!$    
!!$    new_dims_global(1) = newnumprocx / Nx
!!$    IF(rank.GE. 2) &
!!$         new_dims_global(2) = newnumprocy / Ny
!!$    IF(rank .GE. 3) &
!!$         new_dims_global(3) = newnumprocz / Nz
!!$    IF(rank.GE. 4) &
!!$         new_dims_global(4) = newnumprocr / Nr
!!$    
!!$    
!!$    DO iprocr = 1, numprocr - 1
!!$       DO iprocz = 1, numprocz - 1
!!$          DO iprocy = 1, numprocy - 1
!!$             DO iprocx = 1, numprocx - 1
!!$                
!!$                ipro = ipro + 1
!!$                
!!$                xmin = REAL(grid_ranks_x(ipro) * Nx,dp) * dx
!!$                xmax = REAL((grid_ranks_x(ipro) + 1) * Nx,dp) * dx
!!$                
!!$                ymin = REAL(grid_ranks_y(ipro) * Ny,dp) * dy 
!!$                ymax = REAL((grid_ranks_y(ipro) + 1) * Ny,dp ) * dy 
!!$                
!!$                zmin = REAL(grid_ranks_z(ipro) * Nz, dp) * dz
!!$                zmax = REAL((grid_ranks_z(ipro) + 1) * Nz, dp) * dz 
!!$                
!!$                rmin = REAL(grid_ranks_r(ipro) * Nr, dp ) * dr
!!$                rmax = REAL((grid_ranks_r(ipro) + 1) * Nr, dp) * dr 
!!$                
!!$                IF( (x0min.GE.xmin .AND. x0min.LE.xmax .OR. & 
!!$                     x0min.LE.xmin .AND. x0max.GE.xmax .OR. &
!!$                     x0max .GE. xmin .AND. x0max .LE. xmax) .AND. &
!!$                     (y0min.GE.ymin .AND. y0min.LE.ymax .OR. & 
!!$                     y0min.LE.ymin .AND. y0max.GE.ymax .OR. &
!!$                     y0max .GE. ymin .AND. y0max .LE. ymax) .AND. &
!!$                     (z0min.GE.zmin .AND. z0min.LE.zmax .OR. & 
!!$                     z0min.LE.zmin .AND. z0max.GE.zmax .OR. &
!!$                     z0max .GE. zmin .AND. z0max .LE. zmax) .AND. &
!!$                     (r0min.GE.rmin .AND. r0min.LE.rmax .OR. & 
!!$                     r0min.LE.rmin .AND. r0max.GE.rmax .OR. &
!!$                     r0max .GE. rmin .AND. r0max .LE. rmax) ) THEN
!!$                   
!!$                   DO ip = 1, np
!!$                      members(ip) = ipro
!!$                      new_grid_addresses(ip) = ipro                   
!!$                   ENDDO
!!$                   
!!$                ENDIF
!!$                
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    ! Create a global group, used for I/O
!!$    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
!!$    
!!$    CALL MPI_GROUP_INCL(simgroup, np, members, iogroup, ierror)
!!$    
!!$    CALL MPI_COMM_CREATE( globalcomm, iogroup, newcomm, ierror)
!!$    
!!$    CALL MPI_COMM_RANK( newcomm, newrank, ierror)
!!$    
!!$    new_ipgrid = 1
!!$    ipro = -1
!!$    
!!$    DO iprocr = 1, newnumprocr - 1
!!$       DO iprocz = 1, newnumprocz - 1
!!$          DO iprocy = 1, newnumprocy - 1
!!$             DO iprocx = 1, newnumprocx - 1
!!$                ipro = ipro + 1
!!$                IF (newrank .EQ. ipro) THEN
!!$                   new_ipgrid(1) = iprocx
!!$                   IF(rank.GE. 2) &
!!$                        new_ipgrid(2) = iprocy
!!$                   IF(rank.GE. 3) &
!!$                        new_ipgrid(3) = iprocz
!!$                   IF(rank.GE. 4) &
!!$                        new_ipgrid(4) = iprocr
!!$                ENDIF
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    
!!$    DEALLOCATE(members)
!!$    
!!$  END SUBROUTINE get_subset
!!$  
  !----------------------------------------------------!
  
#endif
  
END MODULE patchwork
