MODULE patchwork
  
  USE constants
  USE tools
#if _COM_MPI
 USE MPI
#endif

  IMPLICIT NONE
  
  PRIVATE

  PUBLIC             :: get_subset_cartesian_coordinates
  PUBLIC             :: get_subset_cylindrical_coordinates
  
#if _COM_MPI  
  PUBLIC             :: sum_grid_pieces_xy 
  PUBLIC             :: sum_grid_pieces_xz
  PUBLIC             :: sum_grid_pieces_yz
  PUBLIC             :: sum_grid_pieces_rz
  PUBLIC             :: sum_grid_pieces_xyz
#endif

  INTERFACE get_subset_cartesian_coordinates
     MODULE PROCEDURE get_subset_cartesian_coordinates_serial  
  END INTERFACE get_subset_cartesian_coordinates
  
  INTERFACE get_subset_cylindrical_coordinates
     MODULE PROCEDURE get_subset_cylindrical_coordinates_serial  
  END INTERFACE get_subset_cylindrical_coordinates
  
CONTAINS
  
  !----------------------------------------------------!
  
  
  SUBROUTINE get_subset_cartesian_coordinates_serial(domain, &
       offset,dims,xpts,ypts,zpts,rpts)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: domain(:, :)
    INTEGER, INTENT(OUT)           :: offset(:)
    INTEGER, INTENT(OUT)           :: dims(:)
    REAL(dp), INTENT(IN)           :: xpts(:)
    REAL(dp), INTENT(IN)           :: ypts(:)
    REAL(dp), INTENT(IN), OPTIONAL :: zpts(:)
    REAL(dp), INTENT(IN), OPTIONAL :: rpts(:)

    ! First, find out the offset of the subset
    IF(ALL(xpts.GT.domain(1,1))) THEN
       offset(1) = 0
       IF(ALL(xpts.GT.domain(2,1))) THEN
          dims = 0
       ELSE
          dims(1) = MAXLOC(xpts,1,xpts.LT.domain(2,1))
       ENDIF
    ELSE
       offset(1) = MAXLOC(xpts,1,xpts.LT.domain(1,1))
       dims(1) = MAXLOC(xpts,1,xpts.LT.domain(2,1)) &
            - MAXLOC(xpts,1,xpts.LT.domain(1,1))
    ENDIF
    
    IF(ALL(ypts.GT.domain(1,2))) THEN
       offset(2) = 0
       IF(ALL(ypts.GT.domain(2,2))) THEN
          dims(2) = 0
       ELSE
          dims(2) = MAXLOC(ypts,1,ypts.LT.domain(2,2))
       ENDIF
    ELSE
       offset(2) = MAXLOC(ypts,1,ypts.LT.domain(1,2))
       dims(2) = MAXLOC(ypts,1,ypts.LT.domain(2,2)) &
            - MAXLOC(ypts,1,ypts.LT.domain(1,2))
    ENDIF

    IF(PRESENT(zpts)) THEN
       IF(ALL(zpts.GT.domain(1,3))) THEN
          offset(3) = 0
          IF(ALL(zpts.GT.domain(2,3))) THEN
             dims(3) = 0
          ELSE
             dims(3) = MAXLOC(zpts,1,zpts.LT.domain(2,3))
          ENDIF
       ELSE
          offset(3) = MAXLOC(zpts,1,zpts.LT.domain(1,3))
          dims(3) = MAXLOC(zpts,1,zpts.LT.domain(2,3)) &
               - MAXLOC(zpts,1,zpts.LT.domain(1,3))
       ENDIF
    ENDIF
    
    IF(PRESENT(rpts)) THEN
       IF(ALL(rpts.GT.domain(1,4))) THEN
          offset(4) = 0
          IF(ALL(rpts.GT.domain(2,4))) THEN
             dims(4) = 0
          ELSE
             dims(4) = MAXLOC(rpts,1,rpts.LT.domain(2,4))
          ENDIF
       ELSE
          offset(4) = MAXLOC(rpts,1,rpts.LT.domain(1,4))
          dims(4) = MAXLOC(rpts,1,rpts.LT.domain(2,4)) &
               - MAXLOC(rpts,1,rpts.LT.domain(1,4))
       ENDIF
    ENDIF
    
  END SUBROUTINE get_subset_cartesian_coordinates_serial
  
  !------------------------------------------!

  SUBROUTINE get_subset_cylindrical_coordinates_serial(domain,&
       offset,dims,rhopts,zpts,rpts)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: domain(:, :)
    INTEGER, INTENT(OUT)            :: offset(:)
    INTEGER, INTENT(OUT)            :: dims(:)
    REAL(dp), INTENT(IN)            :: rhopts(:)
    REAL(dp), INTENT(IN)            :: zpts(:)
    REAL(dp), INTENT(IN), OPTIONAL  :: rpts(:)

    !-------------------------------------------------------------!
    
    ! First, find out the offset of the subset
    IF(PRESENT(rpts)) THEN
       IF(ALL(rpts.GT.domain(1,1))) THEN
          offset(1) = 0
          IF(ALL(rpts.GT.domain(2,1))) THEN
             dims(1) = 0
          ELSE
             dims(1) = MAXLOC(rpts,1,rpts.LT.domain(2,1))
          ENDIF
       ELSE
          offset(1) = MAXLOC(rpts,1,rpts.LT.domain(1,1))
          dims(1) = MAXLOC(rpts,1,rpts.LT.domain(2,1)) &
               - MAXLOC(rpts,1,rpts.LT.domain(1,1))
       ENDIF
    ENDIF
    
    IF(ALL(rhopts.GT.domain(1,2))) THEN
       offset(2) = 0
       IF(ALL(rhopts.GT.domain(2,2))) THEN
          dims(2) = 0
       ELSE
          dims(2) = MAXLOC(rhopts,1,rhopts.LT.domain(2,2))
       ENDIF
    ELSE
       offset(2) = MAXLOC(rhopts,1,rhopts.LT.domain(1,2))
       dims(2) = MAXLOC(rhopts,1,rhopts.LT.domain(2,2)) &
            - MAXLOC(rhopts,1,rhopts.LT.domain(1,2))
    ENDIF
    
    IF(ALL(zpts.GT.domain(1,3))) THEN
       offset(3) = 0
       IF(ALL(zpts.GT.domain(2,3))) THEN
          dims(3) = 0
       ELSE
          dims(3) = MAXLOC(zpts,1,zpts.LT.domain(2,3))
       ENDIF
    ELSE
       offset(3) = MAXLOC(zpts,1,zpts.LT.domain(1,3))
       dims(3) = MAXLOC(zpts,1,zpts.LT.domain(2,3)) &
            - MAXLOC(zpts,1,zpts.LT.domain(1,3))
    ENDIF
    
  END SUBROUTINE get_subset_cylindrical_coordinates_serial
  
  !------------------------------------------!
  
  !----------------------------------------------------!
  !----------------------------------------------------!
  
#if _COM_MPI  
  
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
  
#endif
  
END MODULE patchwork
