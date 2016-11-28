MODULE io_subset

  USE HDF5
  USE constants_pop
  USE tools
  USE patchwork
#if _COM_MPI
  USE MPI
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC               :: write_subset
  PUBLIC               :: read_subset

  
  INTERFACE write_subset
#if _COM_MPI     
     MODULE PROCEDURE write_subset_parallel
#endif
  END INTERFACE write_subset
  
  INTERFACE read_subset
     MODULE PROCEDURE dread_subset_serial
#if _COM_MPI
     MODULE PROCEDURE dread_subset_parallel
#endif
  END INTERFACE read_subset
  
CONTAINS

      !-------------------------------------------------------------------!
#if _COM_MPI
  
  SUBROUTINE write_subset_parallel(wave, domain, &
       rank, dims_local, dims_global, &
       globalcomm, grid_rank,grid_addresses, filename, &
       xpts, ypts, zpts, rpts)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)         :: wave(:)
    REAL(dp), INTENT(IN)         :: domain(:, :)
    INTEGER, INTENT(IN)          :: rank
    INTEGER, INTENT(IN)          :: dims_local(:)
    INTEGER, INTENT(IN)          :: dims_global(:)
    INTEGER, INTENT(IN)          :: globalcomm
    INTEGER, INTENT(IN)          :: grid_rank(:)
    INTEGER, INTENT(IN)          :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename 
    REAL(dp), INTENT(IN), OPTIONAL :: xpts(:), ypts(:)
    REAL(dp), INTENT(IN), OPTIONAL :: zpts(:), rpts(:)
    
    !! Shaped arrrays
    REAL(dp), ALLOCATABLE        :: wave2D(:, :)
    REAL(dp), ALLOCATABLE        :: wave3D(:, :, :)
    REAL(dp), ALLOCATABLE        :: wave4D(:, :, :, :)  
    REAL(dp), ALLOCATABLE        :: wave2D_redux(:, :)
    REAL(dp), ALLOCATABLE        :: wave3D_redux(:, :, :)
    REAL(dp), ALLOCATABLE        :: wave4D_redux(:, :, :, :)  
    INTEGER, ALLOCATABLE         :: i_am_subset_local(:, :, :, :)
    INTEGER, ALLOCATABLE         :: i_am_subset(:, :, :, :)
    INTEGER, ALLOCATABLE         :: my_subset_dims(:, :, :, :, :)
    INTEGER, ALLOCATABLE         :: subset_dims_global(:, :, :, :, :)
    INTEGER, ALLOCATABLE         :: local_offset(:)
    INTEGER, ALLOCATABLE         :: subset_dims(:)
    INTEGER, ALLOCATABLE         :: total_dims(:)
    INTEGER                      :: numtotalprocs
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ix,iy,iz,ir
    INTEGER                      :: ierror
    
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Dataspace identifier in file 
    INTEGER(HID_T)                 :: filespace
    ! Property list identifier 
    INTEGER(HID_T)                 :: plist_id
    ! Dataspace identifier in memory
    INTEGER(HID_T)                 :: memspace
    ! Group identifiers
    INTEGER(HID_T)                 :: group_id 
    ! Dataspace identifier
    INTEGER(HID_T)                 :: dspace_id    
    ! Dataset identifier  
    INTEGER(HID_T)                 :: dset_id      
    ! Dataset dimensions  
    INTEGER(HSIZE_T), ALLOCATABLE  :: file_dims(:)
    INTEGER(HSIZE_T), ALLOCATABLE  :: proc_dims(:)
    ! Memory variables, per proc
    INTEGER(HSIZE_T), ALLOCATABLE  :: cont(:)  
    INTEGER(HSSIZE_T), ALLOCATABLE :: offset(:) 
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name
    CHARACTER(LEN = 100)           :: setname
    
    !---------------------------------!
    
    IF(rank.EQ.2) THEN
       ALLOCATE(wave2D(1:dims_local(1),1:dims_local(2)))
       wave2D = RESHAPE(wave,SHAPE(wave2D))
    ELSEIF(rank.EQ.3) THEN
       ALLOCATE(wave3D(1:dims_local(1),1:dims_local(2), &
            1:dims_local(3)))
       wave3D = RESHAPE(wave,SHAPE(wave3D))
    ELSEIF(rank.EQ.4) THEN
       ALLOCATE(wave4D(1:dims_local(1),1:dims_local(2), &
            1:dims_local(3),1:dims_local(4)))
       wave4D = RESHAPE(wave,SHAPE(wave4D))
    ELSE 
       WRITE(*,*) 'Number of dimensions not implemented!'
       STOP
    ENDIF
    
    ALLOCATE(local_offset(1:rank))
    ALLOCATE(subset_dims(1:rank))
    ALLOCATE(total_dims(1:rank)) 
    ALLOCATE(file_dims(1:rank))
    ALLOCATE(proc_dims(1:rank))
    ALLOCATE(offset(1:rank))
    ALLOCATE(cont(1:rank))
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    numtotalprocs = numprocx * numprocy * numprocz * numprocr
    
    ALLOCATE(i_am_subset_local(0:numprocx-1,0:numprocy-1,0:numprocz-1,0:numprocr-1))
    ALLOCATE(i_am_subset(0:numprocx-1,0:numprocy-1,0:numprocz-1,0:numprocr-1))
    i_am_subset_local = 0
    i_am_subset = 0
    ALLOCATE(my_subset_dims(0:numprocx-1,0:numprocy-1,0:numprocz-1,0:numprocr-1,1:rank))
    ALLOCATE(subset_dims_global(0:numprocx-1,0:numprocy-1,0:numprocz-1,0:numprocr-1,1:rank))
    my_subset_dims = 0
    subset_dims_global = 0
    
    ! Check if the processor are within the domain
    IF(rank.EQ.2) THEN
       CALL get_subset_cartesian_coordinates(domain,local_offset,subset_dims,&
            xpts(1:dims_local(1)),ypts(1:dims_local(2)))
    ELSEIF(rank.EQ.3) THEN
       CALL get_subset_cartesian_coordinates(domain,local_offset,subset_dims,&
            xpts(1:dims_local(1)),ypts(1:dims_local(2)), &
            zpts(1:dims_local(3)))
    ELSEIF(rank.EQ.4) THEN
       CALL get_subset_cartesian_coordinates(domain,local_offset,subset_dims,&
            xpts(1:dims_local(1)),ypts(1:dims_local(2)), &
            zpts(1:dims_local(3)),rpts(1:dims_local(4)))
    ELSE
       WRITE(*,*) 'Number of dimensions not implemented!'
       STOP
    ENDIF
    
    IF(ANY(subset_dims.EQ.0)) THEN
       i_am_subset_local(ipx,ipy,ipz,ipr) = 0
    ELSE
       i_am_subset_local(ipx,ipy,ipz,ipr) = 1
    ENDIF
    
    CALL MPI_ALLREDUCE(i_am_subset_local,i_am_subset,numtotalprocs, &
         MPI_INTEGER,MPI_SUM,globalcomm,ierror)
    
    my_subset_dims(ipx,ipy,ipz,ipr,:) = subset_dims(:)
    
    CALL MPI_ALLREDUCE(my_subset_dims,subset_dims_global,numtotalprocs*rank, &
         MPI_INTEGER,MPI_SUM,globalcomm,ierror)
    
    iogroup_size = SUM(i_am_subset)
    ALLOCATE(members(1:iogroup_size))
    
    ipro = 0
    DO iprocr = 0, numprocr - 1
       DO iprocz = 0, numprocz - 1
          DO iprocy = 0, numprocy - 1
             DO iprocx = 0, numprocx - 1
                IF(i_am_subset(iprocx,iprocy,iprocz,iprocr).EQ.1) THEN
                   ipro = ipro + 1
                   members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                        iprocr+1,numprocx,numprocy,numprocz,numprocr))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! Create a global group, used for I/O
    CALL MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    ! Open and write to the file only the procs that contribute
    IF(i_am_subset(ipx,ipy,ipz,ipr).EQ.1) THEN
       
       ! Obtain the piece of the data needed
       IF(rank.EQ.2) THEN
          ALLOCATE(wave2D_redux(1:subset_dims(1),1:subset_dims(2)))
          DO iy = 1, subset_dims(2)
             DO ix = 1, subset_dims(1)
                wave2D_redux(ix,iy) = wave2D(local_offset(1)+ix,local_offset(2)+iy)
             ENDDO
          ENDDO
          
       ELSEIF(rank.EQ.3) THEN
          ALLOCATE(wave3D_redux(1:subset_dims(1),1:subset_dims(2),1:subset_dims(3)))
          
          DO iz = 1, subset_dims(3)
             DO iy = 1, subset_dims(2)
                DO ix = 1, subset_dims(1)
                   wave3D_redux(ix,iy,iz) = wave3D(local_offset(1)+ix,local_offset(2)+iy, &
                        local_offset(3)+iz)
                ENDDO
             ENDDO
          ENDDO
          
       ELSEIF(rank.EQ.4) THEN
          ALLOCATE(wave4D_redux(1:subset_dims(1),1:subset_dims(2),1:subset_dims(3),&
               1:subset_dims(4)))
          
          DO ir = 1, subset_dims(4)
             DO iz = 1, subset_dims(3)
                DO iy = 1, subset_dims(2)
                   DO ix = 1, subset_dims(1)
                      wave4D_redux(ix,iy,iz,ir) = wave4D(local_offset(1)+ix,local_offset(2)+iy, &
                           local_offset(3)+iz,local_offset(4)+ir)
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          
       ENDIF

       ! Initialize fortran interface
       CALL h5open_f(error)
       
       ! Setup file access property list with 
       ! parallel I/O access.
       CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
       CALL h5pset_fapl_mpio_f(plist_id, iocomm, MPI_INFO_NULL, error)
       
       
       ! Create the file collectively.
       name = TRIM(filename) // '.h5'  
       CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
       CALL h5pclose_f(plist_id, error)
       
       ! Create group in the root group.     
       name = '/' // TRIM(filename)
       !write(*,*) 'name: ',name
       CALL h5gcreate_f(file_id, filename, group_id, error)
       
       ! Assignate dims to the file
       file_dims = 0
       proc_dims = 0
       total_dims = 0
       
       IF(rank.GE.2) THEN
          DO iprocx = 0, numprocx - 1
             total_dims(1) = total_dims(1) + subset_dims_global(iprocx,ipy,ipz,ipr,1)
          ENDDO
          DO iprocy = 0, numprocy - 1
             total_dims(2) = total_dims(2) + subset_dims_global(ipx,iprocy,ipz,ipr,2)
          ENDDO
       ENDIF
       IF(rank.GE.3) THEN
          DO iprocz = 0, numprocz - 1
             total_dims(3) = total_dims(3) + subset_dims_global(ipx,ipy,iprocz,ipr,3)
          ENDDO
       ENDIF
       IF(rank.GE.4) THEN
          DO iprocr = 0, numprocr - 1
             total_dims(4) = total_dims(4) + subset_dims_global(ipx,ipy,ipz,iprocr,4)
          ENDDO
       ENDIF

       proc_dims = subset_dims
       file_dims = total_dims
       
       ! Create the data space for the  dataset.
       CALL h5screate_simple_f(rank, file_dims, filespace, error)
       CALL h5screate_simple_f(rank, proc_dims, memspace, error)
       
       
       setname = '/' // TRIM(filename) // '/'// TRIM(filename)
       !write(*,*) 'setname: ',setname
       CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, filespace, &
            dset_id, error)!, plist_id)
       CALL h5sclose_f(filespace, error)
       
       ! Each process defines dataset in memory and writes it to the hyperslab
       ! in the file. 
       cont  = proc_dims
       offset = 0
       
       IF(rank.GE.2) THEN
          DO iprocx = 0, ipx - 1
             offset(1) = offset(1) + subset_dims_global(iprocx,ipy,ipz,ipr,1)
          ENDDO
          
          DO iprocy = 0, ipy - 1
             offset(2) = offset(2) + subset_dims_global(ipx,iprocy,ipz,ipr,2)
          ENDDO
       ENDIF
       
       IF(rank.GE.3) THEN
          DO iprocz = 0, ipz - 1
             offset(3) = offset(3) + subset_dims_global(ipx,ipy,iprocz,ipr,3)
          ENDDO
       ENDIF
       
       IF(rank.GE.4) THEN
          DO iprocr = 0, ipr - 1
             offset(4) = offset(4) + subset_dims_global(ipx,ipy,ipz,iprocr,4)
          ENDDO
       ENDIF
       
       ! Select hyperslab in the file.
       CALL h5dget_space_f(dset_id, filespace, error)
       CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset,&
            cont, error)!, stride, blck)
       
       ! Create property list for collective dataset write
       CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
       CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
       
       ! Write the dataset collectively.
       IF (rank.EQ.2) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave2D_redux, file_dims, &
               error, file_space_id = filespace, mem_space_id = memspace, &
               xfer_prp = plist_id)
       ELSEIF(rank.EQ.3) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave3D_redux, file_dims, &
               error, file_space_id = filespace, mem_space_id = memspace, &
               xfer_prp = plist_id)
       ELSEIF(rank.EQ.4) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave4D_redux, file_dims, &
               error, file_space_id = filespace, mem_space_id = memspace, &
               xfer_prp = plist_id)
       ELSE
          WRITE(*,*) 'Number of dimensions not implemented!'
          STOP
       ENDIF
       
       ! Write the dataset independently. 
       !    CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dimsfi,error, &
       !                     file_space_id = filespace, mem_space_id = memspace)
       
       ! Close dataspaces.
       CALL h5sclose_f(filespace, error)
       CALL h5sclose_f(memspace, error)
       
       ! Close the dataset.
       CALL h5dclose_f(dset_id, error)
       
       ! Close the property list.
       CALL h5pclose_f(plist_id, error)
       
       ! Close the group.
       CALL h5gclose_f(group_id, error)
       
       ! Close the file.
       CALL h5fclose_f(file_id, error)
       
       ! Close FORTRAN interfaces and HDF5 library.
       CALL h5close_f(error)

       IF(rank.EQ.2) THEN
          DEALLOCATE(wave2D_redux)
       ELSEIF(rank.EQ.3) THEN
          DEALLOCATE(wave3D_redux)
       ELSEIF(rank.EQ.4) THEN
          DEALLOCATE(wave4D_redux)
       ENDIF

    ENDIF
    
    ! Release communicator and groups
    IF(i_am_subset(ipx,ipy,ipz,ipr).EQ.1) &
         CALL MPI_COMM_FREE(iocomm, ierror)
    CALL MPI_GROUP_FREE(iogroup, ierror)
    CALL MPI_GROUP_FREE(simgroup, ierror)
    
    ! Deallocate stuff
    DEALLOCATE(i_am_subset_local,i_am_subset)
    DEALLOCATE(my_subset_dims,subset_dims_global)
    DEALLOCATE(members)
    DEALLOCATE(local_offset, offset)
    DEALLOCATE(subset_dims, total_dims, cont)
    DEALLOCATE(file_dims, proc_dims)
    IF(rank.EQ.2) THEN
       DEALLOCATE(wave2D)
    ELSEIF(rank.EQ.3) THEN
       DEALLOCATE(wave3D)
    ELSEIF(rank.EQ.4) THEN
       DEALLOCATE(wave4D)
    ENDIF
    
  END SUBROUTINE write_subset_parallel
#endif
  !-------------------------------------------------------------------!

  SUBROUTINE dread_subset_serial(filename_in, psi, rank, dims, offset)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)      :: filename_in
    REAL(dp), INTENT(OUT)             :: psi(:)
    INTEGER, INTENT(IN)               :: rank
    INTEGER, INTENT(IN)               :: dims(:)
    INTEGER, INTENT(IN), OPTIONAL     :: offset(:)

    REAL(dp), ALLOCATABLE             :: psi2D(:, :)
    REAL(dp), ALLOCATABLE             :: psi3D(:, :, :)
    REAL(dp), ALLOCATABLE             :: psi4D(:, :, :, :)  
    ! File identifier
    INTEGER(HID_T)                    :: file_id
    ! Group identifier
    INTEGER(HID_T)                    :: grp_id
    ! Dataset identifier
    INTEGER(HID_T)                    :: dset_id 
    ! Dataspace identifier
    INTEGER(HID_T)                    :: dataspace 
    ! dset size
    INTEGER(HSIZE_T)                  :: subset_size(rank)
    ! Memspace identifier
    INTEGER(HID_T)                    :: memspace
    ! Hyperslab offset
    INTEGER(HSIZE_T)                  :: subset_offset(rank)
    ! Error flag
    INTEGER                           :: error
    CHARACTER(LEN=100)                :: filename
    CHARACTER(LEN=100)                :: dsetname, groupname
    
    !----------------------------------------------!
    
    ! Initialize FORTRAN interface. 
    CALL h5open_f(error) 

    IF(rank.EQ.2) THEN
       ALLOCATE(psi2D(1:dims(1),1:dims(2)))
    ELSEIF(rank.EQ.3) THEN
       ALLOCATE(psi3D(1:dims(1),1:dims(2),1:dims(3)))
    ELSEIF(rank.EQ.4) THEN
       ALLOCATE(psi4D(1:dims(1),1:dims(2),1:dims(3), &
            1:dims(4)))
    ELSE
       WRITE(*,*) 'No routine for your rank!!'
       STOP
    ENDIF
    
    ! Open the file.
    filename = TRIM(filename_in) // '.h5'
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

    ! Open the group.
    groupname = '/' // TRIM(filename_in)
    CALL h5gopen_f(file_id, groupname, grp_id, error)
    
    ! Open the dataset.
    dsetname = '/' // TRIM(filename_in) // '/' // TRIM(filename_in)
    CALL h5dopen_f(file_id, dsetname, dset_id, error)
    
    ! Get dataset's dataspace identifier and select subset.
    IF(PRESENT(offset)) THEN
       subset_offset = offset
    ELSE
       subset_offset = 0
    ENDIF
    subset_size = dims
    CALL h5dget_space_f(dset_id, dataspace, error)
    CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
         subset_offset, subset_size, error) 
    
    ! Create memspace
    CALL h5screate_simple_f(rank,subset_size,memspace,error)
    
    ! Read the dataset
    IF(rank.EQ.2) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi2D, subset_size, error, &
            memspace, dataspace)
       psi = RESHAPE(psi2D,(/PRODUCT(dims)/))
       DEALLOCATE(psi2D)
    ELSEIF(rank.EQ.3) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi3D, subset_size, error, &
            memspace, dataspace)
       psi = RESHAPE(psi3D,(/PRODUCT(dims)/))
       DEALLOCATE(psi3D)
    ELSEIF(rank.EQ.4) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi4D, subset_size, error, &
            memspace, dataspace)
       psi = RESHAPE(psi4D,(/PRODUCT(dims)/))
       DEALLOCATE(psi4D)
    ENDIF
    
    ! Close everything opened.
    CALL h5sclose_f(dataspace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(grp_id, error)
    CALL h5fclose_f(file_id, error)
    
    ! Close FORTRAN interface.
    CALL h5close_f(error)
    
  END SUBROUTINE dread_subset_serial
  
  !-------------------------------------------------------------------!
  
  !-------------------------------------------------------------------!
#if _COM_MPI
  SUBROUTINE dread_subset_parallel(filename_in, psi, rank, dims, &
       dims_proc, offset, comm)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)      :: filename_in
    REAL(dp), INTENT(OUT)             :: psi(:)
    INTEGER, INTENT(IN)               :: rank
    INTEGER, INTENT(IN)               :: dims(:)
    INTEGER, INTENT(IN)               :: dims_proc(:)
    INTEGER, INTENT(IN), OPTIONAL     :: offset(:)
    INTEGER, INTENT(IN)               :: comm

    REAL(dp), ALLOCATABLE             :: psi2D(:, :)
    REAL(dp), ALLOCATABLE             :: psi3D(:, :, :)
    REAL(dp), ALLOCATABLE             :: psi4D(:, :, :, :)
    ! File identifier
    INTEGER(HID_T)                   :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                    :: plist_id
    ! Group identifier
    INTEGER(HID_T)                    :: grp_id
    ! Dataset identifier
    INTEGER(HID_T)                    :: dset_id 
    ! Dataspace identifier
    INTEGER(HID_T)                    :: dataspace 
    ! dset size
    INTEGER(HSIZE_T)                  :: subset_size(rank)
    INTEGER(HSIZE_T)                  :: subset_size_proc(rank)
    ! Dataspace identifier in file 
    INTEGER(HID_T)                    :: filespace
    ! Memspace identifier
    INTEGER(HID_T)                    :: memspace
    ! Hyperslab offset
    INTEGER(HSIZE_T)                  :: subset_offset(rank)
    ! Error flag
    INTEGER                           :: error
    CHARACTER(LEN=100)                :: filename
    CHARACTER(LEN=100)                :: dsetname, groupname
    
    !----------------------------------------------!
    
    IF(rank.EQ.2) THEN
       ALLOCATE(psi2D(1:dims_proc(1),1:dims_proc(2)))
    ELSEIF(rank.EQ.3) THEN
       ALLOCATE(psi3D(1:dims_proc(1),1:dims_proc(2),1:dims_proc(3)))
    ELSEIF(rank.EQ.4) THEN
       ALLOCATE(psi4D(1:dims_proc(1),1:dims_proc(2),1:dims_proc(3), &
            1:dims_proc(4)))
    ELSE
       WRITE(*,*) 'No routine for your rank!!'
       STOP
    ENDIF
    
    ! Initialize FORTRAN interface. 
    CALL h5open_f(error) 
    
    ! Setup file access property list with 
    ! parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error) 
    
    ! Open the file.
    filename = TRIM(filename_in) // '.h5'
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, &
         error, access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)
    
    ! Open the group.
    groupname = '/' // TRIM(filename_in)
    CALL h5gopen_f(file_id, groupname, grp_id, error)
    
    ! Set dims
    IF(PRESENT(offset)) THEN
       subset_offset = offset
    ELSE
       subset_offset = 0
    ENDIF
    subset_size_proc = dims_proc
    subset_size = dims
    
    ! Create filespace and memspace
    CALL h5screate_simple_f(rank,subset_size,filespace, error)
    CALL h5screate_simple_f(rank,subset_size_proc,memspace,error)
    
    ! Open the dataset.
    dsetname = '/' // TRIM(filename_in) // '/' // TRIM(filename_in)
    CALL h5dopen_f(file_id, dsetname, dset_id, error)
    
    ! Get dataset's dataspace identifier and select subset.
    CALL h5dget_space_f(dset_id, filespace, error)
    CALL h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
         subset_offset, subset_size_proc, error) 
    
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    ! Read the dataset
    IF(rank.EQ.2) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi2D, subset_size, &
            error, file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
       psi = RESHAPE(psi2D,(/PRODUCT(dims_proc)/))
       DEALLOCATE(psi2D)
    ELSEIF(rank.EQ.3) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi3D, subset_size, &
            error, file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
       psi = RESHAPE(psi3D,(/PRODUCT(dims_proc)/))
       DEALLOCATE(psi3D)
    ELSEIF(rank.EQ.4) THEN
       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi4D, subset_size, &
            error, file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
       psi = RESHAPE(psi4D,(/PRODUCT(dims_proc)/))
       DEALLOCATE(psi4D)
    ENDIF

    ! Close everything opened.
    CALL h5sclose_f(filespace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5pclose_f(plist_id, error)
    CALL h5gclose_f(grp_id, error)
    CALL h5fclose_f(file_id, error)
    
    
    ! Close FORTRAN interface.
    CALL h5close_f(error)
    
  END SUBROUTINE dread_subset_parallel
  
  !-------------------------------------------------------------------!
  
#endif
  
END MODULE io_subset
