MODULE io
  
  USE HDF5
  USE constants
#if _COM_MPI
  USE MPI
#endif

  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC               :: write_wave
  
  INTERFACE write_wave
     MODULE PROCEDURE write_wave_serial
#if _COM_MPI
     MODULE PROCEDURE write_wave_parallel
#endif
  END INTERFACE write_wave
  
  
CONTAINS
  
  SUBROUTINE write_wave_serial(wave,rank,dims,filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: wave(:)
    INTEGER, INTENT(IN)            :: rank
    INTEGER, INTENT(IN)            :: dims(:)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    
    REAL(dp), ALLOCATABLE          :: wave2D(:, :)
    REAL(dp), ALLOCATABLE          :: wave3D(:, :, :)
    REAL(dp), ALLOCATABLE          :: wave4D(:, :, :, :)
    
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Dataspace identifier
    INTEGER(HID_T)                 :: dspace_id    
    ! Dataset identifier  
    INTEGER(HID_T)                 :: dset_id      
    ! Dataset dimensions  
    INTEGER(HSIZE_T), ALLOCATABLE  :: file_dims(:)
    ! Error flag
    INTEGER                        :: error     
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name
    CHARACTER(LEN = 100)           :: setname
    
    !------------------------------------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)    
    
    ALLOCATE(file_dims(1:rank))
    
    file_dims = dims
    
    IF(rank.EQ.2) THEN
       ALLOCATE(wave2D(1:dims(1),1:dims(2)))
       wave2D = RESHAPE(wave,SHAPE(wave2D))
    ELSEIF(rank.EQ.3) THEN
       ALLOCATE(wave3D(1:dims(1),1:dims(2),1:dims(3)))
       wave3D = RESHAPE(wave,SHAPE(wave3D))
    ELSEIF(rank.EQ.4) THEN
       ALLOCATE(wave4D(1:dims(1),1:dims(2),1:dims(3), &
            1:dims(4)))
       wave4D = RESHAPE(wave,SHAPE(wave4D))
    ELSE
       WRITE(*,*) 'No routine for your rank!!'
       STOP
    ENDIF
    
    name = TRIM(filename) // '.h5'  
    CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error)
    
    
    ! Create the data spaces for the  datasets. 
    CALL h5screate_simple_f(rank, file_dims, dspace_id, error)
    
    ! Create the dataset.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         dspace_id, dset_id, error)
    
    IF (rank.EQ.2) THEN
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave2D, &
            file_dims, error)
    ELSEIF(rank.EQ.3) THEN
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave3D, &
            file_dims, error)
    ELSEIF(rank.EQ.4) THEN
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave4D, &
            file_dims, error)
    ELSE
       WRITE(*,*) 'No routine for your rank!!'
       STOP
    ENDIF
    
    ! End access to the dataset and release 
    ! resources used by it.
    CALL h5dclose_f(dset_id, error)
    
    ! Terminate access to the data space.
    CALL h5sclose_f(dspace_id, error)
    
    DEALLOCATE(file_dims)
    
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    
    ! Close HDF5 fortran interface
    CALL h5close_f(error)
    
  END SUBROUTINE write_wave_serial
 
  !-----------------------------------------------------!

#if _COM_MPI

  SUBROUTINE write_wave_parallel(wave, rank, dims_local, &
       dims_global, comm, ipgrid, filename )

    IMPLICIT NONE

    REAL(dp), INTENT(IN)         :: wave(:)
    INTEGER, INTENT(IN)          :: rank
    INTEGER, INTENT(IN)          :: dims_local(:)
    INTEGER, INTENT(IN)          :: dims_global(:)
    INTEGER, INTENT(IN)          :: comm
    INTEGER, INTENT(IN)          :: ipgrid(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    !---------------------------------------------!

    !! Shaped arrrays
    REAL(dp), ALLOCATABLE          :: wave2D(:, :)
    REAL(dp), ALLOCATABLE          :: wave3D(:, :, :)
    REAL(dp), ALLOCATABLE          :: wave4D(:, :, :, :)   
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
    INTEGER(HSIZE_T), ALLOCATABLE  :: stride(:)
    INTEGER(HSIZE_T), ALLOCATABLE  :: blck(:)
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name
    CHARACTER(LEN = 100)           :: setname
    
    !---------------------------------------------!
    
    ALLOCATE(file_dims(1:rank))
    ALLOCATE(proc_dims(1:rank))
    ALLOCATE(stride(1:rank))
    ALLOCATE(cont(1:rank))
    ALLOCATE(blck(1:rank))
    ALLOCATE(offset(1:rank))
    
    file_dims = dims_global
    proc_dims = dims_local

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

    ! Initialize fortran interface
    CALL h5open_f(error)
   
    ! Setup file access property list with 
    ! parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error) 
   
    ! Create the file collectively.
    name = TRIM(filename) // '.h5'  
    CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)

    ! Create group in the root group.     
    name = '/' // TRIM(filename)
    !write(*,*) 'name: ',name
    CALL h5gcreate_f(file_id, filename, group_id, error)
    
    
    ! Create the data space for the  dataset.
    CALL h5screate_simple_f(rank, file_dims, filespace, error)
    CALL h5screate_simple_f(rank, proc_dims, memspace, error)

    ! Create chunked dataset.
    !CALL h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    !CALL h5pset_chunk_f(plist_id, rank, proc_dims, error)

    setname = '/' // TRIM(filename) // '/'// TRIM(filename)
    !write(*,*) 'setname: ',setname
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, filespace, &
         dset_id, error)!, plist_id)
    CALL h5sclose_f(filespace, error)
    
    ! Each process defines dataset in memory and writes it to the hyperslab
    ! in the file. 
    stride = 1
    cont  = proc_dims
    !blck = proc_dims(1)
    IF(rank.EQ.2) THEN
       offset(1) = ipgrid(1) * proc_dims(1)
       offset(2) = ipgrid(2) * proc_dims(2)
    ELSEIF (rank.EQ.3) THEN
       offset(1) = ipgrid(1) * proc_dims(1)
       offset(2) = ipgrid(2) * proc_dims(2)
       offset(3) = ipgrid(3) * proc_dims(3)
    ELSEIF (rank.EQ.4) THEN
       offset(1) = ipgrid(1) * proc_dims(1)
       offset(2) = ipgrid(2) * proc_dims(2)
       offset(3) = ipgrid(3) * proc_dims(3)
       offset(4) = ipgrid(4) * proc_dims(4)
    ELSE 
       WRITE(*,*) 'Number of dimensions not implemented!'
       STOP
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
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave2D, file_dims, &
            error, file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
    ELSEIF(rank.EQ.3) THEN
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave3D, file_dims, &
            error, file_space_id = filespace, mem_space_id = memspace, &
            xfer_prp = plist_id)
    ELSEIF(rank.EQ.4) THEN
       CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, wave4D, file_dims, &
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
    
    ! Free arrays
    DEALLOCATE(stride,cont,blck,offset)
    IF(rank.EQ.2) THEN
       DEALLOCATE(wave2D)
    ELSEIF(rank.EQ.3) THEN
       DEALLOCATE(wave3D)
    ELSEIF(rank.EQ.4) THEN
       DEALLOCATE(wave4D)
    ELSE
       WRITE(*,*) 'No routine for those dims!!!'
       STOP
    ENDIF
    
  END SUBROUTINE write_wave_parallel
#endif  

END MODULE io
