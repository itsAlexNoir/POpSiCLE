MODULE io_pop
  
  USE HDF5
  USE constants
  USE tools
  USE coords
#if _COM_MPI
  USE MPI
  USE patchwork
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC               :: write_wave
  PUBLIC               :: open_surface_file
  PUBLIC               :: write_surface_file
  PUBLIC               :: close_surface_file
  PUBLIC               :: read_subset
#if _COM_MPI
  PUBLIC               :: write_slice_xyz
  PUBLIC               :: write_slice_xy
  PUBLIC               :: write_slice_xz
  PUBLIC               :: write_slice_yz
  PUBLIC               :: write_slice_rz
#endif
  
  INTERFACE write_wave
     MODULE PROCEDURE write_wave_serial
#if _COM_MPI
     MODULE PROCEDURE write_wave_parallel
#endif
  END INTERFACE write_wave
  
  INTERFACE open_surface_file
     MODULE PROCEDURE open_surface_file_serial
#if _COM_MPI
     MODULE PROCEDURE open_surface_file_parallel
#endif
  END INTERFACE open_surface_file
  
  INTERFACE write_surface_file
     MODULE PROCEDURE write_surface_file2D_serial
     !MODULE PROCEDURE write_surface_file3D_serial
#if _COM_MPI
     MODULE PROCEDURE write_surface_file2D_parallel
     !MODULE PROCEDURE write_surface_file3D_parallel
#endif
  END INTERFACE write_surface_file
  
  INTERFACE close_surface_file
     MODULE PROCEDURE close_surface_file_serial
#if _COM_MPI
     MODULE PROCEDURE close_surface_file_parallel
#endif
  END INTERFACE close_surface_file
  
  INTERFACE read_subset
     MODULE PROCEDURE read_subset_3Dserial
#if _COM_MPI
     MODULE PROCEDURE read_subset_parallel
#endif
  END INTERFACE read_subset
  
  !--------------------!
  ! Private variables  !
  !--------------------!
  ! File identifier
  INTEGER(HID_T)                 :: file_id
  ! Property list identifier 
  INTEGER(HID_T)                 :: plist_id
  ! Group identifiers
  INTEGER(HID_T)                 :: wave_group
  INTEGER(HID_T)                 :: wavederiv_group 
  INTEGER(HID_T)                 :: field_group
  
  INTEGER                        :: dataset_count
  
  !*************************************************!
  !*************************************************!
  
CONTAINS
  
  SUBROUTINE open_surface_file_serial(filename)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name
    !--------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)
    
    ! Create the file.
    name = TRIM(filename) // '.h5'  
    CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error)
    
    ! Create group in the root group.     
    name = '/wave'     
    CALL h5gcreate_f(file_id, name, wave_group, error)
    name = '/wavederiv'     
    CALL h5gcreate_f(file_id, name, wavederiv_group, error)
    name = '/field'     
    CALL h5gcreate_f(file_id, name, field_group, error)
    
    ! Set the count of datasets to zero. Start 
    dataset_count = 0
    
  END SUBROUTINE open_surface_file_serial
  
  !******************************************************!
#if _COM_MPI  
  SUBROUTINE open_surface_file_parallel(filename, comm)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: comm
    
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name
    !--------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
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
    name = '/wave'     
    CALL h5gcreate_f(file_id, name, wave_group, error)
    name = '/wavederiv'     
    CALL h5gcreate_f(file_id, name, wavederiv_group, error)
    name = '/field'     
    CALL h5gcreate_f(file_id, name, field_group, error)

    ! Set the count of datasets to zero. Start 
    dataset_count = 0
    
  END SUBROUTINE open_surface_file_parallel
#endif
  !***************************************************!
  
  !***************************************************!
  
  SUBROUTINE write_surface_file2D_serial(wave, wavederiv, time, &
       efield, afield, lmax )
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: wave(:)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:)
    REAL(dp), INTENT(IN)           :: time, efield, afield
    INTEGER, INTENT(IN)            :: lmax
    
    INTEGER                        :: wave_rank, field_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :)
    REAL(dp)                       :: field(3)
    INTEGER                        :: numthetapts
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: wavederiv_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: wavederiv_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(2)
    INTEGER(HSIZE_T)               :: field_dim(1)
    ! Error flag
    INTEGER                        :: error
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN = 100)           :: setname
    
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    dataset_count = dataset_count + 1
    WRITE(cstep,'(I6.6)') dataset_count

    ! Number of theta points in a serial run, its
    numthetapts = 2 * lmax + 1
    !Define rank, in this case 2.
    wave_rank = 2
    field_rank = 1
    wave_dims = (/numthetapts, 2/)
    field_dim = 3
    
    ! Write field array
    ALLOCATE(complex_wave(numthetapts,2))
    ALLOCATE(complex_wavederiv(numthetapts,2))
    complex_wave(:,1) = REAL(wave,dp)
    complex_wave(:,2) = AIMAG(wave)
    complex_wavederiv(:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,2) = AIMAG(wavederiv)
    field = (/ time, efield, afield /)
    
    ! Create the data spaces for the  datasets. 
    CALL h5screate_simple_f(wave_rank, wave_dims, wave_dspace_id, error)
    CALL h5screate_simple_f(wave_rank, wave_dims, wavederiv_dspace_id, error)
    CALL h5screate_simple_f(field_rank, field_dim, field_dspace_id, error)
    
    setname = '/wave/'// cstep
    ! Create the dataset.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         wave_dspace_id, wave_dset_id, error)
    CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, &
         wave_dims, error)
    
    setname = '/wavederiv/'// cstep
    ! Create the dataset.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         wavederiv_dspace_id, wavederiv_dset_id, error)
    CALL h5dwrite_f(wavederiv_dset_id, H5T_NATIVE_DOUBLE, complex_wavederiv, &
         wave_dims, error)
    
    setname = '/field/'// cstep
    ! Create the dataset.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         field_dspace_id, field_dset_id, error)
    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, &
         field_dim, error)
    
    ! End access to the dataset and release 
    ! resources used by it.
    CALL h5dclose_f(wave_dset_id, error)
    CALL h5dclose_f(wavederiv_dset_id, error)
    CALL h5dclose_f(field_dset_id, error)
    
    ! Terminate access to the data space.
    CALL h5sclose_f(wave_dspace_id, error)
    CALL h5sclose_f(wavederiv_dspace_id, error)
    CALL h5sclose_f(field_dspace_id, error)
    
    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
  END SUBROUTINE write_surface_file2D_serial
  
  !******************************************************!
#if _COM_MPI
  SUBROUTINE write_surface_file2D_parallel(wave, wavederiv, time, &
       efield, afield, lmax, numthetapts )
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: wave(:)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:)
    REAL(dp), INTENT(IN)           :: time, efield, afield
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: numthetapts
    
    INTEGER                        :: wave_rank, field_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :)
    REAL(dp)                       :: field(3)
    INTEGER                        :: numpivots
    ! Dataspace identifier in file 
    INTEGER(HID_T)                 :: wave_filespace
    INTEGER(HID_T)                 :: wavederiv_filespace
    INTEGER(HID_T)                 :: field_filespace
    ! Dataspace identifier in memory
    INTEGER(HID_T)                 :: wave_memspace
    INTEGER(HID_T)                 :: wavederiv_memspace
    INTEGER(HID_T)                 :: field_memspace
    ! Dataspace identifier
    !INTEGER(HID_T)                 :: wave_dspace_id
    !INTEGER(HID_T)                 :: wavederiv_dspace_id
    !INTEGER(HID_T)                 :: field_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: wavederiv_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_global_dims(2)
    INTEGER(HSIZE_T)               :: wave_dims(2)
    INTEGER(HSIZE_T)               :: field_dim(1)
    ! Memory variables, per proc
    INTEGER(HSIZE_T)               :: cont(2)  
    INTEGER(HSSIZE_T)              :: offset(2)
    INTEGER(HSIZE_T)               :: cont_field(1) 
    INTEGER(HSSIZE_T)              :: offset_field(1)  
    ! Error flag
    INTEGER                        :: error
    
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN = 100)           :: setname
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    dataset_count = dataset_count + 1
    WRITE(cstep,'(I6.6)') dataset_count

    numpivots = 2 * lmax + 1
    !Define rank, in this case 2.
    wave_rank = 2
    field_rank = 1
    wave_dims = (/numthetapts, 2/)
    field_dim = 3
    wave_global_dims = (/numpivots, 2/)
    
    ! Write field array
    ALLOCATE(complex_wave(numthetapts,2))
    ALLOCATE(complex_wavederiv(numthetapts,2))
    complex_wave(:,1) = REAL(wave,dp)
    complex_wave(:,2) = AIMAG(wave)
    complex_wavederiv(:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,2) = AIMAG(wavederiv)
    
    field = (/ time, efield, afield /)
    
    ! Create the data space for the dataset.
    CALL h5screate_simple_f(wave_rank, wave_global_dims, wave_filespace, error)
    CALL h5screate_simple_f(wave_rank, wave_dims, wave_memspace, error)
    CALL h5screate_simple_f(wave_rank, wave_global_dims, wavederiv_filespace, error)
    CALL h5screate_simple_f(wave_rank, wave_dims, wavederiv_memspace, error)
    
    ! Create dataset.
    setname = '/wave/'// cstep
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, wave_filespace, &
         wave_dset_id, error)!, plist_id)
    CALL h5sclose_f(wave_filespace, error)
    setname = '/wavederiv/'// cstep
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, wavederiv_filespace, &
         wavederiv_dset_id, error)!, plist_id)
    CALL h5sclose_f(wavederiv_filespace, error)
    
    cont = wave_dims
    offset = (/ pivot_number(1)-1, 0 /)
    
    ! Select hyperslab in the file.
    CALL h5dget_space_f(wave_dset_id, wave_filespace, error)
    CALL h5sselect_hyperslab_f (wave_filespace, H5S_SELECT_SET_F, offset,&
         cont, error)!, stride, blck)
    CALL h5dget_space_f(wavederiv_dset_id, wavederiv_filespace, error)
    CALL h5sselect_hyperslab_f (wavederiv_filespace, H5S_SELECT_SET_F, offset,&
         cont, error)!, stride, blck)
    
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    ! Write the dataset collectively.
    CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, wave_global_dims, &
         error, file_space_id = wave_filespace, mem_space_id = wave_memspace, &
         xfer_prp = plist_id)
    CALL h5dwrite_f(wavederiv_dset_id, H5T_NATIVE_DOUBLE, complex_wavederiv, wave_global_dims, &
         error, file_space_id = wavederiv_filespace, mem_space_id = wavederiv_memspace, &
         xfer_prp = plist_id)
    
    ! Close dataspaces.
    CALL h5sclose_f(wave_filespace, error)
    CALL h5sclose_f(wave_memspace, error)
    CALL h5sclose_f(wavederiv_filespace, error)
    CALL h5sclose_f(wavederiv_memspace, error)
    
    ! Close the dataset.
    CALL h5dclose_f(wave_dset_id, error)
    CALL h5dclose_f(wavederiv_dset_id, error)
    
    ! For the field, only the first processor in the group
    ! write the field data.

    ! Create the data space for the dataset.
    CALL h5screate_simple_f(field_rank, field_dim, field_filespace, error)
    CALL h5screate_simple_f(field_rank, field_dim, field_memspace, error)
    IF(pivot_number(1).NE.1) &
         CALL h5sselect_none_f(field_memspace,error)
    
    ! Create dataset
    setname = '/field/'// cstep
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, field_filespace, &
         field_dset_id, error)!, plist_id)
    CALL h5sclose_f(field_filespace, error)

    cont_field = field_dim
    offset_field = 0
    ! Select hyperslab in the file.
    CALL h5dget_space_f(field_dset_id, field_filespace, error)
    CALL h5sselect_hyperslab_f (field_filespace, H5S_SELECT_SET_F, offset_field,&
         cont_field, error)!, stride, blck)
    IF(pivot_number(1).NE.1) &
         CALL h5sselect_none_f(field_filespace,error)
       
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    !IF (pivot_number(1) .EQ. 1) &
         CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dim, &
         error, file_space_id = field_filespace, mem_space_id = field_memspace, &
         xfer_prp = plist_id)
    
    ! Close dataspace and dataset
    CALL h5sclose_f(field_filespace, error)
    CALL h5sclose_f(field_memspace, error)
    CALL h5dclose_f(field_dset_id, error)
    
    
!!$    IF (pivot_number(1) .EQ. 1) THEN
!!$       CALL h5screate_simple_f(field_rank, field_dim, field_dspace_id, error)
!!$       setname = '/field/'// cstep
!!$       CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$            field_dspace_id, field_dset_id, error)
!!$       CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, &
!!$            field_dim, error)
!!$       ! End access to the dataset and release 
!!$       ! resources used by it.
!!$       CALL h5dclose_f(field_dset_id, error)
!!$       
!!$       ! Terminate access to the data space.
!!$       CALL h5sclose_f(field_dspace_id, error)
!!$    ENDIF

    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
  END SUBROUTINE write_surface_file2D_parallel
#endif
  
  !***************************************************!
  
  SUBROUTINE close_surface_file_serial( )
    IMPLICIT NONE
    
    ! Error flag
    INTEGER                        :: error, info
    !--------------------------------------------!
    
    ! Close the group.
    CALL h5gclose_f(wave_group, error)
    CALL h5gclose_f(wavederiv_group, error)
    CALL h5gclose_f(field_group, error)
    
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    
    ! Close FORTRAN interfaces and HDF5 library.
    CALL h5close_f(error)
    
  END SUBROUTINE close_surface_file_serial
  
  !***************************************************!
#if _COM_MPI
  SUBROUTINE close_surface_file_parallel( rank )
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: rank
    ! Error flag
    INTEGER                        :: error, info
    !--------------------------------------------!
    
    IF(i_am_surface(rank).EQ.1) THEN
       ! Close the property list.
       CALL h5pclose_f(plist_id, error)
       
       ! Close the group.
       CALL h5gclose_f(wave_group, error)
       CALL h5gclose_f(wavederiv_group, error)
       CALL h5gclose_f(field_group, error)
       
       ! Close the file.
       CALL h5fclose_f(file_id, error)
       
       ! Close FORTRAN interfaces and HDF5 library.
       CALL h5close_f(error)
    ENDIF
    
  END SUBROUTINE close_surface_file_parallel
#endif
  !***************************************************!
  !***************************************************!
  !***************************************************!
 
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
    ! Group identifiers
    INTEGER(HID_T)                 :: group_id 
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
    
    ! Create group in the root group.     
    name = '/' // TRIM(filename)
    !write(*,*) 'name: ',name
    CALL h5gcreate_f(file_id, filename, group_id, error)
    
    ! Create the data spaces for the  datasets. 
    CALL h5screate_simple_f(rank, file_dims, dspace_id, error)
    
    ! Create the dataset.
    setname = '/' // TRIM(filename) // '/'// TRIM(filename)
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
  
  !------------------------------------------------------------!
  
  
  SUBROUTINE write_slice_xyz(density3D, rank, dims_local, dims_global, &
       globalcomm, grid_rank, grid_addresses, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: density3D(:, :, :)
    INTEGER, INTENT(IN)         :: rank
    INTEGER, INTENT(IN)         :: dims_local(:)
    INTEGER, INTENT(IN)         :: dims_global(:)
    INTEGER, INTENT(IN)         :: globalcomm
    INTEGER, INTENT(IN)         :: grid_rank(:)
    INTEGER, INTENT(IN)         :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename 
    
    !--------------------------------------------------!
    
    REAL(dp), ALLOCATABLE        :: densitytotal3D(:, :, :)
    INTEGER                      :: ipgrid(3)
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ierror
    
    !---------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ipgrid  = (/ipx, ipy, ipz/)
    
    ! Allocate global array
    ALLOCATE(densitytotal3D(1:Nx,1:Ny,1:Nz))
    
    ! Sum our slice across processors
    CALL sum_grid_pieces_xyz(density3D, densitytotal3D, rank, &
         dims_local, dims_global, &
         globalcomm, grid_rank, grid_addresses )
    
    ! Get a communicator for those who contribute
    iogroup_size = numprocx * numprocy * numprocz
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    ipro = 0
    DO iprocz = 0, numprocz - 1
       DO iprocy = 0, numprocy - 1
          DO iprocx = 0, numprocx - 1
             ipro = ipro + 1
             members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                  1,numprocx,numprocy,numprocz,numprocr))      
          ENDDO
       ENDDO
    ENDDO
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    IF (ipr.EQ.0) THEN
       CALL write_wave(RESHAPE(densitytotal3D,(/Nx*Ny*Nz/)), 3, &
            (/Nx,Ny,Nz/), (/Nxgl,Nygl,Nzgl/), &
            iocomm, ipgrid, filename )
    ENDIF
    
    DEALLOCATE(densitytotal3D)
    DEALLOCATE(members)
    
  END SUBROUTINE write_slice_xyz
  
  !--------------------------------------------------------!
  
  SUBROUTINE write_slice_xy(density2D, rank, dims_local, dims_global, &
       globalcomm, grid_rank, grid_addresses, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: density2D(:, :)
    INTEGER, INTENT(IN)         :: rank
    INTEGER, INTENT(IN)         :: dims_local(:)
    INTEGER, INTENT(IN)         :: dims_global(:)
    INTEGER, INTENT(IN)         :: globalcomm
    INTEGER, INTENT(IN)         :: grid_rank(:)
    INTEGER, INTENT(IN)         :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename 
    
    !--------------------------------------------------!
    
    REAL(dp), ALLOCATABLE        :: densitytotal2D(:, :)
    INTEGER                      :: ipgrid(2)
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ierror
    
    !---------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ipgrid  = (/ipx, ipy/)
    
    ! Allocate global array
    ALLOCATE(densitytotal2D(1:Nx,1:Ny))
    
    ! Sum our slice across processors
    CALL sum_grid_pieces_xy(density2D, densitytotal2D, rank, &
         dims_local, dims_global, &
         globalcomm, grid_rank, grid_addresses )
    
    ! Get a communicator for those who contribute
    iogroup_size = numprocx * numprocy
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    ipro = 0
    DO iprocy = 0, numprocy - 1
       DO iprocx = 0, numprocx - 1
          ipro = ipro + 1
          members(ipro) = grid_addresses(indx(iprocx+1,iprocy+1,1, &
               1,numprocx,numprocy,numprocz,numprocr))      
       ENDDO
    ENDDO
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    IF (ipz.EQ.0 .AND. ipr.EQ.0) THEN
       CALL write_wave(RESHAPE(densitytotal2D,(/Nx * Ny/)), 2, &
            (/Nx, Ny/), (/Nxgl, Nygl/), &
            iocomm, ipgrid, filename )
    ENDIF
    
    DEALLOCATE(densitytotal2D)
    DEALLOCATE(members)
    
  END SUBROUTINE write_slice_xy
  
  !--------------------------------------------------------!
  
  SUBROUTINE write_slice_xz(density2D, rank, dims_local, dims_global, &
       globalcomm, grid_rank, grid_addresses, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)        :: density2D(:, :)
    INTEGER, INTENT(IN)         :: rank
    INTEGER, INTENT(IN)         :: dims_local(:)
    INTEGER, INTENT(IN)         :: dims_global(:)
    INTEGER, INTENT(IN)         :: globalcomm
    INTEGER, INTENT(IN)         :: grid_rank(:)
    INTEGER, INTENT(IN)         :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    !--------------------------------------------------!
    
    REAL(dp), ALLOCATABLE        :: densitytotal2D(:, :)
    INTEGER                      :: ipgrid(2)
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ierror
    
    !---------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ipgrid  = (/ipx, ipz/)
    
    ! Allocate global array
    ALLOCATE(densitytotal2D(1:Nx,1:Nz))
    
    ! Sum our slice across processors
    CALL sum_grid_pieces_xz(density2D, densitytotal2D, rank, &
         dims_local, dims_global, &
         globalcomm, grid_rank, grid_addresses )
    
    ! Get a communicator for those who contribute
    iogroup_size = numprocx * numprocz
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    ipro = 0
    DO iprocz = 0, numprocz - 1
       DO iprocx = 0, numprocx - 1
          ipro = ipro + 1
          members(ipro) = grid_addresses(indx(iprocx+1,1,iprocz+1, &
               1,numprocx,numprocy,numprocz,numprocr))      
       ENDDO
    ENDDO
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    IF (ipy.EQ.0 .AND. ipr.EQ.0) THEN
       CALL write_wave(RESHAPE(densitytotal2D,(/Nx * Nz/)), 2, &
            (/Nx,Nz/),(/Nxgl,Nzgl/), &          
            iocomm, ipgrid, filename )
    ENDIF
    
    DEALLOCATE(densitytotal2D)
    DEALLOCATE(members)
    
  END SUBROUTINE write_slice_xz

  !--------------------------------------------------------!

   SUBROUTINE write_slice_yz(density2D, rank, dims_local, dims_global, &
       globalcomm, grid_rank, grid_addresses, filename)
    
    IMPLICIT NONE

    REAL(dp), INTENT(IN)        :: density2D(:, :)
    INTEGER, INTENT(IN)         :: rank
    INTEGER, INTENT(IN)         :: dims_local(:)
    INTEGER, INTENT(IN)         :: dims_global(:)
    INTEGER, INTENT(IN)         :: globalcomm
    INTEGER, INTENT(IN)         :: grid_rank(:)
    INTEGER, INTENT(IN)         :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    !--------------------------------------------------!
    
    REAL(dp), ALLOCATABLE        :: densitytotal2D(:, :)
    INTEGER                      :: ipgrid(2)
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ierror
    
    !---------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ipgrid  = (/ipy, ipz/)
    
    ! Allocate global array
    ALLOCATE(densitytotal2D(1:Ny,1:Nz))
    
    ! Sum our slice across processors
    CALL sum_grid_pieces_yz(density2D, densitytotal2D, rank, &
         dims_local, dims_global, &
         globalcomm, grid_rank, grid_addresses )
    
    ! Get a communicator for those who contribute
    iogroup_size = numprocy * numprocz
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    ipro = 0
    DO iprocz = 0, numprocz - 1
       DO iprocy = 0, numprocy - 1
          ipro = ipro + 1
          members(ipro) = grid_addresses(indx(1,iprocy+1,iprocz+1, &
               1,numprocx,numprocy,numprocz,numprocr))      
       ENDDO
    ENDDO
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    IF (ipx.EQ.0 .AND. ipr.EQ.0) THEN
       CALL write_wave(RESHAPE(densitytotal2D,(/Ny * Nz/)), 2, &
            (/Ny, Nz/), (/Nygl, Nzgl/), &
            iocomm, ipgrid, filename )
    ENDIF
    
    DEALLOCATE(densitytotal2D)
    DEALLOCATE(members)
    
  END SUBROUTINE write_slice_yz

  !------------------------------------------------------------------------!

   SUBROUTINE write_slice_rz( density2D, rank, dims_local, dims_global, &
       globalcomm, grid_rank, grid_addresses, filename)
    
    IMPLICIT NONE

    REAL(dp), INTENT(IN)        :: density2D(:, :)
    INTEGER, INTENT(IN)         :: rank
    INTEGER, INTENT(IN)         :: dims_local(:)
    INTEGER, INTENT(IN)         :: dims_global(:)
    INTEGER, INTENT(IN)         :: globalcomm
    INTEGER, INTENT(IN)         :: grid_rank(:)
    INTEGER, INTENT(IN)         :: grid_addresses(:)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    !--------------------------------------------------!
    
    REAL(dp), ALLOCATABLE        :: densitytotal2D(:, :)
    INTEGER                      :: ipgrid(2)
    INTEGER, ALLOCATABLE         :: members(:)
    INTEGER                      :: simgroup, iogroup, iocomm
    INTEGER                      :: iogroup_size, ipro
    INTEGER                      :: numprocx, numprocy, numprocz, numprocr
    INTEGER                      :: Nx, Ny, Nz, Nr
    INTEGER                      :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                      :: ipx, ipy, ipz, ipr
    INTEGER                      :: iprocx, iprocy
    INTEGER                      :: iprocz, iprocr
    INTEGER                      :: ierror
    
    !---------------------------------------------------!
    
    ! Create a global group, used for I/O
    call MPI_COMM_GROUP(globalcomm, simgroup, ierror)
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
    
    ipgrid  = (/ipz, ipr/)
    
    ! Allocate global array
    ALLOCATE(densitytotal2D(1:Nz,1:Nr))
    
    ! Sum our slice across processors
    CALL sum_grid_pieces_rz(density2D, densitytotal2D, rank, &
         dims_local, dims_global, &
         globalcomm, grid_rank, grid_addresses )
    
    ! Get a communicator for those who contribute
    iogroup_size = numprocr * numprocz
    ALLOCATE(members(1:iogroup_size))
    members = 0
    
    ipro = 0
    DO iprocr = 0, numprocr - 1
       DO iprocz = 0, numprocz - 1
          ipro = ipro + 1
          members(ipro) = grid_addresses(indx(1,1,iprocz+1, &
               iprocr+1,numprocx,numprocy,numprocz,numprocr))      
       ENDDO
    ENDDO
    
    CALL MPI_GROUP_INCL(simgroup, iogroup_size, members, iogroup, ierror)
    
    CALL MPI_COMM_CREATE(globalcomm, iogroup, iocomm, ierror)
    
    IF (ipx.EQ.0 .AND. ipy.EQ.0) THEN
       CALL write_wave(RESHAPE(densitytotal2D,(/Nr * Nz/)), 2, &
            (/ Nz, Nr/), (/Nzgl, Nrgl/), &
            iocomm, ipgrid, filename )
    ENDIF
    
    DEALLOCATE(densitytotal2D)
    DEALLOCATE(members)
    
  END SUBROUTINE write_slice_rz
  
  !-------------------------------------------------------------------!

#endif  

  !-------------------------------------------------------------------!

  SUBROUTINE read_subset_3Dserial(filename_in, psi_in, rank, dims, offset)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)      :: filename_in
    REAL(dp), INTENT(IN)              :: psi_in(:, :, :)
    INTEGER, INTENT(IN)               :: rank
    INTEGER, INTENT(IN)               :: dims(:)
    INTEGER, INTENT(IN)               :: offset
    
    ! Dataset identifier
    INTEGER(HID_T)                    :: dset_id 
    ! Dataspace identifier
    INTEGER(HID_T)                    :: dataspace 
    ! Dset size
    INTEGER(HSIZE_T)                  :: subset_size(3)
    ! Hyperslab offset
    INTEGER(HSIZE_T)                  :: subset_offset(3)
    ! Error flag
    INTEGER                           :: error
    CHARACTER(LEN=100)                :: filename, dsetname
    
    !----------------------------------------------!
    
    ! Initialize FORTRAN interface. 
    CALL h5open_f(error) 
    
    ! Open the file.
    filename = TRIM(filename_in) // '.h5'
    CALL h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)
    
    ! Open the  dataset.
    dsetname = '/' // TRIM(filename_in) // '/' // TRIM(filename_in)
    CALL h5dopen_f(file_id, dsetname, dset_id, error)
    
    ! Get dataset's dataspace identifier and select subset.
    subset_ofset = offset
    subset_size = dims
    CALL h5dget_space_f(dset_id, dataspace, error)
    CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
         subset_offset, subset_size, error) 
    
    ! Read the dataset
    CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, psi_in, subset_dims, error)
    
    ! Close everything opened.
    CALL h5sclose_f(dataspace, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5fclose_f(file_id, error)
    
    
    ! Close FORTRAN interface.
    CALL h5close_f(error)
    
  END SUBROUTINE read_subset_serial
  
  !-------------------------------------------------------------------!
  
END MODULE io_pop
