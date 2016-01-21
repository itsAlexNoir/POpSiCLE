MODULE io_surface
  
  USE HDF5
  USE constants
  USE cylboundary
  USE cartboundary
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC               :: open_surface_file
  PUBLIC               :: write_surface_file
  PUBLIC               :: close_surface_file
  PUBLIC               :: get_surface_dims
  PUBLIC               :: read_surface
  
  !-----------------------------!
  
  INTERFACE open_surface_file
     MODULE PROCEDURE open_surface_file_serial
#if _COM_MPI
     MODULE PROCEDURE open_surface_file_parallel
#endif
  END INTERFACE open_surface_file
  
  INTERFACE write_surface_file
     MODULE PROCEDURE write_surface_file2D_serial
     MODULE PROCEDURE write_surface_file3D_serial
#if _COM_MPI
     MODULE PROCEDURE write_surface_file2D_parallel
     MODULE PROCEDURE write_surface_file3D_parallel
#endif
  END INTERFACE write_surface_file
  
  INTERFACE close_surface_file
     MODULE PROCEDURE close_surface_file_serial
#if _COM_MPI
     MODULE PROCEDURE close_surface_file_parallel
#endif
  END INTERFACE close_surface_file
    
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
    WRITE(cstep,'(I6.6)') dataset_count
    dataset_count = dataset_count + 1
    
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
       efield, afield, lmax, surfacerank, numthetaptsperproc )
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: wave(:)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:)
    REAL(dp), INTENT(IN)           :: time, efield, afield
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: surfacerank
    INTEGER, INTENT(IN)            :: numthetaptsperproc
    
    INTEGER                        :: wave_rank, field_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :)
    REAL(dp)                       :: field(3)
    INTEGER                        :: numthetapts
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
    WRITE(cstep,'(I6.6)') dataset_count
    dataset_count = dataset_count + 1
    
    numthetapts = 2 * lmax + 1
    !Define rank, in this case 2.
    wave_rank = 2
    field_rank = 1
    wave_dims = (/numthetaptsperproc, 2/)
    field_dim = 3
    wave_global_dims = (/numthetapts, 2/)
    
    ! Write field array
    ALLOCATE(complex_wave(numthetaptsperproc,2))
    ALLOCATE(complex_wavederiv(numthetaptsperproc,2))
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
    offset = (/ numthetaptsperproc * surfacerank, 0 /)
    
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
    IF(surfacerank.NE.1) &
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
    IF(surfacerank.NE.1) &
         CALL h5sselect_none_f(field_filespace,error)
       
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dim, &
         error, file_space_id = field_filespace, mem_space_id = field_memspace, &
         xfer_prp = plist_id)
    
    ! Close dataspace and dataset
    CALL h5sclose_f(field_filespace, error)
    CALL h5sclose_f(field_memspace, error)
    CALL h5dclose_f(field_dset_id, error)
    
    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
  END SUBROUTINE write_surface_file2D_parallel
#endif
  
  !***************************************************!
  
  !***************************************************!
  
  SUBROUTINE write_surface_file3D_serial(wave, wavederiv, time, &
       efield, afield, lmax )
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
    REAL(dp), INTENT(IN)           :: time, efield, afield
    INTEGER, INTENT(IN)            :: lmax
    
    INTEGER                        :: wave_rank, field_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
    REAL(dp)                       :: field(3)
    INTEGER                        :: numthetapts, numphipts
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: wavederiv_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: wavederiv_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(3)
    INTEGER(HSIZE_T)               :: field_dim(1)
    ! Error flag
    INTEGER                        :: error
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN = 100)           :: setname
    
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    WRITE(cstep,'(I6.6)') dataset_count
    dataset_count = dataset_count + 1
    
    ! Number of theta points in a serial run, its
    numthetapts = 2 * lmax + 1
    numphipts = 2 * lmax + 1
    !Define rank, in this case 2.
    wave_rank = 3
    field_rank = 1
    wave_dims = (/numthetapts, numphipts, 2/)
    field_dim = 3
    
    ! Write field array
    ALLOCATE(complex_wave(numthetapts,numphipts,2))
    ALLOCATE(complex_wavederiv(numthetapts,numphipts,2))
    complex_wave(:,:,1) = REAL(wave,dp)
    complex_wave(:,:,2) = AIMAG(wave)
    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
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
    
  END SUBROUTINE write_surface_file3D_serial
  
  !******************************************************!
#if _COM_MPI
  SUBROUTINE write_surface_file3D_parallel(wave, wavederiv, time, &
       efield, afield, lmax, surfacerank, &
       numthetaptsperproc, numphiptsperproc )
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
    REAL(dp), INTENT(IN)           :: time, efield, afield
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: surfacerank
    INTEGER, INTENT(IN)            :: numthetaptsperproc
    INTEGER, INTENT(IN)            :: numphiptsperproc
    
    INTEGER                        :: wave_rank, field_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
    REAL(dp)                       :: field(3)
    INTEGER                        :: numthetapts, numphipts
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
    INTEGER(HSIZE_T)               :: wave_global_dims(3)
    INTEGER(HSIZE_T)               :: wave_dims(3)
    INTEGER(HSIZE_T)               :: field_dim(1)
    ! Memory variables, per proc
    INTEGER(HSIZE_T)               :: cont(3)  
    INTEGER(HSSIZE_T)              :: offset(3)
    INTEGER(HSIZE_T)               :: cont_field(1) 
    INTEGER(HSSIZE_T)              :: offset_field(1)  
    ! Error flag
    INTEGER                        :: error
    
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN = 100)           :: setname
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    WRITE(cstep,'(I6.6)') dataset_count
    dataset_count = dataset_count + 1
    
    numthetapts = 2 * lmax + 1
    numphipts = 2 * lmax + 1
    wave_rank = 3
    field_rank = 1
    field_dim = 3
    wave_global_dims = (/numthetapts, numphipts, 2/)
    
    wave_dims = (/numthetaptsperproc, numphiptsperproc, 2/)
    ALLOCATE(complex_wave(numthetaptsperproc,numphiptsperproc,2))
    ALLOCATE(complex_wavederiv(numthetaptsperproc,numphiptsperproc,2))
    complex_wave(:,:,1) = REAL(wave,dp)
    complex_wave(:,:,2) = AIMAG(wave)
    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
    
    ! Write field array   
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
    offset = (/ numthetaptsperproc * surfacerank, numphiptsperproc * surfacerank, 0 /)
    
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
    
    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
    !************!
    
    ! For the field, only the first processor in the group
    ! write the field data.
    
    ! Create the data space for the dataset.
    CALL h5screate_simple_f(field_rank, field_dim, field_filespace, error)
    CALL h5screate_simple_f(field_rank, field_dim, field_memspace, error)
    IF(surfacerank.NE.1) &
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
    
    IF(surfacerank.NE.1) &
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
    
  END SUBROUTINE write_surface_file3D_parallel
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
  
  SUBROUTINE get_surface_dims(filename, ntime, &
       numthetapts, numphipts, lmax)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)      :: filename
    INTEGER, INTENT(OUT)              :: ntime
    INTEGER, INTENT(OUT)              :: numthetapts, numphipts
    INTEGER, INTENT(OUT)              :: lmax
    
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Group identifiers
    INTEGER(HID_T)                 :: group_id 
    ! Dataspace identifier
    INTEGER(HID_T)                 :: dspace_id    
    ! Dataset identifier  
    INTEGER(HID_T)                 :: dset_id      
    ! Dataset dimensions  
    INTEGER(HSIZE_T)               :: dspace_dims(3)
    INTEGER(HSIZE_T)               :: dspace_maxdims(3)
    ! Structure info
    INTEGER                        :: storage_type
    INTEGER                        :: nlinks
    INTEGER                        :: maxcorder
    ! Error flag
    INTEGER                        :: error     
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name

    !----------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)    
    
    ! Open file
    name = TRIM(filename) // '.h5'
    CALL h5fopen_f(name,H5F_ACC_RDONLY_F, file_id, error) 
    
    ! Open wave group, to see how many times are on it.
    name = '/wave'
    CALL h5gopen_f(file_id, name, group_id, error)
    ! Get file info
    CALL h5gget_info_f(group_id, storage_type, nlinks, &
         maxcorder, error)

    ! Assign number of time steps
    ntime = nlinks
    
    name = '/wave/000000'
    CALL h5dopen_f(file_id, name, dset_id, error)

    CALL h5dget_space_f(dset_id, dspace_id, error)
    
    CALL h5sget_simple_extent_dims_f(dspace_id, dspace_dims, &
         dspace_maxdims, error)

    IF(dspace_dims(3).EQ.0) THEN
       numthetapts = dspace_dims(1)
       numphipts = 1
       lmax = (numthetapts - 1)/2
    ELSE
       numthetapts = dspace_dims(1)
       numphipts = dspace_dims(2)
       lmax = (numthetapts - 1)/2
    ENDIF
    
    WRITE(*,*) 'Path and name of the file: ',TRIM(filename) // '.h5'
    WRITE(*,*) 'Number of time step in the file: ',ntime
    WRITE(*,*) 'Number of theta points: ',numthetapts
    WRITE(*,*) 'Number of phi points: ',numphipts
    WRITE(*,*) 'Maximun angular momenta: ',lmax
    
    ! Close dataset
    CALL h5dclose_f(dset_id, error)
    ! Close dataspace
    CALL h5sclose_f(dspace_id, error)
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    
    ! Close HDF5 fortran interface
    CALL h5close_f(error)
    
  END SUBROUTINE get_surface_dims


  !***************************************************!
  
  SUBROUTINE read_surface(filename, itime, numthetapts, &
       numphipts, psi_sph, psip_sph, time, efield, afield)

    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)       :: filename
    INTEGER, INTENT(IN)                :: itime
    INTEGER, INTENT(IN)                :: numthetapts
    INTEGER, INTENT(IN)                :: numphipts
    COMPLEX(dp), INTENT(OUT)           :: psi_sph(:, :)
    COMPLEX(dp), INTENT(OUT)           :: psip_sph(:, :)
    REAL(dp), INTENT(OUT)              :: time
    REAL(dp), INTENT(OUT)              :: efield, afield
    
    REAL(dp), ALLOCATABLE         :: psi2D(:, :)
    REAL(dp), ALLOCATABLE         :: psip2D(:, :)
    REAL(dp), ALLOCATABLE         :: psi3D(:, :, :)
    REAL(dp), ALLOCATABLE         :: psip3D(:, :, :)
    REAL(dp)                      :: field(3)
    INTEGER                       :: itheta, iphi
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group_id 
    INTEGER(HID_T)                 :: wavep_group_id 
    INTEGER(HID_T)                 :: field_group_id 
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id    
    INTEGER(HID_T)                 :: wavep_dspace_id    
    INTEGER(HID_T)                 :: field_dspace_id    
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id      
    INTEGER(HID_T)                 :: wavep_dset_id      
    INTEGER(HID_T)                 :: field_dset_id      
    ! Dataset dimensions  
    INTEGER(HSIZE_T), ALLOCATABLE  :: dspace_dims(:)
    INTEGER(HSIZE_T), ALLOCATABLE  :: dspace_offset(:)
    INTEGER(HSIZE_T)               :: field_dim(1)
    INTEGER(HSIZE_T)               :: field_offset(1)
    ! Memspace identifier
    INTEGER(HID_T)                 :: wave_memspace
    INTEGER(HID_T)                 :: wavep_memspace
    INTEGER(HID_T)                 :: field_memspace
    INTEGER                        :: error     
    ! Auxiliary strings
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN = 100)           :: name

    !---------------------------------------------------!
    
    ! Set time step as a string
    WRITE(cstep,'(I6.6)') itime
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)    
    
    ! Open file
    name = TRIM(filename) // '.h5'
    CALL h5fopen_f(name,H5F_ACC_RDONLY_F, file_id, error) 
    
    ! Open wave group. 
    name = '/wave'
    CALL h5gopen_f(file_id, name, wave_group_id, error)
    name = '/wavederiv'
    CALL h5gopen_f(file_id, name, wavep_group_id, error)
    
    ! Open dataset
    name = '/wave/' // cstep
    CALL h5dopen_f(file_id, name, wave_dset_id, error)
    ! Get dataspace identifier
    CALL h5dget_space_f(wave_dset_id, wave_dspace_id, error)
    name = '/wavederiv/' // cstep
    CALL h5dopen_f(file_id, name, wavep_dset_id, error)
    ! Get dataspace identifier
    CALL h5dget_space_f(wavep_dset_id, wavep_dspace_id, error)
    
    IF(numphipts.EQ.1) THEN
       ALLOCATE(dspace_dims(2),dspace_offset(2))
       dspace_dims = (/ numthetapts, 2/)
       ALLOCATE(psi2D(1:numthetapts,1:2))
       ALLOCATE(psip2D(1:numthetapts,1:2))
    ELSE
       ALLOCATE(dspace_dims(3),dspace_offset(3))
       dspace_dims = (/ numthetapts, numphipts, 2/)
       ALLOCATE(psi3D(1:numthetapts,1:numphipts,1:2))
       ALLOCATE(psip3D(1:numthetapts,1:numphipts,1:2))
    ENDIF
    dspace_offset = 0
    
    CALL h5sselect_hyperslab_f(wave_dspace_id, H5S_SELECT_SET_F, &
         dspace_offset, dspace_dims, error)
    CALL h5sselect_hyperslab_f(wavep_dspace_id, H5S_SELECT_SET_F, &
         dspace_offset, dspace_dims, error)
    
    ! Create memspace and read dataset
    IF(numphipts.EQ.1) THEN
       CALL h5screate_simple_f(2,dspace_dims,wave_memspace,error)
       CALL h5screate_simple_f(2,dspace_dims,wavep_memspace,error)
       
       CALL h5dread_f(wave_dset_id, H5T_NATIVE_DOUBLE, psi2D, dspace_dims, error, &
            wave_memspace, wave_dspace_id )
       CALL h5dread_f(wavep_dset_id, H5T_NATIVE_DOUBLE, psip2D, dspace_dims, error, &
            wavep_memspace, wavep_dspace_id )
       
       DO itheta = 1, numthetapts
          psi_sph(itheta,1) = CMPLX(psi2D(itheta,1), psi2D(itheta,2))
          psip_sph(itheta,1) = CMPLX(psip2D(itheta,1), psip2D(itheta,2))
       ENDDO
       
    ELSE
       CALL h5screate_simple_f(3,dspace_dims,wave_memspace,error)
       CALL h5screate_simple_f(3,dspace_dims,wavep_memspace,error)
       
       CALL h5dread_f(wave_dset_id, H5T_NATIVE_DOUBLE, psi2D, dspace_dims, error, &
            wave_memspace, wave_dspace_id )
       CALL h5dread_f(wavep_dset_id, H5T_NATIVE_DOUBLE, psip2D, dspace_dims, error, &
            wavep_memspace, wavep_dspace_id )
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             psi_sph(itheta,iphi) = CMPLX(psi3D(itheta,iphi,1), psi3D(itheta,iphi,2))
             psip_sph(itheta,iphi) = CMPLX(psip3D(itheta,iphi,1), psip3D(itheta,iphi,2))
          ENDDO
       ENDDO
    ENDIF
    
    ! Close (almost) everything
    CALL h5sclose_f(wave_dspace_id, error)
    CALL h5sclose_f(wave_memspace, error)
    CALL h5dclose_f(wave_dset_id, error)
    CALL h5gclose_f(wave_group_id, error)
    
    CALL h5sclose_f(wavep_dspace_id, error)
    CALL h5sclose_f(wavep_memspace, error)
    CALL h5dclose_f(wavep_dset_id, error)
    CALL h5gclose_f(wavep_group_id, error)
    
    DEALLOCATE(dspace_dims, dspace_offset)
    IF(numphipts.EQ.1) THEN
       DEALLOCATE(psi2D, psip2D)
    ELSE
       DEALLOCATE(psi3D, psip3D)
    ENDIF
    
    !
    ! Now, we read the field
    !
    
    ! Open wave group. 
    name = '/field'
    CALL h5gopen_f(file_id, name, field_group_id, error)
    
    ! Open dataset
    name = '/field/' // cstep
    CALL h5dopen_f(file_id, name, field_dset_id, error)
    ! Get dataspace identifier
    CALL h5dget_space_f(field_dset_id, field_dspace_id, error)

    field_dim = 3
    field_offset = 0
    field = 0
    
    CALL h5sselect_hyperslab_f(field_dspace_id, H5S_SELECT_SET_F, &
         field_offset, field_dim, error)
    
    CALL h5screate_simple_f(1,field_dim,field_memspace,error)
    
    CALL h5dread_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dim, error, &
         field_memspace, field_dspace_id )
    
    time = field(1)
    efield = field(2)
    afield = field(3)
    
    CALL h5sclose_f(field_dspace_id, error)
    CALL h5sclose_f(field_memspace, error)
    CALL h5dclose_f(field_dset_id, error)
    CALL h5gclose_f(field_group_id, error)
    
    ! Close the file.
    CALL h5fclose_f(file_id, error)
    
    ! Close HDF5 fortran interface
    CALL h5close_f(error)
    
  END SUBROUTINE read_surface
  
  !***************************************************!
  !***************************************************!
  !***************************************************!
  
END MODULE io_surface
  
