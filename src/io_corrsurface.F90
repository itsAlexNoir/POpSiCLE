MODULE io_corrsurface
  
  USE HDF5
  USE constants
  USE boundary
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC               :: write_correlated_surface_file
  !PUBLIC               :: get_correlated_surface_dims
  !PUBLIC               :: read_correlated_surface
  
  !-----------------------------!
  
  INTERFACE write_correlated_surface_file
     MODULE PROCEDURE write_correlated_surface_file2D_serial
     !MODULE PROCEDURE write_correlated_surface_file3D_serial
#if _COM_MPI
     MODULE PROCEDURE write_correlated_surface_file2D_parallel
     !MODULE PROCEDURE write_correlated_surface_file3D_parallel
#endif
  END INTERFACE write_correlated_surface_file
  
  !--------------------!
  ! Private variables  !
  !--------------------!
  
  INTEGER                        :: dataset_count
  
  !*************************************************!
  !*************************************************!
  
CONTAINS
  
  !***************************************************!
  
  !***************************************************!
  
  SUBROUTINE write_correlated_surface_file2D_serial(filename, wave, &
       wavederiv, time, efield, afield, lmax, basis_order )
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
    REAL(dp), INTENT(IN)           :: time
    REAL(dp), INTENT(IN)           :: efield(:), afield(:)
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: basis_order
    
    INTEGER                        :: wave_rank, field_rank
    INTEGER                        :: time_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
    REAL(dp)                       :: field(3,2)
    INTEGER                        :: numthetapts
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                 :: plist_id
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group
    INTEGER(HID_T)                 :: wavederiv_group 
    INTEGER(HID_T)                 :: field_group
    INTEGER(HID_T)                 :: time_group
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: wavederiv_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id  
    INTEGER(HID_T)                 :: time_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: wavederiv_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    INTEGER(HID_T)                 :: time_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(3)
    INTEGER(HSIZE_T)               :: field_dims(2)
    INTEGER(HSIZE_T)               :: time_dim(1)  
    ! Error flag
    INTEGER                        :: error
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN=100)             :: name
    CHARACTER(LEN = 100)           :: setname
    
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    dataset_count = dataset_count + 1
    WRITE(cstep,'(I6.6)') dataset_count
    
    ! Initialize fortran interface
    CALL h5open_f(error)
    
    ! Open files and groups
    name = TRIM(filename) // '.h5'
    CALL h5fopen_f(name, H5F_ACC_RDWR_F, file_id, error)
    
    name = '/wave'
    CALL h5gopen_f(file_id, name, wave_group, error)
    name = '/wavederiv'
    CALL h5gopen_f(file_id, name, wavederiv_group, error)
    name = '/field'
    CALL h5gopen_f(file_id, name, field_group, error)
    name = '/time'
    CALL h5gopen_f(file_id, name, time_group, error)
    
    ! Number of theta points in a serial run, its
    numthetapts = lmax + 1
    !Define rank, in this case 2.
    wave_rank = 3
    field_rank = 2
    time_rank = 1
    wave_dims = (/basis_order,numthetapts, 2/)
    field_dims = (/ 3,2/)
    time_dim = 1
    
    ! Write field array
    ALLOCATE(complex_wave(basis_order,numthetapts,2))
    ALLOCATE(complex_wavederiv(basis_order,numthetapts,2))
    complex_wave(:,:,1) = REAL(wave,dp)
    complex_wave(:,:,2) = AIMAG(wave)
    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
    field(:,1) = (/ efield(1), efield(2), efield(3) /)
    field(:,2) = (/ afield(1), afield(2), afield(3) /)
    
    ! Create the data spaces for the  datasets. 
    CALL h5screate_simple_f(wave_rank, wave_dims, wave_dspace_id, error)
    CALL h5screate_simple_f(wave_rank, wave_dims, wavederiv_dspace_id, error)
    CALL h5screate_simple_f(field_rank, field_dims, field_dspace_id, error)
    CALL h5screate_simple_f(time_rank, time_dim, time_dspace_id, error)
    
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
         field_dims, error)
    
    setname = '/time/'// cstep
    ! Create the dataset.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         time_dspace_id, time_dset_id, error)
    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, &
         time_dim, error)

    ! End access to the dataset and release 
    ! resources used by it.
    CALL h5dclose_f(wave_dset_id, error)
    CALL h5dclose_f(wavederiv_dset_id, error)
    CALL h5dclose_f(field_dset_id, error)
    CALL h5dclose_f(time_dset_id, error)
    
    ! Terminate access to the data space.
    CALL h5sclose_f(wave_dspace_id, error)
    CALL h5sclose_f(wavederiv_dspace_id, error)
    CALL h5sclose_f(field_dspace_id, error)
    CALL h5sclose_f(time_dspace_id, error)
    
    ! Close the property list, group, file and interface
    CALL h5gclose_f(wave_group, error)
    CALL h5gclose_f(wavederiv_group, error)
    CALL h5gclose_f(field_group, error)
    CALL h5gclose_f(time_group, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
  END SUBROUTINE write_correlated_surface_file2D_serial
  
  !******************************************************!
#if _COM_MPI
  SUBROUTINE write_correlated_surface_file2D_parallel(filename, wave, &
       wavederiv, time, efield, afield, lmax, basis_order, &
       surfacerank, comm )
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
    REAL(dp), INTENT(IN)           :: time
    REAL(dp), INTENT(IN)           :: efield(:), afield(:)
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: basis_order
    INTEGER, INTENT(IN)            :: surfacerank
    INTEGER, INTENT(IN)            :: comm
    
    INTEGER                        :: wave_rank, field_rank
    INTEGER                        :: time_rank
    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
    REAL(dp)                       :: field(3,2)
    INTEGER                        :: numthetaoffset
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                 :: plist_id
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group
    INTEGER(HID_T)                 :: wavederiv_group 
    INTEGER(HID_T)                 :: field_group
    INTEGER(HID_T)                 :: time_group   
    ! Dataspace identifier in file 
    INTEGER(HID_T)                 :: wave_filespace
    INTEGER(HID_T)                 :: wavederiv_filespace
    INTEGER(HID_T)                 :: field_filespace
    INTEGER(HID_T)                 :: time_filespace
    ! Dataspace identifier in memory
    INTEGER(HID_T)                 :: wave_memspace
    INTEGER(HID_T)                 :: wavederiv_memspace
    INTEGER(HID_T)                 :: field_memspace
    INTEGER(HID_T)                 :: time_memspace
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: wavederiv_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id  
    INTEGER(HID_T)                 :: time_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: wavederiv_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    INTEGER(HID_T)                 :: time_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_global_dims(3)
    INTEGER(HSIZE_T)               :: wave_dims(3)
    INTEGER(HSIZE_T)               :: field_dims(2)
    INTEGER(HSIZE_T)               :: time_dim(1)
    ! Memory variables, per proc
    INTEGER(HSIZE_T)               :: cont(3)  
    INTEGER(HSSIZE_T)              :: offset(3)
    INTEGER(HSIZE_T)               :: cont_field(2)
    INTEGER(HSSIZE_T)              :: offset_field(2)  
    INTEGER(HSIZE_T)               :: cont_time(1) 
    INTEGER(HSSIZE_T)              :: offset_time(1)  
    ! Error flag
    INTEGER                        :: error
    
    CHARACTER(LEN=6)               :: cstep
    CHARACTER(LEN=100)             :: name
    CHARACTER(LEN = 100)           :: setname
    
    !------------------------------------------------------!
    
    ! Increase in 1 the count
    dataset_count = dataset_count + 1
    WRITE(cstep,'(I6.6)') dataset_count
    
    ! Initialize fortran interface
    CALL h5open_f(error)
    
    ! Setup file access property list with
    ! parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
    
    ! Open files and groups
    name = TRIM(filename) // '.h5'
    CALL h5fopen_f(name, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)
    
    name = '/wave'
    CALL h5gopen_f(file_id, name, wave_group, error)
    name = '/wavederiv'
    CALL h5gopen_f(file_id, name, wavederiv_group, error)
    name = '/field'
    CALL h5gopen_f(file_id, name, field_group, error)
    name = '/time'
    CALL h5gopen_f(file_id, name, time_group, error)
    
    numthetapts = lmax + 1
    numthetaoffset = INT(numthetapts / numsurfaceprocs)
    !Define rank, in this case 2.
    wave_rank = 3
    field_rank = 2
    time_rank = 1
    wave_dims = (/basis_order,numthetaptsperproc, 2/)
    field_dims = (/ 3, 2/)
    wave_global_dims = (/basis_order,numthetapts, 2/)
    time_dim = 1
    
    ! Write field array
    ALLOCATE(complex_wave(basis_order,numthetaptsperproc,2))
    ALLOCATE(complex_wavederiv(basis_order,numthetaptsperproc,2))
    complex_wave(:,:,1)      = REAL(wave,dp)
    complex_wave(:,:,2)      = AIMAG(wave)
    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
    
    field(:, 1) = (/efield(1), efield(2), efield(3)/)
    field(:, 2) = (/afield(1), afield(2), afield(3)/)
    
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
    offset = (/ 0, numthetaoffset * surfacerank, 0 /)

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
    CALL h5pclose_f(plist_id, error)
    CALL h5sclose_f(wave_filespace, error)
    CALL h5sclose_f(wave_memspace, error)
    CALL h5sclose_f(wavederiv_filespace, error)
    CALL h5sclose_f(wavederiv_memspace, error)
    
    ! Close the dataset.
    CALL h5dclose_f(wave_dset_id, error)
    CALL h5dclose_f(wavederiv_dset_id, error)
    
    ! Deallocate aux arrays for data
    DEALLOCATE(complex_wave,complex_wavederiv)
    
    ! For the field, only the first processor in the group
    ! write the field data.

    ! Create the data space for the dataset.
    CALL h5screate_simple_f(field_rank, field_dims, field_filespace, error)
    CALL h5screate_simple_f(field_rank, field_dims, field_memspace, error)
    IF(surfacerank.NE.1) &
         CALL h5sselect_none_f(field_memspace,error)
    
    ! Create dataset
    setname = '/field/'// cstep
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, field_filespace, &
         field_dset_id, error)!, plist_id)
    CALL h5sclose_f(field_filespace, error)

    cont_field = field_dims
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
    
    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dims, &
         error, file_space_id = field_filespace, mem_space_id = field_memspace, &
         xfer_prp = plist_id)

    CALL h5pclose_f(plist_id, error)
    CALL h5sclose_f(field_filespace, error)
    CALL h5sclose_f(field_memspace, error)
    CALL h5dclose_f(field_dset_id, error)

    !-------------------------------------!

     ! Create the data space for the dataset.
    CALL h5screate_simple_f(time_rank, time_dim, time_filespace, error)
    CALL h5screate_simple_f(time_rank, time_dim, time_memspace, error)
    IF(surfacerank.NE.1) &
         CALL h5sselect_none_f(time_memspace,error)
    
    ! Create dataset
    setname = '/time/'// cstep
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, time_filespace, &
         time_dset_id, error)!, plist_id)
    CALL h5sclose_f(time_filespace, error)

    cont_time = time_dim
    offset_time = 0
    ! Select hyperslab in the file.
    CALL h5dget_space_f(time_dset_id, time_filespace, error)
    CALL h5sselect_hyperslab_f (time_filespace, H5S_SELECT_SET_F, offset_time,&
         cont_time, error)!, stride, blck)
    IF(surfacerank.NE.1) &
         CALL h5sselect_none_f(time_filespace,error)
       
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_dim, &
         error, file_space_id = time_filespace, mem_space_id = time_memspace, &
         xfer_prp = plist_id)
    !---------------!
    
    ! Close dataspace and dataset 
    CALL h5pclose_f(plist_id, error)
    CALL h5sclose_f(time_filespace, error)
    CALL h5sclose_f(time_memspace, error)
    CALL h5dclose_f(time_dset_id, error)

     ! Close the property list, group, file and interface
    CALL h5gclose_f(wave_group, error)
    CALL h5gclose_f(wavederiv_group, error)
    CALL h5gclose_f(field_group, error)
    CALL h5gclose_f(time_group, error)
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
    
  END SUBROUTINE write_correlated_surface_file2D_parallel
#endif
  
  !***************************************************!
  
  !***************************************************!
!!$  
!!$  SUBROUTINE write_surface_file3D_serial(filename, wave, &
!!$       wavederiv, time, efield, afield, lmax )
!!$    IMPLICIT NONE
!!$
!!$    CHARACTER(LEN=*), INTENT(IN)   :: filename
!!$    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
!!$    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
!!$    REAL(dp), INTENT(IN)           :: time
!!$    REAL(dp), INTENT(IN)           :: efield(:), afield(:)
!!$    INTEGER, INTENT(IN)            :: lmax
!!$    
!!$    INTEGER                        :: wave_rank, field_rank
!!$    INTEGER                        :: time_rank
!!$    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
!!$    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
!!$    REAL(dp)                       :: field(3,2)
!!$    INTEGER                        :: numthetapts, numphipts
!!$    ! File identifier
!!$    INTEGER(HID_T)                 :: file_id
!!$    ! Property list identifier 
!!$    INTEGER(HID_T)                 :: plist_id
!!$    ! Group identifiers
!!$    INTEGER(HID_T)                 :: wave_group
!!$    INTEGER(HID_T)                 :: wavederiv_group 
!!$    INTEGER(HID_T)                 :: field_group
!!$    INTEGER(HID_T)                 :: time_group
!!$    ! Dataspace identifier
!!$    INTEGER(HID_T)                 :: wave_dspace_id
!!$    INTEGER(HID_T)                 :: wavederiv_dspace_id
!!$    INTEGER(HID_T)                 :: field_dspace_id  
!!$    INTEGER(HID_T)                 :: time_dspace_id  
!!$    ! Dataset identifier  
!!$    INTEGER(HID_T)                 :: wave_dset_id
!!$    INTEGER(HID_T)                 :: wavederiv_dset_id
!!$    INTEGER(HID_T)                 :: field_dset_id
!!$    INTEGER(HID_T)                 :: time_dset_id
!!$    ! Dataset dimensions 
!!$    INTEGER(HSIZE_T)               :: wave_dims(3)
!!$    INTEGER(HSIZE_T)               :: field_dims(2)
!!$    INTEGER(HSIZE_T)               :: time_dim(1)
!!$    ! Error flag
!!$    INTEGER                        :: error
!!$    CHARACTER(LEN=6)               :: cstep
!!$    CHARACTER(LEN=100)             :: name
!!$    CHARACTER(LEN = 100)           :: setname
!!$    
!!$    !------------------------------------------------------!
!!$    
!!$    ! Increase in 1 the count
!!$    dataset_count = dataset_count + 1
!!$    WRITE(cstep,'(I6.6)') dataset_count
!!$    
!!$    ! Initialize fortran interface
!!$    CALL h5open_f(error)
!!$    
!!$    ! Open files and groups
!!$    name = TRIM(filename) // '.h5'
!!$    CALL h5fopen_f(name, H5F_ACC_RDWR_F, file_id, error)
!!$    
!!$    name = '/wave'
!!$    CALL h5gopen_f(file_id, name, wave_group, error)
!!$    name = '/wavederiv'
!!$    CALL h5gopen_f(file_id, name, wavederiv_group, error)
!!$    name = '/field'
!!$    CALL h5gopen_f(file_id, name, field_group, error)
!!$    name = '/time'
!!$    CALL h5gopen_f(file_id, name, time_group, error)
!!$    
!!$    ! Number of theta points in a serial run, its
!!$    numthetapts = lmax + 1
!!$    numphipts = lmax + 1
!!$    !Define rank, in this case 2.
!!$    wave_rank = 3
!!$    field_rank = 2
!!$    time_rank = 1
!!$    wave_dims = (/numthetapts, numphipts, 2/)
!!$    field_dims = (/3,2/)
!!$    time_dim = 1
!!$    
!!$    ! Write field array
!!$    ALLOCATE(complex_wave(numthetapts,numphipts,2))
!!$    ALLOCATE(complex_wavederiv(numthetapts,numphipts,2))
!!$    complex_wave(:,:,1) = REAL(wave,dp)
!!$    complex_wave(:,:,2) = AIMAG(wave)
!!$    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
!!$    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
!!$    field(:,1) = (/ efield(1), efield(2), efield(3) /)
!!$    field(:,2) = (/ afield(1), afield(2), afield(3) /)
!!$    
!!$    ! Create the data spaces for the  datasets. 
!!$    CALL h5screate_simple_f(wave_rank, wave_dims, wave_dspace_id, error)
!!$    CALL h5screate_simple_f(wave_rank, wave_dims, wavederiv_dspace_id, error)
!!$    CALL h5screate_simple_f(field_rank, field_dims, field_dspace_id, error)
!!$    CALL h5screate_simple_f(time_rank, time_dim, time_dspace_id, error)
!!$    
!!$    setname = '/wave/'// cstep
!!$    ! Create the dataset.
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$         wave_dspace_id, wave_dset_id, error)
!!$    CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, &
!!$         wave_dims, error)
!!$    
!!$    setname = '/wavederiv/'// cstep
!!$    ! Create the dataset.
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$         wavederiv_dspace_id, wavederiv_dset_id, error)
!!$    CALL h5dwrite_f(wavederiv_dset_id, H5T_NATIVE_DOUBLE, complex_wavederiv, &
!!$         wave_dims, error)
!!$    
!!$    setname = '/field/'// cstep
!!$    ! Create the dataset.
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$         field_dspace_id, field_dset_id, error)
!!$    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, &
!!$         field_dims, error)
!!$    
!!$    setname = '/time/'// cstep
!!$    ! Create the dataset.
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$         time_dspace_id, time_dset_id, error)
!!$    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, &
!!$         time_dim, error)
!!$
!!$    ! End access to the dataset and release 
!!$    ! resources used by it.
!!$    CALL h5dclose_f(wave_dset_id, error)
!!$    CALL h5dclose_f(wavederiv_dset_id, error)
!!$    CALL h5dclose_f(field_dset_id, error)
!!$    CALL h5dclose_f(time_dset_id, error)
!!$    
!!$    ! Terminate access to the data space.
!!$    CALL h5sclose_f(wave_dspace_id, error)
!!$    CALL h5sclose_f(wavederiv_dspace_id, error)
!!$    CALL h5sclose_f(field_dspace_id, error)
!!$    CALL h5sclose_f(time_dspace_id, error)
!!$    
!!$    ! Close the property list, group, file and interface
!!$    CALL h5gclose_f(wave_group, error)
!!$    CALL h5gclose_f(wavederiv_group, error)
!!$    CALL h5gclose_f(field_group, error)
!!$    CALL h5gclose_f(time_group, error)
!!$    CALL h5fclose_f(file_id, error)
!!$    CALL h5close_f(error)
!!$    
!!$    ! Deallocate aux arrays for data
!!$    DEALLOCATE(complex_wave,complex_wavederiv)
!!$    
!!$  END SUBROUTINE write_surface_file3D_serial
!!$  
!!$  !******************************************************!
!!$#if _COM_MPI
!!$  SUBROUTINE write_surface_file3D_parallel(filename, wave, &
!!$       wavederiv, time, efield, afield, lmax, &
!!$       surfacerank,numsurfaceprocs, comm, &
!!$       numthetaptsperproc, numphiptsperproc )
!!$    IMPLICIT NONE
!!$    
!!$    CHARACTER(LEN=*), INTENT(IN)   :: filename
!!$    COMPLEX(dp), INTENT(IN)        :: wave(:, :)
!!$    COMPLEX(dp), INTENT(IN)        :: wavederiv(:, :)
!!$    REAL(dp), INTENT(IN)           :: time
!!$    REAL(dp), INTENT(IN)           :: efield(:), afield(:)
!!$    INTEGER, INTENT(IN)            :: lmax
!!$    INTEGER, INTENT(IN)            :: surfacerank
!!$    INTEGER, INTENT(IN)            :: numsurfaceprocs
!!$    INTEGER, INTENT(IN)            :: comm
!!$    INTEGER, INTENT(IN)            :: numthetaptsperproc
!!$    INTEGER, INTENT(IN)            :: numphiptsperproc
!!$    
!!$    INTEGER                        :: wave_rank, field_rank
!!$    INTEGER                        :: time_rank
!!$    REAL(dp), ALLOCATABLE          :: complex_wave(:, :, :)
!!$    REAL(dp), ALLOCATABLE          :: complex_wavederiv(:, :, :)
!!$    REAL(dp)                       :: field(3,2)
!!$    INTEGER                        :: numthetapts, numphipts
!!$    INTEGER                        :: numthetaoffset, numphioffset 
!!$    ! File identifier
!!$    INTEGER(HID_T)                 :: file_id
!!$    ! Property list identifier 
!!$    INTEGER(HID_T)                 :: plist_id
!!$    ! Group identifiers
!!$    INTEGER(HID_T)                 :: wave_group
!!$    INTEGER(HID_T)                 :: wavederiv_group 
!!$    INTEGER(HID_T)                 :: field_group
!!$    INTEGER(HID_T)                 :: time_group
!!$    ! Dataspace identifier in file 
!!$    INTEGER(HID_T)                 :: wave_filespace
!!$    INTEGER(HID_T)                 :: wavederiv_filespace
!!$    INTEGER(HID_T)                 :: field_filespace
!!$    INTEGER(HID_T)                 :: time_filespace
!!$    ! Dataspace identifier in memory
!!$    INTEGER(HID_T)                 :: wave_memspace
!!$    INTEGER(HID_T)                 :: wavederiv_memspace
!!$    INTEGER(HID_T)                 :: field_memspace
!!$    INTEGER(HID_T)                 :: time_memspace
!!$    ! Dataspace identifier
!!$    !INTEGER(HID_T)                 :: wave_dspace_id
!!$    !INTEGER(HID_T)                 :: wavederiv_dspace_id
!!$    !INTEGER(HID_T)                 :: field_dspace_id  
!!$    !INTEGER(HID_T)                 :: time_dspace_id  
!!$    ! Dataset identifier  
!!$    INTEGER(HID_T)                 :: wave_dset_id
!!$    INTEGER(HID_T)                 :: wavederiv_dset_id
!!$    INTEGER(HID_T)                 :: field_dset_id
!!$    INTEGER(HID_T)                 :: time_dset_id
!!$    ! Dataset dimensions 
!!$    INTEGER(HSIZE_T)               :: wave_global_dims(3)
!!$    INTEGER(HSIZE_T)               :: wave_dims(3)
!!$    INTEGER(HSIZE_T)               :: field_dims(2)
!!$    INTEGER(HSIZE_T)               :: time_dim(1)
!!$    ! Memory variables, per proc
!!$    INTEGER(HSIZE_T)               :: cont(3)  
!!$    INTEGER(HSSIZE_T)              :: offset(3)
!!$    INTEGER(HSIZE_T)               :: cont_field(2) 
!!$    INTEGER(HSSIZE_T)              :: offset_field(2)  
!!$    INTEGER(HSIZE_T)               :: cont_time(1) 
!!$    INTEGER(HSSIZE_T)              :: offset_time(1)  
!!$    ! Error flag
!!$    INTEGER                        :: error
!!$    
!!$    CHARACTER(LEN=6)               :: cstep
!!$    CHARACTER(LEN=100)             :: name
!!$    CHARACTER(LEN = 100)           :: setname
!!$
!!$    !------------------------------------------------------!
!!$    
!!$    ! Increase in 1 the count
!!$    dataset_count = dataset_count + 1
!!$    WRITE(cstep,'(I6.6)') dataset_count
!!$    
!!$    ! Initialize fortran interface
!!$    CALL h5open_f(error)
!!$
!!$    ! Setup file access property list with
!!$    ! parallel I/O access.
!!$    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
!!$    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
!!$    
!!$    ! Open files and groups
!!$    name = TRIM(filename) // '.h5'
!!$    CALL h5fopen_f(name, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
!!$    CALL h5pclose_f(plist_id, error)
!!$    
!!$    name = '/wave'
!!$    CALL h5gopen_f(file_id, name, wave_group, error)
!!$    name = '/wavederiv'
!!$    CALL h5gopen_f(file_id, name, wavederiv_group, error)
!!$    name = '/field'
!!$    CALL h5gopen_f(file_id, name, field_group, error)
!!$    name = '/time'
!!$    CALL h5gopen_f(file_id, name, time_group, error)
!!$    
!!$    numthetapts = lmax + 1
!!$    numphipts = lmax + 1
!!$    numthetaoffset = INT(numthetapts / numsurfaceprocs)
!!$    numphioffset = INT(numphipts / numsurfaceprocs)
!!$    wave_rank = 3
!!$    field_rank = 2
!!$    time_rank = 1
!!$    field_dims = (/3,2/)
!!$    time_dim = 1
!!$    wave_global_dims = (/numthetapts, numphipts, 2/)
!!$    
!!$    wave_dims = (/numthetaptsperproc, numphiptsperproc, 2/)
!!$    ALLOCATE(complex_wave(numthetaptsperproc,numphiptsperproc,2))
!!$    ALLOCATE(complex_wavederiv(numthetaptsperproc,numphiptsperproc,2))
!!$    complex_wave(:,:,1) = REAL(wave,dp)
!!$    complex_wave(:,:,2) = AIMAG(wave)
!!$    complex_wavederiv(:,:,1) = REAL(wavederiv,dp)
!!$    complex_wavederiv(:,:,2) = AIMAG(wavederiv)
!!$    
!!$    ! Write field array   
!!$    field(:,1) = (/ efield(1), efield(2), efield(3) /)
!!$    field(:,2) = (/ afield(1), afield(2), afield(3) /)
!!$    
!!$    ! Create the data space for the dataset.
!!$    CALL h5screate_simple_f(wave_rank, wave_global_dims, wave_filespace, error)
!!$    CALL h5screate_simple_f(wave_rank, wave_dims, wave_memspace, error)
!!$    CALL h5screate_simple_f(wave_rank, wave_global_dims, wavederiv_filespace, error)
!!$    CALL h5screate_simple_f(wave_rank, wave_dims, wavederiv_memspace, error)
!!$    
!!$    ! Create dataset.
!!$    setname = '/wave/'// cstep
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, wave_filespace, &
!!$         wave_dset_id, error)!, plist_id)
!!$    CALL h5sclose_f(wave_filespace, error)
!!$    setname = '/wavederiv/'// cstep
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, wavederiv_filespace, &
!!$         wavederiv_dset_id, error)!, plist_id)
!!$    CALL h5sclose_f(wavederiv_filespace, error)
!!$    
!!$    cont = wave_dims
!!$    offset = (/ numthetaoffset * surfacerank, numphioffset * surfacerank, 0 /)
!!$    
!!$    ! Select hyperslab in the file.
!!$    CALL h5dget_space_f(wave_dset_id, wave_filespace, error)
!!$    CALL h5sselect_hyperslab_f (wave_filespace, H5S_SELECT_SET_F, offset,&
!!$         cont, error)!, stride, blck)
!!$    CALL h5dget_space_f(wavederiv_dset_id, wavederiv_filespace, error)
!!$    CALL h5sselect_hyperslab_f (wavederiv_filespace, H5S_SELECT_SET_F, offset,&
!!$         cont, error)!, stride, blck)
!!$    
!!$    ! Create property list for collective dataset write
!!$    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
!!$    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!!$    ! Write the dataset collectively.
!!$    CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, wave_global_dims, &
!!$         error, file_space_id = wave_filespace, mem_space_id = wave_memspace, &
!!$         xfer_prp = plist_id)
!!$    CALL h5dwrite_f(wavederiv_dset_id, H5T_NATIVE_DOUBLE, complex_wavederiv, wave_global_dims, &
!!$         error, file_space_id = wavederiv_filespace, mem_space_id = wavederiv_memspace, &
!!$         xfer_prp = plist_id)
!!$    
!!$    ! Close dataspaces.
!!$    CALL h5pclose_f(plist_id, error)
!!$    CALL h5sclose_f(wave_filespace, error)
!!$    CALL h5sclose_f(wave_memspace, error)
!!$    CALL h5sclose_f(wavederiv_filespace, error)
!!$    CALL h5sclose_f(wavederiv_memspace, error)
!!$    
!!$    ! Close the dataset.
!!$    CALL h5dclose_f(wave_dset_id, error)
!!$    CALL h5dclose_f(wavederiv_dset_id, error)
!!$    
!!$    ! Deallocate aux arrays for data
!!$    DEALLOCATE(complex_wave,complex_wavederiv)
!!$    
!!$    !************!
!!$    
!!$    ! For the field, only the first processor in the group
!!$    ! write the field data.
!!$    
!!$    ! Create the data space for the dataset.
!!$    CALL h5screate_simple_f(field_rank, field_dims, field_filespace, error)
!!$    CALL h5screate_simple_f(field_rank, field_dims, field_memspace, error)
!!$    IF(surfacerank.NE.1) &
!!$         CALL h5sselect_none_f(field_memspace,error)
!!$    
!!$    ! Create dataset
!!$    setname = '/field/'// cstep
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, field_filespace, &
!!$         field_dset_id, error)!, plist_id)
!!$    CALL h5sclose_f(field_filespace, error)
!!$    
!!$    cont_field = field_dims
!!$    offset_field = 0
!!$    ! Select hyperslab in the file.
!!$    CALL h5dget_space_f(field_dset_id, field_filespace, error)
!!$    CALL h5sselect_hyperslab_f (field_filespace, H5S_SELECT_SET_F, offset_field,&
!!$         cont_field, error)!, stride, blck)
!!$    
!!$    IF(surfacerank.NE.1) &
!!$         CALL h5sselect_none_f(field_filespace,error)
!!$    
!!$    ! Create property list for collective dataset write
!!$    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
!!$    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!!$    
!!$    !IF (pivot_number(1) .EQ. 1) &
!!$    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dims, &
!!$         error, file_space_id = field_filespace, mem_space_id = field_memspace, &
!!$         xfer_prp = plist_id)
!!$    
!!$    ! Close dataspace and dataset
!!$    CALL h5pclose_f(plist_id, error)
!!$    CALL h5sclose_f(field_filespace, error)
!!$    CALL h5sclose_f(field_memspace, error)
!!$    CALL h5dclose_f(field_dset_id, error)    
!!$
!!$    !----------------------------------------------------!
!!$
!!$    ! For the field, only the first processor in the group
!!$    ! write the field data.
!!$    
!!$    ! Create the data space for the dataset.
!!$    CALL h5screate_simple_f(time_rank, time_dim, time_filespace, error)
!!$    CALL h5screate_simple_f(time_rank, time_dim, time_memspace, error)
!!$    IF(surfacerank.NE.1) &
!!$         CALL h5sselect_none_f(time_memspace,error)
!!$    
!!$    ! Create dataset
!!$    setname = '/time/'// cstep
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, time_filespace, &
!!$         time_dset_id, error)!, plist_id)
!!$    CALL h5sclose_f(time_filespace, error)
!!$    
!!$    cont_time = time_dim
!!$    offset_time = 0
!!$    ! Select hyperslab in the file.
!!$    CALL h5dget_space_f(time_dset_id, time_filespace, error)
!!$    CALL h5sselect_hyperslab_f (time_filespace, H5S_SELECT_SET_F, offset_time,&
!!$         cont_time, error)!, stride, blck)
!!$    
!!$    IF(surfacerank.NE.1) &
!!$         CALL h5sselect_none_f(time_filespace,error)
!!$    
!!$    ! Create property list for collective dataset write
!!$    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
!!$    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
!!$    
!!$    !IF (pivot_number(1) .EQ. 1) &
!!$    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_dim, &
!!$         error, file_space_id = time_filespace, mem_space_id = time_memspace, &
!!$         xfer_prp = plist_id)
!!$    
!!$    ! Close dataspace and dataset
!!$    CALL h5pclose_f(plist_id, error)
!!$    CALL h5sclose_f(time_filespace, error)
!!$    CALL h5sclose_f(time_memspace, error)
!!$    CALL h5dclose_f(time_dset_id, error)    
!!$    
!!$    !------------------------------------!
!!$    
!!$    ! Close the property list, group, file and interface
!!$    CALL h5gclose_f(wave_group, error)
!!$    CALL h5gclose_f(wavederiv_group, error)
!!$    CALL h5gclose_f(field_group, error)
!!$    CALL h5gclose_f(time_group, error)
!!$    CALL h5fclose_f(file_id, error)
!!$    CALL h5close_f(error)
!!$    
!!$  END SUBROUTINE write_surface_file3D_parallel
!!$#endif
  
  !***************************************************!
  
!!$  SUBROUTINE get_surface_dims(filename, ntime, dt, &
!!$       numthetapts, numphipts, lmax)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    CHARACTER(LEN=*), INTENT(IN)      :: filename
!!$    INTEGER, INTENT(OUT)              :: ntime
!!$    REAL(dp), INTENT(OUT)             :: dt
!!$    INTEGER, INTENT(OUT)              :: numthetapts, numphipts
!!$    INTEGER, INTENT(OUT)              :: lmax
!!$
!!$    REAL(dp)                       :: t1, t2
!!$    ! File identifier
!!$    INTEGER(HID_T)                 :: file_id
!!$    ! Group identifiers
!!$    INTEGER(HID_T)                 :: group_id 
!!$    ! Dataspace identifier
!!$    INTEGER(HID_T)                 :: dspace_id
!!$    INTEGER(HID_T)                 :: memspace
!!$    ! Dataset identifier  
!!$    INTEGER(HID_T)                 :: dset_id      
!!$    ! Dataset dimensions  
!!$    INTEGER(HSIZE_T)               :: dspace_dims(3)
!!$    INTEGER(HSIZE_T)               :: dspace_maxdims(3)
!!$    INTEGER(HSIZE_T)               :: time_dim(1)
!!$    INTEGER(HSIZE_T)               :: time_offset(1)
!!$    ! Structure info
!!$    INTEGER                        :: storage_type
!!$    INTEGER                        :: nlinks
!!$    INTEGER                        :: maxcorder
!!$    ! Error flag
!!$    INTEGER                        :: error     
!!$    ! Auxiliary strings
!!$    CHARACTER(LEN = 100)           :: name
!!$    
!!$    !----------------------------------------------!
!!$    
!!$    ! Initialize FORTRAN interface for hdf5.
!!$    CALL h5open_f(error)    
!!$    
!!$    ! Open file
!!$    name = TRIM(filename) // '.h5'
!!$    CALL h5fopen_f(name,H5F_ACC_RDONLY_F, file_id, error) 
!!$    
!!$    ! Open wave group, to see how many times are on it.
!!$    name = '/wave'
!!$    CALL h5gopen_f(file_id, name, group_id, error)
!!$    ! Get file info
!!$    CALL h5gget_info_f(group_id, storage_type, nlinks, &
!!$         maxcorder, error)
!!$
!!$    ! Assign number of time steps
!!$    ntime = nlinks
!!$    
!!$    name = '/wave/000000'
!!$    CALL h5dopen_f(file_id, name, dset_id, error)
!!$
!!$    CALL h5dget_space_f(dset_id, dspace_id, error)
!!$    
!!$    CALL h5sget_simple_extent_dims_f(dspace_id, dspace_dims, &
!!$         dspace_maxdims, error)
!!$    
!!$    IF(dspace_dims(3).EQ.0) THEN
!!$       numthetapts = dspace_dims(1)
!!$       numphipts = 1
!!$       lmax = numthetapts - 1
!!$    ELSE
!!$       numthetapts = dspace_dims(1)
!!$       numphipts = dspace_dims(2)
!!$       lmax = numthetapts - 1
!!$    ENDIF
!!$    
!!$    ! Close dataspaces and sets opened so far.
!!$    CALL h5sclose_f(dspace_id, error)
!!$    CALL h5dclose_f(dset_id, error)
!!$    CALL h5gclose_f(group_id, error)
!!$
!!$    time_dim = 1
!!$    time_offset = 0
!!$    t1 = 0.0_dp
!!$    t2 = 0.0_dp
!!$    
!!$    IF(ntime.GT.1) THEN
!!$       ! Open time group
!!$       name = '/time'
!!$       CALL h5gopen_f(file_id, name, group_id, error)
!!$       
!!$       ! Calculate t1
!!$       name = '/time/000000'
!!$       CALL h5dopen_f(file_id, name, dset_id, error)
!!$       
!!$       CALL h5dget_space_f(dset_id, dspace_id, error)
!!$       
!!$       CALL h5sget_simple_extent_dims_f(dspace_id, dspace_dims, &
!!$            dspace_maxdims, error)
!!$       
!!$       CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
!!$            time_offset, time_dim, error)
!!$       CALL h5screate_simple_f(1,time_dim,memspace,error)
!!$       
!!$       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, t1, time_dim, error, &
!!$            memspace, dspace_id )
!!$       
!!$       CALL h5sclose_f(dspace_id, error)
!!$       CALL h5sclose_f(memspace, error)
!!$       CALL h5dclose_f(dset_id, error)
!!$       
!!$       ! Calculate t2
!!$       name = '/time/000001'
!!$       CALL h5dopen_f(file_id, name, dset_id, error)
!!$       
!!$       CALL h5dget_space_f(dset_id, dspace_id, error)
!!$       
!!$       CALL h5sget_simple_extent_dims_f(dspace_id, dspace_dims, &
!!$            dspace_maxdims, error)
!!$       
!!$       CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, &
!!$            time_offset, time_dim, error)
!!$       CALL h5screate_simple_f(1,time_dim,memspace,error)
!!$       
!!$       CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, t2, time_dim, error, &
!!$            memspace, dspace_id )
!!$
!!$       ! Close again!
!!$       CALL h5sclose_f(dspace_id, error)
!!$       CALL h5sclose_f(memspace, error)
!!$       CALL h5dclose_f(dset_id, error)
!!$       CALL h5gclose_f(group_id,error)
!!$    ELSE
!!$       WRITE(*,*) 'File contains only 1 or 0 timesteps.'
!!$    ENDIF
!!$
!!$    ! Determine time step
!!$    dt = t2 - t1
!!$    
!!$    WRITE(*,*)
!!$    WRITE(*,*)           '---------------------------------------------------------&
!!$         &----------------------'
!!$    WRITE(*,'(A,A,A)')   'Path and name of the file:      ',TRIM(filename) // '.h5'
!!$    WRITE(*,'(A,I8)')    'Number of time step in the file: ',ntime
!!$    WRITE(*,'(A,F0.6)')  'Time step:                       ',dt 
!!$    WRITE(*,'(A,I4)')    'Number of theta points:          ',numthetapts
!!$    WRITE(*,'(A,I4)')    'Number of phi points:            ',numphipts
!!$    WRITE(*,'(A,I4)')    'Maximum angular momenta:         ',lmax
!!$    
!!$    ! Close the file.
!!$    CALL h5fclose_f(file_id, error)
!!$    
!!$    ! Close HDF5 fortran interface
!!$    CALL h5close_f(error)
!!$    
!!$  END SUBROUTINE get_surface_dims
!!$
!!$
!!$  !***************************************************!
!!$  
!!$  SUBROUTINE read_surface(filename, itime, numthetapts, &
!!$       numphipts, psi_sph, psip_sph, time, efield, afield)
!!$
!!$    IMPLICIT NONE
!!$    
!!$    CHARACTER(LEN=*), INTENT(IN)       :: filename
!!$    INTEGER, INTENT(IN)                :: itime
!!$    INTEGER, INTENT(IN)                :: numthetapts
!!$    INTEGER, INTENT(IN)                :: numphipts
!!$    COMPLEX(dp), INTENT(OUT)           :: psi_sph(:, :)
!!$    COMPLEX(dp), INTENT(OUT)           :: psip_sph(:, :)
!!$    REAL(dp), INTENT(OUT)              :: time
!!$    REAL(dp), INTENT(OUT)              :: efield(:), afield(:)
!!$    
!!$    REAL(dp), ALLOCATABLE         :: psi2D(:, :)
!!$    REAL(dp), ALLOCATABLE         :: psip2D(:, :)
!!$    REAL(dp), ALLOCATABLE         :: psi3D(:, :, :)
!!$    REAL(dp), ALLOCATABLE         :: psip3D(:, :, :)
!!$    REAL(dp)                      :: field(3,2)
!!$    INTEGER                       :: itheta, iphi
!!$    ! File identifier
!!$    INTEGER(HID_T)                 :: file_id
!!$    ! Group identifiers
!!$    INTEGER(HID_T)                 :: wave_group_id 
!!$    INTEGER(HID_T)                 :: wavep_group_id 
!!$    INTEGER(HID_T)                 :: field_group_id 
!!$    INTEGER(HID_T)                 :: time_group_id 
!!$    ! Dataspace identifier
!!$    INTEGER(HID_T)                 :: wave_dspace_id    
!!$    INTEGER(HID_T)                 :: wavep_dspace_id    
!!$    INTEGER(HID_T)                 :: field_dspace_id    
!!$    INTEGER(HID_T)                 :: time_dspace_id    
!!$    ! Dataset identifier  
!!$    INTEGER(HID_T)                 :: wave_dset_id      
!!$    INTEGER(HID_T)                 :: wavep_dset_id      
!!$    INTEGER(HID_T)                 :: field_dset_id      
!!$    INTEGER(HID_T)                 :: time_dset_id      
!!$    ! Dataset dimensions  
!!$    INTEGER(HSIZE_T), ALLOCATABLE  :: dspace_dims(:)
!!$    INTEGER(HSIZE_T), ALLOCATABLE  :: dspace_offset(:)
!!$    INTEGER(HSIZE_T)               :: field_dims(2)
!!$    INTEGER(HSIZE_T)               :: field_offset(2)
!!$    INTEGER(HSIZE_T)               :: time_dim(1)
!!$    INTEGER(HSIZE_T)               :: time_offset(1)
!!$    ! Memspace identifier
!!$    INTEGER(HID_T)                 :: wave_memspace
!!$    INTEGER(HID_T)                 :: wavep_memspace
!!$    INTEGER(HID_T)                 :: field_memspace
!!$    INTEGER(HID_T)                 :: time_memspace
!!$    INTEGER                        :: error     
!!$    ! Auxiliary strings
!!$    CHARACTER(LEN=6)               :: cstep
!!$    CHARACTER(LEN = 100)           :: name
!!$
!!$    !---------------------------------------------------!
!!$    
!!$    ! Set time step as a string
!!$    WRITE(cstep,'(I6.6)') itime
!!$    
!!$    ! Initialize FORTRAN interface for hdf5.
!!$    CALL h5open_f(error)    
!!$    
!!$    ! Open file
!!$    name = TRIM(filename) // '.h5'
!!$    CALL h5fopen_f(name,H5F_ACC_RDONLY_F, file_id, error) 
!!$    
!!$    ! Open wave group. 
!!$    name = '/wave'
!!$    CALL h5gopen_f(file_id, name, wave_group_id, error)
!!$    name = '/wavederiv'
!!$    CALL h5gopen_f(file_id, name, wavep_group_id, error)
!!$    
!!$    ! Open dataset
!!$    name = '/wave/' // cstep
!!$    CALL h5dopen_f(file_id, name, wave_dset_id, error)
!!$    ! Get dataspace identifier
!!$    CALL h5dget_space_f(wave_dset_id, wave_dspace_id, error)
!!$    name = '/wavederiv/' // cstep
!!$    CALL h5dopen_f(file_id, name, wavep_dset_id, error)
!!$    ! Get dataspace identifier
!!$    CALL h5dget_space_f(wavep_dset_id, wavep_dspace_id, error)
!!$    
!!$    IF(numphipts.EQ.1) THEN
!!$       ALLOCATE(dspace_dims(2),dspace_offset(2))
!!$       dspace_dims = (/ numthetapts, 2/)
!!$       ALLOCATE(psi2D(1:numthetapts,1:2))
!!$       ALLOCATE(psip2D(1:numthetapts,1:2))
!!$    ELSE
!!$       ALLOCATE(dspace_dims(3),dspace_offset(3))
!!$       dspace_dims = (/ numthetapts, numphipts, 2/)
!!$       ALLOCATE(psi3D(1:numthetapts,1:numphipts,1:2))
!!$       ALLOCATE(psip3D(1:numthetapts,1:numphipts,1:2))
!!$    ENDIF
!!$    dspace_offset = 0
!!$    
!!$    CALL h5sselect_hyperslab_f(wave_dspace_id, H5S_SELECT_SET_F, &
!!$         dspace_offset, dspace_dims, error)
!!$    CALL h5sselect_hyperslab_f(wavep_dspace_id, H5S_SELECT_SET_F, &
!!$         dspace_offset, dspace_dims, error)
!!$    
!!$    ! Create memspace and read dataset
!!$    IF(numphipts.EQ.1) THEN
!!$       CALL h5screate_simple_f(2,dspace_dims,wave_memspace,error)
!!$       CALL h5screate_simple_f(2,dspace_dims,wavep_memspace,error)
!!$       
!!$       CALL h5dread_f(wave_dset_id, H5T_NATIVE_DOUBLE, psi2D, dspace_dims, error, &
!!$            wave_memspace, wave_dspace_id )
!!$       CALL h5dread_f(wavep_dset_id, H5T_NATIVE_DOUBLE, psip2D, dspace_dims, error, &
!!$            wavep_memspace, wavep_dspace_id )
!!$       
!!$       DO itheta = 1, numthetapts
!!$          psi_sph(itheta,1) = CMPLX(psi2D(itheta,1), psi2D(itheta,2), dp)
!!$          psip_sph(itheta,1) = CMPLX(psip2D(itheta,1), psip2D(itheta,2), dp)
!!$       ENDDO
!!$       
!!$    ELSE
!!$       CALL h5screate_simple_f(3,dspace_dims,wave_memspace,error)
!!$       CALL h5screate_simple_f(3,dspace_dims,wavep_memspace,error)
!!$       
!!$       CALL h5dread_f(wave_dset_id, H5T_NATIVE_DOUBLE, psi2D, dspace_dims, error, &
!!$            wave_memspace, wave_dspace_id )
!!$       CALL h5dread_f(wavep_dset_id, H5T_NATIVE_DOUBLE, psip2D, dspace_dims, error, &
!!$            wavep_memspace, wavep_dspace_id )
!!$       
!!$       DO iphi = 1, numphipts
!!$          DO itheta = 1, numthetapts
!!$             psi_sph(itheta,iphi) = CMPLX(psi3D(itheta,iphi,1), psi3D(itheta,iphi,2), dp)
!!$             psip_sph(itheta,iphi) = CMPLX(psip3D(itheta,iphi,1), psip3D(itheta,iphi,2), dp)
!!$          ENDDO
!!$       ENDDO
!!$    ENDIF
!!$    
!!$    ! Close (almost) everything
!!$    CALL h5sclose_f(wave_dspace_id, error)
!!$    CALL h5sclose_f(wave_memspace, error)
!!$    CALL h5dclose_f(wave_dset_id, error)
!!$    CALL h5gclose_f(wave_group_id, error)
!!$    
!!$    CALL h5sclose_f(wavep_dspace_id, error)
!!$    CALL h5sclose_f(wavep_memspace, error)
!!$    CALL h5dclose_f(wavep_dset_id, error)
!!$    CALL h5gclose_f(wavep_group_id, error)
!!$    
!!$    DEALLOCATE(dspace_dims, dspace_offset)
!!$    IF(numphipts.EQ.1) THEN
!!$       DEALLOCATE(psi2D, psip2D)
!!$    ELSE
!!$       DEALLOCATE(psi3D, psip3D)
!!$    ENDIF
!!$    
!!$    !
!!$    ! Now, we read the field
!!$    !
!!$    
!!$    ! Open field group. 
!!$    name = '/field'
!!$    CALL h5gopen_f(file_id, name, field_group_id, error)
!!$    
!!$    ! Open dataset
!!$    name = '/field/' // cstep
!!$    CALL h5dopen_f(file_id, name, field_dset_id, error)
!!$    ! Get dataspace identifier
!!$    CALL h5dget_space_f(field_dset_id, field_dspace_id, error)
!!$    
!!$    field_dims = (/3,2/)
!!$    field_offset = 0
!!$    field = 0
!!$    
!!$    CALL h5sselect_hyperslab_f(field_dspace_id, H5S_SELECT_SET_F, &
!!$         field_offset, field_dims, error)
!!$    
!!$    CALL h5screate_simple_f(2,field_dims,field_memspace,error)
!!$    
!!$    CALL h5dread_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dims, error, &
!!$         field_memspace, field_dspace_id )
!!$    
!!$    efield(:) = field(:,1)
!!$    afield(:) = field(:,2)
!!$    
!!$    CALL h5sclose_f(field_dspace_id, error)
!!$    CALL h5sclose_f(field_memspace, error)
!!$    CALL h5dclose_f(field_dset_id, error)
!!$    CALL h5gclose_f(field_group_id, error)
!!$    
!!$    !--------!
!!$    
!!$    ! Open field group. 
!!$    name = '/time'
!!$    CALL h5gopen_f(file_id, name, time_group_id, error)
!!$    
!!$    ! Open dataset
!!$    name = '/time/' // cstep
!!$    CALL h5dopen_f(file_id, name, time_dset_id, error)
!!$    ! Get dataspace identifier
!!$    CALL h5dget_space_f(time_dset_id, time_dspace_id, error)
!!$    
!!$    time_dim = 1
!!$    time_offset = 0
!!$    time = 0.0_dp
!!$    
!!$    CALL h5sselect_hyperslab_f(time_dspace_id, H5S_SELECT_SET_F, &
!!$         time_offset, time_dim, error)
!!$    
!!$    CALL h5screate_simple_f(1,time_dim,time_memspace,error)
!!$    
!!$    CALL h5dread_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_dim, error, &
!!$         time_memspace, time_dspace_id )
!!$    
!!$    CALL h5sclose_f(time_dspace_id, error)
!!$    CALL h5sclose_f(time_memspace, error)
!!$    CALL h5dclose_f(time_dset_id, error)
!!$    CALL h5gclose_f(time_group_id, error)
!!$    
!!$    ! Close the file.
!!$    CALL h5fclose_f(file_id, error)
!!$    
!!$    ! Close HDF5 fortran interface
!!$    CALL h5close_f(error)
!!$    
!!$  END SUBROUTINE read_surface
  
  !***************************************************!
  !***************************************************!
  !***************************************************!
  
END MODULE io_corrsurface
