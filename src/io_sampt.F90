MODULE io_sampt

  USE constants
  USE fourier
  USE HDF5
#if _COM_MPI
  USE MPI
#endif
  IMPLICIT NONE

  PRIVATE
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public subroutines                                                        !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  PUBLIC                        :: initialize_detector
  PUBLIC                        :: create_detector_file
  PUBLIC                        :: write_detector_points
  PUBLIC                        :: get_detector_params
  PUBLIC                        :: read_detector
  PUBLIC                        :: get_solid_angles
  PUBLIC                        :: calculate_sampt_pes

  INTERFACE initialize_detector
     MODULE PROCEDURE initialize_detector_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_detector_parallel
#endif
  END INTERFACE initialize_detector

  INTERFACE create_detector_file
     MODULE PROCEDURE create_detector_file_serial
#if _COM_MPI
     MODULE PROCEDURE create_detector_file_parallel
#endif
  END INTERFACE create_detector_file

  INTERFACE write_detector_points
     MODULE PROCEDURE write_detector_points_serial
!#if _COM_MPI
     !MODULE PROCEDURE write_detector_points_parallel
!#endif
  END INTERFACE write_detector_points

  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public variables                                                          !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Private variables                                                         !
  !                                                                           !
  !---------------------------------------------------------------------------!

  INTEGER                        :: numdetectorpts
  REAL(dp), ALLOCATABLE          :: detectorpoints(:, :)
  INTEGER, ALLOCATABLE           :: detectorpoints_indx(:, :)
  INTEGER, ALLOCATABLE           :: i_am_in_local(:)
  INTEGER, ALLOCATABLE           :: i_am_detector(:)
  INTEGER, ALLOCATABLE           :: detector_members(:)
  INTEGER                        :: numdetectorprocs
  INTEGER                        :: dataset_count
  
CONTAINS
  
  SUBROUTINE initialize_detector_serial(nodetectorpts, detectorpts, x_ax, y_ax, z_ax )
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)            :: nodetectorpts
    REAL(dp), INTENT(IN)           :: detectorpts(:, :)
    REAL(dp), INTENT(IN)           :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN), OPTIONAL :: z_ax(:)

    INTEGER                        :: ipts, inum
    
    !------------------------------------------------------------------------!
    
    numdetectorpts = nodetectorpts
    ! Allocate array for indexes
    IF(PRESENT(z_ax)) THEN
       ALLOCATE(detectorpoints_indx(1:numdetectorpts,3))
       ALLOCATE(detectorpoints(1:numdetectorpts,3))
    ELSE
       ALLOCATE(detectorpoints_indx(1:numdetectorpts,2))
       ALLOCATE(detectorpoints(1:numdetectorpts,2))
    ENDIF
    ALLOCATE(i_am_in_local(numdetectorpts))
    
    ! Calculate what array index is for each point
    ! Check if it is out of extent
    DO ipts = 1, numdetectorpts
       IF( (detectorpts(ipts,1).GT.MAXVAL(x_ax)) &
            .OR. (detectorpts(ipts,1).LT.MINVAL(x_ax)) &
            .OR. (detectorpts(ipts,2).GT.MAXVAL(y_ax)) &
            .OR. (detectorpts(ipts,2).LT.MINVAL(y_ax))) THEN
          IF(PRESENT(z_ax)) THEN
             IF( (detectorpts(ipts,3) .GT. MAXVAL(z_ax)) &
                  .OR. (detectorpts(ipts,3) .LT. MINVAL(z_ax)) )&
                  detectorpoints_indx(ipts,:) = 0
          ENDIF
       ELSE
          i_am_in_local(ipts) = 1
          detectorpoints_indx(ipts,1) = MINLOC(x_ax,1,x_ax.GE.detectorpts(ipts,1))
          detectorpoints(ipts,1) = x_ax(detectorpoints_indx(ipts,1))
          
          detectorpoints_indx(ipts,2) = MINLOC(y_ax,1,y_ax.GE.detectorpts(ipts,2))
          detectorpoints(ipts,2) = y_ax(detectorpoints_indx(ipts,2))
          
          IF(PRESENT(z_ax)) THEN
             detectorpoints_indx(ipts,3) = MINLOC(z_ax,1,z_ax.GE.detectorpts(ipts,3))
             detectorpoints(ipts,3) = z_ax(detectorpoints_indx(ipts,3))
          ENDIF
       ENDIF
    ENDDO
    
  END SUBROUTINE initialize_detector_serial
  
  !**************************************************************!
#if _COM_MPI 
  SUBROUTINE initialize_detector_parallel(nodetectorpts, detectorpts, &
       mpi_rank,mpi_size, comm, newcomm, detector_rank, &
       x_ax, y_ax, z_ax )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: nodetectorpts
    REAL(dp), INTENT(IN)           :: detectorpts(:, :)
    INTEGER, INTENT(IN)            :: mpi_rank
    INTEGER, INTENT(IN)            :: mpi_size
    INTEGER, INTENT(IN)            :: comm
    INTEGER, INTENT(OUT)           :: newcomm
    INTEGER, INTENT(OUT)           :: detector_rank
    REAL(dp), INTENT(IN)           :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN), OPTIONAL :: z_ax(:)

    INTEGER, ALLOCATABLE           :: i_am_detector_local(:)
    INTEGER                        :: simgroup, detectorgroup
    INTEGER                        :: ierror
    INTEGER                        :: ipts, inum, ii
    
    !------------------------------------------------------------------------!
    
    ! Allocate arrays
    numdetectorpts = nodetectorpts
    IF(PRESENT(z_ax)) THEN
       ALLOCATE(detectorpoints_indx(1:numdetectorpts,3))
       ALLOCATE(detectorpoints(1:numdetectorpts,3))
    ELSE
       ALLOCATE(detectorpoints_indx(1:numdetectorpts,2))
       ALLOCATE(detectorpoints(1:numdetectorpts,2))
    ENDIF
    ALLOCATE(i_am_in_local(numdetectorpts))
    ALLOCATE(i_am_detector_local(0:mpi_size-1))
    ALLOCATE(i_am_detector(0:mpi_size-1))
    i_am_in_local       = 0
    i_am_detector_local = 0
    i_am_detector       = 0
    
    
    ! Calculate what array index is for each point
    ! Check if it is out of extent
    DO ipts = 1, numdetectorpts
       IF( (detectorpts(ipts,1).GT.MAXVAL(x_ax)) &
            .OR. (detectorpts(ipts,1).LT.MINVAL(x_ax)) &
            .OR. (detectorpts(ipts,2).GT.MAXVAL(y_ax)) &
            .OR. (detectorpts(ipts,2).LT.MINVAL(y_ax))) THEN
          IF(PRESENT(z_ax)) THEN
               IF( (detectorpts(ipts,3) .GT. MAXVAL(z_ax)) &
               .OR. (detectorpts(ipts,3) .LT. MINVAL(z_ax))) &
               detectorpoints_indx(ipts,:) = 0
            ENDIF
       ELSE
          detectorpoints_indx(ipts,1) = MINLOC(x_ax,1,x_ax.GE.detectorpts(ipts,1))
          detectorpoints(ipts,1) = x_ax(detectorpoints_indx(ipts,1))

          detectorpoints_indx(ipts,2) = MINLOC(y_ax,1,y_ax.GE.detectorpts(ipts,2))
          detectorpoints(ipts,2) = y_ax(detectorpoints_indx(ipts,2))
          
          IF(PRESENT(z_ax)) THEN
             detectorpoints_indx(ipts,3) = MINLOC(z_ax,1,z_ax.GE.detectorpts(ipts,3))
             detectorpoints(ipts,3) = z_ax(detectorpoints_indx(ipts,3))
          ENDIF
          i_am_in_local(ipts) = 1
       ENDIF
    ENDDO

    ! Raise this flag to know this processor is participating.
    i_am_detector_local(mpi_rank) = 1
    
    ! Communicate to all processors if they are
    ! at the surface
    CALL MPI_ALLREDUCE(i_am_detector_local, i_am_detector, mpi_size, &
         MPI_INTEGER, MPI_SUM, comm, ierror )
    
    numdetectorprocs = SUM(i_am_detector)
    
    ALLOCATE(detector_members(0:numdetectorprocs-1))
    ii = -1
    DO inum = 0, mpi_size-1
       IF(i_am_detector(inum).EQ.1) THEN
          ii = ii + 1
          detector_members(ii) = inum
       ENDIF
    ENDDO
    
    ! Create a global group
    CALL MPI_COMM_GROUP(comm, simgroup, ierror)
    
    CALL MPI_GROUP_INCL(simgroup, numdetectorprocs, detector_members, &
         detectorgroup, ierror)
    ! Actually, this line creates the communicator
    CALL MPI_COMM_CREATE(comm, detectorgroup, newcomm, ierror)
    
    IF(i_am_detector_local(mpi_rank).EQ.1) &
         CALL MPI_COMM_RANK(newcomm, detector_rank, ierror)
       
  END SUBROUTINE initialize_detector_parallel
#endif
  !**************************************************************!
  
  SUBROUTINE create_detector_file_serial(filename, rank)
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(IN)            :: rank

    INTEGER                        :: wave_rank
    INTEGER                        :: points_rank
    INTEGER                        :: time_rank
    INTEGER                        :: field_rank
    REAL(dp)                       :: complex_wave(2)
    REAL(dp)                       :: field(2,3)
    REAL(dp)                       :: time(1)
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc1 = (/1/)
    INTEGER(HSIZE_T), DIMENSION(2) :: dimsc2 = (/1,1/)
    INTEGER(HSIZE_T), DIMENSION(3) :: dimsc3 = (/1,1,1/)
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                 :: crp_list
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: points_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id  
    INTEGER(HID_T)                 :: time_dspace_id  
    ! Dataset identifier
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: points_dset_id
    INTEGER(HID_T)                 :: field_dset_id  
    INTEGER(HID_T)                 :: time_dset_id 
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(2)
    INTEGER(HSIZE_T)               :: points_dims(2)
    INTEGER(HSIZE_T)               :: field_dims(3)
    INTEGER(HSIZE_T)               :: time_dims(1)
    INTEGER(HSIZE_T)               :: maxdims1D(1)
    INTEGER(HSIZE_T)               :: maxdims2D(2)
    INTEGER(HSIZE_T)               :: maxdims3D(3)
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name, setname
    CHARACTER(LEN=6)               :: cstep, cpts
    INTEGER                        :: ipts
    !--------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)
    
    ! Create the file.
    name = TRIM(filename) // '.h5'  
    CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error)
    
    !-------- Dataset detector points ---------------!
    ! Create the data spaces for points dataset.
    points_rank = rank
    points_dims = (/numdetectorpts,rank/)
    CALL h5screate_simple_f(points_rank, points_dims, points_dspace_id, error)
    
    setname = '/detector_points'
    ! Create the dataset for points.
    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
         points_dspace_id, points_dset_id, error)
    CALL h5dwrite_f(points_dset_id, H5T_NATIVE_DOUBLE, detectorpoints, &
         points_dims, error)
    
    ! End access to the dataset and release 
    ! resources used by it.
    CALL h5dclose_f(points_dset_id, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(points_dspace_id, error)
    !-------------------------------------------------!

    ! Create group in the root group.     
    name = '/wavefunctions'     
    CALL h5gcreate_f(file_id, name, wave_group, error)
    
    ! Create the extended datasets for the wavefunction,
    ! field and time.
    !Create the data space with unlimited dimensions.
    maxdims2D = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
    wave_dims = (/1,2/)
    wave_rank = 2
    !Modify dataset creation properties, i.e. enable chunking
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    CALL h5pset_chunk_f(crp_list, wave_rank, dimsc2, error)
    
    DO ipts = 1, numdetectorpts
       CALL h5screate_simple_f(wave_rank, wave_dims, wave_dspace_id, error, maxdims2D)
       
       !Create a dataset using cparms creation properties
       WRITE(cpts,'(I6.6)') ipts
       name = '/wavefunctions/'// cpts
       CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, wave_dspace_id, &
            wave_dset_id, error, crp_list)
       CALL h5sclose_f(wave_dspace_id, error)
       
       !Write data array to dataset
       complex_wave = (/0.0_dp, 0.0_dp/)
       CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, wave_dims, error)
       
       CALL h5dclose_f(wave_dset_id, error)
    ENDDO
    ! Close property list
    CALL h5pclose_f(crp_list, error)
    CALL h5gclose_f(wave_group, error) 
    !-------------------------------!
    ! Field dataset
    field_dims = (/1,2,3/)
    field_rank = 3
    maxdims3D = (/H5S_UNLIMITED_F,H5S_UNLIMITED_F,H5S_UNLIMITED_F/)
    !Modify dataset creation properties, i.e. enable chunking
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    CALL h5pset_chunk_f(crp_list, field_rank, dimsc3, error)
    
    CALL h5screate_simple_f(field_rank, field_dims, field_dspace_id, error, maxdims3D)
    !Create a dataset using cparms creation properties
    name = '/field'
    CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, field_dspace_id, &
         field_dset_id, error, crp_list )
    CALL h5sclose_f(field_dspace_id, error)
    
    !Write data array to dataset
    field(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    field(2,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dims, error)
    
    CALL h5dclose_f(field_dset_id, error)
    ! Close property list
    CALL h5pclose_f(crp_list, error)   
    !-----------------------------------!
    
    ! Time dataset
    time_dims = 1
    time_rank = 1
    maxdims1D = (/H5S_UNLIMITED_F/)
    !Modify dataset creation properties, i.e. enable chunking
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    CALL h5pset_chunk_f(crp_list, time_rank, dimsc1, error)

    CALL h5screate_simple_f(time_rank, time_dims, time_dspace_id, error, maxdims1D)
    !Create a dataset using cparms creation properties
    name = '/time'
    CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, time_dspace_id, &
         time_dset_id, error, crp_list )
    CALL h5sclose_f(time_dspace_id, error)
    
    !Write data array to dataset
    time = 0.0_dp
    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_dims, error)
    
    CALL h5dclose_f(time_dset_id, error)
    ! Close property list
    CALL h5pclose_f(crp_list, error)   
    !-----------------------------------!
    
    ! Close the property list, group, file and interface
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
    ! Set the count of datasets to 1. It will start at 1.
    dataset_count = 1
    
  END SUBROUTINE create_detector_file_serial
  
  !******************************************************!
#if _COM_MPI
  SUBROUTINE create_detector_file_parallel(filename, rank, comm, &
       detector_rank)
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(IN)            :: rank
    INTEGER, INTENT(IN)            :: comm
    INTEGER, INTENT(IN)            :: detector_rank
    
    INTEGER                        :: wave_rank
    INTEGER                        :: points_rank
    INTEGER                        :: time_rank
    INTEGER                        :: field_rank
    REAL(dp)                       :: complex_wave(2)
    REAL(dp)                       :: field(2,3)
    REAL(dp)                       :: time(1)
    INTEGER(HSIZE_T), DIMENSION(1) :: dimsc1 = (/1/)
    INTEGER(HSIZE_T), DIMENSION(2) :: dimsc2 = (/1,1/)
    INTEGER(HSIZE_T), DIMENSION(3) :: dimsc3 = (/1,1,1/)
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                 :: crp_list
    INTEGER(HID_T)                 :: plist_id
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_filespace
    INTEGER(HID_T)                 :: points_filespace
    INTEGER(HID_T)                 :: field_filespace 
    INTEGER(HID_T)                 :: time_filespace
    INTEGER(HID_T)                 :: wave_memspace
    INTEGER(HID_T)                 :: points_memspace
    INTEGER(HID_T)                 :: field_memspace 
    INTEGER(HID_T)                 :: time_memspace    
    ! Dataset identifier
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: points_dset_id
    INTEGER(HID_T)                 :: field_dset_id  
    INTEGER(HID_T)                 :: time_dset_id 
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(2)
    INTEGER(HSIZE_T)               :: points_dims(2)
    INTEGER(HSIZE_T)               :: field_dims(3)
    INTEGER(HSIZE_T)               :: time_dims(1)
    INTEGER(HSIZE_T)               :: maxdims1D(1)
    INTEGER(HSIZE_T)               :: maxdims2D(2)
    INTEGER(HSIZE_T)               :: maxdims3D(3)
    ! Memory variables
    INTEGER(HSIZE_T)               :: wave_cont(2)
    INTEGER(HSIZE_T)               :: points_cont(2)
    INTEGER(HSIZE_T)               :: field_cont(3)
    INTEGER(HSIZE_T)               :: time_cont(1)
    INTEGER(HSSIZE_T)              :: wave_offset(2)
    INTEGER(HSSIZE_T)              :: points_offset(2)
    INTEGER(HSSIZE_T)              :: field_offset(3)
    INTEGER(HSSIZE_T)              :: time_offset(1)
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group
    ! Error flag
    INTEGER                        :: error, info
    ! Auxiliary strings
    CHARACTER(LEN = 100)           :: name, setname
    CHARACTER(LEN=6)               :: cstep, cpts
    INTEGER                        :: ipts
    !--------------------------------------------!
    
    ! Initialize FORTRAN interface for hdf5.
    CALL h5open_f(error)

    ! Setup file access property list with 
    ! parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
    
    ! Create the file.
    name = TRIM(filename) // '.h5'  
    CALL h5fcreate_f(name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    CALL h5pclose_f(plist_id, error)
    
!!$    !-------- Dataset detector points ---------------!
!!$    ! Create the data spaces for points dataset.
!!$    points_rank = rank
!!$    points_dims = (/numdetectorpts,rank/)
!!$    CALL h5screate_simple_f(points_rank, points_dims, points_dspace_id, error)
!!$    
!!$    setname = '/detector_points'
!!$    ! Create the dataset for points.
!!$    CALL h5dcreate_f(file_id, setname, H5T_NATIVE_DOUBLE, &
!!$         points_dspace_id, points_dset_id, error)
!!$    CALL h5dwrite_f(points_dset_id, H5T_NATIVE_DOUBLE, detectorpoints, &
!!$         points_dims, error)
!!$    
!!$    ! End access to the dataset and release 
!!$    ! resources used by it.
!!$    CALL h5dclose_f(points_dset_id, error)
!!$    ! Terminate access to the data space.
!!$    CALL h5sclose_f(points_dspace_id, error)
!!$    !-------------------------------------------------!
!!$
!!$    ! Create group in the root group.     
!!$    name = '/wavefunctions'     
!!$    CALL h5gcreate_f(file_id, name, wave_group, error)
!!$    
!!$    ! Create the extended datasets for the wavefunction,
!!$    ! field and time.
!!$    !Create the data space with unlimited dimensions.
!!$    maxdims2D = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
!!$    wave_dims = (/1,2/)
!!$    wave_rank = 2
!!$    !Modify dataset creation properties, i.e. enable chunking
!!$    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
!!$    CALL h5pset_chunk_f(crp_list, wave_rank, dimsc2, error)
!!$    
!!$    DO ipts = 1, numdetectorpts
!!$       CALL h5screate_simple_f(wave_rank, wave_dims, wave_dspace_id, error, maxdims2D)
!!$       
!!$       !Create a dataset using cparms creation properties
!!$       WRITE(cpts,'(I6.6)') ipts
!!$       name = '/wavefunctions/'// cpts
!!$       CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, wave_dspace_id, &
!!$            wave_dset_id, error, crp_list)
!!$       CALL h5sclose_f(wave_dspace_id, error)
!!$       
!!$       !Write data array to dataset
!!$       complex_wave = (/0.0_dp, 0.0_dp/)
!!$       CALL h5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, complex_wave, wave_dims, error)
!!$       
!!$       CALL h5dclose_f(wave_dset_id, error)
!!$    ENDDO
!!$    ! Close property list
!!$    CALL h5pclose_f(crp_list, error)
!!$    CALL h5gclose_f(wave_group, error) 
!!$    !-------------------------------!
    ! Field dataset
    field_dims = (/1,2,3/)
    field_rank = 3
    maxdims3D = (/H5S_UNLIMITED_F,H5S_UNLIMITED_F,H5S_UNLIMITED_F/)
    !Modify dataset creation properties, i.e. enable chunking
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    CALL h5pset_chunk_f(crp_list, field_rank, dimsc3, error)

    CALL h5screate_simple_f(field_rank, field_dims, field_filespace, error, maxdims3D)
    CALL h5screate_simple_f(field_rank, field_dims, field_memspace, error, maxdims3D)
    IF(detector_rank.NE.1) &
         CALL h5sselect_none_f(field_memspace,error)
    
    !Create a dataset using cparms creation properties
    name = '/field'
    CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, field_filespace, &
         field_dset_id, error, crp_list )
    CALL h5sclose_f(field_filespace, error)
    
    field_cont   = field_dims
    field_offset = 0
    ! Select hyperslab in the file.
    CALL h5dget_space_f(field_dset_id, field_filespace, error)
    CALL h5sselect_hyperslab_f (field_filespace, H5S_SELECT_SET_F, field_offset,&
         field_cont, error)!, stride, blck)
    IF(detector_rank.NE.1) &
         CALL h5sselect_none_f(field_filespace,error)
       
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    !Write data array to dataset
    field(1,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    field(2,:) = (/0.0_dp, 0.0_dp, 0.0_dp/)
    CALL h5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_dims, &
         error, file_space_id = field_filespace, mem_space_id = field_memspace, &
         xfer_prp = plist_id)
    
    CALL h5sclose_f(field_filespace, error)
    CALL h5sclose_f(field_memspace, error)
    CALL h5dclose_f(field_dset_id, error)
    CALL h5pclose_f(crp_list, error)
    CALL h5pclose_f(plist_id, error)
    !-----------------------------------!
    
    ! Time dataset
    time_dims = 1
    time_rank = 1
    maxdims1D = (/H5S_UNLIMITED_F/)
    !Modify dataset creation properties, i.e. enable chunking
    CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, error)
    CALL h5pset_chunk_f(crp_list, time_rank, dimsc1, error)
    
    CALL h5screate_simple_f(time_rank, time_dims, time_filespace, error, maxdims1D)
    CALL h5screate_simple_f(time_rank, time_dims, time_memspace, error, maxdims1D)
    IF(detector_rank.NE.1) &
         CALL h5sselect_none_f(time_memspace,error)
    !Create a dataset using cparms creation properties
    name = '/time'
    CALL h5dcreate_f(file_id, name, H5T_NATIVE_DOUBLE, time_filespace, &
         time_dset_id, error, crp_list )
    CALL h5sclose_f(time_filespace, error)
    
    time_cont = time_dims
    time_offset = 0
    ! Select hyperslab in the file.
    CALL h5dget_space_f(time_dset_id, time_filespace, error)
    CALL h5sselect_hyperslab_f (time_filespace, H5S_SELECT_SET_F, time_offset,&
         time_cont, error)!, stride, blck)
    IF(detector_rank.NE.1) &
         CALL h5sselect_none_f(time_filespace,error)
       
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error) 
    CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
    !Write data array to dataset
    time = 0.0_dp
    CALL h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_dims, &
         error, file_space_id = time_filespace, mem_space_id = time_memspace, &
         xfer_prp = plist_id)

    CALL h5sclose_f(time_filespace, error)
    CALL h5sclose_f(time_memspace, error)
    CALL h5dclose_f(time_dset_id, error)
    ! Close property list
    CALL h5pclose_f(crp_list, error)   
    CALL h5pclose_f(plist_id, error)
    !-----------------------------------!
    
    ! Close the property list, group, file and interface
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
    ! Set the count of datasets to 1. It will start at 1.
    dataset_count = 1
    
  END SUBROUTINE create_detector_file_parallel
#endif
  !***************************************************!
  
  SUBROUTINE write_detector_points_serial(filename, psi, rank, dims, &
       time, efield, afield )
    
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    COMPLEX(dp), INTENT(IN)        :: psi(:)
    INTEGER, INTENT(IN)            :: rank
    INTEGER, INTENT(IN)            :: dims(:)
    REAL(dp), INTENT(IN)           :: time
    REAL(dp), INTENT(IN)           :: efield(:), afield(:)

    COMPLEX(dp), ALLOCATABLE       :: wave2D(:, :)
    COMPLEX(dp), ALLOCATABLE       :: wave3D(:, :, :)
    INTEGER                        :: wave_rank
    INTEGER                        :: field_rank
    INTEGER                        :: time_rank
    REAL(dp)                       :: cwave(2)
    REAL(dp)                       :: field(2,3)
    ! File identifier
    INTEGER(HID_T)                 :: file_id
    ! Property list identifier 
    INTEGER(HID_T)                 :: cpr_id
    ! Group identifiers
    INTEGER(HID_T)                 :: wave_group
    ! Dataspace identifier
    INTEGER(HID_T)                 :: wave_dspace_id
    INTEGER(HID_T)                 :: field_dspace_id
    INTEGER(HID_T)                 :: time_dspace_id  
    ! Dataset identifier  
    INTEGER(HID_T)                 :: wave_dset_id
    INTEGER(HID_T)                 :: field_dset_id
    INTEGER(HID_T)                 :: time_dset_id
    ! Dataset dimensions 
    INTEGER(HSIZE_T)               :: wave_dims(2)
    INTEGER(HSIZE_T)               :: field_dims(3)
    INTEGER(HSIZE_T)               :: time_dims(1)
    INTEGER(HSSIZE_T)              :: wave_offset(2)
    INTEGER(HSSIZE_T)              :: field_offset(3)
    INTEGER(HSSIZE_T)              :: time_offset(1)
    INTEGER(HSIZE_T)               :: wave_count(2)
    INTEGER(HSIZE_T)               :: field_count(3)
    INTEGER(HSIZE_T)               :: time_count(1)
    ! Memspace
    INTEGER(HID_T)                 :: wave_memspace
    INTEGER(HID_T)                 :: field_memspace
    INTEGER(HID_T)                 :: time_memspace
    ! Error flag
    INTEGER                        :: error, ipts
    CHARACTER(LEN=6)               :: cstep, cpts
    CHARACTER(LEN=100)             :: name
    CHARACTER(LEN=100)             :: setname
    
    !-------------------------------------------------!
    
    ! Increase in 1 the count 
    dataset_count = dataset_count + 1
    WRITE(cstep,'(I6.6)') dataset_count

    
    IF(rank.EQ.2) THEN
       ALLOCATE(wave2D(1:dims(1),1:dims(2)))
       wave2D = RESHAPE(psi,SHAPE(wave2D))
    ELSEIF(rank.EQ.3) THEN
       ALLOCATE(wave3D(1:dims(1),1:dims(2),1:dims(3)))
       wave3D = RESHAPE(psi,SHAPE(wave3D))
    ELSE
       WRITE(*,*) 'Wavefunction dimensions not implemented '&
            &'in write detector points.'
       STOP
    ENDIF
    
    ! Initialize fortran interface
    CALL h5open_f(error)
    
    ! Open files and groups
    name = TRIM(filename) // '.h5'
    CALL h5fopen_f(name, H5F_ACC_RDWR_F, file_id, error)

    !-------------- Write the detector datapoints value
    name = '/wavefunctions'
    CALL h5gopen_f(file_id, name, wave_group, error)
    
    wave_rank   = 2
    wave_dims   = (/dataset_count,2/)
    wave_offset = (/dataset_count -1,0/)
    wave_count  = (/1,2/)
    DO ipts = 1, numdetectorpts
       IF(i_am_in_local(ipts).EQ.1) THEN
          IF(rank.EQ.2) THEN
             cwave = &
                  (/REAL(wave2D(detectorpoints_indx(ipts,1),detectorpoints_indx(ipts,2))),&
                  AIMAG(wave2D(detectorpoints_indx(ipts,1),detectorpoints_indx(ipts,2))) /) 
          ELSEIF(rank.EQ.3) THEN
             cwave = &
                  (/REAL(wave3D(detectorpoints_indx(ipts,1),detectorpoints_indx(ipts,2),&
                  detectorpoints_indx(ipts,3))),&
                  AIMAG(wave3D(detectorpoints_indx(ipts,1),detectorpoints_indx(ipts,2),&
                  detectorpoints_indx(ipts,3)) ) /) 
          ELSE
             WRITE(*,*) 'Rank not implemented in write_detector_file.'
             STOP
          ENDIF
       ELSE
          cwave = ZERO
       ENDIF
          ! Open dataset
          WRITE(cpts,'(I6.6)') ipts
          name = '/wavefunctions/' // cpts
          CALL h5dopen_f(file_id, name, wave_dset_id, error)
          !Extend the dataset.
          CALL h5dset_extent_f(wave_dset_id, wave_dims, error)
          
          CALL h5screate_simple_f (wave_rank, wave_count, wave_memspace, error)
          CALL h5dget_space_f(wave_dset_id, wave_dspace_id, error)
          CALL h5sselect_hyperslab_f(wave_dspace_id, H5S_SELECT_SET_F, &
               wave_offset, wave_count, error)
          CALL H5dwrite_f(wave_dset_id, H5T_NATIVE_DOUBLE, cwave, wave_count, error, &
               wave_memspace, wave_dspace_id)
          CALL h5sclose_f(wave_dspace_id, error)
          CALL h5dclose_f(wave_dset_id, error)
    ENDDO
    ! Close group
    CALL h5gclose_f(wave_group, error)
    !--------------------------------------------------!

    !-------------- Write to the field dataset
    field_rank   = 3
    field_dims   = (/dataset_count,2,3/)
    field_offset = (/dataset_count-1,0,0/)
    field_count  = (/1,2,3/)

    field(1,:) = efield
    field(2,:) = afield
    ! Open dataset
    name = '/field'
    CALL h5dopen_f(file_id, name, field_dset_id, error)
    !Extend the dataset.
    CALL h5dset_extent_f(field_dset_id, field_dims, error)
    
    CALL h5screate_simple_f (field_rank, field_count, field_memspace, error)
    CALL h5dget_space_f(field_dset_id, field_dspace_id, error)
    CALL h5sselect_hyperslab_f(field_dspace_id, H5S_SELECT_SET_F, &
         field_offset, field_count, error)
    CALL H5dwrite_f(field_dset_id, H5T_NATIVE_DOUBLE, field, field_count, error, &
         field_memspace, field_dspace_id)
    CALL h5sclose_f(field_dspace_id, error)
    CALL h5dclose_f(field_dset_id, error)
    !--------------------------------------------------!

    !-------------- Write to the time dataset
    time_rank   = 1
    time_dims   = dataset_count
    time_offset = dataset_count-1
    time_count  = 1

    ! Open dataset
    name = '/time'
    CALL h5dopen_f(file_id, name, time_dset_id, error)
    !Extend the dataset.
    CALL h5dset_extent_f(time_dset_id, time_dims, error)
    
    CALL h5screate_simple_f (time_rank, time_count, time_memspace, error)
    CALL h5dget_space_f(time_dset_id, time_dspace_id, error)
    CALL h5sselect_hyperslab_f(time_dspace_id, H5S_SELECT_SET_F, &
         time_offset, time_count, error)
    CALL H5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time, time_count, error, &
         time_memspace, time_dspace_id)
    CALL h5sclose_f(time_dspace_id, error)
    CALL h5dclose_f(time_dset_id, error)
    !--------------------------------------------------!
    
    
    ! Close file and interface
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
    IF(rank.EQ.2) THEN
       DEALLOCATE(wave2D)
    ELSEIF(rank.EQ.3) THEN
       DEALLOCATE(wave3D)
    ENDIF
    
  END SUBROUTINE write_detector_points_serial
  
  !**************************************************************************!
!!$#if _COM_MPI
!!$  SUBROUTINE write_detector_files_parallel()
!!$
!!$    IMPLICIT NONE
!!$
!!$
!!$
!!$
!!$  END SUBROUTINE write_detector_files_parallel
!!$#endif
  !**************************************************************************!

  SUBROUTINE get_detector_params(filename,nodetectorpts,timesteps, &
       detectorpts_rank)
    
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(OUT)           :: nodetectorpts
    INTEGER, INTENT(OUT)           :: timesteps
    INTEGER, INTENT(OUT)           :: detectorpts_rank
 
    CHARACTER(LEN=100)             :: name
    INTEGER(HID_T)                 :: file_id
    INTEGER(HID_T)                 :: dset_id
    INTEGER(HID_T)                 :: dataspace
    INTEGER(HSIZE_T)               :: dimsr(2), maxdimsr(2) 
    INTEGER                        :: error, rankr
    
    !---------------------------------------------------------!

    ! Initlialize interface
    CALL h5open_f(error)
    
    ! Open file (Read only)
    name = filename // '.h5'
    CALL h5fopen_f (name, H5F_ACC_RDONLY_F, file_id, error)

    !Open the  dataset.
    name = '/field'
    CALL h5dopen_f(file_id, name, dset_id, error)
    
    !Get dataset's dataspace handle for time dataset.
    CALL h5dget_space_f(dset_id, dataspace, error)

    !Get dataspace's dimensions.
    CALL h5sget_simple_extent_dims_f(dataspace, dimsr, maxdimsr, error)

    timesteps = dimsr(1)
    ! Close dataspce and dataset
    CALL h5sclose_f(dataspace, error)
    CALL h5dclose_f(dset_id, error)

    ! To know how many detector points
    !Open the  dataset.
    name = '/detector_points'
    CALL h5dopen_f(file_id, name, dset_id, error)
    
    !Get dataset's dataspace handle for time dataset.
    CALL h5dget_space_f(dset_id, dataspace, error)

    !Get dataspace's dimensions.
    CALL h5sget_simple_extent_dims_f(dataspace, dimsr, maxdimsr, error)

    nodetectorpts = dimsr(1)
    detectorpts_rank = dimsr(2)
    
    ! Close dataspce and dataset
    CALL h5sclose_f(dataspace, error)
    CALL h5dclose_f(dset_id, error)
    
    ! Close file and interface
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
  END SUBROUTINE get_detector_params

  !**************************************************************************!
  
  SUBROUTINE read_detector(filename,points,wavefs,time,field)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    REAL(dp), INTENT(OUT)          :: points(:, :)
    COMPLEX(dp), INTENT(OUT)       :: wavefs(:, :)
    REAL(dp), INTENT(OUT)          :: time(:)
    REAL(dp), INTENT(OUT)          :: field(:, :, :)

    INTEGER                        :: numpts, timesteps
    REAL(dp), ALLOCATABLE          :: cwave(:, :)
    CHARACTER(LEN=100)             :: name
    INTEGER(HID_T)                 :: file_id
    INTEGER(HID_T)                 :: wave_group
    INTEGER(HID_T)                 :: dset_id
    INTEGER(HID_T)                 :: dspace
    INTEGER(HID_T)                 :: memspace
    INTEGER(HSIZE_T)               :: dimsr1(1), maxdimsr1(1)
    INTEGER(HSIZE_T)               :: dimsr2(2), maxdimsr2(2)
    INTEGER(HSIZE_T)               :: dimsr3(3), maxdimsr3(3) 
    INTEGER                        :: error, rankr, ipts
    CHARACTER(LEN=6)               :: cpts
    
    !--------------------------------------------------!
    
    ! Initlialize interface
    CALL h5open_f(error)
    
    !Open the file.
    name = filename // '.h5'
    CALL h5fopen_f (name, H5F_ACC_RDONLY_F, file_id, error)
    
    !--------- Read detector points -----------------!
    !Open the  dataset.
    name = '/detector_points'
    CALL h5dopen_f(file_id, name, dset_id, error)
    
    !Get dataset's dataspace handle.
    CALL h5dget_space_f(dset_id, dspace, error)
    
    !Get dataspace's rank.
    CALL h5sget_simple_extent_ndims_f(dspace, rankr, error)
    
    !Get dataspace's dimensions.
    CALL h5sget_simple_extent_dims_f(dspace, dimsr2, maxdimsr2, error)

    ! Get number of detector points
    numpts = dimsr2(1)
    ! Create memory dataspace
    CALL h5screate_simple_f(rankr, dimsr2, memspace, error)
    
    !Read data 
    points = 0.0_dp
    CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, points, dimsr2, &
         error, memspace, dspace)
    
    !Close the objects that were opened.
    CALL h5sclose_f(dspace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    
    !-------------------------------------------------------------!
    
    !--------- Read time -----------------!
    !Open the  dataset.
    name = '/time'
    CALL h5dopen_f(file_id, name, dset_id, error)
    
    !Get dataset's dataspace handle.
    CALL h5dget_space_f(dset_id, dspace, error)
    
    !Get dataspace's rank.
    CALL h5sget_simple_extent_ndims_f(dspace, rankr, error)
    
    !Get dataspace's dimensions.
    CALL h5sget_simple_extent_dims_f(dspace, dimsr1, maxdimsr1, error)
    
    timesteps = dimsr1(1)
    ! Create memory dataspace
    CALL h5screate_simple_f(rankr, dimsr1, memspace, error)
    
    !Read data 
    time = 0.0_dp
    CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, time, dimsr1, &
         error, memspace, dspace)
    
    !Close the objects that were opened.
    CALL h5sclose_f(dspace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    
    !-------------------------------------------------------------!
    
    ! -------- Read detector points -----------------!
    name = '/wavefunctions'
    CALL h5gopen_f(file_id, name, wave_group, error)
    wavefs = ZIMAGONE
    ! Allocate buffer array for read detector data.
    ALLOCATE(cwave(timesteps,2))
    
    !Open the  dataset.
    DO ipts = 1, numpts
       WRITE(cpts,'(I6.6)') ipts
       name = '/wavefunctions/' //cpts
       CALL h5dopen_f(file_id, name, dset_id, error)
       
       !Get dataset's dataspace handle.
       CALL h5dget_space_f(dset_id, dspace, error)
       
       !Get dataspace's rank.
       CALL h5sget_simple_extent_ndims_f(dspace, rankr, error)
       
       !Get dataspace's dimensions.
       CALL h5sget_simple_extent_dims_f(dspace, dimsr2, maxdimsr2, error)
       
       ! Create memory dataspace
       CALL h5screate_simple_f(rankr, dimsr2, memspace, error)
       
       !Read data 
       CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, cwave, dimsr2, &
            error, memspace, dspace)
       
       wavefs(ipts,:) = CMPLX(cwave(:,1),cwave(:,2))
       cwave = 0.0_dp
       
       !Close the objects that were opened.
       CALL h5sclose_f(dspace, error)
       CALL h5sclose_f(memspace, error)
       CALL h5dclose_f(dset_id, error)
    ENDDO
    
    DEALLOCATE(cwave)
    CALL h5gclose_f(wave_group,error)
    
    !-------------------------------------------------------------!
    
    ! -------- Read detector points -----------------!
    !Open the  dataset.
    name = '/field'
    CALL h5dopen_f(file_id, name, dset_id, error)
    
    !Get dataset's dataspace handle.
    CALL h5dget_space_f(dset_id, dspace, error)
    
    !Get dataspace's rank.
    CALL h5sget_simple_extent_ndims_f(dspace, rankr, error)
    
    !Get dataspace's dimensions.
    CALL h5sget_simple_extent_dims_f(dspace, dimsr3, maxdimsr3, error)
    
    ! Create memory dataspace
    CALL h5screate_simple_f(rankr, dimsr3, memspace, error)
    
    !Read data 
    field = 0.0_dp
    CALL H5dread_f(dset_id, H5T_NATIVE_DOUBLE, field, dimsr3, &
         error, memspace, dspace)
    
    !Close the objects that were opened.
    CALL h5sclose_f(dspace, error)
    CALL h5sclose_f(memspace, error)
    CALL h5dclose_f(dset_id, error)
    
    !-------------------------------------------------------------!
    
    ! Close file and interface
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
    
  END SUBROUTINE read_detector
  
  !**************************************************************************!

  SUBROUTINE get_solid_angles(detectorpoints,solid_angles)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)   :: detectorpoints(:, :)
    REAL(dp), INTENT(OUT)  :: solid_angles(:)

    INTEGER               :: ipt
    INTEGER               :: rank, nopts
    INTEGER               :: dims(2)
    REAL(dp)              :: thetapt1, thetapt2
    REAL(dp)              :: phipt1, phipt2
    REAL(dp)              :: dtheta, dphi, rpt
    !----------------------------------------------!

    solid_angles = 0.0_dp
    ! Get dimensions for detector
    ! (number of points and rank).
    dims = SHAPE(detectorpoints)
    nopts = dims(1)
    rank  = dims(2)

    thetapt1 = 0.0_dp
    phipt1 = 0.0_dp
    DO ipt = 1, nopts
       rpt = SQRT(detectorpoints(ipt,1)**2 + &
            detectorpoints(ipt,2)**2)
       IF(rpt.EQ.0.0_dp) THEN
          thetapt2 = pi / 2.0_dp
       ELSE
          thetapt2 = ACOS(detectorpoints(ipt,2)/rpt)
       ENDIF
       
       solid_angles(ipt) = -COS(thetapt2) + COS(thetapt1)
       dtheta = thetapt2 - thetapt1
       solid_angles(ipt) = solid_angles(ipt) * dtheta
       thetapt1 = thetapt2
    ENDDO
    IF(rank.EQ.2) THEN
       solid_angles = solid_angles * twopi
    ELSEIF(rank.EQ.3) THEN
       phipt2 = pi - ATAN2(detectorpoints(ipt,2),-detectorpoints(ipt,1))
       dphi = phipt2 - phipt1
       solid_angles = solid_angles * dphi
       phipt1 = phipt2
    ELSE
       WRITE(*,*) 'No solid angle for this geometry.'
       STOP
    ENDIF

  END SUBROUTINE get_solid_angles
  
  !**************************************************************************!

  SUBROUTINE calculate_sampt_pes(wavefs,detectorpts,time, &
       pes,w,pad)

    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)          :: wavefs(:, :)
    REAL(dp), INTENT(IN)             :: detectorpts(:, :)
    REAL(dp), INTENT(IN)             :: time(:)
    REAL(dp), INTENT(OUT)            :: pes(:)
    REAL(dp), INTENT(OUT)            :: w(:)
    REAL(dp), INTENT(OUT), OPTIONAL  :: pad(:, :)

    INTEGER                          :: dims(2)
    COMPLEX(dp), ALLOCATABLE         :: trans(:)
    REAL(dp), ALLOCATABLE            :: modsq(:)
    REAL(dp), ALLOCATABLE            :: solid_angle(:)
    REAL(dp)                         :: dt, dw
    INTEGER                          :: numwpts, rank
    INTEGER                          :: nodetectorpts
    INTEGER                          :: timesteps
    INTEGER                          :: ipt, itime
    
    !------------------------------------------------------!

    pes = 0.0_dp
    pad = 0.0_dp
    
    dims = SHAPE(detectorpts)
    nodetectorpts = dims(1)
    rank = dims(2)
    timesteps = SIZE(time)
    dt = time(2) - time(1)
    
    ! Create energy array
    dw = twopi / (timesteps * dt)
    numwpts = INT((timesteps-1)*0.5+1)
    w = 0.0_dp
    DO ipt = 1, numwpts
       w(ipt) = REAL(ipt-1,dp) * dw
    ENDDO
    
    ALLOCATE(trans(timesteps))
    ALLOCATE(modsq(timesteps))
    
    ! Get solid angles for weighting
    ALLOCATE(solid_angle(nodetectorpts))
    CALL get_solid_angles(detectorpts,solid_angle)
    
    DO ipt = 1, nodetectorpts
       CALL FourierTransform(wavefs(ipt,:),trans,1,(/timesteps/))
       trans = trans * dt / SQRT(twopi)
       DO itime = 1, numwpts
          modsq(itime) = REAL(CONJG(trans(itime)) * &
               trans(itime), dp)
          pes(itime) = pes(itime) + solid_angle(ipt) * &
               SQRT(w(itime)) * modsq(itime)
          IF(PRESENT(pad)) &
               pad(ipt,itime) = modsq(itime)
       ENDDO
    ENDDO

    ! Release memory
    DEALLOCATE(trans)
    DEALLOCATE(modsq)
    DEALLOCATE(solid_angle)
    
  END SUBROUTINE calculate_sampt_pes
  
END MODULE io_sampt
