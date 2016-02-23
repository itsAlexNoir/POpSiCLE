!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: surface
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 21/10/2015
!
! DESCRIPTION:
!> \brief Surface module
!> \details This module contains routines to get the surface flux
!> from a wavefunction expanded in a spherical harmonic basis.
!
!------------------------------------------------------------------------------
MODULE surface
  
  USE constants
  USE tools
  USE boundary
  USE io_surface
#if _COM_MPI
  USE MPI
#endif
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: initialize_cylindrical_surface
  PUBLIC       :: initialize_cartesian_surface
  PUBLIC       :: get_cylindrical_surface
  PUBLIC       :: get_cartesian_surface
  PUBLIC       :: delete_surface2D
  PUBLIC       :: delete_surface3D
  
  !-------------------------------------------!
  
  INTERFACE initialize_cylindrical_surface
     MODULE PROCEDURE initialize_cylindrical_surface_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cylindrical_surface_parallel
#endif
  END INTERFACE initialize_cylindrical_surface
  
  INTERFACE initialize_cartesian_surface
     MODULE PROCEDURE initialize_cartesian3D_surface_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cartesian3D_surface_parallel
#endif
  END INTERFACE initialize_cartesian_surface
  
  INTERFACE get_cylindrical_surface
     MODULE PROCEDURE get_cylindrical_surface_serial
#if _COM_MPI
     MODULE PROCEDURE get_cylindrical_surface_parallel
#endif
  END INTERFACE get_cylindrical_surface
  
  INTERFACE get_cartesian_surface
     MODULE PROCEDURE get_cartesian3D_surface_serial
#if _COM_MPI
     MODULE PROCEDURE get_cartesian3D_surface_parallel
#endif
  END INTERFACE get_cartesian_surface
  
  INTERFACE delete_surface2D
     MODULE PROCEDURE delete_surface2D_serial
#if _COM_MPI
     MODULE PROCEDURE delete_surface2D_parallel
#endif
  END INTERFACE delete_surface2D
  
  INTERFACE delete_surface3D
     MODULE PROCEDURE delete_surface3D_serial
#if _COM_MPI
     MODULE PROCEDURE delete_surface3D_parallel
#endif
  END INTERFACE delete_surface3D
  
  !-------------------------------------------!
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_local(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_global(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dr(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dtheta(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_deriv(:)
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_local(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_global(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dtheta(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dphi(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_deriv(:, :)
  
  COMPLEX(dp), ALLOCATABLE  :: psi_l(:, :)
  COMPLEX(dp), ALLOCATABLE  :: psi_lm(:, :, :)
  
  REAL(dp)                  :: rb, deltar
  INTEGER                   :: numpts, numrpts
  INTEGER                   :: numthetapts, numphipts
  INTEGER                   :: numthetaptsperproc
  INTEGER                   :: numphiptsperproc
  INTEGER                   :: surfacecomm
  INTEGER                   :: surfacerank
  INTEGER                   :: numsurfaceprocs
  
CONTAINS
  !********************************************************************!
  !********************************************************************!
  SUBROUTINE initialize_cylindrical_surface_serial(rho_ax, z_ax, dims, &
       Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       write_to_file, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)        :: dims(:)
    REAL(dp), INTENT(IN)       :: Rboundary
    REAL(dp), INTENT(IN)       :: radius_tolerance
    INTEGER, INTENT(IN)        :: fd_rule
    REAL(dp), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)        :: lmax
    LOGICAL, INTENT(IN)        :: write_to_file
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    !------------------------------------------------------------!
    
    ! Initialize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    !! Initialize spherical grid
    !! For scattered interpolation uncomment the paragraph below.
    !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
    !     Rboundary, radius_tolerance, fd_rule, dr, lmax, numpts, &
    !     numrpts, numthetapts )
    
    !For bicubic interpolation, the paragraph below.
    CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
         fd_rule, dr, lmax, numrpts, numthetapts)
    
    
    ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))   
    ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
    ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))
    ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
    
    spherical_wave2D        = ZERO
    spherical_wave2D_dr     = ZERO
    spherical_wave2D_dtheta = ZERO
    spherical_wave2D_deriv  = ZERO
    
    rb = Rboundary
    deltar = dr
    
    IF(write_to_file) &
         CALL create_surface_file(filename)
    
  END SUBROUTINE initialize_cylindrical_surface_serial
  
  !****************************************************************!
#if _COM_MPI
  SUBROUTINE initialize_cylindrical_surface_parallel(rho_ax, z_ax, dims, &
       Rboundary, radius_tolerance, fd_rule, dr, lmax, rank, size, comm, &
       write_to_file, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)        :: dims(:)
    REAL(dp), INTENT(IN)       :: Rboundary
    REAL(dp), INTENT(IN)       :: radius_tolerance
    INTEGER, INTENT(IN)        :: fd_rule
    REAL(dp), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)        :: lmax
    INTEGER, INTENT(IN)        :: rank, size
    INTEGER, INTENT(IN)        :: comm
    LOGICAL, INTENT(IN)        :: write_to_file
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER                    :: numthetaoffset
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    ! Initialize spherical grid
    ! For scattered interpolation
    !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
    !     radius_tolerance, fd_rule, dr, lmax, rank, size, comm, &
    !     numpts, numrpts, numthetapts, surfacerank, numsurfaceprocs, &
    !     surfacecomm, numthetaptsperproc )
    
    ! For bicubic interpolation
    CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
         fd_rule, dr, lmax, rank, size, comm, &
         numrpts, numthetapts, surfacerank, numsurfaceprocs, &
         surfacecomm, numthetaptsperproc )

    numthetaoffset = INT(numthetapts / numsurfaceprocs)
    
    IF(i_am_surface(rank).EQ.1) THEN
       ALLOCATE(spherical_wave2D_local(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_global(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))
       IF(numthetaoffset .EQ. 0) THEN
          ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
       ELSE
          ALLOCATE(spherical_wave2D(1:numrpts,1:numthetaptsperproc))
          ALLOCATE(spherical_wave2D_deriv(1:numthetaptsperproc))
       ENDIF
       
       spherical_wave2D_local  = ZERO
       spherical_wave2D_global = ZERO
       spherical_wave2D_dr     = ZERO
       spherical_wave2D_dtheta = ZERO
       spherical_wave2D        = ZERO
       spherical_wave2D_deriv  = ZERO
       
       rb = Rboundary
       deltar = dr
       
       IF(write_to_file) &
            CALL create_surface_file(filename, surfacecomm)
    ENDIF
    
  END SUBROUTINE initialize_cylindrical_surface_parallel
#endif  
  !****************************************************************!
  
  SUBROUTINE initialize_cartesian3D_surface_serial(x_ax, y_ax, z_ax, &
       dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       write_to_file, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)        :: dims(:)
    REAL(dp), INTENT(IN)       :: Rboundary
    REAL(dp), INTENT(IN)       :: radius_tolerance
    INTEGER, INTENT(IN)        :: fd_rule
    REAL(dp), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)        :: lmax
    LOGICAL, INTENT(IN)        :: write_to_file
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    ! Initialize spherical grid
    ! For scattered interpolation
    !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
    !     Rboundary, radius_tolerance, fd_rule, dr, lmax, numpts, &
    !     numrpts, numthetapts, numphipts )

    ! Tricubic interpolation
    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
         Rboundary, fd_rule, dr, lmax, &
         numrpts, numthetapts, numphipts )
    
    ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_deriv(1:numthetapts,1:numphipts))
    
    spherical_wave3D        = ZERO
    spherical_wave3D_dr     = ZERO
    spherical_wave3D_dtheta = ZERO
    spherical_wave3D_dphi   = ZERO
    spherical_wave3D_deriv  = ZERO
    
    rb = Rboundary
    deltar = dr
    
    IF(write_to_file) &
         CALL create_surface_file(filename)
    
  END SUBROUTINE initialize_cartesian3D_surface_serial
  
  !****************************************************************!
#if _COM_MPI
  SUBROUTINE initialize_cartesian3D_surface_parallel(x_ax, y_ax, z_ax, &
       dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, rank, &
       size, comm, write_to_file, filename)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)        :: dims(:)
    REAL(dp), INTENT(IN)       :: Rboundary
    REAL(dp), INTENT(IN)       :: radius_tolerance
    INTEGER, INTENT(IN)        :: fd_rule
    REAL(dp), INTENT(IN)       :: dr
    INTEGER, INTENT(IN)        :: lmax
    INTEGER, INTENT(IN)        :: rank, size
    INTEGER, INTENT(IN)        :: comm
    LOGICAL, INTENT(IN)        :: write_to_file
    CHARACTER(LEN=*), INTENT(IN) :: filename
    
    INTEGER                    :: numthetaoffset, numphioffset
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    ! Initialize spherical grid
    ! For scattered interpolation
    !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
    !     Rboundary, radius_tolerance, fd_rule, dr, lmax, &
    !     rank, size, comm, numpts, numrpts, numthetapts, numphipts, &
    !     surfacerank, numsurfaceprocs, surfacecomm, &
    !     numthetaptsperproc, numphiptsperproc )
    
    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
         Rboundary, fd_rule, dr, lmax, &
         rank, size, comm, numrpts, numthetapts, numphipts, &
         surfacerank, numsurfaceprocs, surfacecomm, &
         numthetaptsperproc, numphiptsperproc )
    
    numthetaoffset = INT(numthetapts / numsurfaceprocs)
    numphioffset = INT(numphipts / numsurfaceprocs) 
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       ALLOCATE(spherical_wave3D_local(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_global(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
       IF((numthetaoffset.EQ.0) .OR. (numphioffset.EQ.0)) THEN
          ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,&
               1:numphipts))
          ALLOCATE(spherical_wave3D_deriv(1:numthetapts,&
            1:numphipts))
       ELSE
          ALLOCATE(spherical_wave3D(1:numrpts,1:numthetaptsperproc,&
               1:numphiptsperproc))
          ALLOCATE(spherical_wave3D_deriv(1:numthetaptsperproc,&
               1:numphiptsperproc))
       ENDIF
       spherical_wave3D_local  = ZERO
       spherical_wave3D_global = ZERO
       spherical_wave3D_dr     = ZERO
       spherical_wave3D_dtheta = ZERO
       spherical_wave3D_dphi   = ZERO
       spherical_wave3D        = ZERO
       spherical_wave3D_deriv  = ZERO
       
       rb = Rboundary
       deltar = dr
       
       IF(write_to_file) &
            CALL create_surface_file(filename, surfacecomm)
    ENDIF
    
  END SUBROUTINE initialize_cartesian3D_surface_parallel
#endif  
  !****************************************************************!
  !****************************************************************!
  
  SUBROUTINE get_cylindrical_surface_serial(filename, wavefunc, &
       rho_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, write_to_file)
    
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :)
    REAL(dp), INTENT(IN)         :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax
    LOGICAL, INTENT(IN)          :: write_to_file
    
    INTEGER                      :: middle_pt
    !----------------------------------------------------------!
    
    ! For scattered interpolation
    !CALL get_cylindrical_boundary( wavefunc, spherical_wave2D, &
    !     spherical_wave2D_dr, spherical_wave2D_dtheta, &
    !     'quadratic')
    
    ! For bicubic interpolation
    CALL get_cylindrical_boundary( rho_ax, z_ax, dims, &
         wavefunc, fd_rule, spherical_wave2D, &
         spherical_wave2D_dr, spherical_wave2D_dtheta)
    
    middle_pt = fd_rule + 1
    CALL make_wave_boundary_derivative(spherical_wave2D,&
         spherical_wave2D_deriv,fd_rule,deltar,&
         numrpts,numthetapts)
    
    ! Write boundary points to a HDF5 file
    IF (write_to_file) &
         CALL write_surface_file(filename, spherical_wave2D(middle_pt,:), &
         spherical_wave2D_deriv, time, efield, afield, lmax)
        
  END SUBROUTINE get_cylindrical_surface_serial

  !**********************************************************************!
#if _COM_MPI
  SUBROUTINE get_cylindrical_surface_parallel(filename, wavefunc, &
       rho_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, rank, write_to_file)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :)
    REAL(dp), INTENT(IN)         :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax, rank
    LOGICAL, INTENT(IN)          :: write_to_file
    
    INTEGER                      :: numtotalpts
    INTEGER                      :: numthetaoffset
    INTEGER                      :: middle_pt, offset
    INTEGER                      :: ir, itheta, ierror
    
    !----------------------------------------------------------!
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       !! For scattered interpolation
       !CALL get_cylindrical_boundary( wavefunc, spherical_wave2D_local, &
       !     spherical_wave2D_dr, spherical_wave2D_dtheta, &
       !     'quadratic', rank)
       
       CALL get_cylindrical_boundary( rho_ax, z_ax, dims, &
            wavefunc, fd_rule, spherical_wave2D_local, &
            spherical_wave2D_dr, spherical_wave2D_dtheta, rank)
       
       ! Communicate the spherical wavefunction
       numtotalpts = numrpts * numthetapts
       spherical_wave2D_global = ZERO
       CALL MPI_ALLREDUCE(spherical_wave2D_local, spherical_wave2D_global, &
            numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
       
       middle_pt = fd_rule + 1
       numthetaoffset = INT(numthetapts / numsurfaceprocs)

       IF(numthetaoffset.EQ.0) THEN
          CALL make_wave_boundary_derivative(spherical_wave2D,&
               spherical_wave2D_deriv,fd_rule,deltar,&
               numrpts,numthetapts)
          
          ! Write boundary points to a HDF5 file
          IF (write_to_file .AND. (surfacerank.EQ.0)) &
               CALL write_surface_file(filename, spherical_wave2D(middle_pt,:), &
               spherical_wave2D_deriv, time, efield, afield, lmax)
          
       ELSE
          offset = numthetaoffset * surfacerank       
          
          DO itheta = 1, numthetaptsperproc
             DO ir = 1, numrpts
                spherical_wave2D(ir,itheta) = &
                     spherical_wave2D_global(ir,offset+itheta)
             ENDDO
          ENDDO
          
          CALL make_wave_boundary_derivative(spherical_wave2D,&
               spherical_wave2D_deriv,fd_rule,deltar,&
               numrpts,numthetaptsperproc)
       
          ! Write boundary points to a HDF5 file
          IF (write_to_file) &
               CALL write_surface_file(filename, spherical_wave2D(middle_pt,:), &
               spherical_wave2D_deriv, time, efield, afield, lmax, &
               surfacerank, numsurfaceprocs, surfacecomm, numthetaptsperproc)
       ENDIF
    ENDIF
 
  END SUBROUTINE get_cylindrical_surface_parallel
#endif
  !*****************************************************************!
  
  SUBROUTINE get_cartesian3D_surface_serial(filename, wavefunc, &
       x_ax, y_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, write_to_file)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
    REAL(dp), INTENT(IN)         :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax
    LOGICAL, INTENT(IN)          :: write_to_file
    
    INTEGER                      :: middle_pt
    !----------------------------------------------------------!
    
    ! For scattered interpolation
    !CALL get_cartesian_boundary( wavefunc, spherical_wave3D, &
    !     spherical_wave3D_dr, spherical_wave3D_dtheta, &
    !     spherical_wave3D_dphi, 'quadratic')
    
    ! For tricubic interpolation
    CALL get_cartesian_boundary( x_ax, y_ax, z_ax, dims, &
         wavefunc, fd_rule, spherical_wave3D, &
         spherical_wave3D_dr, spherical_wave3D_dtheta, &
         spherical_wave3D_dphi)
    
    middle_pt = fd_rule + 1
    CALL make_wave_boundary_derivative(spherical_wave3D,&
         spherical_wave3D_deriv,fd_rule,deltar,&
         numrpts,numthetapts,numphipts)
    
    ! Write boundary points to a HDF5 file
    IF (write_to_file) &
         CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
         spherical_wave3D_deriv, time, efield, afield, lmax)
    
  END SUBROUTINE get_cartesian3D_surface_serial
  
  !**********************************************************************!
#if _COM_MPI
  SUBROUTINE get_cartesian3D_surface_parallel(filename, wavefunc, &
       x_ax, y_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, rank, write_to_file)
    
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
    REAL(dp), INTENT(IN)         :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax, rank
    LOGICAL, INTENT(IN)          :: write_to_file
    
    INTEGER                      :: numtotalpts
    INTEGER                      :: middle_pt
    INTEGER                      :: thetaoffset, phioffset
    INTEGER                      :: numthetaoffset, numphioffset
    INTEGER                      :: ir, itheta, iphi, ierror
    !----------------------------------------------------------!
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       ! For scattered interpolation
       !CALL get_cartesian_boundary( wavefunc, spherical_wave3D_local, &
       !     spherical_wave3D_dr, spherical_wave3D_dtheta, &
       !     spherical_wave3D_dphi, 'quadratic',rank)

       ! For tricubic interpolation
       CALL get_cartesian_boundary( x_ax, y_ax, z_ax, dims, &
            wavefunc, fd_rule, spherical_wave3D_local, &
            spherical_wave3D_dr, spherical_wave3D_dtheta, &
            spherical_wave3D_dphi, rank)
       
       ! Communicate the spherical wavefunction
       numtotalpts = numrpts * numthetapts * numphipts
       spherical_wave3D_global = ZERO
       CALL MPI_ALLREDUCE(spherical_wave3D_local, spherical_wave3D_global, &
            numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
       
       middle_pt = fd_rule + 1
       numthetaoffset = INT(numthetapts / numsurfaceprocs)
       numphioffset = INT(numphipts / numsurfaceprocs)
       
       IF((numthetaoffset.EQ.0) .OR. (numphioffset.EQ.0)) THEN
          
          CALL make_wave_boundary_derivative(spherical_wave3D,&
               spherical_wave3D_deriv,fd_rule,deltar,&
               numrpts,numthetapts,numphipts)
          
          ! Write boundary points to a HDF5 file
          IF (write_to_file .AND. (surfacerank.EQ.0)) &
               CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
               spherical_wave3D_deriv, time, efield, afield, lmax)
          
       ELSE
          thetaoffset = numthetaoffset * surfacerank       
          phioffset = numphioffset * surfacerank       
          
          DO iphi = 1, numphiptsperproc
             DO itheta = 1, numthetaptsperproc
                DO ir = 1, numrpts
                   spherical_wave3D(ir,itheta, iphi) = &
                        spherical_wave3D_global(ir,thetaoffset+itheta, phioffset+iphi)
                ENDDO
             ENDDO
          ENDDO
          
          CALL make_wave_boundary_derivative(spherical_wave3D,&
               spherical_wave3D_deriv,fd_rule,deltar,&
               numrpts,numthetaptsperproc,numphiptsperproc)
          
          ! Write boundary points to a HDF5 file
          IF (write_to_file) &
               CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
               spherical_wave3D_deriv, time, efield, afield, lmax, &
               surfacerank, numsurfaceprocs, surfacecomm, &
               numthetaptsperproc, numphiptsperproc )
       ENDIF
    ENDIF
    
  END SUBROUTINE get_cartesian3D_surface_parallel
#endif
  
  !**********************************************************************!
  !**********************************************************************!
  
  SUBROUTINE delete_surface2D_serial()
    
    IMPLICIT NONE
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    DEALLOCATE(spherical_wave2D,spherical_wave2D_deriv)
    DEALLOCATE(spherical_wave2D_dr,spherical_wave2D_dtheta)
    
  END SUBROUTINE delete_surface2D_serial
  
  !---------------------------------------------------------------!
#if _COM_MPI
  SUBROUTINE delete_surface2D_parallel(rank)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)        :: rank
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(i_am_surface(rank) .EQ. 1) THEN  
       DEALLOCATE(spherical_wave2D_local,spherical_wave2D_global)
       DEALLOCATE(spherical_wave2D,spherical_wave2D_deriv)
       DEALLOCATE(spherical_wave2D_dr,spherical_wave2D_dtheta)
    ENDIF
    
  END SUBROUTINE delete_surface2D_parallel
#endif
  !----------------------------------------------------------!
  
  SUBROUTINE delete_surface3D_serial()
    
    IMPLICIT NONE

    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
    DEALLOCATE(spherical_wave3D_dr)
    DEALLOCATE(spherical_wave3D_dtheta)
    DEALLOCATE(spherical_wave3D_dphi)
    
  END SUBROUTINE delete_surface3D_serial
  
  !-----------------------------------------------------------!
#if _COM_MPI
  SUBROUTINE delete_surface3D_parallel(rank)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)         :: rank
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       DEALLOCATE(spherical_wave3D_local,spherical_wave3D_global)
       DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
       DEALLOCATE(spherical_wave3D_dr)
       DEALLOCATE(spherical_wave3D_dtheta)
       DEALLOCATE(spherical_wave3D_dphi)
    ENDIF
    
  END SUBROUTINE delete_surface3D_parallel
#endif
  
END MODULE surface
