!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: corrsurface
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 17/03/2016
!
! DESCRIPTION:
!> \brief Correlated surface module
!> \details This module contains routines to get the surface flux
!> from the electronic part of a wavefunction expanded in
!> a spherical harmonic basis.
!
!------------------------------------------------------------------------------
MODULE corrsurface
  
  USE constants_pop
  USE tools
  USE corrtools
  USE boundary
  USE io_surface
  USE io_corrsurface
#if _COM_MPI
  USE MPI
#endif
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: initialize_correlated_cylindrical_surface
  !PUBLIC       :: initialize_correlated_cartesian_surface
  PUBLIC       :: get_correlated_cylindrical_surface
  !PUBLIC       :: get_correlated_cartesian_surface
  PUBLIC       :: delete_correlated_surface2D
  !PUBLIC       :: delete_correlated_surface3D
  
  !-------------------------------------------!
  
!!$  INTERFACE initialize_correlated_cartesian_surface
!!$     MODULE PROCEDURE initialize_correlated_cartesian3D_surface_serial
!!$#if _COM_MPI
!!$     MODULE PROCEDURE initialize_correlated_cartesian3D_surface_parallel
!!$#endif
!!$  END INTERFACE initialize_correlated_cartesian_surface
    
!!$  INTERFACE get_correlated_cartesian_surface
!!$     MODULE PROCEDURE get_correlated_cartesian3D_surface_serial
!!$#if _COM_MPI
!!$     MODULE PROCEDURE get_correlated_cartesian3D_surface_parallel
!!$#endif
!!$  END INTERFACE get_correlated_cartesian_surface
  
  
!!$  INTERFACE delete_correlated_surface3D
!!$     MODULE PROCEDURE delete_correlated_surface3D_serial
!!$#if _COM_MPI
!!$     MODULE PROCEDURE delete_correlated_surface3D_parallel
!!$#endif
!!$  END INTERFACE delete_correlated_surface3D
  
  !-------------------------------------------!
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D1D_local(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D1D_global(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D1D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D1D_dr(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D1D_dtheta(:, :)
  COMPLEX(dp), ALLOCATABLE  :: sphericaln_wave2D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: sphericaln_wave2D_deriv(:, :)
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_local(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_global(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dtheta(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dphi(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_deriv(:, :)
  
  COMPLEX(dp), ALLOCATABLE  :: psi_l(:, :)
  COMPLEX(dp), ALLOCATABLE  :: psi_lm(:, :, :)
  REAL(dp), ALLOCATABLE     :: basis_set(:, :)
  
  REAL(dp)                  :: rb, deltar
  INTEGER                   :: basis_order
  INTEGER                   :: numbigrpts
  
CONTAINS
  !********************************************************************!
  !********************************************************************!
  
  SUBROUTINE initialize_correlated_cylindrical_surface(rho_ax, z_ax,&
       bigr_ax, dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       order_basis, filename, mpi_rank, mpi_size, comm)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)           :: bigr_ax(:)
    INTEGER, INTENT(IN)            :: dims(:)
    REAL(dp), INTENT(IN)           :: Rboundary
    REAL(dp), INTENT(IN)           :: radius_tolerance
    INTEGER, INTENT(IN)            :: fd_rule
    REAL(dp), INTENT(IN)           :: dr
    INTEGER, INTENT(IN)            :: lmax
    INTEGER, INTENT(IN)            :: order_basis
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(IN), OPTIONAL  :: mpi_rank, mpi_size
    INTEGER, INTENT(IN), OPTIONAL  :: comm
    
    INTEGER                        :: numthetaoffset
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    rb = Rboundary
    deltar = dr
    basis_order = order_basis
    numbigrpts = SIZE(bigr_ax)
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       ! Initialize spherical grid
       ! For scattered interpolation
       !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
       !     radius_tolerance, fd_rule, dr, lmax, rank, size, comm )
       
       ! For bicubic interpolation
       CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
            fd_rule, dr, lmax, mpi_rank, mpi_size, comm )
       
       numthetaoffset = INT(numthetapts / numsurfaceprocs)
       
       IF(i_am_surface(mpi_rank).EQ.1) THEN
          
          
          ALLOCATE(spherical_wave2D1D_local(1:numbigrpts,1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D1D_global(1:numbigrpts,1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D1D_dr(1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D1D_dtheta(1:numrpts,1:numthetapts))
          IF(numthetaoffset .EQ. 0) THEN
             ALLOCATE(sphericaln_wave2D(1:basis_order,1:numrpts,1:numthetapts))
             ALLOCATE(sphericaln_wave2D_deriv(1:basis_order,1:numthetapts))
          ELSE
             ALLOCATE(sphericaln_wave2D(1:basis_order,1:numrpts,1:numthetaptsperproc))
             ALLOCATE(sphericaln_wave2D_deriv(1:basis_order,1:numthetaptsperproc))
          ENDIF
          ALLOCATE(basis_set(1:numbigrpts,1:basis_order))
          
          spherical_wave2D1D_local  = ZERO
          spherical_wave2D1D_global = ZERO
          
          ! Create basis set
          CALL create_sine_basis(basis_set,basis_order, bigr_ax)
          
          CALL create_surface_file(filename, surfacecomm)
       ENDIF
    ELSE

       !! Initialize spherical grid
       !! For scattered interpolation uncomment the paragraph below.
       !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
       !     Rboundary, radius_tolerance, fd_rule, dr, lmax )
       
       !For bicubic interpolation, the paragraph below.
       CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
            fd_rule, dr, lmax )
       
       ! Allocate arrays
       ALLOCATE(spherical_wave2D1D(1:numbigrpts,1:numrpts,1:numthetapts))   
       ALLOCATE(spherical_wave2D1D_dr(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D1D_dtheta(1:numrpts,1:numthetapts))
       ALLOCATE(sphericaln_wave2D(1:basis_order,1:numrpts,1:numthetapts))   
       ALLOCATE(sphericaln_wave2D_deriv(1:basis_order,1:numthetapts))
       ALLOCATE(basis_set(1:numbigrpts,1:basis_order))
       
       ! Create basis set
       CALL create_sine_basis(basis_set,order_basis, bigr_ax)
       
       CALL create_surface_file(filename)
    ENDIF

    spherical_wave2D1D        = ZERO
    spherical_wave2D1D_dr     = ZERO
    spherical_wave2D1D_dtheta = ZERO
    sphericaln_wave2D         = ZERO
    sphericaln_wave2D_deriv   = ZERO
    basis_set                 = 0.0_dp
    
    
       
  END SUBROUTINE initialize_correlated_cylindrical_surface

  !****************************************************************!
  
!!$  SUBROUTINE initialize_cartesian3D_surface_serial(x_ax, y_ax, z_ax, &
!!$       dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, &
!!$       write_to_file, filename)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    REAL(dp), INTENT(IN)       :: x_ax(:), y_ax(:), z_ax(:)
!!$    INTEGER, INTENT(IN)        :: dims(:)
!!$    REAL(dp), INTENT(IN)       :: Rboundary
!!$    REAL(dp), INTENT(IN)       :: radius_tolerance
!!$    INTEGER, INTENT(IN)        :: fd_rule
!!$    REAL(dp), INTENT(IN)       :: dr
!!$    INTEGER, INTENT(IN)        :: lmax
!!$    LOGICAL, INTENT(IN)        :: write_to_file
!!$    CHARACTER(LEN=*), INTENT(IN) :: filename
!!$    
!!$    !------------------------------------------------------------!
!!$    
!!$    ! Initilize finite-difference coefficients
!!$    CALL initialize_fd_coeffs(fd_rule)
!!$    
!!$    ! Initialize spherical grid
!!$    ! For scattered interpolation
!!$    !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
!!$    !     Rboundary, radius_tolerance, fd_rule, dr, lmax, numpts, &
!!$    !     numrpts, numthetapts, numphipts )
!!$
!!$    ! Tricubic interpolation
!!$    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
!!$         Rboundary, fd_rule, dr, lmax, &
!!$         numrpts, numthetapts, numphipts )
!!$    
!!$    ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
!!$    ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
!!$    ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
!!$    ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
!!$    ALLOCATE(spherical_wave3D_deriv(1:numthetapts,1:numphipts))
!!$    
!!$    spherical_wave3D        = ZERO
!!$    spherical_wave3D_dr     = ZERO
!!$    spherical_wave3D_dtheta = ZERO
!!$    spherical_wave3D_dphi   = ZERO
!!$    spherical_wave3D_deriv  = ZERO
!!$    
!!$    rb = Rboundary
!!$    deltar = dr
!!$    
!!$    IF(write_to_file) &
!!$         CALL create_surface_file(filename)
!!$    
!!$  END SUBROUTINE initialize_cartesian3D_surface_serial
!!$  
!!$  !****************************************************************!
!!$#if _COM_MPI
!!$  SUBROUTINE initialize_cartesian3D_surface_parallel(x_ax, y_ax, z_ax, &
!!$       dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, rank, &
!!$       size, comm, write_to_file, filename)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    REAL(dp), INTENT(IN)       :: x_ax(:), y_ax(:), z_ax(:)
!!$    INTEGER, INTENT(IN)        :: dims(:)
!!$    REAL(dp), INTENT(IN)       :: Rboundary
!!$    REAL(dp), INTENT(IN)       :: radius_tolerance
!!$    INTEGER, INTENT(IN)        :: fd_rule
!!$    REAL(dp), INTENT(IN)       :: dr
!!$    INTEGER, INTENT(IN)        :: lmax
!!$    INTEGER, INTENT(IN)        :: rank, size
!!$    INTEGER, INTENT(IN)        :: comm
!!$    LOGICAL, INTENT(IN)        :: write_to_file
!!$    CHARACTER(LEN=*), INTENT(IN) :: filename
!!$    
!!$    INTEGER                    :: numthetaoffset, numphioffset
!!$    
!!$    !------------------------------------------------------------!
!!$    
!!$    ! Initilize finite-difference coefficients
!!$    CALL initialize_fd_coeffs(fd_rule)
!!$    
!!$    ! Initialize spherical grid
!!$    ! For scattered interpolation
!!$    !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
!!$    !     Rboundary, radius_tolerance, fd_rule, dr, lmax, &
!!$    !     rank, size, comm, numpts, numrpts, numthetapts, numphipts, &
!!$    !     surfacerank, numsurfaceprocs, surfacecomm, &
!!$    !     numthetaptsperproc, numphiptsperproc )
!!$    
!!$    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
!!$         Rboundary, fd_rule, dr, lmax, &
!!$         rank, size, comm, numrpts, numthetapts, numphipts, &
!!$         surfacerank, numsurfaceprocs, surfacecomm, &
!!$         numthetaptsperproc, numphiptsperproc )
!!$    
!!$    numthetaoffset = INT(numthetapts / numsurfaceprocs)
!!$    numphioffset = INT(numphipts / numsurfaceprocs) 
!!$    
!!$    IF(i_am_surface(rank) .EQ. 1) THEN
!!$       ALLOCATE(spherical_wave3D_local(1:numrpts,1:numthetapts,1:numphipts))
!!$       ALLOCATE(spherical_wave3D_global(1:numrpts,1:numthetapts,1:numphipts))
!!$       ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
!!$       ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
!!$       ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
!!$       IF((numthetaoffset.EQ.0) .OR. (numphioffset.EQ.0)) THEN
!!$          ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,&
!!$               1:numphipts))
!!$          ALLOCATE(spherical_wave3D_deriv(1:numthetapts,&
!!$               1:numphipts))
!!$       ELSE
!!$          ALLOCATE(spherical_wave3D(1:numrpts,1:numthetaptsperproc,&
!!$               1:numphiptsperproc))
!!$          ALLOCATE(spherical_wave3D_deriv(1:numthetaptsperproc,&
!!$               1:numphiptsperproc))
!!$       ENDIF
!!$       spherical_wave3D_local  = ZERO
!!$       spherical_wave3D_global = ZERO
!!$       spherical_wave3D_dr     = ZERO
!!$       spherical_wave3D_dtheta = ZERO
!!$       spherical_wave3D_dphi   = ZERO
!!$       spherical_wave3D        = ZERO
!!$       spherical_wave3D_deriv  = ZERO
!!$       
!!$       rb = Rboundary
!!$       deltar = dr
!!$       
!!$       IF(write_to_file) &
!!$            CALL create_surface_file(filename, surfacecomm)
!!$       
!!$    ENDIF
!!$    
!!$  END SUBROUTINE initialize_cartesian3D_surface_parallel
!!$#endif  
  !**********************************************************************!
  !**********************************************************************!
    
  !**********************************************************************!

  SUBROUTINE get_correlated_cylindrical_surface(filename, wavefunc, &
       rho_ax, z_ax, bigr_ax, dims, fd_rule, time, efield, afield, &
       lmax, mpi_rank, mpi_size )
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
    REAL(dp), INTENT(IN)         :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)         :: bigr_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax
    INTEGER, INTENT(IN), OPTIONAL :: mpi_rank, mpi_size
    
    INTEGER                      :: numtotalpts
    INTEGER                      :: numthetaoffset
    INTEGER                      :: middle_pt, offset
    COMPLEX(dp)                  :: suma
    INTEGER                      :: ir, itheta, ierror
    INTEGER                      :: ibigr, iorder
    
    !----------------------------------------------------------!
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN
          DO ibigr = 1, numbigrpts
             !! For scattered interpolation
             !CALL get_cylindrical_boundary( wavefunc, spherical_wave2D_local, &
             !     spherical_wave2D_dr, spherical_wave2D_dtheta, &
             !     'quadratic', rank)
             
             CALL get_cylindrical_boundary( rho_ax, z_ax, dims, &
                  wavefunc(ibigr,:,:), fd_rule, spherical_wave2D1D_local(ibigr,:,:), &
                  spherical_wave2D1D_dr, spherical_wave2D1D_dtheta, mpi_rank)
          ENDDO
          
          ! Communicate the spherical wavefunction
          numtotalpts = numrpts * numthetapts * numbigrpts
          spherical_wave2D1D_global = ZERO
#if _COM_MPI
          CALL MPI_ALLREDUCE(spherical_wave2D1D_local, spherical_wave2D1D_global, &
               numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
#endif
          middle_pt = fd_rule + 1
          numthetaoffset = INT(numthetapts / numsurfaceprocs)
          
          IF(numthetaoffset.EQ.0) THEN
             
             spherical_wave2D1D = spherical_wave2D1D_global
             
             ! Project the basis set against the wavefunction at the boundary
             DO iorder = 1, basis_order
                DO itheta = 1, numthetapts
                   DO ir = 1, numrpts
                      suma = ZERO
                      DO ibigr = 1, numbigrpts
                         suma = suma + basis_set(ibigr,iorder) * &
                              spherical_wave2D1D(ibigr,ir,itheta)
                      ENDDO
                      sphericaln_wave2D(iorder,ir,itheta) = suma
                   ENDDO
                ENDDO
                
                CALL make_wave_boundary_derivative(sphericaln_wave2D(iorder,:,:),&
                     sphericaln_wave2D_deriv(iorder,:),fd_rule,deltar,&
                     numrpts,numthetapts)
                
             ENDDO
             ! Write boundary points to a HDF5 file
!!$             IF (surfacerank.EQ.0) &
!!$                  CALL write_correlated_surface_file(filename,&
!!$                  sphericaln_wave2D(:,middle_pt,:), &
!!$                  sphericaln_wave2D_deriv, time, efield, afield, lmax, basis_order)
             
          ELSE
             offset = numthetaoffset * surfacerank       
             
             DO itheta = 1, numthetaptsperproc
                DO ir = 1, numrpts
                   DO ibigr = 1, numbigrpts
                      spherical_wave2D1D(ibigr,ir,itheta) = &
                           spherical_wave2D1D_global(ibigr,ir,offset+itheta)
                   ENDDO
                ENDDO
             ENDDO
             
             ! Project the basis set against the wavefunction at the boundary
             DO iorder = 1, basis_order
                DO itheta = 1, numthetapts
                   DO ir = 1, numrpts
                      suma = ZERO
                      DO ibigr = 1, numbigrpts
                         suma = suma + basis_set(ibigr,iorder) * &
                              spherical_wave2D1D(ibigr,ir,itheta)
                      ENDDO
                      sphericaln_wave2D(iorder,ir,itheta) = suma
                   ENDDO
                ENDDO
                
                CALL make_wave_boundary_derivative(sphericaln_wave2D(iorder,:,:),&
                     sphericaln_wave2D_deriv(iorder,:),fd_rule,deltar,&
                     numrpts,numthetaptsperproc)
                
             ENDDO
             ! Write boundary points to a HDF5 file
!!$             CALL write_correlated_surface_file(filename, &
!!$                  sphericaln_wave2D(:,middle_pt,:), &
!!$                  sphericaln_wave2D_deriv, time, efield, afield, lmax, basis_order,&
!!$                  surfacerank, numsurfaceprocs, surfacecomm, numthetaptsperproc)
          ENDIF
       ENDIF
    ELSE
       
       DO ibigr = 1, numbigrpts
          ! For scattered interpolation
          !CALL get_cylindrical_boundary( wavefunc, spherical_wave2D, &
          !     spherical_wave2D_dr, spherical_wave2D_dtheta, &
          !     'quadratic')
          
          ! For bicubic interpolation
          CALL get_cylindrical_boundary( rho_ax, z_ax, dims, &
               wavefunc(ibigr,:,:), fd_rule, spherical_wave2D1D(ibigr,:,:), &
               spherical_wave2D1D_dr, spherical_wave2D1D_dtheta)
       ENDDO
       
       middle_pt = fd_rule + 1
       ! Project the basis set against the wavefunction at the boundary
       DO iorder = 1, basis_order
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                suma = ZERO
                DO ibigr = 1, numbigrpts
                   suma = suma + basis_set(ibigr,iorder) * &
                        spherical_wave2D1D(ibigr,ir,itheta)
                ENDDO
                sphericaln_wave2D(iorder,ir,itheta) = suma
             ENDDO
          ENDDO
          
          ! Calculate the derivative
          CALL make_wave_boundary_derivative(sphericaln_wave2D(iorder,:,:),&
               sphericaln_wave2D_deriv(iorder,:),fd_rule,deltar,&
               numrpts,numthetapts)
       ENDDO
       
       ! Write boundary points to a HDF5 file
       CALL write_correlated_surface_file(filename, sphericaln_wave2D(:,middle_pt,:), &
            sphericaln_wave2D_deriv, time, efield, afield, lmax, basis_order)
       
    ENDIF
    
    
  END SUBROUTINE get_correlated_cylindrical_surface

  !*****************************************************************!
  
!!$  SUBROUTINE get_cartesian3D_surface_serial(filename, wavefunc, &
!!$       x_ax, y_ax, z_ax, dims, fd_rule, time, efield, afield, &
!!$       lmax, write_to_file)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    CHARACTER(LEN=*), INTENT(IN) :: filename
!!$    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
!!$    REAL(dp), INTENT(IN)         :: x_ax(:), y_ax(:), z_ax(:)
!!$    INTEGER, INTENT(IN)          :: dims(:)
!!$    INTEGER, INTENT(IN)          :: fd_rule
!!$    REAL(dp), INTENT(IN)         :: time
!!$    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
!!$    INTEGER, INTENT(IN)          :: lmax
!!$    LOGICAL, INTENT(IN)          :: write_to_file
!!$    
!!$    INTEGER                      :: middle_pt
!!$    !----------------------------------------------------------!
!!$    
!!$    ! For scattered interpolation
!!$    !CALL get_cartesian_boundary( wavefunc, spherical_wave3D, &
!!$    !     spherical_wave3D_dr, spherical_wave3D_dtheta, &
!!$    !     spherical_wave3D_dphi, 'quadratic')
!!$    
!!$    ! For tricubic interpolation
!!$    CALL get_cartesian_boundary( x_ax, y_ax, z_ax, dims, &
!!$         wavefunc, fd_rule, spherical_wave3D, &
!!$         spherical_wave3D_dr, spherical_wave3D_dtheta, &
!!$         spherical_wave3D_dphi)
!!$    
!!$    middle_pt = fd_rule + 1
!!$    CALL make_wave_boundary_derivative(spherical_wave3D,&
!!$         spherical_wave3D_deriv,fd_rule,deltar,&
!!$         numrpts,numthetapts,numphipts)
!!$    
!!$    ! Write boundary points to a HDF5 file
!!$    IF (write_to_file) &
!!$         CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
!!$         spherical_wave3D_deriv, time, efield, afield, lmax)
!!$    
!!$  END SUBROUTINE get_cartesian3D_surface_serial
!!$  
!!$  !**********************************************************************!
!!$#if _COM_MPI
!!$  SUBROUTINE get_cartesian3D_surface_parallel(filename, wavefunc, &
!!$       x_ax, y_ax, z_ax, dims, fd_rule, time, efield, afield, &
!!$       lmax, rank, write_to_file)
!!$    
!!$    IMPLICIT NONE
!!$
!!$    CHARACTER(LEN=*), INTENT(IN) :: filename
!!$    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
!!$    REAL(dp), INTENT(IN)         :: x_ax(:), y_ax(:), z_ax(:)
!!$    INTEGER, INTENT(IN)          :: dims(:)
!!$    INTEGER, INTENT(IN)          :: fd_rule
!!$    REAL(dp), INTENT(IN)         :: time
!!$    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
!!$    INTEGER, INTENT(IN)          :: lmax, rank
!!$    LOGICAL, INTENT(IN)          :: write_to_file
!!$    
!!$    INTEGER                      :: numtotalpts
!!$    INTEGER                      :: middle_pt
!!$    INTEGER                      :: thetaoffset, phioffset
!!$    INTEGER                      :: numthetaoffset, numphioffset
!!$    INTEGER                      :: ir, itheta, iphi, ierror
!!$    !----------------------------------------------------------!
!!$    
!!$    IF(i_am_surface(rank) .EQ. 1) THEN
!!$       ! For scattered interpolation
!!$       !CALL get_cartesian_boundary( wavefunc, spherical_wave3D_local, &
!!$       !     spherical_wave3D_dr, spherical_wave3D_dtheta, &
!!$       !     spherical_wave3D_dphi, 'quadratic',rank)
!!$
!!$       ! For tricubic interpolation
!!$       CALL get_cartesian_boundary( x_ax, y_ax, z_ax, dims, &
!!$            wavefunc, fd_rule, spherical_wave3D_local, &
!!$            spherical_wave3D_dr, spherical_wave3D_dtheta, &
!!$            spherical_wave3D_dphi, rank)
!!$       
!!$       ! Communicate the spherical wavefunction
!!$       numtotalpts = numrpts * numthetapts * numphipts
!!$       spherical_wave3D_global = ZERO
!!$       CALL MPI_ALLREDUCE(spherical_wave3D_local, spherical_wave3D_global, &
!!$            numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
!!$       
!!$       middle_pt = fd_rule + 1
!!$       numthetaoffset = INT(numthetapts / numsurfaceprocs)
!!$       numphioffset = INT(numphipts / numsurfaceprocs)
!!$       
!!$       IF((numthetaoffset.EQ.0) .OR. (numphioffset.EQ.0)) THEN
!!$          
!!$          spherical_wave3D =  spherical_wave3D_global
!!$          CALL make_wave_boundary_derivative(spherical_wave3D,&
!!$               spherical_wave3D_deriv,fd_rule,deltar,&
!!$               numrpts,numthetapts,numphipts)
!!$          
!!$          ! Write boundary points to a HDF5 file
!!$          IF (write_to_file .AND. (surfacerank.EQ.0)) &
!!$               CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
!!$               spherical_wave3D_deriv, time, efield, afield, lmax)
!!$          
!!$       ELSE
!!$          thetaoffset = numthetaoffset * surfacerank       
!!$          phioffset = numphioffset * surfacerank       
!!$          
!!$          DO iphi = 1, numphiptsperproc
!!$             DO itheta = 1, numthetaptsperproc
!!$                DO ir = 1, numrpts
!!$                   spherical_wave3D(ir,itheta, iphi) = &
!!$                        spherical_wave3D_global(ir,thetaoffset+itheta, phioffset+iphi)
!!$                ENDDO
!!$             ENDDO
!!$          ENDDO
!!$          
!!$          CALL make_wave_boundary_derivative(spherical_wave3D,&
!!$               spherical_wave3D_deriv,fd_rule,deltar,&
!!$               numrpts,numthetaptsperproc,numphiptsperproc)
!!$          
!!$          ! Write boundary points to a HDF5 file
!!$          IF (write_to_file) &
!!$               CALL write_surface_file(filename, spherical_wave3D(middle_pt,:, :), &
!!$               spherical_wave3D_deriv, time, efield, afield, lmax, &
!!$               surfacerank, numsurfaceprocs, surfacecomm, &
!!$               numthetaptsperproc, numphiptsperproc )
!!$       ENDIF
!!$    ENDIF
!!$    
!!$  END SUBROUTINE get_cartesian3D_surface_parallel
!!$#endif
!!$  
  !**********************************************************************!
  !**********************************************************************!

  SUBROUTINE delete_correlated_surface2D(mpi_rank)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN), OPTIONAL    :: mpi_rank
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(PRESENT(mpi_rank)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN  
          DEALLOCATE(spherical_wave2D1D_local,spherical_wave2D1D_global)
       ENDIF
    ENDIF
    DEALLOCATE(sphericaln_wave2D,sphericaln_wave2D_deriv)
    DEALLOCATE(spherical_wave2D1D_dr,spherical_wave2D1D_dtheta)
    DEALLOCATE(spherical_wave2D1D)
    DEALLOCATE(basis_set)
    
  END SUBROUTINE delete_correlated_surface2D

  !----------------------------------------------------------!
  
!!$  SUBROUTINE delete_surface3D_serial()
!!$    
!!$    IMPLICIT NONE
!!$
!!$    ! Deallocate fd coeffs
!!$    CALL delete_fd_coeffs()
!!$    
!!$    DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
!!$    DEALLOCATE(spherical_wave3D_dr)
!!$    DEALLOCATE(spherical_wave3D_dtheta)
!!$    DEALLOCATE(spherical_wave3D_dphi)
!!$    
!!$  END SUBROUTINE delete_surface3D_serial
!!$  
!!$  !-----------------------------------------------------------!
!!$#if _COM_MPI
!!$  SUBROUTINE delete_surface3D_parallel(rank)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    INTEGER, INTENT(IN)         :: rank
!!$    
!!$    ! Deallocate fd coeffs
!!$    CALL delete_fd_coeffs()
!!$    
!!$    IF(i_am_surface(rank) .EQ. 1) THEN
!!$       DEALLOCATE(spherical_wave3D_local,spherical_wave3D_global)
!!$       DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
!!$       DEALLOCATE(spherical_wave3D_dr)
!!$       DEALLOCATE(spherical_wave3D_dtheta)
!!$       DEALLOCATE(spherical_wave3D_dphi)
!!$    ENDIF
!!$    
!!$  END SUBROUTINE delete_surface3D_parallel
!!$#endif
  
END MODULE corrsurface
