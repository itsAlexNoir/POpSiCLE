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

  !-------------------------------------------!
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_local(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dr(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dtheta(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_deriv(:)
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_local(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dtheta(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dphi(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_deriv(:, :)
  
  COMPLEX(dp), ALLOCATABLE  :: psi_l(:, :)
  COMPLEX(dp), ALLOCATABLE  :: psi_lm(:, :, :)
  
  REAL(dp)                  :: rb, deltar
  INTEGER                   :: maxpts
  
CONTAINS
  !********************************************************************!
  !********************************************************************!
  
  SUBROUTINE initialize_cylindrical_surface(rho_ax, z_ax, dims, &
       Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       filename, itime, mpi_rank, mpi_size, comm)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)             :: dims(:)
    REAL(dp), INTENT(IN)            :: Rboundary
    REAL(dp), INTENT(IN)            :: radius_tolerance
    INTEGER, INTENT(IN)             :: fd_rule
    REAL(dp), INTENT(IN)            :: dr
    INTEGER, INTENT(IN)             :: lmax
    CHARACTER(LEN=*), INTENT(IN)    :: filename
    INTEGER, INTENT(IN)             :: itime
    INTEGER, INTENT(IN), OPTIONAL   :: mpi_rank, mpi_size
    INTEGER, INTENT(IN), OPTIONAL   :: comm
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    rb = Rboundary
    deltar = dr
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       ! Initialize spherical grid
       ! For scattered interpolation
       !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
       !     radius_tolerance, fd_rule, dr, lmax, rank, size, comm )
       
       ! For bicubic interpolation
       CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
            fd_rule, dr, lmax, mpi_rank, mpi_size, comm )
       
       IF(i_am_surface(mpi_rank).EQ.1) THEN
          ALLOCATE(spherical_wave2D_local(1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
          ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))

          spherical_wave2D_local  = ZERO
          spherical_wave2D        = ZERO
          
          IF(surfacerank .EQ. 0) THEN
             ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
             
             IF(itime .EQ. 0) &
                  CALL create_surface_file(filename)
          ENDIF
          
       ENDIF
    ELSE
       
       !! Initialize spherical grid
       !! For scattered interpolation uncomment the paragraph below.
       !CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
       !     Rboundary, radius_tolerance, fd_rule, dr, lmax )
       
       !For bicubic interpolation, the paragraph below.
       CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
            fd_rule, dr, lmax)
       
       ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))   
       ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
       
       IF(itime .EQ. 0) &
            CALL create_surface_file(filename)
    ENDIF
    
    spherical_wave2D        = ZERO
    spherical_wave2D_dr     = ZERO
    spherical_wave2D_dtheta = ZERO
    spherical_wave2D_deriv  = ZERO
    
  END SUBROUTINE initialize_cylindrical_surface
  
  !****************************************************************!
  
  !****************************************************************!
  
  SUBROUTINE initialize_cartesian_surface(x_ax, y_ax, z_ax, &
       dims, Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       filename, itime, mpi_rank, mpi_size, comm )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)            :: dims(:)
    REAL(dp), INTENT(IN)           :: Rboundary
    REAL(dp), INTENT(IN)           :: radius_tolerance
    INTEGER, INTENT(IN)            :: fd_rule
    REAL(dp), INTENT(IN)           :: dr
    INTEGER, INTENT(IN)            :: lmax
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(IN)            :: itime
    INTEGER, INTENT(IN), OPTIONAL  :: mpi_rank, mpi_size
    INTEGER, INTENT(IN), OPTIONAL  :: comm
       
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule)
    
    rb = Rboundary
    deltar = dr
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       ! Initialize spherical grid
       ! For scattered interpolation
       !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
       !     Rboundary, radius_tolerance, fd_rule, dr, lmax, &
       !     rank, size, comm )
       
       CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
            Rboundary, fd_rule, dr, lmax, &
            mpi_rank, mpi_size, comm )
              
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN
          ALLOCATE(spherical_wave3D_local(1:numrpts,1:numthetapts,1:numphipts))
          ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
          ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
          ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
          ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
          
          spherical_wave3D_local  = ZERO
          spherical_wave3D        = ZERO
          
          IF(surfacerank .EQ. 0) THEN
             ALLOCATE(spherical_wave3D_deriv(1:numthetapts,&
                  1:numphipts))
             
             IF(itime .EQ. 0) &
                  CALL create_surface_file(filename)
          ENDIF
          
       ENDIF
    ELSE
       
       ! Initialize spherical grid
       ! For scattered interpolation
       !CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
       !     Rboundary, radius_tolerance, fd_rule, dr, lmax )
       
       ! Tricubic interpolation
       CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
            Rboundary, fd_rule, dr, lmax )
       
       ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_deriv(1:numthetapts,1:numphipts))
       
       IF(itime .EQ. 0) &
            CALL create_surface_file(filename)
    ENDIF
    
    spherical_wave3D        = ZERO
    spherical_wave3D_dr     = ZERO
    spherical_wave3D_dtheta = ZERO
    spherical_wave3D_dphi   = ZERO
    spherical_wave3D_deriv  = ZERO
    
  END SUBROUTINE initialize_cartesian_surface
  
  !****************************************************************!
  !****************************************************************!
  
  !**********************************************************************!

  SUBROUTINE get_cylindrical_surface(filename, itime, wavefunc, &
       rho_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, mpi_rank, mpi_size)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN)           :: itime
    COMPLEX(dp), INTENT(IN)       :: wavefunc(:, :)
    REAL(dp), INTENT(IN)          :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)           :: dims(:)
    INTEGER, INTENT(IN)           :: fd_rule
    REAL(dp), INTENT(IN)          :: time
    REAL(dp), INTENT(IN)          :: efield(:), afield(:)
    INTEGER, INTENT(IN)           :: lmax
    INTEGER, INTENT(IN), OPTIONAL :: mpi_rank, mpi_size
    
    INTEGER                      :: numtotalpts
    INTEGER                      :: middle_pt
    INTEGER                      :: ir, itheta, ierror
    
    !----------------------------------------------------------!
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN
          !! For scattered interpolation
          !CALL get_cylindrical_boundary( wavefunc, spherical_wave2D_local, &
          !     spherical_wave2D_dr, spherical_wave2D_dtheta, &
          !     'quadratic', rank)
          
          CALL get_cylindrical_boundary( rho_ax, z_ax, dims, &
               wavefunc, fd_rule, spherical_wave2D_local, &
               spherical_wave2D_dr, spherical_wave2D_dtheta, mpi_rank)
          
          ! Communicate the spherical wavefunction
          numtotalpts = numrpts * numthetapts
          spherical_wave2D        = ZERO
          
#if _COM_MPI
          CALL MPI_ALLREDUCE(spherical_wave2D_local, spherical_wave2D, &
               numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
#endif
          
          IF(surfacerank.EQ.0) THEN
             
             middle_pt = fd_rule + 1    
             
             CALL make_wave_boundary_derivative(spherical_wave2D,&
                  spherical_wave2D_deriv,fd_rule,deltar,&
                  numrpts,numthetapts)
             
             ! Write boundary points to a HDF5 file
             CALL write_surface_file(filename, itime, &
                  spherical_wave2D(middle_pt,:), &
                  spherical_wave2D_deriv, time, efield, afield  )

          ENDIF
       ENDIF
       
    ELSE
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
       CALL write_surface_file(filename, itime, &
            spherical_wave2D(middle_pt,:), &
            spherical_wave2D_deriv, time, efield, afield )
    ENDIF
    
  END SUBROUTINE get_cylindrical_surface
  
  !**********************************************************************!
  !**********************************************************************!
  
  SUBROUTINE get_cartesian_surface(filename, itime, wavefunc, &
       x_ax, y_ax, z_ax, dims, fd_rule, time, efield, afield, &
       lmax, mpi_rank, mpi_size)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: itime
    COMPLEX(dp), INTENT(IN)      :: wavefunc(:, :, :)
    REAL(dp), INTENT(IN)         :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)          :: dims(:)
    INTEGER, INTENT(IN)          :: fd_rule
    REAL(dp), INTENT(IN)         :: time
    REAL(dp), INTENT(IN)         :: efield(:), afield(:)
    INTEGER, INTENT(IN)          :: lmax
    INTEGER, INTENT(IN), OPTIONAL :: mpi_rank, mpi_size
    
    INTEGER                      :: numtotalpts
    INTEGER                      :: middle_pt
    INTEGER                      :: ir, itheta, iphi, ierror
    
    !----------------------------------------------------------!
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN
          ! For scattered interpolation
          !CALL get_cartesian_boundary( wavefunc, spherical_wave3D_local, &
          !     spherical_wave3D_dr, spherical_wave3D_dtheta, &
          !     spherical_wave3D_dphi, 'quadratic',rank)
          
          ! For tricubic interpolation
          CALL get_cartesian_boundary( x_ax, y_ax, z_ax, dims, &
               wavefunc, fd_rule, spherical_wave3D_local, &
               spherical_wave3D_dr, spherical_wave3D_dtheta, &
               spherical_wave3D_dphi, mpi_rank)
          
          ! Communicate the spherical wavefunction
          numtotalpts = numrpts * numthetapts * numphipts
          spherical_wave3D = ZERO

#if _COM_MPI
          CALL MPI_ALLREDUCE(spherical_wave3D_local, spherical_wave3D, &
               numtotalpts, MPI_DOUBLE_COMPLEX, MPI_SUM, surfacecomm, ierror)
#endif
          
          IF(surfacerank .EQ. 0) THEN
             
             middle_pt = fd_rule + 1
             
             CALL make_wave_boundary_derivative(spherical_wave3D,&
                  spherical_wave3D_deriv,fd_rule,deltar,&
                  numrpts,numthetapts,numphipts)
             
             ! Write boundary points to a HDF5 file
             CALL write_surface_file(filename, itime, spherical_wave3D(middle_pt,:, :), &
                  spherical_wave3D_deriv, time, efield, afield )
             
          ENDIF
       ENDIF
       
    ELSE
       
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
       CALL write_surface_file(filename, itime, &
            spherical_wave3D(middle_pt,:, :), &
            spherical_wave3D_deriv, time, efield, afield)
    ENDIF
    
  END SUBROUTINE get_cartesian_surface
  
  !**********************************************************************!
  !**********************************************************************!
  
  SUBROUTINE delete_surface2D(mpi_rank)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN), OPTIONAL      :: mpi_rank
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(PRESENT(mpi_rank)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN  
          DEALLOCATE(spherical_wave2D_local,spherical_wave2D)
          DEALLOCATE(spherical_wave2D_dr,spherical_wave2D_dtheta)
       ENDIF
       IF(surfacerank .EQ. 0) &
            DEALLOCATE(spherical_wave2D_deriv)
       
    ELSE
       DEALLOCATE(spherical_wave2D,spherical_wave2D_deriv)
       DEALLOCATE(spherical_wave2D_dr,spherical_wave2D_dtheta)
    ENDIF
    
  END SUBROUTINE delete_surface2D

  !----------------------------------------------------------!
  
  SUBROUTINE delete_surface3D(mpi_rank)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN), OPTIONAL     :: mpi_rank
    
    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(PRESENT(mpi_rank)) THEN
       IF(i_am_surface(mpi_rank) .EQ. 1) THEN
          DEALLOCATE(spherical_wave3D_local,spherical_wave3D)
          DEALLOCATE(spherical_wave3D_dr)
          DEALLOCATE(spherical_wave3D_dtheta)
          DEALLOCATE(spherical_wave3D_dphi)
       ENDIF
       IF(surfacerank .EQ. 0) &
            DEALLOCATE(spherical_wave3D_deriv)
       
    ELSE
       DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
       DEALLOCATE(spherical_wave3D_dr)
       DEALLOCATE(spherical_wave3D_dtheta)
       DEALLOCATE(spherical_wave3D_dphi)
    ENDIF
    
  END SUBROUTINE delete_surface3D
  
  !----------------------------------------------------------!
  !----------------------------------------------------------!
  !----------------------------------------------------------!
  
  
END MODULE surface
