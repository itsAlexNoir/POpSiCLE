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
  USE coords
  USE io_pop
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
  !PUBLIC       :: get_cylindrical_flux
  !PUBLIC       :: get_cartesian_flux
  
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
  
!!$  INTERFACE get_cartesian_flux
!!$     MODULE PROCEDURE get_cartesian2D_flux
!!$     MODULE PROCEDURE get_cartesian3D_flux
!!$  END INTERFACE get_cartesian_flux
  
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

  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dr(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_dtheta(:, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave2D_deriv(:)
  
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dr(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dtheta(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_dphi(:, :, :)
  COMPLEX(dp), ALLOCATABLE  :: spherical_wave3D_deriv(:, :)
  
  COMPLEX(dp), ALLOCATABLE  :: psi_l(:, :)
  COMPLEX(dp), ALLOCATABLE  :: psi_lm(:, :, :)

  REAL(dp)                  :: rb
  INTEGER                   :: numpts, numrpts
  INTEGER                   :: numthetapts, numphipts
  INTEGER                   :: surfacecomm
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

    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule,dr)
    
    ! Initialize spherical grid
    CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
         Rboundary, radius_tolerance, fd_rule, dr, lmax, numpts, &
         numrpts, numthetapts )
    
    ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))
    ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
    ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))
    ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
   
    rb = Rboundary
    
    IF(write_to_file) &
         CALL open_surface_file(filename)
    
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
        
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule,dr)
    
    ! Initialize spherical grid
    CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, Rboundary, &
         radius_tolerance, fd_rule, dr, lmax, rank, size, comm, &
         numpts, numrpts, numthetapts, numsurfaceprocs, surfacecomm )
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       ALLOCATE(spherical_wave2D(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_dr(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_dtheta(1:numrpts,1:numthetapts))
       ALLOCATE(spherical_wave2D_deriv(1:numthetapts))
       
       rb = Rboundary
       
       IF(write_to_file) &
            CALL open_surface_file(filename, surfacecomm)
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
    CALL initialize_fd_coeffs(fd_rule,dr)
    
    ! Initialize spherical grid
    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
         Rboundary, radius_tolerance, fd_rule, dr, lmax, numpts, &
         numrpts, numthetapts, numphipts )
    
    ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(spherical_wave3D_deriv(1:numthetapts,1:numphipts))
    
    rb = Rboundary
    
    IF(write_to_file) &
         CALL open_surface_file(filename)
    
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
    
    !------------------------------------------------------------!
    
    ! Initilize finite-difference coefficients
    CALL initialize_fd_coeffs(fd_rule,dr)
   
    ! Initialize spherical grid
    CALL initialize_cartesian_boundary(x_ax, y_ax, z_ax, dims, &
         Rboundary, radius_tolerance, fd_rule, dr, lmax, &
         rank, size, comm, numpts, numrpts, numthetapts, numphipts, &
         numsurfaceprocs, surfacecomm )
   
    IF(i_am_surface(rank) .EQ. 1) THEN
       ALLOCATE(spherical_wave3D(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dr(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dtheta(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_dphi(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(spherical_wave3D_deriv(1:numthetapts,1:numphipts))
       
       rb = Rboundary
       
       IF(write_to_file) &
            CALL open_surface_file(filename, surfacecomm)
    ENDIF
    
  END SUBROUTINE initialize_cartesian3D_surface_parallel
#endif  
  !****************************************************************!
  !****************************************************************!
  
  SUBROUTINE get_cylindrical_surface_serial(wavefunc, fd_rule, time, &
       efield, afield, lmax, write_to_file)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: wavefunc(:, :)
    INTEGER, INTENT(IN)       :: fd_rule
    REAL(dp), INTENT(IN)      :: time, efield, afield
    INTEGER, INTENT(IN)       :: lmax
    LOGICAL, INTENT(IN)       :: write_to_file
    
    INTEGER                   :: middle_pt
    !----------------------------------------------------------!
    
    CALL get_cylindrical_boundary( wavefunc, spherical_wave2D, &
         spherical_wave2D_dr, spherical_wave2D_dtheta, &
         'quadratic')
    
    middle_pt = fd_rule + 1
    CALL make_wave_boundary_derivative(spherical_wave2D,&
         spherical_wave2D_deriv,fd_rule,numrpts,numthetapts)
    
    ! Write boundary points to a HDF5 file
    IF (write_to_file) &
         CALL write_surface_file(spherical_wave2D(middle_pt,:), &
         spherical_wave2D_deriv, time, efield, afield, lmax)
        
  END SUBROUTINE get_cylindrical_surface_serial

  !**********************************************************************!
#if _COM_MPI
  SUBROUTINE get_cylindrical_surface_parallel(wavefunc, fd_rule, time, &
       efield, afield, lmax, rank, write_to_file)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: wavefunc(:, :)
    INTEGER, INTENT(IN)       :: fd_rule
    REAL(dp), INTENT(IN)      :: time, efield, afield
    INTEGER, INTENT(IN)       :: lmax, rank
    LOGICAL, INTENT(IN)       :: write_to_file
    
    INTEGER                   :: middle_pt
    !----------------------------------------------------------!
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       CALL get_cylindrical_boundary( wavefunc, spherical_wave2D, &
            spherical_wave2D_dr, spherical_wave2D_dtheta, &
            'quadratic')
       
       middle_pt = fd_rule + 1
       CALL make_wave_boundary_derivative(spherical_wave2D,&
            spherical_wave2D_deriv,fd_rule,numrpts,numthetapts)
       
       ! Write boundary points to a HDF5 file
       IF (write_to_file) &
            CALL write_surface_file(spherical_wave2D(middle_pt,:), &
            spherical_wave2D_deriv, time, efield, afield, lmax, numthetapts)
    ENDIF
    
  END SUBROUTINE get_cylindrical_surface_parallel
#endif
  !*****************************************************************!
  
  SUBROUTINE get_cartesian3D_surface_serial(wavefunc, fd_rule, time, &
       efield, afield, lmax, write_to_file)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: wavefunc(:, :, :)
    INTEGER, INTENT(IN)       :: fd_rule
    REAL(dp), INTENT(IN)      :: time, efield, afield
    INTEGER, INTENT(IN)       :: lmax
    LOGICAL, INTENT(IN)       :: write_to_file
    
    INTEGER                   :: middle_pt
    !----------------------------------------------------------!
    
    CALL get_cartesian_boundary( wavefunc, spherical_wave3D, &
         spherical_wave3D_dr, spherical_wave3D_dtheta, &
         spherical_wave3D_dphi, 'quadratic')
    
    middle_pt = fd_rule + 1
    CALL make_wave_boundary_derivative(spherical_wave3D,&
         spherical_wave3D_deriv,fd_rule,numrpts,numthetapts,numphipts)
    
    ! Write boundary points to a HDF5 file
    IF (write_to_file) &
         CALL write_surface_file(spherical_wave3D(middle_pt,:, :), &
         spherical_wave3D_deriv, time, efield, afield, lmax)
    
  END SUBROUTINE get_cartesian3D_surface_serial
  
  !**********************************************************************!
#if _COM_MPI
  SUBROUTINE get_cartesian3D_surface_parallel(wavefunc, fd_rule, time, &
       efield, afield, lmax, rank, write_to_file)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: wavefunc(:, :, :)
    INTEGER, INTENT(IN)       :: fd_rule
    REAL(dp), INTENT(IN)      :: time, efield, afield
    INTEGER, INTENT(IN)       :: lmax, rank
    LOGICAL, INTENT(IN)       :: write_to_file
    
    INTEGER                   :: middle_pt
    !----------------------------------------------------------!
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       CALL get_cartesian_boundary( wavefunc, spherical_wave3D, &
            spherical_wave3D_dr, spherical_wave3D_dtheta, &
            spherical_wave3D_dphi, 'quadratic')
       
       middle_pt = fd_rule + 1
       CALL make_wave_boundary_derivative(spherical_wave3D,&
            spherical_wave3D_deriv,fd_rule,numrpts,numthetapts, numphipts)
       
       ! Write boundary points to a HDF5 file
       IF (write_to_file) &
            CALL write_surface_file(spherical_wave3D(middle_pt,:, :), &
            spherical_wave3D_deriv, time, efield, afield, lmax, &
            numthetapts, numphipts )
    ENDIF
    
  END SUBROUTINE get_cartesian3D_surface_parallel
#endif
  
  !**********************************************************************!
  !**********************************************************************!
  
  SUBROUTINE delete_surface2D_serial()
    
    IMPLICIT NONE
    
    ! Close HDF5 file
    CALL close_surface_file()
    
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
    
    ! Close HDF5 parallel file
    CALL close_surface_file(rank)

    ! Deallocate fd coeffs
    CALL delete_fd_coeffs()
    
    IF(i_am_surface(rank) .EQ. 1) THEN  
       DEALLOCATE(spherical_wave2D,spherical_wave2D_deriv)
       DEALLOCATE(spherical_wave2D_dr,spherical_wave2D_dtheta)
    ENDIF
    
  END SUBROUTINE delete_surface2D_parallel
#endif
  !----------------------------------------------------------!
  
  SUBROUTINE delete_surface3D_serial()
    
    IMPLICIT NONE
    
    ! Close HDF5 file
    CALL close_surface_file()
    
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
    
    ! Close HDF5 parallel file
    CALL close_surface_file(rank)
    
    IF(i_am_surface(rank) .EQ. 1) THEN
       DEALLOCATE(spherical_wave3D,spherical_wave3D_deriv)
       DEALLOCATE(spherical_wave3D_dr)
       DEALLOCATE(spherical_wave3D_dtheta)
       DEALLOCATE(spherical_wave3D_dphi)
    ENDIF
    
  END SUBROUTINE delete_surface3D_parallel
#endif
  !***********************************************************!
  
!!$  SUBROUTINE get_cylindrical_flux( wavefunc )
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    COMPLEX(dp), INTENT(IN)   :: wavefunc(:,:)
!!$    
!!$    !-----------------------------------------!
!!$    
!!$    !ALLOCATE(psi_l(1:numrpts,0:lmax))
!!$    ! Initialize Spherical Harmonics
!!$    !CALL initialize_spherical_harmonics(lmax, costheta_boundary)
!!$
!!$    
!!$    CALL make_sht(sphfunc, costheta_boundary, &
!!$         theta_weights, lmax, psi_l)
!!$
!!$    CALL calculate_flux( psi_l)
!!$    
!!$  END SUBROUTINE get_cylindrical_flux
  
  !---------------------------------------------------------------!
  
  
  !---------------------------------------------------------------!
  
!!$  SUBROUTINE calculate_radialflux2D(fl, vectorpot )
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    COMPLEX(dp), INTENT(IN)      :: fl(:,:)
!!$    REAL(dp), INTENT(IN)         :: vectorpot 
!!$    
!!$    COMPLEX(dp)                  :: term1, term2, term3
!!$    INTEGER                      :: ik
!!$    
!!$    DO ik = 1, numkpts
!!$       
!!$       term1 = ZERO
!!$       term2 = ZERO
!!$       term3 = ZERO
!!$       
!!$       DO il = 0, lmax
!!$          
!!$          CALL sphbessjy(il,x,sphbessel,sy,sphbesselp,syp)
!!$          
!!$          term1 = term1 + krpts(ik) * rb * 0.5_dp * sphbesselp - &
!!$               0.5_dp * sphbessel ) * fl(indxrb,il)
!!$          
!!$          
!!$          term2 = term2 - rb * 0.5_dp * sphbessel * flp 
!!$          
!!$          term3 = term3 - ZIMAG * 0.5_dp / SQRT(pi) * rb * &
!!$               vectorpot * fl(indxrb,il) * 
!!$          
!!$          DO ill = 0, lmax
!!$             
!!$          ENDDO
!!$
!!$       ENDDO
!!$
!!$    ENDDO
!!$          
!!$  END SUBROUTINE calculate_flux

  
END MODULE surface
