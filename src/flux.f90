!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: flux
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 05/11/2015
!
! DESCRIPTION:
!> \brief Flux module
!> \details This module contains routines to calculate the flux passing
!> through a surface, and integrating this surface over time.
!
!------------------------------------------------------------------------------
MODULE flux

  USE constants
  USE gaussleg
  USE bessel
  USE sht
  USE io_surface

  IMPLICIT NONE

  PRIVATE

  ! Public routines
  PUBLIC        :: initialize_tsurff
  PUBLIC        :: get_volkov_phase
  !PUBLIC        :: get_flux
  !PUBLIC        :: calculate_flux

  !-------------------------------------------------------------!
  
  ! Public variables
  INTEGER, PUBLIC               :: numkpts
  REAL(dp), PUBLIC              :: dk, kmax
  REAL(dp), PUBLIC, ALLOCATABLE :: kaxis(:)
  INTEGER, PUBLIC               :: ntime
  
  
  ! Private variables
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :), psip_sph(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_lm(:, :), psip_lm(:, :)
  REAL(dp)                      :: rb
  INTEGER                       :: lmax, lmaxtotal, mmax
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts
  REAL(dp), PUBLIC, ALLOCATABLE :: theta_ax(:), costheta_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: phi_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_weights(:)

  !-------------------------------------------------------------!
  !-------------------------------------------------------------!
  
  
CONTAINS
  
  
  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, ddk, kkmax)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: ddk, kkmax
    
    INTEGER                          :: iphi, ik
    
    ! Create momentum axis
    kmax = kkmax
    numkpts = INT( kmax / dk)
    
    ALLOCATE(kaxis(1:numkpts))
    
    DO ik = 1, numkpts
       kaxis(ik) = REAL(ik-1,dp) * dk
    ENDDO
    
    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename,ntime, &
         numthetapts, numphipts, lmaxtotal)
    
    ! Assign  surface radius
    rb = radb
    ! Assign angular momenta
    IF(lmaxtotal.LT.lmax_desired) THEN
       WRITE(*,*) 'Maximum angular momenta desired is greater than&
            & the one on the file. '
       STOP
    ELSE
       lmax = lmax_desired
       IF(numphipts.EQ.1) THEN
          mmax = 1
       ELSE
          mmax = lmax
       ENDIF
    ENDIF
    
    ! Allocate wavefunction arrays
    ALLOCATE(psi_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psi_lm(0:lmax,0:mmax))
    ALLOCATE(psip_lm(0:lmax,0:mmax))
    
    ! Create spherical coordinates
    ALLOCATE(theta_ax(1:numthetapts),costheta_ax(1:numthetapts))
    ALLOCATE(gauss_weights(1:numthetapts))
    ALLOCATE(phi_ax(1:numphipts))
    
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta_ax, gauss_weights) 
    theta_ax = ACOS(costheta_ax)
    
    DO iphi = 1, numphipts
       phi_ax(iphi) = REAL(iphi,dp) * twopi / REAL(numphipts,dp)
    ENDDO
    
    ! Create spherical harmonics
    CALL initialize_spherical_harmonics(lmax,costheta_ax)
    
    
  END SUBROUTINE initialize_tsurff
  
  !**************************************!
  
  SUBROUTINE get_volkov_phase(phase, kaxis, afield, dt)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(INOUT) :: phase(:)
    REAL(dp), INTENT(IN)    :: kaxis(:)
    REAL(dp), INTENT(IN)    :: afield, dt
    
    REAL(dp)                :: newterm(SIZE(kaxis))
    
    newterm(:) = EXP(-ZIMAGONE * 0.5_dp * (kaxis(:) * afield)**2 * dt)
    
    phase = phase * newterm
    
  END SUBROUTINE get_volkov_phase
  
  !****************************************************************!

!!$  SUBROUTINE get_flux(dt, afield, )
!!$    
!!$    IMPLICIT NONE
!!$
!!$    REAL(dp), INTENT(IN)       :: dt
!!$    REAL(dp), INTENT(IN)       :: afield
!!$
!!$    COMPLEX(dp), ALLOCATABLE  :: intflux(:, :)
!!$    !()
!!$    INTEGER                   :: itime
!!$
!!$    !----------------------------------------!
!!$
!!$    ALLOCATE(intflux(0:lmax,0:mmax))
!!$    intflux = ZERO
!!$    
!!$    
!!$    DO itime = 1, ntime
!!$       ! Read flux from file
!!$       CALL read_surface(filename, itime, psi_sph )
!!$       ! Decompose into spherical harmonics
!!$       CALL make_sht(psitheta, theta_ax, gauss_weights, &
!!$            lmax, psil)
!!$       
!!$       
!!$       ! Calculate flux
!!$       intflux = ZERO
!!$       CALL calculate_flux(psil, psilp, lmax, intflux)
!!$       
!!$       CALL get_volkov_phase(volkov_phase, kaxis, afield, dt)
!!$       
!!$       intflux = intflux * volkov_phase
!!$       
!!$       DO il = 0, lmax
!!$          DO ik = 1, numkpts
!!$             bl(ik,il) = bl(ik,il) + intflux(ik,il)
!!$          ENDDO
!!$       ENDDO
!!$       
!!$    ENDDO
!!$    
!!$    bl = bl * dt
!!$    
!!$  END SUBROUTINE get_flux
!!$
!!$  !***************************************!
!!$
!!$  SUBROUTINE calculate_flux
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    
!!$    
!!$  END SUBROUTINE calculate_flux
  
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

END MODULE flux
