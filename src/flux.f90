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
  PUBLIC        :: get_flux
  !PUBLIC        :: calculate_flux
  
  !-------------------------------------------------------------!
  
  ! Public variables
  REAL(dp), PUBLIC, ALLOCATABLE :: kaxis(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: theta_ax(:), costheta_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: phi_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_weights(:)
  
  
  ! Private variables
  INTEGER                       :: numkpts
  REAL(dp)                      :: dk, kmax
  INTEGER                       :: ntime
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :), psip_sph(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_lm(:, :), psip_lm(:, :)
  REAL(dp)                      :: rb
  INTEGER                       :: lmax, lmaxtotal, mmax
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts
  
  !-------------------------------------------------------------!
  !-------------------------------------------------------------!
  
CONTAINS
  
  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, &
       ddk, kkmax, maxkpts, mmax_out)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: ddk, kkmax
    INTEGER, INTENT(OUT)             :: maxkpts
    INTEGER, INTENT(OUT)             :: mmax_out
    INTEGER                          :: iphi, ik

    !-------------------------------------------------!
    
    ! Create momentum axis
    kmax = kkmax
    numkpts = INT( kmax / dk)
    maxkpts = numkpts
    
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
          mmax = 0
       ELSE
          mmax = lmax
       ENDIF
    ENDIF
    
    mmax_out = mmax
    
    ! Allocate wavefunction arrays
    ALLOCATE(psi_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psi_lm(0:mmax,0:lmax))
    ALLOCATE(psip_lm(0:mmax,0:lmax))
    
    
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
  
  SUBROUTINE get_volkov_phase(phase, kaxis, afield)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(INOUT) :: phase(:)
    REAL(dp), INTENT(IN)       :: kaxis(:)
    REAL(dp), INTENT(IN)       :: afield
    
    COMPLEX(dp)                :: newterm(SIZE(kaxis))
    INTEGER                    :: ik
    !-------------------------------------------------!
    
    DO ik = 1, SIZE(kaxis)
       newterm(:) = EXP(-ZIMAGONE * 0.5_dp * &
            (kaxis(ik) * afield)**2)
    ENDDO
    
    phase = phase * newterm
    
  END SUBROUTINE get_volkov_phase
  
  !****************************************************************!
  
  SUBROUTINE get_flux(filename, blm)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    COMPLEX(dp), INTENT(OUT)     :: blm(:, :, :)

    REAL(dp)                     :: dt
    REAL(dp), ALLOCATABLE        :: time(:)
    REAL(dp), ALLOCATABLE        :: efield(:), afield(:)
    COMPLEX(dp), ALLOCATABLE     :: intflux(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: volkov_phase(:)
    INTEGER                      :: itime, ik, il, im
    
    !----------------------------------------!
    
    ALLOCATE(intflux(1:numkpts,0:lmax,0:mmax))
    intflux = ZERO
    ALLOCATE(volkov_phase(1:numkpts))
    volkov_phase = ZERO
    ALLOCATE(time(ntime))
    ALLOCATE(efield(ntime),afield(ntime))

    WRITE(*,*)
    WRITE(*,*)
    
    ! We begin the time loop!!
    DO itime = 1, ntime

       WRITE(*,'(A,I3)') 'Begin loop number: ',itime
       
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            psi_sph, psip_sph, time(itime), efield(itime), afield(itime) )
       
       ! Decompose into spherical harmonics:
       ! The wavefunction
       CALL make_sht(psi_sph, lmax, gauss_weights, &
            psi_lm)
       ! His first derivative
       !CALL make_sht(psip_sph, lmax, gauss_weights, &
       !     psip_lm)
       
       
       ! Calculate flux
       intflux = ZERO
       !CALL calculate_flux(psi_lm, psip_lm, lmax, intflux)
       
       ! Get the Volkov phase (one per time step)
       CALL get_volkov_phase(volkov_phase, kaxis, afield(itime))
       
       DO im = 0, mmax
          DO il = 0, lmax
             DO ik = 1, numkpts
                blm(ik,il,im) = blm(ik,il,im) + &
                     intflux(ik,il,im) * volkov_phase(ik)
             ENDDO
          ENDDO
       ENDDO
       
    ENDDO
    
    !blm = blm * dt
    
  END SUBROUTINE get_flux

  !***************************************!
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
