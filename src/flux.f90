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
  PUBLIC        :: calculate_time_integrand
  PUBLIC        :: get_angular_resolution
  
  !-------------------------------------------------------------!
  
  ! Public variables
  REAL(dp), PUBLIC, ALLOCATABLE :: k_ax(:)
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
  INTEGER                       :: lmax_total
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts
  
  !-------------------------------------------------------------!
  !-------------------------------------------------------------!
  
CONTAINS
  
  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, &
       ddk, kkmax, maxkpts, lmax_total, mmax)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: ddk, kkmax
    INTEGER, INTENT(OUT)             :: maxkpts
    INTEGER, INTENT(OUT)             :: lmax_total
    INTEGER, INTENT(OUT)             :: mmax

    INTEGER                          :: lmax, mmin
    INTEGER                          :: iphi, ik

    !-------------------------------------------------!
    
    ! Create momentum axis
    kmax = kkmax
    dk = ddk
    numkpts = INT( kmax / dk)
    maxkpts = numkpts
    
    ALLOCATE(k_ax(1:numkpts))
    
    DO ik = 1, numkpts
       k_ax(ik) = REAL(ik,dp) * dk
    ENDDO
    
    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename,ntime, &
         numthetapts, numphipts, lmax_total)

    ! Assign  surface radius
    rb = radb
    ! Assign angular momenta
    IF(lmax_total.LT.lmax_desired) THEN
       WRITE(*,*) 'Maximum angular momenta desired is greater than&
            & the one on the file. '
       STOP
    ELSE
       lmax = lmax_desired
       IF(numphipts.EQ.1) THEN
          mmax = 0
          mmin = 0
       ELSE
          mmax = lmax
          mmin = -lmax
       ENDIF
    ENDIF
    
    ! Allocate wavefunction arrays
    ALLOCATE(psi_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psi_lm(mmin:mmax,0:lmax))
    ALLOCATE(psip_lm(mmin:mmax,0:lmax))
    
    
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
  
  SUBROUTINE get_flux(filename, lmax, mmax, b_lm)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: lmax, mmax
    COMPLEX(dp), INTENT(OUT)     :: b_lm(:,-mmax:, 0:)
    
    REAL(dp)                     :: time1, time2, dt
    REAL(dp), ALLOCATABLE        :: time(:)
    REAL(dp), ALLOCATABLE        :: efield(:), afield(:)
    COMPLEX(dp), ALLOCATABLE     :: intflux(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: volkov_phase(:)
    INTEGER                      :: mmin
    INTEGER                      :: itime, ik, il, im
    
    !----------------------------------------!
    
    mmin = -mmax
    ! Allocate and set to zero
    ALLOCATE(intflux(1:numkpts,:mmax,0:lmax))
    ALLOCATE(volkov_phase(1:numkpts))
    ALLOCATE(time(ntime))
    ALLOCATE(efield(ntime),afield(ntime))
    
    intflux = ZERO
    b_lm = ZERO
    volkov_phase = ZERO
    
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
       CALL make_sht(psi_sph, lmax, mmax, gauss_weights, &
            psi_lm)
       ! His first derivative
       CALL make_sht(psip_sph, lmax, mmax, gauss_weights, &
            psip_lm)
       
       
       ! Calculate flux
       intflux = ZERO
       CALL calculate_time_integrand(psi_lm, psip_lm, &
            lmax, mmax, afield(itime), intflux)
       
       ! Get the Volkov phase (one per time step)
       CALL get_volkov_phase(volkov_phase, k_ax, afield(itime))
       
       DO il = 0, lmax
          DO im = mmin, mmax
             DO ik = 1, numkpts
                intflux(ik,im,il) = intflux(ik,im,il) * &
                     volkov_phase(ik)
             ENDDO
          ENDDO
       ENDDO
       
       IF(itime.EQ.1) THEN
          time1 = time(itime)
          b_lm = b_lm + intflux * 0.5_dp
       ELSEIF(itime.EQ.2) THEN
          time2 = time(itime)
          b_lm = b_lm + intflux
       ELSEIF(itime.EQ.ntime) THEN
          b_lm = b_lm + intflux * 0.5_dp
       ELSE
          b_lm = b_lm + intflux
       ENDIF
       
    ENDDO
    
    dt = time2 - time1
    
    b_lm = b_lm * ZIMAGONE * dt * EXP(dt)
    
  END SUBROUTINE get_flux
  
  !***************************************!
  
  SUBROUTINE calculate_time_integrand( func_lm, funcp_lm, &
       lmax, mmax, afield, integrand )
    
    IMPLICIT NONE
    
    
    COMPLEX(dp), INTENT(IN)       :: func_lm(-mmax:, 0:)
    COMPLEX(dp), INTENT(IN)       :: funcp_lm(-mmax:, 0:)
    INTEGER, INTENT(IN)           :: lmax, mmax
    REAL(dp), INTENT(IN)          :: afield
    COMPLEX(dp), INTENT(OUT)      :: integrand(:, :, :)
    
    COMPLEX(dp)                   :: term1, term2
    COMPLEX(dp)                   :: term3, term4  
    INTEGER                       :: maxfactlog, mmin
    REAL(dp), ALLOCATABLE         :: factlog(:)
    REAL(dp), ALLOCATABLE         :: krb_ax(:)
    REAL(dp), ALLOCATABLE         :: jl(:, :), jlp(:, :)
    REAL(dp), ALLOCATABLE         :: jy(:), jyp(:)
    REAL(dp)                      :: sqrtcoeff
    INTEGER                       :: ik, il, ill, im
    
    !---------------------------------------------------!

    maxfactlog = 3 * lmax + 1
    mmin = -mmax
    
    ALLOCATE(krb_ax(1:numkpts))
    ALLOCATE(jl(0:lmax,1:numkpts),jlp(0:lmax,1:numkpts))
    ALLOCATE(jy(1:numkpts),jyp(1:numkpts))
    ALLOCATE(factlog(1:maxfactlog))
    
    integrand = ZERO
    jl = 0.0_dp
    jlp = 0.0_dp
    krb_ax = k_ax * rb
    
    DO il = 0, lmax       
       CALL sphbessjy(il,krb_ax,jl(il,:),jy,jlp(il,:),jyp)
    ENDDO
    
    CALL make_factlog(factlog, maxfactlog)
    
    DO ik = 1, numkpts  
       DO il = 0, lmax
          DO im = mmin, mmax
             
             term1 = ZERO
             term2 = ZERO
             term3 = ZERO
             term4 = ZERO
             
             term1 = 0.5_dp * (-ZIMAGONE)**il * &
                  (krb_ax(ik) * jlp(il,ik) - jl(il,ik)) * &
                  psi_lm(im,il)
             
             term2 =  - 0.5_dp * (-ZIMAGONE)**il * &
                  rb * jl(il,ik) * psip_lm(im,il)
             
             term3 = 0.5_dp * ZIMAGONE * afield * rb * &
                  psi_lm(im,il)
             
             DO ill = 0, lmax
                write(*,*) 'loop', ill
                write(*,*) ABS(im).GT.il .OR. 0.GT.1 .OR. ABS(im).GT.ill .OR. &
                     (2 * MAX(il, 1, ill)).GT.(il +1 + ill) .OR.              &
                     (0 + 0).NE.0
                
                write(*,*) ABS(im).GT.il
                write(*,*) 0.GT.1
                write(*,*) ABS(im).GT.ill
                write(*,*) (2 * MAX(il, 1, ill)).GT.(il +1 + ill)
                write(*,*) (0 + im).NE.im



                
                sqrtcoeff = SQRT(REAL(2 * il + 1,dp)) / REAL(2 * ill + 1,dp)
                
                term4 = term4 + (-ZIMAGONE)**ill * &
                     sqrtcoeff * jl(ill,ik) * &
                     clebsch(il, 1, ill, im, 0, im, factlog, maxfactlog) !* &
                     !clebsch(il, 1, ill, 0, 0, 0, factlog, maxfactlog)
             ENDDO
             term3 = term3 * term4
             integrand(ik,im,il) = term1 + term2 + term3
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(krb_ax)
    DEALLOCATE(jl,jlp)
    DEALLOCATE(jy,jyp)
    DEALLOCATE(factlog)
    
  END SUBROUTINE calculate_time_integrand
  
  !***********************************************************!
  
  SUBROUTINE get_angular_resolution(b_lm, dPomega, lmax, mmax)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: b_lm(:, -mmax:, 0:)
    COMPLEX(dp), INTENT(OUT)      :: dPomega(:, :, :)
    INTEGER, INTENT(IN)           :: lmax,mmax
    
    INTEGER                       :: maxfactlog, mmin
    REAL(dp), ALLOCATABLE         :: factlog(:)
    REAL(dp), ALLOCATABLE         :: rnormal(:, :)
    REAL(dp), ALLOCATABLE         :: legenpl(:, :, :)
    INTEGER                       :: ik,il,im,itheta,iphi
    !---------------------------------------------------------!
    
    dPomega = ZERO
    maxfactlog = 3 * lmax + 1
    
    ALLOCATE(factlog(1:maxfactlog))
    ALLOCATE(rnormal(mmin:mmax,0:lmax))
    ALLOCATE(legenpl(0:numthetapts-1,mmin:mmax,0:lmax))
    
    
    CALL make_factlog(factlog,maxfactlog)
    CALL make_rnormal(rnormal, factlog, lmax, maxfactlog)
    CALL make_legendre(legenpl, costheta_ax, lmax, numthetapts)
    
    DO ik = 1, numkpts
       DO il = 0, lmax
          DO im = mmin, mmax
             IF(im.LT.0) THEN
                DO iphi = 1, numphipts
                   DO itheta = 1, numthetapts
                      dPomega(ik,itheta,iphi) = b_lm(ik,im,il) * &
                           1.0_dp / k_ax(ik) * (-1.0_dp)**ABS(im) * &
                           rnormal(ABS(im),il) * legenpl(itheta-1,ABS(im),il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ENDDO
                ENDDO
             ELSE
                DO iphi = 1, numphipts
                   DO itheta = 1, numthetapts
                      dPomega(ik,itheta,iphi) = b_lm(ik,im,il) * &
                           1.0_dp / k_ax(ik) * &
                           rnormal(im,il) * legenpl(itheta-1,im,il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(factlog,rnormal,legenpl)
    
    
  END SUBROUTINE get_angular_resolution
    
  !***********************************************************!
  !***********************************************************!
  
  
END MODULE flux
