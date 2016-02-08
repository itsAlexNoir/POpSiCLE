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
  USE omp_lib

  IMPLICIT NONE

  PRIVATE

  ! Public routines
  PUBLIC        :: initialize_tsurff
  PUBLIC        :: get_volkov_phase
  PUBLIC        :: get_flux
  PUBLIC        :: calculate_time_integrand

  !-------------------------------------------------------------!

  ! Public variables
  REAL(dp), PUBLIC, ALLOCATABLE :: k_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: theta_ax(:), costheta_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: phi_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_weights(:)


  ! Private variables
  INTEGER                       :: numkpts
  REAL(dp)                      :: dt, dk, kmax, kmax_th
  INTEGER                       :: ntime
  INTEGER                       :: maxfactlog
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :), psip_sph(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_lm(:, :), psip_lm(:, :)
  REAL(dp)                      :: rb
  REAL(dp), ALLOCATABLE         :: jl(:, :), jlp(:, :)
  REAL(dp), ALLOCATABLE         :: jy(:), jyp(:)
  REAL(dp), ALLOCATABLE         :: krb_ax(:)
  REAL(dp), ALLOCATABLE         :: factlog(:)
  INTEGER                       :: lmax_total
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts

  !-------------------------------------------------------------!
  !-------------------------------------------------------------!

CONTAINS

  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, &
       kmax_input, maxkpts, maxthetapts, maxphipts, &
       lmax_total, mmax)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: kmax_input
    INTEGER, INTENT(OUT)             :: maxkpts
    INTEGER, INTENT(OUT)             :: maxthetapts
    INTEGER, INTENT(OUT)             :: maxphipts
    INTEGER, INTENT(OUT)             :: lmax_total
    INTEGER, INTENT(OUT)             :: mmax

    INTEGER                          :: lmax, mmin
    INTEGER                          :: iphi, ik, il

    !-------------------------------------------------!

    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename, ntime, dt, &
         numthetapts, numphipts, lmax_total)

    maxthetapts = numthetapts
    maxphipts = numphipts

    ! Create momentum axis
    kmax_th = twopi / dt
    dk = twopi / ntime / dt

    kmax = kmax_input
    if(kmax.GT.kmax_th) THEN
       WRITE(*,*) 'Kmax desired large than the theoretical.'
       STOP
    ENDIF

    WRITE(*,'(A,F10.6)')  'Theoretical kmax:                ',kmax_th
    WRITE(*,'(A,F9.6)')  'Desired kmax:                    ',kmax
    WRITE(*,'(A,F9.6)')  'dk:                              ',dk
    WRITE(*,*)           '---------------------------------------------------------&
         &---------------------'
    WRITE(*,*)

    numkpts = INT( kmax / dk)
    maxkpts = numkpts

    ALLOCATE(k_ax(1:numkpts))

    DO ik = 1, numkpts
       k_ax(ik) = REAL(ik,dp) * dk
    ENDDO


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
    maxfactlog = 3 * lmax + 1
    ALLOCATE(factlog(0:maxfactlog))
    CALL make_factlog(factlog, maxfactlog)

    CALL initialize_spherical_harmonics(lmax,costheta_ax)

    ! Initialize Bessel functions
    ALLOCATE(jl(0:lmax,1:numkpts))
    ALLOCATE(jlp(0:lmax,1:numkpts))
    ALLOCATE(jy(1:numkpts),jyp(1:numkpts))
    ALLOCATE(krb_ax(1:numkpts))

    krb_ax = k_ax * rb

    jl  = 0.0_dp
    jlp = 0.0_dp
    jy  = 0.0_dp
    jyp = 0.0_dp

    DO il = 0, lmax
       CALL sphbessjy(il,krb_ax,jl(il,:),jy,jlp(il,:),jyp)
    ENDDO


  END SUBROUTINE initialize_tsurff

  !**************************************!

  SUBROUTINE get_volkov_phase(phase, afield, deltat)

    IMPLICIT NONE

    COMPLEX(dp), INTENT(INOUT) :: phase(:, :, :)
    REAL(dp), INTENT(IN)       :: afield(:)
    REAL(dp), INTENT(IN)       :: deltat

    REAL(dp), ALLOCATABLE      :: newterm(:, :, :)
    REAL(dp)                   :: term1, term2
    REAL(dp)                   :: term3, term4
    REAL(dp)                   :: afieldsq
    INTEGER                    :: ik, itheta, iphi

    !-------------------------------------------------!

    ALLOCATE(newterm(1:numkpts, 1:numthetapts, 1:numphipts))

    afieldsq = afield(1)**2 + afield(2)**2 + &
         afield(3)**2

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
    !$OMP& PRIVATE(term1,term2,term3,term4)

    !$OMP DO COLLAPSE(3)
    DO ik = 1, numkpts
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             term1 =  k_ax(ik) * k_ax(ik)

             term2 = k_ax(ik) * afield(1) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             term3 = k_ax(ik) * afield(2) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             term4 = k_ax(ik) * afield(3) * COS(theta_ax(itheta))
             newterm(ik,itheta,iphi) = term1 + afieldsq - 2.0_dp * (term2 + term3 + term4)

             phase(ik,itheta,iphi) = phase(ik,itheta,iphi) * &
                  EXP(-ZIMAGONE * 0.5_dp * newterm(ik,itheta,iphi) * deltat)
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL

    DEALLOCATE(newterm)

  END SUBROUTINE get_volkov_phase

  !****************************************************************!

  SUBROUTINE get_flux(filename, lmax, mmax, bk)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: lmax, mmax
    COMPLEX(dp), INTENT(OUT)     :: bk(:, :, :)

    REAL(dp), ALLOCATABLE        :: time(:)
    REAL(dp), ALLOCATABLE        :: efield(:, :)
    REAL(dp), ALLOCATABLE        :: afield(:, :)
    COMPLEX(dp), ALLOCATABLE     :: intflux(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: volkov_phase(:, :, :)
    INTEGER                      :: mmin
    INTEGER                      :: itime, il, im
    INTEGER                      :: ik, itheta, iphi

    !----------------------------------------------------!

    mmin = -mmax
    ! Allocate and set to zero
    ALLOCATE(intflux(1:numkpts,1:numthetapts,1:numphipts))
    ALLOCATE(volkov_phase(1:numkpts,1:numthetapts,1:numphipts))
    ALLOCATE(time(ntime))
    ALLOCATE(efield(3,ntime),afield(3,ntime))

    intflux = ZERO
    bk = ZERO
    volkov_phase = ZONE

    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '----------------------'
    WRITE(*,*) 'Begin time integral...'

    ! We begin the time loop!!
    DO itime = 1, ntime
       
       !WRITE(*,'(A,I6)') 'Loop number: ',itime
       
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            psi_sph, psip_sph, time(itime), efield(:,itime), afield(:,itime) )
       
       ! Decompose into spherical harmonics:
       ! The wavefunction
       CALL make_sht(psi_sph, lmax, mmax, gauss_weights, &
            psi_lm)
       ! His first derivative
       CALL make_sht(psip_sph, lmax, mmax, gauss_weights, &
            psip_lm)
       ! Divide by rb
       psi_lm = psi_lm * rb
       psip_lm = psip_lm * rb
       
       
       ! Calculate flux
       intflux = ZERO
       CALL calculate_time_integrand(psi_lm, psip_lm, &
            lmax, mmax, afield(:,itime), intflux)
       
       ! Get the Volkov phase (one per time step)
       CALL get_volkov_phase(volkov_phase, afield(:,itime), dt)
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi)

       !$OMP DO COLLAPSE(3)
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ik = 1, numkpts

                intflux(ik,itheta,iphi) = intflux(ik,itheta,iphi) * &
                     volkov_phase(ik,itheta,iphi)

                IF((itime.EQ.1) .OR. (itime.EQ.ntime)) THEN
                   bk(ik,itheta,iphi) = bk(ik,itheta,iphi) + &
                        intflux(ik,itheta,iphi) * 0.5_dp
                ELSE
                   bk(ik,itheta,iphi) = bk(ik,itheta,iphi) + &
                        intflux(ik,itheta,iphi)
                ENDIF

             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO NOWAIT
       !$OMP END PARALLEL
       ! End loop over coordinates

    ENDDO ! End time loop

    WRITE(*,*)
    WRITE(*,*) 'Time integral finished!'
    WRITE(*,*) '------------------------'
    WRITE(*,*)

    bk = bk * ZIMAGONE * dt

  END SUBROUTINE get_flux

  !***************************************!

  SUBROUTINE calculate_time_integrand( func_lm, funcp_lm, &
       lmax, mmax, afield, integrand )

    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)       :: func_lm(-mmax:, 0:)
    COMPLEX(dp), INTENT(IN)       :: funcp_lm(-mmax:, 0:)
    INTEGER, INTENT(IN)           :: lmax, mmax
    REAL(dp), INTENT(IN)          :: afield(:)
    COMPLEX(dp), INTENT(OUT)      :: integrand(:, :, :)

    COMPLEX(dp)                   :: term1, term2
    COMPLEX(dp)                   :: term3, term4
    INTEGER                       :: mmin
    REAL(dp)                      :: sqrtcoeff, sph_harm
    INTEGER                       :: il, ill, im
    INTEGER                       :: ik, itheta, iphi

    !---------------------------------------------------!

    mmin = -mmax
    sph_harm = 0.0_dp

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
    !$OMP& PRIVATE(term1,term2,term3,term4) &
    !$OMP& PRIVATE(sph_harm, sqrtcoeff)

    !$OMP DO COLLAPSE(3)
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             DO il = 0, lmax
                DO im = mmin, mmax
                   term1 = ZERO
                   term2 = ZERO
                   term3 = ZERO
                   term4 = ZERO
                   sph_harm = 0.0_dp
                   
                   term1 = 0.5_dp * (-ZIMAGONE)**il * &
                        (krb_ax(ik) * jlp(il,ik) - jl(il,ik)) * &
                        psi_lm(im,il)
                   
                   term2 =  - 0.5_dp * (-ZIMAGONE)**il * &
                        rb * jl(il,ik) * psip_lm(im,il)
                   
                   term3 = 0.5_dp * ZIMAGONE / SQRT(pi) * afield(3) * rb * &
                        psi_lm(im,il)
                   
                   DO ill = 0, lmax
                      sqrtcoeff = SQRT(REAL(2 * il + 1,dp)) / REAL(2 * ill + 1,dp)
                      sph_harm = 0.0_dp
                      IF(im.LT.0) THEN
                         sph_harm = (-1.0_dp)**ABS(im) * normfact(ABS(im),ill) * &
                              legenpl(itheta-1,ABS(im),ill) * &
                              EXP(ZIMAGONE * im * phi_ax(iphi))
                      ELSE
                         sph_harm = normfact(im,ill) * &
                              legenpl(itheta-1,im,ill) * &
                              EXP(ZIMAGONE * im * phi_ax(iphi))
                      ENDIF

                      term4 = term4 + (-ZIMAGONE)**ill * &
                           sqrtcoeff * jl(ill,ik) * &
                           clebsch(il, 1, ill, im, 0, im, factlog, maxfactlog) * &
                           clebsch(il, 1, ill, 0, 0, 0, factlog, maxfactlog) * &
                           sph_harm
                   ENDDO

                   term3 = term3 * term4

                   IF(im.LT.0) THEN
                      sph_harm = (-1.0_dp)**ABS(im) * &
                           normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ELSE
                      sph_harm = normfact(im,il) * legenpl(itheta-1,im,il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ENDIF
                   
                   !! Finally, add together the terms
                   integrand(ik,itheta,iphi) = integrand(ik,itheta,iphi) + &
                        (term1 + term2) * sph_harm + term3
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
  END SUBROUTINE calculate_time_integrand
  
  !***********************************************************!
  
END MODULE flux
