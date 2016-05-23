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
  USE spline_interp
  USE tools
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
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_th_weights(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_time(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_time_weights(:)
  REAL(dp), PUBLIC              :: dphi
  
  ! Private variables
  INTEGER                       :: numkpts
  REAL(dp)                      :: dt, dk, kmax, kmax_th
  REAL(dp)                      :: coulomb_exp_ener
  INTEGER                       :: ntime, ntime_gauss
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :), psip_sph(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_sph_v(:, :), psip_sph_v(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_lm(:, :), psip_lm(:, :)
  REAL(dp)                      :: rb, total_time
  REAL(dp), ALLOCATABLE         :: jl(:, :), jlp(:, :)
  REAL(dp), ALLOCATABLE         :: jy(:), jyp(:)
  REAL(dp), ALLOCATABLE         :: krb_ax(:)
  INTEGER                       :: lmax_total
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts
  LOGICAL                       :: gauge_trans_desired
  
  !-------------------------------------------------------------!
  !-------------------------------------------------------------!
  
CONTAINS
  
  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, &
       ddk, kmax_input, numgausstime, maxkpts, maxthetapts, maxphipts, &
       lmax_total, mmax, gauge_trans_on, coulomb_exp )
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: ddk
    REAL(dp), INTENT(IN)             :: kmax_input
    INTEGER, INTENT(IN)              :: numgausstime
    INTEGER, INTENT(OUT)             :: maxkpts
    INTEGER, INTENT(OUT)             :: maxthetapts
    INTEGER, INTENT(OUT)             :: maxphipts
    INTEGER, INTENT(OUT)             :: lmax_total
    INTEGER, INTENT(OUT)             :: mmax
    REAL(dp), INTENT(IN), OPTIONAL   :: coulomb_exp
    LOGICAL, INTENT(IN), OPTIONAL    :: gauge_trans_on
    
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
    dk = ddk

    ! Set how many points will be in the time quadrature
    ntime_gauss = numgausstime
    
    ! This factor shift the energy if we
    ! are dealing with fixed nuclei molecules
    coulomb_exp_ener = 0.0_dp
    IF(PRESENT(coulomb_exp)) &
         coulomb_exp_ener = coulomb_exp

    ! If this option is on, we will transform
    ! from length to velocity gauge the wavefunction
    IF(gauge_trans_on) THEN
       gauge_trans_desired = .TRUE.
    ELSE
       gauge_trans_desired = .FALSE.
    ENDIF
    
    kmax = kmax_input
    if(kmax.GT.kmax_th) THEN
       WRITE(*,*) 'Kmax desired large than the theoretical.'
       STOP
    ENDIF
    
    WRITE(*,'(A,F10.6)')  'Theoretical kmax:                ',kmax_th
    WRITE(*,'(A,F9.6)')  'Desired kmax:                     ',kmax
    WRITE(*,'(A,F9.6)')  'dk:                               ',dk
    WRITE(*,'(A,F0.6)')  'Time in the file (a.u.):          ', ntime * dt
    WRITE(*,'(A,F0.6)')  'Time in the file (fs):            ', ntime * dt * autime_fs
    WRITE(*,'(A,I4.4)')  'Points in time quadrature:        ', ntime_gauss
    !WRITE(*,'(A,F0.6)')  'Total time (a.u.):                 ', ntimetotal * dt
    !WRITE(*,'(A,F0.6)')  'Total time (fs):                   ', ntimetotal * dt * autime_fs
    
    IF(PRESENT(coulomb_exp)) &
         WRITE(*,'(A,F9.6)')  'Coulomb explosion energy (a.u.): ', coulomb_exp_ener
    WRITE(*,*)           '---------------------------------------------------------&
         &---------------------'
    WRITE(*,*)
    
    numkpts = INT( kmax / dk) + 1
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
    
    IF(gauge_trans_desired) THEN
       ALLOCATE(psi_sph_v(1:numthetapts,1:numphipts))
       ALLOCATE(psip_sph_v(1:numthetapts,1:numphipts))
    ENDIF

    ALLOCATE(psi_lm(mmin:mmax,0:lmax))
    ALLOCATE(psip_lm(mmin:mmax,0:lmax))
    
    ! Create spherical coordinates
    ALLOCATE(theta_ax(1:numthetapts),costheta_ax(1:numthetapts))
    ALLOCATE(gauss_th_weights(1:numthetapts))
    ALLOCATE(phi_ax(1:numphipts))
    
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta_ax, gauss_th_weights)
    theta_ax = ACOS(costheta_ax)
    
    dphi =  twopi / REAL(numphipts,dp)
    DO iphi = 1, numphipts
       phi_ax(iphi) = REAL(iphi,dp) * dphi
    ENDDO
    
    ! Create spherical harmonics
    CALL initialize_legendre_stuff(lmax, costheta_ax)
    CALL initialize_spherical_harmonics(lmax,mmax,costheta_ax,phi_ax)
    
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
    
    ! Create time axis for time integration
    ALLOCATE(gauss_time(1:ntime_gauss))
    ALLOCATE(gauss_time_weights(1:ntime_gauss))
    total_time = ntime * dt
    CALL get_gauss_stuff(0.0_dp,total_time,gauss_time,gauss_time_weights)
   
    
  END SUBROUTINE initialize_tsurff

  !**************************************!
  
  SUBROUTINE get_volkov_phase(afield,time,maxtime,phase, &
       coulomb_exp_ener)
    
    IMPLICIT NONE

    REAL(dp), INTENT(IN)           :: afield(:, :)
    REAL(dp), INTENT(IN)           :: time(:)
    REAL(dp), INTENT(IN)           :: maxtime
    COMPLEX(dp), INTENT(OUT)       :: phase(:, :, :)
    REAL(dp), INTENT(IN), OPTIONAL :: coulomb_exp_ener
    
    REAL(dp), ALLOCATABLE          :: time_ax(:)
    REAL(dp), ALLOCATABLE          :: weights(:)
    REAL(dp), ALLOCATABLE          :: newafield(:, :)
    REAL(dp), ALLOCATABLE          :: af2(:)
    REAL(dp)                       :: term1, term2
    REAL(dp)                       :: term3, term4
    REAL(dp)                       :: cou_ener, newterm
    REAL(dp)                       :: fdcoeffs(-2:2,2)
    REAL(dp)                       :: afield_deriv
    INTEGER                        :: nt, left
    INTEGER                        :: lower, upper
    INTEGER                        :: fdrule, rulepts
    INTEGER                        :: ik, itheta, iphi
    INTEGER                        :: itime, i

    !-------------------------------------------------!

    nt      = SIZE(time)
    fdrule  = 2
    rulepts = 2 * fdrule
    
    ALLOCATE(time_ax(1:ntime_gauss))
    ALLOCATE(weights(1:ntime_gauss))
    ALLOCATE(af2(1:nt))
    ALLOCATE(newafield(3,ntime_gauss))
    
    newafield = 0.0_dp
    phase = ZERO
    
    IF(PRESENT(coulomb_exp_ener)) THEN
       cou_ener = coulomb_exp_ener
    ELSE
       cou_ener = 0.0_dp
    ENDIF
    
    ! Create axis for integration
    CALL get_gauss_stuff(0.0_dp,maxtime,time_ax,weights)
    
    ! Get index of the last time value
    left = MAXLOC(time,DIM=1,MASK=(time.LE.maxtime))
    
    IF(left-fdrule.LT.LBOUND(time,DIM=1)) THEN
       lower = left
       upper = left + rulepts
    ELSEIF(left+fdrule+1.GT.UBOUND(time,DIM=1)) THEN
       lower = left - rulepts
       upper = left
    ELSE
       lower = left - fdrule
       upper = left + fdrule
    ENDIF
    ! Calculate fd weights for differentiation in z case             
    CALL fdweights(time(left),time(lower:upper),&
         rulepts+1,2,fdcoeffs)
    
    DO i = 1, 3
       af2 = 0.0_dp
       CALL make_derivative(afield(i,lower:upper),afield_deriv,fdrule, &
            dt,fdcoeffs(:,2))
       CALL get_scnd_derivatives_spline(time,afield(i,:),0.0_dp,afield_deriv,af2)
       
       DO itime = 1, ntime_gauss
          CALL spline_interpolation(time,afield(i,:),af2,time_ax(itime),newafield(i,itime))
       ENDDO
    ENDDO
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi,itime) &
    !$OMP& PRIVATE(term1,term2,term3,term4,newterm)
    
    !$OMP DO COLLAPSE(3)
    DO ik = 1, numkpts
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts

             newterm = 0.0_dp
             DO itime = 1, ntime_gauss
                term1 =  k_ax(ik) * k_ax(ik)
                
                term2 = k_ax(ik) * newafield(1,itime) * &
                     SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                
                term3 = k_ax(ik) * newafield(2,itime) * &
                     SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                
                term4 = k_ax(ik) * newafield(3,itime) * costheta_ax(itheta)
                
                newterm = newterm + &
                     (term1 + 2.0_dp * (term2 + term3 + term4)) * weights(itime)
                
             ENDDO
             phase(ik,itheta,iphi) = EXP( ZIMAGONE * (0.5_dp * newterm + &
                  cou_ener * maxtime) )
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    DEALLOCATE(time_ax,weights)
    DEALLOCATE(af2,newafield)
    
  END SUBROUTINE get_volkov_phase

  !****************************************************************!

  SUBROUTINE get_flux(filename, lmax, mmax, bk)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: lmax, mmax
    COMPLEX(dp), INTENT(OUT)     :: bk(:, :, :)

    REAL(dp), ALLOCATABLE        :: time(:), af2(:)
    REAL(dp), ALLOCATABLE        :: efield(:, :)
    REAL(dp), ALLOCATABLE        :: afield(:, :)
    REAL(dp), ALLOCATABLE        :: gauss_afield(:, :)
    COMPLEX(dp), ALLOCATABLE     :: intflux(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: volkov_phase(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: func_sph(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: funcp_sph(:, :, :)
    REAL(dp), ALLOCATABLE        :: func2_re(:), func2_imag(:)
    REAL(dp), ALLOCATABLE        :: funcp2_re(:), funcp2_imag(:)
    REAL(dp)                     :: new_func_re, new_func_imag
    INTEGER                      :: mmin
    INTEGER                      :: itime, igtime, il, im
    INTEGER                      :: ik, itheta, iphi, i

    !----------------------------------------------------!

    mmin = -mmax
    ! Allocate and set to zero
    ALLOCATE(intflux(1:numkpts,1:numthetapts,1:numphipts))
    ALLOCATE(volkov_phase(1:numkpts,1:numthetapts,1:numphipts))
    ALLOCATE(time(ntime))
    ALLOCATE(efield(3,ntime),afield(3,ntime))
    ALLOCATE(gauss_afield(3,ntime_gauss),af2(ntime))
    ALLOCATE(func_sph(1:numthetapts,1:numphipts,1:ntime))
    ALLOCATE(funcp_sph(1:numthetapts,1:numphipts,1:ntime))
    ALLOCATE(func2_re(1:ntime),func2_imag(1:ntime))
    ALLOCATE(funcp2_re(1:ntime),funcp2_imag(1:ntime))

    intflux = ZERO
    bk = ZERO
    volkov_phase = ZONE
    
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '------------------------'
    WRITE(*,*) 'Begin time integral...'
    
    ! First, get the time axis and vector potential
    DO itime = 1, ntime
       !WRITE(*,'(A,I6)') 'Reading loop number: ',itime
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            psi_sph, psip_sph, time(itime), efield(:,itime), afield(:,itime) )
    ENDDO
    
    ! Get afield on the Gauss time axis.
    DO i = 1, 3
       CALL get_scnd_derivatives_spline(time,afield(i,:),0.0_dp, 0.0_dp, af2)
       
       DO itime = 1, ntime_gauss
          CALL spline_interpolation(time,afield(i,:),af2, &
               gauss_time(itime),gauss_afield(i,itime))
       ENDDO
    ENDDO
    
    !-------------------------------------------!
    
    DO itime = 1, ntime
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            func_sph(:,:,itime), funcp_sph(:,:,itime), time(itime), &
            efield(:,itime), afield(:,itime) )
    ENDDO
    
    ! We begin the time loop!!
    DO igtime  = 1, ntime_gauss
       
       !WRITE(*,*) 'Loop number: ',igtime
       ! Get function (for all times).
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             
             func2_re    = 0.0_dp
             func2_imag  = 0.0_dp
             funcp2_re   = 0.0_dp
             funcp2_imag = 0.0_dp
             
             ! Get second derivatives, for real and imaginary parts
             CALL get_scnd_derivatives_spline(time,REAL(func_sph(itheta,iphi,:)), &
                  0.0_dp,0.0_dp,func2_re)
             CALL get_scnd_derivatives_spline(time,AIMAG(func_sph(itheta,iphi,:)), &
                  0.0_dp,0.0_dp,func2_imag)
             
             ! Get second derivatives, for real and imaginary parts
             CALL get_scnd_derivatives_spline(time,REAL(funcp_sph(itheta,iphi,:)), &
                  0.0_dp,0.0_dp,funcp2_re)
             CALL get_scnd_derivatives_spline(time,AIMAG(funcp_sph(itheta,iphi,:)), &
                  0.0_dp,0.0_dp,funcp2_imag)
             
             CALL spline_interpolation(time,REAL(func_sph(itheta,iphi,:)),&
                  func2_re, gauss_time(igtime),new_func_re)
             CALL spline_interpolation(time,AIMAG(func_sph(itheta,iphi,:)),&
                  func2_imag, gauss_time(igtime),new_func_imag)
             
             psi_sph(itheta,iphi) = CMPLX(new_func_re,new_func_imag)
             
             CALL spline_interpolation(time,REAL(funcp_sph(itheta,iphi,:)), &
                  funcp2_re, gauss_time(igtime),new_func_re)
             CALL spline_interpolation(time,AIMAG(funcp_sph(itheta,iphi,:)), &
                  funcp2_imag, gauss_time(igtime),new_func_imag)
             
             psip_sph(itheta,iphi) = CMPLX(new_func_re,new_func_imag)
             
          ENDDO
       ENDDO
       
       IF(gauge_trans_desired) THEN
          CALL gauge_transform_l2v(psi_sph, psi_sph_v, gauss_afield(:,igtime))
          CALL gauge_transform_l2v(psip_sph, psip_sph_v, gauss_afield(:,igtime))
          psi_sph  = psi_sph_v
          psip_sph = psip_sph_v 
       ENDIF
       
       ! Decompose into spherical harmonics:
       ! The wavefunction
       CALL make_sht(psi_sph, lmax, mmax, gauss_th_weights, &
            psi_lm)
       ! His first derivative
       CALL make_sht(psip_sph, lmax, mmax, gauss_th_weights, &
            psip_lm)
       ! Divide by rb
       psi_lm = psi_lm * rb
       psip_lm = psip_lm * rb + psi_lm / rb
       
       ! Calculate flux
       intflux = ZERO
       CALL calculate_time_integrand(psi_lm, psip_lm, &
            lmax, mmax, gauss_afield(:,igtime), intflux)
       
       ! Get the Volkov phase (one per time step)
       CALL get_volkov_phase(afield, time, gauss_time(igtime), &
            volkov_phase, coulomb_exp_ener )
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi)
       
       !$OMP DO COLLAPSE(3)
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ik = 1, numkpts
                
                intflux(ik,itheta,iphi) = intflux(ik,itheta,iphi) * &
                     volkov_phase(ik,itheta,iphi)
                
                
                bk(ik,itheta,iphi) = bk(ik,itheta,iphi) + &
                     intflux(ik,itheta,iphi) * gauss_time_weights(igtime)
                
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
    
    bk = bk * SQRT(2.0_dp / pi) * ZIMAGONE
    
    ! Deallocate like no other
    DEALLOCATE(time,efield,afield)
    DEALLOCATE(gauss_time,gauss_afield,af2)
    DEALLOCATE(intflux,volkov_phase)
    DEALLOCATE(func_sph,funcp_sph)
    DEALLOCATE(func2_re,func2_imag)
    DEALLOCATE(funcp2_re,funcp2_imag)
    
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
    COMPLEX(dp)                   :: term3
    INTEGER                       :: mmin
    COMPLEX(dp)                   :: suma
    INTEGER                       :: il, im
    INTEGER                       :: ill, imm
    INTEGER                       :: ik, itheta, iphi
    INTEGER                       :: jtheta, jphi
    
    !---------------------------------------------------!

    mmin = -mmax

    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
    !$OMP& PRIVATE(term1,term2,term3,jtheta,jphi,suma) &
    !$OMP& PRIVATE(il,im,ill,imm)
    
    !$OMP DO COLLAPSE(3)
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             DO il = 0, lmax
                DO im = mmin, mmax
                   term1 = ZERO
                   term2 = ZERO
                   term3 = ZERO
                   
                   term1 = 0.5_dp * (-ZIMAGONE)**il * &
                        (krb_ax(ik) * jlp(il,ik) - jl(il,ik)) * &
                        func_lm(im,il)
                   
                   term2 =  - 0.5_dp * (-ZIMAGONE)**il * &
                        rb * jl(il,ik) * funcp_lm(im,il)
                   
                   
                   DO ill = 0, lmax
                      DO imm = mmin, mmax
                         
                         suma = ZERO
                         DO jphi = 1, numphipts
                            DO jtheta = 1, numthetapts
                               suma = suma + CONJG(sph_harmonics(jphi,jtheta,im,il)) * &
                                    SIN(theta_ax(jtheta)) * COS(phi_ax(jphi)) * &
                                    sph_harmonics(jphi,jtheta,imm,ill) * &
                                    gauss_th_weights(jtheta)
                            ENDDO
                         ENDDO
                         IF(numphipts.EQ.1) THEN
                            suma = suma * twopi
                         ELSE
                            suma = suma * dphi
                         ENDIF
                         
                         term3 = term3 + suma * afield(1) * psi_lm(imm,ill)
                         
                         suma = ZERO
                         DO jphi = 1, numphipts
                            DO jtheta = 1, numthetapts
                               suma = suma + CONJG(sph_harmonics(jphi,jtheta,im,il)) * &
                                    costheta_ax(jtheta) * sph_harmonics(jphi,jtheta,imm,il) * &
                                    gauss_th_weights(jtheta)
                            ENDDO
                         ENDDO
                         IF(numphipts.EQ.1) THEN
                            suma = suma * twopi
                         ELSE
                            suma = suma * dphi
                         ENDIF
                         
                         term3 = term3 + suma * afield(3) * psi_lm(imm,ill)
                         
                      ENDDO
                   ENDDO
                   
                   term3 = term3 * (-ZIMAGONE)**il * (-ZIMAGONE) * &
                        jl(il,ik) * rb
                   
                   !! Finally, add together the terms
                   integrand(ik,itheta,iphi) = integrand(ik,itheta,iphi) + &
                        (term1 + term2 + term3) * sph_harmonics(iphi,itheta,im,il)
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
   
  END SUBROUTINE calculate_time_integrand

  !***********************************************************!
  !***********************************************************!
  
  
  !   SUBROUTINE calculate_time_integrand( func_lm, funcp_lm, &
  !      lmax, mmax, afield, integrand )
  
  !   IMPLICIT NONE
  
  !   COMPLEX(dp), INTENT(IN)       :: func_lm(-mmax:, 0:)
  !   COMPLEX(dp), INTENT(IN)       :: funcp_lm(-mmax:, 0:)
  !   INTEGER, INTENT(IN)           :: lmax, mmax
  !   REAL(dp), INTENT(IN)          :: afield(:)
  !   COMPLEX(dp), INTENT(OUT)      :: integrand(:, :, :)
  
  !   COMPLEX(dp)                   :: term1, term2
  !   COMPLEX(dp)                   :: term3, term4
  !   INTEGER                       :: mmin
  !   REAL(dp)                      :: sqrtcoeff
  !   COMPLEX(dp)                   :: sph_harm
  !   INTEGER                       :: il, ill, im
  !   INTEGER                       :: ik, itheta, iphi

  !   !---------------------------------------------------!

  !   mmin = -mmax
  !   sph_harm = ZERO

  !   !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
  !   !$OMP& PRIVATE(term1,term2,term3,term4) &
  !   !$OMP& PRIVATE(sph_harm, sqrtcoeff)

  !   !$OMP DO COLLAPSE(3)
  !   DO iphi = 1, numphipts
  !      DO itheta = 1, numthetapts
  !         DO ik = 1, numkpts
  !            DO il = 0, lmax
  !               DO im = mmin, mmax
  !                  term1 = ZERO
  !                  term2 = ZERO
  !                  term3 = ZERO
  !                  term4 = ZERO
  !                  sph_harm = ZERO
                   
  !                  term1 = 0.5_dp * (-ZIMAGONE)**il * &
  !                       (krb_ax(ik) * jlp(il,ik) + jl(il,ik)) * &
  !                       func_lm(im,il)
                   
  !                  term2 =  - 0.5_dp * (-ZIMAGONE)**il * &
  !                       rb * jl(il,ik) * funcp_lm(im,il)
                   
  !                  term3 = - 0.5_dp * ZIMAGONE / SQRT(pi) * afield(3) * rb * &
  !                       func_lm(im,il)
                   
  !                  DO ill = 0, lmax
  !                     sqrtcoeff = SQRT(REAL(2 * il + 1,dp)) / REAL(2 * ill + 1,dp)
  !                     sph_harm = ZERO
  !                     IF(im.LT.0) THEN
  !                        sph_harm = (-1.0_dp)**ABS(im) * normfact(ABS(im),ill) * &
  !                             legenpl(itheta-1,ABS(im),ill) * &
  !                             EXP(ZIMAGONE * im * phi_ax(iphi))
  !                     ELSE
  !                        sph_harm = normfact(im,ill) * &
  !                             legenpl(itheta-1,im,ill) * &
  !                             EXP(ZIMAGONE * im * phi_ax(iphi))
  !                     ENDIF
                      
  !                     term4 = term4 + (-ZIMAGONE)**ill * &
  !                          sqrtcoeff * jl(ill,ik) * &
  !                          clebsch(il, 1, ill, im, 0, im, factlog, maxfactlog) * &
  !                          clebsch(il, 1, ill, 0, 0, 0, factlog, maxfactlog) * &
  !                          sph_harm
  !                  ENDDO
                   
  !                  term3 = term3 * term4
                   
  !                  IF(im.LT.0) THEN
  !                     sph_harm = (-1.0_dp)**ABS(im) * &
  !                          normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) * &
  !                          EXP(ZIMAGONE * im * phi_ax(iphi))
  !                  ELSE
  !                     sph_harm = normfact(im,il) * legenpl(itheta-1,im,il) * &
  !                          EXP(ZIMAGONE * im * phi_ax(iphi))
  !                  ENDIF
                   
  !                  !! Finally, add together the terms
  !                  integrand(ik,itheta,iphi) = integrand(ik,itheta,iphi) + &
  !                       (term1 + term2) * sph_harm + term3
                   
  !               ENDDO
  !            ENDDO
  !         ENDDO
  !      ENDDO
  !   ENDDO
  !   !$OMP END DO NOWAIT
  !   !$OMP END PARALLEL
    
    
  ! END SUBROUTINE calculate_time_integrand

  !***********************************************************!
  
  SUBROUTINE gauge_transform_l2v(psi_length,psi_vel,afield)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: psi_length(:, :)
    COMPLEX(dp), INTENT(OUT)  :: psi_vel(:, :)
    REAL(dp), INTENT(IN)      :: afield(:)
    
    INTEGER                   :: dims(2)
    INTEGER                   :: numthetapts
    INTEGER                   :: numphipts
    INTEGER                   :: itheta, iphi
    
    !--------------------------------------!
    
    dims = shape(psi_length)
    numthetapts = dims(1)
    numphipts   = dims(2)
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          psi_vel(itheta,iphi) = psi_length(itheta,iphi) * &
               EXP(- ZIMAGONE * ( &
               afield(1) * rb * SIN(theta_ax(itheta)) * COS(phi_ax(iphi)) + &
               afield(2) * rb * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi)) + &
               afield(3) * rb * costheta_ax(itheta) ))
       ENDDO
    ENDDO
    
  END SUBROUTINE gauge_transform_l2v
  
  !***********************************************************!
  !***********************************************************!
  
END MODULE flux
