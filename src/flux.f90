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
  !USE bessel
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
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_phi_weights(:)
  
  ! Private variables
  INTEGER                       :: numkpts
  REAL(dp)                      :: dt, dk, kmax, kmax_th
  INTEGER                       :: ntime
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :), psip_sph(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_sph_v(:, :), psip_sph_v(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_lm(:, :), psip_lm(:, :)
  REAL(dp)                      :: rb, total_time
  REAL(dp), ALLOCATABLE         :: jl(:, :), jlp(:, :)
  REAL(dp), ALLOCATABLE         :: yl(:), ylp(:)
  REAL(dp), ALLOCATABLE         :: krb_ax(:)
  REAL(dp), ALLOCATABLE         :: xcoupling(:, :, :, :)
  REAL(dp), ALLOCATABLE         :: ycoupling(:, :, :, :)
  REAL(dp), ALLOCATABLE         :: zcoupling(:, :, :, :)
  INTEGER                       :: lmax_total
  INTEGER                       :: numpts, numrpts
  INTEGER                       :: numthetapts, numphipts
  LOGICAL                       :: gauge_trans_desired
  
  !-------------------------------------------------------------!
  !-------------------------------------------------------------!
  
CONTAINS
  
  SUBROUTINE initialize_tsurff(filename, radb, lmax_desired, &
       ddk, kmax_input, maxkpts, maxthetapts, maxphipts, &
       lmax_total, gauge_trans_on )
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    INTEGER, INTENT(IN)              :: lmax_desired
    REAL(dp), INTENT(IN)             :: radb
    REAL(dp), INTENT(IN)             :: ddk
    REAL(dp), INTENT(IN)             :: kmax_input
    INTEGER, INTENT(OUT)             :: maxkpts
    INTEGER, INTENT(OUT)             :: maxthetapts
    INTEGER, INTENT(OUT)             :: maxphipts
    INTEGER, INTENT(OUT)             :: lmax_total
    LOGICAL, INTENT(IN), OPTIONAL    :: gauge_trans_on

    COMPLEX(dp)                      :: suma
    COMPLEX(dp)                      :: sumb    
    REAL(dp)                         :: dphi 
    INTEGER                          :: lmax
    INTEGER                          :: iphi, ik
    INTEGER                          :: il, im, ill, imm
    INTEGER                          :: jtheta, jphi
    INTEGER                          :: ifail
    
    !-------------------------------------------------!

    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename, ntime, dt, &
         numthetapts, numphipts, lmax_total)

    maxthetapts = numthetapts
    maxphipts = numphipts
    
    ! Create momentum axis
    kmax_th = twopi / dt
    dk = ddk
    
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
    
    WRITE(*,*)
    
    numkpts = INT( kmax / dk) + 1
    maxkpts = numkpts
    
    ALLOCATE(k_ax(1:numkpts))
    
    DO ik = 1, numkpts
       k_ax(ik) = REAL(ik-1,dp) * dk
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
    ENDIF
    
    ! Allocate wavefunction arrays
    ALLOCATE(psi_sph(1:numthetapts,1:numphipts))
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts))
    
    IF(gauge_trans_desired) THEN
       ALLOCATE(psi_sph_v(1:numthetapts,1:numphipts))
       ALLOCATE(psip_sph_v(1:numthetapts,1:numphipts))
    ENDIF
    
    ALLOCATE(psi_lm(-lmax:lmax,0:lmax))
    ALLOCATE(psip_lm(-lmax:lmax,0:lmax))
    psi_lm  = ZERO
    psip_lm = ZERO
    
    ! Create spherical coordinates
    ALLOCATE(theta_ax(1:numthetapts),costheta_ax(1:numthetapts))
    ALLOCATE(gauss_th_weights(1:numthetapts))
    ALLOCATE(phi_ax(1:numphipts))
    ALLOCATE(gauss_phi_weights(1:numphipts))
    
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta_ax, gauss_th_weights)
    theta_ax = ACOS(costheta_ax)
    
    dphi =  twopi / REAL(numphipts,dp)
    DO iphi = 1, numphipts
       phi_ax(iphi) = REAL(iphi-1,dp) * dphi
       gauss_phi_weights(iphi) = dphi
    ENDDO
    
    ! Create spherical harmonics
    CALL initialize_legendre_stuff(lmax, costheta_ax)
    CALL initialize_spherical_harmonics(lmax,costheta_ax,phi_ax)
    
    IF(numphipts.EQ.1) THEN
       ALLOCATE(xcoupling(0:0,0:lmax,0:0,0:lmax))
       ALLOCATE(ycoupling(0:0,0:lmax,0:0,0:lmax))
       ALLOCATE(zcoupling(0:0,0:lmax,0:0,0:lmax))
       xcoupling = 0.0_dp
       ycoupling = 0.0_dp
       zcoupling = 0.0_dp     
       
       ! Calculate matrix elements for laser-matter couling in
       ! z direction laser polarisation.
       
       DO ill = 0, lmax
          DO il = 0, lmax
             
             IF (ABS(il - ill).EQ.1) THEN
                ! Calculate matrix elements between spherical harmonics
                suma = ZERO
                DO jphi = 1, numphipts
                   DO jtheta = 1, numthetapts
                      suma = suma + &
                           CONJG(sph_harmonics(jphi,jtheta,0,il)) * &
                           sph_harmonics(jphi,jtheta,0,ill) * &
                           costheta_ax(jtheta) * &
                           gauss_th_weights(jtheta) * &
                           gauss_phi_weights(jphi)
                   ENDDO
                ENDDO
                zcoupling(0,il,0,ill) = suma
             ENDIF
             
          ENDDO
       ENDDO
       
    ELSE
       ALLOCATE(xcoupling(-lmax:lmax,0:lmax,-lmax:lmax,0:lmax))
       ALLOCATE(ycoupling(-lmax:lmax,0:lmax,-lmax:lmax,0:lmax))
       ALLOCATE(zcoupling(-lmax:lmax,0:lmax,-lmax:lmax,0:lmax))
       xcoupling = 0.0_dp
       ycoupling = 0.0_dp
       zcoupling = 0.0_dp     
       
       ! Calculate matrix elements for laser-matter couling in
       ! x and y direction laser polarisation.
       
       DO ill = 0, lmax
          DO imm = -ill, ill
             
             DO il = 0, lmax
                DO im = -il, il

                   IF (ABS(il - ill).EQ.1 .AND. ABS(im - imm).EQ.1) THEN
                      
                      ! Calculate matrix elements between spherical harmonics
                      suma = ZERO
                      sumb = ZERO
                      
                      DO jphi = 1, numphipts
                         DO jtheta = 1, numthetapts
                            
                            suma = suma + &
                                 CONJG(sph_harmonics(jphi,jtheta,im,il)) * &
                                 sph_harmonics(jphi,jtheta,imm,ill) * &
                                 SIN(theta_ax(jtheta)) * COS(phi_ax(jphi)) * &
                                 gauss_th_weights(jtheta) * &
                                 gauss_phi_weights(jphi)

                            sumb = sumb + &
                                 CONJG(sph_harmonics(jphi,jtheta,im,il)) * &
                                 sph_harmonics(jphi,jtheta,imm,ill) * &
                                 SIN(theta_ax(jtheta)) * SIN(phi_ax(jphi)) * &
                                 gauss_th_weights(jtheta) * &
                                 gauss_phi_weights(jphi)
                            
                         ENDDO
                      ENDDO
                      
                      xcoupling(im,il,imm,ill) = suma
                      ycoupling(im,il,imm,ill) = AIMAG(sumb)
                      
                   ENDIF
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       
       ! Calculate matrix elements for laser-matter couling in
       ! z direction laser polarisation.
       
       DO ill = 0, lmax
          DO imm = -ill, ill
             DO il = 0, lmax
                DO im = -il, il
                   
                   IF (ABS(il - ill).EQ.1 .AND. ABS(im - imm).EQ.0) THEN
                      ! Calculate matrix elements between spherical harmonics
                      
                      suma = ZERO
                      
                      DO jphi = 1, numphipts
                         DO jtheta = 1, numthetapts
                            
                            suma = suma + &
                                 CONJG(sph_harmonics(jphi,jtheta,im,il)) * &
                                 sph_harmonics(jphi,jtheta,imm,ill) * &
                                 costheta_ax(jtheta) * &
                                 gauss_th_weights(jtheta) * &
                                 gauss_phi_weights(jphi)
                            
                         ENDDO
                      ENDDO
                      
                      zcoupling(im,il,imm,ill) = suma
                      
                   ENDIF
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       
    ENDIF
    
    ! Initialize spherical Bessel functions
    ALLOCATE(jl(0:lmax,1:numkpts))
    ALLOCATE(jlp(0:lmax,1:numkpts))
    ALLOCATE(yl(0:lmax))
    ALLOCATE(ylp(0:lmax))
    ALLOCATE(krb_ax(1:numkpts))
    
    krb_ax = k_ax * rb
    jl  = 0.0_dp
    jlp = 0.0_dp
    yl  = 0.0_dp
    ylp = 0.0_dp
    
    DO ik = 1, numkpts
       CALL sbesjy(krb_ax(ik),lmax,jl(:,ik),yl,jlp(:,ik),ylp,ifail)
    ENDDO
    
  END SUBROUTINE initialize_tsurff

  !**************************************!
  
  SUBROUTINE get_volkov_phase(afield, time, maxtime, phase)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: afield(:, :)
    REAL(dp), INTENT(IN)           :: time(:)
    REAL(dp), INTENT(IN)           :: maxtime
    COMPLEX(dp), INTENT(OUT)       :: phase(:, :, :)
    
    REAL(dp)                       :: integralAx
    REAL(dp)                       :: integralAy
    REAL(dp)                       :: integralAz
    INTEGER                        :: ik, itheta, iphi
    INTEGER                        :: itime, numtimesteps
    
    !-------------------------------------------------!
    
    phase = ZERO
    
    numtimesteps = MAXLOC(time,DIM=1,MASK=(time.LE.maxtime))
    
    integralAx = 0.0_dp
    integralAy = 0.0_dp
    integralAz = 0.0_dp
    
    DO itime = 1, numtimesteps
       IF( (itime.EQ.1) .OR. (itime.EQ.numtimesteps)) THEN
          integralAx = integralAx + &
               afield(1,itime) * 0.5_dp
          
          integralAy = integralAy + &
               afield(2,itime) * 0.5_dp
          
          integralAz = integralAz + &
               afield(3,itime) * 0.5_dp
       ELSE
          integralAx = integralAx + &
               afield(1,itime)
          
          integralAy = integralAy + &
               afield(2,itime)
          
          integralAz = integralAz + &
               afield(3,itime)
       ENDIF
       
    ENDDO
    
    integralAx = integralAx * dt
    integralAy = integralAy * dt
    integralAz = integralAz * dt
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi)
    
    !$OMP DO COLLAPSE(3)    
    DO ik = 1, numkpts
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             phase(ik,itheta,iphi) = EXP( ZIMAGONE * &
                  ( 0.5_dp * k_ax(ik) * k_ax(ik) * maxtime ) + &
                  ( k_ax(ik) * integralAx * SIN(theta_ax(itheta)) * COS(phi_ax(iphi)) ) + &
                  ( k_ax(ik) * integralAy * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi)) ) + &
                  ( k_ax(ik) * integralAz * costheta_ax(itheta) ) )
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    
  END SUBROUTINE get_volkov_phase
  
  !****************************************************************!

  SUBROUTINE get_flux(filename, lmax, bk)

    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN)          :: lmax
    COMPLEX(dp), INTENT(OUT)     :: bk(:, :, :)

    REAL(dp)                     :: popsurface
    REAL(dp), ALLOCATABLE        :: time(:)
    REAL(dp), ALLOCATABLE        :: efield(:, :)
    REAL(dp), ALLOCATABLE        :: afield(:, :)
    COMPLEX(dp), ALLOCATABLE     :: intflux(:, :, :)
    COMPLEX(dp), ALLOCATABLE     :: volkov_phase(:, :, :)
    INTEGER                      :: itime, il, im
    INTEGER                      :: ik, itheta, iphi
     
    !----------------------------------------------------!
    
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
    WRITE(*,*) '------------------------'
    WRITE(*,*) 'Begin time integral...'
    
    ! We begin the time loop!!
    DO itime = 1, ntime
       
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            psi_sph, psip_sph, time(itime), efield(:,itime), afield(:,itime) )
       
       popsurface = 0.0_dp
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts 
             
             popsurface = popsurface +                        &
                  REAL(CONJG(psi_sph(itheta, iphi)) *         &
                  psi_sph(itheta, iphi), dp) *                &
                  gauss_th_weights(itheta) * gauss_phi_weights(iphi)
             
          ENDDO
       ENDDO
       
       WRITE(*,'(A14, I6,A14, F20.16)' ) ' Loop number =', itime,  &
            ' Popsurface =', popsurface
       
       IF(gauge_trans_desired) THEN
          CALL gauge_transform_l2v(psi_sph, psi_sph_v, afield, itime)
          CALL gauge_transform_l2v(psip_sph, psip_sph_v, afield, itime)
          psi_sph  = psi_sph_v
          psip_sph = psip_sph_v 
       ENDIF
       
       ! Decompose into spherical harmonics:
       ! The wavefunction
       CALL make_sht(psi_sph, lmax, gauss_th_weights, &
            gauss_phi_weights, psi_lm)
       ! His first derivative
       CALL make_sht(psip_sph, lmax, gauss_th_weights, &
            gauss_phi_weights, psip_lm)
       
       ! Calculate flux
       intflux = ZERO
       CALL calculate_time_integrand(psi_lm, psip_lm, &
            lmax, afield(:,itime), intflux)
       
       ! Get the Volkov phase (one per time step)
       CALL get_volkov_phase(afield, time, time(itime), &
            volkov_phase )

       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi)
       
       !$OMP DO COLLAPSE(3)       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ik = 1, numkpts
                
                intflux(ik,itheta,iphi) = intflux(ik,itheta,iphi) * &
                     volkov_phase(ik,itheta,iphi)
                                
                IF((itime.EQ.1) .OR. (itime.EQ.ntime)) THEN
                   bk(ik,itheta,iphi) = bk(ik,itheta,iphi) + &
                        intflux(ik,itheta,iphi) * 0.5
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
    
    WRITE(*,*) '-------------------------'
    WRITE(*,*) '                         '
    WRITE(*,*) '    End time integral    '
    WRITE(*,*) '                         '
    WRITE(*,*) '-------------------------'
    
    bk = bk * ZIMAGONE * dt * SQRT(2.0_dp / pi) * rb * rb
    
    ! Deallocate like no other
    DEALLOCATE(time,efield,afield)
    DEALLOCATE(intflux,volkov_phase)
    
  END SUBROUTINE get_flux
  
  !***************************************!
  
  SUBROUTINE calculate_time_integrand( func_lm, funcp_lm, &
       lmax, afield, integrand )
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: func_lm(-lmax:, 0:)
    COMPLEX(dp), INTENT(IN)       :: funcp_lm(-lmax:, 0:)
    INTEGER, INTENT(IN)           :: lmax
    REAL(dp), INTENT(IN)          :: afield(:)
    COMPLEX(dp), INTENT(OUT)      :: integrand(:, :, :)
    
    COMPLEX(dp)                   :: term1, term2
    COMPLEX(dp)                   :: fieldterm, suma
    INTEGER                       :: il, im
    INTEGER                       :: ill, imm
    INTEGER                       :: ik, itheta, iphi
    
    !---------------------------------------------------!
    
    IF(numphipts.EQ.1) THEN
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
       !$OMP& PRIVATE(term1,term2,fieldterm,suma) &
       !$OMP& PRIVATE(il,ill)
       
       !$OMP DO COLLAPSE(3)
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ik = 1, numkpts
                
                suma = ZERO
                
                DO il = 0, lmax
                   
                   term1 = func_lm(0,il) * &
                        sph_harmonics(iphi,itheta,0,il)
                   term2 = funcp_lm(0,il) * &
                        sph_harmonics(iphi,itheta,0,il)
                   
                   fieldterm = ZERO
                   DO ill = 0, lmax
                      
                      IF (ABS(il - ill).EQ.1 ) THEN   
                         fieldterm = fieldterm + &
                              zcoupling(0,il,0,ill) * &
                              afield(3) * func_lm(0,ill)
                      ENDIF
                   ENDDO
                   
                   fieldterm = ZIMAGONE * sph_harmonics(iphi,itheta,0,il) * fieldterm
                   
                   suma = suma + ((-ZIMAGONE)**il) * &
                        (  0.5_dp * k_ax(ik) * jlp(il,ik) * term1 - &
                        0.5_dp * jl(il,ik) * term2 -&
                        jl(il,ik) * fieldterm )
                   
                ENDDO
                integrand(ik,itheta,iphi) = suma
             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO NOWAIT
       !$OMP END PARALLEL
       
    ELSE
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,itheta,iphi) &
       !$OMP& PRIVATE(term1,term2,fieldterm,suma) &
       !$OMP& PRIVATE(il,im,ill,imm)
       
       !$OMP DO COLLAPSE(3)
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ik = 1, numkpts
                
                suma = ZERO
                
                DO il = 0, lmax
                   
                   term1 = ZERO
                   term2 = ZERO
                   
                   DO im = -il, il
                      
                      term1 = term1 + func_lm(im,il) * &
                           sph_harmonics(iphi,itheta,im,il)
                      term2 = term2 + funcp_lm(im,il) * &
                           sph_harmonics(iphi,itheta,im,il)
                      
                      fieldterm = ZERO
                      DO ill = 0, lmax
                         DO imm = -ill, ill
                            
                            IF (ABS(il - ill).EQ.1 .AND.   &
                                 ABS(im - imm).EQ.1) THEN
                               
                               fieldterm = fieldterm +  &
                                    (xcoupling(im, il, imm, ill) * &
                                    afield(1) + &
                                    ycoupling(im, il, imm, ill) * &
                                    ZIMAGONE * afield(2)) * &
                                    func_lm(imm, ill)
                               
                            ENDIF
                            
                            IF (ABS(il - ill).EQ.1 .AND.  &
                                 ABS(im - imm).EQ.0) THEN
                               
                               fieldterm = fieldterm + &
                                    zcoupling(im, il, imm, ill) * &
                                    afield(3) * func_lm(imm, ill)
                               
                            ENDIF
                            
                         ENDDO
                      ENDDO
                      
                      fieldterm = ZIMAGONE * sph_harmonics(iphi,itheta,im,il) * fieldterm
                      
                   ENDDO
                   
                   suma = suma + ( (-ZIMAGONE)**il) * &
                        (  0.5_dp * k_ax(ik) * jlp(il,ik) * term1 - &
                        0.5_dp * jl(il,ik) * term2 - &
                        jl(il,ik) * fieldterm )                   
                   
                ENDDO
                integrand(ik,itheta,iphi) = suma
             ENDDO
          ENDDO
       ENDDO
       !$OMP END DO NOWAIT
       !$OMP END PARALLEL
       
    ENDIF
    
  END SUBROUTINE calculate_time_integrand
  
  !***********************************************************!
  
  SUBROUTINE gauge_transform_l2v(psi_length,psi_vel,afield, currtime )
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: psi_length(:, :)
    COMPLEX(dp), INTENT(OUT)  :: psi_vel(:, :)
    REAL(dp), INTENT(IN)      :: afield(:, :)
    INTEGER, INTENT(IN)       :: currtime
    
    INTEGER                   :: itheta, iphi, itime
    REAL(dp)                  :: adotr, tfac, asq, expfac
    
    !--------------------------------------!
    
    DO itime = 1, currtime
       
       asq  = afield(1, itime) * afield(1, itime) +  &
            afield(2, itime) * afield(2, itime) + &
            afield(3, itime) * afield(3, itime)
       
       asq = 0.5_dp * asq
       
       tfac = 2.0_dp
       
       IF ((itime.EQ.1) .OR. (itime.EQ.currtime)) THEN
          
          tfac = 1.0_dp
          
       ENDIF
       
       expfac = expfac + asq * tfac
       
    ENDDO
    
    expfac = 0.5_dp * expfac
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(itheta,iphi, adotr)
    
    !$OMP DO COLLAPSE(2)
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          
          adotr = afield(1, currtime) * SIN(theta_ax(itheta)) *              &
               COS(phi_ax(iphi)) +                                        &
               afield(2, currtime) * SIN(theta_ax(itheta)) *              &
               SIN(phi_ax(iphi)) +                                        &
               afield(3, currtime) * costheta_ax(itheta) 
          
          adotr = adotr * rb
          
          psi_vel(itheta, iphi) = psi_length(itheta, iphi) *                 &
               EXP(ZIMAGONE * (expfac - adotr))
          
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
  END SUBROUTINE gauge_transform_l2v
  
  !***********************************************************!
  !***********************************************************!
  
END MODULE flux
