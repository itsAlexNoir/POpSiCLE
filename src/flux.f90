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
  
  USE comm_surff
  USE constants_pop
  USE gaussleg
  !USE bessel
  USE tools
  USE sharmonics
  USE io_surface
  
  IMPLICIT NONE
  
  PRIVATE
  
  ! Public routines
  
  PUBLIC        :: initialize_tsurff
  PUBLIC        :: get_volkov_phase
  PUBLIC        :: get_flux
  PUBLIC        :: calculate_time_integrand
  PUBLIC        :: make_kpoints
  
  !-------------------------------------------------------------!
  
  ! Public variables
  
  ! Private variables
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: probk3D(:, :, :)
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: bktime(:, :, :)
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: bktimeold(:, :, :)
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: bk(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: bkaverage(:, :, :)
  
  COMPLEX(dp), ALLOCATABLE         :: sph_harmonics(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE         :: sph_harmonicsk(:, :, :, :)
  
  COMPLEX(dp), ALLOCATABLE         :: psi_sph(:, :)
  COMPLEX(dp), ALLOCATABLE         :: psip_sph(:, :)
  
  COMPLEX(dp), ALLOCATABLE         :: psi_lm(:, :)
  COMPLEX(dp), ALLOCATABLE         :: psip_lm(:, :)
  
  REAL(dp), ALLOCATABLE	           :: krspts(:)
  REAL(dp), ALLOCATABLE            :: jl(:, :)
  REAL(dp), ALLOCATABLE            :: jlp(:, :)
  REAL(dp), ALLOCATABLE            :: yl(:)
  REAL(dp), ALLOCATABLE            :: ylp(:)
  
  REAL(dp), ALLOCATABLE            :: xcoupling(:, :, :, :)
  REAL(dp), ALLOCATABLE            :: ycoupling(:, :, :, :)
  REAL(dp), ALLOCATABLE            :: zcoupling(:, :, :, :)
  
  REAL(dp), ALLOCATABLE	           :: efield(:, :)
  REAL(dp), ALLOCATABLE	           :: afield(:, :)
  REAL(dp), ALLOCATABLE	           :: quiver(:, :)
  REAL(dp), ALLOCATABLE	           :: asqintegral(:)
  
  
CONTAINS
  
  
  SUBROUTINE initialize_tsurff( filename )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    CHARACTER(LEN = *), INTENT(IN)   :: filename

    INTEGER                          :: ik, itheta, iphi
    INTEGER                          :: il, im, ill, imm
    INTEGER                          :: jmaxorder, ifail
    INTEGER                          :: itime
    REAL(dp)                         :: t1, a1, a2, a3
    
    !------------------------------------------------------------------------!
    
    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename, tmesh%numtimes, tmesh%deltat, &
         rmesh%numthetapts, rmesh%numphipts, rmesh%lmax)
    
    ALLOCATE(tmesh%timepts(1:tmesh%numtimes))
    ALLOCATE(efield(1:3, tmesh%numtimes))
    ALLOCATE(afield(1:3, tmesh%numtimes))
        
    IF (.NOT. tmesh%truncate_time_integration) THEN
       
       tmesh%maxtime_integrate = tmesh%numtimes
       
    ELSE
       
       IF (tmesh%maxtime_integrate.GT.tmesh%numtimes) THEN
          
          tmesh%maxtime_integrate = tmesh%numtimes
          
       ENDIF
       
    ENDIF
    
    IF (commsurff%iprocessor.EQ.0) THEN
       
       WRITE(*, *)					                       &
            '--t-SURFF Surface File Information------------------------------'
       
       WRITE(*, *) 
       
       WRITE(*,'(A46, I9)')                                                  &
	    ' Number of time steps in the file:           ', tmesh%numtimes
       WRITE(*,'(A46, F9.3)')                                                &
	    ' Deltat:                                     ', tmesh%deltat 
       WRITE(*,'(A46, I9)')                                                  &
	    ' Number of theta points:                     ', rmesh%numthetapts
       WRITE(*,'(A46, I9)')                                                  &
	    ' Number of phi points:                       ', rmesh%numphipts
       WRITE(*, *) 
       
       WRITE(*, *)					                       &
            '----------------------------------------------------------------'
       
       WRITE(*, *) 
       
    ENDIF
    
    ! Create momentum axis
    
    kmesh%emax_theoretical = 0.5_dp * twopi * twopi / (tmesh%deltat *        &
         tmesh%deltat)
    
    IF (kmesh%emax.GT.kmesh%emax_theoretical) THEN
       
       CALL tsurff_stop( 'Maximum momentum desired too large.' )
       
    ENDIF
    
    IF (commsurff%iprocessor.EQ.0) THEN
       
       WRITE(*, *)					                       &
            '--t-SURFF Miscellaneous Information-----------------------------'
       
       WRITE(*, *)
       
       WRITE(*, '(A46, F9.3)')                                               &
	    ' Theoretical emax:                           ',                   &
            kmesh%emax_theoretical
       
       WRITE(*, '(A46, F9.3)')                                               &
	    ' Time in the file (a.u.):                    ',                   &
            tmesh%numtimes * tmesh%deltat
       WRITE(*, '(A46, F9.3)')                                               &
	    ' Time in the file (fs):                      ',                   &
            tmesh%numtimes * tmesh%deltat * autime_fs
       
       WRITE(*, '(A46, F9.3)')                                               &
	    ' Time to propagate (a.u.):                   ',                   &
            tmesh%maxtime_integrate * tmesh%deltat
       WRITE(*, '(A46, F9.3)')                                               &
	    ' Time to propagate (fs):                     ',                   &
            tmesh%maxtime_integrate * tmesh%deltat * autime_fs
       
       WRITE(*, *) 
       
       WRITE(*, *)					                       &
            '----------------------------------------------------------------'
       
    ENDIF
    
    ! Allocate wavefunction arrays
    
    ALLOCATE(psi_sph(1:rmesh%numthetapts, 1:rmesh%numphipts))
    ALLOCATE(psip_sph(1:rmesh%numthetapts, 1:rmesh%numphipts))

    ALLOCATE(psi_lm(-kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    ALLOCATE(psip_lm(-kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    
    ! Create spherical coordinates
    
    ALLOCATE(rmesh%thetapts(1:rmesh%numthetapts))
    ALLOCATE(rmesh%sinthetapts(1:rmesh%numthetapts))
    ALLOCATE(rmesh%costhetapts(1:rmesh%numthetapts))
    ALLOCATE(rmesh%wtheta(1:rmesh%numthetapts))
    ALLOCATE(rmesh%phipts(1:rmesh%numphipts))
    ALLOCATE(rmesh%sinphipts(1:rmesh%numphipts))
    ALLOCATE(rmesh%cosphipts(1:rmesh%numphipts))
    ALLOCATE(rmesh%wphi(1:rmesh%numphipts))
    
    CALL make_gauss_lobatto( -1.0_dp, 1.0_dp, rmesh%costhetapts, rmesh%wtheta )
    
    rmesh%thetapts    = ACOS(rmesh%costhetapts)
    rmesh%sinthetapts = SIN(rmesh%thetapts)
    rmesh%deltaphi    = twopi / REAL(rmesh%numphipts, dp)
    
    DO iphi = 1, rmesh%numphipts
       
       rmesh%phipts(iphi)    = REAL(iphi - 1, dp) * rmesh%deltaphi 
       rmesh%sinphipts(iphi) = SIN(rmesh%phipts(iphi))
       rmesh%cosphipts(iphi) = COS(rmesh%phipts(iphi))
       rmesh%wphi(iphi)      = rmesh%deltaphi
       
    ENDDO
    
    ! Create spherical harmonics
    ALLOCATE(sph_harmonics(1:rmesh%numphipts, 1:rmesh%numthetapts,          &
         -kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    ALLOCATE(sph_harmonicsk(1:kmesh%numphipts, 1:kmesh%numthetapts,          &
         -kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    
    CALL create_spherical_harmonics( sph_harmonics, kmesh%lmax, &
         rmesh%costhetapts, rmesh%phipts )
    CALL create_spherical_harmonics( sph_harmonicsk, kmesh%lmax, &
         kmesh%costhetapts, kmesh%phipts )
    
    ALLOCATE(xcoupling(-kmesh%lmax:kmesh%lmax,0:kmesh%lmax,&
         -kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    ALLOCATE(ycoupling(-kmesh%lmax:kmesh%lmax,0:kmesh%lmax,&
         -kmesh%lmax:kmesh%lmax,0:kmesh%lmax))
    ALLOCATE(zcoupling(-kmesh%lmax:kmesh%lmax,0:kmesh%lmax,&
         -kmesh%lmax:kmesh%lmax,0:kmesh%lmax))

    CALL create_spherical_harmonics_couplings( xcoupling, ycoupling, zcoupling, &
         sph_harmonics, kmesh%lmax, rmesh%thetapts, rmesh%phipts, &
         rmesh%wtheta, rmesh%wphi )
        
    ! Initialize Bessel functions
    
    ALLOCATE(jl(0:kmesh%lmax, 1:kmesh%maxkpts))
    ALLOCATE(jlp(0:kmesh%lmax, 1:kmesh%maxkpts))
    ALLOCATE(yl(0:kmesh%lmax))
    ALLOCATE(ylp(0:kmesh%lmax))
    ALLOCATE(krspts(1:kmesh%maxkpts))
    
    krspts = kmesh%kpts * rmesh%rsurface
    
    jl  = 0.0_dp
    jlp = 0.0_dp
    yl  = 0.0_dp
    ylp = 0.0_dp
      
    DO ik = 1, kmesh%maxkpts
       
       CALL sbesjy( krspts(ik), kmesh%lmax, jl(0:kmesh%lmax, ik), yl,        &
            jlp(0:kmesh%lmax, ik), ylp, ifail )
       
    ENDDO
    
    ! Allocate momentum amplitude array
    
    ALLOCATE(bktime(1:kmesh%maxkpts, 1:kmesh%maxthetapts,                    &
         1:kmesh%maxphipts))
    ALLOCATE(bktimeold(1:kmesh%maxkpts, 1:kmesh%maxthetapts,                 &
         1:kmesh%maxphipts))
    ALLOCATE(bk(1:kmesh%maxkpts, 1:kmesh%maxthetapts, 1:kmesh%maxphipts,     &
         1:tmesh%numtimes))
    ALLOCATE(bkaverage(1:kmesh%maxkpts, 1:kmesh%maxthetapts,                 &
         1:kmesh%maxphipts))
    ALLOCATE(probk3D(1:kmesh%maxkpts, 1:kmesh%maxthetapts,                   &
         1:kmesh%maxphipts))
    
    ! Allocate and set to zero
    
    ALLOCATE(quiver(1:3, tmesh%numtimes))
    ALLOCATE(asqintegral(tmesh%numtimes))

    CALL read_field_surface( filename, tmesh%timepts, efield, afield )

    CALL make_field_integral_terms(  )
        
  END SUBROUTINE initialize_tsurff
  
  !-----------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------!
  
  SUBROUTINE make_kpoints( lmax, maxkpts, maxthetapts, maxphipts, deltaenergy, &
       eneoffset )
    
    IMPLICIT NONE

    !------------------------------------------------------------------------!

    INTEGER, INTENT(IN)     :: lmax
    INTEGER, INTENT(IN)     :: maxkpts
    INTEGER, INTENT(IN)     :: maxthetapts
    INTEGER, INTENT(IN)     :: maxphipts
    REAL(dp), INTENT(IN)    :: deltaenergy
    REAL(dp), INTENT(IN), OPTIONAL :: eneoffset
    
    !------------------------------------------------------------------------!
    
    INTEGER                 :: ik, itheta, iphi
    INTEGER                 :: ithetag
    REAL(dp), ALLOCATABLE   :: costheta_global(:)
    REAL(dp), ALLOCATABLE   :: wtheta_global(:)
    
    !------------------------------------------------------------------------!

    kmesh%lmax         = lmax
    kmesh%maxkpts      = maxkpts
    kmesh%maxthetapts  = maxthetapts
    kmesh%maxphipts    = maxphipts
    kmesh%deltaenergy  = deltaenergy
    IF(PRESENT(eneoffset)) THEN
       kmesh%eneoffset    = eneoffset
    ELSE
       kmesh%eneoffset    = 0.0_dp
    ENDIF
    
    kmesh%numkpts      = commsurff%numproc1dk * kmesh%maxkpts
    kmesh%numthetapts  = commsurff%numproc1dt * kmesh%maxthetapts
    kmesh%numphipts    = commsurff%numproc1dp * kmesh%maxphipts
      
    ALLOCATE(kmesh%kpts(1:kmesh%maxkpts))
    ALLOCATE(kmesh%energypts(1:kmesh%maxkpts))
    ALLOCATE(kmesh%thetapts(1:kmesh%maxthetapts))
    ALLOCATE(kmesh%sinthetapts(1:kmesh%maxthetapts))
    ALLOCATE(kmesh%costhetapts(1:kmesh%maxthetapts))
    ALLOCATE(kmesh%wtheta(1:kmesh%maxthetapts))
    ALLOCATE(kmesh%phipts(1:kmesh%maxphipts))
    ALLOCATE(kmesh%sinphipts(1:kmesh%maxphipts))
    ALLOCATE(kmesh%cosphipts(1:kmesh%maxphipts))
    ALLOCATE(kmesh%wphi(1:kmesh%maxphipts))
    
    DO ik = 1, kmesh%maxkpts 
       
       kmesh%energypts(ik) = kmesh%eneoffset +             &
            REAL(kmesh%maxkpts * commsurff%ipk + ik, dp) * &
            kmesh%deltaenergy
       kmesh%kpts(ik)      = SQRT(2.0_dp * kmesh%energypts(ik))
       
    ENDDO
    
    kmesh%emin         = kmesh%energypts(1)
    kmesh%emax         = kmesh%energypts(kmesh%maxkpts)
    
    ! Prepare Gauss nodes and weights with a Gauss-Lobatto subroutine
    
    ALLOCATE(costheta_global(1:kmesh%numthetapts))
    ALLOCATE(wtheta_global(1:kmesh%numthetapts))
    
    CALL make_gauss_lobatto( -1.0_dp, 1.0_dp, costheta_global, wtheta_global )
    
    DO itheta = 1, kmesh%maxthetapts
       
       ithetag                   = kmesh%maxthetapts * commsurff%ipt +       &
            itheta
       
       kmesh%costhetapts(itheta) = costheta_global(ithetag)
       
       kmesh%wtheta(itheta)      = wtheta_global(ithetag)
       
       kmesh%thetapts(itheta)    = ACOS(kmesh%costhetapts(itheta))
       
       kmesh%sinthetapts(itheta) = SIN(kmesh%thetapts(itheta))
       
    ENDDO
    
    DEALLOCATE(costheta_global)
    DEALLOCATE(wtheta_global)
    
    kmesh%deltaphi = twopi / REAL(kmesh%numphipts, dp)
    
    DO iphi = 1, kmesh%maxphipts 
       
       kmesh%phipts(iphi) =                                                  &
            REAL(kmesh%maxphipts * commsurff%ipp + iphi - 1, dp) *                &
            kmesh%deltaphi
       
       kmesh%sinphipts(iphi) = SIN(kmesh%phipts(iphi))
       
       kmesh%cosphipts(iphi) = COS(kmesh%phipts(iphi))

       kmesh%wphi(iphi)      = kmesh%deltaphi
       
    ENDDO
    
  END SUBROUTINE make_kpoints
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
  SUBROUTINE get_flux( filename, iorbital )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    CHARACTER(LEN = *), INTENT(IN) :: filename
    INTEGER, INTENT(IN)	           :: iorbital
    
    !------------------------------------------------------------------------!
    
    INTEGER	                     :: ik, itheta, iphi, itime
    INTEGER	                     :: itimeaveragemin
    REAL(dp)	                     :: tfac, tintfac
    COMPLEX(dp)	                     :: fluxintegral
    COMPLEX(dp)	                     :: volkov_phase
    
    !------------------------------------------------------------------------!
    
    bk = ZERO
    
    itimeaveragemin = tmesh%maxtime_integrate -                              &
         INT(tmesh%timetoaverageover / tmesh%deltat)
    
    IF (commsurff%iprocessor.EQ.0) THEN
       
       WRITE(*, *)
       
       WRITE(*, *)					                       &
            '--Begin time integral-------------------------------------------'
       
       WRITE(*, *)
       
       WRITE(*, '(A46, I9)')                                                 &
	    ' Orbital number:                             ', iorbital
       WRITE(*, *)
       
    ENDIF

    itime = 1
    
    CALL read_wave_surface( filename, itime - 1 , psi_sph, psip_sph )
    
    IF (commsurff%iprocessor.EQ.0) THEN
       
       WRITE(*, '(1X, A8, 1X, I6, 1X, A2, 1X, I6)')                          &
            'Timestep', itime, 'of', tmesh%maxtime_integrate
       
    ENDIF
    
    IF (rmesh%gauge_transform) THEN
       
       CALL gauge_transform_l2v( itime )
       
    ENDIF
    
    ! Decompose into spherical harmonics:
    
    CALL make_sht( psi_sph, sph_harmonics, kmesh%lmax, &
         rmesh%wtheta, rmesh%wphi, psi_lm)
    CALL make_sht( psip_sph, sph_harmonics, kmesh%lmax, &
         rmesh%wtheta, rmesh%wphi, psip_lm)
        
    ! Calculate flux
    
    DO iphi = 1, kmesh%maxphipts
       DO itheta = 1, kmesh%maxthetapts
          DO ik = 1, kmesh%maxkpts
             
             CALL calculate_time_integrand( fluxintegral, ik, itheta,        &
                  iphi, itime )
             
             CALL get_volkov_phase( volkov_phase, ik, itheta, iphi,	       &
                  itime )
             
             bktimeold(ik, itheta, iphi) = fluxintegral * volkov_phase
             
          ENDDO
       ENDDO
    ENDDO
    
    DO itime = 2, tmesh%maxtime_integrate
       
       IF (commsurff%iprocessor.EQ.0) THEN
          
          WRITE(*, '(1X, A8, 1X, I6, 1X, A2, 1X, I6)')                       &
	       'Timestep', itime, 'of', tmesh%maxtime_integrate
          
       ENDIF
       
       ! Read wavefunction at Rb from file
       CALL read_wave_surface( filename, itime - 1, psi_sph, psip_sph )
       
       IF (rmesh%gauge_transform) THEN
          
          CALL gauge_transform_l2v( itime )
          
       ENDIF
       
       ! Decompose into spherical harmonics:
       
       CALL make_sht( psi_sph, sph_harmonics, &
            kmesh%lmax, rmesh%wtheta, rmesh%wphi, psi_lm)
       CALL make_sht( psip_sph, sph_harmonics, &
            kmesh%lmax, rmesh%wtheta, rmesh%wphi, psip_lm)
       
       ! Calculate flux
       
       DO iphi = 1, kmesh%maxphipts
          DO itheta = 1, kmesh%maxthetapts
             DO ik = 1, kmesh%maxkpts
                
                CALL calculate_time_integrand( fluxintegral, ik, itheta,     &
                     iphi, itime )
                
                CALL get_volkov_phase( volkov_phase, ik, itheta, iphi,       &
                     itime )

                bktime(ik, itheta, iphi)    = fluxintegral * volkov_phase
                
                bk(ik, itheta, iphi, itime) =                                &
                     bk(ik, itheta, iphi, itime - 1) +                       &
                     bktime(ik, itheta, iphi) + bktimeold(ik, itheta, iphi)
                
                bktimeold(ik, itheta, iphi) = bktime(ik, itheta, iphi)
                
             ENDDO
          ENDDO
       ENDDO
       
    ENDDO
    
    IF (tmesh%timeaverageforrydberg) THEN
       
       DO iphi = 1, kmesh%maxphipts
          DO itheta = 1, kmesh%maxthetapts
             DO ik = 1, kmesh%maxkpts
                
                bkaverage(ik, itheta, iphi) = bk(ik, itheta, iphi,           &
                     itimeaveragemin)
                
             ENDDO
          ENDDO
       ENDDO
       
       DO itime = itimeaveragemin + 1, tmesh%maxtime_integrate
          
          DO iphi = 1, kmesh%maxphipts
             DO itheta = 1, kmesh%maxthetapts
                DO ik = 1, kmesh%maxkpts
                   
                   bkaverage(ik, itheta, iphi) =                             &
                        bkaverage(ik, itheta, iphi) +                             &
                        bk(ik, itheta, iphi, itime - 1) +                         &
                        bk(ik, itheta, iphi, itime)
		  
                ENDDO
             ENDDO
          ENDDO
          
       ENDDO
       
       tintfac = 0.25_dp * tmesh%deltat * tmesh%deltat /                     &
            (tmesh%timepts(tmesh%maxtime_integrate) -                   &
            tmesh%timepts(itimeaveragemin))
       
    ELSE   
       
       tintfac   = 0.5_dp * tmesh%deltat
       
       DO iphi = 1, kmesh%maxphipts
          DO itheta = 1, kmesh%maxthetapts
             DO ik = 1, kmesh%maxkpts
                
                bkaverage(ik, itheta, iphi) =  bk(ik, itheta, iphi,          &
                     tmesh%maxtime_integrate)
                
             ENDDO
          ENDDO
       ENDDO
       
       
    ENDIF
    
    IF (commsurff%iprocessor.EQ.0) THEN
       
       WRITE(*, *)
       
       WRITE(*, *)					                       &
            '--End time integral---------------------------------------------'
       
       WRITE(*, *)
       
    ENDIF
    
    bkaverage = bkaverage * ZIMAGONE * SQRT(2.0_dp / pi) * rmesh%rsurface *     &
         rmesh%rsurface * tintfac
    
    
  END SUBROUTINE get_flux
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
  SUBROUTINE calculate_time_integrand( fluxintegral, ik, itheta, iphi,     &
       itime )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    INTEGER, INTENT(IN)               :: ik, itheta, iphi, itime
    COMPLEX(dp), INTENT(OUT)	      :: fluxintegral
    
    !------------------------------------------------------------------------!
    
    INTEGER                           :: il, im, ill, imm
    COMPLEX(dp)	                      :: term1, term2
    COMPLEX(dp)                       :: suma, fieldterm
    
    !------------------------------------------------------------------------!
    
    IF ( rmesh%numphipts .EQ. 1 ) THEN
       
       suma = ZERO
       
       DO il = 0, kmesh%lmax
          
          term1 = psi_lm(0, il) *			               &
               sph_harmonicsk(iphi, itheta, 0, il)
          term2 = psip_lm(0, il) *			               &
               sph_harmonicsk(iphi, itheta, 0, il)
          
          fieldterm = ZERO
          
          DO ill = 0, kmesh%lmax
             
             IF (ABS(il - ill).EQ.1) THEN
                
                fieldterm = fieldterm +			               &
                     (xcoupling(0, il, 0, ill) *                       &
                     afield(1, itime) +  		               &
                     ycoupling(0, il, 0, ill) *                        &
                     ZIMAGONE * afield(2, itime)) *	               &
                     psi_lm(0, ill)
                
             ENDIF
             
             IF (ABS(il - ill).EQ.1) THEN
                
                fieldterm = fieldterm +			               &
                     zcoupling(0, il, 0, ill) *                        &
                     afield(3, itime) * psi_lm(0, ill)
                
             ENDIF
             
          ENDDO
          
          fieldterm = ZIMAGONE * fieldterm *  		               &
               sph_harmonicsk(iphi, itheta, 0, il)
          
          
          suma = suma + ((-ZIMAGONE) ** il) *			       &
               (0.5_dp * kmesh%kpts(ik) * jlp(il, ik) * term1 -	       &
               0.5_dp * jl(il, ik) * term2 - jl(il, ik) * fieldterm)
          
       ENDDO

       fluxintegral = suma
       
    ELSE
       
       suma = ZERO
       
       DO il = 0, kmesh%lmax
          
          term1 = ZERO
          term2 = ZERO
          
          DO im = -il, il
             
             term1 = term1 + psi_lm(im, il) *	                      &
                  sph_harmonicsk(iphi, itheta, im, il)
             term2 = term2 + psip_lm(im, il) *			      &
                  sph_harmonicsk(iphi, itheta, im, il)
             
             fieldterm = ZERO
             
             DO ill = 0, kmesh%lmax
                DO imm = -ill, ill
                   
                   IF (ABS(il - ill).EQ.1 .AND.  		      &
                        ABS(im - imm).EQ.1) THEN
                      
                      fieldterm = fieldterm +			      &
                           (xcoupling(im, il, imm, ill) *             &
                           afield(1, itime) +  		              &
                           ycoupling(im, il, imm, ill) *              &
                           ZIMAGONE * afield(2, itime)) *	      &
                           psi_lm(imm, ill)
                      
                   ENDIF
                   
                   IF (ABS(il - ill).EQ.1 .AND.  		      &
                        ABS(im - imm).EQ.0) THEN
                      
                      fieldterm = fieldterm +			      &
                           zcoupling(im, il, imm, ill) *              &
                           afield(3, itime) * psi_lm(imm, ill)
                      
                   ENDIF
                   
                ENDDO
             ENDDO
             
             fieldterm = ZIMAGONE * fieldterm *  		      &
                  sph_harmonicsk(iphi, itheta, im, il)
             
          ENDDO
          
          suma = suma + ((-ZIMAGONE) ** il) *			      &
               (0.5_dp * kmesh%kpts(ik) * jlp(il, ik) * term1 -	      &
               0.5_dp * jl(il, ik) * term2 - jl(il, ik) * fieldterm)
          
       ENDDO
       
       fluxintegral = suma
       
    ENDIF
    
  END SUBROUTINE calculate_time_integrand
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
  SUBROUTINE get_volkov_phase( volkov_phase, ik, itheta, iphi, itime )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    INTEGER, INTENT(IN)               :: ik, itheta, iphi, itime
    COMPLEX(dp), INTENT(OUT)          :: volkov_phase
    
    !------------------------------------------------------------------------!
    
    REAL(dp)                          :: kdota, kene
    
    !------------------------------------------------------------------------!
    
    
    kdota = kmesh%kpts(ik) * kmesh%sinthetapts(itheta) *                     &
         kmesh%cosphipts(iphi) * quiver(1, itime) +		       &
         kmesh%kpts(ik) * kmesh%sinthetapts(itheta) *                     &
         kmesh%sinphipts(iphi) * quiver(2, itime) +		       &
         kmesh%kpts(ik) * kmesh%costhetapts(itheta) *                     &
         quiver(3, itime)
    kene  = kmesh%energypts(ik) 
    
    
    volkov_phase = EXP(ZIMAGONE * (kene * tmesh%timepts(itime) + kdota))
    
    
  END SUBROUTINE get_volkov_phase
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
  SUBROUTINE make_field_integral_terms( )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    INTEGER                   :: icoord, itime
    REAL(dp)                  :: asq, asqold
    
    !------------------------------------------------------------------------!
    
    quiver      = 0.0_dp
    asqintegral = 0.0_dp
    
    
    asqold = afield(1, 1) * afield(1, 1) + afield(2, 1) * afield(2, 1) +     &
         afield(3, 1) * afield(3, 1)
    
    DO itime = 2, tmesh%numtimes
       DO icoord = 1, 3
          
          quiver(icoord, itime) = quiver(icoord, itime - 1) +                &
               afield(icoord, itime - 1) +                &
               afield(icoord, itime)
          
       ENDDO
       
       asq = afield(1, itime) * afield(1, itime) +                           &
            afield(2, itime) * afield(2, itime) +                           &
            afield(3, itime) * afield(3, itime)
       
       asqintegral(itime) = asqintegral(itime - 1) + asqold + asq
       
       asqold = asq
       
    ENDDO
    
    quiver      = 0.5_dp * quiver * tmesh%deltat
    asqintegral = 0.5_dp * asqintegral * tmesh%deltat
    
    
  END SUBROUTINE make_field_integral_terms
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
  SUBROUTINE gauge_transform_l2v( itime )
    
    IMPLICIT NONE
    
    !------------------------------------------------------------------------!
    
    INTEGER, INTENT(IN)       :: itime
    
    INTEGER                   :: itheta, iphi
    REAL(dp)                  :: adotr, adotrhat, gterm
    
    !------------------------------------------------------------------------!
    
    DO iphi = 1, rmesh%numphipts
       DO itheta = 1, rmesh%numthetapts
          
          adotrhat = afield(1, itime) * rmesh%sinthetapts(itheta) *       &
               rmesh%cosphipts(iphi) +                                    &
               afield(2, itime) * rmesh%sinthetapts(itheta) *             &
               rmesh%sinphipts(iphi) +                                    &
               afield(3, itime) * rmesh%costhetapts(itheta) 
          
          adotr = adotrhat * rmesh%rsurface
          
          gterm = EXP(ZIMAGONE * (0.5_dp * asqintegral(itime) - adotr))
          
          psip_sph(itheta, iphi) = gterm * (psip_sph(itheta, iphi) -      &
               ZIMAGONE * adotrhat *                                      &
               psi_sph(itheta, iphi))
          psi_sph(itheta, iphi)  = psi_sph(itheta, iphi) * gterm
          
       ENDDO
    ENDDO
    
  END SUBROUTINE gauge_transform_l2v
  
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  !------------------------------------------------------------------------!
  
END MODULE flux
