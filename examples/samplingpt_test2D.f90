!-----------------------------------------------------!
!
!  PROGRAM Sampling point test 2D
!  \brief This example test the detector method.
!  \details This program runs and example for
!  the detector method. In cylindrical coordinates,
!  It tracks the evolution of a gaussian wavepacket.
!  It produces a file with the all the recorder points
!  at the surface, and also the computed spectra.
!
!-----------------------------------------------------!
PROGRAM samplingpt_test2D

  USE popsicle
  
  IMPLICIT NONE
  
  INTEGER                     :: numxpts, numzpts
  INTEGER                     :: dims(2), halfxpts
  INTEGER                     :: rank, lmax
  REAL(dp)                    :: dx, dz, dt
  REAL(dp)                    :: dkx, dkz
  COMPLEX(dp), ALLOCATABLE    :: psi(:, :), psik(:, :)
  REAL(dp), ALLOCATABLE       :: x_ax(:) , z_ax(:)
  REAL(dp), ALLOCATABLE       :: kx_ax(:) , kz_ax(:)
  INTEGER                     :: numdetectorpts
  REAL(dp)                    :: timept, numtimesteps
  INTEGER                     :: timesteps
  INTEGER                     :: no_detectors
  REAL(dp)                    :: rad_detector
  INTEGER                     :: numwpts
  REAL(dp), ALLOCATABLE       :: pes(:), pad(:, :, :)
  REAL(dp)                    :: sumatot, sumaion, rpt
  
  ! Waveapcket parameters
  REAL(dp)                    :: px0, pz0,chi0, mu0
  REAL(dp)                    :: r, r0, s, s0
  REAL(dp)                    :: omega, omega0

  ! Pulse parameters
  REAL(dp)                    :: intensity, intensity0
  REAL(dp)                    :: wavelength, wavelength0
  REAL(dp)                    :: period0, durationpulse0
  INTEGER                     :: nocycles, afterpulse
  REAL(dp)                    :: interactiontime
  REAL(dp)                    :: efield(3), afield(3)
  REAL(dp)                    :: e0, a0, wl
  REAL(dp)                    :: ft, gt, intgt, intgt2
  
  LOGICAL                     :: wave_snapshot
  INTEGER                     :: ix, iz, itime, ipt
  
  !------------------------------------------------!

  WRITE(*,*)
  WRITE(*,*) '------------------------------------'
  WRITE(*,*)
  WRITE(*,*) '    Sampling point method test      '
  WRITE(*,*)
  WRITE(*,*) '------------------------------------'

  WRITE(*,*)
  WRITE(*,*) 'Create arrays and coordiante axis...'
  WRITE(*,*) '------------------------------------'
    
  numxpts   = 401
  numzpts   = 401
  dims      = (/numxpts, numzpts/)
  rank      = 2
  dx        = 0.4_dp
  dz        = 0.4_dp
  dkx       = twopi / ( numxpts * dx )
  dkz       = twopi / ( numzpts * dz )
  halfxpts   = (numxpts - 1) * 0.5
  
  ! Allocate arrays
  ALLOCATE(x_ax(numxpts),z_ax(numzpts))
  ALLOCATE(psi(numxpts,numzpts))
  ALLOCATE(kx_ax(numxpts),kz_ax(numzpts))
  ALLOCATE(psik(numxpts,numzpts))
  
  DO ix = 1, numxpts
     x_ax(ix) = (-REAL(numxpts-1,dp) * 0.5_dp + &
          REAL(ix-1,dp) ) * dx
     kx_ax(ix) = (-REAL(numxpts-1,dp) * 0.5_dp + &
          REAL(ix-1,dp) ) * dkx
  ENDDO
  
  DO iz = 1, numzpts
     z_ax(iz) = (-REAL(numzpts-1,dp) * 0.5_dp + &
          REAL(iz-1,dp) ) * dz 
     kz_ax(iz) = (-REAL(numzpts-1,dp) * 0.5_dp + &
          REAL(iz-1,dp) ) * dkz 
  ENDDO
  
  ! Set the detector
  WRITE(*,*)
  WRITE(*,*) 'Swicth on the detector '
  WRITE(*,*) '(Initialize detector points)...'
  WRITE(*,*) '------------------------------------'
  
  numdetectorpts = 100
  rad_detector = 30.0_dp
  lmax = numdetectorpts - 1
  
  CALL initialize_cylindrical_surface(x_ax(halfxpts:numxpts), &
       z_ax(1:numzpts), dims, rad_detector, &
       radius_tolerance=0.5_dp, fd_rule=2, dr=0.1_dp, &
       lmax=lmax, filename='./results/detector_data')
  
  ! Open file for ionised population
  wave_snapshot = .TRUE.
  OPEN(UNIT=33,FORM='formatted',FILE='./results/population.dat')
  OPEN(UNIT=35,FORM='formatted',FILE='./results/field.dat')
  
  ! Initial parameters for the wavepacket
  px0 = 1.0_dp
  pz0 = 1.0_dp
  r0 = 0.0_dp
  s0 = 0.0_dp
  mu0 = 3.0_dp
  chi0 = 3.0_dp
  omega0 = 0.0_dp

  WRITE(*,*) '------------------------------------'
  WRITE(*,*) 'Start evolution...'
  WRITE(*,*) '------------------------------------'

  DO iz = 1, numzpts
     DO ix = 1, numxpts
        psi(ix,iz) = ( mu0 / (pi * mu0**2))**(0.25_dp) * &
             ( chi0 / (pi * chi0**2))**(0.25_dp) * &
             EXP(ZIMAGONE * ( px0 * x_ax(ix) + pz0 * z_ax(iz) )) * &
             EXP( - (x_ax(ix) - r0)**2 / (2.0_dp * mu0) ) * &
             EXP( - (z_ax(iz) - s0)**2 / (2.0_dp * chi0) ) * &
             EXP( - ZIMAGONE * omega0 )

        psik(ix,iz) = (mu0 / pi)**(0.25_dp) * &
             (chi0 / pi)**(0.25_dp) * &
             EXP(-ZIMAGONE * r0 * (px0 - kx_ax(ix))) * &
             EXP(-ZIMAGONE * s0 * (pz0 - kz_ax(iz))) * &
             EXP(- 0.5_dp * mu0 * (kx_ax(ix) - px0)**2) * &
             EXP(- 0.5_dp * chi0 * (kx_ax(iz) - pz0)**2)
     ENDDO
  ENDDO

  timept = 0.0_dp
  efield = 0.0_dp
  afield = 0.0_dp
  sumatot = 0.0_dp
  sumaion = 0.0_dp
  DO ix = 1, numxpts
     DO iz = 1, numzpts
        sumatot = sumatot + &
             REAL(CONJG(psi(ix,iz)) * psi(ix,iz))
        rpt = SQRT(x_ax(ix)**2 + z_ax(iz)**2)
        IF(rpt.GT.rad_detector) THEN
           sumaion = sumaion + &
                REAL(CONJG(psi(ix,iz)) * psi(ix,iz))
           
        ENDIF
     ENDDO
  ENDDO
  
  sumatot = sumatot * dz * dx
  sumaion = sumaion * dz * dx
  WRITE(*,*) timept, sumatot, sumaion
  WRITE(33,*) timept, sumatot, sumaion
  WRITE(35,*) timept, efield(3), afield(3)
  
  IF(wave_snapshot) THEN
     OPEN(unit=30,form='formatted',file='./results/initial_psi.dat')
     DO ix = 1, numxpts
        DO iz = 1, numzpts
           write(30,*) ABS(psi(ix,iz))**2
        ENDDO
     ENDDO
     close(30)


     OPEN(unit=30,form='formatted',file='./results/initial_psik.dat')
     DO ix = 1, numxpts
        DO iz = 1, numzpts
           write(30,*) ABS(psik(ix,iz))**2
        ENDDO
     ENDDO
     close(30)
  ENDIF
  
  !------------------- Start the evolution. -------------------!
  !---------------- The detector is recording! ----------------!

  ! Time step
  dt = 0.2_dp
  
  ! Set electric field
  intensity = 2.0E14_dp
  wavelength = 23.0_dp
  nocycles   = 5
  afterpulse = 30
  
  wavelength0 = wavelength / aulength_nm
  wl = twopi * speed_light / wavelength0
  period0 = twopi / wl
  e0 = SQRT(intensity / intensity_au)
  a0 = e0 / wl
  durationpulse0  = nocycles * period0
  interactiontime  = ( nocycles + afterpulse ) * period0
  numtimesteps = INT(interactiontime / dt)
  ft = 0.0_dp
  gt = 0.0_dp
  intgt = 0.0_dp
  intgt2 = 0.0_dp
      
  DO itime = 1, numtimesteps
     timept = REAL(itime,dp) * dt
     
     ! Envelope and derivatives
     IF(timept.LE.durationpulse0) THEN
        ft = SIN( pi * timept / durationpulse0 )**2 * &
             SIN(wl * timept)
     ELSE
        ft = 0.0_dp
     ENDIF
     gt = gt + ft * dt
     intgt = intgt + gt * dt
     intgt2 = intgt2 + gt**2 * dt
     
     ! Electric field and vector potential
     efield(3) = e0 * ft
     afield(3) = - a0 * gt
     
     ! Wavepackets phases
     r = r0 + px0 * timept
     s = s0 + pz0 * timept - e0 * intgt
     omega = omega0 + &
          0.5_dp * (px0**2 + pz0**2) * timept - &
          pz0 * e0 * intgt + &
          0.5_dp * e0**2 * intgt2 + &
          0.5_dp * ( ATAN2(timept, mu0) + ATAN2(timept,chi0) )
     
     DO ix = 1, numxpts
        DO iz = 1, numzpts
           psi(ix,iz) = ( mu0 / (pi * (mu0**2 + timept**2)))**(0.25_dp) * &
                ( chi0 / (pi * (chi0**2 + timept**2 ) ))**(0.25_dp) * &
                EXP(ZIMAGONE * ( px0 * x_ax(ix) + pz0 * z_ax(iz) )) * &
                EXP( - (x_ax(ix) - r)**2 / (2.0_dp * (mu0 + ZIMAGONE * timept)) ) * &
                EXP( - (z_ax(iz) - s)**2 / (2.0_dp * (chi0 + ZIMAGONE * timept)) ) * &
                EXP( - ZIMAGONE * omega)
        ENDDO
     ENDDO

     CALL get_cylindrical_surface('./results/detector_data', psi(halfxpts:numxpts,:), &
          x_ax(halfxpts:numxpts), z_ax(1:numzpts),dims, 2, timept, efield, afield, &
          lmax)
     
     sumatot = 0.0_dp
     sumaion = 0.0_dp
     DO ix = 1, numxpts
        DO iz = 1, numzpts
           sumatot = sumatot + &
                REAL(CONJG(psi(ix,iz)) * psi(ix,iz))
           
           rpt = SQRT(x_ax(ix)**2 + z_ax(iz)**2)
           IF(rpt.GT.rad_detector) THEN
              sumaion = sumaion + &
                   REAL(CONJG(psi(ix,iz)) * psi(ix,iz))
              
           ENDIF
        ENDDO
     ENDDO
     
     sumatot = sumatot * dz * dx
     sumaion = sumaion * dz * dx
     WRITE(*,*) timept, sumatot, sumaion
     WRITE(33,*) timept, sumatot, sumaion
     WRITE(35,*) timept, efield(3), afield(3)
     
  ENDDO

  WRITE(*,*) '------------------------------------'
  WRITE(*,*) 'End of evolution.'
  WRITE(*,*) '------------------------------------'
  IF(wave_snapshot) THEN
     OPEN(unit=30,form='formatted',file='./results/final_psi.dat')
     DO ix = 1, numxpts
        DO iz = 1, numzpts
           write(30,*) ABS(psi(ix,iz))**2
        ENDDO
     ENDDO
     CLOSE(30)
  ENDIF
  
  CLOSE(33)
  !-------------------------- Stop! --------------------------!
  

  !-------------------- Read back the data -------------------!
  
  timesteps    = 0
  no_detectors = 0
  
  CALL initialize_sampling_points('./results/detector_data', rad_detector, &
       no_detectors, numthetapts, numphipts, timeafter_fs=0.0_dp, &
       coulomb_exp=0.0_dp, maxwpts=numwpts )

  timesteps = 2 * (numwpts - 1) + 1
  IF(numphipts.EQ.1) THEN
     rank = 2
  ELSE
     rank=3
  ENDIF
  
  WRITE(*,*)
  WRITE(*,*) 'Reading data back...'
  WRITE(*,*) '------------------------------------'
  WRITE(*,*) 'Number of detector ponts: ', no_detectors
  WRITE(*,*) 'Number of timesteps in the simulations: ',timesteps
  WRITE(*,*) 'Dimension of the calculation: ', rank
  WRITE(*,*) '------------------------------------'

  ! Allocate arrays for reading back
  ALLOCATE(pes(numwpts))
  ALLOCATE(pad(numthetapts,numphipts,numwpts))
  
  !-----------------------------------------------------------!
  !---------- Now, from the data, compute the spectra --------!
  
  WRITE(*,*)
  WRITE(*,*) 'Calculating pes and pad...'
  WRITE(*,*) '------------------------------------'
  
  CALL calculate_sampling_pes('./results/detector_data', pes, pad)
  
  ! Save pes and pad
  WRITE(*,*)
  WRITE(*,*) 'Writing results to file...'
  WRITE(*,*) '------------------------------------'
  
  CALL write_sampling_pes('./results/pes',pes)
  
  CALL write_sampling_polar('./results/pad', pad,'pad')
  
  ! Release memory
  DEALLOCATE(x_ax,z_ax)
  DEALLOCATE(psi)
  DEALLOCATE(pes,pad)

  WRITE(*,*)
  WRITE(*,*) '------------------------------------'  
  WRITE(*,*) 'Fin' 
  WRITE(*,*) '------------------------------------'
  WRITE(*,*) '------------------------------------'

    
END PROGRAM samplingpt_test2D
