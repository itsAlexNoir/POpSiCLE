PROGRAM sampling_calculator
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: numwpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts 
  REAL(dp)                   :: coulomb_exp_energy
  REAL(dp), ALLOCATABLE      :: pes(:), pad(:, :, :)
  REAL(dp)                   :: aftertime_fs
  CHARACTER(LEN=100)         :: filename_sampling
  CHARACTER(LEN=100)         :: filename_wavevstime
  CHARACTER(LEN=100)         :: filename_pes
  CHARACTER(LEN=100)         :: filename_polar
  CHARACTER(LEN=100)         :: filename_amplitude
  CHARACTER(LEN=1)           :: answer

  LOGICAL                    :: desired_wavevstime
  LOGICAL                    :: desired_polar
  LOGICAL                    :: desired_amplitude
  
  !--------------------------------------!

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '*******************************'
  WRITE(*,*) '   Sampling calculator app     '
  WRITE(*,*) '*******************************'
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Welcome to the sampling calculator. This little program'
  WRITE(*,*) ' will help to calculate photoelectron spectra from '
  WRITE(*,*) ' a surface file (HDF5), with the sampling point method.'
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Initializing...'
  WRITE(*,*)

  
  WRITE(*, *) 'Path of the surface file:'
  WRITE(*, *) '-------------------------'
  READ(*, '(1A80)') filename_sampling
  WRITE(*, *)
  
  WRITE(*, *) 'Boundary radius (a.u.):'
  WRITE(*, *) '-----------------------'
  READ(*, *) radius_boundary
  WRITE(*, *)

  WRITE(*, *) 'Do you want to add time at the end of the evolution? (y or n)'
  WRITE(*, *) '---------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  aftertime_fs = 0.0_dp
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     WRITE(*, *) 'Time added after the evolution (fs):'
     WRITE(*, *) '--------------------------------'
     READ(*, *) aftertime_fs
     WRITE(*, *)
     
  ENDIF
  
  WRITE(*, *) 'Do you want to add the Coulomb explosion energy? (y or n)'
  WRITE(*, *) '---------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  coulomb_exp_energy = 0.0_dp
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     WRITE(*, *) 'Coulomb explosion energy (a.u.):'
     WRITE(*, *) '--------------------------------'
     READ(*, *) coulomb_exp_energy
     WRITE(*, *)
     
     
  ENDIF

  WRITE(*, *) 'Do you want to plot the probabillity density vs time.'
  WRITE(*, *) '-----------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN

     desired_wavevstime = .TRUE.
     
     WRITE(*, *) 'Name of the wave vs time file.  '
     WRITE(*, *) '--------------------------------'
     READ(*, *) filename_wavevstime
     WRITE(*, *)
     
     
  ENDIF
  
  WRITE(*, *) 'Name of the PES file:'
  WRITE(*, *) '---------------------'
  READ(*,'(1A80)') filename_pes
  WRITE(*, *)
  

  WRITE(*, *) 'Do you want the 3D spherical scattering amplitude? (y or n)'
  WRITE(*, *) '-----------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)

  desired_amplitude = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     desired_amplitude = .TRUE.

     WRITE(*, *) 'Name of the spherical amplitude file:'
     WRITE(*, *) '-------------------------------------'
     READ(*,'(1A80)') filename_amplitude
     WRITE(*, *)
   
  ENDIF
  
  WRITE(*, *) 'Do you want polar scattering amplitude? (y or n)'
  WRITE(*, *) '------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)

  desired_polar = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     desired_polar = .TRUE.

     WRITE(*, *) 'Name of the polar amplitude file:'
     WRITE(*, *) '---------------------------------'
     READ(*,'(1A80)') filename_polar
     WRITE(*, *)
   
  ENDIF
  

  CALL initialize_sampling_points(filename_sampling, radius_boundary, &
       numdetectorpts, numthetapts, numphipts, aftertime_fs, &
       coulomb_exp_energy, numwpts )

  
  ! Allocate momentum amplitude array
  ALLOCATE(pes(1:numwpts))
  ALLOCATE(pad(1:numthetapts,1:numphipts,1:numwpts))

  ! The name is self-explanatory
  CALL calculate_sampling_pes(filename_sampling, pes, pad)
  
  !-------------------!
  !    OUTPUT TIME!   !
  !-------------------!

  IF(desired_wavevstime) &
       CALL get_wave_time_file(filename_wavevstime)
  
  CALL write_sampling_pes(filename_pes,pes)
  
  IF(desired_polar) &
       CALL write_sampling_polar(filename_polar, pad)
  
  IF(desired_amplitude) &
       CALL write_sampling_3Damplitude(filename_amplitude, pad)
  
  ! Deallocate the amplitude
  DEALLOCATE(pes,pad)
  
  WRITE(*,*) 'Our spectra is ready!!!'
  
END PROGRAM sampling_calculator
