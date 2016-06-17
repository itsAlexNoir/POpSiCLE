PROGRAM tsurff_test
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax, mmin, mmax
  INTEGER                    :: lmax_total
  INTEGER                    :: numkpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts 
  REAL(dp)                   :: dk
  REAL(dp)                   :: k_cutoff
  REAL(dp)                   :: coulomb_exp_energy
  COMPLEX(dp), ALLOCATABLE   :: b(:, :, :)
  CHARACTER(LEN=100)         :: filename_surf
  CHARACTER(LEN=100)         :: filename_pes
  CHARACTER(LEN=100)         :: filename_mes
  CHARACTER(LEN=100)         :: filename_polar
  CHARACTER(LEN=100)         :: filename_amplitude
  CHARACTER(LEN=1)           :: answer

  LOGICAL                    :: desired_gauge_trans
  LOGICAL                    :: desired_mes
  LOGICAL                    :: desired_polar
  LOGICAL                    :: desired_amplitude
  
  !--------------------------------------!

  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '****************************'
  WRITE(*,*) '     PES calculator app     '
  WRITE(*,*) '****************************'
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Welcome to the PES calculator. This little program'
  WRITE(*,*) ' will help to calculate photoelectron spectra from '
  WRITE(*,*) ' a surface file (HDF5). '
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Initializing...'
  WRITE(*,*)

  
  WRITE(*, *) 'Path of the surface file:'
  WRITE(*, *) '-------------------------'
  READ(*, '(1A80)') filename_surf
  WRITE(*, *)
  WRITE(*, *) 'Filename: ',filename_surf
  WRITE(*, *)
  
  WRITE(*, *) 'Boundary radius (a.u.):'
  WRITE(*, *) '-----------------------'
  READ(*, *) radius_boundary
  WRITE(*, *)
  WRITE(*, '(A,F9.3)') 'Rb :',radius_boundary
  WRITE(*, *)
  
  WRITE(*, *) 'Momentum grid spacing (dk, in a.u.):'
  WRITE(*, *) '------------------------------------'
  READ(*, *) dk
  WRITE(*, *)
  WRITE(*, '(A,F9.3)') 'dk: ',dk
  WRITE(*, *)

  WRITE(*, *) 'Momentum cut-off (a.u.):'
  WRITE(*, *) '------------------------'
  READ(*, *) k_cutoff
  WRITE(*, *)
  WRITE(*, '(A,F9.3)') 'Cut-off: ',k_cutoff
  WRITE(*, *)

  WRITE(*, *) 'Maximum angular momentum:'
  WRITE(*, *) '-------------------------'
  READ(*, *) lmax
  WRITE(*, *)
  WRITE(*, *) 'lmax :',lmax
  WRITE(*, *)
    
  WRITE(*, *) 'Do you want to add the Coulomb explosion energy? (y or n)'
  WRITE(*, *) '---------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  coulomb_exp_energy = 0.0_dp
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN

     WRITE(*, *)
     WRITE(*, *) 'Yes'
     WRITE(*, *)
     WRITE(*, *) 'Coulomb explosion energy (a.u.):'
     WRITE(*, *) '--------------------------------'
     READ(*, *) coulomb_exp_energy
     WRITE(*, *)
     WRITE(*, '(A,F9.6)') 'Coulomb explosion energy :',coulomb_exp_energy
     WRITE(*, *)
     
  ELSE
     
     WRITE(*, *)
     WRITE(*, *) 'No'
     WRITE(*, *)
     
  ENDIF

  WRITE(*, *) 'Do you want to gauge transform from length to velocity (y or n)'
  WRITE(*, *) '---------------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  desired_gauge_trans = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN

     WRITE(*, *)
     WRITE(*, *) 'Yes'
     WRITE(*, *)
     desired_gauge_trans = .TRUE.
     
  ELSE
     WRITE(*, *)
     WRITE(*, *) 'No'
     WRITE(*, *)
     
  ENDIF
  
  
  WRITE(*, *) 'Name of the PES file:'
  WRITE(*, *) '---------------------'
  READ(*,'(1A80)') filename_pes
  WRITE(*, *)
  WRITE(*, *) 'Filename for PES :',filename_pes
  WRITE(*, *)
  
  WRITE(*, *) 'Do you want momentum electron? (y or n)'
  WRITE(*, *) '---------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)

  desired_mes = .FALSE.

  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN

     WRITE(*, *)
     WRITE(*, *) 'Yes'
     WRITE(*, *)
     
     desired_mes = .TRUE.

     WRITE(*, *) 'Name of the MES file:'
     WRITE(*, *) '---------------------'
     READ(*,'(1A80)') filename_mes
     WRITE(*, *)
     WRITE(*, *) 'Filename for MES :',filename_mes
     WRITE(*, *)

  ELSE

     WRITE(*, *)
     WRITE(*, *) 'No'
     WRITE(*, *)
     
  ENDIF

  WRITE(*, *) 'Do you want the 3D spherical scattering amplitude? (y or n)'
  WRITE(*, *) '-----------------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)

  desired_amplitude = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     WRITE(*, *)
     WRITE(*, *) 'Yes'
     WRITE(*, *)

     desired_amplitude = .TRUE.

     WRITE(*, *) 'Name of the spherical amplitude file:'
     WRITE(*, *) '-------------------------------------'
     READ(*,'(1A80)') filename_amplitude
     WRITE(*, *)
     WRITE(*, *) 'Filename for amplitude: ',filename_amplitude
     WRITE(*, *)

  ELSE

     WRITE(*, *)
     WRITE(*, *) 'No'
     WRITE(*, *)

  ENDIF
  
  WRITE(*, *) 'Do you want polar scattering amplitude? (y or n)'
  WRITE(*, *) '------------------------------------------------'
  READ(*, '(1A1)') answer
  WRITE(*, *)
  
  desired_polar = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN

     WRITE(*, *)
     WRITE(*, *) 'Yes'
     WRITE(*, *)
     
     desired_polar = .TRUE.
     
     WRITE(*, *) 'Name of the polar amplitude file:'
     WRITE(*, *) '---------------------------------'
     READ(*,'(1A80)') filename_polar
     WRITE(*, *)
     WRITE(*, *) 'Filename for polar amplitude :',filename_polar
     WRITE(*, *)

  ELSE

     WRITE(*, *)
     WRITE(*, *) 'No'
     WRITE(*, *)
     
  ENDIF
  
  
  CALL initialize_tsurff(filename_surf, radius_boundary, lmax, &
       dk, k_cutoff, numkpts, numthetapts, numphipts, &
       lmax_total, mmax, desired_gauge_trans, coulomb_exp_energy )
  
  mmin = - mmax
  
  ! Allocate momentum amplitude array
  ALLOCATE(b(1:numkpts,1:numthetapts,1:numphipts))
  
  ! The name is self-explanatory
  CALL get_flux(filename_surf, lmax, mmax, b )

  !-------------------!
  !    OUTPUT TIME!   !
  !-------------------!

  CALL write_pes(b,filename_pes)
  
  IF(desired_polar) &
       CALL write_polar_amplitude(b, filename_polar)
  
  IF(desired_amplitude) &
       CALL write_amplitude(b, filename_amplitude)
  
  IF(desired_mes) &
       CALL write_mes(b,filename_mes)
  
  ! Deallocate the amplitude
  DEALLOCATE(b)
  
  WRITE(*,*) 'Our spectra is ready!!!'

END PROGRAM tsurff_test
