PROGRAM tsurff_test
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax, mmin, mmax
  INTEGER                    :: lmax_total
  INTEGER                    :: numkpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts 
  REAL(dp)                   :: k_cutoff
  COMPLEX(dp), ALLOCATABLE   :: b(:, :, :)
  CHARACTER(LEN=100)         :: filename_surf
  CHARACTER(LEN=100)         :: filename_pes
  CHARACTER(LEN=100)         :: filename_mes
  CHARACTER(LEN=100)         :: filename_polar
  CHARACTER(LEN=100)         :: filename_amplitude
  CHARACTER(LEN=1)           :: answer
  
  LOGICAL                    :: desired_mes
  LOGICAL                    :: desired_polar
  LOGICAL                    :: desired_amplitude
  
  
  !--------------------------------------!
  
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
  READ(*, '(1A80)') filename_surf
  
  WRITE(*, *) 'Name of the PES file:'
  READ(*,'(1A80)') filename_pes

  WRITE(*, *) 'Boundary radius:'
  READ(*, *) radius_boundary

  WRITE(*, *) 'Cut-off frequency k:'
  READ(*, *) k_cutoff

  WRITE(*, *) 'Maximum angular momentum:'
  READ(*, *) lmax
  
  WRITE(*, *) 'Do you want momentum electron? (y or n)'
  READ(*, '(1A1)') answer

  desired_mes = .FALSE.

  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     desired_mes = .TRUE.

     WRITE(*, *) 'Name of the MES file:'
     READ(*,'(1A80)') filename_mes
     
  ENDIF

  WRITE(*, *) 'Do you want the 3D spherical scattering amplitude? (y or n)'
  READ(*, '(1A1)') answer
  
  desired_amplitude = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     desired_amplitude = .TRUE.

     WRITE(*, *) 'Name of the spherical amplitude file:'
     READ(*,'(1A80)') filename_amplitude
     
  ENDIF
  
  WRITE(*, *) 'Do you want polar scattering amplitude? (y or n)'
  READ(*, '(1A1)') answer
  
  desired_polar = .FALSE.
  
  IF ((answer .EQ. 'Y') .OR. (answer .EQ. 'y')) THEN
     
     desired_polar = .TRUE.

     WRITE(*, *) 'Name of the polar amplitude file:'
     READ(*,'(1A80)') filename_polar
     
  ENDIF
  
  
  CALL initialize_tsurff(filename_surf, radius_boundary, lmax, &
       k_cutoff, numkpts, numthetapts, numphipts, &
       lmax_total, mmax )
  
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
