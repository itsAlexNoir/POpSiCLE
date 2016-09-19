PROGRAM tsurff_calculator
  
  USE popsicle
  USE MPI
  
  IMPLICIT NONE
  
  INTEGER                     :: numprock, numproctheta, numprocphi
  INTEGER                     :: maxkpts, maxthetapts, maxphipts
  REAL(dp)                    :: deltaenergy
  INTEGER                     :: numorbitals, iorbitaloffset
  INTEGER                     :: minorbital, maxorbital
  INTEGER                     :: iorbital, itime, ios
  INTEGER                     :: lmax
  
  CHARACTER(LEN=100)          :: filename_surf_part
  CHARACTER(LEN=100)          :: filename_surf
  
  CHARACTER(LEN=100)          :: filename_amplitude_part
  CHARACTER(LEN=100)          :: filename_amplitude
  
  CHARACTER(LEN=3)            :: corbital
  CHARACTER(LEN = 100)        :: inputfile
  CHARACTER(LEN = 100)        :: data_directory
  
  !---------------------------------------------------------------------------!   
  
  !
  ! Read input file
  !
  
  ! Get input file
  CALL get_command_argument(1, inputfile)
  IF (LEN_TRIM(inputfile) == 0) &
       CALL tsurff_stop('No input file specified.')
  
  ! Open input file
  OPEN( UNIT = 12, FORM = 'formatted', FILE = inputfile, &
       STATUS= 'OLD', IOSTAT = ios, ACTION = 'READ')
  IF(ios .NE. 0) &
       CALL tsurff_stop('Error opening input file.')
  
  READ(12, *)
  READ(12, *)
  READ(12, *) data_directory
  READ(12, *)
  READ(12, *) filename_surf_part
  READ(12, *)
  READ(12, *) filename_amplitude_part
  READ(12, *)
  READ(12, *) maxkpts, maxthetapts, maxphipts
  READ(12, *)
  READ(12, *) numprock, numproctheta, numprocphi
  READ(12, *)
  READ(12, *) deltaenergy
  READ(12, *)
  READ(12, *) rmesh%rsurface
  READ(12, *)
  READ(12, *) lmax
  READ(12, *)
  READ(12, *) numorbitals
  READ(12, *)
  READ(12, *) minorbital
  READ(12, *)
  READ(12, *) maxorbital
  READ(12, *)
  READ(12, *) iorbitaloffset
  READ(12, *)
  READ(12, *) rmesh%gauge_transform
  READ(12, *)
  READ(12, *) tmesh%truncate_time_integration
  READ(12, *)
  READ(12, *) tmesh%maxtime_integrate
  READ(12, *)
  READ(12, *) tmesh%timeaverageforrydberg
  READ(12, *)
  READ(12, *) tmesh%timetoaverageover
  
  CLOSE( UNIT = 12 )
  
  CALL initialize_communications( numprock, numproctheta, numprocphi )
  
  CALL output_communications( )
  
  CALL check_processor_count( )
  
  CALL initialize_ip_addresses( )
  
  CALL make_kpoints( lmax, maxkpts, maxthetapts, maxphipts, deltaenergy )
  
  IF (commsurff%iprocessor.EQ.0) THEN
     
     WRITE(*, *)
     WRITE(*, *)
     WRITE(*, *) '****************************'
     WRITE(*, *) '     PES calculator app     '
     WRITE(*, *) '****************************'
     WRITE(*, *)
     WRITE(*, *)
     WRITE(*, *) ' Welcome to the PES calculator. This little program'
     WRITE(*, *) ' will help to calculate photoelectron spectra from '
     WRITE(*, *) ' a surface file (HDF5). '
     WRITE(*, *)
     WRITE(*, *)
     WRITE(*, *) ' Initializing...'
     WRITE(*, *)
     
     WRITE(*, *)
     
     WRITE(*, *)					                       &
          '--t-SURFF Parameters--------------------------------------------'
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, I9)')   ' Number of orbitals:                         ', &
          numorbitals
     WRITE(*, '(A46, I9)')   ' Minimum Orbital number:                     ', &
          minorbital
     WRITE(*, '(A46, I9)')   ' Maximum Orbital number:                     ', &
          maxorbital
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, F9.3)') ' Maximum energy value:                       ', &
          kmesh%emax
     WRITE(*, '(A46, F9.3)') ' Energy spacing:                             ', &
          kmesh%deltaenergy
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, I9)') ' Number of local theta plotting points:      ',   &
          kmesh%maxthetapts 
     WRITE(*, '(A46, I9)') ' Number of global theta plotting points:     ',   &
          kmesh%numthetapts 
     WRITE(*, *) 
     
     WRITE(*, '(A46,I9 )') ' Number of local phi plotting points:        ',   &
          kmesh%maxphipts 
     WRITE(*, '(A46,I9 )') ' Number of global phi plotting points:       ',   &
          kmesh%numphipts 
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, F9.3)') ' Surface radius:                             ', &
          rmesh%rsurface
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, L9)')   ' Are we truncating the time integration?:    ', &
          tmesh%truncate_time_integration
     
     IF (tmesh%truncate_time_integration) THEN
        
        WRITE(*, '(A46, I9)') ' Maximum time point for integration:         ',&
             tmesh%maxtime_integrate 
        
     ENDIF
     
     WRITE(*, *) 
     
     WRITE(*, '(A46, L9)')   ' Transform from length to velocity gauge?:   ', &
          rmesh%gauge_transform
     WRITE(*, *) 
     
     WRITE(*, *)					                       &
          '----------------------------------------------------------------'
     
     WRITE(*, *) 
     
  ENDIF

  iorbital = 1
  WRITE(corbital, '(I3.3)') iorbital
  
  filename_surf        = TRIM(data_directory) // '/' //	               &
       TRIM(filename_surf_part) // '.' //                &
       corbital // '.d'
  
  CALL initialize_tsurff( filename_surf  )

  DO iorbital = minorbital, maxorbital
     
     WRITE(corbital, '(I3.3)') iorbital
     
     filename_surf        = TRIM(data_directory) // '/' //	               &
          TRIM(filename_surf_part) // '.' //                &
          corbital // '.d'
     
     filename_amplitude   = TRIM(data_directory) // '/' //	               &
          TRIM(filename_amplitude_part) // '.' //           &
          corbital // '.' // 		               &
          commsurff%cprocessor // '.dat' 
     
     CALL get_flux( filename_surf, iorbital )
     
     IF (commsurff%iprocessor.EQ.0) THEN
        
        WRITE(*, *)					                       &
             '--Outputting spectral information-------------------------------'
        
        WRITE(*, *)
        
        WRITE(*, '(A46, A100)')                                               &
             ' PES file:                                   ',                 &
             filename_amplitude
        
     ENDIF
     
     OPEN (UNIT = 20, FORM = 'unformatted', FILE = filename_amplitude)
     
     probk3D = REAL(bkaverage, dp)
     
     WRITE(20) probk3D
     
     probk3D = AIMAG(bkaverage)
     
     WRITE(20) probk3D
     
     CLOSE(UNIT = 20)
     
  ENDDO
  
  
  IF (commsurff%iprocessor.EQ.0) THEN
     
     WRITE(*, *) 
     
     WRITE(*, *)					                       &
          '--End of calculation--------------------------------------------'
     
     WRITE(*, *)
     
  ENDIF
  
  ! Finalize the communications
  
  CALL finalize_communications( )
  
  
 END PROGRAM tsurff_calculator
 
