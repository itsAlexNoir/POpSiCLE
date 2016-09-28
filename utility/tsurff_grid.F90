PROGRAM tsurff_grid

  USE popsicle_aux
  USE MPI
  
  IMPLICIT NONE
  
  !--Program variables-----------------------------------------------------!
  
  TYPE(COMMTYPE)             :: communicator      
  
  INTEGER                    :: maxkpts, maxthetapts, maxphipts
  INTEGER                    :: numkpts, lmax
  REAL(dp)                   :: deltaenergy
  INTEGER                    :: numproc1dk        
  INTEGER                    :: numproc1dt       
  INTEGER                    :: numproc1dp         
  INTEGER                    :: maxproc1dk         
  INTEGER                    :: maxproc1dt       
  INTEGER                    :: maxproc1dp         
  
  INTEGER                    :: numorbitals, iorbitaloffset
  INTEGER                    :: iproc, ierror, errorcode    
  
  INTEGER                    :: ik, itheta, iphi
  INTEGER                    :: ikg, ithetag, iphig
  INTEGER                    :: iprock, iproct, iprocp
  INTEGER                    :: ipro, iorbital, ios
  COMPLEX(dp), ALLOCATABLE   :: amplitude(:, :, :)
  COMPLEX(dp), ALLOCATABLE   :: sph_harmonics(:, :, :, :)
  REAL(dp), ALLOCATABLE      :: probktpr(:, :, :)
  REAL(dp), ALLOCATABLE      :: probktpi(:, :, :)
  REAL(dp), ALLOCATABLE      :: probktppror(:, :, :)
  REAL(dp), ALLOCATABLE      :: probktpproi(:, :, :)
  LOGICAL                    :: radamp_desired
  LOGICAL                    :: pes_desired
  LOGICAL                    :: polar_desired
  LOGICAL                    :: radamp_angbasis_desired
  CHARACTER(LEN = 3)         :: corbital
  CHARACTER(LEN = 4)         :: cprocessor
  CHARACTER(LEN = 100)       :: data_directory
  CHARACTER(LEN = 100)       :: filename 
  CHARACTER(LEN = 100)       :: filename_part
  CHARACTER(LEN = 100)       :: inputfile
  CHARACTER(LEN = 100)       :: radamp_filename
  CHARACTER(LEN = 100)       :: pes_filename
  CHARACTER(LEN = 100)       :: polar_filename
  CHARACTER(LEN = 100)       :: radamp_angbasis_filename
  
  !------------------------------------------------------------------------!
  
  ! Start up MPI
  
  CALL MPI_init( ierror )
  
  ! Find out number of processors.
  
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD, communicator%numprocessors,           &
       ierror )
  
  !communicator%maxprocessor = communicator%numprocessors - 1
  
  ! Find out number of the processor we are working on.
  
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, communicator%iprocessor, ierror )
  
  !-----------------------------------------------------------------------!
  !
  ! Read input file
  !
  
  ! Get input file
  CALL get_command_argument(1, inputfile)
  IF (LEN_TRIM(inputfile) == 0) THEN
     IF(communicator%iprocessor .EQ. 0) &
          WRITE(*, *) 'No input file specified.'
     
     CALL MPI_ABORT( MPI_COMM_WORLD, errorcode, ierror )
  ENDIF
  
  ! Open input file
  OPEN( UNIT = 12, FORM = 'formatted', FILE = inputfile, &
       STATUS= 'OLD', IOSTAT = ios, ACTION = 'READ')
  IF(ios .NE. 0) THEN
     IF(communicator%iprocessor .EQ. 0) &
          WRITE(*, *) 'Error opening input file.'
     
     CALL MPI_ABORT( MPI_COMM_WORLD, errorcode, ierror )
  ENDIF
  
  READ(12, *)
  READ(12, *)
  READ(12, *) data_directory
  READ(12, *)
  READ(12, *) !filename_surf_part
  READ(12, *)
  READ(12, *) filename_part
  READ(12, *)
  READ(12, *) maxkpts, maxthetapts, maxphipts
  READ(12, *)
  READ(12, *) numproc1dk, numproc1dt, numproc1dp
  READ(12, *)
  READ(12, *) deltaenergy
  READ(12, *)
  READ(12, *) !rmesh%rsurface
  READ(12, *)
  READ(12, *) lmax
  READ(12, *)
  READ(12, *) !numorbitals
  READ(12, *)
  READ(12, *) !minorbital
  READ(12, *)
  READ(12, *) !maxorbital
  READ(12, *)
  READ(12, *) iorbitaloffset
  READ(12, *)
  READ(12, *) !rmesh%gauge_transform
  READ(12, *)
  READ(12, *) !tmesh%truncate_time_integration
  READ(12, *)
  READ(12, *) !tmesh%maxtime_integrate
  READ(12, *)
  READ(12, *) !tmesh%timeaverageforrydberg
  READ(12, *)
  READ(12, *) !tmesh%timetoaverageover
  READ(12, *)
  READ(12, *) radamp_desired, radamp_filename
  READ(12, *)
  READ(12, *) pes_desired, pes_filename
  READ(12, *)
  READ(12, *) polar_desired, polar_filename
  READ(12, *)
  READ(12, *) radamp_angbasis_desired, radamp_angbasis_filename
  
  
  CLOSE( UNIT = 12 )

  !-----------------------------------------------------------------------!
  
  numkpts        = maxkpts * numproc1dk
  numthetapts    = maxthetapts * numproc1dt
  numphipts      = maxphipts * numproc1dp
  
  maxproc1dk     = numproc1dk - 1
  maxproc1dt     = numproc1dt - 1
  maxproc1dp     = numproc1dp - 1 
  
  numorbitals    = communicator%numprocessors

  ! Initialize dummy communications (for k coordinate parameters)
  commsurff%numproc1dk     = 1
  commsurff%numproc1dt     = 1
  commsurff%numproc1dp     = 1
  
  commsurff%maxproc1dk     = commsurff%numproc1dk - 1
  commsurff%maxproc1dt     = commsurff%numproc1dt - 1
  commsurff%maxproc1dp     = commsurff%numproc1dp - 1 
  
  CALL make_kpoints( lmax, numkpts, numthetapts, numphipts, deltaenergy )

  ALLOCATE(probktppror(1:maxkpts, 1:maxthetapts, 1:maxphipts))
  ALLOCATE(probktpproi(1:maxkpts, 1:maxthetapts, 1:maxphipts))
  ALLOCATE(probktpr(1:numkpts, 1:numthetapts, 1:numphipts))
  ALLOCATE(probktpi(1:numkpts, 1:numthetapts, 1:numphipts))
  
  iorbital = communicator%iprocessor + 1 + iorbitaloffset
  
  WRITE(corbital, '(I3.3)') iorbital
  
  ipro = -1
  
  DO iprocp = 0, maxproc1dp    
     DO iproct = 0, maxproc1dt    
        DO iprock = 0, maxproc1dk	 
           
           ipro = ipro + 1
           
           WRITE(cprocessor, '(I4.4)') ipro
           
           filename = TRIM(data_directory) // '/' //                       &
                TRIM(filename_part) // '.' //  corbital // '.' //    &
                cprocessor // '.dat' 
           
           OPEN(UNIT = 22, FORM = 'unformatted', FILE = filename)
           
           READ(22) probktppror
           READ(22) probktpproi
           
           CLOSE(UNIT = 22)
           
           DO iphi = 1, maxphipts
	          
              iphig = maxphipts * iprocp + iphi
              
              DO itheta = 1, maxthetapts 
                 
                 ithetag = maxthetapts * iproct + itheta
                 
                 DO ik = 1, maxkpts 
                    
                    ikg = maxkpts * iprock + ik
                    
                    probktpr(ikg, ithetag, iphig) =                        &
                         probktppror(ik, itheta, iphi)
                    
                    probktpi(ikg, ithetag, iphig) =                        &
                         probktpproi(ik, itheta, iphi)

                 ENDDO
              ENDDO
           ENDDO
           
        ENDDO
     ENDDO
  ENDDO


  !---------------------------------------------------------------------!
  !
  ! Output observables
  !
  
  ALLOCATE(amplitude(1:numkpts, 1:numthetapts, 1:numphipts))
  amplitude = CMPLX(probktpr,probktpi)

  IF (radamp_desired) THEN
     filename = TRIM(data_directory) // '/' // TRIM(radamp_filename) // '.' //  &
          corbital
     
     CALL write_radial_amplitude(amplitude, filename)

     CALL write_momwave(amplitude, filename)

  ENDIF
  
  IF (radamp_angbasis_desired) THEN
     
     ALLOCATE(sph_harmonics(1:kmesh%maxphipts, kmesh%maxthetapts, &
          -kmesh%lmax:kmesh%lmax, 0:kmesh%lmax))
     
     CALL create_spherical_harmonics(sph_harmonics, kmesh%lmax, &
          kmesh%costhetapts, kmesh%phipts)
     
     filename = TRIM(data_directory) // '/' // TRIM(radamp_angbasis_filename) // '.' //  &
          corbital
     
     CALL write_radial_amp_lm(amplitude, filename, sph_harmonics)

     DEALLOCATE(sph_harmonics)
     
  ENDIF

  IF (pes_desired) THEN
     filename = TRIM(data_directory) // '/' // TRIM(pes_filename) // '.' //  &
          corbital
     
     CALL write_pes(amplitude, filename)
     
  ENDIF
  
  IF (polar_desired) THEN
     filename = TRIM(data_directory) // '/' // TRIM(polar_filename) // '.' //  &
          corbital 
     
     CALL write_polar_amplitude(amplitude, filename)
     
  ENDIF
    
  !--------------------------!
  ! Release resources
  
  DEALLOCATE(probktppror)
  DEALLOCATE(probktpproi)
  DEALLOCATE(probktpr)
  DEALLOCATE(probktpi)
  DEALLOCATE(amplitude)
  
  CALL MPI_finalize( ierror )
  
END PROGRAM tsurff_grid


