MODULE samplingpt

  USE constants
  USE gaussleg
  USE fourier
  USE io_surface
  USE io_pop
  
  IMPLICIT NONE
  
  PRIVATE
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public subroutines                                                        !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  PUBLIC                        :: initialize_sampling_points
  PUBLIC                        :: calculate_sampling_pes
  PUBLIC                        :: get_wave_time_file
  PUBLIC                        :: write_sampling_pes
  PUBLIC                        :: write_sampling_polar
  PUBLIC                        :: write_sampling_3Damplitude
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public variables                                                          !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Private variables                                                         !
  !                                                                           !
  !---------------------------------------------------------------------------!

  ! Public variables
  REAL(dp), PUBLIC, ALLOCATABLE :: w(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: theta_ax(:), costheta_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: phi_ax(:)
  REAL(dp), PUBLIC, ALLOCATABLE :: gauss_th_weights(:)
  
  ! Private variables
  INTEGER                       :: numkpts
  REAL(dp)                      :: dt, dw, wmax
  REAL(dp)                      :: coulomb_exp_ener
  INTEGER                       :: ntime, ntimeafter, ntimetotal
  REAL(dp), ALLOCATABLE         :: time(:), efield(:, :), afield(:, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :, :), psip_sph(:, :)
  REAL(dp)                      :: rb
  INTEGER                       :: numpts, numrpts, numwpts
  INTEGER                       :: numthetapts, numphipts
  INTEGER                       :: numtotalsamplingpts
  LOGICAL                       :: aftercycles
  
CONTAINS
 
  !**************************************************************************!
  !**************************************************************************!
  !**************************************************************************!
  !**************************************************************************!

  SUBROUTINE initialize_sampling_points(filename, radb, maxtotalsamplingpts, &
       maxthetapts, maxphipts, timeafter_fs, coulomb_exp, maxwpts )
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    REAL(dp), INTENT(IN)             :: radb
    INTEGER, INTENT(OUT)             :: maxtotalsamplingpts
    INTEGER, INTENT(OUT)             :: maxthetapts
    INTEGER, INTENT(OUT)             :: maxphipts
    REAL(dp), INTENT(IN), OPTIONAL   :: coulomb_exp
    REAL(dp), INTENT(IN), OPTIONAL   :: timeafter_fs
    INTEGER, INTENT(OUT), OPTIONAL   :: maxwpts
    
    INTEGER                          :: lmax
    INTEGER                          :: iphi, iw, il
    
    !-------------------------------------------------!

    ! Find out the number of time steps in the file
    CALL get_surface_dims(filename, ntime, dt, &
         numthetapts, numphipts, lmax)
    
    maxthetapts = numthetapts
    maxphipts = numphipts
    numtotalsamplingpts = numthetapts * numphipts
    maxtotalsamplingpts = numtotalsamplingpts
    
    ! If we want to padd with zeros at the end
    ! of the time evolution
    ntimeafter = 0.0_dp
    aftercycles = .FALSE.
    IF(PRESENT(timeafter_fs)) THEN
       aftercycles = .TRUE.
       ntimeafter = INT(timeafter_fs / autime_fs / dt)
    ENDIF
    ! Total time for the calculation
    ntimetotal = ntime + ntimeafter
    
    ! This factor shift the energy if we
    ! are dealing with fixed nuclei molecules
    coulomb_exp_ener = 0.0_dp
    IF(PRESENT(coulomb_exp)) &
         coulomb_exp_ener = coulomb_exp

    ! Maximum energy
    dw = twopi / (ntimetotal * dt)
    numwpts = INT((ntimetotal-1)*0.5+1)
    wmax = (numwpts - 1) * dw
    IF(PRESENT(maxwpts)) &
         maxwpts = numwpts
    
    WRITE(*,'(A,I3.3)')  'Total number of sampling points: ',numtotalsamplingpts
    WRITE(*,'(A,I3.3)')  'Number of points in theta:       ',numthetapts
    WRITE(*,'(A,I3.3)')  'Number of points in phi:         ',numphipts
    
    WRITE(*,'(A,F10.6)')  'Maximum energy:                  ',wmax
    WRITE(*,'(A,F9.6)')  'dw:                               ',dw
    WRITE(*,'(A,F0.6)')  'Time in the file (a.u.):          ', ntime * dt
    WRITE(*,'(A,F0.6)')  'Time in the file (fs):            ', ntime * dt * autime_fs
    IF(aftercycles) THEN     
       WRITE(*,'(A,F0.6)')  'Time added after the pulse (a.u.): ', ntimeafter * dt
       WRITE(*,'(A,F0.6)')  'Time added after the pulse (fs):   ', timeafter_fs
    ENDIF
    WRITE(*,'(A,F0.6)')  'Total time (a.u.):                 ', ntimetotal * dt
    WRITE(*,'(A,F0.6)')  'Total time (fs):                   ', ntimetotal * dt * autime_fs
    
    IF(PRESENT(coulomb_exp)) &
         WRITE(*,'(A,F9.6)')  'Coulomb explosion energy (a.u.): ', coulomb_exp_ener
    WRITE(*,*)           '---------------------------------------------------------&
         &---------------------'
    WRITE(*,*)
    
    ! Assign  surface radius
    rb = radb

    ! Create energy array
    ALLOCATE(w(1:numwpts))
    w = 0.0_dp
    DO iw = 1, numwpts
       w(iw) = REAL(iw-1,dp) * dw
    ENDDO
    
    ! Allocate wavefunction arrays and rest of stuff
    ALLOCATE(psi_sph(1:numthetapts,1:numphipts,1:ntime))
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts))
    ALLOCATE(time(1:ntime))
    ALLOCATE(efield(1:3,1:ntime),afield(1:3,1:ntime))
    
    ! Create spherical coordinates
    ALLOCATE(theta_ax(1:numthetapts),costheta_ax(1:numthetapts))
    ALLOCATE(gauss_th_weights(1:numthetapts))
    ALLOCATE(phi_ax(1:numphipts))
    
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta_ax, gauss_th_weights)
    theta_ax = ACOS(costheta_ax)
    
    DO iphi = 1, numphipts
       phi_ax(iphi) = REAL(iphi,dp) * twopi / REAL(numphipts,dp)
    ENDDO
    
    
  END SUBROUTINE initialize_sampling_points
  
    !**************************************************************************!
  
  SUBROUTINE calculate_sampling_pes(filename,pes,pad)
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    REAL(dp), INTENT(OUT)            :: pes(:)
    REAL(dp), INTENT(OUT), OPTIONAL  :: pad(:, :, :)
    
    COMPLEX(dp), ALLOCATABLE         :: before(:), after(:)
    REAL(dp), ALLOCATABLE            :: modsq(:)
    REAL(dp), ALLOCATABLE            :: solid_angle(:, :)
    INTEGER                          :: itheta, iphi
    INTEGER                          :: iw, itime
    
    !------------------------------------------------------!

    pes = 0.0_dp
    pad = 0.0_dp

    ! Read quantities from file
    DO itime = 1, ntime
       
       ! Read wavefunction at Rb from file
       CALL read_surface(filename, itime - 1, numthetapts, numphipts, &
            psi_sph(:,:,itime), psip_sph, &
            time(itime), efield(:,itime), afield(:,itime) )
       
    ENDDO

    ! If we want to pad with zero the end of the array.
    IF(aftercycles) THEN
       DO itime = ntime + 1, ntimetotal
          psi_sph(:,:,itime) = ZERO
       ENDDO
    ENDIF

    !-----------------------------!
    
    ALLOCATE(before(ntimetotal))
    ALLOCATE(after(ntimetotal))
    ALLOCATE(modsq(ntimetotal))

    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts

          before = psi_sph(itheta,iphi,:) * &
               EXP(ZIMAGONE * coulomb_exp_ener * time)

          CALL FourierTransform(before,after,1,&
               (/ntimetotal/),inverse=.TRUE.)
          
          after = after * dt / SQRT(twopi)
          
          DO iw = 1, numwpts
             modsq(iw) = REAL(CONJG(after(iw)) * &
                  after(iw), dp)
             
             IF(PRESENT(pad)) &
                  pad(itheta,iphi,iw) = modsq(iw)
             
             pes(iw) = pes(iw) + gauss_th_weights(itheta) * &
                  SQRT(w(iw)) * modsq(iw)

          ENDDO
       ENDDO
    ENDDO

    IF(numphipts.EQ.1) THEN
       pes = pes * twopi
    ELSE
       pes = pes *  twopi / REAL(numphipts,dp)
    ENDIF
       
    ! Release memory
    DEALLOCATE(before,after)
    DEALLOCATE(modsq)
    
  END SUBROUTINE calculate_sampling_pes

  !****************************************************************!

  SUBROUTINE get_wave_time_file(filename)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)     :: filename

    REAL(dp), ALLOCATABLE            :: wavevstime(:)
    CHARACTER(LEN=100)               :: name
    INTEGER                          :: itime
    
    !----------------------------------------------!

    ALLOCATE(wavevstime(1:ntimetotal))
    wavevstime = REAL( CONJG(psi_sph(1,1,:)) * psi_sph(1,1,:),dp)
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=33,FORM='formatted',FILE=name)
    
    DO itime = 1, ntimetotal
       WRITE(33,*) time(itime), wavevstime(itime)
    ENDDO
    
    CLOSE(33)

    DEALLOCATE(wavevstime)
    
  END SUBROUTINE get_wave_time_file
  
  !**************************************************************************!

  SUBROUTINE write_sampling_pes(filename, pes)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    REAL(dp), INTENT(IN)           :: pes(:)

    CHARACTER(LEN=100)             :: name
    INTEGER                        :: iw
    !---------------------------------------!
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
    DO iw = 1, numwpts
       WRITE(55,*) w(iw), pes(iw)
    ENDDO

    CLOSE(55)
    
  END SUBROUTINE write_sampling_pes

  !**************************************************************************!

  SUBROUTINE write_sampling_polar(filename, pad, groupname)
    
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)   :: filename
    REAL(dp), INTENT(IN)           :: pad(:, :, :)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname

    REAL(dp), ALLOCATABLE          :: prob2D(:, :)
    REAL(dp)                       :: suma
    INTEGER                        :: itheta, iphi, iw
    !------------------------------------------------------!    
    
    ALLOCATE(prob2D(1:numthetapts,1:numwpts))
    
    IF(numphipts.EQ.1) THEN
       prob2D = pad(:,1,:) * twopi
    ELSE

       DO iw = 1, numwpts
          DO itheta = 1, numthetapts
             suma = 0.0_dp
             DO iphi = 1, numphipts
                suma = suma + pad(itheta,iphi,iw)
             ENDDO
             prob2D(itheta,iw) = suma
          ENDDO
       ENDDO
       
       prob2D = prob2D *  twopi / REAL(numphipts,dp)
       
    ENDIF
    
    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(prob2D,(/numthetapts*numwpts/)),2,&
            (/numthetapts,numwpts/),filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(prob2D,(/numthetapts*numwpts/)),2,&
            (/numthetapts,numwpts/),filename)
    ENDIF
    
    DEALLOCATE(prob2D)
    
    
  END SUBROUTINE write_sampling_polar
    
  !**************************************************************************!

  SUBROUTINE write_sampling_3Damplitude(filename, pad, groupname)

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)           :: filename
    REAL(dp), INTENT(IN)                   :: pad(:,:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    !------------------------------------------------------!

    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(pad,(/numthetapts*numphipts*numwpts/)),2,&
            (/numthetapts,numphipts,numwpts/),filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(pad,(/numthetapts*numphipts*numwpts/)),2,&
            (/numthetapts,numphipts,numwpts/),filename)
    ENDIF
    
  END SUBROUTINE write_sampling_3Damplitude

  !**************************************************************************!

END MODULE samplingpt
