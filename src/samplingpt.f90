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
  COMPLEX(dp), ALLOCATABLE      :: psi_sph(:, :, :), psip_sph(:, :, :)
  COMPLEX(dp), ALLOCATABLE      :: psi_sph_v(:, :)
  REAL(dp)                      :: rb
  INTEGER                       :: numpts, numrpts, numwpts
  INTEGER                       :: numthetapts, numphipts
  INTEGER                       :: numtotalsamplingpts
  LOGICAL                       :: aftercycles
  LOGICAL                       :: gauge_trans_desired
  
CONTAINS
 
  !**************************************************************************!
  !**************************************************************************!
  !**************************************************************************!
  !**************************************************************************!

  SUBROUTINE initialize_sampling_points(filename, radb, maxtotalsamplingpts, &
       maxthetapts, maxphipts, timeafter_fs, gauge_trans_on, coulomb_exp, maxwpts )
    
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN)     :: filename
    REAL(dp), INTENT(IN)             :: radb
    INTEGER, INTENT(OUT)             :: maxtotalsamplingpts
    INTEGER, INTENT(OUT)             :: maxthetapts
    INTEGER, INTENT(OUT)             :: maxphipts
    REAL(dp), INTENT(IN), OPTIONAL   :: coulomb_exp
    REAL(dp), INTENT(IN), OPTIONAL   :: timeafter_fs
    LOGICAL, INTENT(IN), OPTIONAL    :: gauge_trans_on
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

    ! If this option is on, we will transform
    ! from length to velocity gauge the wavefunction
    IF(gauge_trans_on) THEN
       gauge_trans_desired = .TRUE.
    ELSE
       gauge_trans_desired = .FALSE.
    ENDIF
    
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
    ALLOCATE(psip_sph(1:numthetapts,1:numphipts,1:ntime))

    IF(gauge_trans_desired) THEN
       ALLOCATE(psi_sph_v(1:numthetapts,1:numphipts))
    ENDIF
    
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
            psi_sph(:,:,itime), psip_sph(:,:,itime), &
            time(itime), efield(:,itime), afield(:,itime) )
       
       IF(gauge_trans_desired) THEN
          CALL gauge_transform_l2v(psi_sph(:,:,itime), psi_sph_v, afield(:,itime))
          psi_sph(:,:,itime)  = psi_sph_v
       ENDIF
       
    ENDDO
    
    !-----------------------------!
    
    ALLOCATE(before(ntimetotal))
    ALLOCATE(after(ntimetotal))
    ALLOCATE(modsq(ntimetotal))
    before = ZERO
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          
          DO itime = 1, ntime
             before(itime) = psi_sph(itheta,iphi,itime) * &
                  EXP(ZIMAGONE * coulomb_exp_ener * time(itime))
          ENDDO
             
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
    REAL(dp), ALLOCATABLE            :: wavederivvstime(:)    
    REAL(dp), ALLOCATABLE            :: current(:)
    CHARACTER(LEN=100)               :: name
    INTEGER                          :: itime
    
    !----------------------------------------------!
    
    ALLOCATE(wavevstime(1:ntime),wavederivvstime(1:ntime),current(1:ntime))
    wavevstime = REAL( CONJG(psi_sph(1,1,:)) * psi_sph(1,1,:),dp)
    wavederivvstime = REAL( CONJG(psip_sph(1,1,:)) * psip_sph(1,1,:),dp)
    current = -ZIMAGONE * 0.5_dp * (psi_sph(1,1,:) * CONJG(psip_sph(1,1,:)) - &
         CONJG(psi_sph(1,1,:)) * psi_sph(1,1,:))
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=23,FORM='formatted',FILE=name)
    name = TRIM(filename) // '_psi.dat'
    OPEN(UNIT=33,FORM='formatted',FILE=name)
    name = TRIM(filename) // '_psip.dat'
    OPEN(UNIT=43,FORM='formatted',FILE=name)
    name = TRIM(filename) // '_current.dat'
    OPEN(UNIT=53,FORM='formatted',FILE=name)
    
    DO itime = 1, ntime
       WRITE(23,*) time(itime), wavevstime, wavederivvstime
       WRITE(33,*) time(itime),REAL(psi_sph(1,1,itime)),AIMAG(psi_sph(1,1,itime))
       WRITE(43,*) time(itime),REAL(psip_sph(1,1,itime)),AIMAG(psip_sph(1,1,itime))
       WRITE(53,*) time(itime), current(itime)
    ENDDO
    
    CLOSE(33)
    
    DEALLOCATE(wavevstime,wavederivvstime)
    
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
  
  SUBROUTINE gauge_transform_l2v(psi_length,psi_vel,afield)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)   :: psi_length(:, :)
    COMPLEX(dp), INTENT(OUT)  :: psi_vel(:, :)
    REAL(dp), INTENT(IN)      :: afield(:)
    
    INTEGER                   :: dims(2)
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
  
END MODULE samplingpt
