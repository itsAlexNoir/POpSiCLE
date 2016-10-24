!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: observables
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 05/02/2016
!
! DESCRIPTION:
!> \brief observables module
!> \details This module contains routines to calculate various observables
!> from the outcome of a photoionisation.
!
!------------------------------------------------------------------------------
MODULE observables
  
  USE constants_pop
  USE comm_surff
  USE gaussleg
  USE sharmonics
  USE flux
  USE io_pop
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC        :: get_angular_resolution
  PUBLIC        :: get_amplitude3D
  PUBLIC        :: get_polar_amplitude
  PUBLIC        :: get_radial_amplitude
  PUBLIC        :: get_pes
  PUBLIC        :: get_momwave
  PUBLIC        :: get_radial_amp_lm
  PUBLIC        :: write_amplitude3D
  PUBLIC        :: write_polar_amplitude
  PUBLIC        :: write_radial_amplitude
  PUBLIC        :: write_pes
  PUBLIC        :: write_momwave
  PUBLIC        :: write_radial_amp_lm
  
CONTAINS
  
  SUBROUTINE get_angular_resolution(b_lm, dPomega, lmax, &
       sph_harmonics)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)           :: lmax
    COMPLEX(dp), INTENT(IN)       :: b_lm(:, -lmax:, 0:)
    COMPLEX(dp), INTENT(OUT)      :: dPomega(:, :, :)
    COMPLEX(dp), INTENT(IN)       :: sph_harmonics(:, :, -lmax:, 0:)
    
    INTEGER                       :: dims(3)
    INTEGER                       :: numkpts, numthetapts
    INTEGER                       :: numphipts
    INTEGER                       :: ik,il,im,itheta,iphi
    
    !---------------------------------------------------------!

    dims = SHAPE(dPomega)
    
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)
    
    dPomega = ZERO
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             DO il = 0, lmax
                DO im = -il, il
                   
                   dPomega(ik,itheta,iphi) = b_lm(ik,im,il) * &
                        sph_harmonics(im,il,iphi,itheta)
                   
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE get_angular_resolution
  
  !***********************************************************!

  SUBROUTINE get_radial_amplitude(bk, bk_rad)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_rad(:)

    REAL(dp)                   :: suma
    INTEGER                    :: ik, itheta, iphi

    !---------------------------------------!
    
    bk_rad = 0.0_dp
    
    DO ik = 1, kmesh%maxkpts
       suma = 0.0_dp
       DO iphi = 1, kmesh%maxphipts
          DO itheta = 1, kmesh%maxthetapts
             suma = suma + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * &
                  bk(ik,itheta,iphi),dp) * &
                  kmesh%wtheta(itheta) * &
                  kmesh%wphi(iphi)
          ENDDO
       ENDDO
       bk_rad(ik) = suma
    ENDDO
    
    
  END SUBROUTINE get_radial_amplitude
  
  !***********************************************************!
  
  SUBROUTINE write_radial_amplitude(bk, filename, groupname)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    REAL(dp), ALLOCATABLE          :: probk1D(:)
    CHARACTER(LEN=100)             :: name
    INTEGER                        :: ik, itheta, iphi
    !--------------------------------------------------!
    
    ALLOCATE(probk1D(1:kmesh%maxkpts))
    
    CALL get_radial_amplitude(bk,probk1D)
    
    ! IF(PRESENT(groupname)) THEN
    !    CALL write_wave(probk1D,1,(/numkpts/),filename,&
    !         groupname,groupname)
    ! ELSE
    !    CALL write_wave(probk1D,1,(/numkpts/),filename)
    ! ENDIF
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
    DO ik = 1, kmesh%maxkpts
       WRITE(55,*) kmesh%kpts(ik), probk1D(ik)
    ENDDO
    
    CLOSE(55)
    
    DEALLOCATE(probk1D)
    
  END SUBROUTINE write_radial_amplitude
  
  !***********************************************************!

  SUBROUTINE get_pes(bk, bk_pes)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)      :: bk(:, :, :)
    REAL(dp), INTENT(OUT)        :: bk_pes(:)
    
    INTEGER                      :: ik, itheta, iphi
    !-------------------------------------------------!

    
    CALL get_radial_amplitude(bk, bk_pes)
    
    DO ik = 1, kmesh%maxkpts
       bk_pes(ik) = bk_pes(ik) * kmesh%kpts(ik)
    ENDDO
    
  END SUBROUTINE get_pes

  !***********************************************************!

  SUBROUTINE write_pes(bk, filename, groupname)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname

    REAL(dp), ALLOCATABLE          :: prob1D(:)
    CHARACTER(LEN=100)             :: name
    INTEGER                        :: ik, itheta, iphi
    !--------------------------------------------------!
    
    ALLOCATE(prob1D(1:kmesh%maxkpts))
    
    CALL get_pes(bk,prob1D)
    
    ! IF(PRESENT(groupname)) THEN
    !    CALL write_wave(probk1D,1,(/numkpts/),filename,&
    !         groupname,groupname)
    ! ELSE
    !    CALL write_wave(probk1D,1,(/numkpts/),filename)
    ! ENDIF
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
    DO ik = 1, kmesh%maxkpts
       WRITE(55,*) kmesh%energypts(ik), prob1D(ik)
    ENDDO
    
    CLOSE(55)
    
    DEALLOCATE(prob1D)

  END SUBROUTINE write_pes

  !***********************************************************!

  SUBROUTINE get_polar_amplitude(bk, bk_polar)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_polar(:, :)

    INTEGER                    :: ik, itheta, iphi

    !----------------------------------------------!
    
    bk_polar = 0.0_dp
    
    DO iphi = 1, kmesh%maxphipts
       DO itheta = 1, kmesh%maxthetapts
          DO ik = 1, kmesh%maxkpts
             bk_polar(ik,itheta) = bk_polar(ik,itheta) + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * &
                  bk(ik,itheta,iphi),dp) * &
                  kmesh%wphi(iphi)
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE get_polar_amplitude
  
  !***********************************************************!
  
  SUBROUTINE write_polar_amplitude(bk, filename, silofile, groupname)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)                :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    LOGICAL, OPTIONAL                      :: silofile
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    REAL(dp), ALLOCATABLE                  :: probk2D(:, :)
    INTEGER                                :: dims(2)
    
    !--------------------------------------------------!
    
    ALLOCATE(probk2D(1:kmesh%maxkpts,1:kmesh%maxthetapts))
    
    CALL get_polar_amplitude(bk,probk2D)
    
    dims = (/ kmesh%maxkpts, kmesh%maxthetapts /)
    
    IF(silofile) THEN
       CALL create_silo_file(filename)
       CALL write_silo_rectilinear_mesh2D(filename,dims, &
            kmesh%kpts, kmesh%thetapts )
       CALL write_silo_function2D(filename, dims, probk2D)
    ENDIF
    
    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(probk2D, &
            (/kmesh%maxkpts*kmesh%maxthetapts/)),2,&
            (/kmesh%maxkpts,kmesh%maxthetapts/), &
            filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(probk2D, &
            (/kmesh%maxkpts*kmesh%maxthetapts/)),2,&
            (/kmesh%maxkpts,kmesh%maxthetapts/),filename)
    ENDIF
    
    DEALLOCATE(probk2D)
    
  END SUBROUTINE write_polar_amplitude
  
  !***********************************************************!
  
  SUBROUTINE get_amplitude3D(bk, bk_real)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_real(:, :, :)
    
    REAL(dp)                   :: dphi
    INTEGER                    :: numkpts, numthetapts
    INTEGER                    :: numphipts, dims(3)
    INTEGER                    :: ik, itheta, iphi
    
    !----------------------------------------------!
    
    bk_real = 0.0_dp
    
    DO iphi = 1, kmesh%maxphipts
       DO itheta = 1, kmesh%maxthetapts
          DO ik = 1, kmesh%maxkpts
             bk_real(ik,itheta,iphi) = bk_real(ik,itheta,iphi) + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * &
                  bk(ik,itheta,iphi),dp)
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE get_amplitude3D
  
  !***********************************************************!
  
  SUBROUTINE write_amplitude3D(bk, filename, silofile, groupname)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)                :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    LOGICAL, OPTIONAL                      :: silofile
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname

    REAL(dp), ALLOCATABLE                  :: probk3D(:, :, :)
    INTEGER                                :: dims(3)
    
    !--------------------------------------------------!
    
    ALLOCATE(probk3D(1:kmesh%maxkpts,1:kmesh%maxthetapts,&
         1:kmesh%maxphipts))
    
    CALL get_amplitude3D(bk,probk3D)
    
    dims = (/ kmesh%maxkpts, kmesh%maxthetapts, kmesh%maxphipts /)
    
    IF(silofile) THEN
       CALL create_silo_file(filename)
       CALL write_silo_rectilinear_mesh3D(filename,dims, &
            kmesh%kpts, kmesh%thetapts, kmesh%phipts)
       CALL write_silo_function3D(filename, dims, probk3D)
    ENDIF
    
    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(probk3D, &
            (/kmesh%maxkpts*kmesh%maxthetapts* &
            kmesh%maxphipts/)),3,&
            (/kmesh%maxkpts,kmesh%maxthetapts, &
            kmesh%maxphipts/), &
            filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(probk3D, &
            (/kmesh%maxkpts*kmesh%maxthetapts*kmesh%maxphipts/)),3,&
            (/kmesh%maxkpts,kmesh%maxthetapts,kmesh%maxphipts/), &
            filename)
    ENDIF

    DEALLOCATE(probk3D)
    
  END SUBROUTINE write_amplitude3D
  
  !***********************************************************!
  
  SUBROUTINE get_momwave(bk, cbk_rad)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    COMPLEX(dp), INTENT(OUT)   :: cbk_rad(:)

    COMPLEX(dp)                :: csuma
    INTEGER                    :: ik, itheta, iphi

    !---------------------------------------!
    
    cbk_rad = ZERO
    
    DO ik = 1, kmesh%numkpts
       csuma = ZERO
       DO iphi = 1, kmesh%numphipts
          DO itheta = 1, kmesh%numthetapts
             csuma = csuma + &
                  bk(ik,itheta,iphi) * &
                  kmesh%wtheta(itheta) * &
                  kmesh%wphi(iphi)
          ENDDO
       ENDDO
       cbk_rad(ik) = csuma
    ENDDO
    
    
  END SUBROUTINE get_momwave
  
  !***********************************************************!
  
  SUBROUTINE write_momwave(bk, filename, groupname)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    COMPLEX(dp), ALLOCATABLE       :: probk1D(:)
    INTEGER                        :: numkpts, numthetapts
    INTEGER                        :: numphipts, dims(3)
    CHARACTER(LEN=100)             :: name
    INTEGER                        :: ik, itheta, iphi
    
    !--------------------------------------------------!
    
    !ALLOCATE(probk1D(1:kmesh%numkpts))
    
    !CALL get_momwave(bk,probk1D)
    
    name = TRIM(filename) // '_psik.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
!!$    DO ik = 1, numkpts
!!$       WRITE(55,*) k_ax(ik), REAL(probk1D(ik),dp),AIMAG(probk1D(ik))
!!$    ENDDO
    
    DO ik = 1, kmesh%numkpts
       WRITE(55,*) kmesh%kpts(ik), REAL(bk(ik,1,1),dp),AIMAG(bk(ik,1,1))
    ENDDO
    
    CLOSE(55)
    
    !DEALLOCATE(probk1D)
    
  END SUBROUTINE write_momwave

  !***********************************************************!
  
  SUBROUTINE get_radial_amp_lm(bk, b_lm, sph_harmonics)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    REAL(dp), INTENT(OUT)          :: b_lm(:, -kmesh%lmax:, 0:)
    COMPLEX(dp), INTENT(IN)       :: sph_harmonics(:, :, &
         -kmesh%lmax:, 0:)
    
    COMPLEX(dp), ALLOCATABLE       :: cb_lm(:, :, :)
    INTEGER                        :: ik, il, im
    
    !----------------------------------------------!
    
    ALLOCATE(cb_lm(1:kmesh%numkpts,-kmesh%lmax:kmesh%lmax, &
         0:kmesh%lmax))
    
    b_lm  = ZERO
    cb_lm = ZERO
    
    DO ik = 1, kmesh%numkpts
       CALL make_sht(bk(ik,:,:),sph_harmonics,kmesh%lmax, &
            kmesh%wtheta,kmesh%wphi, cb_lm(ik,:,:))
    ENDDO
    
    DO ik = 1, kmesh%numkpts
       DO il = 0, kmesh%lmax
          DO im = -il, il
             b_lm(ik,im,il) = REAL(CONJG(cb_lm(ik,im,il)) * &
                  cb_lm(ik,im,il),dp)
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(cb_lm)
    
  END SUBROUTINE get_radial_amp_lm
  
  !***********************************************************!
  
  SUBROUTINE write_radial_amp_lm(bk, filename, sph_harmonics)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    COMPLEX(dp), INTENT(IN)       :: sph_harmonics(:, :, &
         -kmesh%lmax:, 0:)
    
    REAL(dp), ALLOCATABLE          :: b_lm(:, :, :)
    CHARACTER(LEN=3)               :: clmax, cmmax
    CHARACTER(LEN=100)             :: blmformat, name
    INTEGER                        :: ik, il, im
    
    !----------------------------------------------!

    ALLOCATE(b_lm(1:kmesh%numkpts,-kmesh%lmax: &
         kmesh%lmax,0:kmesh%lmax))
    
    CALL get_radial_amp_lm(bk,b_lm, sph_harmonics)
    
    WRITE(cmmax,'(I3.3)') 2 * kmesh%lmax + 1
    
    blmformat = '(1X,F21.16,'// cmmax //'(1X,F21.16))'
    
    DO il = 0, kmesh%lmax
       
       WRITE(clmax,'(I3.3)') il
       name = TRIM(filename) // '_l_' // clmax // '.dat'
       OPEN(UNIT=55,FORM='formatted',FILE=name)
       
       DO ik = 1, kmesh%numkpts
          WRITE(55,blmformat) kmesh%kpts(ik), (b_lm(ik,im,il), im=-il,il)
       ENDDO
       
       CLOSE(55)
       
    ENDDO
    
    DEALLOCATE(b_lm)
    
  END SUBROUTINE write_radial_amp_lm
  
  !***********************************************************!
  !***********************************************************!
  
END MODULE observables
