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
  
  USE constants
  USE gaussleg
  USE sht
  USE flux
  USE io_pop
  
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC        :: get_angular_resolution
  PUBLIC        :: get_amplitude
  PUBLIC        :: get_polar_amplitude
  PUBLIC        :: get_mes
  PUBLIC        :: get_pes
  PUBLIC        :: write_amplitude
  PUBLIC        :: write_polar_amplitude
  PUBLIC        :: write_mes
  PUBLIC        :: write_pes
  
  
CONTAINS
  
  SUBROUTINE get_angular_resolution(b_lm, dPomega, lmax, mmax)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: b_lm(:, -mmax:, 0:)
    COMPLEX(dp), INTENT(OUT)      :: dPomega(:, :, :)
    INTEGER, INTENT(IN)           :: lmax,mmax
    
    INTEGER                       :: mmin, dims(3)
    INTEGER                       :: numkpts, numthetapts
    INTEGER                       :: numphipts
    INTEGER                       :: ik,il,im,itheta,iphi
    
    !---------------------------------------------------------!
    
    dims = SHAPE(dPomega)
    
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)
    
    dPomega = ZERO
    
    DO ik = 1, numkpts
       DO il = 0, lmax
          DO im = mmin, mmax
             IF(im.LT.0) THEN
                DO iphi = 1, numphipts
                   DO itheta = 1, numthetapts
                      dPomega(ik,itheta,iphi) = b_lm(ik,im,il) * &
                           1.0_dp / k_ax(ik) * (-1.0_dp)**ABS(im) * &
                           normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ENDDO
                ENDDO
             ELSE
                DO iphi = 1, numphipts
                   DO itheta = 1, numthetapts
                      dPomega(ik,itheta,iphi) = b_lm(ik,im,il) * &
                           1.0_dp / k_ax(ik) * &
                           normfact(im,il) * legenpl(itheta-1,im,il) * &
                           EXP(ZIMAGONE * im * phi_ax(iphi))
                   ENDDO
                ENDDO
             ENDIF
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE get_angular_resolution

  !***********************************************************!

  SUBROUTINE get_mes(bk, bk_rad)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_rad(:)

    REAL(dp)                   :: dphi
    INTEGER                    :: numkpts, numthetapts
    INTEGER                    :: numphipts, dims(3)
    INTEGER                    :: ik, itheta, iphi

    !---------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)
    
    IF(numphipts.EQ.1) THEN
       dphi = 1.0_dp
    ELSE
       dphi = phi_ax(2) - phi_ax(1)
    ENDIF
    
    bk_rad = 0.0_dp
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             bk_rad(ik) = bk_rad(ik) + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * bk(ik,itheta,iphi)) * &
                  gauss_th_weights(itheta)
          ENDDO
       ENDDO
    ENDDO


    IF(numphipts.EQ.1) THEN
       bk_rad = bk_rad * twopi
    ELSE
       bk_rad = bk_rad * dphi
    ENDIF
    
  END SUBROUTINE get_mes
  
  !***********************************************************!
  
  SUBROUTINE write_mes(bk, filename, groupname)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    REAL(dp), ALLOCATABLE          :: probk1D(:)
    INTEGER                        :: numkpts, numthetapts
    INTEGER                        :: numphipts, dims(3)
    CHARACTER(LEN=100)             :: name
    INTEGER                        :: ik, itheta, iphi
    !--------------------------------------------------!
    
    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    ALLOCATE(probk1D(1:numkpts))

    CALL get_mes(bk,probk1D)

    ! IF(PRESENT(groupname)) THEN
    !    CALL write_wave(probk1D,1,(/numkpts/),filename,&
    !         groupname,groupname)
    ! ELSE
    !    CALL write_wave(probk1D,1,(/numkpts/),filename)
    ! ENDIF
    
    name = TRIM(filename) // '.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
    DO ik = 1, numkpts
       WRITE(55,*) k_ax(ik), probk1D(ik)
    ENDDO

    CLOSE(55)
    
    DEALLOCATE(probk1D)

  END SUBROUTINE write_mes

  !***********************************************************!

  SUBROUTINE get_pes(bk, bk_pes)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)      :: bk(:, :, :)
    REAL(dp), INTENT(OUT)        :: bk_pes(:)

    INTEGER                      :: numkpts, numthetapts
    INTEGER                      :: numphipts, dims(3)
    INTEGER                      :: ik, itheta, iphi
    !-------------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    CALL get_mes(bk, bk_pes)

    DO ik = 1, numkpts
      bk_pes(ik) = bk_pes(ik) * k_ax(ik)
    ENDDO

  END SUBROUTINE get_pes

  !***********************************************************!

  SUBROUTINE write_pes(bk, filename, groupname)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname

    REAL(dp), ALLOCATABLE          :: prob1D(:)
    INTEGER                        :: numkpts, numthetapts
    INTEGER                        :: numphipts, dims(3)
    CHARACTER(LEN=100)             :: name
    INTEGER                        :: ik, itheta, iphi
    !--------------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    ALLOCATE(prob1D(1:numkpts))

    CALL get_pes(bk,prob1D)

    ! IF(PRESENT(groupname)) THEN
    !    CALL write_wave(probk1D,1,(/numkpts/),filename,&
    !         groupname,groupname)
    ! ELSE
    !    CALL write_wave(probk1D,1,(/numkpts/),filename)
    ! ENDIF

    name = TRIM(filename) // '.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)

    DO ik = 1, numkpts
      WRITE(55,*) k_ax(ik)**2*0.5_dp, prob1D(ik)
    ENDDO

    CLOSE(55)
    
    DEALLOCATE(prob1D)

  END SUBROUTINE write_pes

  !***********************************************************!

  SUBROUTINE get_polar_amplitude(bk, bk_polar)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_polar(:, :)

    REAL(dp)                   :: dphi
    INTEGER                    :: numkpts, numthetapts
    INTEGER                    :: numphipts, dims(3)
    INTEGER                    :: ik, itheta, iphi

    !----------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    bk_polar = 0.0_dp
    IF(numphipts.EQ.1) THEN
       dphi = 1.0_dp
    ELSE
       dphi = phi_ax(2) - phi_ax(1)
    ENDIF

    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             bk_polar(ik,itheta) = bk_polar(ik,itheta) + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * bk(ik,itheta,iphi))
          ENDDO
       ENDDO
    ENDDO

    
    IF(numphipts.EQ.1) THEN
       bk_polar = bk_polar * twopi
    ELSE
       bk_polar = bk_polar * dphi
    ENDIF
    
  END SUBROUTINE get_polar_amplitude
  
  !***********************************************************!

  SUBROUTINE write_polar_amplitude(bk, filename, groupname)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)                :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname

    INTEGER                                :: numkpts, numthetapts
    INTEGER                                :: numphipts, dims(3)
    REAL(dp), ALLOCATABLE                  :: probk2D(:, :)

    !--------------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    ALLOCATE(probk2D(1:numkpts,1:numthetapts))

    CALL get_polar_amplitude(bk,probk2D)

    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(probk2D,(/numkpts*numthetapts/)),2,&
            (/numkpts,numthetapts/),filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(probk2D,(/numkpts*numthetapts/)),2,&
            (/numkpts,numthetapts/),filename)
    ENDIF

    DEALLOCATE(probk2D)

  END SUBROUTINE write_polar_amplitude

    !***********************************************************!

  SUBROUTINE get_amplitude(bk, bk_real)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    REAL(dp), INTENT(OUT)      :: bk_real(:, :, :)

    REAL(dp)                   :: dphi
    INTEGER                    :: numkpts, numthetapts
    INTEGER                    :: numphipts, dims(3)
    INTEGER                    :: ik, itheta, iphi

    !----------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    bk_real = 0.0_dp

    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ik = 1, numkpts
             bk_real(ik,itheta,iphi) = bk_real(ik,itheta,iphi) + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * bk(ik,itheta,iphi))
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE get_amplitude

  !***********************************************************!
  
  SUBROUTINE write_amplitude(bk, filename, groupname)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)                :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)           :: filename
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: groupname
    
    INTEGER                                :: numkpts
    INTEGER                                :: numthetapts
    INTEGER                                :: numphipts
    INTEGER                                :: dims(3)
    REAL(dp), ALLOCATABLE                  :: probk3D(:, :, :)
    
    !--------------------------------------------------!
    
    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)
    
    ALLOCATE(probk3D(1:numkpts,1:numthetapts,1:numphipts))
    
    CALL get_amplitude(bk,probk3D)
    
    IF(PRESENT(groupname)) THEN
       CALL write_wave(RESHAPE(probk3D,(/numkpts*numthetapts*numphipts/)),3,&
            (/numkpts,numthetapts,numphipts/),filename,groupname,groupname)
    ELSE
       CALL write_wave(RESHAPE(probk3D,(/numkpts*numthetapts*numphipts/)),3,&
            (/numkpts,numthetapts,numphipts/),filename)
    ENDIF
    
    DEALLOCATE(probk3D)
    
  END SUBROUTINE write_amplitude
  
  !***********************************************************!
  !***********************************************************!

END MODULE observables
