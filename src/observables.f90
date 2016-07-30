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
  USE omp_lib
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC        :: get_angular_resolution
  PUBLIC        :: get_amplitude
  PUBLIC        :: get_polar_amplitude
  PUBLIC        :: get_mes
  PUBLIC        :: get_pes
  PUBLIC        :: get_momwave
  PUBLIC        :: get_mes_angbasis
  PUBLIC        :: write_amplitude
  PUBLIC        :: write_polar_amplitude
  PUBLIC        :: write_mes
  PUBLIC        :: write_pes
  PUBLIC        :: write_momwave
  PUBLIC        :: write_mes_angbasis
  
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

    mmin = -mmax
   
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

    REAL(dp)                   :: dphi, suma
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
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik) &
    !$OMP& PRIVATE(itheta,iphi,suma)
    
    !$OMP DO
    DO ik = 1, numkpts
       suma = 0.0_dp
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             suma = suma + &
                  REAL(CONJG(bk(ik,itheta,iphi)) * bk(ik,itheta,iphi),dp) * &
                  gauss_th_weights(itheta)
          ENDDO
       ENDDO
       bk_rad(ik) = suma
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
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
                  REAL(CONJG(bk(ik,itheta,iphi) * bk(ik,itheta,iphi)),dp)
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
                  REAL(CONJG(bk(ik,itheta,iphi) * bk(ik,itheta,iphi)),dp)
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
  
  SUBROUTINE get_momwave(bk, cbk_rad)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)    :: bk(:, :, :)
    COMPLEX(dp), INTENT(OUT)   :: cbk_rad(:)

    REAL(dp)                   :: dphi
    COMPLEX(dp)                :: csuma
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
    
    cbk_rad = ZERO
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik) &
    !$OMP& PRIVATE(itheta,iphi,csuma)
    
    !$OMP DO
    DO ik = 1, numkpts
       csuma = ZERO
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             csuma = csuma + &
                  bk(ik,itheta,iphi) * gauss_th_weights(itheta)
          ENDDO
       ENDDO
       cbk_rad(ik) = csuma
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    IF(numphipts.EQ.1) THEN
       cbk_rad = cbk_rad * twopi
    ELSE
       cbk_rad = cbk_rad * dphi
    ENDIF
        
    
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
    
    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)

    !ALLOCATE(probk1D(1:numkpts))
    
    !CALL get_momwave(bk,probk1D)
    
    name = TRIM(filename) // '_psik.dat'
    OPEN(UNIT=55,FORM='formatted',FILE=name)
    
!!$    DO ik = 1, numkpts
!!$       WRITE(55,*) k_ax(ik), REAL(probk1D(ik),dp),AIMAG(probk1D(ik))
!!$    ENDDO

    DO ik = 1, numkpts
       WRITE(55,*) k_ax(ik), REAL(bk(ik,1,1),dp),AIMAG(bk(ik,1,1))
    ENDDO
    
    CLOSE(55)
    
    !DEALLOCATE(probk1D)
    
  END SUBROUTINE write_momwave

  !***********************************************************!
  
  SUBROUTINE get_mes_angbasis(bk, b_lm, lmax)
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    INTEGER, INTENT(IN)            :: lmax
    REAL(dp), INTENT(OUT)          :: b_lm(:, -lmax:, 0:)
    
    COMPLEX(dp), ALLOCATABLE       :: cb_lm(:, :, :)
    INTEGER                        :: numkpts, numthetapts
    INTEGER                        :: numphipts
    INTEGER                        :: dims(3)
    INTEGER                        :: ik, il, im
    
    !----------------------------------------------!
    
    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)    
        
    ALLOCATE(cb_lm(1:numkpts,-lmax:lmax,0:lmax))

    DO ik = 1, numkpts
       CALL make_sht(bk(ik,:,:),lmax,gauss_th_weights,&
            cb_lm(ik,:,:))
    ENDDO
    
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ik,il,im)
    
    !$OMP DO COLLAPSE(3)
    DO il = 0, lmax
       DO im = -il, il
          DO ik = 1, numkpts
             b_lm(ik,im,il) = REAL(CONJG(cb_lm(ik,im,il)) * &
                  cb_lm(ik,im,il),dp)
          ENDDO
       ENDDO
    ENDDO
    !$OMP END DO NOWAIT
    !$OMP END PARALLEL
    
    DEALLOCATE(cb_lm)
    
  END SUBROUTINE get_mes_angbasis
  
  !***********************************************************!
  
  SUBROUTINE write_mes_angbasis(bk, filename, lmax)
    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)        :: bk(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: filename
    INTEGER, INTENT(IN)            :: lmax

    REAL(dp), ALLOCATABLE          :: b_lm(:, :, :)
    INTEGER                        :: numkpts, numthetapts
    INTEGER                        :: numphipts
    INTEGER                        :: dims(3)
    CHARACTER(LEN=3)               :: clmax, cmmax
    CHARACTER(LEN=100)             :: blmformat, name
    INTEGER                        :: ik, il, im
    
    !----------------------------------------------!

    dims = SHAPE(bk)
    numkpts     = dims(1)
    numthetapts = dims(2)
    numphipts   = dims(3)
    
    ALLOCATE(b_lm(1:numkpts,-lmax:lmax,0:lmax))
    
    CALL get_mes_angbasis(bk,b_lm,lmax)
    
    WRITE(cmmax,'(I3.3)') 2 * lmax + 1
    
    blmformat = '(1X,F21.16,'// cmmax //'(1X,F21.16))'
    
    DO il = 0, lmax
       
       WRITE(clmax,'(I3.3)') il
       name = TRIM(filename) // '_l_' // clmax // '.dat'
       OPEN(UNIT=55,FORM='formatted',FILE=name)
       
       DO ik = 1, numkpts
          WRITE(55,blmformat) k_ax(ik), (b_lm(ik,im,il), im=-il,il)
       ENDDO
       
       CLOSE(55)
       
    ENDDO
    
    DEALLOCATE(b_lm)
    
  END SUBROUTINE write_mes_angbasis
  
  !***********************************************************!
  !***********************************************************!
  
END MODULE observables
