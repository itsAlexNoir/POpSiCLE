!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: interp
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Interpolating routines.
!> \details This module contains all the subprograms needed for
!> an interpolation of the wavefunction, including routines for
!> multidimensions grid. Many of the routines in this module have
!> been extract form the numerical recipes books.
!
!------------------------------------------------------------------------------

MODULE interp
  
  USE constants
  USE sheppack
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC      :: create_interpolant
  PUBLIC      :: interpolate
  PUBLIC      :: destroy_interpolant
  
  
  INTERFACE create_interpolant
     MODULE PROCEDURE create_interpolant2D
     MODULE PROCEDURE dcreate_interpolant3D
     MODULE PROCEDURE zcreate_interpolant3D
  END INTERFACE create_interpolant
  
  INTERFACE interpolate
     MODULE PROCEDURE interpolate2D
     MODULE PROCEDURE dinterpolate3D
     MODULE PROCEDURE zinterpolate3D
  END INTERFACE interpolate
  
  INTERFACE destroy_interpolant
     MODULE PROCEDURE destroy_interpolant2D
     MODULE PROCEDURE ddestroy_interpolant3D
     MODULE PROCEDURE zdestroy_interpolant3D
  END INTERFACE destroy_interpolant
  
  
  ! Variables
  REAL(dp), ALLOCATABLE        :: func_re(:)
  REAL(dp), ALLOCATABLE        :: func_im(:)

  ! QSHEP2D control parameters
  INTEGER                      :: NQ, NR, NW
  REAL(dp)                     :: DX_RE, DY_RE, RMAX_RE
  REAL(dp)                     :: XMIN_RE, YMIN_RE
  REAL(dp), ALLOCATABLE        :: XYZMIN_RE(:), XYZDEL_RE(:)
  INTEGER, ALLOCATABLE         :: LNEXT_RE(:)
  INTEGER, ALLOCATABLE         :: LCELL_RE(:, :)
  INTEGER, ALLOCATABLE         :: LCELL3_RE(:, :, :)
  REAL(dp), ALLOCATABLE        :: RSQ_RE(:), RW_RE(:)
  REAL(dp), ALLOCATABLE        :: AQ_RE(:, :), AL_RE(:, :)
  
  REAL(dp)                     :: DX_IM, DY_IM, RMAX_IM
  REAL(dp)                     :: XMIN_IM, YMIN_IM
  REAL(dp), ALLOCATABLE        :: XYZMIN_IM(:), XYZDEL_IM(:)
  INTEGER, ALLOCATABLE         :: LNEXT_IM(:)
  INTEGER, ALLOCATABLE         :: LCELL_IM(:, :)
  INTEGER, ALLOCATABLE         :: LCELL3_IM(:, :, :)
  REAL(dp), ALLOCATABLE        :: RSQ_IM(:), RW_IM(:)
  REAL(dp), ALLOCATABLE        :: AQ_IM(:, :), AL_IM(:, :)
  
  
CONTAINS
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE create_interpolant
  !
  !> \brief Create the interpolant for scattered interpolantion
  !> \details Create and interpolant for scattered interpolation.
  !> The routine can choose between different orders.
  !
  !
  !> \param[in] xx Array of tabulated values
  !> \param[in] x Given value to be search
  !> \return locate Array index of xx in which x is centered.
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE create_interpolant2D(numpts, x1, x2, func, method )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    COMPLEX(dp), INTENT(IN)         :: func(:)
    CHARACTER(LEN=*), INTENT(IN)    :: method
    INTEGER                         :: ierror
    
    
    IF(method.EQ.'quadratic') THEN
       
       ! Assign control parameters for interpolation
       NR = 3
       NQ = 13
       NW = 19
       
       ! Allocate arrays for the routine
       ! The real ones
       ALLOCATE(LNEXT_RE(numpts))
       ALLOCATE(LCELL_RE(NR, NR))
       ALLOCATE(RSQ_RE(numpts))
       ALLOCATE(RW_RE(numpts))
       ALLOCATE(AQ_RE(5,numpts))
       ALLOCATE(AL_RE(2,numpts))
       ! The imaginary ones
       ALLOCATE(LNEXT_IM(numpts))
       ALLOCATE(LCELL_IM(NR, NR))
       ALLOCATE(RSQ_IM(numpts))
       ALLOCATE(RW_IM(numpts))
       ALLOCATE(AQ_IM(5,numpts))
       ALLOCATE(AL_IM(2,numpts))
       
       ALLOCATE(func_re(1:numpts))
       ALLOCATE(func_im(1:numpts))
       
       func_re = REAL(func,dp)
       func_im = AIMAG(func)

       !WRITE(*,*) 'Creating the interpolant...'
       CALL QSHEP2(numpts, x1, x2, func_re, NQ, NW, NR, LCELL_RE, &
            LNEXT_RE, XMIN_RE, YMIN_RE, DX_RE, DY_RE, &
            RMAX_RE, RSQ_RE, AQ_RE, ierror)
       IF (ierror /= 0) THEN
          WRITE(*,*) 'Error calculating the (real) interpolant at QSHEP2.'
          WRITE(*,*) 'ierror: ',ierror
          STOP
       ENDIF
       
       CALL QSHEP2(numpts, x1, x2, func_im, NQ, NW, NR, LCELL_IM, &
            LNEXT_IM, XMIN_IM, YMIN_IM, DX_IM, DY_IM, &
            RMAX_IM, RSQ_IM, AQ_IM, ierror)
       IF (ierror /= 0) THEN
          WRITE(*,*) 'Error calculating the (imag) interpolant at QSHEP2.'
          WRITE(*,*) 'ierror: ',ierror
          STOP
       ENDIF
       
    ELSEIF(method.eq.'cubic') THEN
       
       
    ELSEIF(method.eq.'linear') THEN
       
    ELSE
       
       WRITE(*,*) 'The input method has been implemented.'
       STOP
       
    END IF
    
  END SUBROUTINE create_interpolant2D
  
  !--------------------------------------------------------------!
  
  SUBROUTINE dcreate_interpolant3D(numpts, x1, x2, x3, func, method )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    REAL(dp), INTENT(IN)            :: x3(:)
    REAL(dp), INTENT(IN)            :: func(:)
    CHARACTER(LEN=*), INTENT(IN)    :: method
    INTEGER                         :: ierror
    
    IF(method.EQ.'quadratic') THEN
       
       ! Assign control parameters for interpolation
       NR = 3
       NQ = 17
       NW = 32
       
       ! Allocate arrays for the routine
       ! The real ones
       ALLOCATE(LNEXT_RE(numpts))
       ALLOCATE(LCELL3_RE(NR, NR, NR))
       ALLOCATE(RSQ_RE(numpts))
       ALLOCATE(RW_RE(numpts))
       ALLOCATE(XYZMIN_RE(3))
       ALLOCATE(XYZDEL_RE(3))
       ALLOCATE(AQ_RE(9,numpts))
       ALLOCATE(AL_RE(2,numpts))
       
       !WRITE(*,*) 'Creating the interpolant...'
       CALL QSHEP3(numpts, x1, x2, x3, func, NQ, NW, NR, LCELL3_RE, &
            LNEXT_RE, XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, &
            AQ_RE, ierror)
       IF (ierror /= 0) THEN
          WRITE(*,*) 'Error calculating the (real) interpolant at QSHEP3.'
          WRITE(*,*) 'ierror: ',ierror
          STOP
       ENDIF
       
    ELSEIF(method.eq.'linear') THEN
       
    ELSE
       
       WRITE(*,*) 'The input method has been implemented.'
       STOP
       
    END IF
  END SUBROUTINE dcreate_interpolant3D
  
  !-------------------------------------------------!
  
  SUBROUTINE zcreate_interpolant3D(numpts, x1, x2, x3, func, method )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    REAL(dp), INTENT(IN)            :: x3(:)
    COMPLEX(dp), INTENT(IN)         :: func(:)
    CHARACTER(LEN=*), INTENT(IN)    :: method
    INTEGER                         :: ierror
    
    
    IF(method.EQ.'quadratic') THEN
       
       ! Assign control parameters for interpolation
       NR = 3
       NQ = 17
       NW = 32
       
       ! Allocate arrays for the routine
       ! The real ones
       ALLOCATE(LNEXT_RE(numpts))
       ALLOCATE(LCELL3_RE(NR, NR, NR))
       ALLOCATE(RSQ_RE(numpts))
       ALLOCATE(RW_RE(numpts))
       ALLOCATE(XYZMIN_RE(3))
       ALLOCATE(XYZDEL_RE(3))
       ALLOCATE(AQ_RE(9,numpts))
       ALLOCATE(AL_RE(2,numpts))
       ! The imaginary ones
       ALLOCATE(LNEXT_IM(numpts))
       ALLOCATE(LCELL3_IM(NR, NR, NR))
       ALLOCATE(RSQ_IM(numpts))
       ALLOCATE(RW_IM(numpts))
       ALLOCATE(XYZMIN_IM(3))
       ALLOCATE(XYZDEL_IM(3))
       ALLOCATE(AQ_IM(9,numpts))
       ALLOCATE(AL_IM(2,numpts))
       
       ALLOCATE(func_re(1:numpts))
       ALLOCATE(func_im(1:numpts))
       
       func_re = REAL(func,dp)
       func_im = AIMAG(func)
       
       !WRITE(*,*) 'Creating the interpolant...'
       CALL QSHEP3(numpts, x1, x2, x3, func_re, NQ, NW, NR, LCELL3_RE, &
            LNEXT_RE, XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, &
            AQ_RE, ierror)
       IF (ierror /= 0) THEN
          WRITE(*,*) 'Error calculating the (real) interpolant at QSHEP3.'
          WRITE(*,*) 'ierror: ',ierror
          STOP
       ENDIF
       
       CALL QSHEP3(numpts, x1, x2, x3, func_im, NQ, NW, NR, LCELL3_IM, &
            LNEXT_IM, XYZMIN_IM, XYZDEL_IM, RMAX_IM, RSQ_IM, &
            AQ_IM, ierror)
       IF (ierror /= 0) THEN
          WRITE(*,*) 'Error calculating the (imag) interpolant at QSHEP3.'
          WRITE(*,*) 'ierror: ',ierror
          STOP
       ENDIF
       
       
    ELSEIF(method.eq.'linear') THEN
       
    ELSE
       
       WRITE(*,*) 'The input method has been implemented.'
       STOP
       
    END IF
    
  END SUBROUTINE zcreate_interpolant3D
  
  !--------------------------------------------------------------!
  
  SUBROUTINE interpolate2D(numpts, y1, y2, x1, x2, func, method, &
       interp_val, interp_val_dx, interp_val_dy)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                  :: numpts
    REAL(dp), INTENT(IN)                 :: y1
    REAL(dp), INTENT(IN)                 :: y2
    REAL(dp), INTENT(IN)                 :: x1(:)
    REAL(dp), INTENT(IN)                 :: x2(:)
    COMPLEX(dp), INTENT(IN)              :: func(:)
    CHARACTER(LEN=*), INTENT(IN)         :: method
    COMPLEX(dp), INTENT(OUT)             :: interp_val
    COMPLEX(dp), INTENT(OUT), OPTIONAL   :: interp_val_dx
    COMPLEX(dp), INTENT(OUT), OPTIONAL   :: interp_val_dy
    
    REAL(dp)                     :: interp_val_re
    REAL(dp)                     :: interp_val_im
    REAL(dp)                     :: interp_val_dx_re
    REAL(dp)                     :: interp_val_dx_im
    REAL(dp)                     :: interp_val_dy_re
    REAL(dp)                     :: interp_val_dy_im  
    INTEGER                      :: ierror
    
    !-----------------------------------------------!
    
    IF (PRESENT(interp_val_dx).AND.PRESENT(interp_val_dy) ) THEN
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          CALL QS2GRD(y1, y2, numpts, x1, x2, func_re, NR, LCELL_RE, LNEXT_RE,&
               XMIN_RE, YMIN_RE, DX_RE, DY_RE, RMAX_RE, RSQ_RE, AQ_RE, interp_val_re, &
               interp_val_dx_re, interp_val_dy_re, ierror)
          IF (ierror /= 0) THEN
             WRITE (*, *) 'QS2GRD - ERROR!'
             WRITE (*, *) 'Error in QS2GRD (REAL), IER = ', ierror
             STOP
          END IF
          
          
          ! Interpolate the imaginary part
          CALL QS2GRD(y1, y2, numpts, x1, x2, func_im, NR, LCELL_IM, LNEXT_IM,&
               XMIN_IM, YMIN_IM, DX_IM, DY_IM, RMAX_IM, RSQ_IM, AQ_IM, interp_val_im, &
               interp_val_dx_im, interp_val_dy_im, ierror)
          IF (ierror /= 0) THEN
             WRITE (*, *) 'QS2GRD - ERROR!'
             WRITE (*, *) 'Error in QS2GRD (IMAG), IER = ', ierror
             STOP
          END IF
          
          !---------------------------------!
          
       ENDIF
       
       interp_val    = CMPLX(interp_val_re, interp_val_im, dp)
       interp_val_dx = CMPLX(interp_val_dx_re, interp_val_dx_im, dp)
       interp_val_dy = CMPLX(interp_val_dy_re, interp_val_dy_im, dp)
       
    
    ELSEIF (.NOT.PRESENT(interp_val_dx).AND. .NOT.PRESENT(interp_val_dy) ) THEN
       
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          interp_val_re = QS2VAL(y1, y2, numpts, x1, x2, func_re, NR, LCELL3_RE, LNEXT_RE, &
               XMIN_RE, YMIN_RE, DX_RE, DY_RE, RMAX_RE, RSQ_RE, AQ_RE )
          IF (interp_val_re .EQ. 0) THEN
             WRITE (*, *) 'Q2VAL - ERROR!'
             WRITE (*, *) 'Error in QS2VAL (REAL), IER = ', interp_val_re
             STOP
          END IF
          
          
          ! Interpolate the imaginary part
          interp_val_im = QS2VAL(y1, y2, numpts, x1, x2, func_im, NR, LCELL3_IM, LNEXT_IM, &
               XMIN_IM, YMIN_IM, DX_IM, DY_IM, RMAX_IM, RSQ_IM, AQ_IM )
          IF (interp_val_im .EQ. 0) THEN
             WRITE (*, *) 'Q2VAL - ERROR!'
             WRITE (*, *) 'Error in Q23VAL (IMAG), IER = ', interp_val_im
             STOP
          END IF
          
          !---------------------------------!
       ENDIF
       
       interp_val    = CMPLX(interp_val_re, interp_val_im, dp)

    ELSE
       WRITE(*,*) 'Missing derivatives arrays!'
       RETURN
    ENDIF
    
    
  END SUBROUTINE interpolate2D

  !-----------------------------------------------------------!

    
  SUBROUTINE dinterpolate3D(numpts, y1, y2, y3, x1, x2, x3, func, method,&
       interp_val, interp_val_dx, interp_val_dy, interp_val_dz )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                 :: numpts
    REAL(dp), INTENT(IN)                :: y1
    REAL(dp), INTENT(IN)                :: y2
    REAL(dp), INTENT(IN)                :: y3  
    REAL(dp), INTENT(IN)                :: x1(:)
    REAL(dp), INTENT(IN)                :: x2(:)
    REAL(dp), INTENT(IN)                :: x3(:)
    REAL(dp), INTENT(IN)                :: func(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    REAL(dp), INTENT(OUT)               :: interp_val
    REAL(dp), INTENT(OUT), OPTIONAL     :: interp_val_dx
    REAL(dp), INTENT(OUT), OPTIONAL     :: interp_val_dy
    REAL(dp), INTENT(OUT), OPTIONAL     :: interp_val_dz
    INTEGER                             :: ierror
    
    !-----------------------------------------------!
    
    
    IF (PRESENT(interp_val_dx).AND.PRESENT(interp_val_dy) &
         .AND.PRESENT(interp_val_dz) ) THEN
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          CALL QS3GRD(y1, y2, y3, numpts, x1, x2, x3, func, NR, LCELL3_RE, LNEXT_RE,&
               XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, AQ_RE, interp_val, &
               interp_val_dx, interp_val_dy, interp_val_dz, ierror)
          IF ((ierror /= 0) .AND. (ierror /=2)) THEN
             WRITE (*, *) 'Q3GRD - ERROR!'
             WRITE (*, *) 'Error in QS3GRD, IER = ', ierror
             STOP
          END IF
          
          !---------------------------------!
          
       ENDIF
       
    ELSEIF (.NOT.PRESENT(interp_val_dx).AND. .NOT.PRESENT(interp_val_dy) &
         .AND. .NOT.PRESENT(interp_val_dz) ) THEN
       
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          interp_val = QS3VAL(y1, y2, y3, numpts, x1, x2, x3, func, NR, LCELL3_RE, LNEXT_RE,&
               XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, AQ_RE )
!!$          IF (interp_val .EQ. 0) THEN
!!$             WRITE (*, *) 'Q3VAL - ERROR!'
!!$             WRITE (*, *) 'Error in QS3VAL (REAL), IER = ', interp_val
!!$             STOP
!!$          END IF
          
          !---------------------------------!
       ENDIF
       
    ELSE
       WRITE(*,*) 'Missing derivatives arrays!'
       RETURN
    ENDIF
    
  END SUBROUTINE dinterpolate3D
  
  !---------------------------------------------------------------------------!
  
  SUBROUTINE zinterpolate3D(numpts, y1, y2, y3, x1, x2, x3, func, method, &
       interp_val, interp_val_dx, interp_val_dy, interp_val_dz )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)                 :: numpts
    REAL(dp), INTENT(IN)                :: y1
    REAL(dp), INTENT(IN)                :: y2
    REAL(dp), INTENT(IN)                :: y3  
    REAL(dp), INTENT(IN)                :: x1(:)
    REAL(dp), INTENT(IN)                :: x2(:)
    REAL(dp), INTENT(IN)                :: x3(:)
    COMPLEX(dp), INTENT(IN)             :: func(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT)            :: interp_val
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: interp_val_dx
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: interp_val_dy
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: interp_val_dz
    
    REAL(dp)                     :: interp_val_re
    REAL(dp)                     :: interp_val_im
    REAL(dp)                     :: interp_val_dx_re
    REAL(dp)                     :: interp_val_dx_im
    REAL(dp)                     :: interp_val_dy_re
    REAL(dp)                     :: interp_val_dy_im  
    REAL(dp)                     :: interp_val_dz_re
    REAL(dp)                     :: interp_val_dz_im
    INTEGER                      :: ierror
    
    !-----------------------------------------------!
    
    
    IF (PRESENT(interp_val_dx).AND.PRESENT(interp_val_dy) &
         .AND.PRESENT(interp_val_dz) ) THEN
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          CALL QS3GRD(y1, y2, y3, numpts, x1, x2, x3, func_re, NR, LCELL3_RE, LNEXT_RE,&
               XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, AQ_RE, interp_val_re, &
               interp_val_dx_re, interp_val_dy_re, interp_val_dz_re, ierror)
          IF ((ierror /= 0).AND.(ierror /=2)) THEN
             WRITE (*, *) 'QSHEP3 - ERROR!'
             WRITE (*, *) 'Error in QS3GRD (REAL), IER = ', ierror
             !STOP
          END IF
          
          
          ! Interpolate the imaginary part
          CALL QS3GRD(y1, y2, y3, numpts, x1, x2, x3, func_im, NR, LCELL3_IM, LNEXT_IM,&
               XYZMIN_IM, XYZDEL_IM, RMAX_IM, RSQ_IM, AQ_IM, interp_val_im, &
               interp_val_dx_im, interp_val_dy_im, interp_val_dz_im, ierror)
          IF ((ierror /= 0).AND. (ierror /=2)) THEN
             WRITE (*, *) 'QSHEP3 - ERROR!'
             WRITE (*, *) 'Error in QS3GRD (IMAG), IER = ', ierror
             !STOP
          END IF
          
       !---------------------------------!
          
       ENDIF
       
       interp_val    = CMPLX(interp_val_re, interp_val_im)
       interp_val_dx = CMPLX(interp_val_dx_re, interp_val_dx_im)
       interp_val_dy = CMPLX(interp_val_dy_re, interp_val_dy_im)
       interp_val_dz = CMPLX(interp_val_dz_re, interp_val_dz_im)
       
    ELSEIF (.NOT.PRESENT(interp_val_dx).AND. .NOT.PRESENT(interp_val_dy) &
         .AND. .NOT.PRESENT(interp_val_dz) ) THEN
       
       IF (method .EQ. 'quadratic') THEN
          ! Interpolate the real part
          interp_val_re = QS3VAL(y1, y2, y3, numpts, x1, x2, x3, func_re, NR, LCELL3_RE, LNEXT_RE,&
               XYZMIN_RE, XYZDEL_RE, RMAX_RE, RSQ_RE, AQ_RE )
!!$          IF (interp_val_re .EQ. 0) THEN
!!$             WRITE (*, *) 'Q3VAL - ERROR!'
!!$             WRITE (*, *) 'Error in QS3VAL (REAL), IER = ', interp_val_re
!!$             STOP
!!$          END IF
          
          
          ! Interpolate the imaginary part
          interp_val_im = QS3VAL(y1, y2, y3, numpts, x1, x2, x3, func_im, NR, LCELL3_IM, LNEXT_IM,&
               XYZMIN_IM, XYZDEL_IM, RMAX_IM, RSQ_IM, AQ_IM )
!!$          IF (interp_val_im .EQ. 0) THEN
!!$             WRITE (*, *) 'Q3VAL - ERROR!'
!!$             WRITE (*, *) 'Error in QS3VAL (IMAG), IER = ', interp_val_im
!!$             STOP
!!$          END IF
          
          !---------------------------------!
       ENDIF
       
    ELSE
       WRITE(*,*) 'Missing derivatives arrays!'
       RETURN
    ENDIF
    
  END SUBROUTINE zinterpolate3D
  
  !-------------------------------------!
  
  SUBROUTINE destroy_interpolant2D(numpts, x1, x2, func )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    COMPLEX(dp), INTENT(IN)         :: func(:)
    
    
    ! Allocate arrays for the routine
    ! The real ones
    DEALLOCATE(LNEXT_RE)
    DEALLOCATE(LCELL_RE)
    DEALLOCATE(RSQ_RE)
    DEALLOCATE(RW_RE)
    DEALLOCATE(AQ_RE)
    DEALLOCATE(AL_RE)
    ! The imaginary ones
    DEALLOCATE(LNEXT_IM)
    DEALLOCATE(LCELL_IM)
    DEALLOCATE(RSQ_IM)
    DEALLOCATE(RW_IM)
    DEALLOCATE(AQ_IM)
    DEALLOCATE(AL_IM)

    DEALLOCATE(func_re,func_im)
    
  END SUBROUTINE destroy_interpolant2D

  !-------------------------------------!
  
  SUBROUTINE ddestroy_interpolant3D( numpts, x1, x2, x3, func )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    REAL(dp), INTENT(IN)            :: x3(:)
    REAL(dp), INTENT(IN)            :: func(:)
    
    ! Allocate arrays for the routine
    ! The real ones
    DEALLOCATE(LNEXT_RE)
    DEALLOCATE(LCELL3_RE)
    DEALLOCATE(RSQ_RE)
    DEALLOCATE(RW_RE)
    DEALLOCATE(XYZMIN_RE)
    DEALLOCATE(XYZDEL_RE)
    DEALLOCATE(AQ_RE)
    DEALLOCATE(AL_RE)
    
  END SUBROUTINE ddestroy_interpolant3D
  
  !--------------------------------------!
  SUBROUTINE zdestroy_interpolant3D( numpts, x1, x2, x3, func )
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)             :: numpts
    REAL(dp), INTENT(IN)            :: x1(:)
    REAL(dp), INTENT(IN)            :: x2(:)
    REAL(dp), INTENT(IN)            :: x3(:)
    COMPLEX(dp), INTENT(IN)         :: func(:)
    
    ! Allocate arrays for the routine
    ! The real ones
    DEALLOCATE(LNEXT_RE)
    DEALLOCATE(LCELL3_RE)
    DEALLOCATE(RSQ_RE)
    DEALLOCATE(RW_RE)
    DEALLOCATE(XYZMIN_RE)
    DEALLOCATE(XYZDEL_RE)
    DEALLOCATE(AQ_RE)
    DEALLOCATE(AL_RE)
    ! The imaginary ones
    DEALLOCATE(LNEXT_IM)
    DEALLOCATE(LCELL3_IM)
    DEALLOCATE(RSQ_IM)
    DEALLOCATE(RW_IM)
    DEALLOCATE(XYZMIN_IM)
    DEALLOCATE(XYZDEL_IM)
    DEALLOCATE(AQ_IM)
    DEALLOCATE(AL_IM)
    
    DEALLOCATE(func_re,func_im)
    
  END SUBROUTINE zdestroy_interpolant3D
  
  !------------------------------------------!
  !------------------------------------------!
  
END MODULE interp
