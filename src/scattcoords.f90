

! .d8888.  .o88b.  .d8b.  d888888b d888888b
! 88'  YP d8P  Y8 d8' `8b `~~88~~' `~~88~~'
! `8bo.   8P      88ooo88    88       88
!   `Y8b. 8b      88~~~88    88       88
! db   8D Y8b  d8 88   88    88       88
! `8888Y'  `Y88P' YP   YP    YP       YP
!
!
!  .o88b.  .d88b.   .d88b.  d8888b. d8888b. .d8888.
! d8P  Y8 .8P  Y8. .8P  Y8. 88  `8D 88  `8D 88'  YP
! 8P      88    88 88    88 88oobY' 88   88 `8bo.
! 8b      88    88 88    88 88`8b   88   88   `Y8b.
! Y8b  d8 `8b  d8' `8b  d8' 88 `88. 88  .8D db   8D
!  `Y88P'  `Y88P'   `Y88P'  88   YD Y8888D' `8888Y'


MODULE scattcoords

  USE constants
  USE scatt_interp

  IMPLICIT NONE

  PRIVATE

  PUBLIC                  :: dscatt_cartesian2spherical2D
  PUBLIC                  :: zscatt_cartesian2spherical2D
  PUBLIC                  :: dscatt_cartesian2spherical3D
  PUBLIC                  :: zscatt_cartesian2spherical3D
  PUBLIC                  :: dscatt_cylindrical2spherical2D
  PUBLIC                  :: zscatt_cylindrical2spherical2D

  INTEGER                 :: numpts, numrpts
  INTEGER                 :: numthetapts, numphipts

CONTAINS

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !               _          _           ___         _            _         _
  !   __ __ _ _ _| |_ ___ __(_)__ _ _ _ |_  )____ __| |_  ___ _ _(_)__ __ _| |
  !  / _/ _` | '_|  _/ -_|_-< / _` | ' \ / /(_-< '_ \ ' \/ -_) '_| / _/ _` | |
  !  \__\__,_|_|  \__\___/__/_\__,_|_||_/___/__/ .__/_||_\___|_| |_\__\__,_|_|
  !                                            |_|

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  !                          ___ ___
  !                         |_  )   \
  !                          / /| |) |
  !                         /___|___/
  !

  SUBROUTINE dscatt_cartesian2spherical2D(psi_cart,psi_sph, x_ax, y_ax, &
       r_ax, theta_ax, method, psi_sph_dr, psi_sph_dth,rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)                :: psi_cart(:, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dr(:, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dth(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(:, :)
    INTEGER                             :: ix, iy, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt, xpt, ypt
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    REAL(dp), ALLOCATABLE               :: psi_inp(:)

    !-------------------------------------------------------------------!

    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)

    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in2D(dims_out(1),dims_out(2)))

    psi_inp = ZERO
    psi_sph = ZERO
    i_am_in2D = 0

    ii = 0

    DO iy = 1, dims_in(2)
       DO ix = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 )

          x_inp(ii) = x_ax(ix)
          y_inp(ii) = y_ax(iy)
          psi_inp(ii) = psi_cart(ix, iy)
       ENDDO
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method),rank)

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta),rank )
                ENDIF
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank )
                ENDIF
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN
       ENDIF

    ELSE

       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method))

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta) )
                ENDIF
             ENDDO
          ENDDO


       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)

                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta) )
                ENDIF
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN
       ENDIF

    ENDIF

    DEALLOCATE(i_am_in2D)
    DEALLOCATE(x_inp,y_inp,psi_inp)

  END SUBROUTINE dscatt_cartesian2spherical2D

  !-------------------------------------------------------------------!

  SUBROUTINE zscatt_cartesian2spherical2D(psi_cart,psi_sph, x_ax, y_ax, &
       r_ax, theta_ax, method, psi_sph_dr, psi_sph_dth,rank)

    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(:, :)
    INTEGER                             :: ix, iy, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt, xpt, ypt
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)

    !-------------------------------------------------------------------!

    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)

    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in2D(dims_out(1),dims_out(2)))

    psi_inp = ZERO
    psi_sph = ZERO
    i_am_in2D = 0

    ii = 0

    DO iy = 1, dims_in(2)
       DO ix = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 )

          x_inp(ii) = x_ax(ix)
          y_inp(ii) = y_ax(iy)
          psi_inp(ii) = psi_cart(ix, iy)
       ENDDO
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method),rank)

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta),rank )
                ENDIF
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank )
                ENDIF
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN
       ENDIF

    ELSE

       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method))

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta) )
                ENDIF
             ENDDO
          ENDDO


       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN

          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)

                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   xpt = r_ax(ir) * COS(theta_ax(itheta))
                   ypt = r_ax(ir) * SIN(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,xpt, ypt, &
                        x_inp, y_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta) )
                ENDIF
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN
       ENDIF

    ENDIF

    DEALLOCATE(i_am_in2D)
    DEALLOCATE(x_inp,y_inp,psi_inp)

  END SUBROUTINE zscatt_cartesian2spherical2D

  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !
  !                          _______
  !                         |__ /   \
  !                          |_ \ |) |
  !                         |___/___/
  !
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!

  SUBROUTINE dscatt_cartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)                :: psi_cart(:, :, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dr(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dth(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    REAL(dp), ALLOCATABLE               :: psi_inp(:)

    !-----------------------------------------------------------------!

    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)

    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in3D(dims_out(1),dims_out(2),dims_out(3)))

    psi_inp = ZERO
    psi_sph = ZERO
    i_am_in3D = 0

    ii = 0
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )

             x_inp(ii) = x_ax(ix)
             y_inp(ii) = y_ax(iy)
             z_inp(ii) = z_ax(iz)
             psi_inp(ii) = psi_cart(ix,iy,iz)

          ENDDO
       ENDDO
    ENDDO

    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))

             IF ((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in3D(ir,itheta,iphi) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO


    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,z_inp, &
            psi_inp,TRIM(method), rank)

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                           psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                           psi_sph_dphi(ir, itheta, iphi), rank )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth).AND. &
            .NOT. PRESENT(psi_sph_dphi) ) THEN

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi), rank=rank )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN
       ENDIF

    ELSE

       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp, z_inp, &
            psi_inp,TRIM(method) )

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                           psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                           psi_sph_dphi(ir, itheta, iphi) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth).AND. &
            .NOT. PRESENT(psi_sph_dphi) ) THEN

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN
       ENDIF

    ENDIF

    DEALLOCATE(i_am_in3D)
    DEALLOCATE(x_inp,y_inp,z_inp,psi_inp)

  END SUBROUTINE dscatt_cartesian2spherical3D

  !-------------------------------------------------------------------------!
  !-------------------------------------------------------------------------!

  SUBROUTINE zscatt_cartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)

    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)

    !---------------------------------------------------------------------!

    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)

    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in3D(dims_out(1),dims_out(2),dims_out(3)))

    psi_inp = ZERO
    psi_sph = ZERO
    i_am_in3D = 0

    ii = 0
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )

             x_inp(ii) = x_ax(ix)
             y_inp(ii) = y_ax(iy)
             z_inp(ii) = z_ax(iz)
             psi_inp(ii) = psi_cart(ix, iy, iz)

          ENDDO
       ENDDO
    ENDDO

    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))
             IF ((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in3D(ir,itheta,iphi) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp, z_inp, &
            psi_inp,TRIM(method), rank)

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                           psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                           psi_sph_dphi(ir, itheta, iphi), rank )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth).AND. &
            .NOT. PRESENT(psi_sph_dphi) ) THEN

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi), rank=rank )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN
       ENDIF

    ELSE

       CALL create_scatt_interpolant(numpts_in,x_inp,y_inp,z_inp, &
            psi_inp,TRIM(method))

       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                           psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                           psi_sph_dphi(ir, itheta, iphi) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth).AND. &
            .NOT. PRESENT(psi_sph_dphi) ) THEN

          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))

                      CALL scatt_interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                   ENDIF
                ENDDO
             ENDDO
          ENDDO

       ELSE
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN
       ENDIF

    ENDIF

    DEALLOCATE(i_am_in3D)
    DEALLOCATE(x_inp,y_inp,z_inp,psi_inp)

  END SUBROUTINE zscatt_cartesian2spherical3D

  !-----------------------------------------------------------------------------!

  !          _ _         _     _         _ ___         _            _         _
  !  __ _  _| (_)_ _  __| |_ _(_)__ __ _| |_  )____ __| |_  ___ _ _(_)__ __ _| |
  ! / _| || | | | ' \/ _` | '_| / _/ _` | |/ /(_-< '_ \ ' \/ -_) '_| / _/ _` | |
  ! \__|\_, |_|_|_||_\__,_|_| |_\__\__,_|_/___/__/ .__/_||_\___|_| |_\__\__,_|_|
  !     |__/                                     |_|

  !-----------------------------------------------------------------------------!

  !                          ___ ___
  !                         |_  )   \
  !                          / /| |) |
  !                         /___|___/
  !

  SUBROUTINE dscatt_cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, method, psi_sph_dx, psi_sph_dy, rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)                :: psi_cyl(:, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dx(:, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(:, :)
    INTEGER                             :: irho, iz, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rhopt, zpt, rpt
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    REAL(dp), ALLOCATABLE               :: rho_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    REAL(dp), ALLOCATABLE               :: psi_inp(:)

    !--------------------------------------------------------------!

    dims_in = SHAPE(psi_cyl)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)

    ALLOCATE(rho_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in2D(dims_out(1),dims_out(2)))

    psi_inp = 0.0_dp
    psi_sph = 0.0_dp
    i_am_in2D = 0

    ii = 0

    DO iz = 1, dims_in(2)
       DO irho = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(rho_ax(irho)**2 + z_ax(iz)**2 )

          rho_inp(ii) = rho_ax(irho)
          z_inp(ii)   = z_ax(iz)

          psi_inp(ii) = psi_cyl(irho, iz)
       ENDDO
    ENDDO

    ! Check if the new axis are within the cartesian domain.
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          rhopt = r_ax(ir) * SIN(theta_ax(itheta))
          zpt = r_ax(ir) * COS(theta_ax(itheta))

          IF ((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             i_am_in2D(ir,itheta) = 1
          ENDIF
       ENDDO
    ENDDO
    
    
    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,rho_inp,z_inp,&
            psi_inp,TRIM(method),rank)
       
       IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
          psi_sph_dx = 0.0_dp
          psi_sph_dy = 0.0_dp

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))

                   CALL scatt_interpolate(numpts_in,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta),rank)
                ENDIF
             ENDDO
          ENDDO

       ELSE

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank)
                ENDIF
             ENDDO
          ENDDO

       ENDIF

    ELSE
       CALL create_scatt_interpolant(numpts_in,rho_inp,z_inp,&
            psi_inp,TRIM(method))

       IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
          psi_sph_dx = 0.0_dp
          psi_sph_dy = 0.0_dp

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts_in,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta))
                ENDIF
             ENDDO
          ENDDO

       ELSE

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta))
                ENDIF
             ENDDO
          ENDDO

       ENDIF

    ENDIF

    DEALLOCATE(i_am_in2D)
    DEALLOCATE(rho_inp, z_inp,psi_inp)


  END SUBROUTINE dscatt_cylindrical2spherical2D

  !----------------------------------------------------------------!

   SUBROUTINE zscatt_cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, method, psi_sph_dx, psi_sph_dy, rank)

    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)             :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(:, :)
    INTEGER                             :: irho, iz, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rhopt, zpt, rpt
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    REAL(dp), ALLOCATABLE               :: rho_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)

    !--------------------------------------------------------------!

    dims_in = SHAPE(psi_cyl)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)

    ALLOCATE(rho_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    ALLOCATE(i_am_in2D(dims_out(1),dims_out(2)))

    psi_inp = ZERO
    psi_sph = ZERO
    i_am_in2D = 0

    ii = 0

    DO iz = 1, dims_in(2)
       DO irho = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(rho_ax(irho)**2 + z_ax(iz)**2 )

          rho_inp(ii) = rho_ax(irho)
          z_inp(ii)   = z_ax(iz)

          psi_inp(ii) = psi_cyl(irho, iz)
       ENDDO
    ENDDO

    ! Check if the new axis are within the cartesian domain.
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          rhopt = r_ax(ir) * SIN(theta_ax(itheta))
          zpt = r_ax(ir) * COS(theta_ax(itheta))

          IF ((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             i_am_in2D(ir,itheta) = 1
          ENDIF
       ENDDO
    ENDDO

    
    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts_in,rho_inp,z_inp,&
            psi_inp,TRIM(method),rank)
       
       IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
          psi_sph_dx = ZERO
          psi_sph_dy = ZERO
          
          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN
                   
                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   
                   CALL scatt_interpolate(numpts_in,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta),rank)
                ENDIF
             ENDDO
          ENDDO

       ELSE

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank)
                ENDIF
             ENDDO
          ENDDO

       ENDIF

    ELSE
       CALL create_scatt_interpolant(numpts_in,rho_inp,z_inp,&
            psi_inp,TRIM(method))
       
       IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
          psi_sph_dx = ZERO
          psi_sph_dy = ZERO
          
          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts_in,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta), &
                        psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta))
                ENDIF
             ENDDO
          ENDDO

       ELSE

          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                IF(i_am_in2D(ir,itheta).EQ.1) THEN

                   rhopt = r_ax(ir) * SIN(theta_ax(itheta))
                   zpt = r_ax(ir) * COS(theta_ax(itheta))
                   CALL scatt_interpolate(numpts,rhopt, zpt, rho_inp, z_inp, &
                        psi_inp, TRIM(method), psi_sph(ir,itheta))
                ENDIF
             ENDDO
          ENDDO

       ENDIF

    ENDIF

    DEALLOCATE(i_am_in2D)
    DEALLOCATE(rho_inp, z_inp,psi_inp)


  END SUBROUTINE zscatt_cylindrical2spherical2D

  !----------------------------------------------------------------!

END MODULE scattcoords
