!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: sht
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 05/10/2015
!
! DESCRIPTION:
!> \brief Spherical harmonic transform
!> \details This module contains routines to
!> decompose a function in spherical coordinates into spherical
!> harmonics. It contains two types: For functions that depends only
!> in theta, and also for functions in 3D, for theta and phi.
!
!------------------------------------------------------------------------------
MODULE sht
  
  USE constants
  USE gaussleg
  USE fourier
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC       :: initialize_legendre_stuff
  PUBLIC       :: initialize_spherical_harmonics
  PUBLIC       :: make_sht
    
  ! Public variables
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: legenpl(:, :, :)
  REAL(dp), ALLOCATABLE, PUBLIC    :: normfact(:, :)
  COMPLEX(dp), ALLOCATABLE, PUBLIC :: sph_harmonics(:, :, :, :)
  
CONTAINS
  
  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE initialize_legendre_stuff
  !
  !>  \brief This subroutine initializes the Legendre polynomial and
  !>  his correpondent normalization factors for the maximum
  !>  angular momentum required.
  !>  The routines that compute the polynomial and the factor are located
  !>  at the module gaussleg. This is only a wrapper to call them once at
  !>  the beginning of the simulation.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lmax Maximum angular momenta
  !> \param[in] axis The abscissa (cos theta) for Legendre polynomials
  !> \param[out] legenpl The array containing all the polynomials needed
  !> \param[out] normfact The array containing all the normalization factors
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE initialize_legendre_stuff(lmax, axis)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)     :: lmax
    REAL(dp), INTENT(IN)    :: axis(:)
    
    REAL(dp), ALLOCATABLE   :: factlog(:)
    INTEGER                 :: maxthetapts
    INTEGER                 :: maxfactlog
    INTEGER                 :: il
    
    ! Number of points in (theta) axis
    ! The -1 is becasue it assumes that
    ! the array start at 0, not 1.
    maxthetapts = SIZE(axis) - 1
    ! Maximum factorial needed, 3*lmax + 1
    maxfactlog = 3 * lmax + 1
    
    ! Allocate arrays
    ALLOCATE(factlog(0:maxfactlog))
    ALLOCATE(legenpl(0:maxthetapts, 0:lmax, 0:lmax))
    ALLOCATE(normfact(0:lmax, 0:lmax))
    
    ! Calculate the Legendre polynomial
    CALL make_legendre(legenpl, axis, lmax, maxthetapts)
    
    ! Calculate the factorials needed for the normalization
    ! factors.
    CALL make_factlog(factlog, maxfactlog)
    
    ! Calculate the Normalization factor
    CALL make_rnormal(normfact, factlog, lmax, maxfactlog)
    
    ! Deallocate
    DEALLOCATE(factlog)
    
  END SUBROUTINE initialize_legendre_stuff

  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE initialize_spherical_harmonics
  !
  !>  \brief This subroutine initializes the spherical harmonic matrix.
  !>  The routines that compute the polynomial and the factor are located
  !>  at the module gaussleg. This is only a wrapper to call them once at
  !>  the beginning of the simulation.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lmax Maximum angular momenta
  !> \param[in] mmax Maximum magnetic angular momenta
  !> \param[in] costheta_axis The abscissa (cos theta) for Legendre polynomials
  !> \param[in] phi_axis Phi angle axis
  !> \param[out] sph_harmonics The spherical harmonics matrix
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE initialize_spherical_harmonics(lmax, mmax, costheta_axis, phi_axis)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)     :: lmax
    INTEGER, INTENT(IN)     :: mmax    
    REAL(dp), INTENT(IN)    :: costheta_axis(:)
    REAL(dp), INTENT(IN)    :: phi_axis(:)

    REAL(dp), ALLOCATABLE   :: legendre_polynomial(:, :, :)
    REAL(dp), ALLOCATABLE   :: normfactor(:, :)
    REAL(dp), ALLOCATABLE   :: factlog(:)
    INTEGER                 :: maxthetapts
    INTEGER                 :: maxphipts
    INTEGER                 :: maxfactlog
    INTEGER                 :: mmin
    INTEGER                 :: il, im
    INTEGER                 :: itheta, iphi
    
    !-------------------------------------------!
    
    ! Number of points in (theta) axis
    ! The -1 is becasue it assumes that
    ! the array start at 0, not 1.
    maxthetapts = SIZE(costheta_axis)
    maxphipts = SIZE(phi_axis)
    ! Maximum factorial needed, 3*lmax + 1
    maxfactlog = 3 * lmax + 1
    ! Set minimum magnetic angular momentum
    mmin = - mmax
    
    ! Allocate arrays
    IF(mmax.EQ.0) THEN
    ALLOCATE(sph_harmonics(1:maxphipts,1:maxthetapts,0:1,0:lmax))
    ELSE
       ALLOCATE(sph_harmonics(1:maxphipts,1:maxthetapts,mmin:mmax,0:lmax))
    ENDIF
    ALLOCATE(factlog(0:maxfactlog))
    ALLOCATE(legendre_polynomial(0:maxthetapts-1, 0:lmax, 0:lmax))
    ALLOCATE(normfactor(0:lmax, 0:lmax))
    sph_harmonics = ZERO
    
    ! Calculate the Legendre polynomial
    CALL make_legendre(legendre_polynomial, costheta_axis, lmax, maxthetapts-1)
    
    ! Calculate the factorials needed for the normalization
    ! factors.
    CALL make_factlog(factlog, maxfactlog)
    
    ! Calculate the Normalization factor
    CALL make_rnormal(normfactor, factlog, lmax, maxfactlog)

    DO il = 0, lmax
       DO im = mmin, mmax
          DO itheta = 1, maxthetapts
             DO iphi = 1, maxphipts
                
                IF(im.LT.0) THEN
                   sph_harmonics(iphi,itheta,im,il) = (-1.0_dp)**ABS(im) &
                        * normfactor(ABS(im),il) * &
                        legendre_polynomial(itheta-1,ABS(im),il) * &
                        EXP(ZIMAGONE * im * phi_axis(iphi))
                ELSE
                   sph_harmonics(iphi,itheta,im,il) = normfactor(im,il) * &
                        legendre_polynomial(itheta-1,im,il) * &
                        EXP(ZIMAGONE * im * phi_axis(iphi))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! Deallocate
    DEALLOCATE(legendre_polynomial,normfactor)
    DEALLOCATE(factlog)
    
  END SUBROUTINE initialize_spherical_harmonics

  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE make_sht
  !
  !>  \brief This subroutine performs the (fast) Spherical
  !>    Harmonics Transform (SHT). It uses a FFT to decompose
  !>    into the magnetic quantum number m, and a Gauss-Legendre
  !>    quadrature to decompose into angular momentum number l.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] func Function to be transformed
  !> \param[in] lmax Maximum angular momentum
  !> \param[in] weights Gauss weigths needed to perform the quadrature
  !> \param[out] func_lm Transformed function, decomposed into l,m
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE make_sht(func, lmax, mmax, weights, func_lm)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: func(:, :)
    INTEGER, INTENT(IN)            :: lmax, mmax
    REAL(dp), INTENT(IN)           :: weights(:)
    COMPLEX(dp), INTENT(OUT)       :: func_lm(-mmax:, 0:)
    
    INTEGER                        :: dims(2)
    COMPLEX(dp), ALLOCATABLE       :: coeffm(:), gm(:, :)
    COMPLEX(dp)                    :: sum
    REAL(dp)                       :: dphi
    INTEGER                        :: numthetapts, numphipts
    INTEGER                        :: mmin
    INTEGER                        :: il, im, itheta
    
    !----------------------------------------------------!
    
    ! What is the maximum number of theta and phi points?
    !numthetapts = SIZE(axis)
    dims = SHAPE(func)
    numthetapts = dims(1)
    numphipts   = dims(2)
    
    mmin = -mmax
    ! Allocate auxiliary arrays
    IF(mmax.NE.0) THEN
       ALLOCATE(coeffm(mmin:mmax))
       dphi = twopi / REAL(numphipts,dp)
    ENDIF
    
    ALLOCATE(gm(1:numthetapts,mmin:mmax))
    
    ! Initialize array to zero
    func_lm = ZERO
    gm = ZERO
    
    ! Do first the FFT (if necessary)
    IF(numphipts.EQ.1) THEN
       DO itheta = 1, numthetapts
          gm(itheta,0) = func(itheta,1)
       ENDDO
    ELSE
       DO itheta = 1, numthetapts
          coeffm = ZERO
          CALL FourierTransform(func(itheta,:),coeffm,1,(/numphipts/))
          
          ! Reorder coeffm array after fourier transform
          CALL fftshift(coeffm(mmin:mmax),gm(itheta,mmin:mmax))
          
       ENDDO
       gm = gm * dphi
    ENDIF
    
    ! Now, the Gauss-Legendre quadrature
    DO il = 0, lmax
       DO im = mmin, mmax
          sum = ZERO
          IF(im.LT.0) THEN
             DO itheta = 1, numthetapts
                sum = sum + gm(itheta, im) * &
                     normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) &
                     * (-1.0_dp)**ABS(im) * weights(itheta)
             ENDDO
          ELSE
             DO itheta = 1, numthetapts
                sum = sum + gm(itheta, im) * &
                     normfact(im,il) * legenpl(itheta-1,im,il) &
                     * weights(itheta)
             ENDDO
          ENDIF
          func_lm(im,il) = sum
       ENDDO
    ENDDO
    
    IF(mmax.EQ.0) &
         func_lm = func_lm * twopi
    
    ! Deallocate stuff
    DEALLOCATE(gm)
    IF(mmax.NE.0) &
         DEALLOCATE(coeffm)
    
  END SUBROUTINE make_sht

  !*****************************************************************!
END MODULE sht
