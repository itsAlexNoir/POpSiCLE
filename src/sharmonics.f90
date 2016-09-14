!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: sharmonics
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 05/10/2015
!
! DESCRIPTION:
!> \brief Spherical harmonics module
!> \details This module contains routines to
!> decompose a function in spherical coordinates into spherical
!> harmonics. It contains two types: For functions that depends only
!> in theta, and also for functions in 3D, for theta and phi.
!
!------------------------------------------------------------------------------
MODULE sharmonics
  
  USE constants_pop
  USE gaussleg
  !USE fourier
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: initialize_legendre_stuff
  PUBLIC       :: create_spherical_harmonics
  PUBLIC       :: create_spherical_harmonics_couplings
  PUBLIC       :: make_sht
  
  ! Public variables
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: legenpl(:, :, :)
  REAL(dp), ALLOCATABLE, PUBLIC    :: normfact(:, :)
  
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
    ALLOCATE(legenpl(0:maxthetapts,0:lmax,0:lmax))
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
  !> \param[in] costheta_axis The abscissa (cos theta) for Legendre polynomials
  !> \param[in] phi_axis Phi angle axis
  !> \param[out] sph_harmonics The spherical harmonics matrix
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE create_spherical_harmonics(sph_harmonics, lmax, costheta_axis, phi_axis)
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)      :: lmax
    COMPLEX(dp), INTENT(OUT) :: sph_harmonics(:, :, -lmax:, 0:)
    REAL(dp), INTENT(IN)     :: costheta_axis(:)
    REAL(dp), INTENT(IN)     :: phi_axis(:)
    
    REAL(dp), ALLOCATABLE    :: legendre_polynomial(:, :, :)
    REAL(dp), ALLOCATABLE    :: normfactor(:, :)
    REAL(dp), ALLOCATABLE    :: factlog(:)
    INTEGER                  :: maxthetapts
    INTEGER                  :: maxphipts
    INTEGER                  :: maxfactlog
    INTEGER                  :: il, im, mm
    INTEGER                  :: itheta, iphi

    !-------------------------------------------!
    
    ! Number of points in (theta) axis
    ! The -1 is becasue it assumes that
    ! the array start at 0, not 1.
    maxthetapts = SIZE(costheta_axis)
    maxphipts = SIZE(phi_axis)
    ! Maximum factorial needed, 3*lmax + 1
    maxfactlog = 20 * lmax
    
    ! Allocate arrays
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
       DO im = 0, il
          
          DO itheta = 1, maxthetapts
             DO iphi = 1, maxphipts
                
                sph_harmonics(iphi,itheta,im,il) = &
                     ((-1.0_dp) ** ABS(im)) * normfactor(im,il) * &
                     legendre_polynomial(itheta - 1,im,il) * &
                     EXP(ZIMAGONE * im * phi_axis(iphi))
                
                IF ( im .GT. 0) THEN
                   sph_harmonics(iphi,itheta,-im,il) = &
                        ((-1.0_dp) ** ABS(im)) * &
                        CONJG(sph_harmonics(iphi, itheta, im, il))
                ENDIF
                
             ENDDO
          ENDDO
          
       ENDDO
    ENDDO
    
    
    ! Deallocate
    DEALLOCATE(legendre_polynomial,normfactor)
    DEALLOCATE(factlog)
    
  END SUBROUTINE create_spherical_harmonics
  
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
  !> \param[in] costheta_axis The abscissa (cos theta) for Legendre polynomials
  !> \param[in] phi_axis Phi angle axis
  !> \param[out] sph_harmonics The spherical harmonics matrix
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE create_spherical_harmonics_couplings( xcoupling, ycoupling, zcoupling, &
       sph_harmonics, lmax, theta_axis, phi_axis, theta_w, phi_w)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)           :: lmax
    REAL(dp), INTENT(OUT)         :: xcoupling(-lmax:, 0:, -lmax:, 0:)
    REAL(dp), INTENT(OUT)         :: ycoupling(-lmax:, 0:, -lmax:, 0:)
    REAL(dp), INTENT(OUT)         :: zcoupling(-lmax:, 0:, -lmax:, 0:)
    COMPLEX(dp), INTENT(IN)       :: sph_harmonics(:, :, -lmax:, 0:)
    REAL(dp), INTENT(IN)          :: theta_axis(:)
    REAL(dp), INTENT(IN)          :: phi_axis(:)
    REAL(dp), INTENT(IN)          :: theta_w(:)
    REAL(dp), INTENT(IN)          :: phi_w(:)

    COMPLEX(dp)                   :: suma, sumb
    INTEGER                       :: numthetapts, numphipts
    INTEGER                       :: il, im, ill, imm
    INTEGER                       :: itheta, iphi
    
    !------------------------------------------------------!
    
    numthetapts = SIZE(theta_axis)
    numphipts   = SIZE(phi_axis)
    
    xcoupling = 0.0_dp
    ycoupling = 0.0_dp
    zcoupling = 0.0_dp    
    
    ! Calculate matrix elements for laser-matter couling in
    ! x and y direction laser polarisation.
    
    DO ill = 0, lmax
       DO imm = -ill, ill
          
          DO il = 0, lmax
             DO im = -il, il
                
                IF (ABS(il - ill).EQ.1 .AND. ABS(im - imm).EQ.1) THEN
                   
                   ! Calculate matrix elements between spherical harmonics
                   suma = ZERO
                   sumb = ZERO
                   
                   DO iphi = 1, numphipts
                      DO itheta = 1, numthetapts
                         
                         suma = suma + &
                              CONJG(sph_harmonics(iphi,itheta,im,il)) * &
                              sph_harmonics(iphi,itheta,imm,ill) * &
                              SIN(theta_axis(itheta)) * COS(phi_axis(iphi)) * &
                              theta_w(itheta) * phi_w(iphi)
                         
                         sumb = sumb + &
                              CONJG(sph_harmonics(iphi,itheta,im,il)) * &
                              sph_harmonics(iphi,itheta,imm,ill) * &
                              SIN(theta_axis(itheta)) * SIN(phi_axis(iphi)) * &
                              theta_w(itheta) * phi_w(iphi)
                         
                      ENDDO
                   ENDDO
                   
                   xcoupling(im,il,imm,ill) = suma
                   ycoupling(im,il,imm,ill) = AIMAG(sumb)
                   
                ENDIF
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! Calculate matrix elements for laser-matter couling in
    ! z direction laser polarisation.
    
    DO ill = 0, lmax
       DO imm = -ill, ill
          DO il = 0, lmax
             DO im = -il, il
                
                IF (ABS(il - ill).EQ.1 .AND. ABS(im - imm).EQ.0) THEN
                   ! Calculate matrix elements between spherical harmonics
                   
                   suma = ZERO
                   
                   DO iphi = 1, numphipts
                      DO itheta = 1, numthetapts
                         
                         suma = suma + &
                              CONJG(sph_harmonics(iphi,itheta,im,il)) * &
                              sph_harmonics(iphi,itheta,imm,ill) * &
                              COS(theta_axis(itheta)) * &
                              theta_w(itheta) * phi_w(iphi)
                         
                      ENDDO
                   ENDDO
                   
                   zcoupling(im,il,imm,ill) = suma
                   
                ENDIF
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    
  END SUBROUTINE create_spherical_harmonics_couplings
  
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
  
  SUBROUTINE make_sht(func, sph_harmonics, lmax, th_weights, phi_weights, func_lm)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)        :: func(:, :)
    INTEGER, INTENT(IN)            :: lmax
    COMPLEX(dp), INTENT(IN)        :: sph_harmonics(:, :, -lmax:, 0:)
    REAL(dp), INTENT(IN)           :: th_weights(:)
    REAL(dp), INTENT(IN)           :: phi_weights(:)
    COMPLEX(dp), INTENT(OUT)       :: func_lm(-lmax:, 0:)
    
    INTEGER                        :: dims(2)
    COMPLEX(dp), ALLOCATABLE       :: coeffm(:), gm(:, :)
    COMPLEX(dp)                    :: sum
    INTEGER                        :: maxthetapts, maxphipts
    INTEGER                        :: il, im, itheta, iphi
    
    !----------------------------------------------------!
    
    ! What is the maximum number of theta and phi points?
    !numthetapts = SIZE(axis)
    dims = SHAPE(func)
    maxthetapts = dims(1)
    maxphipts   = dims(2)
    
    ! Initialize array to zero
    func_lm = ZERO
    
    DO il = 0, lmax
       DO im = -il, il
          
          sum = ZERO
          
          DO iphi = 1, maxphipts
             DO itheta = 1, maxthetapts
                
                sum = sum + func(itheta, iphi) *                               &
                     CONJG(sph_harmonics(iphi, itheta, im, il)) *             &
                     th_weights(itheta) * phi_weights(iphi)
                
             ENDDO
          ENDDO
          func_lm(im,il) = sum
       ENDDO
    ENDDO
    
       
!!$       ALLOCATE(coeffm(-lmax:lmax))
!!$       ALLOCATE(gm(1:maxthetapts,-lmax:lmax))
!!$       coeffm = ZERO
!!$       gm     = ZERO
!!$       dphi = twopi / REAL(maxphipts,dp)
!!$       
!!$       DO itheta = 1, maxthetapts
!!$          coeffm = ZERO
!!$          CALL FourierTransform(func(itheta,:),coeffm,1,(/maxphipts/))
!!$          
!!$          ! Reorder coeffm array after fourier transform
!!$          CALL fftshift(coeffm(-lmax:lmax),gm(itheta,-lmax:lmax))
!!$          
!!$       ENDDO
!!$       gm = gm * dphi
!!$       
!!$       ! Now, the Gauss-Legendre quadrature
!!$       DO il = 0, lmax
!!$          DO im = -il, il
!!$             sum = ZERO
!!$             
!!$             DO itheta = 1, maxthetapts
!!$                sum = sum + gm(itheta, im) * &
!!$                     normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) &
!!$                     * (-1.0_dp)**ABS(im) * weights(itheta)
!!$             ENDDO
!!$             
!!$             func_lm(im,il) = sum
!!$          ENDDO
!!$       ENDDO
!!$       
!!$       ! Deallocate stuff
!!$       DEALLOCATE(gm)
!!$       DEALLOCATE(coeffm)
!!$       
       
   
  END SUBROUTINE make_sht
  
  !*****************************************************************!
END MODULE sharmonics
