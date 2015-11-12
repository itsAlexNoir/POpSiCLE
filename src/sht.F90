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
  ! USE fftw
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: initialize_spherical_harmonics
  PUBLIC       :: make_sht
  
  ! Interfaces
  
  INTERFACE make_sht
     MODULE PROCEDURE make_sht_2D
     !MODULE PROCEDURE make_sht_3D
  END INTERFACE make_sht
  
  ! Public variables
  
  REAL(dp), ALLOCATABLE, PUBLIC  :: legenpl(:, :, :)
  REAL(dp), ALLOCATABLE, PUBLIC  :: normfact(:, :)
  
CONTAINS
  
  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE initialize_spherical_harmonics
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
  
  SUBROUTINE initialize_spherical_harmonics(lmax, axis)
    
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
    CALL make_rnormal(normfact, factlog, lmax, maxthetapts)
    
    ! Deallocate
    DEALLOCATE(factlog)
    
  END SUBROUTINE initialize_spherical_harmonics
  
  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE make_sht_2D
  !
  !>  \brief This subroutine returns the Gauss nodes of a particular
  !>  range of values, and the correspondent Gauss weights
  !>  \brief This subroutine make the SHT for a function that depends
  !>  only on theta. It uses a Gauss-Legendre algorithm, in other words,
  !>  it computes a Gaussian quadrature with gauss nodes and weights. 
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lower_bound Minimum value of the passed array
  !> \param[in] upper_bound Maximum value of the passed array
  !> \param[out] axis Returned axis with gaussian nodes.
  !> \param[out] weights Returned weights array.
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE make_sht_2D(func, axis, weights, lmax, rad_coeffl0)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)    :: func(:, :)
    REAL(dp), INTENT(IN)       :: axis(:)
    REAL(dp), INTENT(IN)       :: weights(:)
    INTEGER, INTENT(IN)        :: lmax
    COMPLEX(dp), INTENT(OUT)   :: rad_coeffl0(1:5, 0:lmax)
    
    INTEGER                    :: numrpts, numthetapts
    REAL(dp)                   :: sum
    INTEGER                    :: il, ir, itheta
    
    ! What is the maximum number of points in r and theta?
    numrpts = SIZE(func,1)
    numthetapts = SIZE(axis)
    
    ! Initialize array to zero
    rad_coeffl0 = ZERO
    
    DO ir = 1, numrpts
       DO il = 0, lmax
          sum = 0.0_dp
          DO itheta = 1, numthetapts
             sum = sum + func(ir,itheta) * &
                  normfact(0,il) * legenpl(itheta-1,0,il) * weights(itheta)
          ENDDO
          rad_coeffl0(ir,il) = sum
       ENDDO
    ENDDO
           
  END SUBROUTINE make_sht_2D
  
  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE make_sht_3D
  !
  !>  \brief This subroutine returns the Gauss nodes of a particular
  !>  range of values, and the correspondent Gauss weights
  !>  \brief This subroutine make the SHT for a function that depends
  !>  on theta and phi. It uses a Gauss-Legendre algorithm, in other words,
  !>  it computes a Gaussian quadrature with gauss nodes and weights. 
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lower_bound Minimum value of the passed array
  !> \param[in] upper_bound Maximum value of the passed array
  !> \param[out] axis Returned axis with gaussian nodes.
  !> \param[out] weights Returned weights array.
  !
  !=======================================================================
  !=======================================================================
  
!!$  SUBROUTINE make_sht_3D(func, th_axis, weights, phi_axis, lmax, rad_coefflm)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    COMPLEX(dp), INTENT(IN)     :: func(:, :, :)
!!$    REAL(dp), INTENT(IN)        :: th_axis(:)
!!$    REAL(dp), INTENT(IN)        :: weights(:)
!!$    REAL(dp), INTENT(IN)        :: phi_axis(:)
!!$    INTEGER, INTENT(IN)         :: lmax
!!$    COMPLEX(dp), INTENT(OUT)    :: rad_coefflm(: ,:, :)
!!$    
!!$    COMPLEX(dp), ALLOCATABLE   :: coeffm(:)
!!$    COMPLEX(dp), ALLOCATABLE   :: coeffmtheta(:, :)
!!$    INTEGER                    :: numrpts, numthetapts
!!$    INTEGER                    :: numphipts
!!$    INTEGER                    :: il, ir, itheta, iphi
!!$    
!!$    ! What is the maximum number of points in r and theta?
!!$    numrpts = SIZE(func,1)
!!$    numthetapts = SIZE(th_axis)
!!$    numphipts = SIZE(phi_axis)
!!$    ! Initialize array to zero
!!$    rad_coefflm = ZERO
!!$    
!!$    ALLOCATE(coeffm(1:numphipts))
!!$    ALLOCATE(coeffmtheta(1:numthetapts,1:numphipts))
!!$    
!!$    DO ir = 1, numrpts
!!$       
!!$       ! FFT to calculate m number
!!$       DO itheta = 1, numthetapts
!!$          coeffm = ZERO
!!$          CALL FourierTransform(func(ir,itheta,:),coeffm,1,(/numphipts/))
!!$          ! Reorder coeffm array after fourier transform
!!$          shift = numphipts - (numphipts+1)/2
!!$          DO im = 1, shift
!!$             coeffmtheta(itheta,im) = coeffm(im+shift)
!!$          ENDDO
!!$          
!!$          DO im = shift+1, numphipts
!!$             coeffmtheta(itheta,im) = coeffm(im-shift)
!!$          ENDDO
!!$          
!!$       ENDDO
!!$       
!!$       ! Now it is the turn of l
!!$       DO il = 0, lmax
!!$          DO im = 
!!$             sum = 0.0_dp
!!$             DO itheta = 1, numthetapts
!!$                sum = sum + coeffmtheta(itheta,im) * &
!!$                     normfact(im,il) * legenpl(itheta-1,im,il) * weights(itheta)
!!$             ENDDO
!!$             rad_coefflm(ir,il,im) = sum
!!$          ENDDO
!!$       ENDDO
!!$       
!!$    ENDDO ! End of r loop
!!$    
!!$    ! Deallocate auxiliary arrays
!!$    DEALLOCATE(coeffm,coeffmtheta)
!!$ 
!!$  END SUBROUTINE make_sht_3D
  
END MODULE sht
