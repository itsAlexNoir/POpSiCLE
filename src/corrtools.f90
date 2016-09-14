!-----------------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------------
!                                                                              
! Module Correlated tools                                                             
!
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 01/03/2016
!
! DESCRIPTION:
!> \brief Auxiliary functions needed to calculate correlated spectra.
!> \details This module contains routines related to the correlated spectra,
!>          such as creating sine basis and a Runge-Kutta integrator.
! 
!------------------------------------------------------------------------------!
MODULE corrtools

  USE constants_pop
  
  IMPLICIT NONE
  
  PUBLIC             :: create_sine_basis

  !------------------------------------------! 
  
CONTAINS
  
  SUBROUTINE  create_sine_basis(sineb, order, r_axis)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(OUT)         :: sineb(:, :)
    INTEGER, INTENT(IN)           :: order
    REAL(dp), INTENT(IN)          :: r_axis(:)
    
    REAL(dp)                      :: maxr
    INTEGER                       :: numrpts
    INTEGER                       :: ir, iorder
    
    !-----------------------------------------------!
    
    maxr = MAXVAL(r_axis)
    numrpts = SIZE(r_axis)
    
    DO iorder = 1, order
       DO ir = 1, numrpts 
          sineb(ir,iorder) = SIN(REAL(iorder,dp) * pi * &
               r_axis(ir) / maxr)
       ENDDO
    ENDDO
    
  END SUBROUTINE create_sine_basis
  
  !-------------------------------------!

!!$  SUBROUTINE rungekutta4(sinecoeff, sinecoeffdt, energy_sineb, elec_flux, &
!!$       elec_flux_h2, elec_flux_h)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    COMPLEX(dp), INTENT(IN)      :: sinecoeff(:, :, :, :)
!!$    COMPLEX(dp), INTENT(OUT)     :: sinecoeffdt(:, :, :, :)
!!$    COMPLEX(dp), INTENT(IN)      :: energy_sineb(:, :)
!!$    COMPLEX(dp), INTENT(IN)      :: elec_flux(:, :, :, :)
!!$    COMPLEX(dp), INTENT(IN)      :: elec_flux_h2(:, :, :, :)
!!$    COMPLEX(dp), INTENT(IN)      :: elec_flux_h(:, :, :, :)
!!$    
!!$    INTEGER                      :: coeff_dims(4)
!!$    COMPLEX(dp), ALLOCATABLE     :: k1(:, :, :, :)
!!$    COMPLEX(dp), ALLOCATABLE     :: k2(:, :, :, :)
!!$    COMPLEX(dp), ALLOCATABLE     :: k3(:, :, :, :)
!!$    COMPLEX(dp), ALLOCATABLE     :: k4(:, :, :, :)
!!$    INTEGER                      :: numkpts, numthetapts
!!$    INTEGER                      :: numphipts,numorder
!!$    
!!$    !--------------------------------------!
!!$    
!!$    coeffs_dims = SHAPE(sinecoeff)
!!$
!!$    numkpts     = coeff_dims(1)
!!$    numthetapts = coeff_dims(2)
!!$    numphipts   = coeff_dims(3)
!!$    numorder    = coeff_dims(4)
!!$    
!!$    ! Allocate auxiliary arrays (RK4 k's)
!!$    ALLOCATE(k1(1:numkepts,1:numthetapts,1:numphipts,1:numorder))
!!$    ALLOCATE(k2(1:numkepts,1:numthetapts,1:numphipts,1:numorder))
!!$    ALLOCATE(k3(1:numkepts,1:numthetapts,1:numphipts,1:numorder))
!!$    ALLOCATE(k4(1:numkepts,1:numthetapts,1:numphipts,1:numorder))
!!$    k1 = ZERO
!!$    k2 = ZERO
!!$    k3 = ZERO
!!$    k4 = ZERO
!!$
!!$    DO im = 1, numorder
!!$       DO iphi = 1, numphipts
!!$          DO itheta = 1, numthetapts
!!$             DO ik = 1, numkpts
!!$                
!!$                DO imm = 1, numorder
!!$                   k1(ik,itheta,iphi,im) = k1(ik,itheta,iphi,im) + &
!!$                        energy_sineb(imm,im) * sinecoeff(ik,itheta,iphi,imm)
!!$                ENDDO
!!$                
!!$                k1(ik,itheta,iphi,im) = k1(ik,itheta,iphi,im) - &
!!$                     elec_flux(ik,itheta,iphi,im)
!!$                
!!$             ENDDO
!!$          ENDDO
!!$       ENDDO
!!$    ENDDO
!!$    k1 = k1 * (-ZIMAGONE)
!!$    
!!$    
!!$    
!!$    ! Deallocate
!!$    DEALLOCATE(k1,k2,k3,k4)
!!$    
!!$  END SUBROUTINE rungekutta4
    
END MODULE corrtools
