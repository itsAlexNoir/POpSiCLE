!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: cubic_interp
!> \author Alejandro de la Calle, Queen's University Belfast,
!> \author Daniel Dundas, Queen's University Belfast,
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Interpolating routines,
!> \details This module contains all the subprograms needed for
!> an interpolation of the wavefunction, including routines for
!> multidimensions grid, Many of the routines in this module have
!> been extract form the numerical recipes books,
!
!------------------------------------------------------------------------------
MODULE cubic_interp
  
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC   :: initialize_bicubic_matrix
  PUBLIC   :: initialize_tricubic_matrix
  PUBLIC   :: bicubic_interpolation
  PUBLIC   :: tricubic_interpolation
  PUBLIC   :: delete_bicubic_matrix
  PUBLIC   :: delete_tricubic_matrix
  
  REAL(dp), ALLOCATABLE, PUBLIC  :: wbi(:, :)
  REAL(dp), ALLOCATABLE, PUBLIC  :: wtri(:, :)
  
  
CONTAINS
  
  !**********************************************************************!
  
  SUBROUTINE initialize_bicubic_matrix( )
    IMPLICIT NONE
    
    ALLOCATE(wbi(16,16))
    
    wbi(:,1)  = (/1,0,-3,2,0,0,0,0,-3,0,9,-6,2,0,-6,4/)
    wbi(:,2)  = (/0,0,0,0,0,0,0,0,3,0,-9,6,-2,0,6,-4/)
    wbi(:,3)  = (/0,0,0,0,0,0,0,0,0,0,9,-6,0,0,-6,4/)
    wbi(:,4)  = (/0,0,3,-2,0,0,0,0,0,0,-9,6,0,0,6,-4/)
    wbi(:,5)  = (/0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2/)
    wbi(:,6)  = (/0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2/)
    wbi(:,7)  = (/0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2/)
    wbi(:,8)  = (/0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2/)
    wbi(:,9)  = (/0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2/)
    wbi(:,10) = (/0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2/)
    wbi(:,11) = (/0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2/)
    wbi(:,12) = (/0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2/)
    wbi(:,13) = (/0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1/)
    wbi(:,14) = (/0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1/)
    wbi(:,15) = (/0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1/)
    wbi(:,16) = (/0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1/)
    
  END SUBROUTINE initialize_bicubic_matrix
  
  !*************************************************!
  
  SUBROUTINE initialize_tricubic_matrix( )
    IMPLICIT NONE
    
    ALLOCATE(wtri(64,64))

    wtri(:,1) = (/1,0,-3,2,0,0,0,0,-3,0,9,-6,2,0,-6,4,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,-3,0,9,-6,0,0,0,0,9,0,-27,18,-6,0,18,&
         -12,2,0,-6,4,0,0,0,0,-6,0,18,-12,4,0,-12,8/)
    wtri(:,2) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,3,0,-9,6,0,0,0,0,-9,0,27,-18,6,0,-18,12,-2,0,6,&
         -4,0,0,0,0,6,0,-18,12,-4,0,12,-8/)
    wtri(:,3) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,0,-27,18,-6,0,18,-12,0,0,0,0,0,&
         0,0,0,-6,0,18,-12,4,0,-12,8/)
    wtri(:,4) = (/0,0,0,0,0,0,0,0,3,0,-9,6,-2,0,6,-4,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-9,0,27,-18,6,0,-18,12,0,0,0,0,0,&
         0,0,0,6,0,-18,12,-4,0,12,-8/)
    wtri(:,5) = (/0,0,3,-2,0,0,0,0,0,0,-9,6,0,0,6,-4,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,-9,6,0,0,0,0,0,0,27,-18,0,0,-18,12,0,0,6,-4,0,0,&
         0,0,0,0,-18,12,0,0,12,-8/)
    wtri(:,6) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,9,-6,0,0,0,0,0,0,-27,18,0,0,18,-12,0,0,-6,4,0,0,0,0,0,&
         0,18,-12,0,0,-12,8/)
    wtri(:,7) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,0,27,-18,0,0,-18,12,0,0,0,0,0,0,0,0,0,0,&
         -18,12,0,0,12,-8/)
    wtri(:,8) = (/0,0,0,0,0,0,0,0,0,0,9,-6,0,0,-6,4,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,0,0,-27,18,0,0,18,-12,0,0,0,0,0,0,0,0,0,0,&
         18,-12,0,0,-12,8/)
    wtri(:,9) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-3,2,0,0,0,0,-3,0,9,&
         -6,2,0,-6,4,-2,0,6,-4,0,0,0,0,6,0,-18,12,-4,0,12,-8,1,0,-3,2,0,0,0,&
         0,-3,0,9,-6,2,0,-6,4/)
    wtri(:,10) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,-1,0,3,-2,0,0,0,0,3,0,-9,6,-2,0,6,-4,1,0,-3,2,0,0,0,0,-3,0,9,&
         -6,2,0,-6,4/)
    wtri(:,11) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,-3,0,9,-6,2,0,-6,4,0,0,0,0,0,0,0,0,3,0,-9,6,-2,&
         0,6,-4/)
    wtri(:,12) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,-9,6,&
         -2,0,6,-4,0,0,0,0,0,0,0,0,-6,0,18,-12,4,0,-12,8,0,0,0,0,0,0,0,0,3,0,&
         -9,6,-2,0,6,-4/)
    wtri(:,13) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,-2,0,0,0,0,0,0,-9,6,&
         0,0,6,-4,0,0,-6,4,0,0,0,0,0,0,18,-12,0,0,-12,8,0,0,3,-2,0,0,0,0,0,0,&
         -9,6,0,0,6,-4/)
    wtri(:,14) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,-3,2,0,0,0,0,0,0,9,-6,0,0,-6,4,0,0,3,-2,0,0,0,0,0,0,-9,6,&
         0,0,6,-4/)
    wtri(:,15) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,0,-9,6,0,0,6,-4,0,0,0,0,0,0,0,0,0,0,9,-6,0,0,&
         -6,4/)
    wtri(:,16) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,9,-6,0,&
         0,-6,4,0,0,0,0,0,0,0,0,0,0,-18,12,0,0,12,-8,0,0,0,0,0,0,0,0,0,0,9,-6,&
         0,0,-6,4/)
    wtri(:,17) = (/0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,-3,0,9,-6,6,0,-18,12,-3,0,9,-6,0,0,0,0,2,0,-6,4,-4,0,&
         12,-8,2,0,-6,4/)
    wtri(:,18) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,3,0,-9,6,-6,0,18,-12,3,0,-9,6,0,0,0,0,-2,0,6,-4,4,0,-12,8,&
         -2,0,6,-4/)
    wtri(:,19) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,-3,0,9,-6,3,0,-9,6,0,0,0,0,0,0,0,0,2,0,-6,4,-2,0,6,-4/)
    wtri(:,20) = (/0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,3,0,-9,6,-3,0,9,-6,0,0,0,0,0,0,0,0,-2,0,6,-4,2,0,-6,4/)
    wtri(:,21) = (/0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,-9,6,0,0,18,-12,0,0,-9,6,0,0,0,0,0,0,6,-4,0,0,-12,8,0,0,6,-4/)
    wtri(:,22) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,9,-6,0,0,-18,12,0,0,9,-6,0,0,0,0,0,0,-6,4,0,0,12,-8,0,0,-6,4/)
    wtri(:,23) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,-9,6,0,0,9,-6,0,0,0,0,0,0,0,0,0,0,6,-4,0,0,-6,4/)
    wtri(:,24) = (/0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,9,-6,0,0,-9,6,0,0,0,0,0,0,0,0,0,0,-6,4,0,0,6,-4/)
    wtri(:,25) = (/0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,-3,6,-3,0,0,0,0,0,9,-18,9,0,-6,12,-6,0,2,-4,2,0,0,0,0,0,-6,12,-6,&
         0,4,-8,4/)
    wtri(:,26) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,3,-6,3,0,0,0,0,0,-9,18,-9,0,6,-12,6,0,-2,4,-2,0,0,0,0,0,6,-12,6,0,-4,8,-4/)
    wtri(:,27) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,9,-18,9,0,-6,12,-6,0,0,0,0,0,0,0,0,0,-6,12,-6,0,4,-8,4/)
    wtri(:,28) = (/0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,-9,18,-9,0,6,-12,6,0,0,0,0,0,0,0,0,0,6,-12,6,0,-4,8,-4/)
    wtri(:,29) = (/0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,3,-3,0,0,0,0,0,0,-9,9,0,0,6,-6,0,0,-2,2,0,0,0,0,0,0,6,-6,0,0,-4,4/)
    wtri(:,30) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,-3,3,0,0,0,0,0,0,9,-9,0,0,-6,6,0,0,2,-2,0,0,0,0,0,0,-6,6,0,0,4,-4/)
    wtri(:,31) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,-9,9,0,0,6,-6,0,0,0,0,0,0,0,0,0,0,6,-6,0,0,-4,4/)
    wtri(:,32) = (/0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,9,-9,0,0,-6,6,0,0,0,0,0,0,0,0,0,0,-6,6,0,0,4,-4/)
    wtri(:,33) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,&
         2,0,0,0,0,-2,0,6,-4,4,0,-12,8,-2,0,6,-4,0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2/)
    wtri(:,34) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,-1,0,3,-2,2,0,-6,4,-1,0,3,-2,0,0,0,0,1,0,-3,2,-2,0,6,-4,1,0,-3,2/)
    wtri(:,35) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,1,0,-3,2,-1,0,3,-2,0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2/)
    wtri(:,36) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2,0,&
         0,0,0,0,0,0,0,2,0,-6,4,-2,0,6,-4,0,0,0,0,0,0,0,0,-1,0,3,-2,1,0,-3,2/)
    wtri(:,37) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2,0,&
         0,0,0,0,0,-6,4,0,0,12,-8,0,0,-6,4,0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2/)
    wtri(:,38) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,-3,2,0,0,6,-4,0,0,-3,2,0,0,0,0,0,0,3,-2,0,0,-6,4,0,0,3,-2/)
    wtri(:,39) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,3,-2,0,0,-3,2,0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2/)
    wtri(:,40) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,&
         -2,0,0,0,0,0,0,0,0,0,0,6,-4,0,0,-6,4,0,0,0,0,0,0,0,0,0,0,-3,2,0,0,3,-2/)
    wtri(:,41) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,&
         -4,2,0,-2,4,-2,0,0,0,0,0,6,-12,6,0,-4,8,-4,0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2/)
    wtri(:,42) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         -1,2,-1,0,0,0,0,0,3,-6,3,0,-2,4,-2,0,1,-2,1,0,0,0,0,0,-3,6,-3,0,2,-4,2/)
    wtri(:,43) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,-3,6,-3,0,2,-4,2,0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2/)
    wtri(:,44) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,&
         -2,0,0,0,0,0,0,0,0,0,-6,12,-6,0,4,-8,4,0,0,0,0,0,0,0,0,0,3,-6,3,0,-2,4,-2/)
    wtri(:,45) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,&
         2,0,0,2,-2,0,0,0,0,0,0,-6,6,0,0,4,-4,0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2/)
    wtri(:,46) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,1,-1,0,0,0,0,0,0,-3,3,0,0,2,-2,0,0,-1,1,0,0,0,0,0,0,3,-3,0,0,-2,2/)
    wtri(:,47) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,3,-3,0,0,-2,2,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2/)
    wtri(:,48) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,&
         2,-2,0,0,0,0,0,0,0,0,0,0,6,-6,0,0,-4,4,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,2,-2/)
    wtri(:,49) = (/0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,-3,6,-3,0,6,-12,6,0,-3,6,-3,0,0,0,0,0,2,-4,2,0,-4,8,-4,0,2,-4,2/)
    wtri(:,50) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,3,-6,3,0,-6,12,-6,0,3,-6,3,0,0,0,0,0,-2,4,-2,0,4,-8,4,0,-2,4,-2/)
    wtri(:,51) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,-3,6,-3,0,3,-6,3,0,0,0,0,0,0,0,0,0,2,-4,2,0,-2,4,-2/)
    wtri(:,52) = (/0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,3,-6,3,0,-3,6,-3,0,0,0,0,0,0,0,0,0,-2,4,-2,0,2,-4,2/)
    wtri(:,53) = (/0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,3,-3,0,0,-6,6,0,0,3,-3,0,0,0,0,0,0,-2,2,0,0,4,-4,0,0,-2,2/)
    wtri(:,54) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,-3,3,0,0,6,-6,0,0,-3,3,0,0,0,0,0,0,2,-2,0,0,-4,4,0,0,2,-2/)
    wtri(:,55) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,3,-3,0,0,-3,3,0,0,0,0,0,0,0,0,0,0,-2,2,0,0,2,-2/)
    wtri(:,56) = (/0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,0,-3,3,0,0,3,-3,0,0,0,0,0,0,0,0,0,0,2,-2,0,0,-2,2/)
    wtri(:,57) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,&
         -2,1,0,0,0,0,0,-2,4,-2,0,4,-8,4,0,-2,4,-2,0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1/)
    wtri(:,58) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,-1,2,-1,0,2,-4,2,0,-1,2,-1,0,0,0,0,0,1,-2,1,0,-2,4,-2,0,1,-2,1/)
    wtri(:,59) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,1,-2,1,0,-1,2,-1,0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1/)
    wtri(:,60) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,&
         -2,1,0,0,0,0,0,0,0,0,0,2,-4,2,0,-2,4,-2,0,0,0,0,0,0,0,0,0,-1,2,-1,0,1,-2,1/)
    wtri(:,61) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,&
         -1,1,0,0,0,0,0,0,2,-2,0,0,-4,4,0,0,2,-2,0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1/)
    wtri(:,62) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,1,-1,0,0,-2,2,0,0,1,-1,0,0,0,0,0,0,-1,1,0,0,2,-2,0,0,-1,1/)
    wtri(:,63) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,&
         0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,1,-1,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1/)
    wtri(:,64) = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,&
         0,-1,1,0,0,0,0,0,0,0,0,0,0,-2,2,0,0,2,-2,0,0,0,0,0,0,0,0,0,0,1,-1,0,0,-1,1/)
    
  END SUBROUTINE initialize_tricubic_matrix
  
  !*************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE get_bicubic_coeffs
  !
  !> \brief Calculate coefficients used in bicubic interpolation,
  !> \details Given arrays y, y1, y2, and y12, each of length 4,
  !> containing the function, gradients, and cross derivative at the
  !> four grid points of a rectangular grid cell (numbered counterclockwise
  !> from the lower left), and given d1 and d2, the length of the grid cell
  !> in the 1- and 2- directions, this routine returns the 4 × 4 table c
  !> that is used by routine bcuint for bicubic interpolation,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE get_bicubic_coeffs(y,y1,y2,y12,d1,d2,wt,c)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)  :: d1,d2
    REAL(dp), INTENT(IN)  :: y(4)
    REAL(dp), INTENT(IN)  :: y1(4),y2(4)
    REAL(dp), INTENT(IN)  :: y12(4)
    REAL(dp), INTENT(IN)  :: wt(16,16)
    REAL(dp), INTENT(OUT) :: c(4,4)
    
    REAL(dp)  :: x(16)
    
    x(1:4)=y
    x(5:8)=y1*d1
    x(9:12)=y2*d2
    x(13:16)=y12*d1*d2
    x=matmul(wt,x)
    c=reshape(x,(/4,4/),order=(/2,1/))
    
  END SUBROUTINE get_bicubic_coeffs
  
  !********************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE get_tricubic_coeffs
  !
  !> \brief Calculate coefficients used in bicubic interpolation,
  !> \details Given arrays y, y1, y2, y3 and y12, y13, y23 and y123
  !> each of length 8, containing the function, gradients, and cross derivative
  !> at the eight grid points of a cubic grid cell (numbered counterclockwise
  !> from the lower left), and given d1, d2 and d3, the length of the grid cell
  !> in the 1-, 2- and 3-directions, this routine returns the 8 × 8 table c
  !> that is used by routine bcuint for bicubic interpolation,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y3  Value of the derivative of the function in the third
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] y13 Value of the cross-derivative at known points
  !> \param[in] y23 Value of the cross-derivative at known points
  !> \param[in] y123 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE get_tricubic_coeffs(y,y1,y2,y3,y12,y13,y23,y123,d1,d2,d3,wt,c)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)  :: d1,d2,d3
    REAL(dp), INTENT(IN)  :: y(8)
    REAL(dp), INTENT(IN)  :: y1(8),y2(8),y3(8)
    REAL(dp), INTENT(IN)  :: y12(8),y13(8),y23(8)
    REAL(dp), INTENT(IN)  :: y123(8)
    REAL(dp), INTENT(IN)  :: wt(64,64)
    REAL(dp), INTENT(OUT) :: c(4,4,4)
    
    REAL(dp)  :: x(64)
    
    !----------------------------------------------!
    
    x(1:8)   = y
    x(9:16)  = y1 * d1
    x(17:24) = y2 * d2
    x(25:32) = y3 * d3
    x(33:40) = y12 * d1 * d2
    x(41:48) = y13 * d1 * d3
    x(49:56) = y23 * d2 * d3
    x(57:64) = y123 * d1 * d2 * d3
    
    x=matmul(wt,x)
    c=reshape(x,(/4,4,4/),order=(/3,2,1/))
    
  END SUBROUTINE get_tricubic_coeffs
  
  !********************************************************!
  
  !---------------------------------------------------------------------------
  !
  !  SUBROUTINE bicubic_interpolation
  !
  !> \brief Get interpolated values from bicubic interpolation
  !> \details  Bicubic interpolation within a grid square,
  !> Input quantities are y,y1,y2,y12, x1l and x1u, the lower and upper
  !> coordinates of the grid square in the 1- direction; x2l and x2u
  !> likewise for the 2-direction; and x1,x2, the coordinates of the
  !> desired point for the interpolation, The interpolated function value
  !> is returned as ansy, and the interpolated gradient values as
  !> ansy1 and ansy2,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE bicubic_interpolation(y,y1,y2,y12,x1l,x1u,x2l,x2u,&
       x1,x2,ansy,ansy1,ansy2)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: y(4)
    REAL(dp), INTENT(IN) :: y1(4),y2(4),y12(4)
    REAL(dp), INTENT(IN) :: x1l,x1u,x2l,x2u
    REAL(dp), INTENT(IN) :: x1,x2
    REAL(dp), INTENT(OUT):: ansy,ansy1,ansy2
    
    REAL(dp) :: t,u
    REAL(dp) :: d1, d2
    REAL(dp) :: c(4,4)
    INTEGER  :: i
    
    !---------------------------------------------------------------!
    
    d1 = x1u - x1l
    d2 = x2u - x2l
    
    IF ( (d1.EQ.0.0_dp).OR.(d2.EQ.0.0_dp)) THEN
       WRITE(*,*) 'Problem with input values in bicubic interpolate,'
       WRITE(*,*) 'Are the boundary pair equal?'
       STOP
    ENDIF
    
    CALL get_bicubic_coeffs(y,y1,y2,y12,d1,d2,wbi,c)
    
    t = (x1-x1l) / d1
    u = (x2-x2l) / d2
    
    ansy = 0.0_dp
    ansy2 = 0.0_dp
    ansy1 = 0.0_dp
    
    DO i = 4, 1, -1
       ansy = t * ansy + ((c(i,4) * u + c(i,3)) * u + c(i,2)) * u + c(i,1)
       ansy2 = t * ansy2 + (3.0_sp * c(i,4) * u + 2.0_dp * c(i,3)) * u + c(i,2)
       ansy1 = u * ansy1 + (3.0_sp * c(4,i) * t + 2.0_dp * c(3,i)) * t + c(2,i)
    ENDDO
    
    ansy1 = ansy1 / d1
    ansy2 = ansy2 / d2
    
  END SUBROUTINE bicubic_interpolation
  
  !**********************************************************************!
  
  !---------------------------------------------------------------------------
  !
  !  SUBROUTINE tricubic_interpolation
  !
  !> \brief Get interpolated values from tricubic interpolation
  !> \details  Tricubic interpolation within a grid cube,
  !> Input quantities are y,y1,y2,y3,y12,y13,y23,y123,
  !> x1l and x1u, the lower and upper coordinates of the grid square
  !> in the 1- direction; x2l, x2u and x3l and x3u likewise for the
  !> 2 and 3-direction; and x1,x2,x3 the coordinates of the
  !> desired point for the interpolation, The interpolated function value
  !> is returned as ansy, and the interpolated gradient values as
  !> ansy1 and ansy2,
  !
  !> \param[in] y Function value at known points
  !> \param[in] y1 Value of the derivative of the function in the first
  !>  direction at known points
  !> \param[in] y2  Value of the derivative of the function in the second
  !>  direction at known points
  !> \param[in] y3  Value of the derivative of the function in the thrird
  !>  direction at known points
  !> \param[in] y12 Value of the cross-derivative at known points x and y
  !> \param[in] y13 Value of the cross-derivative at known points x and z
  !> \param[in] y23 Value of the cross-derivative at known points y and z
  !> \param[in] y123 Value of the cross-derivative at known points x, y and z
  !> \param[in] wt Matrix with tabulated coefficients from linear transformation
  !> \return c Matrix with the desired coefficients
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE tricubic_interpolation(y,y1,y2,y3,y12,y13,y23,y123,&
       x1l,x1u,x2l,x2u,x3l,x3u,&
       x1,x2,x3,ansy,ansy1,ansy2,ansy3)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: y(8)
    REAL(dp), INTENT(IN) :: y1(8),y2(8),y3(8)
    REAL(dp), INTENT(IN) :: y12(8), y13(8), y23(8)
    REAL(dp), INTENT(IN) :: y123(8)
    REAL(dp), INTENT(IN) :: x1l,x1u,x2l,x2u,x3l,x3u
    REAL(dp), INTENT(IN) :: x1,x2,x3
    REAL(dp), INTENT(OUT):: ansy,ansy1,ansy2,ansy3
    
    REAL(dp) :: t, u, v
    REAL(dp) :: d1, d2, d3
    REAL(dp) :: c(4,4,4)
    INTEGER  :: i, j ,k
    
    !---------------------------------------------------------------!
    
    d1 = x1u - x1l
    d2 = x2u - x2l
    d3 = x3u - x3l
    
    IF ( (d1.EQ.0.0_dp).OR.(d2.EQ.0.0_dp).OR.(d3.EQ.0.0_dp)) THEN
       WRITE(*,*) 'Problem with input values in bicubic interpolate,'
       WRITE(*,*) 'Are the boundary pair equal?'
       STOP
    ENDIF
    
    CALL get_tricubic_coeffs(y,y1,y2,y3,y12,y13,y23,y123,&
         d1,d2,d3,wtri,c)
    
    t = (x1-x1l) / d1
    u = (x2-x2l) / d2
    v = (x3-x3l) / d3
    
    ansy  = 0.0_dp
    ansy3 = 0.0_dp
    ansy2 = 0.0_dp
    ansy1 = 0.0_dp
    
    DO i = 1, 4
       DO j = 1, 4
          DO k = 1, 4
             ansy  = ansy  + c(i,j,k) * t**(i-1) * u**(j-1) * v**(k-1)
             ansy1 = ansy1 + c(i,j,k) * t**(i-2) * u**(j-1) * v**(k-1)
             ansy2 = ansy2 + c(i,j,k) * t**(i-1) * u**(j-2) * v**(k-1)
             ansy3 = ansy3 + c(i,j,k) * t**(i-1) * u**(j-1) * v**(k-2)
          ENDDO
       ENDDO
    ENDDO
    
    ansy1 = ansy1 / d1
    ansy2 = ansy2 / d2
    ansy3 = ansy3 / d3
    
  END SUBROUTINE tricubic_interpolation
  
  !**********************************************************************!

  SUBROUTINE delete_bicubic_matrix()
    
    IMPLICIT NONE
    
    DEALLOCATE(wbi)
    
  END SUBROUTINE delete_bicubic_matrix
  
  !********************************************************************!
  
  SUBROUTINE delete_tricubic_matrix
    
    IMPLICIT NONE
    
    DEALLOCATE(wtri)
    
  END SUBROUTINE delete_tricubic_matrix
  
  !********************************************************************!
  !*********************************************!
END MODULE cubic_interp

