!-----------------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------------
!                                                                              
! Module Constants                                                             
!
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Constants numbers and data types
!> \details As the name says, this module contains physical constants and other
!> quantities of interest. Also, it constains data types definitions and kinds
! 
!------------------------------------------------------------------------------!

MODULE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  !-----------------------------------
  ! Kinds
  !-----------------------------------
  !> Double precision type length.
  INTEGER, PARAMETER, PUBLIC             :: dp         = KIND(1.0D0)
  !> Single precision type length.
  INTEGER, PARAMETER, PUBLIC             :: sp         = KIND(1.0)
  !> \f$\sqrt{2}\f$
  REAL(dp), PARAMETER, PUBLIC            :: SQRTTWO    = 1.4142135623730951_dp
  ! Complex numbers
  !> Imaginary unit.
  COMPLEX(dp), PARAMETER, PUBLIC         :: ZIMAGONE   = (0.0_dp, 1.0_dp)
  !> Complex number one.
  COMPLEX(dp), PARAMETER, PUBLIC         :: ZONE       = (1.0_dp, 0.0_dp)
  !> Zero, in complex type.
  COMPLEX(dp), PARAMETER, PUBLIC         :: ZERO       = (0.0_dp, 0.0_dp)
  
  !--------------------------------------
  ! Physical and mathematical constants
  !--------------------------------------

  !> \todo document this.
  REAL(dp), PARAMETER, PUBLIC :: amu_to_internal  = 103.6426867_dp
  !> in eV / K units
  REAL(dp), PARAMETER, PUBLIC :: boltzmann_k      = 8.61734215E-5_dp
  !> One eV in hartree (atomic units of energy).
  REAL(dp), PARAMETER, PUBLIC :: ev_to_hartree    = 0.03674932600434263_dp
  !> One hartree (atomic units of energy) in eV.
  REAL(dp), PARAMETER, PUBLIC :: hatree_to_ev     = 27.211383411_dp
  !> \f$\hbar\f$
  REAL(dp), PARAMETER, PUBLIC :: hbar             = 0.65821188926_dp !in eV fs
  !> \f$\pi\f$
  REAL(dp), PARAMETER, PUBLIC :: pi               = 3.141592653589793_dp
  !> \f$2\pi\f$
  REAL(dp), PARAMETER, PUBLIC :: twopi            = 2.0_dp * pi
  !> \f$\frac{1}{\pi}\f$
  REAL(dp), PARAMETER, PUBLIC :: oneovertwopi     = 1.0_dp / twopi
  !> \f$\sqrt{\frac{1}{2\pi}}\f$
  REAL(dp), PARAMETER, PUBLIC :: sqrtoneovertwopi = SQRT(oneovertwopi)
  !> \f$\sqrt{\frac{1}{(2\pi)^3}}\f$
  REAL(dp), PARAMETER, PUBLIC :: oneovertwopi2over3 = SQRT(oneovertwopi) &
       * SQRT(oneovertwopi) * SQRT(oneovertwopi)
  !> One atomic unit of time in femtoseconds
  REAL(dp), PARAMETER, PUBLIC :: autime_to_fs     = 2.418884326505E-2 ! fs
  !> One bohr (one atomic unit of length) in angstroms
  REAL(dp), PARAMETER, PUBLIC :: bohr_to_ang      = 0.5291772083_dp
  !> Unit of current intensity.
  REAL(dp), PARAMETER, PUBLIC :: efs_to_microamps = 160.2176462_dp
  !> Dielectric constant in vacuum
  REAL(dp), PARAMETER, PUBLIC :: epsilon0         = 5.526349954E-3_dp
  !> Conversion factor for converting between 
  !! wavelength and frequency in atomic units. 
  REAL(dp), PARAMETER, PUBLIC :: waveconvert      = 45.7720588235D0
  !> Atomic unit of intensity.
  REAL(dp), PARAMETER, PUBLIC :: intensityau      = 3.51D16		  
  !> Mass factor for the electron to use in kinetic energy. In atomis units
  REAL(dp), PARAMETER         :: massfacelectron  = 0.5D0
  !> Mass factor for the field term. In atomis units
  REAL(dp), PARAMETER         :: massfacfield     = 1.0D0
  
END MODULE constants
