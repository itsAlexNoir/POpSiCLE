!
!
!     .o88b.  .d88b.   .d88b.  d8888b. d8888b. .d8888.
!    d8P  Y8 .8P  Y8. .8P  Y8. 88  `8D 88  `8D 88'  YP
!    8P      88    88 88    88 88oobY' 88   88 `8bo.
!    8b      88    88 88    88 88`8b   88   88   `Y8b.
!    Y8b  d8 `8b  d8' `8b  d8' 88 `88. 88  .8D db   8D
!     `Y88P'  `Y88P'   `Y88P'  88   YD Y8888D' `8888Y'
!
MODULE coords

  USE scattcoords
  USE cubcoords


  IMPLICIT NONE

  PUBLIC           cartesian2spherical
  PUBLIC           cylindrical2spherical
  
  
  INTERFACE cartesian2spherical
     MODULE PROCEDURE dscatt_cartesian2spherical2D
     MODULE PROCEDURE zscatt_cartesian2spherical2D
     MODULE PROCEDURE dscatt_cartesian2spherical3D
     MODULE PROCEDURE zscatt_cartesian2spherical3D
     MODULE PROCEDURE dcubic_cartesian2spherical3D
     MODULE PROCEDURE zcubic_cartesian2spherical3D
  END INTERFACE cartesian2spherical
  
  INTERFACE cylindrical2spherical
     MODULE PROCEDURE dscatt_cylindrical2spherical2D
     MODULE PROCEDURE zscatt_cylindrical2spherical2D
     MODULE PROCEDURE dcubic_cylindrical2spherical2D
     MODULE PROCEDURE zcubic_cylindrical2spherical2D
  END INTERFACE cylindrical2spherical
  
  
END MODULE coords
