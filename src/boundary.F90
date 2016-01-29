MODULE boundary
  
  USE scattboundcyl
  USE cubboundcyl
  USE scattboundcart
  !USE cubboundcart
  
  IMPLICIT NONE
  
  !PRIVATE
  
  PUBLIC           initialize_cylindrical_boundary
  PUBLIC           get_cylindrical_boundary
  PUBLIC           delete_cylindrical_boundary
  PUBLIC           initialize_cartesian_boundary
  PUBLIC           get_cartesian_boundary
  PUBLIC           delete_cartesian_boundary
  
  
  INTERFACE initialize_cylindrical_boundary
     MODULE PROCEDURE initialize_scatt_cylindrical_boundary2D_serial
     MODULE PROCEDURE initialize_bicubic_cylindrical_boundary2D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_scatt_cylindrical_boundary2D_parallel
     MODULE PROCEDURE initialize_bicubic_cylindrical_boundary2D_parallel
#endif
  END INTERFACE initialize_cylindrical_boundary
  
  INTERFACE initialize_cartesian_boundary
     MODULE PROCEDURE initialize_scatt_cartesian_boundary2D
     MODULE PROCEDURE initialize_scatt_cartesian_boundary3D_serial
     !MODULE PROCEDURE initialize_bicubic_cartesian_boundary3D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_scatt_cartesian_boundary3D_parallel
     !MODULE PROCEDURE initialize_bicubic_cartesian_boundary3D_parallel
#endif
  END INTERFACE initialize_cartesian_boundary
  
  INTERFACE get_cylindrical_boundary
     MODULE PROCEDURE dget_scatt_cylindrical_boundary2D
     MODULE PROCEDURE dget_bicubic_cylindrical_boundary2D
     MODULE PROCEDURE zget_scatt_cylindrical_boundary2D
     MODULE PROCEDURE zget_bicubic_cylindrical_boundary2D
  END INTERFACE get_cylindrical_boundary
  
  INTERFACE get_cartesian_boundary
     MODULE PROCEDURE get_scatt_cartesian_boundary2D
     MODULE PROCEDURE dget_scatt_cartesian_boundary3D
     MODULE PROCEDURE zget_scatt_cartesian_boundary3D
     !MODULE PROCEDURE dget_bicubic_cartesian_boundary3D
     !MODULE PROCEDURE zget_bicubic_catesian_boundary3D
  END INTERFACE get_cartesian_boundary
  
  INTERFACE delete_cylindrical_boundary
     MODULE PROCEDURE delete_scatt_cylindrical_boundary2D
     MODULE PROCEDURE delete_bicubic_cylindrical_boundary2D
  END INTERFACE delete_cylindrical_boundary
  
  INTERFACE delete_cartesian_boundary
     MODULE PROCEDURE delete_scatt_cartesian_boundary2D
     MODULE PROCEDURE delete_scatt_cartesian_boundary3D
     !MODULE PROCEDURE delete_bicubic_cylindrical_boundary3D
  END INTERFACE delete_cartesian_boundary
  
END MODULE boundary
