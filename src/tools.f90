MODULE tools

  USE constants
  
  IMPLICIT NONE
  
  PRIVATE

  ! Public routines
  PUBLIC        :: indx
  PUBLIC        :: get_simulation_parameters
  PUBLIC        :: initialize_fd_coeffs
  PUBLIC        :: delete_fd_coeffs
  PUBLIC        :: make_wave_boundary_derivative

  ! Interfaces
  INTERFACE make_wave_boundary_derivative
     MODULE PROCEDURE make_wave_boundary_derivative2D
     MODULE PROCEDURE make_wave_boundary_derivative3D
  END INTERFACE make_wave_boundary_derivative

  ! Private variables
  REAL(dp), ALLOCATABLE        :: fd_coeff1(:)
  REAL(dp), ALLOCATABLE        :: fd_coeff2(:)
  
CONTAINS
  
  FUNCTION indx(ix, iy, iz, ir, Nx, Ny, Nz, Nr)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)      :: ix, iy, iz, ir
    INTEGER, INTENT(IN)      :: Nx, Ny, Nz, Nr
    INTEGER                  :: indx
    
    indx = Nx*Ny*Nz*(ir-1) + Nx*Ny*(iz-1) &
         + Nx*(iy-1) + ix
    
  END FUNCTION indx
  
  !------------------------------------------!
  
  SUBROUTINE get_simulation_parameters( rank, &
       dims_local, dims_global, grid_rank, &
       Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
       ipx, ipy, ipz, ipr, &
       numprocx, numprocy, numprocz, numprocr )
    
    IMPLICIT NONE

    INTEGER, INTENT(IN)     :: rank
    INTEGER, INTENT(IN)     :: dims_local(:)
    INTEGER, INTENT(IN)     :: dims_global(:)
    INTEGER, INTENT(IN)     :: grid_rank(:)
    INTEGER, INTENT(OUT)    :: Nx, Ny, Nz, Nr
    INTEGER, INTENT(OUT)    :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER, INTENT(OUT)    :: ipx, ipy, ipz, ipr
    INTEGER, INTENT(OUT)    :: numprocx, numprocy
    INTEGER, INTENT(OUT)    :: numprocz, numprocr
    
    !-----------------------------------------------!
    
    IF (rank .EQ. 1) THEN
       Nx = dims_local(1)
       Ny = 1
       Nz = 1
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = 1
       Nzgl = 1
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = 0
       ipz = 0
       ipr = 0
       
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = 1
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = 1
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = 0
       ipr = 0
       
    ELSEIF (rank .EQ. 3) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = 1
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = 1
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = grid_rank(3)
       ipr = 0
       
    ELSEIF (rank .EQ. 4) THEN
       
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = dims_local(4)
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = dims_global(4)
       
       ipx = grid_rank(1)
       ipy = grid_rank(2)
       ipz = grid_rank(3)
       ipr = grid_rank(4)
       
    ELSE
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF
    
    ! Work out the number of processors per dimension
    numprocx = Nxgl / Nx
    numprocy = Nygl / Ny
    numprocz = Nzgl / Nz
    numprocr = Nrgl / Nr
    
  END SUBROUTINE get_simulation_parameters

  !****************************************!
  
  SUBROUTINE initialize_fd_coeffs(fd_rule, deltar)

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: fd_rule
    REAL(dp), INTENT(IN)       :: deltar 
    
    ALLOCATE(fd_coeff1(-fd_rule:fd_rule))
    ALLOCATE(fd_coeff2(-fd_rule:fd_rule))
    
    IF (fd_rule.EQ.6) THEN
       
       fd_coeff1(-6)  =      150.0_dp / 831600.0_dp
       fd_coeff1(-5)  = -   2160.0_dp / 831600.0_dp
       fd_coeff1(-4)  =    14850.0_dp / 831600.0_dp
       fd_coeff1(-3)  = -  66000.0_dp / 831600.0_dp
       fd_coeff1(-2)  =   222750.0_dp / 831600.0_dp
       fd_coeff1(-1)  = - 712800.0_dp / 831600.0_dp
       fd_coeff1(0)  =       0.0_dp
       fd_coeff1(1)  =   712800.0_dp / 831600.0_dp
       fd_coeff1(2)  = - 222750.0_dp / 831600.0_dp
       fd_coeff1(3)  =    66000.0_dp / 831600.0_dp
       fd_coeff1(4)  = -  14850.0_dp / 831600.0_dp
       fd_coeff1(5)  =     2160.0_dp / 831600.0_dp
       fd_coeff1(6)  = -    150.0_dp / 831600.0_dp
       
       fd_coeff2(-6)  = -    50.0_dp / 831600.0_dp
       fd_coeff2(-5)  =      864.0_dp / 831600.0_dp
       fd_coeff2(-4)  = -   7425.0_dp / 831600.0_dp
       fd_coeff2(-3)  =    44000.0_dp / 831600.0_dp
       fd_coeff2(-2)  = - 222750.0_dp / 831600.0_dp
       fd_coeff2(-1)  =  1425600.0_dp / 831600.0_dp
       fd_coeff2(0)  = -2480478.0_dp / 831600.0_dp
       fd_coeff2(1)  =  1425600.0_dp / 831600.0_dp
       fd_coeff2(2)  = - 222750.0_dp / 831600.0_dp
       fd_coeff2(3)  =    44000.0_dp / 831600.0_dp
       fd_coeff2(4)  = -   7425.0_dp / 831600.0_dp
       fd_coeff2(5)  =      864.0_dp / 831600.0_dp
       fd_coeff2(6)  = -     50.0_dp / 831600.0_dp
       
    ELSEIF (fd_rule.EQ.5) THEN
       
       fd_coeff1(-5)  = -   20.0_dp / 25200.0_dp
       fd_coeff1(-4)  =    250.0_dp / 25200.0_dp
       fd_coeff1(-3)  = - 1500.0_dp / 25200.0_dp
       fd_coeff1(-2)  =   6000.0_dp / 25200.0_dp
       fd_coeff1(-1)  = -21000.0_dp / 25200.0_dp
       fd_coeff1(0)  =      0.0_dp 
       fd_coeff1(1)  =  21000.0_dp / 25200.0_dp
       fd_coeff1(2)  = - 6000.0_dp / 25200.0_dp
       fd_coeff1(3)  =   1500.0_dp / 25200.0_dp
       fd_coeff1(4)  = -  250.0_dp / 25200.0_dp
       fd_coeff1(5)  =     20.0_dp / 25200.0_dp
       
       fd_coeff2(-5)  =      8.0_dp / 25200.0_dp
       fd_coeff2(-4)  = -  125.0_dp / 25200.0_dp
       fd_coeff2(-3)  =   1000.0_dp / 25200.0_dp
       fd_coeff2(-2)  = - 6000.0_dp / 25200.0_dp
       fd_coeff2(-1)  =  42000.0_dp / 25200.0_dp
       fd_coeff2(0)  = -73766.0_dp / 25200.0_dp
       fd_coeff2(1)  =  42000.0_dp / 25200.0_dp
       fd_coeff2(2)  = - 6000.0_dp / 25200.0_dp
       fd_coeff2(3)  =   1000.0_dp / 25200.0_dp
       fd_coeff2(4)  = -  125.0_dp / 25200.0_dp
       fd_coeff2(5)  =      8.0_dp / 25200.0_dp
       
    ELSEIF (fd_rule.EQ.4) THEN
       
       fd_coeff1(-4)  =     18.0_dp / 5040.0_dp
       fd_coeff1(-3)  = -  192.0_dp / 5040.0_dp
       fd_coeff1(-2)  =   1008.0_dp / 5040.0_dp
       fd_coeff1(-1)  = - 4032.0_dp / 5040.0_dp
       fd_coeff1(0)  =      0.0_dp 
       fd_coeff1(1)  =   4032.0_dp / 5040.0_dp
       fd_coeff1(2)  = - 1008.0_dp / 5040.0_dp
       fd_coeff1(3)  =    192.0_dp / 5040.0_dp
       fd_coeff1(4)  = -   18.0_dp / 5040.0_dp
       
       fd_coeff2(-4)  = -    9.0_dp / 5040.0_dp
       fd_coeff2(-3)  =    128.0_dp / 5040.0_dp
       fd_coeff2(-2)  = - 1008.0_dp / 5040.0_dp
       fd_coeff2(-1)  =   8064.0_dp / 5040.0_dp
       fd_coeff2(0)  = -14350.0_dp / 5040.0_dp
       fd_coeff2(1)  =   8064.0_dp / 5040.0_dp
       fd_coeff2(2)  = - 1008.0_dp / 5040.0_dp
       fd_coeff2(3)  =    128.0_dp / 5040.0_dp
       fd_coeff2(4)  = -    9.0_dp / 5040.0_dp
       
    ELSEIF (fd_rule.EQ.3) THEN
       
       fd_coeff1(-3)  = -  3.0_dp / 180.0_dp
       fd_coeff1(-2)  =   27.0_dp / 180.0_dp
       fd_coeff1(-1)  = -135.0_dp / 180.0_dp
       fd_coeff1(0)  =    0.0_dp 
       fd_coeff1(1)  =  135.0_dp / 180.0_dp
       fd_coeff1(2)  = - 27.0_dp / 180.0_dp
       fd_coeff1(3)  =    3.0_dp / 180.0_dp
       
       fd_coeff2(-3)  =    2.0_dp / 180.0_dp
       fd_coeff2(-2)  = - 27.0_dp / 180.0_dp
       fd_coeff2(-1)  =  270.0_dp / 180.0_dp
       fd_coeff2(0)  = -490.0_dp / 180.0_dp
       fd_coeff2(1)  =  270.0_dp / 180.0_dp
       fd_coeff2(2)  = - 27.0_dp / 180.0_dp
       fd_coeff2(3)  =    2.0_dp / 180.0_dp
       
    ELSEIF (fd_rule.EQ.2) THEN
       
       fd_coeff1(-2)    =  1.0_dp / 12.0_dp
       fd_coeff1(-1)    = -8.0_dp / 12.0_dp
       fd_coeff1(0)	  =  0.0_dp 
       fd_coeff1(1)	  =  8.0_dp / 12.0_dp
       fd_coeff1(2)	  = -1.0_dp / 12.0_dp
       
       fd_coeff2(-2)    = - 1.0_dp / 12.0_dp
       fd_coeff2(-1)    =  16.0_dp / 12.0_dp
       fd_coeff2(0)	  = -30.0_dp / 12.0_dp
       fd_coeff2(1)	  =  16.0_dp / 12.0_dp
       fd_coeff2(2)	  = - 1.0_dp / 12.0_dp
       
    ELSEIF (fd_rule.EQ.1) THEN
       
       fd_coeff1(-1)  = -1.0_dp / 2.0_dp
       fd_coeff1(0)  =  0.0_dp 
       fd_coeff1(1)  =  1.0_dp / 2.0_dp
       
       fd_coeff2(-1)  =  1.0_dp / 1.0_dp
       fd_coeff2(0)  = -2.0_dp / 1.0_dp
       fd_coeff2(1)  =  1.0_dp / 1.0_dp
       
    ELSE
    	WRITE(*,*) 'We do not have the coefficients for &
    	& that rule. Really sorry!'
    	STOP
    ENDIF
    
    fd_coeff1        =  fd_coeff1 / (deltar)
    fd_coeff2        =  fd_coeff2 / (deltar * deltar)
    
  END SUBROUTINE initialize_fd_coeffs
  
  !********************************************************!
  
  SUBROUTINE delete_fd_coeffs()
    IMPLICIT NONE

    DEALLOCATE(fd_coeff1, fd_coeff2)

  END SUBROUTINE delete_fd_coeffs
  
  !********************************************************!
  SUBROUTINE make_wave_boundary_derivative2D(wave,wave_deriv, &
       fd_rule, numrpts, numthetapts)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)     :: wave(:, :)
    COMPLEX(dp), INTENT(OUT)    :: wave_deriv(:)
    INTEGER, INTENT(IN)         :: fd_rule
    INTEGER, INTENT(IN)         :: numrpts
    INTEGER, INTENT(IN)         :: numthetapts
    
    COMPLEX(dp)                 :: sum
    INTEGER                      :: ifd, itheta
    !-------------------------------------------!
    
    IF(numrpts.LE.fd_rule) THEN
       WRITE(*,*) 'There are less points in the function &
            &than points in the rule!'
       STOP
    ENDIF
    
    DO itheta = 1, numthetapts
       sum = ZERO
       DO ifd = -fd_rule , fd_rule
          sum = sum + &
               fd_coeff1(ifd) * wave(fd_rule+1+ifd,itheta)
       ENDDO
       wave_deriv(itheta) = sum
    ENDDO
    
    
  END SUBROUTINE make_wave_boundary_derivative2D

  !**************************************************!

    SUBROUTINE make_wave_boundary_derivative3D(wave,wave_deriv, &
       fd_rule, numrpts, numthetapts, numphipts)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)     :: wave(:, :, :)
    COMPLEX(dp), INTENT(OUT)    :: wave_deriv(:, :)
    INTEGER, INTENT(IN)         :: fd_rule
    INTEGER, INTENT(IN)         :: numrpts
    INTEGER, INTENT(IN)         :: numthetapts
    INTEGER, INTENT(IN)         :: numphipts

    COMPLEX(dp)                 :: sum
    INTEGER                     :: ifd, itheta, iphi
    !-------------------------------------------!
    
    IF(numrpts.LE.fd_rule) THEN
       WRITE(*,*) 'There are less points in the function &
            &than points in the rule!'
       STOP
    ENDIF
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          sum = ZERO
          DO ifd = -fd_rule , fd_rule
             sum = sum + &
                  fd_coeff1(ifd) * wave(fd_rule+1+ifd,itheta,iphi) 
          ENDDO
          wave_deriv(itheta,iphi) = sum
       ENDDO
    ENDDO
    
  END SUBROUTINE make_wave_boundary_derivative3D
  
END MODULE tools
