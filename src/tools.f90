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
  PUBLIC        :: make_derivative
  PUBLIC        :: make_cross_derivative
  PUBLIC        :: fdweights
  PUBLIC        :: interv
  PUBLIC        :: prime_sieve
  PUBLIC        :: trial_division
  PUBLIC        :: reshape_numptsperproc

  ! Interfaces
  INTERFACE make_wave_boundary_derivative
     MODULE PROCEDURE make_wave_boundary_derivative2D
     MODULE PROCEDURE make_wave_boundary_derivative3D
  END INTERFACE make_wave_boundary_derivative

  INTERFACE make_derivative
     MODULE PROCEDURE dmake_derivative
     MODULE PROCEDURE zmake_derivative
  END INTERFACE make_derivative

  INTERFACE make_cross_derivative
     MODULE PROCEDURE dmake_cross_derivative2D
     MODULE PROCEDURE zmake_cross_derivative2D
     MODULE PROCEDURE dmake_cross_derivative3D
     MODULE PROCEDURE zmake_cross_derivative3D
  END INTERFACE make_cross_derivative
  
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
  
  SUBROUTINE initialize_fd_coeffs(fd_rule)

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: fd_rule
    
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
    
  END SUBROUTINE initialize_fd_coeffs
  
  !********************************************************!
  
  SUBROUTINE delete_fd_coeffs()
    IMPLICIT NONE

    DEALLOCATE(fd_coeff1, fd_coeff2)

  END SUBROUTINE delete_fd_coeffs
  
  !********************************************************!
  SUBROUTINE make_wave_boundary_derivative2D(wave,wave_deriv, &
       fd_rule, deltar, numrpts, numthetapts)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)     :: wave(:, :)
    COMPLEX(dp), INTENT(OUT)    :: wave_deriv(:)
    INTEGER, INTENT(IN)         :: fd_rule
    REAL(dp), INTENT(IN)        :: deltar
    INTEGER, INTENT(IN)         :: numrpts
    INTEGER, INTENT(IN)         :: numthetapts
    
    COMPLEX(dp)                 :: sum
    INTEGER                     :: ifd, itheta
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
               fd_coeff1(ifd) / deltar * &
               wave(fd_rule+1+ifd,itheta)
       ENDDO
       wave_deriv(itheta) = sum
    ENDDO
    
    
  END SUBROUTINE make_wave_boundary_derivative2D

  !**************************************************!

    SUBROUTINE make_wave_boundary_derivative3D(wave,wave_deriv, &
       fd_rule, deltar, numrpts, numthetapts, numphipts)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)     :: wave(:, :, :)
    COMPLEX(dp), INTENT(OUT)    :: wave_deriv(:, :)
    INTEGER, INTENT(IN)         :: fd_rule
    REAL(dp), INTENT(IN)        :: deltar
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
                  fd_coeff1(ifd) / deltar * &
                  wave(fd_rule+1+ifd,itheta,iphi) 
          ENDDO
          wave_deriv(itheta,iphi) = sum
       ENDDO
    ENDDO
    
  END SUBROUTINE make_wave_boundary_derivative3D

  !*******************************************************!
  
  SUBROUTINE dmake_derivative(func,func_dx,&
       fd_rule, dx, fdcoeffs )
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule    
    REAL(dp), INTENT(IN)           :: func(-fd_rule:)
    REAL(dp), INTENT(OUT)          :: func_dx
    REAL(dp), INTENT(IN)           :: dx
    REAL(dp), INTENT(IN), OPTIONAL :: fdcoeffs(-fd_rule:fd_rule)
    
    INTEGER                        :: ifd
    
    !------------------------------------------------!
    
    func_dx = 0.0_dp
    
    IF(PRESENT(fdcoeffs)) THEN
       DO ifd = -fd_rule, fd_rule
          func_dx = func_dx + &
               fdcoeffs(ifd) * func(ifd)
       ENDDO
    ELSE
       DO ifd = -fd_rule, fd_rule
          func_dx = func_dx + &
               fd_coeff1(ifd) / dx * func(ifd)
       ENDDO
    ENDIF
    
  END SUBROUTINE dmake_derivative

   !*******************************************************!
  
  SUBROUTINE zmake_derivative(func,func_dx,&
       fd_rule, dx, fdcoeffs)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule  
    COMPLEX(dp), INTENT(IN)        :: func(-fd_rule:)
    COMPLEX(dp), INTENT(OUT)       :: func_dx
    REAL(dp), INTENT(IN)           :: dx
    REAL(dp), INTENT(IN), OPTIONAL :: fdcoeffs(-fd_rule:fd_rule)

    INTEGER                        :: ifd
    
    !------------------------------------------------!
    
    func_dx = ZERO
    
    IF(PRESENT(fdcoeffs)) THEN
       DO ifd = -fd_rule, fd_rule
          func_dx = func_dx + &
               fdcoeffs(ifd) * func(ifd)
          
       ENDDO
    ELSE
       DO ifd = -fd_rule, fd_rule
          func_dx = func_dx + &
               fd_coeff1(ifd) / dx * func(ifd)
       ENDDO
    ENDIF
    
  END SUBROUTINE zmake_derivative
  
  !*******************************************************!
  
  SUBROUTINE dmake_cross_derivative2D(func,func_dxdy,&
       fd_rule, dx, dy, xcoeffs, ycoeffs)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule
    REAL(dp), INTENT(IN)           :: func(-fd_rule:,&
         -fd_rule:)
    REAL(dp), INTENT(OUT)          :: func_dxdy
    REAL(dp), INTENT(IN)           :: dx, dy
    REAL(dp), INTENT(IN), OPTIONAL :: xcoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: ycoeffs(-fd_rule:fd_rule)
    
    REAL(dp)                   :: func_dx(-fd_rule:fd_rule)
    INTEGER                    :: ifdx, ifdy
    
    !------------------------------------------------!

    func_dx   = 0.0_dp
    func_dxdy = 0.0_dp

    IF(PRESENT(xcoeffs)) THEN
       DO ifdy = -fd_rule, fd_rule
          DO ifdx = -fd_rule, fd_rule
             func_dx(ifdy) = func_dx(ifdy) + &
                  xcoeffs(ifdx) * func(ifdx,ifdy)
          ENDDO
       ENDDO
    ELSE
       DO ifdy = -fd_rule, fd_rule
          DO ifdx = -fd_rule, fd_rule
             func_dx(ifdy) = func_dx(ifdy) + &
                  fd_coeff1(ifdx) / dx * func(ifdx,ifdy)
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(ycoeffs)) THEN
       DO ifdy = -fd_rule, fd_rule
          func_dxdy = func_dxdy + &
               ycoeffs(ifdy) * func_dx(ifdy)
       ENDDO
    ELSE
       DO ifdy = -fd_rule, fd_rule
          func_dxdy = func_dxdy + &
               fd_coeff1(ifdy) / dy * func_dx(ifdy)
       ENDDO
    ENDIF
    
  END SUBROUTINE dmake_cross_derivative2D

  !*******************************************************!
  
  SUBROUTINE zmake_cross_derivative2D(func,func_dxdy,&
       fd_rule, dx, dy, xcoeffs, ycoeffs)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule
    COMPLEX(dp), INTENT(IN)        :: func(-fd_rule:,&
         -fd_rule:)
    COMPLEX(dp), INTENT(OUT)       :: func_dxdy
    REAL(dp), INTENT(IN)           :: dx, dy
    REAL(dp), INTENT(IN), OPTIONAL :: xcoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: ycoeffs(-fd_rule:fd_rule)
    
    COMPLEX(dp)                    :: func_dx(-fd_rule:fd_rule)
    INTEGER                        :: ifdx, ifdy
    
    !------------------------------------------------!

    func_dx   = ZERO
    func_dxdy = ZERO
    
    IF(PRESENT(xcoeffs)) THEN
       DO ifdy = -fd_rule, fd_rule
          DO ifdx = -fd_rule, fd_rule
             func_dx(ifdy) = func_dx(ifdy) + &
                  xcoeffs(ifdx) * func(ifdx,ifdy)
          ENDDO
       ENDDO
    ELSE
       DO ifdy = -fd_rule, fd_rule
          DO ifdx = -fd_rule, fd_rule
             func_dx(ifdy) = func_dx(ifdy) + &
                  fd_coeff1(ifdx) / dx * func(ifdx,ifdy)
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(ycoeffs)) THEN
       DO ifdy = -fd_rule, fd_rule
          func_dxdy = func_dxdy + &
               ycoeffs(ifdy) * func_dx(ifdy)
       ENDDO
    ELSE
       DO ifdy = -fd_rule, fd_rule
          func_dxdy = func_dxdy + &
               fd_coeff1(ifdy) / dy * func_dx(ifdy)
       ENDDO
    ENDIF
    
  END SUBROUTINE zmake_cross_derivative2D

    !*******************************************************!
  
  SUBROUTINE dmake_cross_derivative3D(func,func_dxdydz,&
       fd_rule, dx, dy, dz, xcoeffs, ycoeffs, zcoeffs )
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule
    REAL(dp), INTENT(IN)           :: func(-fd_rule:,&
         -fd_rule:,-fd_rule:)
    REAL(dp), INTENT(OUT)          :: func_dxdydz
    REAL(dp), INTENT(IN)           :: dx, dy, dz
    REAL(dp), INTENT(IN), OPTIONAL :: xcoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: ycoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: zcoeffs(-fd_rule:fd_rule)
    
    REAL(dp)                       :: func_dx(-fd_rule:fd_rule,&
         -fd_rule:fd_rule)
    REAL(dp)                       :: func_dxdy(-fd_rule:fd_rule)
    INTEGER                        :: ifdx, ifdy, ifdz
    
    !------------------------------------------------!
    
    func_dx     = 0.0_dp
    func_dxdy   = 0.0_dp
    func_dxdydz = 0.0_dp
    
    IF(PRESENT(xcoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             DO ifdx = -fd_rule, fd_rule
                func_dx(ifdy,ifdz) = func_dx(ifdy,ifdz) + &
                     xcoeffs(ifdx) * func(ifdx,ifdy,ifdz)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             DO ifdx = -fd_rule, fd_rule
                func_dx(ifdy,ifdz) = func_dx(ifdy,ifdz) + &
                     fd_coeff1(ifdx) / dx * func(ifdx,ifdy,ifdz)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(ycoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             func_dxdy(ifdz) = func_dxdy(ifdz) + &
                  ycoeffs(ifdy) * func_dx(ifdy,ifdz)
          ENDDO
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             func_dxdy(ifdz) = func_dxdy(ifdz) + &
                  fd_coeff1(ifdy) / dy * func_dx(ifdy,ifdz)
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(zcoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          func_dxdydz = func_dxdydz + &
               zcoeffs(ifdz) * func_dxdy(ifdz)
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          func_dxdydz = func_dxdydz + &
               fd_coeff1(ifdz) / dz * func_dxdy(ifdz)
       ENDDO
    ENDIF
        
  END SUBROUTINE dmake_cross_derivative3D

  !*********************************************************!
  
  SUBROUTINE zmake_cross_derivative3D(func,func_dxdydz,&
       fd_rule, dx, dy, dz, xcoeffs, ycoeffs, zcoeffs )
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)            :: fd_rule
    COMPLEX(dp), INTENT(IN)        :: func(-fd_rule:,&
         -fd_rule:,-fd_rule:)
    COMPLEX(dp), INTENT(OUT)       :: func_dxdydz
    REAL(dp), INTENT(IN)           :: dx, dy, dz
    REAL(dp), INTENT(IN), OPTIONAL :: xcoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: ycoeffs(-fd_rule:fd_rule)
    REAL(dp), INTENT(IN), OPTIONAL :: zcoeffs(-fd_rule:fd_rule)
    
    COMPLEX(dp)                    :: func_dx(-fd_rule:fd_rule,&
         -fd_rule:fd_rule)
    COMPLEX(dp)                    :: func_dxdy(-fd_rule:fd_rule)
    INTEGER                        :: ifdx, ifdy, ifdz
    
    !------------------------------------------------!
    
    func_dx     = ZERO
    func_dxdy   = ZERO
    func_dxdydz = ZERO
    
    IF(PRESENT(xcoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             DO ifdx = -fd_rule, fd_rule
                func_dx(ifdy,ifdz) = func_dx(ifdy,ifdz) + &
                     xcoeffs(ifdx) * func(ifdx,ifdy,ifdz)
             ENDDO
          ENDDO
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             DO ifdx = -fd_rule, fd_rule
                func_dx(ifdy,ifdz) = func_dx(ifdy,ifdz) + &
                     fd_coeff1(ifdx) / dx * func(ifdx,ifdy,ifdz)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(ycoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             func_dxdy(ifdz) = func_dxdy(ifdz) + &
                  ycoeffs(ifdy) * func_dx(ifdy,ifdz)
          ENDDO
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          DO ifdy = -fd_rule, fd_rule
             func_dxdy(ifdz) = func_dxdy(ifdz) + &
                  fd_coeff1(ifdy) / dy * func_dx(ifdy,ifdz)
          ENDDO
       ENDDO
    ENDIF

    IF(PRESENT(zcoeffs)) THEN
       DO ifdz = -fd_rule, fd_rule
          func_dxdydz = func_dxdydz + &
               zcoeffs(ifdz) * func_dxdy(ifdz)
       ENDDO
    ELSE
       DO ifdz = -fd_rule, fd_rule
          func_dxdydz = func_dxdydz + &
               fd_coeff1(ifdz) / dz * func_dxdy(ifdz)
       ENDDO
    ENDIF
        
  END SUBROUTINE zmake_cross_derivative3D
  
  !*******************************************************!
  !*******************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE fdweights
  !
  !> \brief Create finite-difference weights.
  !> \details Create coefficients to be used for
  !> finite-difference rules.
  !
  !> \param[in] xi Central grid point.
  !> \param[in] x Grid points.
  !> \param[in] n Order of the finite-difference rule.
  !> \param[in] m Order of the highest derivative + 1 (This includes 0th).
  !> \param[out] Finite-difference coefficients.
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE fdweights(xi, x, n, m, c)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)     :: xi
    REAL(dp), INTENT(IN)     :: x(:)
    INTEGER, INTENT(IN)      :: n
    INTEGER, INTENT(IN)      :: m
    REAL(dp), INTENT(OUT)    :: c(:,:)
    
    REAL(dp)                 :: c1,c2,c3,c4,c5
    INTEGER                  :: mn
    INTEGER                  :: i,j,k

    !---------------------------------------------------!
    
    c1 = 1.0_dp
    c4 = x(1) - xi
    
    DO k = 1, m
       DO j = 1, n
          c(j,k) = 0.0_dp
       ENDDO
    ENDDO
    
    c(1,1) = 1.0_dp
    DO i = 2, n
       mn = MIN(i,m)
       c2 = 1.0_dp
       c5 = c4
       c4 = x(i) - xi
       DO j = 1, i-1
          c3 = x(i) - x(j)
          c2 = c2 * c3
          DO k = mn, 2, -1
             c(i,k) = c1 * ((k-1) * c(i-1,k-1) - c5 * c(i-1,k)) / c2
          ENDDO
          c(i,1) = -c1 * c5 * c(i-1,1) / c2
          DO k = mn, 2, -1
             c(j,k) = (c4 * c(j,k) - (k-1) * c(j,k-1)) / c3
          ENDDO
          c(j,1) = c4 * c(j,1) / c3
       ENDDO
       c1 = c2
    ENDDO
    
  END SUBROUTINE fdweights
  
  !********************************************************!
  
  !----------------------------------------------------------------------------
  !
  !  SUBROUTINE interv
  !
  !> \brief FInd the interval within T containing x.
  !> \details  The program is designed to be efficient in the common
  !> situation that it is called repeatedly, with  x  taken from an
  !> increasing or decreasing sequence. This will happen, e.g., when
  !> a pp function is to be graphed. The first guess for  left  is
  !> therefore taken to be the value returned at the previous call
  !> and stored in the
  !> local variable  ilo . A first check ascertains that  ilo .lt. lxt
  !> (this is necessary since the present call may have nothing to do
  !> with the previous call). Then, if  xt(ilo) .le. x .lt. xt(ilo+1),
  !> we set  left = ilo  and are done after just three comparisons.
  !> Otherwise, we repeatedly double the difference  istep = ihi - ilo
  !> while also moving  ilo  and  ihi  in the direction of  x , until
  !> xt(ilo) .le. x .lt. xt(ihi), after which we use bisection to get,
  !> in addition, ilo+1 = ihi .left = ilo  is then returned.
  !>                   ******  o u t p u t  ******
  !> left, mflag.....both integers, whose value is
  !>
  !>   1     -1      if               x .lt.  xt(1)
  !>   i      0      if   xt(i)  .le. x .lt. xt(i+1)
  !> i      0      if   xt(i)  .lt. x .eq. xt(i+1) .eq. xt(lxt)
  !> i      1      if   xt(i)  .lt.        xt(i+1) .eq. xt(lxt) .lt. x
  !>
  !> In particular,  mflag = 0  is the 'usual' case.  mflag .ne. 0
  !> indicates that  x  lies outside the CLOSED interval
  !> xt(1) .le. y .le. xt(lxt) . The asymmetric treatment of the
  !> intervals is due to the decision to make all pp functions cont-
  !> inuous from the right, but, by returning  mflag = 0  even if
  !> x = xt(lxt), there is the option of having the computed pp function
  !> continuous from the left at  xt(lxt) .
  !
  !> \param[in] xt .a real sequence, of length  lxt ,
  !>   assumed to be nondecreasing.
  !> \param[in] lxt .number of terms in the sequence  xt.
  !> \param[in] x the point whose location with respect to the sequence
  !>   xt is to be determined.
  !> \param[out] left Left index on the interval.
  !> \param[out] mflag Control parameter.
  !
  !----------------------------------------------------------------------------
  
  SUBROUTINE interv(xt,lxt,x,left,mflag)
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)      :: lxt
    REAL(dp), INTENT(IN)     :: xt(lxt)
    REAL(dp), INTENT(IN)     :: x
    INTEGER, INTENT(OUT)     :: left, mflag
    
    INTEGER                  :: ihi,ilo,istep
    INTEGER                  :: middle, i, n
    
    !------------------------------------------------!
    
    ilo=1
    ihi=lxt
    left=0
    mflag=0
    
    IF((x.LT.xt(ilo)).OR.(lxt.LE.1)) THEN
       left = 1
       mflag = -1
    ELSE
       IF(x.GT.xt(ihi)) THEN
          left = lxt
          mflag = 1
       ELSE
          N=INT(LOG(REAL(lxt-1))/LOG(2.0_dp)+2.0_dp)
          DO i = 1, n
             middle=(ihi+ilo)/2
             IF(middle.EQ.ilo) left = ilo
             IF(x.LT.xt(middle)) ihi = middle
             IF(x.GE.xt(middle)) ilo = middle
          ENDDO
       ENDIF
    ENDIF
    IF(left.EQ.0) THEN
       WRITE(*,*) 'Possible problem in interv'
       STOP
    ENDIF
    
  END SUBROUTINE interv

  !********************************************************!
  
  SUBROUTINE prime_sieve(last_number, number_of_primes, &
       numbers )
    
    !  Find prime numbers using array processing.
    !  We use the sieve of Eratosthenes to calculate
    !  the first prime numbers.
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)              :: last_number
    INTEGER, INTENT(INOUT)           :: number_of_primes
    INTEGER, INTENT(OUT)             :: numbers(:)
    
    INTEGER                          :: i, ac
  
    !--------------------------------------------------!
    
    ! Check dimensions for the output array
    IF(SIZE(numbers).ne.last_number) THEN
       WRITE(*,*) 'Error in prime sieve. The dimension '&
            &'of output array is not equal to last number.'
       STOP
    ENDIF
    
    !
    !  Initialize numbers array to 0, 2, 3, ..., last_number--
    !  Zero instead of 1 because 1 is a special case.
    numbers = (/ 0, (ac, ac = 2, last_number) /)
    
    DO i = 2, last_number
       
       ! if this number is prime
       ! eliminate all multiples
       
       IF (numbers(i) /= 0) THEN
          numbers(2*i : last_number : i) = 0
       ENDIF
    ENDDO
    
    !  Count the primes.
    number_of_primes = COUNT (numbers /= 0)
    
    !  Gather them into the front of the array.
    numbers(1:number_of_primes) = PACK(numbers, numbers /= 0)
    
    !  Sample output:
    ! There are   25 prime numbers less than   100
    !     2      3      5      7     11
    !    13     17     19     23     29
    !    31     37     41     43     47
    !    53     59     61     67     71
    !    73     79     83     89     97
    
  END SUBROUTINE prime_sieve

  !********************************************************!
  
  ! Return a list of the prime factors for 
  ! a natural number.
  
  SUBROUTINE trial_division(nn, numfactors, &
       prime_factors )
    
    IMPLICIT NONE
    
    integer, intent(in)        :: nn
    integer, intent(inout)     :: numfactors
    integer, intent(out)       :: prime_factors(:)
    
    integer                    :: sqrtn, n, p
    integer, allocatable       :: prime_list(:)
    integer                    :: numprimes
    integer                    :: ip
    
    !------------------------------------------------!
    
    ! First make sure you are not dealing with one.
    IF (nn .LT. 2) return
    
    ! Get list of primes
    n = nn
    sqrtn = INT(SQRT(REAL(n,dp)))
    
    ALLOCATE(prime_list(1:sqrtn))
    
    
    CALL prime_sieve(sqrtn, numprimes, &
         prime_list )
    
    numfactors    = 0
    prime_factors = 0
    
    DO ip = 1, numprimes
       
       p = prime_list(ip)
       
       IF( p*p .GT. n) EXIT
       
       DO
          
          IF( MOD(n,p) .EQ. 0) THEN
             numfactors = numfactors + 1
             prime_factors(numfactors) = p
             n = n / p
          ELSE
             EXIT
          ENDIF
       ENDDO
       
    ENDDO
    
    ! Check if the number from last division
    ! is a prime number
    IF ( n .GT. 1) THEN
       numfactors = numfactors + 1
       prime_factors(numfactors) = n
    ENDIF

  END SUBROUTINE trial_division
  
  !********************************************************!
  
  SUBROUTINE reshape_numptsperproc(numxpts, numprocs, &
       numprocsio, numxptsperproc, numypts, numyptsperproc)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)              :: numxpts
    INTEGER, INTENT(IN)              :: numprocs
    INTEGER, INTENT(INOUT)           :: numprocsio
    INTEGER, INTENT(INOUT)           :: numxptsperproc
    INTEGER, INTENT(IN), OPTIONAL    :: numypts
    INTEGER, INTENT(INOUT), OPTIONAL :: numyptsperproc
    
    INTEGER                          :: numfactx
    INTEGER                          :: numfacty
    INTEGER, ALLOCATABLE             :: factorsx(:)
    INTEGER, ALLOCATABLE             :: factorsy(:)
    INTEGER                          :: aux, ix, iy
    
    !-------------------------------------------------!
    
    ! First, allocate arrays for prime factors of x
    ALLOCATE(factorsx(1:numxpts))
    
    ! Find out the prime factors that decompose both x
    CALL trial_division(numxpts, numfactx, factorsx)

    ix = numfactx + 1
    numxptsperproc = numxpts

    ! We repeat the same if we have a second dimension
    IF(PRESENT(numypts)) THEN
       ALLOCATE(factorsy(1:numypts))
       
       CALL trial_division(numypts, numfacty, factorsy)
       iy = numfacty + 1
       numyptsperproc = numypts
       
    ENDIF
    
    aux = 1
    DO
       ix = ix - 1
       iy = ix - 1
       
       aux = aux * factorsx(ix)
       IF(aux.GT.numprocs) EXIT
       
       numprocsio = aux
       numxptsperproc = numxptsperproc / factorsx(ix)
       
       IF(PRESENT(numypts)) THEN
          aux = aux * factorsy(iy)
          IF(aux.GT.numprocs) EXIT
          
          numprocsio = aux
          numyptsperproc = numyptsperproc / factorsy(iy)
       ENDIF
       
    ENDDO

    ! Deallocate 
    DEALLOCATE(factorsx)
    IF(PRESENT(numypts)) DEALLOCATE(factorsy)
    
  END SUBROUTINE reshape_numptsperproc
  
  !********************************************************!
  !********************************************************!
  
END MODULE tools
