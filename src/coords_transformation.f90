MODULE coords_transform

  USE constants
  USE interp

  IMPLICIT NONE

  PRIVATE
  
  PUBLIC           initialize_cylindrical_boundary
  PUBLIC           get_cylindrical_boundary
  PUBLIC           cartesian2spherical
  PUBLIC           cylindrical2spherical
  
!!$  INTERFACE initialize_cartesian_grid
!!$     MODULE PROCEDURE initialize_cartesian_grid2D
!!$     MODULE PROCEDURE initialize_cartesian_grid3D
!!$  END INTERFACE initialize_cartesian_grid
  
  INTERFACE initialize_cylindrical_boundary
     MODULE PROCEDURE initialize_cylindrical_boundary2D
     !MODULE PROCEDURE initialize_cylindrical_boundary3D
  END INTERFACE initialize_cylindrical_boundary

  INTERFACE get_cylindrical_boundary
     MODULE PROCEDURE get_cylindrical_boundary2D
     !MODULE PROCEDURE get_cylindrical_boundary3D
  END INTERFACE get_cylindrical_boundary
  
  INTERFACE cartesian2spherical
     !MODULE PROCEDURE cartesian2spherical2D
     MODULE PROCEDURE cartesian2spherical3D
  END INTERFACE cartesian2spherical
  
  INTERFACE cylindrical2spherical
     MODULE PROCEDURE cylindrical2spherical2D
     !MODULE PROCEDURE cylindrical2spherical3D
  END INTERFACE cylindrical2spherical
  
  INTEGER                          :: numpts
  INTEGER                          :: numrpts, numthetapts
  REAL(dp), ALLOCATABLE            :: rpts_scatt(:)
  REAL(dp), ALLOCATABLE            :: theta_scatt(:)
  COMPLEX(dp), ALLOCATABLE         :: psi_scatt(:)
  INTEGER, ALLOCATABLE             :: index_x1(:)
  INTEGER, ALLOCATABLE             :: index_x2(:)
  INTEGER, ALLOCATABLE             :: index_x3(:)
  REAL(dp), ALLOCATABLE            :: rpts(:)
  REAL(dp), ALLOCATABLE            :: theta(:)
  COMPLEX(dp), ALLOCATABLE         :: psi_sph(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi_sph_dx(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi_sph_dy(:,:)
  
CONTAINS
  
  !----------------------------------------------------!
  
  SUBROUTINE initialize_cylindrical_boundary2D(rho_ax, z_ax, dims, &
       halolims, Rs, radtol, fdpts, deltar, deltatheta)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: rho_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: halolims(:,:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    REAL(dp), INTENT(IN)      :: deltatheta
    
    REAL(dp)                  :: minrho, maxrho
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    REAL(dp)                  :: rpt, thetapt
    INTEGER                   :: irho, iz, inum
    INTEGER                   :: ir, itheta
    
    !--------------------------------------------------!
    
    numpts = 0
    
    ! Assign max and min values for the axes
    minrho = rho_ax(1)
    maxrho = rho_ax(dims(1))
    
    minz = z_ax(1)
    maxz = z_ax(dims(2))
    
    minr = SQRT(minrho**2 + minz**2)
    maxr = SQRT(maxrho**2 + maxz**2)
    
    mintheta = ATAN2(rho_ax(halolims(1,1)) , z_ax(halolims(2,1)))
    maxtheta = ATAN2(rho_ax(halolims(1,2)) , z_ax(halolims(2,2)))
    
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts(1))
       ALLOCATE(theta(1))
       ALLOCATE(psi_sph(1,1))
       ALLOCATE(psi_sph_dx(1,1))
       ALLOCATE(psi_sph_dy(1,1))
       
       RETURN
    ENDIF
    
    ! Work out how many points will build the interpolant
    DO iz = halolims(2,1), halolims(2,2)
       DO irho = halolims(1,1), halolims(1,2)
          ! Calculate the point in spherical coordinates
          rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
          !thetapt = ATAN2(rho_ax, z_ax)
          
          ! See if this points lies within the radius of
          ! influence of the boundary
          IF (ABS(rpt-Rs) .LE. radtol ) &
               numpts = numpts + 1
       ENDDO
    ENDDO
    
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))  
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts(1))
       ALLOCATE(theta(1))
       ALLOCATE(psi_sph(1,1))
       ALLOCATE(psi_sph_dx(1,1))
       ALLOCATE(psi_sph_dy(1,1))

    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       
       ALLOCATE(rpts(1:numrpts))
       ALLOCATE(theta(1:numthetapts))
       ALLOCATE(psi_sph(1:numrpts,1:numthetapts))
       ALLOCATE(psi_sph_dx(1:numrpts,1:numthetapts))
       ALLOCATE(psi_sph_dy(1:numrpts,1:numthetapts))

       inum = 0
       
       DO iz = halolims(2,1), halolims(2,2)
          DO irho = halolims(1,1), halolims(1,2)
             
             rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
             thetapt = ATAN2(rho_ax(irho), z_ax(iz))
             
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                inum = inum + 1
                rpts_scatt(inum)  = rpt
                theta_scatt(inum) = thetapt
                index_x1(inum) = irho
                index_x2(inum) = iz
             ENDIF
          ENDDO
       ENDDO
       
       DO ir = 1, numrpts
          rpts(ir) = minr + REAL( ir * deltar, dp )
       ENDDO
       
       DO itheta = 1, numthetapts
          theta(itheta) = mintheta + REAL( itheta * deltatheta )
       ENDDO

    ENDIF
    
  END SUBROUTINE initialize_cylindrical_boundary2D
  
  !-----------------------------------------------------------------!
  
  SUBROUTINE get_cylindrical_boundary2D(psi_cyl,psi_sph, &
       psi_sph_dx, psi_sph_dy, method)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi_sph(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi_sph_dy(:, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    
    INTEGER                           :: inum, ir, itheta
    
    
    psi_scatt = ZERO
    psi_sph = ZERO
    psi_sph_dx = ZERO
    psi_sph_dy = ZERO
    
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cyl(index_x1(inum),index_x2(inum))       
    ENDDO
    
    CALL create_interpolant(numpts,rpts_scatt,theta_scatt,psi_scatt,TRIM(method))
    
    psi_sph_dx = ZERO
    psi_sph_dy = ZERO
    
    DO ir = 1, numrpts
       DO itheta = 1, numthetapts
          CALL interpolate(numpts,rpts(ir), theta(itheta), rpts_scatt, theta_scatt, &
               psi_scatt, TRIM(method), psi_sph(ir,itheta), &
               psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta))
       ENDDO
    ENDDO
    
    
  END SUBROUTINE get_cylindrical_boundary2D
  
  !----------------------------------------------------------------!
  
  
!!$  SUBROUTINE cartesian2spherical2D( psi_cart, dims, coords, rb, &
!!$       psi_sph, d1psi_sph, d2psi_sph )
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    COMPLEX(dp), INTENT(IN)             :: psi_cart(:)
!!$    INTEGER, INTENT(IN)                 :: dims(:)
!!$    REAL(dp), INTENT(IN)                :: coords(:, :, :)
!!$    REAL(dp), INTENT(IN)                :: rb
!!$    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:)
!!$
!!$    INTEGER                             :: maxxpts, maxypts, maxzpts
!!$    
!!$
!!$
!!$    maxxpts = dims(1)
!!$    maxypts = dims(2)
!!$    maxzpts = dims(3)
!!$    
!!$    ! Check if the wavefunction and the boundary are in the same
!!$    ! processor
!!$    !IF ()
!!$    
!!$    
!!$    ! Transform the wavefunction to spherical.
!!$    
!!$    
!!$    
!!$    ! Interpolate the the coordinates-transformed wavefunction onto
!!$    ! a regular grid.
!!$
!!$    
!!$    ! Project the wavefunction onto the spherical harmonics basis.
!!$    
!!$    
!!$  END SUBROUTINE cartesian2spherical2D
  
  
  !----------------------------------------------------------------!
  
  SUBROUTINE cartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dphi(:, :, :)
    
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp), ALLOCATABLE               :: r_inp(:)
    REAL(dp), ALLOCATABLE               :: theta_inp(:)
    REAL(dp), ALLOCATABLE               :: phi_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)
 
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(r_inp(1:numpts_in))
    ALLOCATE(theta_inp(1:numpts_in))
    ALLOCATE(phi_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
    
    ii = 0
    
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             
             r_inp(ii) = rpt
             theta_inp(ii) = ACOS( z_ax(iz) / (rpt) )
             phi_inp(ii) = ATAN2( y_ax(iy) , x_ax(ix) )
             psi_inp(ii) = psi_cart(ix, iy, iz)
          ENDDO
       ENDDO
    ENDDO
    
    CALL create_interpolant(numpts_in,r_inp,theta_inp, phi_inp, &
         psi_inp,TRIM(method))
    
    
    IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
       psi_sph_dr = ZERO
       psi_sph_dth = ZERO
       psi_sph_dphi = ZERO
       
       DO iphi = 1, dims_out(3)
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                     r_inp, theta_inp, phi_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                     psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                     psi_sph_dphi(ir, itheta, iphi) )
             ENDDO
          ENDDO
       ENDDO
       
    ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
         .NOT.PRESENT(psi_sph_dth).AND. &
         .NOT. PRESENT(psi_sph_dphi) ) THEN
       
       DO iphi = 1, dims_out(3)
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                     r_inp, theta_inp, phi_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
             ENDDO
          ENDDO
       ENDDO
       
    ELSE
       WRITE(*,*) 'You must input and array for get the derivatives!'
       RETURN   
    ENDIF
    
    
  END SUBROUTINE cartesian2spherical3D
  
  !---------------------------------------------------------------!
  
  SUBROUTINE cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, method, psi_sph_dx, psi_sph_dy)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dy(:, :)
    
    INTEGER                             :: irho, iz, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    REAL(dp), ALLOCATABLE               :: r_inp(:)
    REAL(dp), ALLOCATABLE               :: theta_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)
    
    dims_in = SHAPE(psi_cyl)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    ALLOCATE(r_inp(1:numpts_in))
    ALLOCATE(theta_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
    
    ii = 0
    
    DO iz = 1, dims_in(2)
       DO irho = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(rho_ax(irho)**2 + z_ax(iz)**2 )
          
          r_inp(ii) = rpt
          theta_inp(ii) = ATAN2( z_ax(iz) , rho_ax(irho) )
          psi_inp(ii) = psi_cyl(irho, iz)
       ENDDO
    ENDDO
    
    
    CALL create_interpolant(numpts_in,r_inp,theta_inp,psi_inp,TRIM(method))
    
    IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
       psi_sph_dx = ZERO
       psi_sph_dy = ZERO
       
       DO ir = 1, numrpts
          DO itheta = 1, numthetapts
             CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), r_inp, theta_inp, &
                  psi_inp, TRIM(method), psi_sph(ir,itheta), &
                  psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta))
          ENDDO
       ENDDO
       
    ELSE
       
       DO ir = 1, numrpts
          DO itheta = 1, numthetapts
             CALL interpolate(numpts,r_ax(ir), theta_ax(itheta), r_inp, theta_inp, &
                  psi_inp, TRIM(method), psi_sph(ir,itheta))
          ENDDO
       ENDDO
       
    ENDIF
    
  END SUBROUTINE cylindrical2spherical2D
  
  
  !----------------------------------------------------------------!

END MODULE coords_transform
