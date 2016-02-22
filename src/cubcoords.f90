
!  .o88b. db    db d8888b. d888888b  .o88b.
! d8P  Y8 88    88 88  `8D   `88'   d8P  Y8
! 8P      88    88 88oooY'    88    8P
! 8b      88    88 88~~~b.    88    8b
! Y8b  d8 88b  d88 88   8D   .88.   Y8b  d8
!  `Y88P' ~Y8888P' Y8888P' Y888888P  `Y88P'
!
!
!
!  .o88b.  .d88b.   .d88b.  d8888b. d8888b. .d8888.
! d8P  Y8 .8P  Y8. .8P  Y8. 88  `8D 88  `8D 88'  YP
! 8P      88    88 88    88 88oobY' 88   88 `8bo.
! 8b      88    88 88    88 88`8b   88   88   `Y8b.
! Y8b  d8 `8b  d8' `8b  d8' 88 `88. 88  .8D db   8D
!  `Y88P'  `Y88P'   `Y88P'  88   YD Y8888D' `8888Y'


MODULE cubcoords
  
  USE constants
  USE tools
  USE cubic_interp
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           dcubic_cartesian2spherical3D
  PUBLIC           zcubic_cartesian2spherical3D
  PUBLIC           dcubic_cylindrical2spherical2D
  PUBLIC           zcubic_cylindrical2spherical2D
    
CONTAINS
  
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  !               _          _           ___         _            _         _
  !   __ __ _ _ _| |_ ___ __(_)__ _ _ _ |_  )____ __| |_  ___ _ _(_)__ __ _| |
  !  / _/ _` | '_|  _/ -_|_-< / _` | ' \ / /(_-< '_ \ ' \/ -_) '_| / _/ _` | |
  !  \__\__,_|_|  \__\___/__/_\__,_|_||_/___/__/ .__/_||_\___|_| |_\__\__,_|_|
  !                                            |_|
  
  !----------------------------------------------------------------!
  
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  !
  !                          _______
  !                         |__ /   \
  !                          |_ \ |) |
  !                         |___/___/
  !
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  SUBROUTINE dcubic_cartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, fdrule, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)                :: psi_cart(:, :, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dr(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dth(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER, ALLOCATABLE                :: cellindex3D(: ,: ,:, :)
    REAL(dp)                            :: psi3D_sph_dr, psi3D_sph_dth
    REAL(dp)                            :: psi3D_sph_dphi
    INTEGER                             :: ix, iy, iz, ii, ipt
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi
    REAL(dp)                            :: psi_in(8)
    REAL(dp)                            :: psi_in_dx(8)
    REAL(dp)                            :: psi_in_dy(8) 
    REAL(dp)                            :: psi_in_dz(8)
    REAL(dp)                            :: psi_in_dxdy(8)
    REAL(dp)                            :: psi_in_dxdz(8)
    REAL(dp)                            :: psi_in_dydz(8)
    REAL(dp)                            :: psi_in_dxdydz(8)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_x(:, :)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_y(:, :)  
    REAL(dp), ALLOCATABLE               :: psi_in_vol_z(:, :)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_xy(:, :, :) 
    REAL(dp), ALLOCATABLE               :: psi_in_vol_xz(:, :, :) 
    REAL(dp), ALLOCATABLE               :: psi_in_vol_yz(:, :, :)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_xyz(:, :, :, :) 
    REAL(dp), ALLOCATABLE               :: xfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: yfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: zfdcoeffs(:, :)  
    REAL(dp)                            :: deriv, dx, dy, dz
    INTEGER                             :: left, mflag
    INTEGER                             :: left_x, left_y, left_z
    INTEGER                             :: rulepts, centralpt
    INTEGER                             :: lower_x, lower_y, lower_z
    INTEGER                             :: upper_x, upper_y, upper_z
    
    !------------------------------------------------------------------!
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(i_am_in3D(1:dims_out(1),1:dims_out(2),1:dims_out(3)))  
    ALLOCATE(cellindex3D(1:dims_out(1),1:dims_out(2),1:dims_out(3),3))
    i_am_in3D = 0
    cellindex3D = 0
    
    psi_sph = 0.0_dp
    
    ! Also, initialize the bicubic matrix
    CALL initialize_tricubic_matrix( )
    
    
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))
             
             IF ((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in3D(ir,itheta,iphi) = 1
                
                ! Look for cell indexes
                CALL interv(x_ax,SIZE(x_ax),xpt,left,mflag)
                cellindex3D(ir,itheta,iphi,1) = left
                CALL interv(y_ax,SIZE(y_ax),ypt,left,mflag)
                cellindex3D(ir,itheta,iphi,2) = left
                CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
                cellindex3D(ir,itheta,iphi,3) = left
               
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Begin with the interpolation
    ! First, if we are parallel
    
    rulepts = 2 * fdrule + 1
    ALLOCATE(xfdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(yfdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_x(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_y(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xy(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xz(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_yz(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xyz(-fdrule:fdrule,-fdrule:fdrule,&
         -fdrule:fdrule,8))
    
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN
                psi_in       = 0.0_dp
                psi_in_dx    = 0.0_dp
                psi_in_dy    = 0.0_dp
                psi_in_dz    = 0.0_dp
                psi_in_dxdy  = 0.0_dp
                psi_in_dxdz  = 0.0_dp
                psi_in_dydz  = 0.0_dp
                psi_in_dxdydz  = 0.0_dp
                
                xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                zpt = r_ax(ir) * COS(theta_ax(itheta))
                
                left_x = cellindex3D(ir,itheta,iphi,1)
                left_y = cellindex3D(ir,itheta,iphi,2)
                left_z = cellindex3D(ir,itheta,iphi,3)
                
                dx = x_ax(left_x+1) - x_ax(left_x)
                dx = y_ax(left_y+1) - y_ax(left_y)
                dz = z_ax(left_z+1) - z_ax(left_z)
                
                psi_in(1) = psi_cart(left_x,left_y,left_z)
                psi_in(2) = psi_cart(left_x+1,left_y,left_z)
                psi_in(3) = psi_cart(left_x+1,left_y+1,left_z)
                psi_in(4) = psi_cart(left_x,left_y+1,left_z)
                psi_in(5) = psi_cart(left_x,left_y,left_z+1)
                psi_in(6) = psi_cart(left_x+1,left_y,left_z+1)
                psi_in(7) = psi_cart(left_x+1,left_y+1,left_z+1)
                psi_in(8) = psi_cart(left_x,left_y+1,left_z+1)
                
                IF(left_x-fdrule.LT.LBOUND(x_ax,DIM=1)) THEN
                   lower_x = left_x
                   upper_x = left_x + rulepts
                ELSEIF(left_x+fdrule+1.GT.UBOUND(x_ax,DIM=1)) THEN
                   lower_x = left_x - rulepts
                   upper_x = left_x
                ELSE
                   lower_x = left_x - fdrule
                   upper_x = left_x + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(x_ax(left_x),x_ax(lower_x:upper_x),&
                     rulepts,2,xfdcoeffs)
                
                IF(left_y-fdrule.LT.LBOUND(y_ax,DIM=1)) THEN
                   lower_y = left_y
                   upper_y = left_y + rulepts
                ELSEIF(left_y+fdrule+1.GT.UBOUND(y_ax,DIM=1)) THEN
                   lower_y = left_y - rulepts
                   upper_y = left_y
                ELSE
                   lower_y = left_y - fdrule
                   upper_y = left_y + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(y_ax(left_y),y_ax(lower_y:upper_y),&
                     rulepts,2,yfdcoeffs)
                
                IF(left_z-fdrule.LT.LBOUND(z_ax,DIM=1)) THEN
                   lower_z = left_z
                   upper_z = left_z + rulepts
                ELSEIF(left_z+fdrule+1.GT.UBOUND(z_ax,DIM=1)) THEN
                   lower_z = left_z - rulepts
                   upper_z = left_z
                ELSE
                   lower_z = left_z - fdrule
                   upper_z = left_z + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in z case             
                CALL fdweights(z_ax(left_z),z_ax(lower_z:upper_z),&
                     rulepts,2,zfdcoeffs)
                
                psi_in_vol_x(:,1) = psi_cart(lower_x:upper_x,left_y,left_z)
                psi_in_vol_x(:,2) = psi_cart(lower_x+1:upper_x+1,left_y,left_z)
                psi_in_vol_x(:,3) = psi_cart(lower_x+1:upper_x+1,left_y+1,left_z)
                psi_in_vol_x(:,4) = psi_cart(lower_x:upper_x,left_y+1,left_z)
                psi_in_vol_x(:,5) = psi_cart(lower_x:upper_x,left_y,left_z+1)
                psi_in_vol_x(:,6) = psi_cart(lower_x+1:upper_x+1,left_y,left_z+1)
                psi_in_vol_x(:,7) = psi_cart(lower_x+1:upper_x+1,left_y+1,left_z+1)
                psi_in_vol_x(:,8) = psi_cart(lower_x:upper_x,left_y+1,left_z+1)

                psi_in_vol_y(:,1) = psi_cart(left_x,lower_y:upper_y,left_z)
                psi_in_vol_y(:,2) = psi_cart(left_x+1,lower_y:upper_y,left_z)
                psi_in_vol_y(:,3) = psi_cart(left_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,4) = psi_cart(left_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,5) = psi_cart(left_x,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,6) = psi_cart(left_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,7) = psi_cart(left_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_y(:,8) = psi_cart(left_x,lower_y+1:upper_y+1,left_z+1)
                
                psi_in_vol_z(:,1) = psi_cart(left_x,left_y,lower_z:upper_z)
                psi_in_vol_z(:,2) = psi_cart(left_x+1,left_y,lower_z:upper_z)
                psi_in_vol_z(:,3) = psi_cart(left_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,4) = psi_cart(left_x,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,5) = psi_cart(left_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,6) = psi_cart(left_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,7) = psi_cart(left_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_z(:,8) = psi_cart(left_x,left_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_xy(:,:,1) = psi_cart(lower_x:upper_x,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,2) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,3) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,4) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,5) = psi_cart(lower_x:upper_x,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,6) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,7) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_xy(:,:,8) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_xz(:,:,1) = psi_cart(lower_x:upper_x,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,2) = psi_cart(lower_x+1:upper_x+1,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,3) = psi_cart(lower_x+1:upper_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,4) = psi_cart(lower_x:upper_x,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,5) = psi_cart(lower_x:upper_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,6) = psi_cart(lower_x+1:upper_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,7) = psi_cart(lower_x+1:upper_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,8) = psi_cart(lower_x:upper_x,left_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_yz(:,:,1) = psi_cart(left_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,2) = psi_cart(left_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,3) = psi_cart(left_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,4) = psi_cart(left_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,5) = psi_cart(left_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,6) = psi_cart(left_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,7) = psi_cart(left_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,8) = psi_cart(left_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_xyz(:,:,:,1) = psi_cart(lower_x:upper_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,2) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,3) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,4) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,5) = psi_cart(lower_x:upper_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,6) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,7) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,8) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                
                
                DO ipt = 1, 8
                   CALL make_derivative(psi_in_vol_x(:,ipt),deriv,fdrule,&
                        dx,xfdcoeffs(:,2))
                   psi_in_dx(ipt) = deriv
                   CALL make_derivative(psi_in_vol_y(:,ipt),deriv,fdrule,&
                        dy,yfdcoeffs(:,2))
                   psi_in_dy(ipt) = deriv
                   CALL make_derivative(psi_in_vol_z(:,ipt),deriv,fdrule,&
                        dz,zfdcoeffs(:,2))
                   psi_in_dz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xy(:,:,ipt),deriv,&
                        fdrule,dx,dy,xfdcoeffs(:,2),yfdcoeffs(:,2))
                   psi_in_dxdy(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xz(:,:,ipt),deriv,&
                        fdrule,dx,dz,xfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dxdz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_yz(:,:,ipt),deriv,&
                        fdrule,dy,dz,yfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dydz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xyz(:,:,:,ipt),deriv,&
                        fdrule,dx,dy,dz,xfdcoeffs(:,2),yfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dxdydz(ipt) = deriv
                ENDDO

                IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
                   CALL tricubic_interpolation(psi_in,psi_in_dx,psi_in_dy, psi_in_dz, &
                        psi_in_dxdy,psi_in_dxdz,psi_in_dydz,psi_in_dxdydz,&
                        x_ax(left_x),x_ax(left_x+1),&
                        y_ax(left_y),y_ax(left_y+1),&
                        z_ax(left_z),z_ax(left_z+1),&
                        xpt,ypt,zpt,&
                        psi_sph(ir,itheta,iphi),psi_sph_dr(ir,itheta,iphi),&
                        psi_sph_dth(ir,itheta,iphi),psi_sph_dphi(ir,itheta,iphi))
                   
                ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
                     .NOT.PRESENT(psi_sph_dth).AND. &
                     .NOT. PRESENT(psi_sph_dphi) ) THEN
                   
                   CALL tricubic_interpolation(psi_in,psi_in_dx,psi_in_dy, psi_in_dz, &
                        psi_in_dxdy,psi_in_dxdz,psi_in_dydz,psi_in_dxdydz,&
                        x_ax(left_x),x_ax(left_x+1),&
                        y_ax(left_y),y_ax(left_y+1),&
                        z_ax(left_z),z_ax(left_z+1),&
                        xpt,ypt,zpt,&
                        psi_sph(ir,itheta,iphi),psi3D_sph_dr,&
                        psi3D_sph_dth,psi3D_sph_dphi)
                   
                ELSE
                   WRITE(*,*) 'You must input an array for getting the derivatives!'
                   RETURN
                ENDIF
                
                !IF(PRESENT(rank)) THEN
                !   WRITE(*,'(A,I3)') '(Real) bicubic interpolate error &
                !        &in rank: ',rank
                !   STOP
                !ELSE
                !   WRITE(*,*) '(Real) bicubic interpolate error. '
                !ENDIF
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    

    DEALLOCATE(i_am_in3D, cellindex3D)
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    
    
  END SUBROUTINE dcubic_cartesian2spherical3D
  
  !-------------------------------------------------------------------------!
  !-------------------------------------------------------------------------!
  
  SUBROUTINE zcubic_cartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, fdrule, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER, ALLOCATABLE                :: cellindex3D(: ,: ,:, :)
    COMPLEX(dp)                         :: psi3D_sph_dr, psi3D_sph_dth
    COMPLEX(dp)                         :: psi3D_sph_dphi
    INTEGER                             :: ix, iy, iz, ii, ipt
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi
    COMPLEX(dp)                         :: psi_in(8)
    COMPLEX(dp)                         :: psi_in_dx(8)
    COMPLEX(dp)                         :: psi_in_dy(8) 
    COMPLEX(dp)                         :: psi_in_dz(8)
    COMPLEX(dp)                         :: psi_in_dxdy(8)
    COMPLEX(dp)                         :: psi_in_dxdz(8)
    COMPLEX(dp)                         :: psi_in_dydz(8)
    COMPLEX(dp)                         :: psi_in_dxdydz(8)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_x(:, :)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_y(:, :)  
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_z(:, :)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_xy(:, :, :) 
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_xz(:, :, :) 
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_yz(:, :, :) 
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_xyz(:, :, :, :)
    REAL(dp)                            :: realpsi_in(8), imagpsi_in(8)
    REAL(dp)                            :: realpsi_in_dx(8)
    REAL(dp)                            :: imagpsi_in_dx(8)
    REAL(dp)                            :: realpsi_in_dy(8)
    REAL(dp)                            :: imagpsi_in_dy(8)
    REAL(dp)                            :: realpsi_in_dz(8)
    REAL(dp)                            :: imagpsi_in_dz(8)
    REAL(dp)                            :: realpsi_in_dxdy(8)
    REAL(dp)                            :: imagpsi_in_dxdy(8)
    REAL(dp)                            :: realpsi_in_dxdz(8)
    REAL(dp)                            :: imagpsi_in_dxdz(8)
    REAL(dp)                            :: realpsi_in_dydz(8)
    REAL(dp)                            :: imagpsi_in_dydz(8)
    REAL(dp)                            :: realpsi_in_dxdydz(8)
    REAL(dp)                            :: imagpsi_in_dxdydz(8)
    REAL(dp)                            :: realpsi3D_sph, imagpsi3D_sph
    REAL(dp)                            :: realpsi3D_sph_dx, imagpsi3D_sph_dx
    REAL(dp)                            :: realpsi3D_sph_dy, imagpsi3D_sph_dy
    REAL(dp)                            :: realpsi3D_sph_dz, imagpsi3D_sph_dz
    REAL(dp), ALLOCATABLE               :: xfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: yfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: zfdcoeffs(:, :)
    COMPLEX(dp)                         :: deriv
    REAL(dp)                            :: dx, dy, dz
    INTEGER                             :: left, mflag
    INTEGER                             :: left_x, left_y, left_z
    INTEGER                             :: rulepts, centralpt
    INTEGER                             :: lower_x, lower_y, lower_z
    INTEGER                             :: upper_x, upper_y, upper_z
    
    !------------------------------------------------------------------!
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(i_am_in3D(1:dims_out(1),1:dims_out(2),1:dims_out(3)))  
    ALLOCATE(cellindex3D(1:dims_out(1),1:dims_out(2),1:dims_out(3),3))
    
    i_am_in3D = 0
    cellindex3D = 0
    
    psi_sph = ZERO
    
    ! Also, initialize the bicubic matrix
    CALL initialize_tricubic_matrix( )
    
    
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))
             
             IF ((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in3D(ir,itheta,iphi) = 1
                
                ! Look for cell indexes
                CALL interv(x_ax,SIZE(x_ax),xpt,left,mflag)
                cellindex3D(ir,itheta,iphi,1) = left
                CALL interv(y_ax,SIZE(y_ax),ypt,left,mflag)
                cellindex3D(ir,itheta,iphi,2) = left
                CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
                cellindex3D(ir,itheta,iphi,3) = left
                
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    ! Begin with the interpolation
    ! First, if we are parallel
    
    rulepts = 2 * fdrule + 1
    ALLOCATE(xfdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(yfdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_x(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_y(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xy(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xz(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_yz(-fdrule:fdrule,-fdrule:fdrule,8))
    ALLOCATE(psi_in_vol_xyz(-fdrule:fdrule,-fdrule:fdrule,&
         -fdrule:fdrule,8))
    
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             IF(i_am_in3D(ir,itheta,iphi).EQ.1) THEN
                psi_in       = 0.0_dp
                psi_in_dx    = 0.0_dp
                psi_in_dy    = 0.0_dp
                psi_in_dz    = 0.0_dp
                psi_in_dxdy  = 0.0_dp
                psi_in_dxdz  = 0.0_dp
                psi_in_dydz  = 0.0_dp
                psi_in_dxdydz  = 0.0_dp
                
                xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                zpt = r_ax(ir) * COS(theta_ax(itheta))
                
                left_x = cellindex3D(ir,itheta,iphi,1)
                left_y = cellindex3D(ir,itheta,iphi,2)
                left_z = cellindex3D(ir,itheta,iphi,3)
                
                dx = x_ax(left_x+1) - x_ax(left_x)
                dx = y_ax(left_y+1) - y_ax(left_y)
                dz = z_ax(left_z+1) - z_ax(left_z)
                
                psi_in(1) = psi_cart(left_x,left_y,left_z)
                psi_in(2) = psi_cart(left_x+1,left_y,left_z)
                psi_in(3) = psi_cart(left_x+1,left_y+1,left_z)
                psi_in(4) = psi_cart(left_x,left_y+1,left_z)
                psi_in(5) = psi_cart(left_x,left_y,left_z+1)
                psi_in(6) = psi_cart(left_x+1,left_y,left_z+1)
                psi_in(7) = psi_cart(left_x+1,left_y+1,left_z+1)
                psi_in(8) = psi_cart(left_x,left_y+1,left_z+1)
                
                IF(left_x-fdrule.LT.LBOUND(x_ax,DIM=1)) THEN
                   lower_x = left_x
                   upper_x = left_x + rulepts
                ELSEIF(left_x+fdrule+1.GT.UBOUND(x_ax,DIM=1)) THEN
                   lower_x = left_x - rulepts
                   upper_x = left_x
                ELSE
                   lower_x = left_x - fdrule
                   upper_x = left_x + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(x_ax(left_x),x_ax(lower_x:upper_x),&
                     rulepts,2,xfdcoeffs)
                
                IF(left_y-fdrule.LT.LBOUND(y_ax,DIM=1)) THEN
                   lower_y = left_y
                   upper_y = left_y + rulepts
                ELSEIF(left_y+fdrule+1.GT.UBOUND(y_ax,DIM=1)) THEN
                   lower_y = left_y - rulepts
                   upper_y = left_y
                ELSE
                   lower_y = left_y - fdrule
                   upper_y = left_y + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(y_ax(left_y),y_ax(lower_y:upper_y),&
                     rulepts,2,yfdcoeffs)
                
                IF(left_z-fdrule.LT.LBOUND(z_ax,DIM=1)) THEN
                   lower_z = left_z
                   upper_z = left_z + rulepts
                ELSEIF(left_z+fdrule+1.GT.UBOUND(z_ax,DIM=1)) THEN
                   lower_z = left_z - rulepts
                   upper_z = left_z
                ELSE
                   lower_z = left_z - fdrule
                   upper_z = left_z + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in z case             
                CALL fdweights(z_ax(left_z),z_ax(lower_z:upper_z),&
                     rulepts,2,zfdcoeffs)
                
                psi_in_vol_x(:,1) = psi_cart(lower_x:upper_x,left_y,left_z)
                psi_in_vol_x(:,2) = psi_cart(lower_x+1:upper_x+1,left_y,left_z)
                psi_in_vol_x(:,3) = psi_cart(lower_x+1:upper_x+1,left_y+1,left_z)
                psi_in_vol_x(:,4) = psi_cart(lower_x:upper_x,left_y+1,left_z)
                psi_in_vol_x(:,5) = psi_cart(lower_x:upper_x,left_y,left_z+1)
                psi_in_vol_x(:,6) = psi_cart(lower_x+1:upper_x+1,left_y,left_z+1)
                psi_in_vol_x(:,7) = psi_cart(lower_x+1:upper_x+1,left_y+1,left_z+1)
                psi_in_vol_x(:,8) = psi_cart(lower_x:upper_x,left_y+1,left_z+1)

                psi_in_vol_y(:,1) = psi_cart(left_x,lower_y:upper_y,left_z)
                psi_in_vol_y(:,2) = psi_cart(left_x+1,lower_y:upper_y,left_z)
                psi_in_vol_y(:,3) = psi_cart(left_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,4) = psi_cart(left_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,5) = psi_cart(left_x,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,6) = psi_cart(left_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,7) = psi_cart(left_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_y(:,8) = psi_cart(left_x,lower_y+1:upper_y+1,left_z+1)
                
                psi_in_vol_z(:,1) = psi_cart(left_x,left_y,lower_z:upper_z)
                psi_in_vol_z(:,2) = psi_cart(left_x+1,left_y,lower_z:upper_z)
                psi_in_vol_z(:,3) = psi_cart(left_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,4) = psi_cart(left_x,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,5) = psi_cart(left_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,6) = psi_cart(left_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,7) = psi_cart(left_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_z(:,8) = psi_cart(left_x,left_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_xy(:,:,1) = psi_cart(lower_x:upper_x,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,2) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,3) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,4) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,5) = psi_cart(lower_x:upper_x,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,6) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,7) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_xy(:,:,8) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_xz(:,:,1) = psi_cart(lower_x:upper_x,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,2) = psi_cart(lower_x+1:upper_x+1,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,3) = psi_cart(lower_x+1:upper_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,4) = psi_cart(lower_x:upper_x,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,5) = psi_cart(lower_x:upper_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,6) = psi_cart(lower_x+1:upper_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,7) = psi_cart(lower_x+1:upper_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,8) = psi_cart(lower_x:upper_x,left_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_yz(:,:,1) = psi_cart(left_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,2) = psi_cart(left_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,3) = psi_cart(left_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,4) = psi_cart(left_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,5) = psi_cart(left_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,6) = psi_cart(left_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,7) = psi_cart(left_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,8) = psi_cart(left_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                
                psi_in_vol_xyz(:,:,:,1) = psi_cart(lower_x:upper_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,2) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,3) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,4) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,5) = psi_cart(lower_x:upper_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,6) = psi_cart(lower_x+1:upper_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,7) = psi_cart(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,8) = psi_cart(lower_x:upper_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                
                
                DO ipt = 1, 8
                   CALL make_derivative(psi_in_vol_x(:,ipt),deriv,fdrule,&
                        dx,xfdcoeffs(:,2))
                   psi_in_dx(ipt) = deriv
                   CALL make_derivative(psi_in_vol_y(:,ipt),deriv,fdrule,&
                        dy,yfdcoeffs(:,2))
                   psi_in_dy(ipt) = deriv
                   CALL make_derivative(psi_in_vol_z(:,ipt),deriv,fdrule,&
                        dz,zfdcoeffs(:,2))
                   psi_in_dz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xy(:,:,ipt),deriv,&
                        fdrule,dx,dy,xfdcoeffs(:,2),yfdcoeffs(:,2))
                   psi_in_dxdy(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xz(:,:,ipt),deriv,&
                        fdrule,dx,dz,xfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dxdz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_yz(:,:,ipt),deriv,&
                        fdrule,dy,dz,yfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dydz(ipt) = deriv
                   CALL make_cross_derivative(psi_in_vol_xyz(:,:,:,ipt),deriv,&
                        fdrule,dx,dy,dz,xfdcoeffs(:,2),yfdcoeffs(:,2),zfdcoeffs(:,2))
                   psi_in_dxdydz(ipt) = deriv
                ENDDO
                
                realpsi_in = REAL(psi_in)
                imagpsi_in = AIMAG(psi_in)
                realpsi_in_dx = REAL(psi_in_dx)
                imagpsi_in_dx = AIMAG(psi_in_dx)
                realpsi_in_dx = REAL(psi_in_dy)
                imagpsi_in_dx = AIMAG(psi_in_dy)
                realpsi_in_dz = REAL(psi_in_dz)
                imagpsi_in_dz = AIMAG(psi_in_dz)
                realpsi_in_dxdy = REAL(psi_in_dxdy)
                imagpsi_in_dxdy = AIMAG(psi_in_dxdy)
                realpsi_in_dxdz = REAL(psi_in_dxdz)
                imagpsi_in_dxdz = AIMAG(psi_in_dxdz)
                realpsi_in_dydz = REAL(psi_in_dydz)
                imagpsi_in_dydz = AIMAG(psi_in_dydz)
                realpsi_in_dxdydz = REAL(psi_in_dxdydz)
                imagpsi_in_dxdydz = AIMAG(psi_in_dxdydz)
                
                CALL tricubic_interpolation(realpsi_in,realpsi_in_dx,realpsi_in_dy,&
                     realpsi_in_dz,realpsi_in_dxdy,realpsi_in_dxdz,realpsi_in_dydz,&
                     realpsi_in_dxdydz,&
                     x_ax(left_x),x_ax(left_x+1),&
                     y_ax(left_y),y_ax(left_y+1),&
                     z_ax(left_z),z_ax(left_z+1),&
                     xpt,ypt,zpt,&
                     realpsi3D_sph,realpsi3D_sph_dx,realpsi3D_sph_dy,realpsi3D_sph_dz)
                
                CALL tricubic_interpolation(imagpsi_in,imagpsi_in_dx,imagpsi_in_dy,&
                     imagpsi_in_dz,imagpsi_in_dxdy,imagpsi_in_dxdz,imagpsi_in_dydz,&
                     imagpsi_in_dxdydz,&
                     x_ax(left_x),x_ax(left_x+1),&
                     y_ax(left_y),y_ax(left_y+1),&
                     z_ax(left_z),z_ax(left_z+1),&
                     xpt,ypt,zpt,&
                     imagpsi3D_sph,imagpsi3D_sph_dx,imagpsi3D_sph_dy,imagpsi3D_sph_dz)
                
                
                psi_sph(ir,itheta,iphi)    = CMPLX(realpsi3D_sph,imagpsi3D_sph)
                
                IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
                   psi_sph_dr(ir,itheta,iphi)   = CMPLX(realpsi3D_sph_dx,imagpsi3D_sph_dx)
                   psi_sph_dth(ir,itheta,iphi)  = CMPLX(realpsi3D_sph_dy,imagpsi3D_sph_dy)
                   psi_sph_dphi(ir,itheta,iphi) = CMPLX(realpsi3D_sph_dz,imagpsi3D_sph_dz)
                   
                ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
                     .NOT.PRESENT(psi_sph_dth).AND. &
                     .NOT. PRESENT(psi_sph_dphi) ) THEN
                   
                ELSE
                   WRITE(*,*) 'You must input an array for getting the derivatives!'
                   RETURN
                ENDIF
                
                !IF(PRESENT(rank)) THEN
                !   WRITE(*,'(A,I3)') '(Real) bicubic interpolate error &
                !        &in rank: ',rank
                !   STOP
                !ELSE
                !   WRITE(*,*) '(Real) bicubic interpolate error. '
                !ENDIF
                
             ENDIF
             
          ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(i_am_in3D, cellindex3D)
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    
    
  END SUBROUTINE zcubic_cartesian2spherical3D
  
  !---------------------------------------------------------------!
  
  !          _ _         _     _         _ ___         _            _         _
  !  __ _  _| (_)_ _  __| |_ _(_)__ __ _| |_  )____ __| |_  ___ _ _(_)__ __ _| |
  ! / _| || | | | ' \/ _` | '_| / _/ _` | |/ /(_-< '_ \ ' \/ -_) '_| / _/ _` | |
  ! \__|\_, |_|_|_||_\__,_|_| |_\__\__,_|_/___/__/ .__/_||_\___|_| |_\__\__,_|_|
  !     |__/                                     |_|

  !---------------------------------------------------------------!
  
  !                          ___ ___
  !                         |_  )   \
  !                          / /| |) |
  !                         /___|___/
  !
  
  SUBROUTINE dcubic_cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, fdrule, psi_sph_dx, psi_sph_dy, rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)                :: psi_cyl(:, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dx(:, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(: ,:)
    INTEGER, ALLOCATABLE                :: cellindex2D(:, :, :)
    REAL(dp)                            :: psi2D_sph_dx, psi2D_sph_dy
    INTEGER                             :: irho, iz, ipt
    INTEGER                             :: ir, itheta
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rhopt, zpt, rpt
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    REAL(dp)                            :: psi_in(4), psi_in_drhodz(4)
    REAL(dp)                            :: psi_in_drho(4), psi_in_dz(4)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_rho(:, :)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_z(:, :)
    REAL(dp), ALLOCATABLE               :: psi_in_vol_rhoz(:, :, :)
    REAL(dp), ALLOCATABLE               :: rhofdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: zfdcoeffs(:, :)  
    REAL(dp)                            :: deriv
    REAL(dp)                            :: drho, dz
    INTEGER                             :: left_rho, left_z
    INTEGER                             :: rulepts, centralpt
    INTEGER                             :: lower_rho, lower_z
    INTEGER                             :: upper_rho, upper_z
    INTEGER                             :: left, mflag
       
    !-----------------------------------------------------------!
    
    dims_in = SHAPE(psi_cyl)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    psi_sph = 0.0_dp
    IF(PRESENT(psi_sph_dx) .AND. PRESENT(psi_sph_dy)) THEN
       psi_sph_dx = 0.0_dp
       psi_sph_dy = 0.0_dp
    ENDIF
    
    ALLOCATE(i_am_in2D(1:dims_out(1),1:dims_out(2)))  
    ALLOCATE(cellindex2D(1:dims_out(1),1:dims_out(2),2))
    i_am_in2D = 0
    cellindex2D = 0
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )

    
    ! Check if the new axis are within the cartesian domain.
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          rhopt = r_ax(ir) * SIN(theta_ax(itheta))
          zpt   = r_ax(ir) * COS(theta_ax(itheta))
          
          IF ((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             i_am_in2D(ir,itheta) = 1
             
             ! Look for cell indexes
             CALL interv(rho_ax,SIZE(rho_ax),rhopt,left,mflag)
             cellindex2D(ir,itheta,1) = left
             CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
             cellindex2D(ir,itheta,2) = left
             
          ENDIF
       ENDDO
    ENDDO
    
    rulepts = 2 * fdrule + 1
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          IF(i_am_in2D(ir,itheta).EQ.1) THEN
             psi_in      = 0.0_dp
             psi_in_drho = 0.0_dp
             psi_in_dz   = 0.0_dp
             psi_in_drhodz  = 0.0_dp
             
             rhopt = r_ax(ir) * SIN(theta_ax(itheta))
             zpt   = r_ax(ir) * COS(theta_ax(itheta))
             
             left_rho = cellindex2D(ir,itheta,1)
             left_z   = cellindex2D(ir,itheta,2)
             
             drho = rho_ax(left_rho+1) - rho_ax(left_rho)
             dz   = z_ax(left_z+1) - z_ax(left_z)
             
             psi_in = (/psi_cyl(left_rho,left_z), psi_cyl(left_rho+1,left_z), &
                  psi_cyl(left_rho+1,left_z+1), psi_cyl(left_rho,left_z+1)/)
             
             
             IF(left_rho-fdrule.LT.LBOUND(rho_ax,DIM=1)) THEN
                lower_rho = left_rho
                upper_rho = left_rho + rulepts
             ELSEIF(left_rho+fdrule+1.GT.UBOUND(rho_ax,DIM=1)) THEN
                lower_rho = left_rho - rulepts
                upper_rho = left_rho
             ELSE
                lower_rho = left_rho - fdrule
                upper_rho = left_rho + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in rho case
             CALL fdweights(rho_ax(left_rho),rho_ax(lower_rho:upper_rho),&
                  rulepts,2,rhofdcoeffs)
             
             IF(left_z-fdrule.LT.LBOUND(z_ax,DIM=1)) THEN
                lower_z = left_z
                upper_z = left_z + rulepts
             ELSEIF(left_z+fdrule+1.GT.UBOUND(z_ax,DIM=1)) THEN
                lower_z = left_z - rulepts
                upper_z = left_z
             ELSE
                lower_z = left_z - fdrule
                upper_z = left_z + fdrule
             ENDIF
             
             ! Calculate fd weights for differentiation in z case             
             CALL fdweights(z_ax(left_z),z_ax(lower_z:upper_z),&
                  rulepts,2,zfdcoeffs)
             
             psi_in_vol_rho(:,1) = psi_cyl(lower_rho:upper_rho,left_z)
             psi_in_vol_rho(:,2) = psi_cyl(lower_rho+1:upper_rho+1,left_z)
             psi_in_vol_rho(:,3) = psi_cyl(lower_rho+1:upper_rho+1,left_z+1)
             psi_in_vol_rho(:,4) = psi_cyl(lower_rho:upper_rho,left_z+1)
             
             psi_in_vol_z(:,1) = psi_cyl(left_rho,lower_z:upper_z)
             psi_in_vol_z(:,2) = psi_cyl(left_rho+1,lower_z:upper_z)
             psi_in_vol_z(:,3) = psi_cyl(left_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_z(:,4) = psi_cyl(left_rho,lower_z+1:upper_z+1)

             psi_in_vol_rhoz(:,:,1) = psi_cyl(lower_rho:upper_rho,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,2) = psi_cyl(lower_rho+1:upper_rho+1,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,3) = psi_cyl(lower_rho+1:upper_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_rhoz(:,:,4) = psi_cyl(lower_rho:upper_rho,lower_z+1:upper_z+1)
             
             DO ipt = 1, 4
                CALL make_derivative(psi_in_vol_rho(:,ipt),deriv,fdrule,&
                     drho,rhofdcoeffs(:,2))
                psi_in_drho(ipt) = deriv
                CALL make_derivative(psi_in_vol_z(:,ipt),deriv,fdrule,&
                     dz,zfdcoeffs(:,2))
                psi_in_dz(ipt) = deriv
                CALL make_cross_derivative(psi_in_vol_rhoz(:,:,ipt),deriv,&
                     fdrule,drho, dz, rhofdcoeffs(:,2),zfdcoeffs(:,2))
                psi_in_drhodz(ipt) = deriv
             ENDDO
             
             IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
                
                CALL bicubic_interpolation(psi_in,psi_in_drho, psi_in_dz, &
                     psi_in_drhodz, rho_ax(left_rho),rho_ax(left_rho+1),&
                     z_ax(left_z),z_ax(left_z+1),rhopt,zpt,&
                     psi_sph(ir,itheta),psi_sph_dx(ir,itheta),psi_sph_dy(ir,itheta))
             ELSE
                
                CALL bicubic_interpolation(psi_in,psi_in_drho, psi_in_dz, &
                     psi_in_drhodz, rho_ax(left_rho),rho_ax(left_rho+1),&
                     z_ax(left_z),z_ax(left_z+1),rhopt,zpt,&
                     psi_sph(ir,itheta),psi2D_sph_dx,psi2D_sph_dy)
             ENDIF
             
             
             !IF(PRESENT(rank)) THEN
             !   WRITE(*,'(A,I3)') '(Complex) bicubic interpolate error &
             !        &in rank: ',rank
             !   STOP
             !ELSE
             !   WRITE(*,*) '(Complex) bicubic interpolate error. '
             !ENDIF
             
          ENDIF
          
       ENDDO
    ENDDO
    
    DEALLOCATE(i_am_in2D, cellindex2D)
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    
    
  END SUBROUTINE dcubic_cylindrical2spherical2D
  
  !---------------------------------------------------------------!
  
  SUBROUTINE zcubic_cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, fdrule, psi_sph_dx, psi_sph_dy, rank)

    IMPLICIT NONE

    COMPLEX(dp), INTENT(IN)             :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in2D(: ,:)
    INTEGER, ALLOCATABLE                :: cellindex2D(:, :, :)
    COMPLEX(dp)                         :: psi2D_sph_dx, psi2D_sph_dy
    INTEGER                             :: irho, iz, ipt
    INTEGER                             :: ir, itheta
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rhopt, zpt, rpt
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    COMPLEX(dp)                         :: psi_in(4), psi_in_drhodz(4)
    COMPLEX(dp)                         :: psi_in_drho(4), psi_in_dz(4)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_rho(:, :)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_z(:, :)
    COMPLEX(dp), ALLOCATABLE            :: psi_in_vol_rhoz(:, :, :)
    REAL(dp)                            :: realpsi_in(4), imagpsi_in(4)
    REAL(dp)                            :: realpsi_in_drho(4)
    REAL(dp)                            :: imagpsi_in_drho(4)
    REAL(dp)                            :: realpsi_in_dz(4)
    REAL(dp)                            :: imagpsi_in_dz(4)
    REAL(dp)                            :: realpsi_in_drhodz(4)
    REAL(dp)                            :: imagpsi_in_drhodz(4)
    REAL(dp)                            :: realpsi2D_sph, imagpsi2D_sph
    REAL(dp)                            :: realpsi2D_sph_dx, imagpsi2D_sph_dx
    REAL(dp)                            :: realpsi2D_sph_dy, imagpsi2D_sph_dy
    REAL(dp), ALLOCATABLE               :: rhofdcoeffs(:, :)
    REAL(dp), ALLOCATABLE               :: zfdcoeffs(:, :)  
    COMPLEX(dp)                         :: deriv
    REAL(dp)                            :: drho, dz
    INTEGER                             :: left_rho, left_z
    INTEGER                             :: rulepts, centralpt
    INTEGER                             :: lower_rho, lower_z
    INTEGER                             :: upper_rho, upper_z
    INTEGER                             :: left, mflag
    
    !-----------------------------------------------------------!
    
    dims_in = SHAPE(psi_cyl)
    dims_out = SHAPE(psi_sph)

    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    psi_sph = ZERO
    IF(PRESENT(psi_sph_dx) .AND. PRESENT(psi_sph_dy)) THEN
       psi_sph_dx = ZERO
       psi_sph_dy = ZERO
    ENDIF
    
    ALLOCATE(i_am_in2D(1:dims_out(1),1:dims_out(2)))  
    ALLOCATE(cellindex2D(1:dims_out(1),1:dims_out(2),2))
    i_am_in2D = 0
    cellindex2D = 0
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )

    
    ! Check if the new axis are within the cartesian domain.
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          rhopt = r_ax(ir) * SIN(theta_ax(itheta))
          zpt   = r_ax(ir) * COS(theta_ax(itheta))
          
          IF ((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             i_am_in2D(ir,itheta) = 1
             
             ! Look for cell indexes
             CALL interv(rho_ax,SIZE(rho_ax),rhopt,left,mflag)
             cellindex2D(ir,itheta,1) = left
             CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
             cellindex2D(ir,itheta,2) = left
             
          ENDIF
       ENDDO
    ENDDO
    
    rulepts = 2 * fdrule + 1
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO itheta = 1, dims_out(2)
       DO ir = 1, dims_out(1)
          IF(i_am_in2D(ir,itheta).EQ.1) THEN
             psi_in      = ZERO
             psi_in_drho = ZERO
             psi_in_dz   = ZERO
             psi_in_drhodz  = ZERO
             
             rhopt = r_ax(ir) * SIN(theta_ax(itheta))
             zpt   = r_ax(ir) * COS(theta_ax(itheta))
             
             left_rho = cellindex2D(ir,itheta,1)
             left_z   = cellindex2D(ir,itheta,2)
             
             drho = rho_ax(left_rho+1) - rho_ax(left_rho)
             dz   = z_ax(left_z+1) - z_ax(left_z)
             
             psi_in = (/psi_cyl(left_rho,left_z), psi_cyl(left_rho+1,left_z), &
                  psi_cyl(left_rho+1,left_z+1), psi_cyl(left_rho,left_z+1)/)
             
             
             IF(left_rho-fdrule.LT.LBOUND(rho_ax,DIM=1)) THEN
                lower_rho = left_rho
                upper_rho = left_rho + rulepts
             ELSEIF(left_rho+fdrule+1.GT.UBOUND(rho_ax,DIM=1)) THEN
                lower_rho = left_rho - rulepts
                upper_rho = left_rho
             ELSE
                lower_rho = left_rho - fdrule
                upper_rho = left_rho + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in rho case
             CALL fdweights(rho_ax(left_rho),rho_ax(lower_rho:upper_rho),&
                  rulepts,2,rhofdcoeffs)
             
             IF(left_z-fdrule.LT.LBOUND(z_ax,DIM=1)) THEN
                lower_z = left_z
                upper_z = left_z + rulepts
             ELSEIF(left_z+fdrule+1.GT.UBOUND(z_ax,DIM=1)) THEN
                lower_z = left_z - rulepts
                upper_z = left_z
             ELSE
                lower_z = left_z - fdrule
                upper_z = left_z + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in z case             
             CALL fdweights(z_ax(left_z),z_ax(lower_z:upper_z),&
                  rulepts,2,zfdcoeffs)

             psi_in_vol_rho(:,1) = psi_cyl(lower_rho:upper_rho,left_z)
             psi_in_vol_rho(:,2) = psi_cyl(lower_rho+1:upper_rho+1,left_z)
             psi_in_vol_rho(:,3) = psi_cyl(lower_rho+1:upper_rho+1,left_z+1)
             psi_in_vol_rho(:,4) = psi_cyl(lower_rho:upper_rho,left_z+1)
             
             psi_in_vol_z(:,1) = psi_cyl(left_rho,lower_z:upper_z)
             psi_in_vol_z(:,2) = psi_cyl(left_rho+1,lower_z:upper_z)
             psi_in_vol_z(:,3) = psi_cyl(left_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_z(:,4) = psi_cyl(left_rho,lower_z+1:upper_z+1)

             psi_in_vol_rhoz(:,:,1) = psi_cyl(lower_rho:upper_rho,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,2) = psi_cyl(lower_rho+1:upper_rho+1,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,3) = psi_cyl(lower_rho+1:upper_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_rhoz(:,:,4) = psi_cyl(lower_rho:upper_rho,lower_z+1:upper_z+1)
             
             DO ipt = 1, 4
                CALL make_derivative(psi_in_vol_rho(:,ipt),deriv,fdrule,&
                     drho,rhofdcoeffs(:,2))
                psi_in_drho(ipt) = deriv
                CALL make_derivative(psi_in_vol_z(:,ipt),deriv,fdrule,&
                     dz,zfdcoeffs(:,2))
                psi_in_dz(ipt) = deriv
                CALL make_cross_derivative(psi_in_vol_rhoz(:,:,ipt),deriv,&
                     fdrule,drho, dz, rhofdcoeffs(:,2),zfdcoeffs(:,2))
                psi_in_drhodz(ipt) = deriv
             ENDDO

             realpsi_in = REAL(psi_in)
             imagpsi_in = AIMAG(psi_in)
             realpsi_in_drho = REAL(psi_in_drho)
             imagpsi_in_drho = AIMAG(psi_in_drho)
             realpsi_in_dz = REAL(psi_in_dz)
             imagpsi_in_dz = AIMAG(psi_in_dz)
             realpsi_in_drhodz = REAL(psi_in_drhodz)
             imagpsi_in_drhodz = AIMAG(psi_in_drhodz)
             
             
             ! Real part
             CALL bicubic_interpolation(realpsi_in,realpsi_in_drho, realpsi_in_dz, &
                  realpsi_in_drhodz, rho_ax(left_rho),rho_ax(left_rho+1),&
                  z_ax(left_z),z_ax(left_z+1),rhopt,zpt,&
                  realpsi2D_sph,realpsi2D_sph_dx,realpsi2D_sph_dy)
             
             ! Imaginary part
             CALL bicubic_interpolation(imagpsi_in,imagpsi_in_drho, imagpsi_in_dz, &
                  imagpsi_in_drhodz, rho_ax(left_rho),rho_ax(left_rho+1),&
                  z_ax(left_z),z_ax(left_z+1),rhopt,zpt,&
                  imagpsi2D_sph,imagpsi2D_sph_dx,imagpsi2D_sph_dy)
             

             psi_sph(ir,itheta)    = CMPLX(realpsi2D_sph,imagpsi2D_sph)
             
             
             IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
                
                psi_sph_dx(ir,itheta) = CMPLX(realpsi2D_sph_dx,imagpsi2D_sph_dx)
                psi_sph_dy(ir,itheta) = CMPLX(realpsi2D_sph_dy,imagpsi2D_sph_dy)
                
             ENDIF
             
             !IF(PRESENT(rank)) THEN
             !   WRITE(*,'(A,I3)') '(Complex) bicubic interpolate error &
             !        &in rank: ',rank
             !   STOP
             !ELSE
             !   WRITE(*,*) '(Complex) bicubic interpolate error. '
             !ENDIF
             
          ENDIF
          
       ENDDO
    ENDDO
    
    DEALLOCATE(i_am_in2D, cellindex2D)
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    

  END SUBROUTINE zcubic_cylindrical2spherical2D
  
  !----------------------------------------------------------------!
  
END MODULE cubcoords
