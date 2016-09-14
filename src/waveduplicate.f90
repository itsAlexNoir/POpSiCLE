MODULE waveduplicate

  USE constants_pop
  USE tools
  USE cubic_interp
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           duplicate_cartesian_wave
  PUBLIC           duplicate_cylindrical_wave
  
  INTERFACE duplicate_cartesian_wave
     MODULE PROCEDURE dduplicate_cartesian_wave3D
     MODULE PROCEDURE zduplicate_cartesian_wave3D
  END INTERFACE duplicate_cartesian_wave
  
  INTERFACE duplicate_cylindrical_wave
     MODULE PROCEDURE dduplicate_cylindrical_wave2D
     MODULE PROCEDURE zduplicate_cylindrical_wave2D
  END INTERFACE duplicate_cylindrical_wave
  
CONTAINS
  
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  
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
  
  SUBROUTINE dduplicate_cartesian_wave3D(psi_original,psi_duplicated, &
       x1_ax, y1_ax, z1_ax, x2_ax, y2_ax, z2_ax, fdrule, &
       psi_dup_dx, psi_dup_dy, psi_dup_dz, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)                :: psi_original(:, :, :)
    REAL(dp), INTENT(OUT)               :: psi_duplicated(:, :, :)
    REAL(dp), INTENT(IN)                :: x1_ax(:), y1_ax(:), z1_ax(:)
    REAL(dp), INTENT(IN)                :: x2_ax(:), y2_ax(:), z2_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_dup_dx(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_dup_dy(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_dup_dz(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER, ALLOCATABLE                :: cellindex3D(: ,: ,:, :)
    REAL(dp)                            :: psi3D_dup_dx, psi3D_dup_dy
    REAL(dp)                            :: psi3D_dup_dz
    INTEGER                             :: ix, iy, iz, ipt
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: xpt, ypt, zpt
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
    
    dims_in = SHAPE(psi_original)
    dims_out = SHAPE(psi_duplicated)

    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(i_am_in3D(1:dims_out(1),1:dims_out(2),1:dims_out(3)))
    ALLOCATE(cellindex3D(1:dims_out(1),1:dims_out(2),1:dims_out(3),3))
    i_am_in3D = 0
    cellindex3D = 0
    
    psi_duplicated = 0.0_dp
    
    ! Also, initialize the bicubic matrix
    CALL initialize_tricubic_matrix( )
    
    ! Check if the new axis are within the cartesian domain.
    DO iz = 1, dims_out(3)
       DO iy = 1, dims_out(2)
          DO ix = 1, dims_out(1)
             xpt = x2_ax(ix)
             ypt = y2_ax(iy)
             zpt = z2_ax(iz)
             
             IF ((xpt.GE.MINVAL(x1_ax)) .AND. (xpt.LE.MAXVAL(x1_ax)) .AND. &
                  (ypt.GE.MINVAL(y1_ax)) .AND. (ypt.LE.MAXVAL(y1_ax)) .AND. &
                  (zpt.GE.MINVAL(z1_ax)) .AND. (zpt.LE.MAXVAL(z1_ax))) THEN
                i_am_in3D(ix,iy,iz) = 1
                
                ! Look for cell indexes
                CALL interv(x1_ax,SIZE(x1_ax),xpt,left,mflag)
                cellindex3D(ix,iy,iz,1) = left
                CALL interv(y1_ax,SIZE(y1_ax),ypt,left,mflag)
                cellindex3D(ix,iy,iz,2) = left
                CALL interv(z1_ax,SIZE(z1_ax),zpt,left,mflag)
                cellindex3D(ix,iy,iz,3) = left
                
             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Begin with the interpolation
    ! First, if we are parallel

    rulepts = 2 * fdrule
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
    
    DO iz = 1, dims_out(3)
       DO iy = 1, dims_out(2)
          DO ix = 1, dims_out(1)
             IF(i_am_in3D(ix,iy,iz).EQ.1) THEN
                psi_in       = 0.0_dp
                psi_in_dx    = 0.0_dp
                psi_in_dy    = 0.0_dp
                psi_in_dz    = 0.0_dp
                psi_in_dxdy  = 0.0_dp
                psi_in_dxdz  = 0.0_dp
                psi_in_dydz  = 0.0_dp
                psi_in_dxdydz  = 0.0_dp

                xpt = x2_ax(ix)
                ypt = y2_ax(iy)
                zpt = z2_ax(iz)

                left_x = cellindex3D(ix,iy,iz,1)
                left_y = cellindex3D(ix,iy,iz,2)
                left_z = cellindex3D(ix,iy,iz,3)

                dx = x1_ax(left_x+1) - x1_ax(left_x)
                dy = y1_ax(left_y+1) - y1_ax(left_y)
                dz = z1_ax(left_z+1) - z1_ax(left_z)

                psi_in(1) = psi_original(left_x,left_y,left_z)
                psi_in(2) = psi_original(left_x+1,left_y,left_z)
                psi_in(3) = psi_original(left_x+1,left_y+1,left_z)
                psi_in(4) = psi_original(left_x,left_y+1,left_z)
                psi_in(5) = psi_original(left_x,left_y,left_z+1)
                psi_in(6) = psi_original(left_x+1,left_y,left_z+1)
                psi_in(7) = psi_original(left_x+1,left_y+1,left_z+1)
                psi_in(8) = psi_original(left_x,left_y+1,left_z+1)

                IF(left_x-fdrule.LT.LBOUND(x1_ax,DIM=1)) THEN
                   lower_x = left_x
                   upper_x = left_x + rulepts
                ELSEIF(left_x+fdrule+1.GT.UBOUND(x1_ax,DIM=1)) THEN
                   lower_x = left_x - rulepts
                   upper_x = left_x
                ELSE
                   lower_x = left_x - fdrule
                   upper_x = left_x + fdrule
                ENDIF

                ! Calculate fd weights for differentiation in x case
                CALL fdweights(x1_ax(left_x),x1_ax(lower_x:upper_x),&
                     rulepts+1,2,xfdcoeffs)

                IF(left_y-fdrule.LT.LBOUND(y1_ax,DIM=1)) THEN
                   lower_y = left_y
                   upper_y = left_y + rulepts
                ELSEIF(left_y+fdrule+1.GT.UBOUND(y1_ax,DIM=1)) THEN
                   lower_y = left_y - rulepts
                   upper_y = left_y
                ELSE
                   lower_y = left_y - fdrule
                   upper_y = left_y + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(y1_ax(left_y),y1_ax(lower_y:upper_y),&
                     rulepts+1,2,yfdcoeffs)

                IF(left_z-fdrule.LT.LBOUND(z1_ax,DIM=1)) THEN
                   lower_z = left_z
                   upper_z = left_z + rulepts
                ELSEIF(left_z+fdrule+1.GT.UBOUND(z1_ax,DIM=1)) THEN
                   lower_z = left_z - rulepts
                   upper_z = left_z
                ELSE
                   lower_z = left_z - fdrule
                   upper_z = left_z + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in z case
                CALL fdweights(z1_ax(left_z),z1_ax(lower_z:upper_z),&
                     rulepts+1,2,zfdcoeffs)

                psi_in_vol_x(:,1) = psi_original(lower_x:upper_x,left_y,left_z)
                psi_in_vol_x(:,2) = psi_original(lower_x+1:upper_x+1,left_y,left_z)
                psi_in_vol_x(:,3) = psi_original(lower_x+1:upper_x+1,left_y+1,left_z)
                psi_in_vol_x(:,4) = psi_original(lower_x:upper_x,left_y+1,left_z)
                psi_in_vol_x(:,5) = psi_original(lower_x:upper_x,left_y,left_z+1)
                psi_in_vol_x(:,6) = psi_original(lower_x+1:upper_x+1,left_y,left_z+1)
                psi_in_vol_x(:,7) = psi_original(lower_x+1:upper_x+1,left_y+1,left_z+1)
                psi_in_vol_x(:,8) = psi_original(lower_x:upper_x,left_y+1,left_z+1)

                psi_in_vol_y(:,1) = psi_original(left_x,lower_y:upper_y,left_z)
                psi_in_vol_y(:,2) = psi_original(left_x+1,lower_y:upper_y,left_z)
                psi_in_vol_y(:,3) = psi_original(left_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,4) = psi_original(left_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,5) = psi_original(left_x,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,6) = psi_original(left_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,7) = psi_original(left_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_y(:,8) = psi_original(left_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_z(:,1) = psi_original(left_x,left_y,lower_z:upper_z)
                psi_in_vol_z(:,2) = psi_original(left_x+1,left_y,lower_z:upper_z)
                psi_in_vol_z(:,3) = psi_original(left_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,4) = psi_original(left_x,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,5) = psi_original(left_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,6) = psi_original(left_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,7) = psi_original(left_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_z(:,8) = psi_original(left_x,left_y+1,lower_z+1:upper_z+1)

                psi_in_vol_xy(:,:,1) = psi_original(lower_x:upper_x,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,2) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,3) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,4) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,5) = psi_original(lower_x:upper_x,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,6) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,7) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_xy(:,:,8) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_xz(:,:,1) = psi_original(lower_x:upper_x,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,2) = psi_original(lower_x+1:upper_x+1,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,3) = psi_original(lower_x+1:upper_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,4) = psi_original(lower_x:upper_x,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,5) = psi_original(lower_x:upper_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,6) = psi_original(lower_x+1:upper_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,7) = psi_original(lower_x+1:upper_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,8) = psi_original(lower_x:upper_x,left_y+1,lower_z+1:upper_z+1)

                psi_in_vol_yz(:,:,1) = psi_original(left_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,2) = psi_original(left_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,3) = psi_original(left_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,4) = psi_original(left_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,5) = psi_original(left_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,6) = psi_original(left_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,7) = psi_original(left_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,8) = psi_original(left_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)

                psi_in_vol_xyz(:,:,:,1) = psi_original(lower_x:upper_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,2) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,3) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,4) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,5) = psi_original(lower_x:upper_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,6) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,7) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,8) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)


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

                IF (PRESENT(psi_dup_dx).AND.PRESENT(psi_dup_dy).AND.PRESENT(psi_dup_dz) ) THEN
                   CALL tricubic_interpolation(psi_in,psi_in_dx,psi_in_dy, psi_in_dz, &
                        psi_in_dxdy,psi_in_dxdz,psi_in_dydz,psi_in_dxdydz,&
                        x1_ax(left_x),x1_ax(left_x+1),&
                        y1_ax(left_y),y1_ax(left_y+1),&
                        z1_ax(left_z),z1_ax(left_z+1),&
                        xpt,ypt,zpt,&
                        psi_duplicated(ix,iy,iz),psi_dup_dx(ix,iy,iz),&
                        psi_dup_dy(ix,iy,iz),psi_dup_dz(ix,iy,iz))

                ELSEIF (.NOT.PRESENT(psi_dup_dx).AND. &
                     .NOT.PRESENT(psi_dup_dy).AND. &
                     .NOT. PRESENT(psi_dup_dz) ) THEN
                   
                   CALL tricubic_interpolation(psi_in,psi_in_dx,psi_in_dy, psi_in_dz, &
                        psi_in_dxdy,psi_in_dxdz,psi_in_dydz,psi_in_dxdydz,&
                        x1_ax(left_x),x1_ax(left_x+1),&
                        y1_ax(left_y),y1_ax(left_y+1),&
                        z1_ax(left_z),z1_ax(left_z+1),&
                        xpt,ypt,zpt,&
                        psi_duplicated(ix,iy,iz),psi3D_dup_dx,&
                        psi3D_dup_dy,psi3D_dup_dz)
                   
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
    
    CALL delete_tricubic_matrix( )
    DEALLOCATE(i_am_in3D, cellindex3D)
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    

  END SUBROUTINE dduplicate_cartesian_wave3D

  !-------------------------------------------------------------------------!
  !-------------------------------------------------------------------------!
  
  SUBROUTINE zduplicate_cartesian_wave3D(psi_original, psi_duplicated, &
       x1_ax, y1_ax, z1_ax, x2_ax, y2_ax, z2_ax, fdrule, psi_dup_dx, psi_dup_dy, &
       psi_dup_dz, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_original(:, :, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_duplicated(:, :, :)
    REAL(dp), INTENT(IN)                :: x1_ax(:), y1_ax(:), z1_ax(:)
    REAL(dp), INTENT(IN)                :: x2_ax(:), y2_ax(:), z2_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_dup_dx(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_dup_dy(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_dup_dz(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank

    INTEGER, ALLOCATABLE                :: i_am_in3D(:, :, :)
    INTEGER, ALLOCATABLE                :: cellindex3D(: ,: ,:, :)
    COMPLEX(dp)                         :: psi3D_dup_dx, psi3D_dup_dy
    COMPLEX(dp)                         :: psi3D_dup_dz
    INTEGER                             :: ix, iy, iz, ipt
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: xpt, ypt, zpt
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
    REAL(dp)                            :: realpsi3D_dup, imagpsi3D_dup
    REAL(dp)                            :: realpsi3D_dup_dx, imagpsi3D_dup_dx
    REAL(dp)                            :: realpsi3D_dup_dy, imagpsi3D_dup_dy
    REAL(dp)                            :: realpsi3D_dup_dz, imagpsi3D_dup_dz
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

    dims_in = SHAPE(psi_original)
    dims_out = SHAPE(psi_duplicated)

    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)

    ALLOCATE(i_am_in3D(1:dims_out(1),1:dims_out(2),1:dims_out(3)))
    ALLOCATE(cellindex3D(1:dims_out(1),1:dims_out(2),1:dims_out(3),3))
    
    i_am_in3D = 0
    cellindex3D = 0

    psi_duplicated = ZERO
    
    ! Also, initialize the bicubic matrix
    CALL initialize_tricubic_matrix( )

    ! Check if the new axis are within the cartesian domain.
    DO iz = 1, dims_out(3)
       DO iy = 1, dims_out(2)
          DO ix = 1, dims_out(1)

             xpt = x2_ax(ix)
             ypt = y2_ax(iy)
             zpt = z2_ax(iz)

             IF ((xpt.GE.MINVAL(x1_ax)) .AND. (xpt.LE.MAXVAL(x1_ax)) .AND. &
                  (ypt.GE.MINVAL(y1_ax)) .AND. (ypt.LE.MAXVAL(y1_ax)) .AND. &
                  (zpt.GE.MINVAL(z1_ax)) .AND. (zpt.LE.MAXVAL(z1_ax))) THEN
                i_am_in3D(ix,iy,iz) = 1

                ! Look for cell indexes
                CALL interv(x1_ax,SIZE(x1_ax),xpt,left,mflag)
                cellindex3D(ix,iy,iz,1) = left
                CALL interv(y1_ax,SIZE(y1_ax),ypt,left,mflag)
                cellindex3D(ix,iy,iz,2) = left
                CALL interv(z1_ax,SIZE(z1_ax),zpt,left,mflag)
                cellindex3D(ix,iy,iz,3) = left

             ENDIF
          ENDDO
       ENDDO
    ENDDO

    ! Begin with the interpolation
    ! First, if we are parallel

    rulepts = 2 * fdrule
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
    
    DO iz = 1, dims_out(3)
       DO iy = 1, dims_out(2)
          DO ix = 1, dims_out(1)
             IF(i_am_in3D(ix,iy,iz).EQ.1) THEN
                psi_in       = 0.0_dp
                psi_in_dx    = 0.0_dp
                psi_in_dy    = 0.0_dp
                psi_in_dz    = 0.0_dp
                psi_in_dxdy  = 0.0_dp
                psi_in_dxdz  = 0.0_dp
                psi_in_dydz  = 0.0_dp
                psi_in_dxdydz  = 0.0_dp
                
                xpt = x2_ax(ix)
                ypt = y2_ax(iy)
                zpt = z2_ax(iz)
                
                left_x = cellindex3D(ix,iy,iz,1)
                left_y = cellindex3D(ix,iy,iz,2)
                left_z = cellindex3D(ix,iy,iz,3)
                
                dx = x1_ax(left_x+1) - x1_ax(left_x)
                dy = y1_ax(left_y+1) - y1_ax(left_y)
                dz = z1_ax(left_z+1) - z1_ax(left_z)
                
                psi_in(1) = psi_original(left_x,left_y,left_z)
                psi_in(2) = psi_original(left_x+1,left_y,left_z)
                psi_in(3) = psi_original(left_x+1,left_y+1,left_z)
                psi_in(4) = psi_original(left_x,left_y+1,left_z)
                psi_in(5) = psi_original(left_x,left_y,left_z+1)
                psi_in(6) = psi_original(left_x+1,left_y,left_z+1)
                psi_in(7) = psi_original(left_x+1,left_y+1,left_z+1)
                psi_in(8) = psi_original(left_x,left_y+1,left_z+1)

                IF(left_x-fdrule.LT.LBOUND(x1_ax,DIM=1)) THEN
                   lower_x = left_x
                   upper_x = left_x + rulepts
                ELSEIF(left_x+fdrule+1.GT.UBOUND(x1_ax,DIM=1)) THEN
                   lower_x = left_x - rulepts
                   upper_x = left_x
                ELSE
                   lower_x = left_x - fdrule
                   upper_x = left_x + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(x1_ax(left_x),x1_ax(lower_x:upper_x),&
                     rulepts+1,2,xfdcoeffs)

                IF(left_y-fdrule.LT.LBOUND(y1_ax,DIM=1)) THEN
                   lower_y = left_y
                   upper_y = left_y + rulepts
                ELSEIF(left_y+fdrule+1.GT.UBOUND(y1_ax,DIM=1)) THEN
                   lower_y = left_y - rulepts
                   upper_y = left_y
                ELSE
                   lower_y = left_y - fdrule
                   upper_y = left_y + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in x case
                CALL fdweights(y1_ax(left_y),y1_ax(lower_y:upper_y),&
                     rulepts+1,2,yfdcoeffs)

                IF(left_z-fdrule.LT.LBOUND(z1_ax,DIM=1)) THEN
                   lower_z = left_z
                   upper_z = left_z + rulepts
                ELSEIF(left_z+fdrule+1.GT.UBOUND(z1_ax,DIM=1)) THEN
                   lower_z = left_z - rulepts
                   upper_z = left_z
                ELSE
                   lower_z = left_z - fdrule
                   upper_z = left_z + fdrule
                ENDIF
                ! Calculate fd weights for differentiation in z case
                CALL fdweights(z1_ax(left_z),z1_ax(lower_z:upper_z),&
                     rulepts+1,2,zfdcoeffs)

                psi_in_vol_x(:,1) = psi_original(lower_x:upper_x,left_y,left_z)
                psi_in_vol_x(:,2) = psi_original(lower_x+1:upper_x+1,left_y,left_z)
                psi_in_vol_x(:,3) = psi_original(lower_x+1:upper_x+1,left_y+1,left_z)
                psi_in_vol_x(:,4) = psi_original(lower_x:upper_x,left_y+1,left_z)
                psi_in_vol_x(:,5) = psi_original(lower_x:upper_x,left_y,left_z+1)
                psi_in_vol_x(:,6) = psi_original(lower_x+1:upper_x+1,left_y,left_z+1)
                psi_in_vol_x(:,7) = psi_original(lower_x+1:upper_x+1,left_y+1,left_z+1)
                psi_in_vol_x(:,8) = psi_original(lower_x:upper_x,left_y+1,left_z+1)

                psi_in_vol_y(:,1) = psi_original(left_x,lower_y:upper_y,left_z)
                psi_in_vol_y(:,2) = psi_original(left_x+1,lower_y:upper_y,left_z)
                psi_in_vol_y(:,3) = psi_original(left_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,4) = psi_original(left_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_y(:,5) = psi_original(left_x,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,6) = psi_original(left_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_y(:,7) = psi_original(left_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_y(:,8) = psi_original(left_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_z(:,1) = psi_original(left_x,left_y,lower_z:upper_z)
                psi_in_vol_z(:,2) = psi_original(left_x+1,left_y,lower_z:upper_z)
                psi_in_vol_z(:,3) = psi_original(left_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,4) = psi_original(left_x,left_y+1,lower_z:upper_z)
                psi_in_vol_z(:,5) = psi_original(left_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,6) = psi_original(left_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_z(:,7) = psi_original(left_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_z(:,8) = psi_original(left_x,left_y+1,lower_z+1:upper_z+1)

                psi_in_vol_xy(:,:,1) = psi_original(lower_x:upper_x,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,2) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,left_z)
                psi_in_vol_xy(:,:,3) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,4) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,left_z)
                psi_in_vol_xy(:,:,5) = psi_original(lower_x:upper_x,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,6) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,left_z+1)
                psi_in_vol_xy(:,:,7) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,left_z+1)
                psi_in_vol_xy(:,:,8) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,left_z+1)

                psi_in_vol_xz(:,:,1) = psi_original(lower_x:upper_x,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,2) = psi_original(lower_x+1:upper_x+1,left_y,lower_z:upper_z)
                psi_in_vol_xz(:,:,3) = psi_original(lower_x+1:upper_x+1,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,4) = psi_original(lower_x:upper_x,left_y+1,lower_z:upper_z)
                psi_in_vol_xz(:,:,5) = psi_original(lower_x:upper_x,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,6) = psi_original(lower_x+1:upper_x+1,left_y,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,7) = psi_original(lower_x+1:upper_x+1,left_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xz(:,:,8) = psi_original(lower_x:upper_x,left_y+1,lower_z+1:upper_z+1)

                psi_in_vol_yz(:,:,1) = psi_original(left_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,2) = psi_original(left_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_yz(:,:,3) = psi_original(left_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,4) = psi_original(left_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_yz(:,:,5) = psi_original(left_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,6) = psi_original(left_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,7) = psi_original(left_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_yz(:,:,8) = psi_original(left_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)

                psi_in_vol_xyz(:,:,:,1) = psi_original(lower_x:upper_x,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,2) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,3) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,4) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,lower_z:upper_z)
                psi_in_vol_xyz(:,:,:,5) = psi_original(lower_x:upper_x,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,6) = psi_original(lower_x+1:upper_x+1,lower_y:upper_y,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,7) = psi_original(lower_x+1:upper_x+1,lower_y+1:upper_y+1,lower_z+1:upper_z+1)
                psi_in_vol_xyz(:,:,:,8) = psi_original(lower_x:upper_x,lower_y+1:upper_y+1,lower_z+1:upper_z+1)


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

                realpsi_in = REAL(psi_in, dp)
                imagpsi_in = AIMAG(psi_in)
                realpsi_in_dx = REAL(psi_in_dx, dp)
                imagpsi_in_dx = AIMAG(psi_in_dx)
                realpsi_in_dx = REAL(psi_in_dy, dp)
                imagpsi_in_dx = AIMAG(psi_in_dy)
                realpsi_in_dz = REAL(psi_in_dz, dp)
                imagpsi_in_dz = AIMAG(psi_in_dz)
                realpsi_in_dxdy = REAL(psi_in_dxdy, dp)
                imagpsi_in_dxdy = AIMAG(psi_in_dxdy)
                realpsi_in_dxdz = REAL(psi_in_dxdz, dp)
                imagpsi_in_dxdz = AIMAG(psi_in_dxdz)
                realpsi_in_dydz = REAL(psi_in_dydz, dp)
                imagpsi_in_dydz = AIMAG(psi_in_dydz)
                realpsi_in_dxdydz = REAL(psi_in_dxdydz, dp)
                imagpsi_in_dxdydz = AIMAG(psi_in_dxdydz)

                CALL tricubic_interpolation(realpsi_in,realpsi_in_dx,realpsi_in_dy,&
                     realpsi_in_dz,realpsi_in_dxdy,realpsi_in_dxdz,realpsi_in_dydz,&
                     realpsi_in_dxdydz,&
                     x1_ax(left_x),x1_ax(left_x+1),&
                     y1_ax(left_y),y1_ax(left_y+1),&
                     z1_ax(left_z),z1_ax(left_z+1),&
                     xpt,ypt,zpt,&
                     realpsi3D_dup,realpsi3D_dup_dx,realpsi3D_dup_dy,realpsi3D_dup_dz)

                CALL tricubic_interpolation(imagpsi_in,imagpsi_in_dx,imagpsi_in_dy,&
                     imagpsi_in_dz,imagpsi_in_dxdy,imagpsi_in_dxdz,imagpsi_in_dydz,&
                     imagpsi_in_dxdydz,&
                     x1_ax(left_x),x1_ax(left_x+1),&
                     y1_ax(left_y),y1_ax(left_y+1),&
                     z1_ax(left_z),z1_ax(left_z+1),&
                     xpt,ypt,zpt,&
                     imagpsi3D_dup,imagpsi3D_dup_dx,imagpsi3D_dup_dy,imagpsi3D_dup_dz)
                
                
                psi_duplicated(ix,iy,iz)    = CMPLX(realpsi3D_dup,imagpsi3D_dup, dp)
                
                IF (PRESENT(psi_dup_dx).AND.PRESENT(psi_dup_dy).AND.PRESENT(psi_dup_dz) ) THEN
                   psi_dup_dx(ix,iy,iz)   = CMPLX(realpsi3D_dup_dx,imagpsi3D_dup_dx, dp)
                   psi_dup_dy(ix,iy,iz)   = CMPLX(realpsi3D_dup_dy,imagpsi3D_dup_dy, dp)
                   psi_dup_dz(ix,iy,iz)   = CMPLX(realpsi3D_dup_dz,imagpsi3D_dup_dz, dp)

                ELSEIF (.NOT.PRESENT(psi_dup_dx).AND. &
                     .NOT.PRESENT(psi_dup_dy).AND. &
                     .NOT. PRESENT(psi_dup_dz) ) THEN
                   
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


    CALL delete_tricubic_matrix( )
    DEALLOCATE(i_am_in3D, cellindex3D)
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    
    
  END SUBROUTINE zduplicate_cartesian_wave3D
  
  !---------------------------------------------------------------!
  !
  !          _ _         _     _         _ ___         _            _         _
  !  __ _  _| (_)_ _  __| |_ _(_)__ __ _| |_  )____ __| |_  ___ _ _(_)__ __ _| |
  ! / _| || | | | ' \/ _` | '_| / _/ _` | |/ /(_-< '_ \ ' \/ -_) '_| / _/ _` | |
  ! \__|\_, |_|_|_||_\__,_|_| |_\__\__,_|_/___/__/ .__/_||_\___|_| |_\__\__,_|_|
  !     |__/                                     |_|
  !
  !---------------------------------------------------------------!

  !                          ___ ___
  !                         |_  )   \
  !                          / /| |) |
  !                         /___|___/
  !
  
  SUBROUTINE dduplicate_cylindrical_wave2D(psi_original,psi_duplicated, &
       rho1_ax, z1_ax, rho2_ax, z2_ax, fdrule, psi_dup_drho, psi_dup_dz, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)                :: psi_original(:, :)
    REAL(dp), INTENT(OUT)               :: psi_duplicated(:, :)
    REAL(dp), INTENT(IN)                :: rho1_ax(:), z1_ax(:)
    REAL(dp), INTENT(IN)                :: rho2_ax(:), z2_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_dup_drho(:, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_dup_dz(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER, ALLOCATABLE                :: i_am_in2D(: ,:)
    INTEGER, ALLOCATABLE                :: cellindex2D(:, :, :)
    REAL(dp)                            :: psi2D_dup_drho, psi2D_dup_dz
    INTEGER                             :: irho, iz, ipt
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
    
    dims_in = SHAPE(psi_original)
    dims_out = SHAPE(psi_duplicated)
    
    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    psi_duplicated = 0.0_dp
    IF(PRESENT(psi_dup_drho) .AND. PRESENT(psi_dup_dz)) THEN
       psi_dup_drho = 0.0_dp
       psi_dup_dz = 0.0_dp
    ENDIF
    
    ALLOCATE(i_am_in2D(1:dims_out(1),1:dims_out(2)))
    ALLOCATE(cellindex2D(1:dims_out(1),1:dims_out(2),2))
    i_am_in2D = 0
    cellindex2D = 0
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )
    
    ! Check if the new axis are within the originalesian domain.
    DO iz = 1, dims_out(2)
       DO irho = 1, dims_out(1)
          rhopt = rho2_ax(irho)
          zpt   = z2_ax(iz)
          
          IF ((rhopt.GE.MINVAL(rho1_ax)) .AND. (rhopt.LE.MAXVAL(rho1_ax)) .AND. &
               (zpt.GE.MINVAL(z1_ax)) .AND. (zpt.LE.MAXVAL(z1_ax))) THEN
             i_am_in2D(irho,iz) = 1
             
             ! Look for cell indexes
             CALL interv(rho1_ax,SIZE(rho1_ax),rhopt,left,mflag)
             cellindex2D(irho,iz,1) = left
             CALL interv(z1_ax,SIZE(z1_ax),zpt,left,mflag)
             cellindex2D(irho,iz,2) = left
             
          ENDIF
       ENDDO
    ENDDO
    
    rulepts = 2 * fdrule
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO iz = 1, dims_out(2)
       DO irho = 1, dims_out(1)
          IF(i_am_in2D(irho,iz).EQ.1) THEN
             psi_in      = 0.0_dp
             psi_in_drho = 0.0_dp
             psi_in_dz   = 0.0_dp
             psi_in_drhodz  = 0.0_dp
             
             rhopt = rho2_ax(irho)
             zpt   = z2_ax(iz)
             
             left_rho = cellindex2D(irho,iz,1)
             left_z   = cellindex2D(irho,iz,2)
             
             drho = rho1_ax(left_rho+1) - rho1_ax(left_rho)
             dz   = z1_ax(left_z+1) - z1_ax(left_z)
             
             psi_in = (/psi_original(left_rho,left_z), psi_original(left_rho+1,left_z), &
                  psi_original(left_rho+1,left_z+1), psi_original(left_rho,left_z+1)/)
             

             IF(left_rho-fdrule.LT.LBOUND(rho1_ax,DIM=1)) THEN
                lower_rho = left_rho
                upper_rho = left_rho + rulepts
             ELSEIF(left_rho+fdrule+1.GT.UBOUND(rho1_ax,DIM=1)) THEN
                lower_rho = left_rho - rulepts
                upper_rho = left_rho
             ELSE
                lower_rho = left_rho - fdrule
                upper_rho = left_rho + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in rho case
             CALL fdweights(rho1_ax(left_rho),rho1_ax(lower_rho:upper_rho),&
                  rulepts,2,rhofdcoeffs)
             
             IF(left_z-fdrule.LT.LBOUND(z1_ax,DIM=1)) THEN
                lower_z = left_z
                upper_z = left_z + rulepts
             ELSEIF(left_z+fdrule+1.GT.UBOUND(z1_ax,DIM=1)) THEN
                lower_z = left_z - rulepts
                upper_z = left_z
             ELSE
                lower_z = left_z - fdrule
                upper_z = left_z + fdrule
             ENDIF
             
             ! Calculate fd weights for differentiation in z case
             CALL fdweights(z1_ax(left_z),z1_ax(lower_z:upper_z),&
                  rulepts,2,zfdcoeffs)
             
             psi_in_vol_rho(:,1) = psi_original(lower_rho:upper_rho,left_z)
             psi_in_vol_rho(:,2) = psi_original(lower_rho+1:upper_rho+1,left_z)
             psi_in_vol_rho(:,3) = psi_original(lower_rho+1:upper_rho+1,left_z+1)
             psi_in_vol_rho(:,4) = psi_original(lower_rho:upper_rho,left_z+1)

             psi_in_vol_z(:,1) = psi_original(left_rho,lower_z:upper_z)
             psi_in_vol_z(:,2) = psi_original(left_rho+1,lower_z:upper_z)
             psi_in_vol_z(:,3) = psi_original(left_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_z(:,4) = psi_original(left_rho,lower_z+1:upper_z+1)

             psi_in_vol_rhoz(:,:,1) = psi_original(lower_rho:upper_rho,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,2) = psi_original(lower_rho+1:upper_rho+1,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,3) = psi_original(lower_rho+1:upper_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_rhoz(:,:,4) = psi_original(lower_rho:upper_rho,lower_z+1:upper_z+1)

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

             IF (PRESENT(psi_dup_drho).AND.PRESENT(psi_dup_dz)) THEN
                
                CALL bicubic_interpolation(psi_in,psi_in_drho, psi_in_dz, &
                     psi_in_drhodz, rho1_ax(left_rho),rho1_ax(left_rho+1),&
                     z1_ax(left_z),z1_ax(left_z+1),rhopt,zpt,&
                     psi_duplicated(irho,iz),psi_dup_drho(irho,iz),psi_dup_dz(irho,iz))
             ELSE
                
                CALL bicubic_interpolation(psi_in,psi_in_drho, psi_in_dz, &
                     psi_in_drhodz, rho1_ax(left_rho),rho1_ax(left_rho+1),&
                     z1_ax(left_z),z1_ax(left_z+1),rhopt,zpt,&
                     psi_duplicated(irho,iz),psi2D_dup_drho,psi2D_dup_dz)
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

    CALL delete_bicubic_matrix( )
    DEALLOCATE(i_am_in2D, cellindex2D)
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    
    
  END SUBROUTINE dduplicate_cylindrical_wave2D
  
  !---------------------------------------------------------------!
  
  SUBROUTINE zduplicate_cylindrical_wave2D(psi_original,psi_duplicated, &
       rho1_ax, z1_ax, rho2_ax, z2_ax, fdrule, psi_dup_drho, psi_dup_dz, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_original(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_duplicated(:, :)
    REAL(dp), INTENT(IN)                :: rho1_ax(:), z1_ax(:)
    REAL(dp), INTENT(IN)                :: rho2_ax(:), z2_ax(:)
    INTEGER, INTENT(IN)                 :: fdrule
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_dup_drho(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_dup_dz(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER, ALLOCATABLE                :: i_am_in2D(: ,:)
    INTEGER, ALLOCATABLE                :: cellindex2D(:, :, :)
    COMPLEX(dp)                         :: psi2D_dup_drho, psi2D_dup_dz
    INTEGER                             :: irho, iz, ipt
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rhopt, zpt
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
    REAL(dp)                            :: realpsi2D_dup, imagpsi2D_dup
    REAL(dp)                            :: realpsi2D_dup_drho, imagpsi2D_dup_drho
    REAL(dp)                            :: realpsi2D_dup_dz, imagpsi2D_dup_dz
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
    
    dims_in = SHAPE(psi_original)
    dims_out = SHAPE(psi_duplicated)
    
    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    psi_duplicated = ZERO
    IF(PRESENT(psi_dup_drho) .AND. PRESENT(psi_dup_dz)) THEN
       psi_dup_drho = ZERO
       psi_dup_dz = ZERO
    ENDIF
    
    ALLOCATE(i_am_in2D(1:dims_out(1),1:dims_out(2)))
    ALLOCATE(cellindex2D(1:dims_out(1),1:dims_out(2),2))
    i_am_in2D = 0
    cellindex2D = 0
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )
    
    ! Check if the new axis are within the cartesian domain.
    DO iz = 1, dims_out(2)
       DO irho = 1, dims_out(1)
          rhopt = rho2_ax(irho)
          zpt   = z2_ax(iz)
          
          IF ((rhopt.GE.MINVAL(rho1_ax)) .AND. (rhopt.LE.MAXVAL(rho1_ax)) .AND. &
               (zpt.GE.MINVAL(z1_ax)) .AND. (zpt.LE.MAXVAL(z1_ax))) THEN
             i_am_in2D(irho,iz) = 1
             
             ! Look for cell indexes
             CALL interv(rho1_ax,SIZE(rho1_ax),rhopt,left,mflag)
             cellindex2D(irho,iz,1) = left
             CALL interv(z1_ax,SIZE(z1_ax),zpt,left,mflag)
             cellindex2D(irho,iz,2) = left
             
          ENDIF
       ENDDO
    ENDDO
    
    rulepts = 2 * fdrule
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO iz = 1, dims_out(2)
       DO irho = 1, dims_out(1)
          IF(i_am_in2D(irho,iz).EQ.1) THEN
             psi_in      = ZERO
             psi_in_drho = ZERO
             psi_in_dz   = ZERO
             psi_in_drhodz  = ZERO
             
             rhopt = rho2_ax(irho)
             zpt   = z2_ax(iz)
             
             left_rho = cellindex2D(irho,iz,1)
             left_z   = cellindex2D(irho,iz,2)
             
             drho = rho1_ax(left_rho+1) - rho1_ax(left_rho)
             dz   = z1_ax(left_z+1) - z1_ax(left_z)
             
             psi_in = (/psi_original(left_rho,left_z), psi_original(left_rho+1,left_z), &
                  psi_original(left_rho+1,left_z+1), psi_original(left_rho,left_z+1)/)
             
             
             IF(left_rho-fdrule.LT.LBOUND(rho1_ax,DIM=1)) THEN
                lower_rho = left_rho
                upper_rho = left_rho + rulepts
             ELSEIF(left_rho+fdrule+1.GT.UBOUND(rho1_ax,DIM=1)) THEN
                lower_rho = left_rho - rulepts
                upper_rho = left_rho
             ELSE
                lower_rho = left_rho - fdrule
                upper_rho = left_rho + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in rho case
             CALL fdweights(rho1_ax(left_rho),rho1_ax(lower_rho:upper_rho),&
                  rulepts,2,rhofdcoeffs)

             IF(left_z-fdrule.LT.LBOUND(z1_ax,DIM=1)) THEN
                lower_z = left_z
                upper_z = left_z + rulepts
             ELSEIF(left_z+fdrule+1.GT.UBOUND(z1_ax,DIM=1)) THEN
                lower_z = left_z - rulepts
                upper_z = left_z
             ELSE
                lower_z = left_z - fdrule
                upper_z = left_z + fdrule
             ENDIF
             ! Calculate fd weights for differentiation in z case
             CALL fdweights(z1_ax(left_z),z1_ax(lower_z:upper_z),&
                  rulepts,2,zfdcoeffs)

             psi_in_vol_rho(:,1) = psi_original(lower_rho:upper_rho,left_z)
             psi_in_vol_rho(:,2) = psi_original(lower_rho+1:upper_rho+1,left_z)
             psi_in_vol_rho(:,3) = psi_original(lower_rho+1:upper_rho+1,left_z+1)
             psi_in_vol_rho(:,4) = psi_original(lower_rho:upper_rho,left_z+1)

             psi_in_vol_z(:,1) = psi_original(left_rho,lower_z:upper_z)
             psi_in_vol_z(:,2) = psi_original(left_rho+1,lower_z:upper_z)
             psi_in_vol_z(:,3) = psi_original(left_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_z(:,4) = psi_original(left_rho,lower_z+1:upper_z+1)

             psi_in_vol_rhoz(:,:,1) = psi_original(lower_rho:upper_rho,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,2) = psi_original(lower_rho+1:upper_rho+1,lower_z:upper_z)
             psi_in_vol_rhoz(:,:,3) = psi_original(lower_rho+1:upper_rho+1,lower_z+1:upper_z+1)
             psi_in_vol_rhoz(:,:,4) = psi_original(lower_rho:upper_rho,lower_z+1:upper_z+1)

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
             
             realpsi_in = REAL(psi_in, dp)
             imagpsi_in = AIMAG(psi_in)
             realpsi_in_drho = REAL(psi_in_drho, dp)
             imagpsi_in_drho = AIMAG(psi_in_drho)
             realpsi_in_dz = REAL(psi_in_dz, dp)
             imagpsi_in_dz = AIMAG(psi_in_dz)
             realpsi_in_drhodz = REAL(psi_in_drhodz, dp)
             imagpsi_in_drhodz = AIMAG(psi_in_drhodz)

             
             ! Real part
             CALL bicubic_interpolation(realpsi_in,realpsi_in_drho, realpsi_in_dz, &
                  realpsi_in_drhodz, rho1_ax(left_rho),rho1_ax(left_rho+1),&
                  z1_ax(left_z),z1_ax(left_z+1),rhopt,zpt,&
                  realpsi2D_dup,realpsi2D_dup_drho,realpsi2D_dup_dz)
             
             ! Imaginary part
             CALL bicubic_interpolation(imagpsi_in,imagpsi_in_drho, imagpsi_in_dz, &
                  imagpsi_in_drhodz, rho1_ax(left_rho),rho1_ax(left_rho+1),&
                  z1_ax(left_z),z1_ax(left_z+1),rhopt,zpt,&
                  imagpsi2D_dup,imagpsi2D_dup_drho,imagpsi2D_dup_dz)


             psi_duplicated(irho,iz)    = CMPLX(realpsi2D_dup,imagpsi2D_dup, dp)
             
             IF (PRESENT(psi_dup_drho).AND.PRESENT(psi_dup_dz)) THEN
                
                psi_dup_drho(irho,iz) = CMPLX(realpsi2D_dup_drho,imagpsi2D_dup_drho, dp)
                psi_dup_dz(irho,iz) = CMPLX(realpsi2D_dup_dz,imagpsi2D_dup_dz, dp)

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
    
    CALL delete_bicubic_matrix( )
    DEALLOCATE(i_am_in2D, cellindex2D)
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    
    
  END SUBROUTINE zduplicate_cylindrical_wave2D
  
  !----------------------------------------------------------------!

END MODULE waveduplicate
