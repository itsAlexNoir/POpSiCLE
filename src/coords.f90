MODULE coords
  
  USE constants
  USE interp
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           cartesian2spherical
  PUBLIC           cylindrical2spherical
  

  INTERFACE cartesian2spherical
     MODULE PROCEDURE cartesian2spherical2D
     MODULE PROCEDURE dcartesian2spherical3D
     MODULE PROCEDURE zcartesian2spherical3D
  END INTERFACE cartesian2spherical
  
  INTERFACE cylindrical2spherical
     MODULE PROCEDURE cylindrical2spherical2D
  END INTERFACE cylindrical2spherical
  
  INTEGER                          :: numpts, numrpts
  INTEGER                          :: numthetapts, numphipts
  
CONTAINS
 
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  SUBROUTINE cartesian2spherical2D(psi_cart,psi_sph, x_ax, y_ax, &
       r_ax, theta_ax, method, psi_sph_dr, psi_sph_dth,rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER                             :: ix, iy, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt, xpt, ypt
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)
    
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
    
    ii = 0
    
    DO iy = 1, dims_in(2)
       DO ix = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 )
          
          x_inp(ii) = x_ax(ix)
          y_inp(ii) = y_ax(iy)
          psi_inp(ii) = psi_cart(ix, iy)
       ENDDO
    ENDDO
    
    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method),rank)
       
       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                
                xpt = r_ax(ir) * COS(theta_ax(itheta))
                ypt = r_ax(ir) * SIN(theta_ax(itheta))
                
                CALL interpolate(numpts_in,xpt, ypt, &
                     x_inp, y_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta), &
                     psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta),rank )
             ENDDO
          ENDDO
          
          
       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN
          
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)

                xpt = r_ax(ir) * COS(theta_ax(itheta))
                ypt = r_ax(ir) * SIN(theta_ax(itheta))
                
                CALL interpolate(numpts_in,xpt, ypt, &
                     x_inp, y_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank )
             ENDDO
          ENDDO
          
       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN   
       ENDIF
       
    ELSE
       
       CALL create_interpolant(numpts_in,x_inp,y_inp,psi_inp,TRIM(method))
       
       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)

                xpt = r_ax(ir) * COS(theta_ax(itheta))
                ypt = r_ax(ir) * SIN(theta_ax(itheta))
                
                CALL interpolate(numpts_in,xpt, ypt, &
                     x_inp, y_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta), &
                     psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta) )
             ENDDO
          ENDDO
          
          
       ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
            .NOT.PRESENT(psi_sph_dth) ) THEN
          
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)

                xpt = r_ax(ir) * COS(theta_ax(itheta))
                ypt = r_ax(ir) * SIN(theta_ax(itheta))
                    
                CALL interpolate(numpts_in,xpt, ypt, &
                     x_inp, y_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta) )
             ENDDO
          ENDDO
          
       ELSE
          WRITE(*,*) 'You must input and array for get the derivatives!'
          RETURN   
       ENDIF
       
    ENDIF
    
  END SUBROUTINE cartesian2spherical2D
  
  !----------------------------------------------------------------!
  
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  SUBROUTINE dcartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, i_am_in, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)                :: psi_cart(:, :, :)
    REAL(dp), INTENT(OUT)               :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    INTEGER, INTENT(OUT),OPTIONAL       :: i_am_in(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dr(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dth(:, :, :)
    REAL(dp), INTENT(OUT), OPTIONAL     :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi 
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    REAL(dp), ALLOCATABLE               :: psi_inp(:)
    
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
        
    ii = 0
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             
             x_inp(ii) = x_ax(ix)
             y_inp(ii) = y_ax(iy)
             z_inp(ii) = z_ax(iz)
             psi_inp(ii) = psi_cart(ix,iy,iz)
             
          ENDDO
       ENDDO
    ENDDO
    
    IF(PRESENT(i_am_in)) THEN
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
                   i_am_in(ir,itheta,iphi) = 1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts_in,x_inp,y_inp,z_inp, &
            psi_inp,TRIM(method), rank)
       
       IF(PRESENT(i_am_in)) THEN
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                              psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                              psi_sph_dphi(ir, itheta, iphi), rank )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
               .NOT.PRESENT(psi_sph_dth).AND. &
               .NOT. PRESENT(psi_sph_dphi) ) THEN
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt, ypt, zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta,iphi), rank=rank )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
          
       ELSE
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp,TRIM(method), psi_sph(ir,itheta, iphi), &
                           psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                           psi_sph_dphi(ir, itheta, iphi), rank )
                   ENDDO
                ENDDO
             ENDDO
             
          ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
               .NOT.PRESENT(psi_sph_dth).AND. &
               .NOT. PRESENT(psi_sph_dphi) ) THEN
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp,TRIM(method), psi_sph(ir,itheta,iphi), rank=rank )
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
          
       ENDIF
    ELSE
       
       CALL create_interpolant(numpts_in,x_inp,y_inp, z_inp, &
            psi_inp,TRIM(method) )
       
       IF(PRESENT(i_am_in)) THEN
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt, ypt, zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                              psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                              psi_sph_dphi(ir, itheta, iphi) )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
               .NOT.PRESENT(psi_sph_dth).AND. &
               .NOT. PRESENT(psi_sph_dphi) ) THEN
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
          
       ELSE
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      
                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp,TRIM(method), psi_sph(ir,itheta, iphi), &
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

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt, ypt, zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp,TRIM(method), psi_sph(ir,itheta,iphi) )
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
       ENDIF
    ENDIF
    
  END SUBROUTINE dcartesian2spherical3D
  
  !-------------------------------------------------------------------------!
  
  SUBROUTINE zcartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, i_am_in, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:), phi_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    INTEGER, INTENT(OUT),OPTIONAL       :: i_am_in(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dphi(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi 
    REAL(dp), ALLOCATABLE               :: x_inp(:)
    REAL(dp), ALLOCATABLE               :: y_inp(:)
    REAL(dp), ALLOCATABLE               :: z_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)
    
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2) *  dims_in(3)
    numpts_out = dims_out(1) *  dims_out(2) *  dims_out(3)
    
    ALLOCATE(x_inp(1:numpts_in))
    ALLOCATE(y_inp(1:numpts_in))
    ALLOCATE(z_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
    
    ii = 0
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             
             x_inp(ii) = x_ax(ix)
             y_inp(ii) = y_ax(iy)
             z_inp(ii) = z_ax(iz)
             psi_inp(ii) = psi_cart(ix, iy, iz)
             
          ENDDO
       ENDDO
    ENDDO
    
    ! Check if the new axis are within the cartesian domain.
    IF(PRESENT(i_am_in)) THEN
       DO iphi = 1, dims_out(3)
          DO itheta = 1, dims_out(2)
             DO ir = 1, dims_out(1)
                xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                zpt = r_ax(ir) * COS(theta_ax(itheta))
                IF ((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                     (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                     (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                   i_am_in(ir,itheta,iphi) = 1
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts_in,x_inp,y_inp, z_inp, &
            psi_inp,TRIM(method), rank)
       
       IF(PRESENT(i_am_in)) THEN    
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN

                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                              psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                              psi_sph_dphi(ir, itheta, iphi), rank )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
               .NOT.PRESENT(psi_sph_dth).AND. &
               .NOT. PRESENT(psi_sph_dphi) ) THEN
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta,iphi), rank=rank )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
          
       ELSE
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      
                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
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

                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
       ENDIF
       
    ELSE
       
       CALL create_interpolant(numpts_in,x_inp,y_inp,z_inp, &
            psi_inp,TRIM(method))
       
       IF(PRESENT(i_am_in)) THEN    
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                         
                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta, iphi), &
                              psi_sph_dr(ir,itheta, iphi), psi_sph_dth(ir,itheta, iphi),&
                              psi_sph_dphi(ir, itheta, iphi) )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
               .NOT.PRESENT(psi_sph_dth).AND. &
               .NOT. PRESENT(psi_sph_dphi) ) THEN
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      IF(i_am_in(ir,itheta,iphi).EQ.1) THEN

                         xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                         ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                         zpt = r_ax(ir) * COS(theta_ax(itheta))
                         
                         CALL interpolate(numpts_in,xpt,ypt,zpt, &
                              x_inp, y_inp, z_inp, &
                              psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                      ENDIF
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
          
       ELSE
          IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
             psi_sph_dr = ZERO
             psi_sph_dth = ZERO
             psi_sph_dphi = ZERO
             
             DO iphi = 1, dims_out(3)
                DO itheta = 1, dims_out(2)
                   DO ir = 1, dims_out(1)
                      
                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
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
                      
                      xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
                      ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
                      zpt = r_ax(ir) * COS(theta_ax(itheta))
                      
                      CALL interpolate(numpts_in,xpt,ypt,zpt, &
                           x_inp, y_inp, z_inp, &
                           psi_inp, TRIM(method), psi_sph(ir,itheta,iphi) )
                   ENDDO
                ENDDO
             ENDDO
             
          ELSE
             WRITE(*,*) 'You must input an array for getting the derivatives!'
             RETURN   
          ENDIF
       ENDIF
       
    ENDIF
    
  END SUBROUTINE zcartesian2spherical3D
  
  !---------------------------------------------------------------!
  
  SUBROUTINE cylindrical2spherical2D(psi_cyl,psi_sph, rho_ax, z_ax, &
       r_ax, theta_ax, method, psi_sph_dx, psi_sph_dy, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: rho_ax(:), z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL       :: rank
    
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
          IF(rpt.EQ.0.0_dp) THEN
             theta_inp(ii) = pi / 2.0_dp
          ELSE
             theta_inp(ii) = ACOS( z_ax(iz) / rpt )
          ENDIF
          psi_inp(ii) = psi_cyl(irho, iz)
       ENDDO
    ENDDO
    
    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts_in,r_inp,theta_inp,psi_inp,TRIM(method),rank)
       
       IF (PRESENT(psi_sph_dx).AND.PRESENT(psi_sph_dy)) THEN
          psi_sph_dx = ZERO
          psi_sph_dy = ZERO
          
          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), r_inp, theta_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta), &
                     psi_sph_dx(ir,itheta), psi_sph_dy(ir,itheta),rank)
             ENDDO
          ENDDO
          
       ELSE
          
          DO ir = 1, numrpts
             DO itheta = 1, numthetapts
                CALL interpolate(numpts,r_ax(ir), theta_ax(itheta), r_inp, theta_inp, &
                     psi_inp, TRIM(method), psi_sph(ir,itheta), rank=rank)
             ENDDO
          ENDDO
          
       ENDIF
       
    ELSE
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
       
    ENDIF
    
  END SUBROUTINE cylindrical2spherical2D
  
  
  !----------------------------------------------------------------!

END MODULE coords
