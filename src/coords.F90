MODULE coords

  USE constants
  USE gaussleg
  USE interp
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE

  PRIVATE

  PUBLIC           initialize_cartesian_boundary
  PUBLIC           initialize_cylindrical_boundary
  PUBLIC           get_cartesian_boundary
  PUBLIC           get_cylindrical_boundary
  PUBLIC           cartesian2spherical
  PUBLIC           cylindrical2spherical
  
  INTERFACE initialize_cartesian_boundary
     MODULE PROCEDURE initialize_cartesian_boundary2D
     MODULE PROCEDURE initialize_cartesian_boundary3D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cartesian_boundary3D_parallel
#endif
  END INTERFACE initialize_cartesian_boundary
  
  INTERFACE initialize_cylindrical_boundary
     MODULE PROCEDURE initialize_cylindrical_boundary2D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cylindrical_boundary2D_parallel
#endif
  END INTERFACE initialize_cylindrical_boundary
  
  INTERFACE get_cartesian_boundary
     MODULE PROCEDURE get_cartesian_boundary2D
     MODULE PROCEDURE get_cartesian_boundary3D
  END INTERFACE get_cartesian_boundary
  
  INTERFACE get_cylindrical_boundary
     MODULE PROCEDURE get_cylindrical_boundary2D
  END INTERFACE get_cylindrical_boundary
  
  INTERFACE cartesian2spherical
     MODULE PROCEDURE cartesian2spherical2D
     MODULE PROCEDURE dcartesian2spherical3D
     MODULE PROCEDURE zcartesian2spherical3D
  END INTERFACE cartesian2spherical
  
  INTERFACE cylindrical2spherical
     MODULE PROCEDURE cylindrical2spherical2D
  END INTERFACE cylindrical2spherical
  
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: rpts_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: theta_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: costheta_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: phi_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: theta_weights(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: pivots(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: phinodes(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: pivot_number(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: phinode_number(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_surface(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: surface_members(:)
  
  INTEGER                          :: numpts, numrpts
  INTEGER                          :: numthetapts, numphipts
  INTEGER                          :: numpivots
  INTEGER                          :: numphinodes
  INTEGER                          :: numsurfaceprocs
  REAL(dp), ALLOCATABLE            :: rpts_scatt(:)
  REAL(dp), ALLOCATABLE            :: theta_scatt(:)
  REAL(dp), ALLOCATABLE            :: phi_scatt(:) 
  COMPLEX(dp), ALLOCATABLE         :: psi_scatt(:)
  INTEGER, ALLOCATABLE             :: index_x1(:)
  INTEGER, ALLOCATABLE             :: index_x2(:)
  INTEGER, ALLOCATABLE             :: index_x3(:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph_dx(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph_dy(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi3D_sph(:,:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi3D_sph_dx(:,:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi3D_sph_dy(:,:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi3D_sph_dz(:,:,:)
  
CONTAINS
  
  !----------------------------------------------------!
  !  INITIALIZE CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE initialize_cartesian_boundary2D(x_ax, y_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, &
       maxpts, maxrpts, maxthetapts )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: x_ax(:)
    REAL(dp), INTENT(IN)      :: y_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(OUT)      :: maxpts
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts

    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: ix, iy, inum
    INTEGER                   :: ir, itheta
    
    !--------------------------------------------------!
    
    numpts = 0
    
    halolims(1,:) = (/ lbound(x_ax), ubound(x_ax) /)
    halolims(2,:) = (/ lbound(y_ax), ubound(y_ax) /)
    
    ! Assign max and min values for the axes
    minx = MINVAL(ABS(x_ax))
    maxx = MAXVAL(ABS(x_ax))
    
    miny = MINVAL(ABS(y_ax))
    maxy = MAXVAL(ABS(y_ax))
    
    minr = SQRT(minx**2 + miny**2)
    maxr = SQRT(maxy**2 + maxy**2)
            
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1)) 
       ALLOCATE(theta_weights(1))
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
       
       RETURN
    ENDIF
    
    mintheta =  twopi
    maxtheta = 0.0_dp
    
    ! Work out how many points will build the interpolant
    DO iy = halolims(2,1), halolims(2,2)
       DO ix = halolims(1,1), halolims(1,2)
          ! Calculate the point in spherical coordinates
          rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 )
          thetapt = pi - ATAN2(y_ax(iy) , -x_ax(ix) )
          
          ! Check the extend of the grid
          minr = MIN(rpt,minr)
          maxr = MAX(rpt,maxr)
          
          ! See if this points lies within the radius of
          ! influence of the boundary
          IF (ABS(rpt-Rs) .LE. radtol ) THEN
             numpts = numpts + 1
             mintheta = MIN(mintheta,thetapt)
             maxtheta = MAX(maxtheta,thetapt)
          ENDIF
       ENDDO
    ENDDO
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))  
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
       
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numpivots = 2 * lmax + 1
       numthetapts = numpivots
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       !numpivots = numthetapts
       
       ALLOCATE(pivots(1:numpivots))
       ALLOCATE(pivot_number(1:numthetapts))
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(theta_weights(1:numthetapts))
       ALLOCATE(psi2D_sph(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dx(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dy(1:numrpts,1:numthetapts))
       
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
!!$       CALL get_gauss_stuff(-1.0_dp, 1.0_dp, &
!!$            pivots, theta_weights)
!!$       
!!$       costheta_boundary = pivots
!!$       theta_boundary = ACOS(costheta_boundary)

       DO itheta = 1, numthetapts 
          theta_boundary(itheta) = REAL(itheta,dp) * twopi / REAL(numthetapts,dp)
       ENDDO
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
       ! Check limits on theta
       IF (mintheta .GT. MINVAL(pivots) .AND. &
            maxtheta .LT. MAXVAL(pivots)) THEN
          WRITE(*,*) 'Pivots extent is bigger than &
               & theta coordinate limits.'
          STOP
       ENDIF
       
       DO itheta = 1, numpivots
          pivot_number(itheta) = itheta
       ENDDO
       
       inum = 0       
       DO iy = halolims(2,1), halolims(2,2)
          DO ix = halolims(1,1), halolims(1,2)
             
             rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 )
             thetapt = pi - ATAN2(y_ax(iy) , -x_ax(ix))
             
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                inum = inum + 1
                rpts_scatt(inum)  = rpt
                theta_scatt(inum) = thetapt
                index_x1(inum) = ix
                index_x2(inum) = iy         
             ENDIF
          ENDDO
       ENDDO
       
    ENDIF
    
    maxpts  = numpts
    maxrpts = numrpts
    maxthetapts = numthetapts
    
  END SUBROUTINE initialize_cartesian_boundary2D
  
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !  INITIALIZE CARTESIAN BOUNDARY 3D
  !----------------------------------------------------!

  
  SUBROUTINE initialize_cartesian_boundary3D_serial(x_ax, y_ax, z_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, &
       maxpts, maxrpts, maxthetapts, maxphipts )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: x_ax(:)
    REAL(dp), INTENT(IN)      :: y_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(OUT)      :: maxpts
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts
    INTEGER, INTENT(OUT)      :: maxphipts
    
    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    REAL(dp)                  :: minphi, maxphi    
    REAL(dp)                  :: deltaphi
    INTEGER                   :: halolims(3,2)
    REAL(dp)                  :: rpt, thetapt, phipt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: ix, iy, iz, inum
    INTEGER                   :: ir, itheta, iphi
    
    !--------------------------------------------------!
    
    numpts = 0
   
    halolims(1,:) = (/ lbound(x_ax), ubound(x_ax) /)
    halolims(2,:) = (/ lbound(y_ax), ubound(y_ax) /)
    halolims(3,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minx = MINVAL(ABS(x_ax))
    maxx = MAXVAL(ABS(x_ax))

    miny = MINVAL(ABS(y_ax))
    maxy = MAXVAL(ABS(y_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minx**2 + miny**2 + minz**2)
    maxr = SQRT(maxx**2 + maxy**2 + maxz**2)
    
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(phi_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(index_x3(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(phinodes(1))       
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(phinode_number(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(psi3D_sph(1,1,1))
       ALLOCATE(psi3D_sph_dx(1,1,1))
       ALLOCATE(psi3D_sph_dy(1,1,1))
       
       RETURN
    ENDIF
    
    mintheta = pi
    maxtheta = 0.0_dp
    minphi   = twopi
    maxphi   = 0.0_dp
    
    ! Work out how many points will build the interpolant
    DO iz = halolims(3,1), halolims(3,2)
       DO iy = halolims(2,1), halolims(2,2)
          DO ix = halolims(1,1), halolims(1,2)
             ! Calculate the point in spherical coordinates
             rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             IF(rpt.EQ.0.0_dp) THEN
                thetapt = pi / 2.0_dp
             ELSE
                thetapt = ACOS(z_ax(iz) / rpt )
             ENDIF
             phipt = pi - ATAN2( y_ax(iy) , - x_ax(ix) )
             
             ! Check the extend of the grid
             minr = MIN(minr,rpt)
             maxr = MAX(maxr,rpt)
             
             ! See if this points lies within the radius of
             ! influence of the boundary
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                numpts = numpts + 1
                mintheta = MIN(mintheta,thetapt)
                maxtheta = MAX(maxtheta,thetapt)
                minphi = MIN(minphi,phipt)
                maxphi = MAX(maxphi,phipt)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))
       ALLOCATE(phi_scatt(1))
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(phinodes(1))       
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(phinode_number(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(psi3D_sph(1,1,1))
       ALLOCATE(psi3D_sph_dx(1,1,1))
       ALLOCATE(psi3D_sph_dy(1,1,1))
       
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(phi_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       ALLOCATE(index_x3(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numpivots = 2 * lmax + 1
       numphinodes = 2 * lmax + 1
       numthetapts = numpivots
       numphipts = numphinodes
       deltaphi = twopi / REAL(numphinodes, dp)
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       !numphipts = INT( (maxphi - minphi) / deltaphi )
       
       ALLOCATE(pivots(1:numpivots))
       ALLOCATE(phinodes(1:numphinodes))
       ALLOCATE(pivot_number(1:numthetapts))
       ALLOCATE(phinode_number(1:numphipts))
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(phi_boundary(1:numphipts))
       ALLOCATE(theta_weights(1:numthetapts))
       ALLOCATE(psi3D_sph(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dx(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dy(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dz(1:numrpts,1:numthetapts,1:numphipts))
       
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
       CALL get_gauss_stuff(-1.0_dp, 1.0_dp, pivots, &
            theta_weights)

       costheta_boundary = pivots
       theta_boundary = ACOS(costheta_boundary)
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
       ! Check limits on theta
       IF (mintheta .GT. MINVAL(ACOS(pivots)) .AND. &
            maxtheta .LT. MAXVAL(ACOS(pivots))) THEN
          WRITE(*,*) 'Pivots extent is bigger than &
               & theta coordinate limits.'
          STOP
       ENDIF
       
       DO iphi = 1, numphinodes
          phinodes(iphi) = REAL( iphi,dp ) * deltaphi
       ENDDO
       
       phi_boundary = phinodes
       
       ! Check limits on phi
       IF (minphi .GT. MINVAL(phinodes) .AND. &
            maxphi .LT. MAXVAL(phinodes)) THEN
          WRITE(*,*) 'Phi nodes extent is bigger than &
               & phi coordinate limits.'
          STOP
       ENDIF
       
       DO itheta = 1, numpivots
          pivot_number(itheta) = itheta
       ENDDO
       
       DO iphi = 1, numphinodes
          phinode_number(iphi) = iphi
       ENDDO
       
       inum = 0
       DO iz = halolims(3,1), halolims(3,2)
          DO iy = halolims(2,1), halolims(2,2)
             DO ix = halolims(1,1), halolims(1,2)
                
                rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
                IF(rpt.EQ.0.0_dp) THEN
                   thetapt = 0.0_dp
                ELSE
                   thetapt = ACOS( z_ax(iz) / rpt)
                ENDIF
                phipt = pi - ATAN2( y_ax(iy),-x_ax(ix) )
                
                IF (ABS(rpt-Rs) .LE. radtol ) THEN
                   inum = inum + 1
                   rpts_scatt(inum)  = rpt
                   theta_scatt(inum) = thetapt
                   phi_scatt(inum) = phipt
                   index_x1(inum) = ix
                   index_x2(inum) = iy
                   index_x3(inum) = iz
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
    ENDIF
    
    maxpts = numpts
    maxrpts = numrpts
    maxthetapts = numthetapts
    maxphipts = numphipts
    
  END SUBROUTINE initialize_cartesian_boundary3D_serial
  
  !-----------------------------------------------------------------!
#if _COM_MPI
  SUBROUTINE initialize_cartesian_boundary3D_parallel(x_ax, y_ax, z_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, rank, size, comm, &
       maxpts, maxrpts, maxthetapts, maxphipts, maxsurfaceprocs, newcomm )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: x_ax(:)
    REAL(dp), INTENT(IN)      :: y_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(IN)       :: rank, size
    INTEGER, INTENT(IN)       :: comm
    INTEGER, INTENT(OUT)      :: maxpts
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts
    INTEGER, INTENT(OUT)      :: maxphipts
    INTEGER, INTENT(OUT)      :: maxsurfaceprocs
    INTEGER, INTENT(OUT)      :: newcomm

    INTEGER, ALLOCATABLE      :: i_am_surface_local(:)
    INTEGER                   :: simgroup, surfacegroup
    INTEGER                   :: ierror
    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    REAL(dp)                  :: minphi, maxphi    
    REAL(dp)                  :: deltaphi
    INTEGER                   :: halolims(3,2)
    REAL(dp)                  :: rpt, thetapt, phipt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: ix, iy, iz, inum
    INTEGER                   :: ir, itheta, iphi, ii
    REAL(dp)                  :: minrr, maxrr
    !--------------------------------------------------!
    
    numpts = 0
   
    halolims(1,:) = (/ lbound(x_ax), ubound(x_ax) /)
    halolims(2,:) = (/ lbound(y_ax), ubound(y_ax) /)
    halolims(3,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minx = MINVAL(ABS(x_ax))
    maxx = MAXVAL(ABS(x_ax))

    miny = MINVAL(ABS(y_ax))
    maxy = MAXVAL(ABS(y_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minx**2 + miny**2 + minz**2)
    maxr = SQRT(maxx**2 + maxy**2 + maxz**2)
    
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(phi_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(index_x3(1))
       ALLOCATE(pivots(1))
       ALLOCATE(phinodes(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(phinode_number(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(psi3D_sph(1,1,1))
       ALLOCATE(psi3D_sph_dx(1,1,1))
       ALLOCATE(psi3D_sph_dy(1,1,1))
       
       RETURN
    ENDIF
    
    mintheta = pi
    maxtheta = 0.0_dp
    minphi   = twopi
    maxphi   = 0.0_dp

    maxrr = 0.0_dp
    minrr = Huge(minrr)

    ! Work out how many points will build the interpolant
    DO iz = halolims(3,1), halolims(3,2)
       DO iy = halolims(2,1), halolims(2,2)
          DO ix = halolims(1,1), halolims(1,2)
             ! Calculate the point in spherical coordinates
             rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             IF(rpt.EQ.0.0_dp) THEN
                thetapt = pi / 2.0_dp
             ELSE
                thetapt = ACOS(z_ax(iz) / rpt )
             ENDIF
             phipt = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
             
             ! Check the extend of the grid
             minr = MIN(minr,rpt)
             maxr = MAX(maxr,rpt)
             
             ! See if this points lies within the radius of
             ! influence of the boundary
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                numpts = numpts + 1
                mintheta = MIN(mintheta,thetapt)
                maxtheta = MAX(maxtheta,thetapt)

                minphi = MIN(minphi,phipt)
                maxphi = MAX(maxphi,phipt)

                minrr = MIN(minrr,rpt)
                maxrr = MAX(maxrr,rpt)
                 
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    IF((ANY(y_ax.GT.0.0_dp)).AND.(ANY(y_ax.LT.0.0_dp)) &
         .AND. (ALL(x_ax.GT.0))) THEN
       minphi = pi - ATAN2(MINVAL(y_ax),-MINVAL(x_ax))
       maxphi = pi - ATAN2(MAXVAL(y_ax),-MINVAL(x_ax))
    ENDIF
    
    ! Create communicator and surface members
    !-----------------------------------------!
    ALLOCATE(i_am_surface_local(0:size-1))
    ALLOCATE(i_am_surface(0:size-1))
    i_am_surface_local = 0
    i_am_surface = 0
    
    IF(numpts.NE.0) &
         i_am_surface_local(rank) = 1
    
    ! Communicate to all processors if they are
    ! at the surface
    CALL MPI_ALLREDUCE(i_am_surface_local, i_am_surface, size, &
         MPI_INTEGER, MPI_SUM, comm, ierror )
    
    DO inum = 0, size - 1
       IF(i_am_surface(inum).EQ.1) &
            numsurfaceprocs = numsurfaceprocs + 1
    ENDDO
    
    ALLOCATE(surface_members(0:numsurfaceprocs-1))
    ii = -1
    DO inum = 0, size-1
       IF(i_am_surface(inum).EQ.1) THEN
          ii = ii + 1
          surface_members(ii) = inum
       ENDIF
    ENDDO
    
    ! Create a global group
    CALL MPI_COMM_GROUP(comm, simgroup, ierror)
    
    CALL MPI_GROUP_INCL(simgroup, numsurfaceprocs, surface_members, &
         surfacegroup, ierror)
    ! Actually, this line creates the communicator
    CALL MPI_COMM_CREATE(comm, surfacegroup, newcomm, ierror)
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))
       ALLOCATE(phi_scatt(1))
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(pivots(1))
       ALLOCATE(phinodes(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(phinode_number(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(psi3D_sph(1,1,1))
       ALLOCATE(psi3D_sph_dx(1,1,1))
       ALLOCATE(psi3D_sph_dy(1,1,1))
       
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(phi_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       ALLOCATE(index_x3(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numpivots = 2 * lmax + 1
       numphinodes = 2 * lmax + 1
       deltaphi = twopi / REAL(numphinodes, dp)
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       !numphipts = INT( (maxphi - minphi) / deltaphi )
       
       ALLOCATE(pivots(1:numpivots))
       ALLOCATE(theta_weights(1:numpivots))
       ALLOCATE(phinodes(1:numphinodes))
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(psi3D_sph(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dx(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dy(1:numrpts,1:numthetapts,1:numphipts))
       ALLOCATE(psi3D_sph_dz(1:numrpts,1:numthetapts,1:numphipts))
       
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
       CALL get_gauss_stuff(-1.0_dp, 1.0_dp, pivots, &
            theta_weights)
       
       ! Find out how many pivots lie in our processor 
       numthetapts = 0
       DO itheta = 1, numpivots
          IF(ACOS(pivots(itheta)).GE. mintheta .AND. &
               ACOS(pivots(itheta)).LE. maxtheta) THEN
             numthetapts = numthetapts + 1
          ENDIF
       ENDDO
       
       DO iphi = 1, numphinodes
          phinodes(iphi) = REAL( iphi,dp ) * deltaphi
       ENDDO
       
       ! Find out how many phi nodes lie in our processor 
       ! Remenber, if we are in the 4th quadrant of the imagninary
       ! plane, minphi and max phi are reversed.
       numphipts = 0
       DO iphi = 1, numphinodes
          IF(minphi.LT.maxphi) THEN
             IF( phinodes(iphi).GE.minphi .AND. &
                  phinodes(iphi).LE.maxphi) THEN
                numphipts = numphipts + 1
             ENDIF
          ELSE
             IF( phinodes(iphi).GE.minphi .OR. &
                  phinodes(iphi).LE.maxphi) THEN
                numphipts = numphipts + 1
             ENDIF
          ENDIF
       ENDDO
       
       ! Allocate arrays once we know the length
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(phi_boundary(1:numphipts))
       ALLOCATE(pivot_number(1:numthetapts))
       ALLOCATE(phinode_number(1:numphipts))
       
       ! Now, assignate those theta points
       inum = 0
       DO itheta = 1, numpivots
          IF(ACOS(pivots(itheta)).GE.mintheta .AND. &
               ACOS(pivots(itheta)).LE. maxtheta) THEN
             inum = inum + 1
             costheta_boundary(inum) = pivots(itheta)
             pivot_number(inum) = itheta             
          ENDIF
       ENDDO
       
       theta_boundary = ACOS(costheta_boundary)
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
       ! Don't forget to assignate phi points
       inum = 0
       DO iphi = 1, numphinodes
          IF(minphi.LT.maxphi) THEN
             IF(phinodes(iphi).GE.minphi .AND. &
                  phinodes(iphi).LE. maxphi) THEN
                inum = inum + 1
                phi_boundary(inum) = phinodes(iphi)
                phinode_number(inum) = iphi
             ENDIF
          ELSE
             IF(phinodes(iphi).GE.minphi .OR. &
                  phinodes(iphi).LE.maxphi) THEN
                inum = inum + 1
                phi_boundary(inum) = phinodes(iphi)
                phinode_number(inum) = iphi
             ENDIF
          ENDIF
       ENDDO
       
!!$       write(*,*) 'x',rank,MINVAL(x_ax),MAXVAL(x_ax)
!!$       !write(*,*) 'y',rank,MINVAL(y_ax),MAXVAL(y_ax)
!!$
!!$       write(*,*) rank,'min/max phi',minphi,maxphi
!!$       write(*,*) rank,'num',numphipts
!!$       
!!$       write(*,*) rank,'min/max theta',mintheta,maxtheta
!!$       write(*,*) rank,'num',numthetapts
!!$       
!!$       write(*,*) rank,'min/max r',minrr,maxrr
!!$       write(*,*) rank,'num',numrpts
!!$       
!!$       
!!$       if(rank.eq.2) then
!!$          write(*,*) 'theta', pivot_number
!!$          write(*,*) 'phi node', phi_boundary
!!$       endif
!!$       
       ! Get points to participate in the interpolant
       inum = 0   
       DO iz = halolims(3,1), halolims(3,2)
          DO iy = halolims(2,1), halolims(2,2)
             DO ix = halolims(1,1), halolims(1,2)
                
                rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
                IF(rpt.EQ.0.0_dp) THEN
                   thetapt = 0.0_dp
                ELSE
                   thetapt = ACOS( z_ax(iz) / rpt)
                ENDIF
                phipt = pi - ATAN2( y_ax(iy),-x_ax(ix) )
                
                IF (ABS(rpt-Rs) .LE. radtol ) THEN
                   inum = inum + 1
                   rpts_scatt(inum)  = rpt
                   theta_scatt(inum) = thetapt
                   phi_scatt(inum) = phipt
                   index_x1(inum) = ix
                   index_x2(inum) = iy
                   index_x3(inum) = iz
                   
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
    ENDIF
    
    maxpts = numpts
    maxrpts = numrpts
    maxthetapts = numthetapts
    maxphipts = numphipts
    maxsurfaceprocs = numsurfaceprocs
    
  END SUBROUTINE initialize_cartesian_boundary3D_parallel
#endif  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  INITIALIZE CYLINDRICAL BOUNDARY 2D
  !----------------------------------------------------!
  
  
  SUBROUTINE initialize_cylindrical_boundary2D_serial(rho_ax, z_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, &
       maxpts, maxrpts, maxthetapts )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: rho_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(OUT)      :: maxpts
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts
    
    REAL(dp)                  :: minrho, maxrho
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: irho, iz, inum
    INTEGER                   :: ir, itheta
    
    !--------------------------------------------------!
    
    numpts = 0
    
    halolims(1,:) = (/ lbound(rho_ax), ubound(rho_ax) /)
    halolims(2,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minrho = MINVAL(ABS(rho_ax))
    maxrho = MAXVAL(ABS(rho_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minrho**2 + minz**2)
    maxr = SQRT(maxrho**2 + maxz**2)
            
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
       
       RETURN
    ENDIF
    
    mintheta = pi
    maxtheta = 0.0_dp
    
    ! Work out how many points will build the interpolant
    DO iz = halolims(2,1), halolims(2,2)
       DO irho = halolims(1,1), halolims(1,2)
          ! Calculate the point in spherical coordinates
          rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
          IF (rpt.EQ.0.0_dp) THEN
             thetapt = pi / 2.0_dp
          ELSE
             thetapt = ACOS(z_ax(iz) / rpt )
          ENDIF
          ! Check the extend of the grid
          minr = MIN(rpt,minr)
          maxr = MAX(rpt,maxr)

          ! See if this points lies within the radius of
          ! influence of the boundary
          IF (ABS(rpt-Rs) .LE. radtol ) THEN
             numpts = numpts + 1
             mintheta = MIN(mintheta,thetapt)
             maxtheta = MAX(maxtheta,thetapt)
          ENDIF
       ENDDO
    ENDDO
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))  
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(pivots(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
       
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numpivots = 2 * lmax + 1
       numthetapts = numpivots
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       
       ALLOCATE(pivots(1:numpivots))
       ALLOCATE(pivot_number(1:numpivots))
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(theta_weights(1:numpivots))
       ALLOCATE(psi2D_sph(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dx(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dy(1:numrpts,1:numthetapts))
       
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
       CALL get_gauss_stuff(-1.0_dp, 1.0_dp, &
            pivots, theta_weights)
       
       costheta_boundary = pivots
       theta_boundary = ACOS(costheta_boundary)
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
       ! Check limits on theta
       IF (mintheta .GT. MINVAL(ACOS(pivots)) .AND. &
            maxtheta .LT. MAXVAL(ACOS(pivots))) THEN
          WRITE(*,*) 'Pivots extent is bigger than &
               & theta coordinate limits.'
          STOP
       ENDIF
       
       DO itheta = 1, numpivots
          pivot_number(itheta) = itheta
       ENDDO
       
       inum = 0
       DO iz = halolims(2,1), halolims(2,2)
          DO irho = halolims(1,1), halolims(1,2)
             
             rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
             IF(rpt.EQ.0.0_dp) THEN
                thetapt = pi / 2.0_dp
             ELSE
                thetapt = ACOS( z_ax(iz) / rpt )
             ENDIF
             
             IF (ABS(rpt-Rs) .LE. radtol) THEN
                inum = inum + 1
                rpts_scatt(inum)  = rpt
                theta_scatt(inum) = thetapt
                index_x1(inum) = irho
                index_x2(inum) = iz
                
             ENDIF
          ENDDO
       ENDDO
       
    ENDIF
    
    maxpts = numpts
    maxrpts = numrpts
    maxthetapts = numthetapts
    
  END SUBROUTINE initialize_cylindrical_boundary2D_serial

  !*****************************************************************!
  !*****************************************************************!
#if _COM_MPI
  
  SUBROUTINE initialize_cylindrical_boundary2D_parallel(rho_ax, z_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, rank, size, comm, &
       maxpts, maxrpts, maxthetapts, maxsurfaceprocs, newcomm )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: rho_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rs
    REAL(dp), INTENT(IN)      :: radtol
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(IN)       :: rank, size
    INTEGER, INTENT(IN)       :: comm
    INTEGER, INTENT(OUT)      :: maxpts
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts
    INTEGER, INTENT(OUT)      :: maxsurfaceprocs
    INTEGER, INTENT(OUT)      :: newcomm
    
    INTEGER, ALLOCATABLE      :: i_am_surface_local(:)
    INTEGER                   :: simgroup, surfacegroup
    INTEGER                   :: ierror
    REAL(dp)                  :: minrho, maxrho
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: irho, iz, inum
    INTEGER                   :: ir, itheta, ii
    
    !--------------------------------------------------!

    numpts = 0
    
    halolims(1,:) = (/ lbound(rho_ax), ubound(rho_ax) /)
    halolims(2,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minrho = MINVAL(ABS(rho_ax))
    maxrho = MAXVAL(ABS(rho_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minrho**2 + minz**2)
    maxr = SQRT(maxrho**2 + maxz**2)
            
    IF ( (Rs.LT.minr).AND.(Rs.GT.maxr) ) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
      
       RETURN
    ENDIF
    
    mintheta = pi
    maxtheta = 0.0_dp
    
    ! Work out how many points will build the interpolant
    DO iz = halolims(2,1), halolims(2,2)
       DO irho = halolims(1,1), halolims(1,2)
          ! Calculate the point in spherical coordinates
          rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
          IF (rpt.EQ.0.0_dp) THEN
             thetapt = pi / 2.0_dp
          ELSE
             thetapt = ACOS( z_ax(iz) / rpt )
          ENDIF
          ! Check the extend of the grid
          minr = MIN(rpt,minr)
          maxr = MAX(rpt,maxr)
          
          ! See if this points lies within the radius of
          ! influence of the boundary
          IF (ABS(rpt-Rs) .LE. radtol ) THEN
             numpts = numpts + 1
             mintheta = MIN(mintheta,thetapt)
             maxtheta = MAX(maxtheta,thetapt)
          ENDIF
       ENDDO
    ENDDO
    
    ! Create communicator and surface members
    !-----------------------------------------!
    ALLOCATE(i_am_surface_local(0:size-1))
    ALLOCATE(i_am_surface(0:size-1))
    i_am_surface_local = 0
    i_am_surface = 0
    
    IF(numpts.NE.0) &
         i_am_surface_local(rank) = 1
    
    ! Communicate to all processors if they are
    ! at the surface
    CALL MPI_ALLREDUCE(i_am_surface_local, i_am_surface, size, &
         MPI_INTEGER, MPI_SUM, comm, ierror )
    
    DO inum = 0, size - 1
       IF(i_am_surface(inum).EQ.1) &
            numsurfaceprocs = numsurfaceprocs + 1
    ENDDO
    
    ALLOCATE(surface_members(0:numsurfaceprocs-1))
    ii = -1
    DO inum = 0, size-1
       IF(i_am_surface(inum).EQ.1) THEN
          ii = ii + 1
          surface_members(ii) = inum
       ENDIF
    ENDDO
    
    ! Create a global group
    CALL MPI_COMM_GROUP(comm, simgroup, ierror)
    
    CALL MPI_GROUP_INCL(simgroup, numsurfaceprocs, surface_members, &
         surfacegroup, ierror)
    ! Actually, this line creates the communicator
    CALL MPI_COMM_CREATE(comm, surfacegroup, newcomm, ierror)
    
    ! Create grids for interpolation
    IF(numpts .EQ. 0) THEN
       ALLOCATE(rpts_scatt(1))
       ALLOCATE(theta_scatt(1))  
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(pivots(1))
       ALLOCATE(pivot_number(1))
       ALLOCATE(psi2D_sph(1,1))
       ALLOCATE(psi2D_sph_dx(1,1))
       ALLOCATE(psi2D_sph_dy(1,1))
       
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numpivots = 2 * lmax + 1
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       
       ALLOCATE(pivots(1:numpivots))
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(theta_weights(1:numpivots))
       ALLOCATE(psi2D_sph(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dx(1:numrpts,1:numthetapts))
       ALLOCATE(psi2D_sph_dy(1:numrpts,1:numthetapts))
       
       ! Create desired axis for interpolation
       ! Radial axis
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
       ! Now theta axis
       CALL get_gauss_stuff(-1.0_dp, 1.0_dp, &
            pivots, theta_weights)
       
       ! Find out how many pivot points lie in our processor 
       numthetapts = 0
       DO itheta = 1, numpivots
          IF(ACOS(pivots(itheta)).GE. mintheta .AND. &
               ACOS(pivots(itheta)).LE. maxtheta) THEN
             numthetapts = numthetapts + 1
          ENDIF
       ENDDO
       
       ! Allocate arrays once we know the length
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(pivot_number(1:numthetapts))
       
       ! Now, assignate those points
       inum = 0
       DO itheta = 1, numpivots
          IF(ACOS(pivots(itheta)).GE.mintheta .AND. &
               ACOS(pivots(itheta)).LE. maxtheta) THEN
             inum = inum + 1
             costheta_boundary(inum) = pivots(itheta)
             pivot_number(inum) = itheta             
          ENDIF
       ENDDO
       
       theta_boundary = ACOS(costheta_boundary)
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
!!$       ! Check limits on theta
!!$       IF (mintheta .GT. min(ACOS(pivots)) .AND. &
!!$            maxtheta .LT. max(ACOS(pivots))) THEN
!!$          WRITE(*,*) 'Pivots extent is bigger than &
!!$               & theta coordinate limits.'
!!$          STOP
!!$       ENDIF
       
       ! Now save theta points from the original grid
       inum = 0
       DO iz = halolims(2,1), halolims(2,2)
          DO irho = halolims(1,1), halolims(1,2)
             
             rpt = SQRT( rho_ax(irho)**2 + z_ax(iz)**2 )
             IF(rpt.EQ.0.0_dp) THEN
                thetapt = pi / 2.0_dp
             ELSE
                thetapt = ACOS( z_ax(iz) / rpt )
             ENDIF
             
             IF (ABS(rpt-Rs) .LE. radtol) THEN
                inum = inum + 1
                rpts_scatt(inum)  = rpt
                theta_scatt(inum) = thetapt
                index_x1(inum) = irho
                index_x2(inum) = iz
                
             ENDIF
          ENDDO
       ENDDO
       
    ENDIF
    
    ! Assign points to return
    maxpts = numpts
    maxrpts = numrpts
    maxthetapts = numthetapts
    maxsurfaceprocs = numsurfaceprocs
    
  END SUBROUTINE initialize_cylindrical_boundary2D_parallel
#endif
  !-----------------------------------------------------------------!

  !-----------------------------------------------------------------!

  !-----------------------------------------------------------------!

  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE get_cartesian_boundary2D(psi_cart,psi2D_sph, &
       psi2D_sph_dx, psi2D_sph_dy, method)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dy(:, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    
    INTEGER                           :: inum, ir, itheta
    
    psi_scatt = ZERO
    psi2D_sph = ZERO
    psi2D_sph_dx = ZERO
    psi2D_sph_dy = ZERO
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum))       
    ENDDO
    
    CALL create_interpolant(numpts,rpts_scatt,theta_scatt,psi_scatt,TRIM(method))
    
    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          CALL interpolate(numpts,rpts_boundary(ir), theta_boundary(itheta), &
               rpts_scatt, theta_scatt, psi_scatt, TRIM(method), &
               psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), psi2D_sph_dy(ir,itheta))
       ENDDO
    ENDDO
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, psi_scatt)
    
  END SUBROUTINE get_cartesian_boundary2D
  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 3D
  !----------------------------------------------------!
  
  SUBROUTINE get_cartesian_boundary3D(psi_cart,psi3D_sph, &
       psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, method)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dx(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dy(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dz(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    
    INTEGER                           :: inum, ir, itheta, iphi
    
    psi_scatt = ZERO
    psi3D_sph = ZERO
    psi3D_sph_dx = ZERO
    psi3D_sph_dy = ZERO
    psi3D_sph_dz = ZERO
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum),index_x3(inum))       
    ENDDO

    CALL create_interpolant(numpts,rpts_scatt,theta_scatt,phi_scatt, &
         psi_scatt,TRIM(method))
 
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             CALL interpolate(numpts,rpts_boundary(ir), theta_boundary(itheta), &
                  phi_boundary(iphi),rpts_scatt, theta_scatt, phi_scatt, &
                  psi_scatt, TRIM(method), psi3D_sph(ir,itheta,iphi), &
                  psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                  psi3D_sph_dz(ir,itheta,iphi))
          ENDDO
       ENDDO
    ENDDO
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, phi_scatt, &
         psi_scatt)
 
  END SUBROUTINE get_cartesian_boundary3D
  
  !----------------------------------------------------!
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CYLINDRICAL BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE get_cylindrical_boundary2D(psi_cyl,psi2D_sph, &
       psi2D_sph_dx, psi2D_sph_dy, method)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cyl(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dy(:, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    
    INTEGER                           :: inum, ir, itheta
    
    psi_scatt = ZERO
    psi2D_sph = ZERO
    psi2D_sph_dx = ZERO
    psi2D_sph_dy = ZERO
    
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cyl(index_x1(inum),index_x2(inum))       
    ENDDO
    
    CALL create_interpolant(numpts,rpts_scatt,theta_scatt,psi_scatt,TRIM(method))

    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          CALL interpolate(numpts,rpts_boundary(ir), theta_boundary(itheta), &
               rpts_scatt, theta_scatt, psi_scatt, TRIM(method), &
               psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), psi2D_sph_dy(ir,itheta))
       ENDDO
    ENDDO
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, psi_scatt)
    
  END SUBROUTINE get_cylindrical_boundary2D
  
  !----------------------------------------------------------------!
  
  SUBROUTINE cartesian2spherical2D(psi_cart,psi_sph, x_ax, y_ax, &
       r_ax, theta_ax, method, psi_sph_dr, psi_sph_dth)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)             :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)            :: psi_sph(:, :)
    REAL(dp), INTENT(IN)                :: x_ax(:), y_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:), theta_ax(:)
    CHARACTER(LEN=*), INTENT(IN)        :: method
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dr(:, :)
    COMPLEX(dp), INTENT(OUT), OPTIONAL  :: psi_sph_dth(:, :)
    
    INTEGER                             :: ix, iy, ii
    INTEGER                             :: ir, itheta, inum
    INTEGER, DIMENSION(2)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp), ALLOCATABLE               :: r_inp(:)
    REAL(dp), ALLOCATABLE               :: theta_inp(:)
    COMPLEX(dp), ALLOCATABLE            :: psi_inp(:)
    
    
    dims_in = SHAPE(psi_cart)
    dims_out = SHAPE(psi_sph)
    
    numpts_in = dims_in(1) *  dims_in(2)
    numpts_out = dims_out(1) *  dims_out(2)
    
    ALLOCATE(r_inp(1:numpts_in))
    ALLOCATE(theta_inp(1:numpts_in))
    ALLOCATE(psi_inp(1:numpts_in))
    
    psi_inp = ZERO
    psi_sph = ZERO
    
    ii = 0
    
    DO iy = 1, dims_in(2)
       DO ix = 1, dims_in(1)
          ii = ii + 1
          rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 )
          
          r_inp(ii) = rpt
          theta_inp(ii) = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
          psi_inp(ii) = psi_cart(ix, iy)
       ENDDO
    ENDDO
    
    CALL create_interpolant(numpts_in,r_inp,theta_inp,psi_inp,TRIM(method))
    
    IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth)) THEN
       psi_sph_dr = ZERO
       psi_sph_dth = ZERO
       
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), &
                  r_inp, theta_inp, &
                  psi_inp, TRIM(method), psi_sph(ir,itheta), &
                  psi_sph_dr(ir,itheta), psi_sph_dth(ir,itheta) )
          ENDDO
       ENDDO
       
       
    ELSEIF (.NOT.PRESENT(psi_sph_dr).AND. &
         .NOT.PRESENT(psi_sph_dth) ) THEN
       
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), &
                  r_inp, theta_inp, &
                  psi_inp, TRIM(method), psi_sph(ir,itheta) )
          ENDDO
       ENDDO
       
    ELSE
       WRITE(*,*) 'You must input and array for get the derivatives!'
       RETURN   
    ENDIF
    
    
  END SUBROUTINE cartesian2spherical2D
  
  !----------------------------------------------------------------!
  
  !----------------------------------------------------------------!
  !----------------------------------------------------------------!
  
  SUBROUTINE dcartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, i_am_in, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi)
    
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
    
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi 
    REAL(dp), ALLOCATABLE               :: r_inp(:)
    REAL(dp), ALLOCATABLE               :: theta_inp(:)
    REAL(dp), ALLOCATABLE               :: phi_inp(:)
    REAL(dp), ALLOCATABLE               :: psi_inp(:)
    
    
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
    
    minr = HUGE(minr) !1.0E9_dp
    maxr = 0.0_dp
    
    mintheta = pi
    maxtheta = 0.0_dp
    
    minphi = pi
    maxphi = -pi
    
    ii = 0
    DO iz = 1, dims_in(3)
       DO iy = 1, dims_in(2)
          DO ix = 1, dims_in(1)
             ii = ii + 1
             rpt = SQRT(x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
             
             r_inp(ii) = rpt
             IF(rpt.EQ.0.0_dp) THEN
                theta_inp(ii) = pi / 2.0_dp
             ELSE
                theta_inp(ii) = ACOS( z_ax(iz) / rpt )
             ENDIF
             phi_inp(ii) = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
             psi_inp(ii) = psi_cart(ix, iy, iz)
             
             minr = MIN(minr,r_inp(ii))
             mintheta = MIN(mintheta,theta_inp(ii))
             minphi = MIN(minphi,phi_inp(ii))
             
             maxr = MAX(maxr,r_inp(ii))
             maxtheta = MAX(maxtheta,theta_inp(ii))
             maxphi = MAX(maxphi,phi_inp(ii))
          ENDDO
       ENDDO
    ENDDO
    
    CALL create_interpolant(numpts_in,r_inp,theta_inp, phi_inp, &
         psi_inp,TRIM(method))
    
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))
             IF ((r_ax(ir).GE.minr) .AND. (r_ax(ir).LE.maxr) .AND. &
                  (theta_ax(itheta).GE.mintheta) .AND. &
                  (theta_ax(itheta).LE.maxtheta) .AND. &
                  (phi_ax(iphi).GE.minphi) .AND. (phi_ax(iphi).LE.maxphi) .AND. &
                  (xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in(ir,itheta,iphi) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    IF(PRESENT(i_am_in)) THEN
       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO
          
          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                      CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                           r_inp, theta_inp, phi_inp, &
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
                      CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                           r_inp, theta_inp, phi_inp, &
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
                   CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                        r_inp, theta_inp, phi_inp, &
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
                   CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                        r_inp, theta_inp, phi_inp, &
                        psi_inp,TRIM(method), psi_sph(ir,itheta,iphi) )
                ENDDO
             ENDDO
          ENDDO
          
       ELSE
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN   
       ENDIF
       
    ENDIF
    
  END SUBROUTINE dcartesian2spherical3D

  !-------------------------------------------------------------------------!
  
  SUBROUTINE zcartesian2spherical3D(psi_cart,psi_sph, x_ax, y_ax, z_ax, &
       r_ax, theta_ax, phi_ax, method, i_am_in, psi_sph_dr, psi_sph_dth, &
       psi_sph_dphi)
    
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
    
    INTEGER                             :: ix, iy, iz, ii
    INTEGER                             :: ir, itheta, iphi, inum
    INTEGER, DIMENSION(3)               :: dims_in, dims_out
    INTEGER                             :: numpts_in, numpts_out
    REAL(dp)                            :: rpt
    REAL(dp)                            :: xpt, ypt, zpt
    REAL(dp)                            :: minr, maxr
    REAL(dp)                            :: mintheta, maxtheta
    REAL(dp)                            :: minphi, maxphi 
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
             IF(rpt.EQ.0.0_dp) THEN
                theta_inp(ii) = pi / 2.0_dp
             ELSE
                theta_inp(ii) = ACOS( z_ax(iz) / rpt )
             ENDIF
             phi_inp(ii) = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
             psi_inp(ii) = psi_cart(ix, iy, iz)
          ENDDO
       ENDDO
    ENDDO
    
    CALL create_interpolant(numpts_in,r_inp,theta_inp, phi_inp, &
         psi_inp,TRIM(method))
    
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, dims_out(3)
       DO itheta = 1, dims_out(2)
          DO ir = 1, dims_out(1)
             xpt = r_ax(ir) * SIN(theta_ax(itheta)) * COS(phi_ax(iphi))
             ypt = r_ax(ir) * SIN(theta_ax(itheta)) * SIN(phi_ax(iphi))
             zpt = r_ax(ir) * COS(theta_ax(itheta))
             IF ((r_ax(ir).GE.minr) .AND. (r_ax(ir).LE.maxr) .AND. &
                  (theta_ax(itheta).GE.mintheta) .AND. &
                  (theta_ax(itheta).LE.maxtheta) .AND. &
                  (phi_ax(iphi).GE.minphi) .AND. (phi_ax(iphi).LE.maxphi) .AND. &
                  (xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in(ir,itheta,iphi) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    IF(PRESENT(i_am_in)) THEN    
       IF (PRESENT(psi_sph_dr).AND.PRESENT(psi_sph_dth).AND.PRESENT(psi_sph_dphi) ) THEN
          psi_sph_dr = ZERO
          psi_sph_dth = ZERO
          psi_sph_dphi = ZERO
          
          DO iphi = 1, dims_out(3)
             DO itheta = 1, dims_out(2)
                DO ir = 1, dims_out(1)
                   IF(i_am_in(ir,itheta,iphi).EQ.1) THEN
                      CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                           r_inp, theta_inp, phi_inp, &
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
                      CALL interpolate(numpts_in,r_ax(ir), theta_ax(itheta), phi_ax(iphi), &
                           r_inp, theta_inp, phi_inp, &
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
          WRITE(*,*) 'You must input an array for getting the derivatives!'
          RETURN   
       ENDIF
    ENDIF
    
  END SUBROUTINE zcartesian2spherical3D
  
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
          IF(rpt.EQ.0.0_dp) THEN
             theta_inp(ii) = pi / 2.0_dp
          ELSE
             theta_inp(ii) = ACOS( z_ax(iz) / rpt )
          ENDIF
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

END MODULE coords
