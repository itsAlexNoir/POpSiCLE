MODULE cartboundary
  
  USE constants
  USE gaussleg
  USE interp
  USE cylboundary
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           initialize_cartesian_boundary
  PUBLIC           get_cartesian_boundary
  PUBLIC           delete_cartesian_boundary
  
  
  INTERFACE initialize_cartesian_boundary
     MODULE PROCEDURE initialize_cartesian_boundary2D
     MODULE PROCEDURE initialize_cartesian_boundary3D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cartesian_boundary3D_parallel
#endif
  END INTERFACE initialize_cartesian_boundary
  
  INTERFACE get_cartesian_boundary
     MODULE PROCEDURE get_cartesian_boundary2D
     MODULE PROCEDURE dget_cartesian_boundary3D
     MODULE PROCEDURE zget_cartesian_boundary3D
  END INTERFACE get_cartesian_boundary
  
  INTERFACE delete_cartesian_boundary
     MODULE PROCEDURE delete_cartesian_boundary2D
     MODULE PROCEDURE delete_cartesian_boundary3D
  END INTERFACE delete_cartesian_boundary
  
  ! Module variables
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: phi_boundary(:)
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: phinodes(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: phinode_number(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_local3D(:, :, :)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_global3D(:, :, :)
  
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
       ALLOCATE(psi3D_sph_dz(1,1,1))
       
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
       ALLOCATE(index_x3(1))
       ALLOCATE(index_x3(1))
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
       ALLOCATE(psi3D_sph_dz(1,1,1))
       
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

    INTEGER                   :: simgroup, surfacegroup
    INTEGER                   :: ierror
    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    REAL(dp)                  :: minphi, maxphi
    REAL(dp)                  :: min_angles_local(2)
    REAL(dp)                  :: max_angles_local(2)
    REAL(dp)                  :: min_angles_global(2)
    REAL(dp)                  :: max_angles_global(2)
    REAL(dp)                  :: deltaphi
    INTEGER                   :: halolims(3,2)
    REAL(dp)                  :: rpt, thetapt, phipt
    REAL(dp)                  :: xpt, ypt, zpt
    REAL(dp)                  :: Rs_start
    INTEGER                   :: ix, iy, iz, inum
    INTEGER                   :: ir, itheta, iphi, ii
    
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
       ALLOCATE(theta_pivots(1))
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
       ALLOCATE(psi3D_sph_dz(1,1,1))
       ALLOCATE(i_am_in_local3D(1,1,1))
       ALLOCATE(i_am_in_global3D(1,1,1))
       ALLOCATE(i_am_surface_local(1))
       ALLOCATE(i_am_surface(1))
       
       RETURN
    ENDIF

    !-----------------!
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
             phipt = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
             
             ! See if this points lies within the radius of
             ! influence of the boundary
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                numpts = numpts + 1
                
                ! Check the extend of the grid
                minr = MIN(minr,rpt)
                maxr = MAX(maxr,rpt)
                
                mintheta = MIN(mintheta,thetapt)
                maxtheta = MAX(maxtheta,thetapt)
                
                minphi = MIN(minphi,phipt)
                maxphi = MAX(maxphi,phipt)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    ! Communicate min and max phi values
    max_angles_local = (/ maxtheta, maxphi /)
    min_angles_local = (/ mintheta, minphi /)
    
    CALL MPI_ALLREDUCE(max_angles_local, max_angles_global, 2, &
         MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierror )
    CALL MPI_ALLREDUCE(min_angles_local, min_angles_global, 2, &
         MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierror )

    numrpts = 2 * fdpts + 1
    numpivots = 2 * lmax + 1
    numphinodes = 2 * lmax + 1
    deltaphi = (max_angles_global(2) - min_angles_global(2)) / REAL(numphinodes, dp)

    ALLOCATE(pivot_number(1))
    ALLOCATE(phinode_number(1))
    ALLOCATE(theta_boundary(1),costheta_boundary(1))
    ALLOCATE(phi_boundary(1))
    
    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(pivots(1:numpivots))
    ALLOCATE(theta_pivots(1:numpivots))
    ALLOCATE(theta_weights(1:numpivots))
    ALLOCATE(phinodes(1:numphinodes))
    ALLOCATE(i_am_in_local3D(1:numrpts,1:numpivots,1:numphinodes))
    ALLOCATE(i_am_in_global3D(1:numrpts,1:numpivots,1:numphinodes))
    i_am_in_local3D = 0
    i_am_in_global3D = 0
    
    ! RADIAL AXIS
    Rs_start = (Rs - deltar * (fdpts + 1))
    DO ir = 1, numrpts
       rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
    ENDDO
    
    ! THETA AXIS
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, pivots, &
         theta_weights)
    
    theta_pivots = ACOS(pivots)
    ! Check the extent
    IF((MINVAL(theta_pivots).LT.min_angles_global(1)) .OR. &
         MAXVAL(theta_pivots).GT.max_angles_global(1)) THEN
       IF(rank.EQ.0) &
            WRITE(*,*) 'Theta pivots are out of extent. '
       STOP
    ENDIF
    
    ! PHI AXIS
    DO iphi = 1, numphinodes
       phinodes(iphi) = REAL( iphi,dp ) * deltaphi
    ENDDO
    !------------!
    
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, numphinodes
       DO itheta = 1, numpivots
          DO ir = 1, numrpts
             xpt = rpts_boundary(ir) * SIN(theta_pivots(itheta)) * COS(phinodes(iphi))
             ypt = rpts_boundary(ir) * SIN(theta_pivots(itheta)) * SIN(phinodes(iphi))
             zpt = rpts_boundary(ir) * COS(theta_pivots(itheta))
             IF ((rpts_boundary(ir).GE.minr) .AND. (rpts_boundary(ir).LE.maxr) .AND. &
                  (theta_pivots(itheta).GE.mintheta) .AND. &
                  (theta_pivots(itheta).LE.maxtheta) .AND. &
                  (phinodes(iphi).GE.minphi) .AND. (phinodes(iphi).LE.maxphi) .AND. &
                  (xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                i_am_in_local3D(ir,itheta,iphi) = 1
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
!!$    ! Get I am in global array if needed.
!!$    CALL MPI_ALLREDUCE(i_am_in_local, i_am_in_global, numrpts*numpivots*numphinodes, &
!!$         MPI_INTEGER, MPI_SUM, comm, ierror )
    
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

    numsurfaceprocs = SUM(i_am_surface)
    
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
       ALLOCATE(index_x3(1))
       ALLOCATE(psi3D_sph(1,1,1))
       ALLOCATE(psi3D_sph_dx(1,1,1))
       ALLOCATE(psi3D_sph_dy(1,1,1))
       ALLOCATE(psi3D_sph_dz(1,1,1))
     
    ELSE
       ALLOCATE(rpts_scatt(1:numpts))
       ALLOCATE(theta_scatt(1:numpts))
       ALLOCATE(phi_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       ALLOCATE(index_x3(1:numpts))
       
       ALLOCATE(psi3D_sph(1:numrpts,1:numpivots,1:numphinodes))
       ALLOCATE(psi3D_sph_dx(1:numrpts,1:numpivots,1:numphinodes))
       ALLOCATE(psi3D_sph_dy(1:numrpts,1:numpivots,1:numphinodes))
       ALLOCATE(psi3D_sph_dz(1:numrpts,1:numpivots,1:numphinodes))
       
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
!!$    maxthetapts = numthetapts
!!$    maxphipts = numphipts
    
    maxthetapts = numpivots
    maxphipts = numphinodes
    
    maxsurfaceprocs = numsurfaceprocs
    
  END SUBROUTINE initialize_cartesian_boundary3D_parallel
#endif  
  
  !----------------------------------------------------!
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE get_cartesian_boundary2D(psi_cart,psi2D_sph, &
       psi2D_sph_dx, psi2D_sph_dy, method, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dy(:, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    INTEGER, INTENT(IN), OPTIONAL     :: rank
    
    INTEGER                           :: inum, ir, itheta
    
    psi_scatt = ZERO
    psi2D_sph = ZERO
    psi2D_sph_dx = ZERO
    psi2D_sph_dy = ZERO
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum))       
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,psi_scatt,&
            TRIM(method),rank)
       
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             CALL interpolate(numpts,rpts_boundary(ir), theta_boundary(itheta), &
                  rpts_scatt, theta_scatt, psi_scatt, TRIM(method), &
                  psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), &
                  psi2D_sph_dy(ir,itheta), rank)
          ENDDO
       ENDDO
       
    ELSE
       
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,psi_scatt,&
            TRIM(method),rank)
       
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             CALL interpolate(numpts,rpts_boundary(ir), theta_boundary(itheta), &
                  rpts_scatt, theta_scatt, psi_scatt, TRIM(method), &
                  psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), &
                  psi2D_sph_dy(ir,itheta), rank)
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, psi_scatt)
    
  END SUBROUTINE get_cartesian_boundary2D
  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 3D
  !----------------------------------------------------!
  
  SUBROUTINE dget_cartesian_boundary3D(psi_cart,psi3D_sph, &
       psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, method, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: psi_cart(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dx(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dy(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dz(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: method
    INTEGER, INTENT(IN), OPTIONAL  :: rank
    
    INTEGER                        :: inum, ir, itheta, iphi

 
    psi_scatt = ZERO
    psi3D_sph = 0.0_dp
    psi3D_sph_dx = 0.0_dp
    psi3D_sph_dy = 0.0_dp
    psi3D_sph_dz = 0.0_dp
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum),&
            index_x3(inum))       
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,phi_scatt, &
            REAL(psi_scatt),TRIM(method),rank)
       
       DO iphi = 1, numphinodes
          DO itheta = 1, numpivots
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   CALL interpolate(numpts,rpts_boundary(ir), theta_pivots(itheta), &
                        phinodes(iphi),rpts_scatt, theta_scatt, phi_scatt, &
                        REAL(psi_scatt), TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ELSE
       
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,phi_scatt, &
            REAL(psi_scatt),TRIM(method),rank)
       
       DO iphi = 1, numphinodes
          DO itheta = 1, numpivots
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   CALL interpolate(numpts,rpts_boundary(ir), theta_pivots(itheta), &
                        phinodes(iphi),rpts_scatt, theta_scatt, phi_scatt, &
                        REAL(psi_scatt), TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, phi_scatt, &
         REAL(psi_scatt))    
    
  END SUBROUTINE dget_cartesian_boundary3D

  !-----------------------------------------!

    SUBROUTINE zget_cartesian_boundary3D(psi_cart,psi3D_sph, &
       psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, method, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dx(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dy(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dz(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    INTEGER, INTENT(IN), OPTIONAL     :: rank
    
    INTEGER                           :: inum, ir, itheta, iphi
    
    psi_scatt = ZERO
    psi3D_sph = ZERO
    psi3D_sph_dx = ZERO
    psi3D_sph_dy = ZERO
    psi3D_sph_dz = ZERO
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum),&
            index_x3(inum))       
    ENDDO

    IF(PRESENT(rank)) THEN
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,phi_scatt, &
            psi_scatt,TRIM(method),rank)
       
       DO iphi = 1, numphinodes
          DO itheta = 1, numpivots
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   CALL interpolate(numpts,rpts_boundary(ir), theta_pivots(itheta), &
                        phinodes(iphi),rpts_scatt, theta_scatt, phi_scatt, &
                        psi_scatt, TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO

    ELSE
       
       CALL create_interpolant(numpts,rpts_scatt,theta_scatt,phi_scatt, &
            psi_scatt,TRIM(method))
       
       DO iphi = 1, numphinodes
          DO itheta = 1, numpivots
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   CALL interpolate(numpts,rpts_boundary(ir), theta_pivots(itheta), &
                        phinodes(iphi),rpts_scatt, theta_scatt, phi_scatt, &
                        psi_scatt, TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_interpolant(numpts, rpts_scatt, theta_scatt, phi_scatt, &
         psi_scatt)
    
  END SUBROUTINE zget_cartesian_boundary3D
  
  !----------------------------------------------------!
  !----------------------------------------------------!
  
  SUBROUTINE delete_cartesian_boundary2D(rpts,theta,rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rpts(:), theta(:)   
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(rpts_scatt,theta_scatt)
    DEALLOCATE(psi_scatt) 
    DEALLOCATE(index_x1,index_x2)
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary)
    DEALLOCATE(psi2D_sph)
    DEALLOCATE(psi2D_sph_dx,psi2D_sph_dy)
    DEALLOCATE(i_am_in_local2D)
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global2D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_cartesian_boundary2D
  
  !----------------------------------------------------------------!
  
  SUBROUTINE delete_cartesian_boundary3D(rpts,theta,phi,rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rpts(:), theta(:), phi(:)
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(rpts_scatt,theta_scatt)
    DEALLOCATE(phi_scatt,psi_scatt) 
    DEALLOCATE(index_x1,index_x2,index_x3)
    DEALLOCATE(pivots)
    DEALLOCATE(theta_pivots)
    DEALLOCATE(phinodes)
    DEALLOCATE(pivot_number)
    DEALLOCATE(phinode_number)
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary,costheta_boundary)
    DEALLOCATE(theta_weights)
    DEALLOCATE(phi_boundary)
    DEALLOCATE(psi3D_sph)
    DEALLOCATE(psi3D_sph_dx,psi3D_sph_dy,psi3D_sph_dz)
    DEALLOCATE(i_am_in_local3D)
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global3D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_cartesian_boundary3D
  
  !---------------------------------------------!
  
END MODULE cartboundary
