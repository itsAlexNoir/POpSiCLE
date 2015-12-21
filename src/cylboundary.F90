MODULE cylboundary
  
  USE constants
  USE gaussleg
  USE interp
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           initialize_cylindrical_boundary
  PUBLIC           get_cylindrical_boundary
  
  
  INTERFACE initialize_cylindrical_boundary
     MODULE PROCEDURE initialize_cylindrical_boundary2D_serial
#if _COM_MPI
     MODULE PROCEDURE initialize_cylindrical_boundary2D_parallel
#endif
  END INTERFACE initialize_cylindrical_boundary
  
  INTERFACE get_cylindrical_boundary
     MODULE PROCEDURE get_cylindrical_boundary2D
  END INTERFACE get_cylindrical_boundary
  
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: rpts_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: theta_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: costheta_boundary(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: theta_weights(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: pivots(:)
  REAL(dp), ALLOCATABLE, PUBLIC    :: theta_pivots(:)
  
  INTEGER, ALLOCATABLE, PUBLIC     :: pivot_number(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_surface_local(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_surface(:) 
  INTEGER, ALLOCATABLE, PUBLIC     :: surface_members(:)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_local2D(:, :)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_global2D(:, :)
  
  INTEGER                          :: numpts, numrpts
  INTEGER                          :: numthetapts
  INTEGER                          :: numpivots
  INTEGER                          :: numsurfaceprocs
  REAL(dp), ALLOCATABLE            :: rpts_scatt(:)
  REAL(dp), ALLOCATABLE            :: theta_scatt(:)
  COMPLEX(dp), ALLOCATABLE         :: psi_scatt(:)
  INTEGER, ALLOCATABLE             :: index_x1(:)
  INTEGER, ALLOCATABLE             :: index_x2(:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph_dx(:,:)
  COMPLEX(dp), ALLOCATABLE         :: psi2D_sph_dy(:,:)
  
CONTAINS
  
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
  
END MODULE cylboundary
