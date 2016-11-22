MODULE scattboundcart
  
  USE constants_pop
  USE gaussleg
  USE scatt_interp
  USE tools
  USE scattboundcyl
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           initialize_scatt_cartesian_boundary2D
  PUBLIC           initialize_scatt_cartesian_boundary3D
  PUBLIC           get_scatt_cartesian_boundary2D
  PUBLIC           dget_scatt_cartesian_boundary3D
  PUBLIC           zget_scatt_cartesian_boundary3D
  
  PUBLIC           delete_scatt_cartesian_boundary2D
  PUBLIC           delete_scatt_cartesian_boundary3D
  
  
  ! Module variables
  
  REAL(dp), ALLOCATABLE, PUBLIC    :: phi_boundary(:)
  
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_local3D(:, :, :)
  INTEGER, ALLOCATABLE, PUBLIC     :: i_am_in_global3D(:, :, :)
  INTEGER, PUBLIC                  :: numphipts
  INTEGER, PUBLIC                  :: numphiptsperproc
  
  INTEGER                          :: numpts
  REAL(dp), ALLOCATABLE            :: xpts_scatt(:)
  REAL(dp), ALLOCATABLE            :: ypts_scatt(:)
  REAL(dp), ALLOCATABLE            :: zpts_scatt(:) 
  COMPLEX(dp), ALLOCATABLE         :: psi_scatt(:)
  INTEGER, ALLOCATABLE             :: index_x1(:)
  INTEGER, ALLOCATABLE             :: index_x2(:)
  INTEGER, ALLOCATABLE             :: index_x3(:)
  
CONTAINS
  
  !----------------------------------------------------!
  !  INITIALIZE CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE initialize_scatt_cartesian_boundary2D(x_ax, y_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, &
       maxpts)
    
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
    
    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: xpt, ypt
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
       ALLOCATE(xpts_scatt(1))
       ALLOCATE(ypts_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1)) 
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local2D(1,1))

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
       ALLOCATE(xpts_scatt(1))
       ALLOCATE(ypts_scatt(1))  
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local2D(1,1))
       
    ELSE
       ALLOCATE(xpts_scatt(1:numpts))
       ALLOCATE(ypts_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       
       numrpts = 2 * fdpts + 1
       numthetapts = lmax + 1
       !numthetapts = INT( (maxtheta - mintheta) / deltatheta )
       
       ALLOCATE(rpts_boundary(1:numrpts))
       ALLOCATE(theta_boundary(1:numthetapts))
       ALLOCATE(costheta_boundary(1:numthetapts))
       ALLOCATE(theta_weights(1:numthetapts))
       ALLOCATE(i_am_in_local2D(1:numrpts,1:numthetapts))
       
       Rs_start = (Rs - deltar * (fdpts + 1))
       DO ir = 1, numrpts
          rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
       ENDDO
       
       DO itheta = 1, numthetapts 
          theta_boundary(itheta) = REAL(itheta,dp) * twopi / REAL(numthetapts,dp)
       ENDDO
       
       !DO itheta = 1, numthetapts
       !   theta_boundary(itheta) = mintheta + &
       !REAL( itheta * deltatheta,dp )
       !ENDDO
       
       ! Check limits on theta
       IF (mintheta .GT. MINVAL(theta_boundary) .AND. &
            maxtheta .LT. MAXVAL(theta_boundary)) THEN
          WRITE(*,*) 'Pivots extent is bigger than &
               & theta coordinate limits.'
          STOP
       ENDIF

       i_am_in_local2D = 0
       ! Check if the new axis are within the cartesian domain.
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             
             xpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
             
             IF((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax))) THEN
                i_am_in_local2D(ir,itheta) = 1
                
             ENDIF
          ENDDO
       ENDDO
       
       inum = 0       
       DO iy = halolims(2,1), halolims(2,2)
          DO ix = halolims(1,1), halolims(1,2)
             
             rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 )
             
             IF (ABS(rpt-Rs) .LE. radtol ) THEN
                inum = inum + 1
                xpts_scatt(inum)  = x_ax(ix)
                ypts_scatt(inum) = y_ax(iy)
                index_x1(inum) = ix
                index_x2(inum) = iy         
             ENDIF
          ENDDO
       ENDDO
       
    ENDIF
    
    maxpts  = numpts
    
    WRITE(*,*)
    WRITE(*,*) '*************************************'
    WRITE(*,*) '  Surface parameters (Cartesian).    '
    WRITE(*,*) '*************************************'
    WRITE(*,*)
    WRITE(*,'(A, F9.3)') ' Radius boundary:                       ',&
         Rs
    WRITE(*,'(A,I3)')    ' Maximum angular momentum:              ',&
         lmax
    WRITE(*,'(A,F9.3)')  ' Radius tolerance :                     ',&
         radtol
    WRITE(*,'(A, F9.3)') ' Delta R :                              ',&
         deltar
    WRITE(*,'(A,I3)')    ' Finite difference rule used :          ',&
         fdpts
    WRITE(*,'(A,I3)')    ' Total number of theta points:          ',&
         numthetapts
    WRITE(*,*)
    WRITE(*,*)
    
    
  END SUBROUTINE initialize_scatt_cartesian_boundary2D
  
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !  INITIALIZE CARTESIAN BOUNDARY 3D
  !----------------------------------------------------!
  
  !-----------------------------------------------------------------!
  
  SUBROUTINE initialize_scatt_cartesian_boundary3D(x_ax, y_ax, z_ax, dims, &
       Rs, radtol, fdpts, deltar, lmax, mpi_rank, mpi_size, comm )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: x_ax(:)
    REAL(dp), INTENT(IN)            :: y_ax(:)
    REAL(dp), INTENT(IN)            :: z_ax(:)
    INTEGER, INTENT(IN)             :: dims(:)
    REAL(dp), INTENT(IN)            :: Rs
    REAL(dp), INTENT(IN)            :: radtol
    INTEGER, INTENT(IN)             :: fdpts
    REAL(dp), INTENT(IN)            :: deltar
    INTEGER, INTENT(IN)             :: lmax
    INTEGER, INTENT(IN), OPTIONAL   :: mpi_rank, mpi_size
    INTEGER, INTENT(IN), OPTIONAL   :: comm
    
    INTEGER                   :: simgroup
    INTEGER                   :: surfacegroup
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
       ALLOCATE(xpts_scatt(1))
       ALLOCATE(ypts_scatt(1)) 
       ALLOCATE(zpts_scatt(1)) 
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(index_x3(1))
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(phi_boundary(1))
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
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       ! Communicate min and max phi values
       max_angles_local = (/ maxtheta, maxphi /)
       min_angles_local = (/ mintheta, minphi /)
       
#if _COM_MPI   
       CALL MPI_ALLREDUCE(max_angles_local, max_angles_global, 2, &
            MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierror )
       CALL MPI_ALLREDUCE(min_angles_local, min_angles_global, 2, &
            MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierror )
#endif
       IF(ALL(x_ax.GT.0.0)) THEN
          minphi = pi - ATAN2(MINVAL(y_ax),-MINVAL(x_ax))
          maxphi = pi - ATAN2(MAXVAL(y_ax),-MINVAL(x_ax))
       ENDIF
    ENDIF
    
    numrpts = 2 * fdpts + 1
    numthetapts = lmax + 1
    numphipts = 2 * lmax + 1
    deltaphi = (max_angles_global(2) - min_angles_global(2)) / REAL(numphipts, dp)

    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(theta_boundary(1:numthetapts))
    ALLOCATE(costheta_boundary(1:numthetapts))
    ALLOCATE(theta_weights(1:numthetapts))
    ALLOCATE(phi_boundary(1:numphipts))
    ALLOCATE(i_am_in_local3D(1:numrpts,1:numthetapts,1:numphipts))
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) &
         ALLOCATE(i_am_in_global3D(1:numrpts,1:numthetapts,1:numphipts))
    
    ! RADIAL AXIS
    Rs_start = (Rs - deltar * (fdpts + 1))
    DO ir = 1, numrpts
       rpts_boundary(ir) = Rs_start + REAL( ir * deltar, dp )
    ENDDO
    
    ! THETA AXIS
    CALL make_gauss_legendre(-1.0_dp, 1.0_dp, costheta_boundary, &
         theta_weights)
    
    theta_boundary = ACOS(costheta_boundary)

#if _COM_MPI
    ! Check the extent
    IF((MINVAL(theta_boundary).LT.min_angles_global(1)) .OR. &
         MAXVAL(theta_boundary).GT.max_angles_global(1)) THEN
       IF(mpi_rank.EQ.0) &
            WRITE(*,*) 'Theta pivots are out of extent. '
       STOP
    ENDIF
#endif
    
    ! PHI AXIS
    DO iphi = 1, numphipts
       phi_boundary(iphi) = REAL( iphi,dp ) * deltaphi
    ENDDO
    
    !------------!
    
    i_am_in_local3D = 0
    ! Check if the new axis are within the cartesian domain.
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             
             xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
             ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
             zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             
             IF((xpt.GE.MINVAL(x_ax)) .AND. (xpt.LE.MAXVAL(x_ax)) .AND. &
                  (ypt.GE.MINVAL(y_ax)) .AND. (ypt.LE.MAXVAL(y_ax)) .AND. &
                  (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
                
                i_am_in_local3D(ir,itheta,iphi) = 1
                
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    
#if _COM_MPI
    !
    ! Set parameters for parallel calculation
    !
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN
       i_am_in_global3D = 0
       ! Get I am in global array if needed.
       CALL MPI_ALLREDUCE(i_am_in_local3D, i_am_in_global3D, numrpts*numthetapts*numphipts, &
            MPI_INTEGER, MPI_SUM, comm, ierror )
       
       ! Check if there is a overlap of points between processors
       IF(ANY(i_am_in_global3D.GT.1)) THEN
          WRITE(*,*) 'There is an overlap of points beetween processors'
          IF(mpi_rank.EQ.0) THEN
             OPEN(UNIT=101,FORM='formatted',FILE='i_am_in.dat')
             DO iphi = 1, numphipts
                DO itheta = 1, numthetapts
                   DO ir = 1, numrpts
                      WRITE(101,*) rpts_boundary(ir), theta_boundary(itheta), &
                           phi_boundary(iphi), i_am_in_global3D(ir,itheta,iphi)
                   ENDDO
                ENDDO
             ENDDO
             CLOSE(101)
          ENDIF
          STOP
       ENDIF
       
       ! Create communicator and surface members
       !-----------------------------------------!
       ALLOCATE(i_am_surface_local(0:mpi_size-1))
       ALLOCATE(i_am_surface(0:mpi_size-1))
       i_am_surface_local = 0
       i_am_surface = 0
       
       IF(numpts.NE.0) &
            i_am_surface_local(mpi_rank) = 1
       
       ! Communicate to all processors if they are
       ! at the surface
       CALL MPI_ALLREDUCE(i_am_surface_local, i_am_surface, mpi_size, &
            MPI_INTEGER, MPI_SUM, comm, ierror )
       
       numsurfaceprocs = SUM(i_am_surface)

       ALLOCATE(surface_members(0:numsurfaceprocs-1))
       ii = -1
       DO inum = 0, mpi_size-1
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
       CALL MPI_COMM_CREATE(comm, surfacegroup, surfacecomm, ierror)
       
       !Assign ranks for those who are on the surface
       IF(i_am_surface_local(mpi_rank).EQ.1) &
            CALL MPI_COMM_RANK( surfacecomm, surfacerank, ierror)
       
    ENDIF
#endif
    
    IF(numpts .EQ. 0) THEN
       ALLOCATE(xpts_scatt(1))
       ALLOCATE(ypts_scatt(1))
       ALLOCATE(zpts_scatt(1))
       ALLOCATE(psi_scatt(1))
       ALLOCATE(index_x1(1))
       ALLOCATE(index_x2(1))
       ALLOCATE(index_x3(1))
       
    ELSE
       ALLOCATE(xpts_scatt(1:numpts))
       ALLOCATE(ypts_scatt(1:numpts))
       ALLOCATE(zpts_scatt(1:numpts))
       ALLOCATE(psi_scatt(1:numpts))
       ALLOCATE(index_x1(1:numpts))
       ALLOCATE(index_x2(1:numpts))
       ALLOCATE(index_x3(1:numpts))
       
       ! Get points to participate in the interpolant
       inum = 0   
       DO iz = halolims(3,1), halolims(3,2)
          DO iy = halolims(2,1), halolims(2,2)
             DO ix = halolims(1,1), halolims(1,2)
                
                rpt = SQRT( x_ax(ix)**2 + y_ax(iy)**2 + z_ax(iz)**2 )
                
                IF (ABS(rpt-Rs) .LE. radtol ) THEN
                   inum = inum + 1
                   
                   xpts_scatt(inum) = x_ax(ix)
                   ypts_scatt(inum) = y_ax(iy)
                   zpts_scatt(inum) = z_ax(iz)
                   
                   index_x1(inum) = ix
                   index_x2(inum) = iy
                   index_x3(inum) = iz
                   
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
    ENDIF
    
    IF(PRESENT(mpi_size) .AND. (mpi_size.GT.1)) THEN    
       IF(mpi_rank.EQ.0) THEN
          WRITE(*,*)
          WRITE(*,*) '*************************************'
          WRITE(*,*) '  Surface parameters (Cartesian).    '
          WRITE(*,*) '*************************************'
          WRITE(*,*)
          WRITE(*,'(A, F9.3)') ' Radius boundary:                              ',&
               Rs
          WRITE(*,'(A,I3)')    ' Maximum angular momentum:                     ',&
               lmax
          WRITE(*,'(A,F9.3)')  ' Radius tolerance :                            ',&
               radtol
          WRITE(*,'(A, F9.3)') ' Delta R :                                     ',&
               deltar
          WRITE(*,'(A,I3)')    ' Finite difference rule used :                 ',&
               fdpts
          WRITE(*,'(A,I3)')    ' Number of processors on the surface:          ',&
               numsurfaceprocs
          WRITE(*,'(A,I3)')    ' Total number of theta points:                 ',&
               numthetapts
          WRITE(*,'(A,I3)')    ' Total number of phi points:                   ',&
               numphipts
          WRITE(*,'(A,F9.3)')  ' Delta phi :                                   ',&
               deltaphi
          WRITE(*,*)
          WRITE(*,*)
       ENDIF
    ELSE
       WRITE(*,*)
       WRITE(*,*) '*************************************'
       WRITE(*,*) '  Surface parameters (Cartesian).    '
       WRITE(*,*) '*************************************'
       WRITE(*,*)
       WRITE(*,'(A, F9.3)') ' Radius boundary:                              ',&
            Rs
       WRITE(*,'(A,I3)')    ' Maximum angular momentum:                     ',&
            lmax
       WRITE(*,'(A,F9.3)')  ' Radius tolerance :                            ',&
            radtol
       WRITE(*,'(A, F9.3)') ' Delta R :                                     ',&
            deltar
       WRITE(*,'(A,I3)')    ' Finite difference rule used :                 ',&
            fdpts
       WRITE(*,'(A,I3)')    ' Total number of theta points:                 ',&
            numthetapts
       WRITE(*,'(A,I3)')    ' Total number of phi points:                   ',&
            numphipts
       WRITE(*,'(A,F9.3)')  ' Delta phi :                                   ',&
            deltaphi
       WRITE(*,*)
       WRITE(*,*)
    ENDIF


       
  END SUBROUTINE initialize_scatt_cartesian_boundary3D
  
  !----------------------------------------------------!
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE get_scatt_cartesian_boundary2D(psi_cart,psi2D_sph, &
       psi2D_sph_dx, psi2D_sph_dy, method, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)          :: psi2D_sph_dy(:, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    INTEGER, INTENT(IN), OPTIONAL     :: rank

    REAL(dp)                          :: xpt, ypt
    INTEGER                           :: inum, ir, itheta
    
    psi_scatt = ZERO
    psi2D_sph = ZERO
    psi2D_sph_dx = ZERO
    psi2D_sph_dy = ZERO
    
    DO inum = 1, numpts
       psi_scatt(inum) = psi_cart(index_x1(inum),index_x2(inum))       
    ENDDO
    
    IF(PRESENT(rank)) THEN
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,psi_scatt,&
            TRIM(method),rank)
       
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             
             xpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
             
             CALL scatt_interpolate(numpts,xpt, ypt, &
                  xpts_scatt, ypts_scatt, psi_scatt, TRIM(method), &
                  psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), &
                  psi2D_sph_dy(ir,itheta), rank)
          ENDDO
       ENDDO
       
    ELSE
       
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,psi_scatt,&
            TRIM(method),rank)
       
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             
             xpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
             
             CALL scatt_interpolate(numpts,xpt, ypt, &
                  xpts_scatt, ypts_scatt, psi_scatt, TRIM(method), &
                  psi2D_sph(ir,itheta), psi2D_sph_dx(ir,itheta), &
                  psi2D_sph_dy(ir,itheta), rank)
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_scatt_interpolant(numpts, xpts_scatt, ypts_scatt, psi_scatt)
    
  END SUBROUTINE get_scatt_cartesian_boundary2D
  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  GET CARTESIAN BOUNDARY 3D
  !----------------------------------------------------!
  
  SUBROUTINE dget_scatt_cartesian_boundary3D(psi_cart,psi3D_sph, &
       psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, method, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)           :: psi_cart(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dx(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dy(:, :, :)
    REAL(dp), INTENT(OUT)          :: psi3D_sph_dz(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)   :: method
    INTEGER, INTENT(IN), OPTIONAL  :: rank

    REAL(dp)                       :: xpt, ypt, zpt
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
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,zpts_scatt, &
            REAL(psi_scatt),TRIM(method),rank)
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   
                   xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                   ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                   zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                   
                   CALL scatt_interpolate(numpts,xpt, ypt, zpt, &
                        xpts_scatt, ypts_scatt, zpts_scatt, &
                        REAL(psi_scatt), TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
    ELSE
       
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,zpts_scatt, &
            REAL(psi_scatt),TRIM(method),rank)
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   
                   xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                   ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                   zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                   
                   CALL scatt_interpolate(numpts,xpt, ypt, zpt, &
                        xpts_scatt, ypts_scatt, zpts_scatt, &
                        REAL(psi_scatt), TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_scatt_interpolant(numpts, xpts_scatt, ypts_scatt, zpts_scatt, &
         REAL(psi_scatt))    
    
  END SUBROUTINE dget_scatt_cartesian_boundary3D

  !-----------------------------------------!

    SUBROUTINE zget_scatt_cartesian_boundary3D(psi_cart,psi3D_sph, &
       psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, method, rank)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)           :: psi_cart(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dx(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dy(:, :, :)
    COMPLEX(dp), INTENT(OUT)          :: psi3D_sph_dz(:, :, :)
    CHARACTER(LEN=*), INTENT(IN)      :: method
    INTEGER, INTENT(IN), OPTIONAL     :: rank

    REAL(dp)                          :: xpt, ypt, zpt
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
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,zpts_scatt, &
            psi_scatt,TRIM(method),rank)
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                   
                   xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                   ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                   zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                   
                   CALL scatt_interpolate(numpts,xpt, ypt, zpt, &
                        xpts_scatt, ypts_scatt, zpts_scatt, &
                        psi_scatt, TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi),rank)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       
    ELSE
       
       CALL create_scatt_interpolant(numpts,xpts_scatt,ypts_scatt,zpts_scatt, &
            psi_scatt,TRIM(method))
       
       DO iphi = 1, numphipts
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN

                   xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                   ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                   zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                   
                   CALL scatt_interpolate(numpts,xpt, ypt, zpt, &
                        xpts_scatt, ypts_scatt, zpts_scatt, &
                        psi_scatt, TRIM(method), psi3D_sph(ir,itheta,iphi), &
                        psi3D_sph_dx(ir,itheta,iphi), psi3D_sph_dy(ir,itheta,iphi), &
                        psi3D_sph_dz(ir,itheta,iphi))
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    
    CALL destroy_scatt_interpolant(numpts, xpts_scatt, ypts_scatt, zpts_scatt, &
         psi_scatt)
    
  END SUBROUTINE zget_scatt_cartesian_boundary3D
  
  !----------------------------------------------------!
  !----------------------------------------------------!
  
  SUBROUTINE delete_scatt_cartesian_boundary2D(rpts,theta,rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rpts(:), theta(:)   
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(xpts_scatt,ypts_scatt)
    DEALLOCATE(psi_scatt) 
    DEALLOCATE(index_x1,index_x2)
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary)
    DEALLOCATE(i_am_in_local2D)
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global2D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_scatt_cartesian_boundary2D
  
  !----------------------------------------------------------------!
  
  SUBROUTINE delete_scatt_cartesian_boundary3D(rpts,theta,phi,rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rpts(:), theta(:), phi(:)
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(xpts_scatt,ypts_scatt)
    DEALLOCATE(zpts_scatt,psi_scatt) 
    DEALLOCATE(index_x1,index_x2,index_x3)
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary,costheta_boundary)
    DEALLOCATE(theta_weights)
    DEALLOCATE(phi_boundary)
    DEALLOCATE(i_am_in_local3D)
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global3D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_scatt_cartesian_boundary3D
  
  !---------------------------------------------!
  
END MODULE scattboundcart
