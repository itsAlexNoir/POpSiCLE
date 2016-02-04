MODULE cubboundcart
  
  USE constants
  USE tools
  USE gaussleg
  USE cubic_interp
  USE scattboundcyl
  USE scattboundcart
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           initialize_tricubic_cartesian_boundary3D_serial
#if _COM_MPI
  PUBLIC           initialize_tricubic_cartesian_boundary3D_parallel
#endif
  PUBLIC           dget_tricubic_cartesian_boundary3D  
  PUBLIC           zget_tricubic_cartesian_boundary3D
  PUBLIC           delete_tricubic_cartesian_boundary3D
  
  ! Module variables
  
  INTEGER, ALLOCATABLE, PUBLIC     :: cellindex3D(: ,: ,:, :)
  
  INTEGER                          :: numrpts
  INTEGER                          :: numthetapts
  INTEGER                          :: numphipts
  INTEGER                          :: numsurfaceprocs
  
CONTAINS
  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  INITIALIZE CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  
  SUBROUTINE initialize_tricubic_cartesian_boundary3D_serial(x_ax, y_ax, z_ax, &
       dims, Rb, fdpts, deltar, lmax, maxrpts, maxthetapts, maxphipts )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: x_ax(:)
    REAL(dp), INTENT(IN)      :: y_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rb
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(OUT)      :: maxrpts
    INTEGER, INTENT(OUT)      :: maxthetapts
    INTEGER, INTENT(OUT)      :: maxphipts
    
    REAL(dp)                  :: minx, maxx
    REAL(dp)                  :: miny, maxy
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    REAL(dp)                  :: minphi, maxphi
    INTEGER                   :: halolims(3,2)
    REAL(dp)                  :: Rb_start
    REAL(dp)                  :: deltaphi
    REAL(dp)                  :: rpt, thetapt, phipt
    REAL(dp)                  :: xpt, ypt, zpt
    INTEGER                   :: ix, iy, iz
    INTEGER                   :: ir, itheta, iphi
    INTEGER                   :: left, mflag
    
    !--------------------------------------------------!
    
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
    
    IF ( (Rb.LT.minr).AND.(Rb.GT.maxr) ) THEN
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local3D(1,1,1))
       ALLOCATE(cellindex3D(1,1,1,1))
       
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
             mintheta = MIN(mintheta,thetapt)
             maxtheta = MAX(maxtheta,thetapt)
             minphi = MIN(minphi,phipt)
             maxphi = MAX(maxphi,phipt)
          ENDDO
       ENDDO
    ENDDO
    
    numrpts = 2 * fdpts + 1
    numthetapts = 2 * lmax + 1
    numphipts = 2 * lmax + 1
    deltaphi = twopi / REAL(numphipts, dp)
    
    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(theta_boundary(1:numthetapts))
    ALLOCATE(phi_boundary(1:numphipts))
    ALLOCATE(costheta_boundary(1:numthetapts))
    ALLOCATE(theta_weights(1:numthetapts))
    ALLOCATE(i_am_in_local3D(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(cellindex3D(1:numrpts,1:numthetapts,1:numphipts,3))
    
    Rb_start = (Rb - deltar * (fdpts + 1))
    DO ir = 1, numrpts
       rpts_boundary(ir) = Rb_start + REAL( ir * deltar, dp )
    ENDDO
    
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, &
         costheta_boundary, theta_weights)
    
    theta_boundary = ACOS(costheta_boundary)
    
    ! Check limits on theta
    IF (mintheta .GT. MINVAL(theta_boundary) .AND. &
         maxtheta .LT. MAXVAL(theta_boundary)) THEN
       WRITE(*,*) 'Pivots extent is bigger than &
            & theta coordinate limits.'
       STOP
    ENDIF
    
    
    DO iphi = 1, numphipts
       phi_boundary(iphi) = REAL( iphi,dp ) * deltaphi
    ENDDO
    
    ! Check limits on phi
    IF (minphi .GT. MINVAL(phi_boundary) .AND. &
         maxphi .LT. MAXVAL(phi_boundary)) THEN
       WRITE(*,*) 'Phi nodes extent is bigger than &
            & phi coordinate limits.'
       STOP
    ENDIF
    
    !------------!
    
    i_am_in_local3D = 0
    cellindex3D = 0
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
    
    ! Also, initialize the tricubic matrix
    CALL initialize_tricubic_matrix( )
    
    maxrpts = numrpts
    maxthetapts = numthetapts
    maxphipts = numphipts
    
    WRITE(*,*)
    WRITE(*,*) '*************************************'
    WRITE(*,*) '  Surface parameters (Cartesian).    '
    WRITE(*,*) '*************************************'
    WRITE(*,*)
    WRITE(*,'(A, F9.3)') ' Radius boundary:                       ',&
         Rb
    WRITE(*,'(A, F9.3)') ' Delta R :                              ',&
         deltar
    WRITE(*,'(A,I3)')    ' Finite difference rule used :          ',&
         fdpts
    WRITE(*,'(A,I3)')    ' Maximum angular momentum:              ',&
         lmax
    WRITE(*,'(A,I3)')    ' Total number of theta points:          ',&
         numthetapts
    WRITE(*,'(A,I3)')    ' Total number of phi points:            ',&
         numphipts
    WRITE(*,*)
    WRITE(*,*)
    
    
  END SUBROUTINE initialize_tricubic_cartesian_boundary3D_serial
  
  !*****************************************************************!
  !*****************************************************************!
#if _COM_MPI
  
  SUBROUTINE initialize_tricubic_cartesian_boundary3D_parallel(x_ax, y_ax, z_ax, &
       dims, Rb, fdpts, deltar, lmax, mpi_rank, mpi_size, comm, &
       maxrpts, maxthetapts, maxphipts, surfacerank, maxsurfaceprocs, newcomm, &
       maxthetaptsperproc, maxphiptsperproc)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: x_ax(:)
    REAL(dp), INTENT(IN)      :: y_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rb
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(IN)       :: mpi_rank, mpi_size
    INTEGER, INTENT(IN)       :: comm
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts, maxphipts
    INTEGER, INTENT(OUT)      :: surfacerank
    INTEGER, INTENT(OUT)      :: maxsurfaceprocs
    INTEGER, INTENT(OUT)      :: newcomm
    INTEGER, INTENT(OUT)      :: maxthetaptsperproc, maxphiptsperproc
    
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
    REAL(dp)                  :: Rb_start    
    REAL(dp)                  :: rpt, thetapt, phipt
    REAL(dp)                  :: xpt, ypt, zpt
    INTEGER                   :: ix, iy, iz, inum
    INTEGER                   :: ir, itheta, iphi, ii
    INTEGER                   :: left, mflag
    
    !--------------------------------------------------!

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
    
    IF ( (Rb.LT.minr).AND.(Rb.GT.maxr) ) THEN
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(phi_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local3D(1,1,1))
       ALLOCATE(i_am_in_global3D(1,1,1))
       ALLOCATE(cellindex3D(1,1,1,1))
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
             phipt = pi - ATAN2( y_ax(iy) , -x_ax(ix) )
             
             ! See if this points lies within the radius of
             ! influence of the boundary
             
             ! Check the extend of the grid
             minr = MIN(minr,rpt)
             maxr = MAX(maxr,rpt)
             
             mintheta = MIN(mintheta,thetapt)
             maxtheta = MAX(maxtheta,thetapt)
             
             minphi = MIN(minphi,phipt)
             maxphi = MAX(maxphi,phipt)
             
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
    
    IF(ALL(x_ax.GT.0.0)) THEN
       minphi = pi - ATAN2(MINVAL(y_ax),-MINVAL(x_ax))
       maxphi = pi - ATAN2(MAXVAL(y_ax),-MINVAL(x_ax))
    ENDIF
    
    numrpts = 2 * fdpts + 1
    numthetapts = 2 * lmax + 1
    numphipts = 2 * lmax + 1
    deltaphi = (max_angles_global(2) - min_angles_global(2)) / REAL(numphipts, dp)
    
    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(theta_boundary(1:numthetapts))
    ALLOCATE(phi_boundary(1:numphipts))
    ALLOCATE(costheta_boundary(1:numthetapts))
    ALLOCATE(theta_weights(1:numthetapts))
    ALLOCATE(i_am_in_local3D(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(i_am_in_global3D(1:numrpts,1:numthetapts,1:numphipts))
    ALLOCATE(cellindex3D(1:numrpts,1:numthetapts,1:numphipts,3))
    i_am_in_local3D = 0
    i_am_in_global3D = 0
    
    ! Create desired axis for interpolation
    ! Radial axis
    Rb_start = (Rb - deltar * (fdpts + 1))
    
    DO ir = 1, numrpts
       rpts_boundary(ir) = Rb_start + REAL( ir * deltar, dp )
    ENDDO
    
    ! Now theta axis
    CALL get_gauss_stuff(-1.0_dp, 1.0_dp, &
         costheta_boundary, theta_weights)
    
    theta_boundary = ACOS(costheta_boundary)

    ! PHI AXIS
    DO iphi = 1, numphipts
       phi_boundary(iphi) = REAL( iphi,dp ) * deltaphi
    ENDDO
    
    !------------!
    
    i_am_in_local3D = 0
    cellindex3D = 0
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
                
                ! Look for cell indexes
                CALL interv(x_ax,SIZE(x_ax),ypt,left,mflag)
                cellindex3D(ir,itheta,iphi,1) = left
                CALL interv(y_ax,SIZE(y_ax),ypt,left,mflag)
                cellindex3D(ir,itheta,iphi,2) = left
                CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
                cellindex3D(ir,itheta,iphi,3) = left
                
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    
    i_am_in_global3D = 0
    ! Get I am in global array if needed.
    CALL MPI_ALLREDUCE(i_am_in_local3D, i_am_in_global3D, &
         numrpts*numthetapts*numphipts, MPI_INTEGER, MPI_SUM, comm, ierror )
    
    ! Check if there is a overlap of points between processors
    IF(ANY(i_am_in_global3D.GT.1)) THEN
       WRITE(*,*) 'There is an overlap of points beetween processors'
       IF(mpi_rank.EQ.0) THEN
          OPEN(UNIT=101,FORM='formatted',FILE='i_am_in.dat')
          DO iphi = 1, numphipts
             DO itheta = 1, numthetapts
                DO ir = 1, numrpts
                   WRITE(101,*) rpts_boundary(ir), &
                        theta_boundary(itheta), phi_boundary(iphi), &
                        i_am_in_global3D(ir,itheta,iphi)
                ENDDO
             ENDDO
          ENDDO
          CLOSE(101)
       ENDIF
       STOP
    ENDIF
    
    ! Also, initialize the bicubic matrix
    CALL initialize_tricubic_matrix( )
    
    ! Create communicator and surface members
    !-----------------------------------------!
    ALLOCATE(i_am_surface_local(0:mpi_size-1))
    ALLOCATE(i_am_surface(0:mpi_size-1))
    i_am_surface_local = 0
    i_am_surface = 0
    
    IF(SUM(i_am_in_local3D).NE.0) &
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
    CALL MPI_COMM_CREATE(comm, surfacegroup, newcomm, ierror)
    
    !Assign ranks for those who are on the surface
    IF(i_am_surface_local(mpi_rank).EQ.1) &
         CALL MPI_COMM_RANK( newcomm, surfacerank, ierror)
    
    ! Assign points to return
    maxrpts = numrpts
    maxthetapts = numthetapts
    maxphipts = numphipts
    maxsurfaceprocs = numsurfaceprocs
    
    ! Set number of points per proc
    ! Take into account if the number of points is not
    ! a multiplier of the number of processors
    maxthetaptsperproc = INT(numthetapts / numsurfaceprocs)
    
    IF(MOD(numthetapts,numsurfaceprocs).NE.0) THEN
       IF((surfacerank-1).EQ.numsurfaceprocs) THEN
          maxthetaptsperproc = numthetapts / numsurfaceprocs + &
               MOD(numthetapts,numsurfaceprocs)
       ENDIF
    ENDIF

    maxphiptsperproc = INT(numphipts / numsurfaceprocs)
    
    IF(MOD(numphipts,numsurfaceprocs).NE.0) THEN
       IF((surfacerank-1).EQ.numsurfaceprocs) THEN
          maxphiptsperproc = numthetapts / numsurfaceprocs + &
               MOD(numphipts,numsurfaceprocs)
       ENDIF
    ENDIF
    
    IF(mpi_rank.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '*************************************'
       WRITE(*,*) '  Surface parameters (Cartesian).    '
       WRITE(*,*) '*************************************'
       WRITE(*,*)
       WRITE(*,'(A, F9.3)') ' Radius boundary:                              ',&
            Rb
       WRITE(*,'(A,I3)')    ' Maximum angular momentum:                     ',&
            lmax
       WRITE(*,'(A, F9.3)') ' Delta R :                                     ',&
            deltar
       WRITE(*,'(A,I3)')    ' Finite difference rule used :                 ',&
            fdpts
       WRITE(*,'(A,I3)')    ' Number of processors on the surface:          ',&
            numsurfaceprocs
       WRITE(*,'(A,I3)')    ' Total number of theta points:                 ',&
            numthetapts
       WRITE(*,'(A,I3)')    ' Number of theta points per processors:        ',&
            maxthetaptsperproc
       IF(MOD(numthetapts,numsurfaceprocs).NE.0) THEN
          WRITE(*,'(A,I3)') ' Number of theta points on the last processor: ',&
               maxthetaptsperproc + MOD(numthetapts,numsurfaceprocs)
       ENDIF
       WRITE(*,'(A,I3)')    ' Total number of phi points:                   ',&
            numphipts
       WRITE(*,'(A,I3)')    ' Number of phi points per processors:          ',&
            maxphiptsperproc
       IF(MOD(numphipts,numsurfaceprocs).NE.0) THEN
          WRITE(*,'(A,I3)') ' Number of phi points on the last processor:   ',&
               maxphiptsperproc + MOD(numphipts,numsurfaceprocs)
       ENDIF
       WRITE(*,'(A,F9.3)')  ' Delta phi :                                   ',&
            deltaphi
       WRITE(*,*)
       WRITE(*,*)
    ENDIF
    
    
  END SUBROUTINE initialize_tricubic_cartesian_boundary3D_parallel
#endif
  !-----------------------------------------------------------------!
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !         GET BICUBIC CARTESIAN BOUNDARY 2D
  !----------------------------------------------------!
  
  SUBROUTINE dget_tricubic_cartesian_boundary3D(x_ax, y_ax, z_ax, dims, psi_cart, &
       fdrule, psi3D_sph, psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)             :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)              :: dims(:)
    REAL(dp), INTENT(IN)             :: psi_cart(:, :, :)
    INTEGER, INTENT(IN)              :: fdrule
    REAL(dp), INTENT(OUT)            :: psi3D_sph(:, :, :)
    REAL(dp), INTENT(OUT)            :: psi3D_sph_dx(:, :, :)
    REAL(dp), INTENT(OUT)            :: psi3D_sph_dy(:, :, :)
    REAL(dp), INTENT(OUT)            :: psi3D_sph_dz(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL    :: rank

    REAL(dp)                         :: psi_in(8)
    REAL(dp)                         :: psi_in_dx(8)
    REAL(dp)                         :: psi_in_dy(8) 
    REAL(dp)                         :: psi_in_dz(8)
    REAL(dp)                         :: psi_in_dxdy(8)
    REAL(dp)                         :: psi_in_dxdz(8)
    REAL(dp)                         :: psi_in_dydz(8)
    REAL(dp)                         :: psi_in_dxdydz(8)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_x(:, :)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_y(:, :)  
    REAL(dp), ALLOCATABLE            :: psi_in_vol_z(:, :)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_xy(:, :, :) 
    REAL(dp), ALLOCATABLE            :: psi_in_vol_xz(:, :, :) 
    REAL(dp), ALLOCATABLE            :: psi_in_vol_yz(:, :, :)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_xyz(:, :, :, :) 
    REAL(dp), ALLOCATABLE            :: xfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: yfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: zfdcoeffs(:, :)  
    REAL(dp)                         :: deriv, dx, dy, dz
    INTEGER                          :: left_x, left_y, left_z
    REAL(dp)                         :: xpt, ypt, zpt
    INTEGER                          :: rulepts, centralpt
    INTEGER                          :: lower_x, lower_y, lower_z
    INTEGER                          :: upper_x, upper_y, upper_z
    INTEGER                          :: ix, iy, iz, ipt
    INTEGER                          :: ir, itheta, iphi
    
    !---------------------------------------------------------------!
    
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
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                psi_in       = 0.0_dp
                psi_in_dx    = 0.0_dp
                psi_in_dy    = 0.0_dp
                psi_in_dz    = 0.0_dp
                psi_in_dxdy  = 0.0_dp
                psi_in_dxdz  = 0.0_dp
                psi_in_dydz  = 0.0_dp
                psi_in_dxdydz  = 0.0_dp
                
                xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                
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
                   upper_y = left_x
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
                
                CALL tricubic_interpolation(psi_in,psi_in_dx,psi_in_dy, psi_in_dz, &
                     psi_in_dxdy,psi_in_dxdz,psi_in_dydz,psi_in_dxdydz,&
                     x_ax(left_x),x_ax(left_x+1),&
                     y_ax(left_y),y_ax(left_y+1),&
                     z_ax(left_z),z_ax(left_z+1),&
                     xpt,ypt,zpt,&
                     psi3D_sph(ir,itheta,iphi),psi3D_sph_dx(ir,itheta,iphi),&
                     psi3D_sph_dy(ir,itheta,iphi),psi3D_sph_dz(ir,itheta,iphi))
                
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
    
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    
  END SUBROUTINE dget_tricubic_cartesian_boundary3D
  
  !******************************************************************!
  
  SUBROUTINE zget_tricubic_cartesian_boundary3D(x_ax, y_ax, z_ax, dims, psi_cart, &
       fdrule, psi3D_sph, psi3D_sph_dx, psi3D_sph_dy, psi3D_sph_dz, rank)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)             :: x_ax(:), y_ax(:), z_ax(:)
    INTEGER, INTENT(IN)              :: dims(:)
    COMPLEX(dp), INTENT(IN)          :: psi_cart(:, :, :)
    INTEGER, INTENT(IN)              :: fdrule
    COMPLEX(dp), INTENT(OUT)         :: psi3D_sph(:, :, :)
    COMPLEX(dp), INTENT(OUT)         :: psi3D_sph_dx(:, :, :)
    COMPLEX(dp), INTENT(OUT)         :: psi3D_sph_dy(:, :, :)
    COMPLEX(dp), INTENT(OUT)         :: psi3D_sph_dz(:, :, :)
    INTEGER, INTENT(IN), OPTIONAL    :: rank
    
    COMPLEX(dp)                      :: psi_in(8)
    COMPLEX(dp)                      :: psi_in_dx(8)
    COMPLEX(dp)                      :: psi_in_dy(8) 
    COMPLEX(dp)                      :: psi_in_dz(8)
    COMPLEX(dp)                      :: psi_in_dxdy(8)
    COMPLEX(dp)                      :: psi_in_dxdz(8)
    COMPLEX(dp)                      :: psi_in_dydz(8)
    COMPLEX(dp)                      :: psi_in_dxdydz(8)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_x(:, :)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_y(:, :)  
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_z(:, :)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_xy(:, :, :) 
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_xz(:, :, :) 
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_yz(:, :, :) 
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_xyz(:, :, :, :)
    REAL(dp)                         :: realpsi_in(8), imagpsi_in(8)
    REAL(dp)                         :: realpsi_in_dx(8)
    REAL(dp)                         :: imagpsi_in_dx(8)
    REAL(dp)                         :: realpsi_in_dy(8)
    REAL(dp)                         :: imagpsi_in_dy(8)
    REAL(dp)                         :: realpsi_in_dz(8)
    REAL(dp)                         :: imagpsi_in_dz(8)
    REAL(dp)                         :: realpsi_in_dxdy(8)
    REAL(dp)                         :: imagpsi_in_dxdy(8)
    REAL(dp)                         :: realpsi_in_dxdz(8)
    REAL(dp)                         :: imagpsi_in_dxdz(8)
    REAL(dp)                         :: realpsi_in_dydz(8)
    REAL(dp)                         :: imagpsi_in_dydz(8)
    REAL(dp)                         :: realpsi_in_dxdydz(8)
    REAL(dp)                         :: imagpsi_in_dxdydz(8)
    REAL(dp)                         :: realpsi3D_sph, imagpsi3D_sph
    REAL(dp)                         :: realpsi3D_sph_dx, imagpsi3D_sph_dx
    REAL(dp)                         :: realpsi3D_sph_dy, imagpsi3D_sph_dy
    REAL(dp)                         :: realpsi3D_sph_dz, imagpsi3D_sph_dz
    REAL(dp), ALLOCATABLE            :: xfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: yfdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: zfdcoeffs(:, :)  
    COMPLEX(dp)                      :: deriv
    REAL(dp)                         :: dx, dy, dz
    INTEGER                          :: left_x, left_y, left_z
    REAL(dp)                         :: xpt, ypt, zpt
    INTEGER                          :: rulepts, centralpt
    INTEGER                          :: lower_x, lower_y, lower_z
    INTEGER                          :: upper_x, upper_y, upper_z
    INTEGER                          :: ix, iy, iz, ipt
    INTEGER                          :: ir, itheta, iphi
    
    !---------------------------------------------------------------!
    
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
    
    DO iphi = 1, numphipts
       DO itheta = 1, numthetapts
          DO ir = 1, numrpts
             IF(i_am_in_local3D(ir,itheta,iphi).EQ.1) THEN
                psi_in       = ZERO
                psi_in_dx    = ZERO
                psi_in_dy    = ZERO
                psi_in_dz    = ZERO
                psi_in_dxdy  = ZERO
                psi_in_dxdz  = ZERO
                psi_in_dydz  = ZERO
                psi_in_dxdydz  = ZERO
                
                xpt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * COS(phi_boundary(iphi))
                ypt = rpts_boundary(ir) * SIN(theta_boundary(itheta)) * SIN(phi_boundary(iphi))
                zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
                
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
                   upper_y = left_x
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
                
                
                psi3D_sph(ir,itheta,iphi)    = CMPLX(realpsi3D_sph,imagpsi3D_sph)
                psi3D_sph_dx(ir,itheta,iphi) = CMPLX(realpsi3D_sph_dx,imagpsi3D_sph_dx)
                psi3D_sph_dy(ir,itheta,iphi) = CMPLX(realpsi3D_sph_dy,imagpsi3D_sph_dy)
                psi3D_sph_dz(ir,itheta,iphi) = CMPLX(realpsi3D_sph_dz,imagpsi3D_sph_dz)
                
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
    
    DEALLOCATE(xfdcoeffs,yfdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_x,psi_in_vol_y,psi_in_vol_z)
    DEALLOCATE(psi_in_vol_xy,psi_in_vol_xz,psi_in_vol_yz,psi_in_vol_xyz)
    
  END SUBROUTINE zget_tricubic_cartesian_boundary3D
  
  !****************************************************************!
  
  SUBROUTINE delete_tricubic_cartesian_boundary3D(x,y,z,w, &
       rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: x(:), y(:), z(:)
    REAL(dp), INTENT(IN)            :: w(:, :)
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary)
    DEALLOCATE(phi_boundary)
    DEALLOCATE(costheta_boundary)
    DEALLOCATE(theta_weights)
    DEALLOCATE(i_am_in_local3D)
    DEALLOCATE(cellindex3D)
    CALL delete_tricubic_matrix( )
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global3D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_tricubic_cartesian_boundary3D
  
  !----------------------------------------------------------------!
  !******************************************************************!
  !******************************************************************!
  
END MODULE cubboundcart
