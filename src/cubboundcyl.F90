MODULE cubboundcyl
  
  USE constants
  USE tools
  USE gaussleg
  USE cubic_interp
  USE scattboundcyl
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC           initialize_bicubic_cylindrical_boundary2D_serial
#if _COM_MPI
  PUBLIC           initialize_bicubic_cylindrical_boundary2D_parallel
#endif
  PUBLIC           dget_bicubic_cylindrical_boundary2D  
  PUBLIC           zget_bicubic_cylindrical_boundary2D
  PUBLIC           delete_bicubic_cylindrical_boundary2D
  
  ! Module variables
  
  INTEGER, ALLOCATABLE, PUBLIC     :: cellindex2D(: ,: ,:)
  
  INTEGER                          :: numrpts
  INTEGER                          :: numthetapts
  INTEGER                          :: numsurfaceprocs
  
CONTAINS
  
  !----------------------------------------------------!
  
  !----------------------------------------------------!
  !  INITIALIZE CYLINDRICAL BOUNDARY 2D
  !----------------------------------------------------!
  
  
  SUBROUTINE initialize_bicubic_cylindrical_boundary2D_serial(rho_ax, z_ax, &
       dims, Rb, fdpts, deltar, lmax, maxrpts, maxthetapts )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: rho_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rb
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(OUT)      :: maxrpts
    INTEGER, INTENT(OUT)      :: maxthetapts
    
    REAL(dp)                  :: minrho, maxrho
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: Rb_start
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: rhopt, zpt
    INTEGER                   :: irho, iz
    INTEGER                   :: ir, itheta
    INTEGER                   :: left, mflag
    
    !--------------------------------------------------!
    
    halolims(1,:) = (/ lbound(rho_ax), ubound(rho_ax) /)
    halolims(2,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minrho = MINVAL(ABS(rho_ax))
    maxrho = MAXVAL(ABS(rho_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minrho**2 + minz**2)
    maxr = SQRT(maxrho**2 + maxz**2)
    
    IF ( (Rb.LT.minr).AND.(Rb.GT.maxr) ) THEN
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local2D(1,1))
       ALLOCATE(cellindex2D(1,1,1))
       
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
          
          mintheta = MIN(mintheta,thetapt)
          maxtheta = MAX(maxtheta,thetapt)
          
       ENDDO
    ENDDO
    
    numrpts = 2 * fdpts + 1
    numthetapts = 2 * lmax + 1
    
    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(theta_boundary(1:numthetapts))
    ALLOCATE(costheta_boundary(1:numthetapts))
    ALLOCATE(theta_weights(1:numthetapts))
    ALLOCATE(i_am_in_local2D(1:numrpts,1:numthetapts))
    ALLOCATE(cellindex2D(1:numrpts,1:numthetapts,2))
    
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
      
    !------------!
    
    i_am_in_local2D = 0
    cellindex2D = 0
    ! Check if the new axis are within the cylindrical domain.
    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          
          rhopt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
          zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
          
          IF((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             
             i_am_in_local2D(ir,itheta) = 1
             
             ! Look for cell indexes
             CALL interv(rho_ax,SIZE(rho_ax),rhopt,left,mflag)
             cellindex2D(ir,itheta,1) = left
             CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
             cellindex2D(ir,itheta,2) = left
          ENDIF
          
       ENDDO
    ENDDO
    
    
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )
    
    maxrpts = numrpts
    maxthetapts = numthetapts
    
    WRITE(*,*)
    WRITE(*,*) '*************************************'
    WRITE(*,*) '  Surface parameters (Cylindrical).  '
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
    WRITE(*,*)
    WRITE(*,*)
    
    
  END SUBROUTINE initialize_bicubic_cylindrical_boundary2D_serial
  
  !*****************************************************************!
  !*****************************************************************!
#if _COM_MPI
  
  SUBROUTINE initialize_bicubic_cylindrical_boundary2D_parallel(rho_ax, z_ax, dims, &
       Rb, fdpts, deltar, lmax, mpi_rank, mpi_size, comm, &
       maxrpts, maxthetapts, surfacerank, maxsurfaceprocs, newcomm, &
       maxthetaptsperproc)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)      :: rho_ax(:)
    REAL(dp), INTENT(IN)      :: z_ax(:)
    INTEGER, INTENT(IN)       :: dims(:)
    REAL(dp), INTENT(IN)      :: Rb
    INTEGER, INTENT(IN)       :: fdpts
    REAL(dp), INTENT(IN)      :: deltar
    INTEGER, INTENT(IN)       :: lmax
    INTEGER, INTENT(IN)       :: mpi_rank, mpi_size
    INTEGER, INTENT(IN)       :: comm
    INTEGER, INTENT(OUT)      :: maxrpts, maxthetapts
    INTEGER, INTENT(OUT)      :: surfacerank
    INTEGER, INTENT(OUT)      :: maxsurfaceprocs
    INTEGER, INTENT(OUT)      :: newcomm
    INTEGER, INTENT(OUT)      :: maxthetaptsperproc
    
    INTEGER                   :: simgroup, surfacegroup
    INTEGER                   :: ierror
    REAL(dp)                  :: minrho, maxrho
    REAL(dp)                  :: minz, maxz
    REAL(dp)                  :: minr, maxr
    REAL(dp)                  :: mintheta, maxtheta
    INTEGER                   :: halolims(2,2)
    REAL(dp)                  :: Rb_start    
    REAL(dp)                  :: rpt, thetapt
    REAL(dp)                  :: rhopt, zpt
    INTEGER                   :: irho, iz, inum
    INTEGER                   :: ir, itheta, ii
    INTEGER                   :: left, mflag
 
    !--------------------------------------------------!
    
    halolims(1,:) = (/ lbound(rho_ax), ubound(rho_ax) /)
    halolims(2,:) = (/ lbound(z_ax), ubound(z_ax) /)
    
    ! Assign max and min values for the axes
    minrho = MINVAL(ABS(rho_ax))
    maxrho = MAXVAL(ABS(rho_ax))
    
    minz = MINVAL(ABS(z_ax))
    maxz = MAXVAL(ABS(z_ax))
    
    minr = SQRT(minrho**2 + minz**2)
    maxr = SQRT(maxrho**2 + maxz**2)
    
    IF ( (Rb.LT.minr).AND.(Rb.GT.maxr) ) THEN
       ALLOCATE(rpts_boundary(1))
       ALLOCATE(theta_boundary(1))
       ALLOCATE(costheta_boundary(1))
       ALLOCATE(theta_weights(1))
       ALLOCATE(i_am_in_local2D(1,1))
       ALLOCATE(i_am_in_global2D(1,1))
       ALLOCATE(cellindex2D(1,1,1))
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
          
          mintheta = MIN(mintheta,thetapt)
          maxtheta = MAX(maxtheta,thetapt)
       ENDDO
    ENDDO
    
    numrpts = 2 * fdpts + 1
    numthetapts = 2 * lmax + 1
    
    ALLOCATE(rpts_boundary(1:numrpts))
    ALLOCATE(theta_boundary(1:numthetapts))
    ALLOCATE(costheta_boundary(1:numthetapts))
    ALLOCATE(theta_weights(1:numthetapts))
    ALLOCATE(i_am_in_local2D(1:numrpts,1:numthetapts))
    ALLOCATE(i_am_in_global2D(1:numrpts,1:numthetapts))
    ALLOCATE(cellindex2D(1:numrpts,1:numthetapts,2))
    
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
    
    !------------!
    
    i_am_in_local2D = 0
    cellindex2D = 0
    ! Check if the new axis are within the cartesian domain.
    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          
          rhopt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
          zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
          
          IF((rhopt.GE.MINVAL(rho_ax)) .AND. (rhopt.LE.MAXVAL(rho_ax)) .AND. &
               (zpt.GE.MINVAL(z_ax)) .AND. (zpt.LE.MAXVAL(z_ax))) THEN
             
             i_am_in_local2D(ir,itheta) = 1

             ! Look for cell indexes
             CALL interv(rho_ax,SIZE(rho_ax),rhopt,left,mflag)
             cellindex2D(ir,itheta,1) = left
             CALL interv(z_ax,SIZE(z_ax),zpt,left,mflag)
             cellindex2D(ir,itheta,2) = left
             
          ENDIF
       ENDDO
    ENDDO
    
    i_am_in_global2D = 0
    ! Get I am in global array if needed.
    CALL MPI_ALLREDUCE(i_am_in_local2D, i_am_in_global2D, numrpts*numthetapts, &
         MPI_INTEGER, MPI_SUM, comm, ierror )
    
    ! Check if there is a overlap of points between processors
    IF(ANY(i_am_in_global2D.GT.1)) THEN
       WRITE(*,*) 'There is an overlap of points beetween processors'
       IF(mpi_rank.EQ.0) THEN
          OPEN(UNIT=101,FORM='formatted',FILE='i_am_in.dat')
          DO itheta = 1, numthetapts
             DO ir = 1, numrpts
                WRITE(101,*) rpts_boundary(ir), &
                     theta_boundary(itheta), i_am_in_global2D(ir,itheta)
             ENDDO
          ENDDO
          CLOSE(101)
       ENDIF
       STOP
    ENDIF
    
    ! Also, initialize the bicubic matrix
    CALL initialize_bicubic_matrix( )
    
    
    ! Create communicator and surface members
    !-----------------------------------------!
    ALLOCATE(i_am_surface_local(0:mpi_size-1))
    ALLOCATE(i_am_surface(0:mpi_size-1))
    i_am_surface_local = 0
    i_am_surface = 0

    
    IF(SUM(i_am_in_local2D).NE.0) &
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
    maxsurfaceprocs = numsurfaceprocs
    
    ! Set number of points per proc
    ! Take into account if the number of points is not
    ! a multiplier of the number of processors
    IF(numsurfaceprocs.GT.numthetapts) THEN
       maxthetaptsperproc = 1
    ELSE
       maxthetaptsperproc = INT(numthetapts / numsurfaceprocs)
       IF(MOD(numthetapts,numsurfaceprocs).NE.0) THEN
          IF(surfacerank.EQ.(numsurfaceprocs-1)) THEN
             maxthetaptsperproc = numthetapts / numsurfaceprocs + &
                  MOD(numthetapts,numsurfaceprocs)
          ENDIF
       ENDIF
    ENDIF
    
    IF(mpi_rank.EQ.0) THEN
       WRITE(*,*)
       WRITE(*,*) '*************************************'
       WRITE(*,*) '  Surface parameters (Cylindrical).  '
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
       WRITE(*,*)
       WRITE(*,*)
    ENDIF
    
    
  END SUBROUTINE initialize_bicubic_cylindrical_boundary2D_parallel
#endif
  !-----------------------------------------------------------------!
  !-----------------------------------------------------------------!
  
  !----------------------------------------------------!
  !         GET BICUBIC CYLINDRICAL BOUNDARY 2D
  !----------------------------------------------------!

  SUBROUTINE dget_bicubic_cylindrical_boundary2D(rho_ax, z_ax, dims, psi_cyl, &
       fdrule, psi2D_sph, psi2D_sph_dx, psi2D_sph_dy, rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)             :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)              :: dims(:)
    REAL(dp), INTENT(IN)             :: psi_cyl(:, :)
    INTEGER, INTENT(IN)              :: fdrule
    REAL(dp), INTENT(OUT)            :: psi2D_sph(:, :)
    REAL(dp), INTENT(OUT)            :: psi2D_sph_dx(:, :)
    REAL(dp), INTENT(OUT)            :: psi2D_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL    :: rank

    REAL(dp)                         :: psi_in(4), psi_in_drhodz(4)
    REAL(dp)                         :: psi_in_drho(4), psi_in_dz(4)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_rho(:, :)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_z(:, :)
    REAL(dp), ALLOCATABLE            :: psi_in_vol_rhoz(:, :, :) 
    REAL(dp), ALLOCATABLE            :: rhofdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: zfdcoeffs(:, :)  
    REAL(dp)                         :: deriv, drho, dz
    INTEGER                          :: left_rho, left_z
    REAL(dp)                         :: rhopt, zpt
    INTEGER                          :: rulepts, centralpt
    INTEGER                          :: lower_rho, lower_z
    INTEGER                          :: upper_rho, upper_z
    INTEGER                          :: irho, iz, ipt
    INTEGER                          :: ir, itheta
    
    !---------------------------------------------------------------!
    
    rulepts = 2 * fdrule
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          IF(i_am_in_local2D(ir,itheta).EQ.1) THEN
             psi_in      = 0.0_dp
             psi_in_drho = 0.0_dp
             psi_in_dz   = 0.0_dp
             psi_in_drhodz  = 0.0_dp
             
             rhopt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
             zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             
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
                  rulepts+1,2,rhofdcoeffs)
             
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
                  rulepts+1,2,zfdcoeffs)
          
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
             
             CALL bicubic_interpolation(psi_in,psi_in_drho, psi_in_dz, &
                  psi_in_drhodz, rho_ax(left_rho),rho_ax(left_rho+1),&
                  z_ax(left_z),z_ax(left_z+1),rhopt,zpt,&
                  psi2D_sph(ir,itheta),psi2D_sph_dx(ir,itheta),&
                  psi2D_sph_dy(ir,itheta))
             
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
    
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    
  END SUBROUTINE dget_bicubic_cylindrical_boundary2D
  
  !******************************************************************!
  
  SUBROUTINE zget_bicubic_cylindrical_boundary2D(rho_ax, z_ax, dims, psi_cyl, &
       fdrule, psi2D_sph, psi2D_sph_dx, psi2D_sph_dy, rank)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)             :: rho_ax(:), z_ax(:)
    INTEGER, INTENT(IN)              :: dims(:)
    COMPLEX(dp), INTENT(IN)          :: psi_cyl(:, :)
    INTEGER, INTENT(IN)              :: fdrule
    COMPLEX(dp), INTENT(OUT)         :: psi2D_sph(:, :)
    COMPLEX(dp), INTENT(OUT)         :: psi2D_sph_dx(:, :)
    COMPLEX(dp), INTENT(OUT)         :: psi2D_sph_dy(:, :)
    INTEGER, INTENT(IN), OPTIONAL    :: rank

    COMPLEX(dp)                      :: psi_in(4), psi_in_drhodz(4)
    COMPLEX(dp)                      :: psi_in_drho(4), psi_in_dz(4)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_rho(:, :)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_z(:, :)
    COMPLEX(dp), ALLOCATABLE         :: psi_in_vol_rhoz(:, :, :)
    REAL(dp)                         :: realpsi_in(4), imagpsi_in(4)
    REAL(dp)                         :: realpsi_in_drho(4)
    REAL(dp)                         :: imagpsi_in_drho(4)
    REAL(dp)                         :: realpsi_in_dz(4)
    REAL(dp)                         :: imagpsi_in_dz(4)
    REAL(dp)                         :: realpsi_in_drhodz(4)
    REAL(dp)                         :: imagpsi_in_drhodz(4)
    REAL(dp)                         :: realpsi2D_sph, imagpsi2D_sph
    REAL(dp)                         :: realpsi2D_sph_dx, imagpsi2D_sph_dx
    REAL(dp)                         :: realpsi2D_sph_dy, imagpsi2D_sph_dy
    REAL(dp), ALLOCATABLE            :: rhofdcoeffs(:, :)
    REAL(dp), ALLOCATABLE            :: zfdcoeffs(:, :)  
    COMPLEX(dp)                      :: deriv
    REAL(dp)                         :: drho, dz
    INTEGER                          :: left_rho, left_z
    REAL(dp)                         :: rhopt, zpt
    INTEGER                          :: rulepts, centralpt
    INTEGER                          :: lower_rho, lower_z
    INTEGER                          :: upper_rho, upper_z
    INTEGER                          :: irho, iz, ipt
    INTEGER                          :: ir, itheta
    
    !---------------------------------------------------------------!
    
    rulepts = 2 * fdrule
    ALLOCATE(rhofdcoeffs(-fdrule:fdrule,2))
    ALLOCATE(zfdcoeffs(-fdrule:fdrule,2))
    
    ALLOCATE(psi_in_vol_rho(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_z(-fdrule:fdrule,4))
    ALLOCATE(psi_in_vol_rhoz(-fdrule:fdrule,-fdrule:fdrule,4))
    
    DO itheta = 1, numthetapts
       DO ir = 1, numrpts
          IF(i_am_in_local2D(ir,itheta).EQ.1) THEN
             psi_in      = ZERO
             psi_in_drho = ZERO
             psi_in_dz   = ZERO
             psi_in_drhodz  = ZERO
             
             rhopt = rpts_boundary(ir) * SIN(theta_boundary(itheta))
             zpt = rpts_boundary(ir) * COS(theta_boundary(itheta))
             
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
                  rulepts+1,2,rhofdcoeffs)
             
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
                  rulepts+1,2,zfdcoeffs)

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
             
             psi2D_sph(ir,itheta)    = CMPLX(realpsi2D_sph,imagpsi2D_sph)
             psi2D_sph_dx(ir,itheta) = CMPLX(realpsi2D_sph_dx,imagpsi2D_sph_dx)
             psi2D_sph_dy(ir,itheta) = CMPLX(realpsi2D_sph_dy,imagpsi2D_sph_dy)
             
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
    
    DEALLOCATE(rhofdcoeffs,zfdcoeffs)
    DEALLOCATE(psi_in_vol_rho,psi_in_vol_z,psi_in_vol_rhoz)
    
  END SUBROUTINE zget_bicubic_cylindrical_boundary2D
  
  !****************************************************************!
  
  SUBROUTINE delete_bicubic_cylindrical_boundary2D(rpts,theta,w, &
       rank)
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)            :: rpts(:), theta(:)   
    REAL(dp), INTENT(IN)             :: w(: ,:)
    INTEGER, INTENT(IN), OPTIONAL   :: rank
    
    DEALLOCATE(rpts_boundary)
    DEALLOCATE(theta_boundary)
    DEALLOCATE(costheta_boundary)
    DEALLOCATE(theta_weights)
    DEALLOCATE(i_am_in_local2D)
    DEALLOCATE(cellindex2D)
    CALL delete_bicubic_matrix( )
    
    ! Delete arrays used for communications
    IF(PRESENT(rank)) THEN
       DEALLOCATE(i_am_in_global2D)
       DEALLOCATE(i_am_surface_local,i_am_surface)
       DEALLOCATE(surface_members)
    ENDIF
    
  END SUBROUTINE delete_bicubic_cylindrical_boundary2D
  
  !----------------------------------------------------------------!
  !******************************************************************!
  !******************************************************************!
  
END MODULE cubboundcyl
