MODULE fourier
  
#if _COM_MPI
  
  USE MPI
  
#endif
  
  USE omp_lib
  USE constants
  USE tools
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC                         :: FourierTransform
  PUBLIC                         :: FFT
  PUBLIC                         :: HankelTransform
  PUBLIC                         :: create_mask
  PUBLIC                         :: create_correlated_mask
  
  INTERFACE create_mask
     MODULE PROCEDURE create_mask4D
     MODULE PROCEDURE create_mask3D
     MODULE PROCEDURE create_mask2D
  END INTERFACE create_mask

  INTERFACE create_correlated_mask
     MODULE PROCEDURE create_correlated_mask3D1D
  END INTERFACE create_correlated_mask
  
  INTERFACE FourierTransform
     MODULE PROCEDURE FourierTransform_serial
#if _COM_MPI
     MODULE PROCEDURE FourierTransform_parallel
#endif
  END INTERFACE FourierTransform
  
CONTAINS
  
  SUBROUTINE FourierTransform_serial( psi, psik, rank, dims_local )
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: psi(:)
    COMPLEX(dp), INTENT(OUT)      :: psik(:)
    INTEGER, INTENT(IN)           :: rank
    INTEGER, INTENT(IN)           :: dims_local(:)
    
    !***************************************************!
    
    !---------------------------!
    ! Copy the original array   !
    !---------------------------!
    psik = psi
    
    !-------------------------------------------------!
    ! Transform the wavefunction into momentum space. !
    !-------------------------------------------------!
    
    WRITE(*,*) 'Transforming to momentum space...'
    
    ! Compute forward transform
    CALL FFT(psik,rank,dims_local)
    
    
  END SUBROUTINE FourierTransform_serial
  
  !----------------------------------------------------!
  
#if _COM_MPI
  
  SUBROUTINE FourierTransform_parallel(psi, psik, rank, dims_local, &
       dims_global, comm_gl, grid_rank, grid_addresses )
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: psi(:)
    COMPLEX(dp), INTENT(OUT)      :: psik(:)
    INTEGER, INTENT(IN)           :: rank
    INTEGER, INTENT(IN)           :: dims_local(:)
    INTEGER, INTENT(IN)           :: dims_global(:)
    INTEGER, INTENT(IN)           :: comm_gl
    INTEGER, INTENT(IN)           :: grid_rank(:)
    INTEGER, INTENT(IN)           :: grid_addresses(:)
    
    !--------------------------------------------------!

    COMPLEX(dp), ALLOCATABLE      :: psi_sd(:)
    COMPLEX(dp), ALLOCATABLE      :: psi_rc(:)
    INTEGER, ALLOCATABLE          :: fftcomm(:)
    INTEGER                       :: group_gl, group_lc
    INTEGER                       :: group_size
    INTEGER, ALLOCATABLE          :: members(:)
    INTEGER                       :: dims_phony(4)
    INTEGER                       :: Nx, Ny, Nz, Nr
    INTEGER                       :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                       :: ipx, ipy, ipz, ipr
    INTEGER                       :: numprocx, numprocy
    INTEGER                       :: numprocz, numprocr
    INTEGER                       :: ix, iy, iz, ir
    INTEGER                       :: igl, shift
    INTEGER                       :: iprocx, iprocy
    INTEGER                       :: iprocz, iprocr
    INTEGER                       :: ipro
    INTEGER                       :: ierror
    
    !***************************************************!
    
    ! Copy the original array
    psik = psi
    
    
    ! Assign dimensions
    CALL get_simulation_parameters( rank, &
         dims_local, dims_global, grid_rank, &
         Nx, Ny, Nz, Nr, Nxgl, Nygl, Nzgl, Nrgl, &
         ipx, ipy, ipz, ipr, &
         numprocx, numprocy, numprocz, numprocr )
        
    ! Create global group (for all the processors)
    CALL MPI_COMM_GROUP(comm_gl,group_gl,ierror)
    
    dims_phony = 1
    
    !********************
    ! FFT in x direction
    !********************
    
    ! In case there is only one processor in this dimension
    IF (numprocx .EQ. 1) THEN
       ALLOCATE(psi_sd(1:Nxgl))
       ALLOCATE(psi_rc(1:Nxgl))
       dims_phony(1) = Nxgl
       
       DO ir = 1, Nr
          DO iz = 1, Nz
             DO iy = 1, Ny
                
                ! Set arrays to zero at the begining
                psi_sd = ZERO
                psi_rc = ZERO
                
                DO ix =1, Nx
                   psi_rc(ix) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                ENDDO
                
                CALL FFT(psi_rc, 1, dims_phony)
                
                ! fftshift the array
                psi_sd = ZERO
                shift = Nxgl - (Nxgl+1)/2
                DO ix = 1, shift
                   psi_sd(ix) = psi_rc(ix+shift)
                ENDDO
                
                DO ix = shift+1, Nxgl
                   psi_sd(ix) = psi_rc(ix-shift)
                ENDDO
                
                ! Copy the result back to the local psik
                DO ix = 1, Nx
                   psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_sd(ix)
                ENDDO
                
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE(psi_sd,psi_rc) 
       
    ELSE
       
       group_size = numprocy * numprocz * numprocr
       ALLOCATE(members(0:numprocx-1))
       members = 0
       ALLOCATE(fftcomm(0:group_size-1))
       fftcomm = 0
       
       ipro = -1
       DO iprocr = 0, numprocr - 1
          DO iprocz = 0, numprocz - 1
             DO iprocy = 0, numprocy - 1
                
                ipro = ipro + 1
                DO iprocx = 0, numprocx - 1
                   ! Look out! We add +1 in the indx routine. Subroutines do not know 
                   ! lower and upper bounds for assumed-shape arrays 
                   members(iprocx)  = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                        iprocr+1,numprocx,numprocy,numprocz,numprocr))
                ENDDO
                
                CALL MPI_GROUP_INCL(group_gl, numprocx, members, group_lc, ierror)
                
                CALL MPI_COMM_CREATE(comm_gl, group_lc, fftcomm(ipro), ierror)
                
             ENDDO
          ENDDO
       ENDDO
       
       ALLOCATE(psi_sd(1:Nxgl))
       ALLOCATE(psi_rc(1:Nxgl))
       dims_phony(1) = Nxgl
       
       ! Work out in which group we are
       ipro = ipy + numprocy * ipz + &
            numprocy * numprocz * ipr
       
       DO ir = 1, Nr
          DO iz = 1, Nz
             DO iy = 1, Ny
                
                ! Set arrays to zero at the begining
                psi_sd = ZERO
                psi_rc = ZERO
                
                DO ix = 1, Nx
                   igl = ix + ipx * Nx
                   psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                ENDDO
                
                ! Send all the elements to all the processors
                CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nxgl, MPI_COMPLEX16,&
                     MPI_SUM, fftcomm(ipro), ierror)
                
                CALL FFT(psi_rc, 1, dims_phony)
                
                ! fftshift the array
                psi_sd = ZERO
                shift = Nxgl - (Nxgl+1)/2
                DO ix = 1, shift
                   psi_sd(ix) = psi_rc(ix+shift)
                ENDDO
                
                DO ix = shift+1, Nxgl
                   psi_sd(ix) = psi_rc(ix-shift)
                ENDDO
                
                
                DO ix = 1, Nx
                   igl = ix + ipx * Nx
                   psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_sd(igl)
                ENDDO
                
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE(psi_sd,psi_rc)
       DEALLOCATE(members,fftcomm)
       
    ENDIF
    
    !********************
    ! FFT in y direction
    !********************    
    
    IF(rank.GE.2) THEN
       
       IF (numprocy .EQ. 1) THEN
          ALLOCATE(psi_sd(1:Nygl))
          ALLOCATE(psi_rc(1:Nygl))
          dims_phony(1) = Nygl
          
          DO ir = 1, Nr
             DO iz = 1, Nz
                DO ix = 1, Nx
                   
                   ! Set arrays to zero at the begining
                   psi_sd = ZERO
                   psi_rc = ZERO
                   
                   DO iy =1, Ny
                      psi_rc(iy) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   ! fftshift the array
                   psi_sd = ZERO
                   shift = Nygl - (Nygl+1)/2
                   DO iy = 1, shift
                      psi_sd(iy) = psi_rc(iy+shift)
                   ENDDO
                   
                   DO iy = shift+1, Nygl
                      psi_sd(iy) = psi_rc(iy-shift)
                   ENDDO
                   
                   ! Copy the result back to the local psik
                   DO iy = 1, Ny
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_sd(iy)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc) 
          
       ELSE
          
          group_size = numprocx * numprocz * numprocr
          ALLOCATE(members(0:numprocy-1))
          members = 0
          ALLOCATE(fftcomm(0:group_size-1))
          fftcomm = 0
          
          ipro = -1
          DO iprocr = 0, numprocr - 1
             DO iprocz = 0, numprocz - 1
                DO iprocx = 0, numprocx - 1
                   
                   ipro = ipro + 1
                   DO iprocy = 0, numprocy - 1
                      ! Look out! We add +1 in the indx routine. Subroutines do not know 
                      ! lower and upper bounds for assumed-shape arrays 
                      members(iprocy)  = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                           iprocr+1,numprocx,numprocy,numprocz,numprocr))
                   ENDDO
                   
                   CALL MPI_GROUP_INCL(group_gl, numprocy, members, group_lc, ierror)
                   
                   CALL MPI_COMM_CREATE(comm_gl, group_lc, fftcomm(ipro), ierror)
                   
                ENDDO
             ENDDO
          ENDDO
          
          ! Work out in which group we are
          ipro = ipx + numprocx * ipz + &
               numprocx * numprocz * ipr
          dims_phony(1) = Nygl
                 
          ALLOCATE(psi_sd(1:Nygl))
          ALLOCATE(psi_rc(1:Nygl))
          
          DO ir = 1, Nr
             DO iz = 1, Nz
                DO ix = 1, Nx
                   
                   ! Set arrays to zero at the begining
                   psi_sd = ZERO
                   psi_rc = ZERO
                   
                   DO iy = 1, Ny
                      igl = iy + ipy * Ny
                      psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   ! Send all the elements to all the processors
                   CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nygl, MPI_COMPLEX16,&
                        MPI_SUM, fftcomm(ipro), ierror)
                   
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   ! fftshift the array
                   psi_sd = ZERO
                   shift = Nygl - (Nygl+1)/2
                   DO iy = 1, shift
                      psi_sd(iy) = psi_rc(iy+shift)
                   ENDDO
                   
                   DO iy = shift+1, Nygl
                      psi_sd(iy) = psi_rc(iy-shift)
                   ENDDO
                   
                   ! Copy the result back to the local psik
                   DO iy = 1, Ny
                      igl = iy + ipy * Ny
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_sd(igl)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc)
          DEALLOCATE(members,fftcomm)
          
       ENDIF
    ENDIF
    
    !********************
    ! FFT in z direction
    !********************    
    
    IF(rank.GE.3) THEN
       
       IF (numprocz .EQ. 1) THEN
          ALLOCATE(psi_sd(1:Nzgl))
          ALLOCATE(psi_rc(1:Nzgl))
          dims_phony(1) = Nzgl
          
          DO ir = 1, Nr
             DO iy = 1, Ny
                DO ix = 1, Nx
                   
                   !Set arrays to zero
                   psi_sd = ZERO
                   psi_rc = ZERO
                   
                   DO iz =1, Nz
                      psi_rc(iz) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   ! fftshift the array
                   psi_sd = ZERO
                   shift = Nzgl - (Nzgl+1)/2
                   DO iz = 1, shift
                      psi_sd(iz) = psi_rc(iz+shift+1)
                   ENDDO
                   
                   DO iz = shift+1, Nzgl
                      psi_sd(iz) = psi_rc(iz-shift)
                   ENDDO
                   
                   ! Copy the result back to the local psik                   
                   DO iz = 1, Nz
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(iz)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc) 
          
       ELSE
          
          group_size = numprocx * numprocy * numprocr
          ALLOCATE(members(0:numprocz-1))
          members = 0
          ALLOCATE(fftcomm(0:group_size-1))
          fftcomm = 0
          
          ipro = -1
          DO iprocr = 0, numprocr - 1
             DO iprocy = 0, numprocy - 1
                DO iprocx = 0, numprocx - 1
                   
                   ipro = ipro + 1
                   DO iprocz = 0, numprocz - 1
                      ! Look out! We add +1 in the indx routine. Subroutines do not know 
                      ! lower and upper bounds for assumed-shape arrays 
                      members(iprocz)  = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                           iprocr+1,numprocx,numprocy,numprocz,numprocr))
                   ENDDO
                   
                   CALL MPI_GROUP_INCL(group_gl, numprocz, members, group_lc, ierror)
                   
                   CALL MPI_COMM_CREATE(comm_gl, group_lc, fftcomm(ipro), ierror)
                ENDDO
             ENDDO
          ENDDO
          
          ! Work out in which group we are
          ipro = ipx + numprocx * ipy + &
               numprocx * numprocy * ipr
          dims_phony(1) = Nzgl
          
          ALLOCATE(psi_sd(1:Nzgl))
          ALLOCATE(psi_rc(1:Nzgl))
          
          DO ir = 1, Nr
             DO iy = 1, Ny
                DO ix = 1, Nx
                   
                   ! Set arrays to zero at the begining
                   psi_sd = ZERO
                   psi_rc = ZERO
                   
                   DO iz = 1, Nz
                      igl = iz + ipz * Nz
                      psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   ! Send all the elements to all the processors
                   CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nzgl, MPI_COMPLEX16,&
                        MPI_SUM, fftcomm(ipro), ierror)
               
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   ! fftshift the array
                   psi_sd = ZERO
                   shift = Nzgl - (Nzgl+1)/2
                   DO iz = 1, shift
                      psi_sd(iz) = psi_rc(iz+shift)
                   ENDDO
                   
                   DO iz = shift+1, Nzgl
                      psi_sd(iz) = psi_rc(iz-shift)
                   ENDDO
                   
                   ! Copy the result back to the local psik                   
                   DO iz = 1, Nz
                      igl = iz + ipz * Nz
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_sd(igl)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc)
          DEALLOCATE(members,fftcomm)
          
       ENDIF
    ENDIF
    
    !********************
    ! FFT in r direction
    !********************
    
    IF(rank.GE.4) THEN
       
       IF (numprocr .EQ. 1) THEN
          ALLOCATE(psi_sd(1:Nrgl))
          ALLOCATE(psi_rc(1:Nrgl))
          dims_phony(1) = Nrgl
          
          DO iz = 1, Nz
             DO iy = 1, Ny
                DO ix = 1, Nx
                   
                   DO ir =1, Nr
                      psi_rc(ir) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   !! fftshift the array
                   !psi_sd = ZERO
                   !shift = Nrgl - (Nrgl+1)/2
                   !DO ir = 1, shift
                   !   psi_sd(ir) = psi_rc(ir+shift+1)
                   !ENDDO
                   !
                   !DO ir = shift+1, Nrgl
                   !   psi_sd(ir) = psi_rc(ir-shift)
                   !ENDDO
                   
                   ! Copy the result back to the local psik
                   DO ir = 1, Nr
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(ir)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc) 
          
       ELSE
          
          group_size = numprocx * numprocy * numprocz
          ALLOCATE(members(0:numprocr-1))
          members = 0
          ALLOCATE(fftcomm(0:group_size-1))
          fftcomm = 0
          
          ipro = -1
          DO iprocz = 0, numprocz - 1
             DO iprocy = 0, numprocy - 1
                DO iprocx = 0, numprocx - 1
                   
                   ipro = ipro + 1
                   DO iprocr = 0, numprocr - 1
                      ! Look out! We add +1 in the indx routine. Subroutines do not know 
                      ! lower and upper bounds for assumed-shape arrays 
                      members(iprocr)  = grid_addresses(indx(iprocx+1,iprocy+1,iprocz+1, &
                           iprocr+1,numprocx,numprocy,numprocz,numprocr))
                   ENDDO
                   
                   CALL MPI_GROUP_INCL(group_gl, numprocr, members, group_lc, ierror)
                   
                   CALL MPI_COMM_CREATE(comm_gl, group_lc, fftcomm(ipro), ierror)
                   
                ENDDO
             ENDDO
          ENDDO
          
          ! Work out in which group we are
          ipro = ipx + numprocx * ipy + &
               numprocx * numprocy * ipz
          
          ALLOCATE(psi_sd(1:Nrgl))
          ALLOCATE(psi_rc(1:Nrgl))
          
          DO iz = 1, Nz
             DO iy = 1, Ny
                DO ix = 1, Nx
                   
                   DO ir = 1, Nr
                      igl = ir + ipr * Nr
                      psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                   ENDDO
                   
                   ! Send all the elements to all the processors
                   CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nrgl, MPI_COMPLEX16,&
                        MPI_SUM, fftcomm(ipro), ierror)
                   
                   dims_phony(1) = Nrgl
                   CALL FFT(psi_rc, 1, dims_phony)
                   
                   !! fftshift the array
                   !psi_sd = ZERO
                   !shift = Nrgl - (Nrgl+1)/2
                   !DO ir = 1, shift
                   !   psi_sd(ir) = psi_rc(ir+shift+1)
                   !ENDDO
                   !
                   !DO ir = shift+1, Nrgl
                   !   psi_sd(ir) = psi_rc(ir-shift)
                   !ENDDO
                   
                   ! Copy the result back to the local psik                   
                   DO ir = 1, Nr
                      igl = ir + ipr * Nr
                      psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(igl)
                   ENDDO
                   
                ENDDO
             ENDDO
          ENDDO
          
          DEALLOCATE(psi_sd,psi_rc)
          DEALLOCATE(members,fftcomm)
          
       ENDIF
    ENDIF
    
  END SUBROUTINE FourierTransform_Parallel
  
#endif
  
  !----------------------------------------------------!
  
#if _USE_MKL
  
  SUBROUTINE FFT(in, rank, dims)
    
    USE MKL_DFTI
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    COMPLEX(dp), INTENT(INOUT)       :: in(:)
    INTEGER, INTENT(IN)              :: rank
    INTEGER, INTENT(IN)              :: dims(:)
    
    ! Execution status
    COMPLEX(dp), ALLOCATABLE         :: out(:, :, :, :)
    integer                          :: status = 0, ignored_status
    type(DFTI_DESCRIPTOR), POINTER   :: hand
    INTEGER                          :: Nx, Ny, Nz, Nr
    
    !------------------------------------------------------------------------!

    ! Create DFTI descriptor
    hand => null()
    
    IF (rank .EQ. 1) THEN
       Nx = dims(1)
       Ny = 1
       Nz = 1
       Nr = 1       
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = 1
       Nr = 1
    ELSEIF (rank .EQ. 3) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = dims(3)
       Nr = 1
    ELSEIF (rank .EQ. 4) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = dims(3)
       Nr = dims(4)
    ELSE
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF
    
    
    ALLOCATE(out(Nx,Ny,Nz,Nr))
    out  = RESHAPE(in,(/ Nx, Ny, Nz, Nr /))
    
    status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 4,[Nx,Ny,Nz,Nr])
    IF (0 /= status) THEN
       WRITE(*,*) 'Error creating descriptor!' 
       STOP
    ENDIF
    
    ! Commit DFTI descriptor
    status = DftiCommitDescriptor(hand)
    IF (0 /= status) THEN
       WRITE(*,*) 'Error committing the descriptor!'
       STOP
    ENDIF
    
    status = DftiComputeForward(hand, out(:,1,1,1))
    
    status = DftiFreeDescriptor(hand)
    
    in  = RESHAPE(out,(/ Nx*Ny*Nz*Nr /) )
    
    DEALLOCATE(out)
    
  END SUBROUTINE FFT
  
  !----------------------------------------------!
  
#elif _USE_FFTW
  
  SUBROUTINE FFT(in, rank, dims)

    USE FFTW
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    COMPLEX(dp), INTENT(INOUT)          :: in(:)
    INTEGER, INTENT(IN)                 :: rank
    INTEGER, INTENT(IN)                 :: dims(:)
    
    ! Execution status
    type(C_PTR)                         :: plan
    type(C_PTR)                         :: p_in, p_out, p_dims
    complex(C_DOUBLE_COMPLEX), pointer  :: in4d(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer  :: out(:, :, :, :)
    INTEGER(C_INT), pointer             :: dims_for_fftw(:)    			      
    INTEGER                             :: Nx, Ny, Nz, Nr
    
    !------------------------------------------------------------------------!
    
    IF (rank .EQ. 1) THEN
       Nx = dims(1)
       Ny = 1
       Nz = 1
       Nr = 1       
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = 1
       Nr = 1
    ELSEIF (rank .EQ. 3) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = dims(3)
       Nr = 1
    ELSEIF (rank .EQ. 4) THEN
       Nx = dims(1)
       Ny = dims(2)
       Nz = dims(3)
       Nr = dims(4)
    ELSE
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF
    
    ! Allocate arrays
    p_in = fftw_alloc_complex(int(Nx * Ny * Nz * Nr, C_SIZE_T))  
    p_out = fftw_alloc_complex(int(Nx * Ny * Nz * Nr, C_SIZE_T))
    p_dims = fftw_alloc_complex(int(4,C_SIZE_T))
    
    CALL c_f_pointer(p_in, in4d, [Nx, Ny, Nz, Nr])
    CALL c_f_pointer(p_out, out, [Nx, Ny, Nz, Nr])
    CALL c_f_pointer(p_dims,dims_for_fftw,[4])
    
    in4d  = RESHAPE(in,(/ Nx, Ny, Nz, Nr /))
    dims_for_fftw = dims(4:1:-1)
    
    plan = fftw_plan_dft(4,dims_for_fftw,in4d,out,FFTW_FORWARD,FFTW_ESTIMATE)
    
    ! Execute it!!
    CALL fftw_execute_dft(plan, in4d, out)
    
    ! Free plan
    CALL fftw_destroy_plan(plan)
    
    in  = RESHAPE(out,(/ Nx*Ny*Nz*Nr /) )
    
    CALL fftw_free(p_in)
    CALL fftw_free(p_out)
    CALL fftw_free(p_dims)
    
  END SUBROUTINE FFT
  
  !----------------------------------------------!  
  
#else
  
  !----------------------------------------------!  
  
  SUBROUTINE FFT
    
    IMPLICIT NONE
    
    WRITE(*,*) 'There is no FFT library linked.'
    
  END SUBROUTINE FFT
  
  !----------------------------------------------!  
  
#endif
    
  !----------------------------------------------!
  
  FUNCTION HankelTransform( in, Nrho, rpts, krpts, jacobian )
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    COMPLEX(dp), INTENT(IN)          :: in(:)
    COMPLEX(dp)                      :: HankelTransform
    INTEGER, INTENT(IN)              :: Nrho
    REAL(dp), INTENT(IN)             :: rpts(:)
    REAL(dp), INTENT(IN)             :: krpts
    REAL(dp), INTENT(IN)             :: jacobian(:)
    INTEGER                          :: irho
    
    !------------------------------------------------------------------------!
    
    HankelTransform = 0.0_dp
    
    DO irho = 1, Nrho
       HankelTransform  = HankelTransform + jacobian(irho) * &
            in(irho) * BESSEL_J0(krpts * rpts(irho))
    ENDDO
    
  END FUNCTION HankelTransform
  
  !----------------------------------------------!
  
  
  !----------------------------------------------!
  
  SUBROUTINE create_mask4D( mask, x_ax, y_ax, z_ax, r_ax, &
       inner_rad, outer_rad )
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    REAL(dp), INTENT(OUT)               :: mask(:, :, :, :)
    REAL(dp), INTENT(IN)                :: inner_rad
    REAL(dp), INTENT(IN)                :: outer_rad
    REAL(dp), INTENT(IN)                :: x_ax(:)
    REAL(dp), INTENT(IN)                :: y_ax(:)
    REAL(dp), INTENT(IN)                :: z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:)    

    INTEGER, ALLOCATABLE                :: dims(:) 
    REAL(dp)                            :: sigma
    REAL(dp)                            :: rad
    INTEGER                             :: ix, iy, iz, ir
    
    !------------------------------------------------------------------------!
    
    ALLOCATE(dims(1:4))
    
    dims = SHAPE(mask)
    
    IF (SIZE(x_ax).NE.dims(1) &
         .OR. SIZE(y_ax).NE.dims(2) &
         .OR. SIZE(z_ax).NE.dims(3) &
         .OR. SIZE(r_ax).NE.dims(4) ) THEN
       WRITE(*,*) 'Dimensions of the input arrays do not match.'
       STOP
    ENDIF
    
    ! Calculate the sigma
    sigma    = (outer_rad - inner_rad) / SQRT( - LOG( 1E-8 ) )
    
    ! Create the mask
    DO ir = 1, dims(4)
       DO iz = 1, dims(3)
          DO iy = 1, dims(2)
             DO ix = 1, dims(1)
                
                rad = SQRT(x_ax(ix) * x_ax(ix) + &
                     y_ax(iy) * y_ax(iy) + &
                     z_ax(iz) * z_ax(iz) + &
                     r_ax(ir) * r_ax(ir) )
                
                IF (rad.LE.inner_rad) THEN 
                   mask(ix,iy,iz,ir) = 0.0E-15_dp
                ELSEIF (rad.GE.inner_rad .AND. rad.LT.outer_rad) THEN
                   mask(ix,iy,iz,ir) = 1.0 - EXP( - ( (rad - inner_rad) / sigma )**2 )
                ELSE
                   mask(ix,iy,iz,ir) = 1.0_dp
                ENDIF
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(dims)
    
  END SUBROUTINE create_mask4D

  !----------------------------------------!
  
  SUBROUTINE create_mask3D( mask, x_ax, y_ax, z_ax, &
       inner_rad, outer_rad )
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    REAL(dp), INTENT(INOUT)             :: mask(:, :, :)
    REAL(dp), INTENT(IN)                :: inner_rad
    REAL(dp), INTENT(IN)                :: outer_rad
    REAL(dp), INTENT(IN)                :: x_ax(:)
    REAL(dp), INTENT(IN)                :: y_ax(:)
    REAL(dp), INTENT(IN)                :: z_ax(:)
    
    INTEGER, ALLOCATABLE                :: dims(:)
    REAL(dp)                            :: sigma
    REAL(dp)                            :: rad
    INTEGER                             :: ix, iy, iz
    
    !------------------------------------------------------------------------!
    
    ALLOCATE(dims(1:3))
    
    dims = SHAPE(mask)
    IF (SIZE(x_ax).NE.dims(1) &
         .OR. SIZE(y_ax).NE.dims(2) &
         .OR. SIZE(z_ax).NE.dims(3)) THEN
       WRITE(*,*) 'Dimensions of the input arrays do not match.' 
       STOP
    ENDIF
    
    ! Calculate the sigma
    sigma    = (outer_rad - inner_rad) / SQRT( - LOG( 1E-8 ) )
    
    ! Create the mask
    DO iz = 1, dims(3)
       DO iy = 1, dims(2)
          DO ix = 1, dims(1)
             
             rad = SQRT(x_ax(ix) * x_ax(ix) + &
                  y_ax(iy) * y_ax(iy) + &
                  z_ax(iz) * z_ax(iz) )
             
             IF (rad.LE.inner_rad) THEN 
                mask(ix,iy,iz) = 0.0_dp
             ELSEIF (rad.GE.inner_rad .AND. rad.LT.outer_rad) THEN
                mask(ix,iy,iz) = 1.0 - EXP( - ( (rad - inner_rad) / sigma )**2 )
             ELSE
                mask(ix,iy,iz) = 1.0_dp
             ENDIF
             
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(dims)
    
    
  END SUBROUTINE create_mask3D
  
  !----------------------------------------!
  
  SUBROUTINE create_mask2D( mask, x_ax, y_ax, &
       inner_rad, outer_rad )
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    REAL(dp), INTENT(INOUT)             :: mask(:, :)
    REAL(dp), INTENT(IN)                :: inner_rad
    REAL(dp), INTENT(IN)                :: outer_rad
    REAL(dp), INTENT(IN)                :: x_ax(:)
    REAL(dp), INTENT(IN)                :: y_ax(:)
    
    INTEGER, ALLOCATABLE                :: dims(:)
    REAL(dp)                            :: sigma
    REAL(dp)                            :: rad
    INTEGER                             :: ix, iy, iz
    
    !------------------------------------------------------------------------!
    
    ALLOCATE(dims(1:2))
    
    dims = SHAPE(mask)
    IF (SIZE(x_ax).NE.dims(1) &
         .OR. SIZE(y_ax).NE.dims(2)) THEN
       WRITE(*,*) 'Dimensions of the input arrays do not match.'
       STOP
    ENDIF
    
    ! Calculate the sigma
    sigma    = (outer_rad - inner_rad) / SQRT( - LOG( 1E-8 ) )
    
    ! Create the mask
    DO iy = 1, dims(2)
       DO ix = 1, dims(1)
          
          rad = SQRT(x_ax(ix) * x_ax(ix) + &
               y_ax(iy) * y_ax(iy))
          
          IF (rad.LE.inner_rad) THEN 
             mask(ix,iy) = 0.0_dp
          ELSEIF (rad.GE.inner_rad .AND. rad.LT.outer_rad) THEN
             mask(ix,iy) = 1.0 - EXP( - ( (rad - inner_rad) / sigma )**2 )
          ELSE
             mask(ix,iy) = 1.0_dp
          ENDIF
          
       ENDDO
    ENDDO
    
    DEALLOCATE(dims)
    
    
  END SUBROUTINE create_mask2D
  
  !----------------------------------------!
  
  SUBROUTINE create_correlated_mask3D1D(mask, x_ax, y_ax, &
       z_ax, r_ax, inner_erad, outer_erad, inner_nrad, outer_nrad)
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    REAL(dp), INTENT(INOUT)             :: mask(:, :, :, :)
    REAL(dp), INTENT(IN)                :: inner_erad
    REAL(dp), INTENT(IN)                :: outer_erad
    REAL(dp), INTENT(IN)                :: inner_nrad
    REAL(dp), INTENT(IN)                :: outer_nrad
    REAL(dp), INTENT(IN)                :: x_ax(:)
    REAL(dp), INTENT(IN)                :: y_ax(:)
    REAL(dp), INTENT(IN)                :: z_ax(:)
    REAL(dp), INTENT(IN)                :: r_ax(:) 
    
    INTEGER, ALLOCATABLE                :: dims(:)
    REAL(dp)                            :: sigma
    REAL(dp)                            :: erad, nrad
    INTEGER                             :: ix, iy, iz, ir
    
    !------------------------------------------------------------------------!
    
    ALLOCATE(dims(1:4))
    
    dims = SHAPE(mask)
    
    IF (SIZE(x_ax).NE.dims(1) &
         .OR. SIZE(y_ax).NE.dims(2) &
         .OR. SIZE(z_ax).NE.dims(3) &
         .OR. SIZE(r_ax).NE.dims(4) ) THEN
       WRITE(*,*) 'Dimensions of the input arrays do not match.'
       STOP
    ENDIF
    
    ! Calculate the sigma
    sigma = ( (outer_erad / inner_erad)**2 + (outer_nrad/inner_nrad)**2 &
         - 1.0_dp ) / SQRT( - LOG( 1.0E-50_dp ) )
    
    ! Create the mask
    DO ir = 1, dims(4)
       DO iz = 1, dims(3)
          DO iy = 1, dims(2)
             DO ix = 1, dims(1)
                
                erad = SQRT(x_ax(ix) * x_ax(ix) + &
                     y_ax(iy) * y_ax(iy) + &
                     z_ax(iz) * z_ax(iz)  )
                
                nrad = r_ax(ir)
                
                IF( (erad/inner_erad)**2 + (nrad/inner_nrad)**2 .LT. 1.0_dp ) THEN
                   mask(ix,iy,iz,ir) = 0.0E-15_dp
                   
                ELSEIF( ( (erad/inner_erad)**2 + (nrad/inner_nrad)**2 .GT. 1.0_dp ) &
                     .AND. &
                     ( (erad/outer_erad)**2 + (nrad/outer_nrad)**2 .LT. 1.0_dp ) ) THEN
                   mask(ix,iy,iz,ir) = 1.0_dp - EXP( - ( ( (erad/inner_erad)**2 + (nrad/inner_nrad)**2 -1.0_dp )/sigma )**2 )
                   
                ELSE
                   mask(ix,iy,iz,ir) = 1.0_dp
                ENDIF
                
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(dims)
    
    
  END SUBROUTINE create_correlated_mask3D1D
  
  !----------------------------------------!
  
END MODULE fourier
