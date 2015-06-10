MODULE fourier
  
#if _COM_MPI

  USE MPI
  
#endif
  
  USE omp_lib
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC                         :: FourierTransform
  PUBLIC                         :: FFT
  PUBLIC                         :: HankelTransform
  PUBLIC                         :: create_mask
  PUBLIC                         :: indx
  
  INTERFACE create_mask
     MODULE PROCEDURE create_mask4D
     MODULE PROCEDURE create_mask3D
     MODULE PROCEDURE create_mask2D
  END INTERFACE create_mask
  
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
       dims_global, mpi_rank, communicators, ipgrid )

    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: psi(:)
    COMPLEX(dp), INTENT(OUT)      :: psik(:)
    INTEGER, INTENT(IN)           :: rank
    INTEGER, INTENT(IN)           :: dims_local(:)
    INTEGER, INTENT(IN)           :: dims_global(:)
    INTEGER, INTENT(IN)           :: mpi_rank
    INTEGER, INTENT(IN)           :: communicators(:)
    INTEGER, INTENT(IN)           :: ipgrid(:)

    COMPLEX(dp), ALLOCATABLE      :: psi_sd(:)
    COMPLEX(dp), ALLOCATABLE      :: psi_rc(:)
    INTEGER                       :: dims_phony(4)
    INTEGER                       :: ierror
    INTEGER                       :: Nx, Ny, Nz, Nr
    INTEGER                       :: Nxgl, Nygl, Nzgl, Nrgl
    INTEGER                       :: ix, iy, iz, ir
    INTEGER                       :: igl
    
    !***************************************************!

    ! Copy the original array
    psik = psi
    
    ! Assign dimensions
    IF (rank .EQ. 1) THEN
       Nx = dims_local(1)
       Ny = 1
       Nz = 1
       Nr = 1

       Nxgl = dims_global(1)
       Nygl = 1
       Nzgl = 1
       Nrgl = 1
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = 1
       Nr = 1

       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = 1
       Nrgl = 1
    ELSEIF (rank .EQ. 3) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = 1

       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = 1
    ELSEIF (rank .EQ. 4) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = dims_local(4)
       
       Nxgl = dims_global(1)
       Nygl = dims_global(2)
       Nzgl = dims_global(3)
       Nrgl = dims_global(4)
    ELSE
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF
    
    dims_phony = 1

    ALLOCATE(psi_sd(1:Nxgl))
    ALLOCATE(psi_rc(1:Nxgl))

    DO ir = 1, Nr
       DO iz = 1, Nz
          DO iy = 1, Ny
             DO ix = 1, Nx
                igl = ix + ipgrid(1) * Nx
                psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
             ENDDO
             
             ! Send all the elements to all the processors
             CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nxgl, MPI_COMPLEX16,&
                  MPI_SUM, communicators(1), ierror)
             
             dims_phony(1) = Nxgl
             
             CALL FFT(psi_rc, 1, dims_phony)
             
             DO ix = 1, Nx
                igl = ix + ipgrid(1) * Nx
                psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(igl)
             ENDDO
             
          ENDDO
       ENDDO
    ENDDO
    
    DEALLOCATE(psi_sd)
    DEALLOCATE(psi_rc)
    

    IF(rank.GE.2) THEN
       
       ALLOCATE(psi_sd(1:Nygl))
       ALLOCATE(psi_rc(1:Nygl))
       
       DO ir = 1, Nr
          DO iz = 1, Nz
             DO ix = 1, Nx
                
                DO iy = 1, Ny
                   igl = iy + ipgrid(2) * Ny
                   psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                ENDDO
                
                ! Send all the elements to all the processors
                CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nygl, MPI_COMPLEX16,&
                     MPI_SUM, communicators(2), ierror)
                
                dims_phony(1) = Nygl
                CALL FFT(psi_rc, 1, dims_phony)
                
                DO iy = 1, Ny
                   igl = iy + ipgrid(2) * Ny
                   psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(igl)
                ENDDO
                
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE(psi_sd)
       DEALLOCATE(psi_rc)
       
    ENDIF
    
    
    IF(rank.GE.3) THEN
       
       ALLOCATE(psi_sd(1:Nzgl))
       ALLOCATE(psi_rc(1:Nzgl))
       
       DO ir = 1, Nr
          DO iy = 1, Ny
             DO ix = 1, Nx
                
                DO iz = 1, Nz
                   igl = iz + ipgrid(3) * Nz
                   psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                ENDDO
                
                ! Send all the elements to all the processors
                CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nzgl, MPI_COMPLEX16,&
                     MPI_SUM, communicators(3), ierror)
                
                dims_phony(1) = Nzgl
                CALL FFT(psi_rc, 1, dims_phony)
                
                DO iz = 1, Nz
                   igl = iz + ipgrid(3) * Nz
                   psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(igl)
                ENDDO
                
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE(psi_sd)
       DEALLOCATE(psi_rc)
       
    ENDIF
    
    
    IF(rank.GE.4) THEN
       
       ALLOCATE(psi_sd(1:Nrgl))
       ALLOCATE(psi_rc(1:Nrgl))
       
       DO iz = 1, Nz
          DO iy = 1, Ny
             DO ix = 1, Nx
                
                DO ir = 1, Nr
                   igl = ir + ipgrid(4) * Nr
                   psi_sd(igl) = psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr))
                ENDDO
                
                ! Send all the elements to all the processors
                CALL MPI_ALLREDUCE(psi_sd, psi_rc, Nrgl, MPI_COMPLEX16,&
                     MPI_SUM, communicators(4), ierror)
                
                dims_phony(1) = Nrgl
                CALL FFT(psi_rc, 1, dims_phony)
                
                DO ir = 1, Nr
                   igl = ir + ipgrid(4) * Nr
                   psik(indx(ix,iy,iz,ir,Nx,Ny,Nz,Nr)) = psi_rc(igl)
                ENDDO
                
             ENDDO
          ENDDO
       ENDDO
       
       DEALLOCATE(psi_sd)
       DEALLOCATE(psi_rc)
       
    ENDIF
    
    
  END SUBROUTINE FourierTransform_Parallel
  
#endif
  
  !----------------------------------------------------!
  
  FUNCTION indx(ix, iy, iz, ir, Nx, Ny, Nz, Nr)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)      :: ix, iy, iz, ir
    INTEGER, INTENT(IN)      :: Nx, Ny, Nz, Nr
    INTEGER                  :: indx
    
    indx = Nx*Ny*Nz*(ir-1) + Nx*Ny*(iz-1) &
         + Nx*(iy-1) + ix
    
  END FUNCTION indx
  
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
    type(C_PTR)                         :: plan, p
    complex(C_DOUBLE_COMPLEX), pointer  :: in4d(:, :, :, :)
    complex(C_DOUBLE_COMPLEX), pointer  :: out(:, :, :, :)
    int(C_INT)
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
    dims_for_fftw = dims

    plan = fftw_plan_dft(4,dims_for_fftw,in4d,out,FFTW_FORWARD,FFTW_ESTIMATE)
    
    ! Execute it!!
    CALL fftw_execute_dft(plan, in4d, out)
    
    ! Free plan
    CALL fftw_destroy_plan(plan)
    
    in  = RESHAPE(out,(/ Nx*Ny*Nz*Nr /) )
    
    fftw_free(p_in)
    fftw_free(p_out)
    fftw_free(p_dims)
    
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
  
!!$#if _COMM_MPI  
!!$  
!!$  SUBROUTINE FourierTransform_parallel(in, rank, dims_local, dims_global)
!!$    
!!$    IMPLICIT NONE
!!$    
!!$    !--Program variables-----------------------------------------------------!
!!$    
!!$    COMPLEX(dp), INTENT(INOUT)          :: in(:)
!!$    INTEGER, INTENT(IN)                 :: rank
!!$    INTEGER, INTENT(IN)                 :: dims_local(:)
!!$    INTEGER, INTENT(IN)                 :: dims_global(:)
!!$    
!!$    ! Execution status
!!$    COMPLEX(dp), ALLOCATABLE            :: out1D(:)
!!$    COMPLEX(dp), ALLOCATABLE            :: out2D(:, :)
!!$    COMPLEX(dp), ALLOCATABLE            :: out3D(:, :, :) 
!!$    COMPLEX(dp), ALLOCATABLE            :: out4D(:, :, :, :)
!!$    COMPLEX(dp), ALLOCATABLE            :: work(:)
!!$    integer                             :: status = 0, ignored_status
!!$    integer                             :: size
!!$    type(DFTI_DESCRIPTOR_DM), POINTER   :: desc
!!$    INTEGER                             :: Nx, Ny, Nz, Nr
!!$    INTEGER              ::  NX_OUT,START_X,START_X_OUT
!!$    
!!$    !------------------------------------------------------------------------!
!!$    
!!$    write(*,*) 'Hola!'
!!$    write(*,*) 'rank!: ',rank
!!$
!!$    IF (rank .EQ. 1) THEN
!!$       Nx = dims_local(1)
!!$       ALLOCATE(out1D(1:Nx))
!!$       out1D = RESHAPE(in,(/ Nx /))
!!$    ELSEIF (rank .EQ. 2) THEN
!!$       Nx = dims_local(1)
!!$       Ny = dims_local(2)
!!$       ALLOCATE(out2D(1:Nx,1:Ny))
!!$       out2D = RESHAPE(in,(/ Nx, Ny /))
!!$    ELSEIF (rank .EQ. 3) THEN
!!$       write(*,*) 'dentro del if!'
!!$       Nx = dims_local(1)
!!$       Ny = dims_local(2)
!!$       Nz = dims_local(3)
!!$       ALLOCATE(out3D(1:Nx,1:Ny,1:Nz))
!!$       out3D = RESHAPE(in,(/ Nx, Ny, Nz /))
!!$       write(*,*) 'Despues'
!!$   ELSEIF (rank .EQ. 4) THEN
!!$       Nx = dims_local(1)
!!$       Ny = dims_local(2)
!!$       Nz = dims_local(3)
!!$       Nr = dims_local(4)
!!$       ALLOCATE(out4D(1:Nx,1:Ny,1:Nz,1:Nr))
!!$       out4D = RESHAPE(in,(/ Nx, Ny, Nz, Nr /))
!!$    ELSE
!!$       write(*,*) 'Am I here?'
!!$       write(*,*) 'rank',rank
!!$       WRITE(*,*) 'No option implemented for these dimensions'
!!$       STOP
!!$    ENDIF
!!$
!!$    ! Allocate memory for the descriptor by calling DftiCreateDescriptorDM
!!$    
!!$    status = DftiCreateDescriptorDM(MPI_COMM_WORLD ,desc,DFTI_DOUBLE,        &
!!$         &                                DFTI_COMPLEX,rank,dims_global) 
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Create parallel descriptor failed')
!!$    END IF
!!$    
!!$    ! Obtain some values of configuration parameters by calls to
!!$    !  DftiGetValueDM
!!$    
!!$    status = DftiGetValueDM(desc,CDFT_LOCAL_SIZE,size)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Get CDFT local size failed')
!!$    END IF
!!$    
!!$    status = DftiGetValueDM(desc,CDFT_LOCAL_NX,NX)
!!$    IF (STATUS .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Get CDFT local maxpts failed')
!!$    ENDIF
!!$    
!!$    STATUS = DftiGetValueDM(desc,CDFT_LOCAL_X_START,START_X)
!!$    IF (STATUS .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Get CDFT local x start failed')      
!!$    END IF
!!$    
!!$    status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_NX,NX_OUT)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Get CDFT local out failed')
!!$    END IF
!!$    
!!$    status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_X_START,START_X_OUT)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Get CDFT local out x start failed')
!!$    ENDIF
!!$    
!!$
!!$    write(*,*) 'size: ',size
!!$    write(*,*) 'Nx: ',NX
!!$    write(*,*) 'Start_x: ',start_x
!!$    write(*,*) 'Nx_out: ',nx_out
!!$    write(*,*) 'start_x_out: ',start_x_out
!!$
!!$    ! Allocate auxiliary array
!!$    ALLOCATE(work(size))
!!$    
!!$    
!!$    ! Specify a value(s) of configuration parameters by a call(s) to
!!$    !  DftiSetValueDM
!!$    
!!$    status = DftiSetValueDM(desc,CDFT_WORKSPACE,work)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Set CDFT workspace failed')
!!$    END IF
!!$    
!!$    !     Perform initialization that facilitates DFT computation by a call to
!!$    !        DftiCommitDescriptorDM
!!$    
!!$    status = DftiCommitDescriptorDM(desc)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Commit descriptor failed')
!!$    END IF
!!$    
!!$    ! Create arrays for local parts of input and output data
!!$    !  (if it is needed) and fill the local part of input data with
!!$    !  values (for more information, see Distributing Data among Processes)
!!$    
!!$      ROOTRANK=0
!!$      STATUS = MKL_CDFT_SCATTERDATA_D(COMM,ROOTRANK,ELEMENTSIZE,2,      &
!!$     &                                LENGTHS,X_IN,NX,START_X,LOCAL)
!!$      IF (ADVANCED_DATA_PRINT .AND. (MPI_RANK .EQ. 0)) THEN
!!$         PRINT *, 'Scatter=',STATUS
!!$      END IF
!!$      IF (STATUS .NE. DFTI_NO_ERROR) THEN
!!$         IF (MPI_RANK .EQ. 0) THEN
!!$            CALL Dfti_Example_Status_Print(STATUS)
!!$            PRINT *, 'TEST FAILED'
!!$         END IF
!!$         FAILURE = .TRUE.
!!$         GOTO 101
!!$      END IF
!!$    
!!$    
!!$    ! Compute the transform by calling
!!$    !  DftiComputeForwardDM or DftiComputeBackward
!!$    !  (Compute Forward transform)
!!$    
!!$    STATUS = DftiComputeForwardDM(desc,in)
!!$    IF (STATUS .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Compute forward cdft failed')
!!$    END IF
!!$    
!!$    
!!$    ! Release memory allocated for a descriptor by a call to
!!$    !  DftiFreeDescriptorDM (Free DftiDM descriptor)
!!$    
!!$    status = DftiFreeDescriptorDM(desc)
!!$    IF (status .NE. DFTI_NO_ERROR) THEN
!!$       CALL tddm_stop('Free descriptor failed')
!!$    END IF
!!$    
!!$    
!!$    DEALLOCATE(work)
!!$    
!!$  END SUBROUTINE FourierTransform_parallel
!!$  
!!$#endif
  
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
    
    REAL(dp), INTENT(INOUT)             :: mask(:, :, :, :)
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
                   mask(ix,iy,iz,ir) = 0.0_dp
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
  
  
END MODULE fourier
