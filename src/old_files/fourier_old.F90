MODULE fourier
  
#if _COMM_MPI

  USE MPI
  
#if _USE_MKL_CLUSTER
  USE MKL_CDFT
#endif
  
#if _USE_FFTW_CLUSTER
  USE MKL_CDFT
#endif

#endif
  
#if _USE_MKL

  USE MKL_DFTI

#elif _USE_FFTW

  USE FFTW

#endif
  
  USE omp_lib
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC                         :: transform_to_momentum_space
  PUBLIC                         :: FourierTransform
  PUBLIC                         :: HankelTransform
  PUBLIC                         :: create_mask
  
#if _COMM_MPI
  PUBLIC                         :: FourierTransform_cluster
#endif
  
  INTERFACE create_mask
     MODULE PROCEDURE create_mask4D
     MODULE PROCEDURE create_mask3D
     MODULE PROCEDURE create_mask2D
  END INTERFACE create_mask
  
  
CONTAINS
  
  SUBROUTINE transform_to_momentum_space(psi, psik, rank, dims_local, dims_global )
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)       :: psi(:)
    COMPLEX(dp), INTENT(OUT)      :: psik(:)
    INTEGER, INTENT(IN)           :: rank
    INTEGER, INTENT(IN)           :: dims_local(:)
    INTEGER, INTENT(IN)           :: dims_global(:)
    
    !***************************************************!
    
    !---------------------------!
    ! Copy the original array   !
    !---------------------------!
    write(*,*) 'Lllegamos!'
    psik = psi
    
    !-------------------------------------------------!
    ! Transform the wavefunction into momentum space. !
    !-------------------------------------------------!
    
    WRITE(*,*) 'Transforming to momentum space...'
    
    ! Compute forward transform
    !CALL FourierTransform_cluster(psik,rank,dims_local,dims_global)
    
   
  END SUBROUTINE transform_to_momentum_space
  
  !----------------------------------------------------!
  
  SUBROUTINE FourierTransform(in, rank, dims)
    
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
    
    write(*,*) 'Estamos en serial'
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
       write(*,*) 'Serial'
       write(*,*) 'rank: ',rank
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
    
  END SUBROUTINE FourierTransform
  
  !----------------------------------------------!

#if _COMM_MPI  
  
  SUBROUTINE FourierTransform_parallel(in, rank, dims_local, dims_global)
    
    IMPLICIT NONE
    
    !--Program variables-----------------------------------------------------!
    
    COMPLEX(dp), INTENT(INOUT)          :: in(:)
    INTEGER, INTENT(IN)                 :: rank
    INTEGER, INTENT(IN)                 :: dims_local(:)
    INTEGER, INTENT(IN)                 :: dims_global(:)
    
    ! Execution status
    COMPLEX(dp), ALLOCATABLE            :: out1D(:)
    COMPLEX(dp), ALLOCATABLE            :: out2D(:, :)
    COMPLEX(dp), ALLOCATABLE            :: out3D(:, :, :) 
    COMPLEX(dp), ALLOCATABLE            :: out4D(:, :, :, :)
    COMPLEX(dp), ALLOCATABLE            :: work(:)
    integer                             :: status = 0, ignored_status
    integer                             :: size
    type(DFTI_DESCRIPTOR_DM), POINTER   :: desc
    INTEGER                             :: Nx, Ny, Nz, Nr
    INTEGER              ::  NX_OUT,START_X,START_X_OUT
    
    !------------------------------------------------------------------------!
    
    write(*,*) 'Hola!'
    write(*,*) 'rank!: ',rank

    IF (rank .EQ. 1) THEN
       Nx = dims_local(1)
       ALLOCATE(out1D(1:Nx))
       out1D = RESHAPE(in,(/ Nx /))
    ELSEIF (rank .EQ. 2) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       ALLOCATE(out2D(1:Nx,1:Ny))
       out2D = RESHAPE(in,(/ Nx, Ny /))
    ELSEIF (rank .EQ. 3) THEN
       write(*,*) 'dentro del if!'
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       ALLOCATE(out3D(1:Nx,1:Ny,1:Nz))
       out3D = RESHAPE(in,(/ Nx, Ny, Nz /))
       write(*,*) 'Despues'
   ELSEIF (rank .EQ. 4) THEN
       Nx = dims_local(1)
       Ny = dims_local(2)
       Nz = dims_local(3)
       Nr = dims_local(4)
       ALLOCATE(out4D(1:Nx,1:Ny,1:Nz,1:Nr))
       out4D = RESHAPE(in,(/ Nx, Ny, Nz, Nr /))
    ELSE
       write(*,*) 'Am I here?'
       write(*,*) 'rank',rank
       WRITE(*,*) 'No option implemented for these dimensions'
       STOP
    ENDIF

    ! Allocate memory for the descriptor by calling DftiCreateDescriptorDM
    
    status = DftiCreateDescriptorDM(MPI_COMM_WORLD ,desc,DFTI_DOUBLE,        &
         &                                DFTI_COMPLEX,rank,dims_global) 
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Create parallel descriptor failed')
    END IF
    
    ! Obtain some values of configuration parameters by calls to
    !  DftiGetValueDM
    
    status = DftiGetValueDM(desc,CDFT_LOCAL_SIZE,size)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Get CDFT local size failed')
    END IF
    
    status = DftiGetValueDM(desc,CDFT_LOCAL_NX,NX)
    IF (STATUS .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Get CDFT local maxpts failed')
    ENDIF
    
    STATUS = DftiGetValueDM(desc,CDFT_LOCAL_X_START,START_X)
    IF (STATUS .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Get CDFT local x start failed')      
    END IF
    
    status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_NX,NX_OUT)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Get CDFT local out failed')
    END IF
    
    status = DftiGetValueDM(desc,CDFT_LOCAL_OUT_X_START,START_X_OUT)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Get CDFT local out x start failed')
    ENDIF
    

    write(*,*) 'size: ',size
    write(*,*) 'Nx: ',NX
    write(*,*) 'Start_x: ',start_x
    write(*,*) 'Nx_out: ',nx_out
    write(*,*) 'start_x_out: ',start_x_out

    ! Allocate auxiliary array
    ALLOCATE(work(size))
    
    
    ! Specify a value(s) of configuration parameters by a call(s) to
    !  DftiSetValueDM
    
    status = DftiSetValueDM(desc,CDFT_WORKSPACE,work)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Set CDFT workspace failed')
    END IF
    
    !     Perform initialization that facilitates DFT computation by a call to
    !        DftiCommitDescriptorDM
    
    status = DftiCommitDescriptorDM(desc)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Commit descriptor failed')
    END IF
    
    ! Create arrays for local parts of input and output data
    !  (if it is needed) and fill the local part of input data with
    !  values (for more information, see Distributing Data among Processes)
    
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
    
    
    ! Compute the transform by calling
    !  DftiComputeForwardDM or DftiComputeBackward
    !  (Compute Forward transform)
    
    STATUS = DftiComputeForwardDM(desc,in)
    IF (STATUS .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Compute forward cdft failed')
    END IF
    
    
    ! Release memory allocated for a descriptor by a call to
    !  DftiFreeDescriptorDM (Free DftiDM descriptor)
    
    status = DftiFreeDescriptorDM(desc)
    IF (status .NE. DFTI_NO_ERROR) THEN
       CALL tddm_stop('Free descriptor failed')
    END IF
    
    
    DEALLOCATE(work)
    
  END SUBROUTINE FourierTransform_parallel
  
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
