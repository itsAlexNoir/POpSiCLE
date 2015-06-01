PROGRAM fourier3D_test
  
  USE fdf
  USE omp_lib
  USE constants
  USE fdmesh_singleproc
  USE io_aux
  USE momprob_tools

  
  IMPLICIT NONE
  
  
  !--Program variables-----------------------------------------------------!
  
  INTEGER                    :: maxxpts, maxypts, maxzpts, maxrpts
  INTEGER                    :: numxpts, numypts, numzpts, numrpts
  INTEGER                    :: maxkxpts, maxkypts, maxkzpts, maxkrpts
  INTEGER                    :: numkxpts, numkypts, numkzpts, numkrpts
  INTEGER                    :: maxkelecpts, maxktotpts
  INTEGER                    :: numkelecpts, numktotpts
  INTEGER                    :: numproc1drx
    INTEGER                    :: numproc1dry
  INTEGER                    :: numproc1dz
    INTEGER                    :: numproc1dr
  INTEGER                    :: maxproc1drx      
  INTEGER                    :: maxproc1dry       
  INTEGER                    :: maxproc1dz         
  INTEGER                    :: maxproc1dr
  
  INTEGER                    :: iproc, ierror       
  
  INTEGER                    :: ix, iy, iz, ir
  INTEGER                    :: ikx, iky, ikz, ikr 
  INTEGER                    :: ixg, iyg, izg, irg
  INTEGER                    :: iprocx, iprocy, iprocz, iprocr
  INTEGER                    :: ipro
  INTEGER                    :: mintime, maxtime, itimestep
  COMPLEX(dp), ALLOCATABLE   :: psipro(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE   :: psi(:, :, :, :)
  COMPLEX(dp), ALLOCATABLE   :: psik(:, :, :, :)
  REAL(dp), ALLOCATABLE      :: jacobian(:)
  REAL(dp), ALLOCATABLE      :: mask(:, :, :, :)
  REAL(dp)                   :: M1, M2
  
  ! DFT arrays 
  COMPLEX(dp), ALLOCATABLE   :: inrz(:, :)
  COMPLEX(dp), ALLOCATABLE   :: inrho(:)
  COMPLEX(dp), ALLOCATABLE   :: outrho(:)
  INTEGER                    :: num_threads
  CHARACTER(LEN = 200)       :: filename
  CHARACTER(LEN = 6)         :: ctime
  CHARACTER(LEN = 4)         :: cprocessor

  REAL(dp)                   :: rad_inner
  REAL(dp)                   :: rad_outer

  LOGICAL                    :: output_rhoz
  LOGICAL                    :: output_rz
  LOGICAL                    :: output_krkrhokz
  LOGICAL                    :: output_krhokz
  LOGICAL                    :: output_krkz
  LOGICAL                    :: output_knke
  LOGICAL                    :: output_ktot
  LOGICAL                    :: mask_on
  LOGICAL                    :: save_psik_real
  LOGICAL                    :: save_psik_complex
  CHARACTER(LEN=100)         :: data_directory

  !------------------------------------------------------------------------!
  
  WRITE(*,*) '********************************'
  WRITE(*,*) ' Fourier3D test.'
  WRITE(*,*) '********************************'
  WRITE(*,*) 
  
  WRITE(*,*) 'Opening...'
  
 
  CALL fdf_init( 'momprob.fdf', 'momprob.dat' )

  
  mintime	   = fdf_integer( 'MinTime', 0 )
  maxtime	   = fdf_integer( 'MaxTime', 0 )
  itimestep	   = fdf_integer( 'TimeStep', 0 )
  maxrpts	   = fdf_integer( 'MaxptsR', 51 )
  maxrhopts        = fdf_integer( 'MaxptsRHO', 51 )
  maxzpts	   = fdf_integer( 'MaxptsZ', 51 )
  
  maxkrpts	   = fdf_integer( 'MaxptsKR', 51 )
  maxkrhopts       = fdf_integer( 'MaxptsKRHO', 51 )
  maxkzpts	   = fdf_integer( 'MaxptsKZ', 51 )
  maxkelecpts      = fdf_integer( 'MaxptsKE', 51)
  maxktotpts       = fdf_integer( 'MaxptsKTOT', 51)
  
  numproc1dr       = fdf_integer( 'NumProcessorsR', 1 )
  numproc1drho     = fdf_integer( 'NumProcessorsRHO', 1 )
  numproc1dz       = fdf_integer( 'NumProcessorsZ', 1 )
  
  maxproc1dr   = numproc1dr - 1
  maxproc1drho = numproc1drho - 1
  maxproc1dz   = numproc1dz - 1 
  
  numrpts	   = numproc1dr   * maxrpts
  numrhopts        = numproc1drho * maxrhopts
  numzpts	   = numproc1dz   * maxzpts
  
  numkrpts	   = numproc1dr   * maxkrpts
  numkrhopts       = numproc1drho * maxkrhopts
  numkzpts	   = numproc1dz   * maxkzpts
  
  numkelecpts	   = numproc1dz   * maxkelecpts
  numktotpts       = numproc1dz   * maxktotpts
  
  mask_on         = fdf_boolean( 'MaskOn', .FALSE.)
  
  M1     = fdf_double( 'M1', 1836.0_dp )
  M2     = fdf_double( 'M2', 1836.0_dp )
  
  rad_inner       = fdf_double( 'InnerRadius', 15.0_dp )
  rad_outer       = fdf_double( 'OuterRadius', 50.0_dp )
  
  output_rhoz     = fdf_boolean( 'GraphicsRhoZ', .FALSE. )
  output_rz       = fdf_boolean( 'GraphicsRZ', .FALSE. )
  output_krkrhokz = fdf_boolean( 'GraphicsKRKRhoKZ', .FALSE. )
  output_krhokz   = fdf_boolean( 'GraphicsKRhoKZ', .FALSE. )
  output_krkz     = fdf_boolean( 'GraphicsKRKZ', .FALSE. )
  output_knke     = fdf_boolean( 'GraphicsKNKE', .FALSE. )
  output_ktot     = fdf_boolean( 'GraphicsKTOT', .FALSE. )
  save_psik_real  = fdf_boolean( 'SavePsikReal', .TRUE. )
  save_psik_complex = fdf_boolean( 'SavePsikComplex', .FALSE. )


  data_directory   = fdf_string('DataDirectory', './000/psi/')
  
  WRITE(ctime, '(I6.6)') maxtime
  
  
  ! Allocate arrays
  WRITE(*,*) 'Allocating arrays...'
  ALLOCATE(psipro(1:maxrpts, 1:maxrhopts, 1:maxzpts))
  ALLOCATE(psi(1:numrpts, 1:numrhopts, 1:numzpts))
  ALLOCATE(psik(1:numkrpts, 1:numkrhopts, 1:numkzpts))
  ALLOCATE(jacobian(1:numrhopts))
  
  ! Initialize and make meshes
  WRITE(*,*) 'Making meshes...'
  
  CALL initialize_mesh_arrays( )
  CALL initialize_momentum_mesh_arrays( )
  CALL make_meshpoints( )
  CALL make_kmeshpoints( )
  
  ! The jacobian for the Hankel transform
  jacobian = SQRT(gpts * gp)

  !-----------------------------------------------!
  ! Create the mask for eliminating bound states  !
  !-----------------------------------------------!

  IF (mask_on) THEN
     ALLOCATE(mask(1:numrpts,1:numrhopts,1:numzpts))
     mask = 0.0_dp
     
     CALL create_mask( mask, numrpts, numrhopts, numzpts, &
          rad_inner, rad_outer )
  ENDIF
  
  !-----------------------------------------------!
  ! Read data from files, across all proccessors. !
  !-----------------------------------------------!
  
  WRITE(*,*) 'Read data...'
  
  ipro = -1
  
  DO iprocz = 0, maxproc1dz     
     DO iprocrho = 0, maxproc1drho    
        DO iprocr = 0, maxproc1dr    
           
           ipro = ipro + 1
           
           WRITE(cprocessor, '(I4.4)') ipro
           
           filename = TRIM(data_directory) // cprocessor //             &
                '/psi.' // cprocessor // '.' // ctime // '.dat'
           
           OPEN(UNIT = 10, FORM = 'unformatted', FILE = filename)
           
           READ(10) psipro
           
           CLOSE(UNIT = 10)
           
           DO iz = 1, maxzpts
              DO irho = 1, maxrhopts
                 DO ir = 1, maxrpts 
                    
                    irg = iprocr * maxrpts + ir
                    irhog = iprocrho * maxrhopts + irho
                    izg = iprocz * maxzpts + iz
                    
                    psi(irg, irhog, izg) = psi(irg, irhog, izg) + psipro(ir, irho, iz)
                    
                 ENDDO
              ENDDO
           ENDDO
           
        ENDDO
     ENDDO
  ENDDO
  
  !----------------------------------------------------!
  ! Apply the mask in order to eliminate bound states. !
  !----------------------------------------------------!
  
  IF (mask_on) psi = psi * mask
  
  !-------------------------------------------------!
  ! Transform the wavefunction into momentum space. !
  !-------------------------------------------------!
  
  WRITE(*,*) 'Transforming to momentum space...'
   
  !---------------------!
  ! R and Z coordinates !
  !---------------------!
  
  WRITE(*,*) '...R and Z coordinates...'
  
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ir,irho,iz,ikrho,inrz) COPYIN(inrz)
  WRITE(*,*) 'Opening thread ',OMP_GET_THREAD_NUM(),' for FFT...'
  ALLOCATE(inrz(1:numrpts,numzpts))
  
  
  !$OMP DO SCHEDULE(STATIC)
  DO irho = 1, numrhopts
     !write(*,*) 'Thread num',OMP_GET_THREAD_NUM()
     
     DO iz = 1, numzpts
        DO ir = 1, numrpts
           inrz(ir,iz) = psi(ir,irho,iz) * SQRT(fp(ir) * hp(iz))
        ENDDO
     ENDDO
     
     ! Compute forward transform
     CALL FourierTransform2D(inrz, numrpts, numzpts)
     
     DO iz = 1, numzpts
        DO ir = 1, numkrpts
           psik(ir,irho,iz) = inrz(ir,iz)
        ENDDO
     ENDDO
     
  ENDDO
  !$OMP END DO NOWAIT
  
  DEALLOCATE(inrz)
  
  !$OMP END PARALLEL
  
  
  !-----------------!
  ! RHO coordinate  !
  !-----------------!
  IF (numrhopts.NE.1) THEN
     WRITE(*,*) '...RHO coordinate...'
     
     !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ir,iz,irho,ikrho,inrho,outrho)  COPYIN(inrho,outrho)
     ALLOCATE(inrho(1:numrhopts))
     ALLOCATE(outrho(1:numrhopts))
     
     !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
     DO iz = 1, numzpts
        DO ir = 1, numrpts
           
           inrho(:) = psik(ir,:,iz)
           
           DO ikrho = 1, numkrhopts
              outrho(ikrho) = HankelTransform(inrho, numrhopts, gpts, kgpts(ikrho), jacobian)
           ENDDO
           
           psik(ir,:,iz) = outrho(:)
           
        ENDDO
     ENDDO
     
     !$OMP END DO NOWAIT
     
     DEALLOCATE(inrho)
     DEALLOCATE(outrho)
     
     !$OMP END PARALLEL
  ENDIF
  

  ! Multiply by momentum space volume element
  psik = psik * (oneovertwopi**3) * &
       mesh%deltar * mesh%deltaz * mesh%deltarho
  
  
  !-----------------!
  ! Save to a file
  !-----------------!
  WRITE(*,*) 'Save to a data file...'
  
  !---------------------------------!
  
  IF (output_rhoz) THEN
     WRITE(*,*) 'Saving psi 2D (rho vs. z) to a file...'
     filename = 'probrhoz'
     CALL write_psi2D(psi, numrpts, numrhopts, numzpts, filename, 'rhoz')
  ENDIF
  
  IF (output_rz) THEN
     WRITE(*,*) 'Saving psi 2D (r vs. z) to a file...'
     filename = 'probrz'
     CALL write_psi2D(psi, numrpts, numrhopts, numzpts, filename, 'rz')
  ENDIF
  
  !--------------------------------!

  
  !---------------------------------!
  
  IF (output_krkrhokz) THEN
     WRITE(*,*) 'Saving 3D psik to a file...'
     filename = 'probkrkrhokz'
     CALL write_psik(psik, numkrpts, numkrhopts, numkzpts, filename, &
          save_psik_real, save_psik_complex)
  ENDIF
  
  !---------------------------------!
  
  IF (output_krhokz) THEN
     WRITE(*,*) 'Saving psik 2D (krho vs. kz) to a file...'
     filename = 'probkrhokz'
     CALL write_psik2D(psik, numkrpts, numkrhopts, numkzpts, filename, 'krhokz')
  ENDIF
  
  IF (output_krkz) THEN
     WRITE(*,*) 'Saving psik 2D (kr vs. kz) to a file...'
     filename = 'probkrkz'
     CALL write_psik2D(psik, numkrpts, numkrhopts, numkzpts, filename, 'krkz')
  ENDIF
  
  !--------------------------------!
  
  IF (output_knke) THEN
     WRITE(*,*) 'Saving psik 2D (kn vs. ke) to a file...'
     filename = 'probknke'
     CALL write_psik2D(psik, numkrpts, numkrhopts, numkzpts, filename, 'knke', numkelecpts)
  ENDIF

  !---------------------------------!
  
  IF (output_ktot) THEN
     WRITE(*,*) 'Saving psik 1D (k total) to a file...'
     filename = 'probktotal'
     CALL write_ktotal(psik, numkrpts, numkrhopts, numkzpts, &
          numkelecpts, numktotpts, M1, M2, filename)
  ENDIF
  
  !---------------------------------!
  
  !---------------------------------!
  
  
  !-------------------------------!
  ! Free memory and descriptors
  !-------------------------------!
  
  WRITE(*,*) 'Free memory...'
  
  DEALLOCATE(psipro)
  DEALLOCATE(psi)
  DEALLOCATE(psik)
  DEALLOCATE(jacobian)

  
  write(*,*) 'Fin!!'
  
END PROGRAM	fourier3D_test
