PROGRAM tsurff_test
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax_desired, mmax
  REAL(dp)                   :: maximumk, deltak
  COMPLEX(dp), ALLOCATABLE   :: blm(:, :, :)
  CHARACTER(LEN=100)         :: filename
  
  !--------------------------------------!
  
  WRITE(*,*) '****************************'
  WRITE(*,*) '       T-SURFF test         '
  WRITE(*,*) '****************************'
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Welcome to the T-SURFF test. This little program'
  WRITE(*,*) ' will help to understand and use a surface file to'
  WRITE(*,*) ' extract the photoelectron spectra from it.'
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) ' Initializing...'
  WRITE(*,*)
  
  
  deltak = 0.01_dp
  
  maximumk = 2.0_dp
  
  radius_boundary = 40.0_dp
  lmax_desired = 10

  
  filename = './data/h2p/surfaces/sphfunc.rb40.000.lmax010'
  
  CALL initialize_tsurff(filename, radius_boundary, lmax_desired, &
       deltak, maximumk, numkpts, mmax )
  
  ALLOCATE(blm(1:numkpts,0:mmax,0:lmax_desired))
  
  CALL get_flux(filename, blm )
  
  DEALLOCATE(blm)
  
END PROGRAM tsurff_test
