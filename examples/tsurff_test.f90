PROGRAM tsurff_test
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax, mmin, mmax
  INTEGER                    :: lmax_total
  INTEGER                    :: numkpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numphipts 
  REAL(dp)                   :: k_cutoff, deltak
  COMPLEX(dp), ALLOCATABLE   :: b(:, :, :)
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
  
  k_cutoff = 2.0_dp
  
  radius_boundary = 50.0_dp
  lmax = 10

  
  filename = 'data/h2p/surfaces/sphfunc.rb50.000.lmax010'
  
  CALL initialize_tsurff(filename, radius_boundary, lmax, &
       deltak, k_cutoff, numkpts, numthetapts, numphipts, &
       lmax_total, mmax )
  
  mmin = - mmax
  
  ALLOCATE(b(1:numkpts,1:numthetapts,1:numphipts))
  
  CALL get_flux(filename, lmax, mmax, b )
  
  CALL write_polar_amplitude(b,'results/probkrktheta',&
       'probkrtheta')
  
  !CALL write_momentum_amplitude(b,'results/probk')
  
  WRITE(*,*) 'Our spectra is ready!!!'
  DEALLOCATE(b)
  
END PROGRAM tsurff_test
