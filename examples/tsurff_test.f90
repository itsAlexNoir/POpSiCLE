PROGRAM tsurff_test
  
  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax, mmin, mmax
  INTEGER                    :: lmax_total
  INTEGER                    :: numkpts
  REAL(dp)                   :: k_cutoff, dk
  REAL(dp)                   :: coulomb_exp
  COMPLEX(dp), ALLOCATABLE   :: b(:, :, :)
  CHARACTER(LEN=100)         :: filename
  LOGICAL                    :: desired_gauge
  
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
  
  
  k_cutoff = 2.5_dp
  dk = 0.005_dp
  radius_boundary = 50.0_dp
  lmax = 10
  coulomb_exp = 0.5_dp
  desired_gauge = .FALSE.
  
  filename = 'data/h2p/surfaces/sphfunc.rb50.000.lmax010'
  
  CALL initialize_tsurff(filename, radius_boundary, lmax, dk, &
       k_cutoff, numkpts, numthetapts, numphipts, &
       lmax_total ,desired_gauge, coulomb_exp )
  
  ALLOCATE(b(1:numkpts,1:numthetapts,1:numphipts))
  
  CALL get_flux(filename, lmax, b )
  
  CALL write_polar_amplitude(b,'results/probkrktheta',&
       'probkrtheta')

  CALL write_mes(b,'results/mes')

  CALL write_pes(b,'results/pes')
  
  WRITE(*,*) 'Our spectra is ready!!!'
  DEALLOCATE(b)
  
END PROGRAM tsurff_test
