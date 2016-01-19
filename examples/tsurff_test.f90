PROGRAM tsurff_test

  USE popsicle
  
  REAL(dp)                   :: radius_boundary
  INTEGER                    :: lmax_desired
  REAL(dp)                   :: maximumk, deltak
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
  
  filename = './data/h2p/surfaces/sphfunc.rb40.000.lmax010'
  
  CALL initialize_tsurff(filename, radius_boundary, lmax_desired, &
       deltak, maximumk )

  
  
END PROGRAM tsurff_test
