PROGRAM legendre_test

  USE popsicle

  IMPLICIT NONE

  INTEGER                :: numrpts
  INTEGER                :: numthetapts
  INTEGER                :: lmax, maxfactlog
  REAL(dp), ALLOCATABLE  :: theta(:), costheta(:)
  REAL(dp), ALLOCATABLE  :: weights(:)
  REAL(dp), ALLOCATABLE  :: legendre(:,:,:), nfactor(:,:)
  REAL(dp), ALLOCATABLE  :: flog(:)
  REAL(dp), ALLOCATABLE  :: trial(:), shcoeff(:)
  REAL(dp)               :: sum, error
  INTEGER                :: itheta, ir
  INTEGER                :: il, im, ill, imm
  
  numthetapts = 21
  lmax = 10
  maxfactlog = 3 * lmax + 1
  
  ALLOCATE(theta(1:numthetapts),costheta(1:numthetapts))
  ALLOCATE(weights(1:numthetapts))
  ALLOCATE(legendre(0:numthetapts-1,0:lmax,0:lmax))
  ALLOCATE(nfactor(0:lmax,0:lmax))
  ALLOCATE(flog(0:maxfactlog))
  ALLOCATE(trial(1:numthetapts))
  ALLOCATE(shcoeff(0:lmax))
  
  CALL get_gauss_stuff(1.0_dp, 3.0_dp, theta, weights)
  
  sum = 0.0_dp
  DO itheta = 1, numthetapts
     sum = sum + 1.0_dp / theta(itheta) * weights(itheta)
  ENDDO
  
  error = abs(sum - LOG(3.0_dp))
  write(*,*)
  write(*,'(A50)') 'Function: 1/x'
  write(*,'(A50)') 'Limits of the integral: [1,3]'
  write(*,'(A50,I4)') 'Number of quadrature points: ',numthetapts
  write(*,'(A50,F15.8)') 'Result:  ',sum
  write(*,'(A50,F15.8)') 'Error: ', error
  write(*,*)
  write(*,*)

  
  theta = 0.0_dp
  costheta = 0.0_dp
  weights = 0.0_dp
  CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta, weights)
  !CALL get_gauss_stuff(0.0_dp, pi, theta, weights)
  theta = acos(costheta)
  
  write(*,*)
  write(*,*) 'our no / theta / costheta / weights...'
  do itheta = 1, numthetapts
     write(*,'(I4,3(F15.8))') itheta, theta(itheta),costheta(itheta), weights(itheta)
  enddo
  write(*,*)
  
  CALL make_legendre(legendre, costheta, lmax, numthetapts-1)
  CALL make_factlog(flog, maxfactlog)
  CALL make_rnormal(nfactor, flog, lmax, numthetapts-1)
  
!!$  im = 0
!!$  imm = 0
!!$  DO il = 0, lmax
!!$     ! DO im = -il, il
!!$     DO ill = 0, lmax
!!$        !       DO imm = -ill, ill        
!!$        sum = 0.0_dp
!!$        DO itheta = 1, numthetapts
!!$           sum = sum + legendre(itheta-1,im,il) * &
!!$                legendre(itheta-1,imm,ill) * &
!!$                weights(itheta)
!!$        ENDDO
!!$        WRITE(*,*) 'il, im, ill, imm, sum'
!!$        WRITE(*,*) il, im, ill, imm, sum
!!$        !        ENDDO
!!$        !      ENDDO
!!$     ENDDO
!!$     write(*,*)
!!$  ENDDO

  
!!$  theta = 0.0_dp
!!$  weights = 0.0_dp
!!$  legendre = 0.0_dp
!!$  CALL get_gauss_stuff(-1.0_dp, 1.0_dp, theta, weights)
!!$  CALL make_legendre(legendre, theta, lmax, numthetapts-1)
  
  !trial = sqrt(1.0/4.0/pi)
  trial = sqrt(3.0 / 4.0 / pi ) * cos(theta)
  !trial = sqrt(5.0 / 16.0  / pi ) * (3.0 * cos(theta)**2 - 1.0 )
  shcoeff = 0.0_dp
  
  DO il = 0, lmax
     DO itheta = 1, numthetapts
        shcoeff(il) = shcoeff(il) + trial(itheta) * legendre(itheta-1,0,il) * &
             nfactor(0,il) * weights(itheta)
     ENDDO
  ENDDO
  
  DO il = 0, lmax
     write(*,*) 'il, shcoeff'
     write(*,*) il, shcoeff(il)
  enddo
  
  deallocate(theta,weights)
  deallocate(legendre,nfactor,flog)
  
  WRITE(*,*) 'The end!'
  
END PROGRAM legendre_test
