PROGRAM legendre_test

  USE popsicle

  IMPLICIT NONE

  INTEGER                  :: numrpts
  INTEGER                  :: numthetapts, numphipts
  INTEGER                  :: lmax, maxfactlog, mmax
  REAL(dp), ALLOCATABLE    :: theta(:), costheta(:)
  REAL(dp), ALLOCATABLE    :: weights(:), phi(:)
  REAL(dp), ALLOCATABLE    :: legendre(:,:,:), nfactor(:,:)
  REAL(dp), ALLOCATABLE    :: flog(:)
  COMPLEX(dp), ALLOCATABLE :: trial(:, :), shcoeff(:, :)
  REAL(dp)                 :: sum, error, dphi
  COMPLEX(dp)              :: csum
  INTEGER                  :: itheta, iphi, ir
  INTEGER                  :: il, im, ill, imm
  
  lmax = 7
  mmax = lmax
  maxfactlog = 3 * lmax + 1
  numthetapts = lmax + 1
  numphipts = 2 * lmax + 1
  
  
  ALLOCATE(theta(1:numthetapts),costheta(1:numthetapts))
  ALLOCATE(weights(1:numthetapts))
  ALLOCATE(phi(1:numphipts))
  ALLOCATE(legendre(0:numthetapts-1,0:lmax,0:lmax))
  ALLOCATE(nfactor(0:lmax,0:lmax))
  ALLOCATE(flog(0:maxfactlog))
  ALLOCATE(trial(1:numthetapts,1:numphipts))
  ALLOCATE(shcoeff(-mmax:mmax,0:lmax))
  
  ! Axis parade!
  CALL get_gauss_stuff(1.0_dp, 3.0_dp, theta, weights)
  
  dphi = twopi / REAL(numphipts,dp)
  
  DO iphi = 1, numphipts
     phi(iphi) = REAL(iphi-1,dp) * dphi
  ENDDO

  
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
  
!!$  write(*,*)
!!$  write(*,*) 'our no / theta / costheta / weights...'
!!$  do itheta = 1, numthetapts
!!$     write(*,'(I4,3(F15.8))') itheta, theta(itheta),costheta(itheta), weights(itheta)
!!$  enddo
!!$  write(*,*)


  CALL make_legendre(legendre, costheta, lmax, numthetapts-1)
  CALL make_factlog(flog, maxfactlog)
  CALL make_rnormal(nfactor, flog, lmax, maxfactlog)

  WRITE(*,*) 'Try the orthonormality of spherical harmonics.'
  
  CALL initialize_spherical_harmonics(lmax,costheta)
  
  shcoeff = ZERO
  sum = 0.0_dp
  csum = ZERO
    ! Trial functions parade!!
  !
  trial = ZERO
  ! l=0, m=0
  !trial(:,1) = sqrt(1.0/4.0/pi)
  !
  ! l=1, m=0
  !trial(:,1) = CMPLX(1.0_dp,1.0_dp) * sqrt(3.0 / 4.0 / pi ) * cos(theta)
  !
  ! l=2, m=0
  !trial(:,1) = sqrt(5.0 / 16.0  / pi ) * (3.0 * cos(theta)**2 - 1.0 )
  ! l=2, m=1
  !DO itheta = 1, numthetapts
  !   DO iphi = 1, numphipts
  !      trial(itheta,iphi) = - SQRT(15.0 / 8.0 / pi) * &
  !           SIN(theta(itheta)) * costheta(itheta) * EXP(ZIMAGONE * phi(iphi))
  !   ENDDO
  !ENDDO
  
  il = 7
  im = -7
  DO itheta = 1, numthetapts
     DO iphi = 1, numphipts
        IF(im.LT.0) THEN
           trial(itheta,iphi) = (-1.0_dp)**ABS(im) * nfactor(ABS(im),il) * legendre(itheta-1,ABS(im),il) * &
                EXP(ZIMAGONE * im * phi(iphi))
        ELSE
           trial(itheta,iphi) = nfactor(im,il) * legendre(itheta-1,im,il) * &
                EXP(ZIMAGONE * im * phi(iphi))
        ENDIF
     ENDDO
  ENDDO
  
  DO iphi = 1, numphipts
     DO itheta = 1, numthetapts
        csum = csum + CONJG(trial(itheta,iphi)) * trial(itheta,iphi) * weights(itheta)
     ENDDO
  ENDDO
  
  csum = csum * dphi
  ! We add 2 pi factor, because we are not integrating in phi.
  !sum = sum * twopi
  
  WRITE(*,*) 'The norm for SH with angular momentum ',il, 'is: ',csum 
  WRITE(*,*)
  
  ! Now test the decomposition into spherical harmonics, A.K.A. SHT
  ! (Spherical Harmonics tranform)
  
  WRITE(*,*) ' Now test the decomposition into spherical harmonics,'
  WRITE(*,*) ' A.K.A. SHT (Spherical Harmonics tranform).'
  
  CALL make_sht(trial,lmax,mmax,weights,shcoeff)

  ! Print the outcome
  DO il = 0, lmax
     DO im = -il,il
        write(*,*) 'il, im, shcoeff'
        write(*,*) il, im, shcoeff(im,il)
     ENDDO
  ENDDO
  
  DEALLOCATE(theta,weights)
  DEALLOCATE(legendre,nfactor,flog)
  
  WRITE(*,*) 'The end!'
  
END PROGRAM legendre_test
