PROGRAM legendre_test

  USE popsicle

  IMPLICIT NONE

  INTEGER                  :: numrpts
  INTEGER                  :: numthetapts, numphipts
  INTEGER                  :: lmax, maxfactlog, mmax
  REAL(dp), ALLOCATABLE    :: theta(:), costheta(:)
  REAL(dp), ALLOCATABLE    :: weights(:), phi(:)
  REAL(dp), ALLOCATABLE    :: r_ax(:), radialfunc(:)
  COMPLEX(dp), ALLOCATABLE :: trial(:, :), shcoeff(:, :)
  COMPLEX(dp), ALLOCATABLE :: func(:, :, :), shfunc(:, :, :)
  REAL(dp)                 :: rsum, error, deltaphi, dr
  COMPLEX(dp)              :: csum
  INTEGER                  :: itheta, iphi, ir
  INTEGER                  :: il, im, ill, imm
  COMPLEX(dp), ALLOCATABLE              :: psi(:)
  
  lmax = 6
  mmax = lmax
  maxfactlog = 3 * lmax + 1
  numthetapts = lmax + 1
  numphipts = 2 * lmax + 1
  numrpts  = 1!500
  
  
  ALLOCATE(theta(1:numthetapts),costheta(1:numthetapts))
  ALLOCATE(weights(1:numthetapts))
  ALLOCATE(phi(1:numphipts))
  ALLOCATE(r_ax(1:numrpts))
  
  ALLOCATE(trial(1:numthetapts,1:numphipts))
  ALLOCATE(shcoeff(-mmax:mmax,0:lmax))
  ALLOCATE(func(1:numrpts,1:numthetapts,1:numphipts))
  ALLOCATE(radialfunc(1:numrpts))
  ALLOCATE(shfunc(1:numrpts,-mmax:mmax,0:lmax))
  ALLOCATE(psi(1:numthetapts))
  
  ! Axis parade!  
  deltaphi = twopi / REAL(numphipts,dp)
  
  DO iphi = 1, numphipts
     phi(iphi) = REAL(iphi-1,dp) * deltaphi
  ENDDO

  dr = 0.05

  DO ir = 1, numrpts
     r_ax(ir) = REAL(ir-1,dp) * dr
  ENDDO

  CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta, weights)
  theta = acos(costheta)
  CALL initialize_legendre_stuff(lmax,costheta)
  CALL initialize_spherical_harmonics(lmax,mmax,costheta,phi)
  
  !******************************************************!

!!$  CALL get_gauss_stuff(1.0_dp, 3.0_dp, theta, weights)
!!$  
!!$  rsum = 0.0_dp
!!$  DO itheta = 1, numthetapts
!!$     rsum = rsum + 1.0_dp / theta(itheta) * weights(itheta)
!!$  ENDDO
!!$  
!!$  error = abs(rsum - LOG(3.0_dp))
!!$  write(*,*)
!!$  write(*,'(A50)') 'Function: 1/x'
!!$  write(*,'(A50)') 'Limits of the integral: [1,3]'
!!$  write(*,'(A50,I4)') 'Number of quadrature points: ',numthetapts
!!$  write(*,'(A50,F15.8)') 'Result:  ',rsum
!!$  write(*,'(A50,F15.8)') 'Error: ', error
!!$  write(*,*)
!!$  write(*,*)
  

  !******************************************************!
  
!!$  theta = 0.0_dp
!!$  costheta = 0.0_dp
!!$  weights = 0.0_dp
!!$  CALL get_gauss_stuff(-1.0_dp, 1.0_dp, costheta, weights)
!!$  theta = acos(costheta)
!!$
!!$  CALL initialize_spherical_harmonics(lmax,costheta)
!!$  
!!$  !write(*,*)
!!$  !write(*,*) 'our no / theta / costheta / weights...'
!!$  !do itheta = 1, numthetapts
!!$  !   write(*,'(I4,3(F15.8))') itheta, theta(itheta),costheta(itheta), weights(itheta)
!!$  !enddo
!!$  !write(*,*)
!!$  
!!$  func = ZERO
!!$  DO iphi = 1, numphipts  
!!$     DO itheta = 1, numthetapts
!!$        DO ir = 1, numrpts
!!$           func(ir,itheta,iphi) = r_ax(ir) * EXP(-r_ax(ir)) * sin(theta(itheta))**4 * &
!!$                costheta(itheta)**2 * sin(phi(iphi))**2
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$  
!!$
!!$  ! First, check the integrals
!!$  rsum = 0.0
!!$
!!$  rsum = SUM(sin(phi)**2) * dphi
!!$
!!$  WRITE(*,*)
!!$  WRITE(*,*) 'Int (sin^2 (phi)) = ',rsum
!!$
!!$  rsum = 0.0
!!$
!!$  rsum = SUM(sin(theta)**4 * costheta**2 * weights)
!!$
!!$  WRITE(*,*)
!!$  WRITE(*,*) 'Int (sin^4(theta) * cos^2(theta)) = ',rsum
!!$  
!!$  WRITE(*,*)
!!$
!!$  !WRITE(*,*) 'Now, try to decompose this wavefunction.'
!!$  !CALL make_sht(trial,lmax,mmax,weights,shcoeff)
!!$  
!!$  radialfunc = 0.0_dp
!!$  DO iphi = 1, numphipts
!!$     DO itheta = 1, numthetapts
!!$        DO ir = 1, numrpts
!!$           radialfunc(ir) = radialfunc(ir) + &
!!$                REAL(func(ir,itheta,iphi),dp) * weights(itheta)
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  radialfunc = radialfunc * dphi
!!$
!!$  open(unit=33,form='formatted',file='results/radialfunc.dat')
!!$
!!$  DO ir = 1, numrpts
!!$     write(33,*) radialfunc(ir)
!!$  ENDDO
!!$
!!$  CLOSE(33)

  !**********************************************************
  
!!$  WRITE(*,*) 'Try the orthonormality of spherical harmonics.'
!!$    
!!$  shcoeff = ZERO
!!$  rsum = 0.0_dp
!!$  csum = ZERO
!!$  ! Trial functions parade!!
!!$  !
!!$  trial = ZERO
!!$  ! l=0, m=0
!!$  !trial(:,1) = sqrt(1.0/4.0/pi)
!!$  !
!!$  ! l=1, m=0
!!$  !trial(:,1) = CMPLX(1.0_dp,1.0_dp, dp) * sqrt(3.0 / 4.0 / pi ) * cos(theta)
!!$  !
!!$  ! l=2, m=0
!!$  !trial(:,1) = sqrt(5.0 / 16.0  / pi ) * (3.0 * cos(theta)**2 - 1.0 )
!!$  ! l=2, m=1
!!$  !DO itheta = 1, numthetapts
!!$  !   DO iphi = 1, numphipts
!!$  !      trial(itheta,iphi) = - SQRT(15.0 / 8.0 / pi) * &
!!$  !           SIN(theta(itheta)) * costheta(itheta) * EXP(ZIMAGONE * phi(iphi))
!!$  !   ENDDO
!!$  !ENDDO
!!$  
!!$  il = 6
!!$  im = 6
!!$  DO itheta = 1, numthetapts
!!$     DO iphi = 1, numphipts
!!$        IF(im.LT.0) THEN
!!$           trial(itheta,iphi) = (-1.0_dp)**ABS(im) * normfact(ABS(im),il) * legenpl(itheta-1,ABS(im),il) * &
!!$                EXP(ZIMAGONE * im * phi(iphi))
!!$        ELSE
!!$           trial(itheta,iphi) = normfact(im,il) * legenpl(itheta-1,im,il) * &
!!$                EXP(ZIMAGONE * im * phi(iphi))
!!$        ENDIF
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  DO iphi = 1, numphipts
!!$     DO itheta = 1, numthetapts
!!$        csum = csum + CONJG(trial(itheta,iphi)) * trial(itheta,iphi) * weights(itheta)
!!$     ENDDO
!!$  ENDDO
!!$  
!!$  csum = csum * dphi
!!$  ! We add 2 pi factor, because we are not integrating in phi.
!!$  !sum = sum * twopi
!!$  
!!$  WRITE(*,*) 'The norm for SH with angular momentum ',il, 'is: ',csum 
!!$  WRITE(*,*)
!!$  
!!$  ! Now test the decomposition into spherical harmonics, A.K.A. SHT
!!$  ! (Spherical Harmonics tranform)
!!$  
!!$  WRITE(*,*) ' Now test the decomposition into spherical harmonics,'
!!$  WRITE(*,*) ' A.K.A. SHT (Spherical Harmonics tranform).'
!!$  
!!$  CALL make_sht(trial,lmax,mmax,weights,shcoeff)
!!$
!!$  ! Print the outcome
!!$  DO il = 0, lmax
!!$     DO im = -il,il
!!$        write(*,*) 'il, im, shcoeff'
!!$        write(*,*) il, im, shcoeff(im,il)
!!$     ENDDO
!!$  ENDDO

  !*************************************************************!

  ! l=5, m=-5
  func = ZERO
  DO iphi = 1, numphipts
     DO itheta = 1, numthetapts
        DO ir = 1, numrpts
           
!!$           func(ir,itheta,iphi) = r_ax(ir) * EXP(-r_ax(ir)) * &
!!$                3.0_dp / 32.0_dp * SQRT(77.0_dp / pi ) * &
!!$                EXP(-ZIMAGONE * 5 * phi(iphi)) * &
!!$                sin(theta(itheta))**5
           func(ir,itheta,iphi) = ZIMAGONE * (3.0_dp / 32.0_dp) * SQRT(77.0_dp / pi ) * &
                EXP(-ZIMAGONE * 5.0_dp * phi(iphi)) * &
                sin(theta(itheta))**5
!!$           func(ir,itheta,iphi) = (1.0_dp / 4.0_dp) * SQRT(5.0_dp / pi) * &
!!$                (3.0_dp * costheta(itheta)**2 - 1.0_dp)
           
!!$           ! l=0, m=0
!!$           func(ir,itheta,iphi) = 0.5_dp / SQRT(pi)

           !func(ir,itheta,iphi) = exp(ZIMAGONE * theta(itheta))
        ENDDO
     ENDDO
  ENDDO
  
  shfunc = ZERO
  
  DO ir = 1, 1!numrpts
     CALL make_sht(func(ir,:,:),lmax,mmax,weights,shfunc(ir,:,:))
  ENDDO
  
!!$  DO il = 0, lmax
!!$     DO im = -mmax, mmax
!!$        DO ir = 1, numrpts
!!$           shfunc(ir,im,il) = shfunc(ir,im,il) * r_ax(ir)
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
      
  !! Print the outcome
  DO il = 0, lmax
     DO im = -il,il
        write(*,*) 'il, im, shcoeff'
        write(*,*) il, im, shfunc(1,im,il)
     ENDDO
  ENDDO
  
  
!!$  OPEN(unit=33,form='formatted',file='results/angfunc.dat')
!!$
!!$  psi = ZERO
!!$  do iphi = 1, numphipts
!!$     do itheta = 1, numthetapts
!!$        csum = ZERO
!!$        DO il = 0, lmax
!!$           DO im = -mmax, mmax
!!$              DO ir = 1, numrpts
!!$                 csum = csum + shfunc(ir,im,il) * &
!!$                      sph_harmonics(iphi,itheta,im,il)
!!$              ENDDO
!!$           ENDDO
!!$        ENDDO
!!$        psi(itheta) = csum
!!$     ENDDO
!!$  enddo
!!$
!!$  DO itheta = 1, numthetapts
!!$     write(33,*) REAL(psi(itheta), dp),AIMAG(psi(itheta))
!!$  ENDDO
!!$     
!!$  DO il = 0, lmax
!!$     DO im = -mmax, mmax
!!$        DO ir = 1, numrpts
!!$           WRITE(33,*) rpts(ir), il, im, REAL(shfunc(ir,im,il),dp), &
!!$                AIMAG(shfunc(ir,im,il))
!!$        ENDDO
!!$     ENDDO
!!$  ENDDO
!!$
!!$  CLOSE(33)
  
  !********************************************!
  
  DEALLOCATE(theta,weights)
  DEALLOCATE(r_ax,phi)

  DEALLOCATE(func,radialfunc)
  DEALLOCATE(trial,shcoeff,shfunc)


  WRITE(*,*) 'The end!'
  
END PROGRAM legendre_test
