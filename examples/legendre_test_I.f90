PROGRAM legendre_test

  USE popsicle_aux

  IMPLICIT NONE

  INTEGER                  :: maxfactlog, maxrpts
  REAL(dp), ALLOCATABLE    :: theta(:), costheta(:)
  REAL(dp), ALLOCATABLE    :: weights(:), phi(:)
  REAL(dp), ALLOCATABLE    :: r_ax(:)
  COMPLEX(dp), ALLOCATABLE :: func(:, :, :), funcp(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: func_lm(:, :, :), funcx_lm(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: funcy_lm(:, :, :), funcz_lm(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: funcp_lm(:, :, :)
  COMPLEX(dp), ALLOCATABLE :: sph_harmonics(:, :, :, :)
  REAL(dp)                 :: rsum, myerror, dr
  COMPLEX(dp)              :: csum
  INTEGER                  :: itheta, iphi, ir
  INTEGER                  :: il, im, ill, imm
  
  rmesh%lmax = 6
  maxfactlog = 3 * rmesh%lmax + 1
  rmesh%numthetapts = rmesh%lmax + 1
  rmesh%numphipts = 2 * rmesh%lmax + 1
  maxrpts  = 1!500
  
  ALLOCATE(rmesh%thetapts(1:rmesh%numthetapts))
  ALLOCATE(rmesh%sinthetapts(1:rmesh%numthetapts))
  ALLOCATE(rmesh%costhetapts(1:rmesh%numthetapts))
  ALLOCATE(rmesh%wtheta(1:rmesh%numthetapts))
  ALLOCATE(rmesh%phipts(1:rmesh%numphipts))
  ALLOCATE(rmesh%sinphipts(1:rmesh%numphipts))
  ALLOCATE(rmesh%cosphipts(1:rmesh%numphipts))
  ALLOCATE(rmesh%wphi(1:rmesh%numphipts))
  
  ALLOCATE(r_ax(1:maxrpts))

  ALLOCATE(sph_harmonics(1:rmesh%numphipts, 1:rmesh%numthetapts,          &
       -rmesh%lmax:rmesh%lmax,0:rmesh%lmax))

  ALLOCATE(func(1:maxrpts,1:rmesh%numthetapts,1:rmesh%numphipts))
  ALLOCATE(funcp(1:maxrpts,1:rmesh%numthetapts,1:rmesh%numphipts))
  
  
  ALLOCATE(func_lm(1:maxrpts,-rmesh%lmax:rmesh%lmax, &
       0:rmesh%lmax))
  ALLOCATE(funcx_lm(1:maxrpts,-rmesh%lmax:rmesh%lmax, &
       0:rmesh%lmax))
  ALLOCATE(funcy_lm(1:maxrpts,-rmesh%lmax:rmesh%lmax, &
       0:rmesh%lmax))
  ALLOCATE(funcz_lm(1:maxrpts,-rmesh%lmax:rmesh%lmax, &
       0:rmesh%lmax))
  ALLOCATE(funcp_lm(1:maxrpts,-rmesh%lmax:rmesh%lmax, &
       0:rmesh%lmax))
  
  dr = 0.05
  DO ir = 1, maxrpts
     r_ax(ir) = REAL(ir-1,dp) * dr
  ENDDO
  
  CALL make_gauss_legendre(-1.0_dp, 1.0_dp, rmesh%costhetapts, rmesh%wtheta)
  
  rmesh%thetapts    = ACOS(rmesh%costhetapts)
  rmesh%sinthetapts = SIN(rmesh%thetapts)
  rmesh%deltaphi    = twopi / REAL(rmesh%numphipts, dp)
  
  DO iphi = 1, rmesh%numphipts
     
     rmesh%phipts(iphi)    = REAL(iphi - 1, dp) * rmesh%deltaphi 
     rmesh%sinphipts(iphi) = SIN(rmesh%phipts(iphi))
     rmesh%cosphipts(iphi) = COS(rmesh%phipts(iphi))
     rmesh%wphi(iphi)      = rmesh%deltaphi
     
  ENDDO
  
  CALL create_spherical_harmonics( sph_harmonics, rmesh%lmax, &
       rmesh%costhetapts, rmesh%phipts )
  
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
  
  func  = ZERO
  funcp = ZERO

  DO iphi = 1, rmesh%numphipts
     DO itheta = 1, rmesh%numthetapts
        DO ir = 1, maxrpts
           
!!$           func(ir,itheta,iphi) = r_ax(ir) * EXP(-r_ax(ir)) * &
!!$                3.0_dp / 32.0_dp * SQRT(77.0_dp / pi ) * &
!!$                EXP(-ZIMAGONE * 5 * phi(iphi)) * &
!!$                sin(theta(itheta))**5

           ! l=5, m=-5
!!$           func(ir,itheta,iphi) = (3.0_dp / 32.0_dp) * SQRT(77.0_dp / pi ) * &
!!$                EXP(-ZIMAGONE * 5.0_dp * rmesh%phipts(iphi) ) * &
!!$                rmesh%sinthetapts(itheta)**5
!!$
!!$           funcp(ir,itheta,iphi) = - (3.0_dp / 32.0_dp) * SQRT(77.0_dp / pi ) * &
!!$                EXP(-ZIMAGONE * 5.0_dp * rmesh%phipts(iphi) ) * &
!!$                rmesh%sinthetapts(itheta)**5
           ! l = 3 , m = 0
           func(ir,itheta,iphi) = (1.0_dp / 4.0_dp) * SQRT(7.0_dp / pi ) * &
                (5.0_dp * rmesh%costhetapts(itheta)**3 - &
                3.0_dp * rmesh%costhetapts(itheta))
           
           funcp(ir,itheta,iphi) = - (1.0_dp / 4.0_dp) * SQRT(7.0_dp / pi ) * &
                (5.0_dp * rmesh%costhetapts(itheta)**3 - &
                3.0_dp * rmesh%costhetapts(itheta))
           
!!$           func(ir,itheta,iphi) = (1.0_dp / 4.0_dp) * SQRT(5.0_dp / pi) * &
!!$                (3.0_dp * costheta(itheta)**2 - 1.0_dp)
           
           ! l=0, m=0
!!$           func(ir,itheta,iphi) = 0.5_dp / SQRT(pi)
!!$           funcp(ir,itheta,iphi) = 0.0_dp
        ENDDO
     ENDDO
  ENDDO
  
  
  DO ir = 1, 1!numrpts
     CALL make_sht(func(ir,:,:), funcp(ir, :, :), sph_harmonics, rmesh%lmax, &
          func_lm(ir,:,:), funcx_lm(ir, :, :), funcy_lm(ir, :, :),  &
          funcz_lm(ir, :, :), funcp_lm(ir, :, :))
  ENDDO
  
  !! Print the outcome
  DO il = 0, rmesh%lmax
     DO im = -il,il
        write(*,*) 'il, im, coeff_lm'
        write(*,*) il, im, func_lm(1,im,il)
        write(*,*) il, im, funcx_lm(1,im,il)
        write(*,*) il, im, funcy_lm(1,im,il)
        write(*,*) il, im, funcz_lm(1,im,il)
        write(*,*) il, im, funcp_lm(1,im,il)
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
  
  WRITE(*,*) 'The end!'
  
END PROGRAM legendre_test
