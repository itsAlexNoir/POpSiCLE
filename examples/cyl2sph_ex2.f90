PROGRAM cyl2sph_ex2

  USE popsicle
  
  IMPLICIT NONE

  INTEGER                    :: numrhopts
  INTEGER                    :: numzpts
  INTEGER                    :: numrpts
  INTEGER                    :: numthetapts
  INTEGER                    :: numpts
  REAL(dp), ALLOCATABLE      :: rho_ax(:)
  REAL(dp), ALLOCATABLE      :: z_ax(:)
  INTEGER                    :: dims(2)
  
  COMPLEX(dp), ALLOCATABLE   :: Y_lm(:, :)
  COMPLEX(dp), ALLOCATABLE   :: R_nl(:, :)
  COMPLEX(dp), ALLOCATABLE   :: cylfunc(:, :)
  REAL(dp)                   :: drho, dz
  COMPLEX(dp), ALLOCATABLE   :: sphfunc(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dr(:, :)
  COMPLEX(dp), ALLOCATABLE   :: sphfunc_dth(:, :)
  COMPLEX(dp), ALLOCATABLE   :: radfunc(:, :)
  COMPLEX(dp)                :: ref_value
  REAL(dp)                   :: dr
  INTEGER                    :: lmax, il
  REAL(dp)                   :: Rboundary, tolerance
  REAL(dp)                   :: thetapt, rpt
  
  REAL(dp)                   :: start_time, end_time
  REAL(dp)                   :: interp_time
  REAL(dp)                   :: maxerror
  INTEGER                    :: irho,iz
  INTEGER                    :: ir, itheta

  !-------------------------------------------------!
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,*) '*******************************'
  WRITE(*,*) '*     cyl2sph_ex2 example     *'
  WRITE(*,*) '*******************************'
  WRITE(*,*)
  
  ! Set number of points
  numrhopts = 200
  numzpts   = 500
  
  ! Set grid spacing
  drho = 0.01
  dz   = 0.02
  
  ! Create axes
  ALLOCATE(rho_ax(1:numrhopts))
  ALLOCATE(z_ax(1:numzpts))
  
  DO irho =1, numrhopts
     rho_ax(irho) = REAL(irho,dp) * drho
  ENDDO
  
  DO iz = 1, numzpts
     z_ax(iz) = REAL(-numzpts / 2 + iz,dp) * dz  
  ENDDO
  
  
  ! We are going to set the 3d hydrogen function
  ! n=3, l=2, m=0
  ALLOCATE(Y_lm(1:numrhopts,1:numzpts))
  ALLOCATE(R_nl(1:numrhopts,1:numzpts))
  ALLOCATE(cylfunc(1:numrhopts,1:numzpts))
  
  DO iz = 1, numzpts
     DO irho = 1, numrhopts
        
        rpt = SQRT(rho_ax(irho)**2 + z_ax(iz)**2)
        thetapt = ACOS(z_ax(iz) / rpt)

        ! Spherical harmonic l=2, m=0
        Y_lm(irho,iz) = SQRT(5.0_dp / 4.0_dp / pi) * &
             ( 3.0_dp * (COS(thetapt))**2 - 1.0_dp)
        
        ! Spherical harmonic l=0, m=0
!!$        Y_lm(irho,iz) = SQRT(1.0_dp / 4.0_dp / pi)
        
        R_nl(irho,iz) = 1.0_dp
        
        cylfunc(irho,iz) = EXP(-rpt / 2.0_dp) * R_nl(irho,iz) * &
             Y_lm(irho,iz)
        
     ENDDO
  ENDDO
  
  ! Initialize the boundary
  Rboundary = 1.0_dp
  tolerance = 0.05_dp
  dr = 0.01_dp
  lmax      = 5
  dims      = (/numrhopts, numzpts/)
  
  WRITE(*,*)
  WRITE(*,*)            '------------------------------------'
  WRITE(*,'(A30,I4)') 'Number of points in rho: ',numrhopts
  WRITE(*,'(A30,I4)') 'Number of points in z: ',numzpts
  
  WRITE(*,'(A30,F9.4)') 'Grid spacing in rho: ',drho
  WRITE(*,'(A30,F9.4)') 'Grid spacing in z: ',dz
  
  WRITE(*,'(A30,F9.4)') 'Radius boundary: ',Rboundary
  WRITE(*,'(A30,F9.4)') 'Radius tolerance: ',tolerance
  WRITE(*,*)            '------------------------------------'
  
  ! Build interpolant
  WRITE(*,*) 'Creating interpolant...'
  CALL cpu_time(start_time)
  
  CALL initialize_cylindrical_boundary(rho_ax, z_ax, dims, &
       Rboundary, tolerance, 2, dr, lmax, numpts, &
       numrpts, numthetapts )
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time

  WRITE(*,'(A40,F9.4)') 'Interpolant time (seconds): ', interp_time
  WRITE(*,*)            '------------------------------------'
  
  WRITE(*,'(A40,I4)') 'Number of interpolation points: ',numpts  
  WRITE(*,'(A40,I4)') 'Number of radial boundary points: ',numrpts
  WRITE(*,'(A40,I4)') 'Number of polar boundary points: ',numthetapts
  WRITE(*,'(A40,F9.4)') 'Grid spacing in r: ',dr
  WRITE(*,'(A40,I4)') 'Maximum angular momenta: ',lmax
  WRITE(*,*)            '------------------------------------'
  
  ALLOCATE(sphfunc(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dr(1:numrpts,1:numthetapts))
  ALLOCATE(sphfunc_dth(1:numrpts,1:numthetapts))
  
  ! Interpolate!!
  WRITE(*,*) 'Interpolating boundary...'
  CALL cpu_time(start_time)
  
  CALL get_cylindrical_boundary(cylfunc, sphfunc, sphfunc_dr, &
       sphfunc_dth, 'quadratic')
  
  CALL cpu_time(end_time)
  
  interp_time = end_time - start_time
  
  WRITE(*,'(A40,F9.4)') 'Interpolation time (seconds): ', interp_time
  
  
  DO itheta = 1, numthetapts
     DO ir = 1, numrpts
!!$        ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) *  0.25_dp * SQRT(5.0_dp / pi) * &
!!$             ( 3.0_dp * (costheta_boundary(itheta))**2 - 1.0_dp)
!!$        ref_value =  EXP(-rpts_boundary(ir) / 2.0_dp) * 0.5_dp * SQRT(1.0_dp / pi)	
        ref_value = SQRT(1.0_dp / 4.0_dp / pi)	
        
!!$        WRITE(*,*) sphfunc(ir,itheta)
!!$        WRITE(*,*) ref_value 
!!$        WRITE(*,*) 'Diff: ', ABS(sphfunc(ir,itheta) - ref_value)
        maxerror = MAX(maxerror, ABS(sphfunc(ir,itheta) - ref_value))
     ENDDO
  ENDDO
  
  WRITE(*,*) 
  WRITE(*,*) '------------RESULTS--------------'
  WRITE(*,*)
  WRITE(*,'(A40,E15.8)') 'Maximum difference in value: ',maxerror
  WRITE(*,*)
  WRITE(*,*)          '----------------------------------------'
  WRITE(*,*)
  WRITE(*,*) 'Now the Spherical Harmonic transform! (SHT)...'
  WRITE(*,*)
  WRITE(*,*)


  ALLOCATE(radfunc(1:numrpts,0:lmax))
  
  ! Initialize Spherical Harmonics
  CALL initialize_spherical_harmonics(lmax, costheta_boundary)
  
  ! Make SHT
  CALL make_sht(sphfunc, costheta_boundary, theta_weights, lmax, radfunc)
  
  
!!!!!!!!!!!!!!!!!!!
!!!! Save axes !!!!
!!!!!!!!!!!!!!!!!!!
  
  ! Rho
  OPEN(unit=33,form='formatted',file='./results/rho_ax.dat')
  DO irho = 1, numrhopts
     WRITE(33,*) rho_ax(irho) 
  ENDDO
  CLOSE(33)
  
  ! z
  OPEN(unit=33,form='formatted',file='./results/z_ax.dat')
  DO iz = 1, numzpts
     WRITE(33,*) z_ax(iz)
  ENDDO
  CLOSE(33)
  
  ! R
  OPEN(unit=33,form='formatted',file='./results/r_ax.dat')
  DO ir = 1, numrpts
     WRITE(33,*) rpts_boundary(ir)
  ENDDO
  CLOSE(33)
  
  ! Theta
  OPEN(unit=33,form='formatted',file='./results/theta_ax.dat')
  DO itheta = 1, numthetapts
     WRITE(33,*) theta_boundary(itheta), costheta_boundary(itheta)
  ENDDO
  CLOSE(33)
  
  ! Save original wavefunction
  OPEN(unit=33,form='formatted',file='./results/cylfunction.dat')
  
  DO iz = 1, numzpts
     DO irho = 1, numrhopts
        WRITE(33,*) ABS(cylfunc(irho,iz))**2 
     ENDDO
  ENDDO
  
  CLOSE(33)
  
  ! Save spherical wavefunction
  OPEN(unit=33,form='formatted',file='./results/sphfunction.dat')
  
  DO itheta = 1, numthetapts
     DO ir = 1, numrpts
        WRITE(33,*) ABS(sphfunc(ir,itheta))**2 
     ENDDO
  ENDDO
  
  CLOSE(33)

  ! Save Legendre polynomials
  OPEN(unit=33,form='formatted',file='./results/legendrepoly.dat')
  
  DO il = 0, lmax
     DO itheta = 0, numthetapts-1
        WRITE(33,*) il, legenpl(itheta, 0, il)
     ENDDO
  ENDDO

  CLOSE(33)

  ! Save spherical wavefunction in sph harmonics basis
  OPEN(unit=33,form='formatted',file='./results/radfunc.dat')
  
  DO il = 0, lmax
     DO ir = 1, numrpts
        WRITE(33,*) REAL(radfunc(ir,il),dp),AIMAG(radfunc(ir,il)) 
     ENDDO
  ENDDO
  
  CLOSE(33)
  
  ! Free memory
  DEALLOCATE(rho_ax,z_ax)
  DEALLOCATE(Y_lm, R_nl)
  DEALLOCATE(cylfunc)
  DEALLOCATE(sphfunc, sphfunc_dr, sphfunc_dth)

END PROGRAM
