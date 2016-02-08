MODULE bessel
  
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: bessj0
  PUBLIC       :: bessj1
  PUBLIC       :: bessjy
  PUBLIC       :: sphbessjy
  
  ! Interfaces
  
  INTERFACE  bessj0
     MODULE PROCEDURE bessj0_s
     MODULE PROCEDURE bessj0_v
  END INTERFACE bessj0
  
  INTERFACE  bessj1
     MODULE PROCEDURE bessj1_s
     MODULE PROCEDURE bessj1_v
  END INTERFACE bessj1
  
  INTERFACE  sphbessjy
     MODULE PROCEDURE sphbessjy_s
     MODULE PROCEDURE sphbessjy_v
  END INTERFACE sphbessjy
  
  INTERFACE beschb
     MODULE PROCEDURE beschb_s
     MODULE PROCEDURE beschb_v
  END INTERFACE beschb
  
  INTERFACE  chebev
     MODULE PROCEDURE chebev_s
     MODULE PROCEDURE chebev_v
  END INTERFACE chebev
  
  INTERFACE poly
     MODULE PROCEDURE poly_dd,poly_ddv,poly_rc
     MODULE PROCEDURE poly_cc,poly_msk_ddv
  END INTERFACE poly
  
  INTEGER, PARAMETER      ::  NPAR_POLY=8
  
  !-------------------------------------------------------------------------!
  
CONTAINS

  ! Returns the Bessel function J0(x) for any real x.
  
  FUNCTION bessj0_s(x)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: bessj0_s
    REAL(dp) :: ax,xx,z
    REAL(dp) :: y
    REAL(dp), DIMENSION(5) :: p = (/1.0_dp,-0.1098628627e-2_dp,&
         0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
    REAL(dp), DIMENSION(5) :: q = (/-0.1562499995e-1_dp,&
         0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
         -0.934945152e-7_dp/)
    REAL(dp), DIMENSION(6) :: r = (/57568490574.0_dp,-13362590354.0_dp,&
         651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
         -184.9052456_dp/)
    REAL(dp), DIMENSION(6) :: s = (/57568490411.0_dp,1029532985.0_dp,&
         9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
    
    IF (ABS(x) < 8.0) THEN
       y=x**2
       bessj0_s=poly(y,r)/poly(y,s)
    ELSE
       ax=ABS(x)
       z=8.0_dp/ax
       y=z**2
       xx=ax-0.785398164_dp
       bessj0_s=SQRT(0.636619772_dp/ax)*(COS(xx) * &
            poly(y,p)-z*SIN(xx)*poly(y,q))
    ENDIF
  END FUNCTION bessj0_s
  
  !----------------------------------------------------!
  
  FUNCTION bessj0_v(x)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)          :: x(:)
    REAL(dp), DIMENSION(SIZE(x))  :: bessj0_v
    
    REAL(dp), DIMENSION(SIZE(x))  :: ax,xx,z
    REAL(dp), DIMENSION(SIZE(x))  :: y
    LOGICAL, DIMENSION(SIZE(x))   :: mask
    REAL(dp), DIMENSION(5)        :: p = (/1.0_dp,-0.1098628627e-2_dp,&
         0.2734510407e-4_dp,-0.2073370639e-5_dp,0.2093887211e-6_dp/)
    REAL(dp), DIMENSION(5)        :: q = (/-0.1562499995e-1_dp,&
         0.1430488765e-3_dp,-0.6911147651e-5_dp,0.7621095161e-6_dp,&
         -0.934945152e-7_dp/)
    REAL(dp), DIMENSION(6)        :: r = (/57568490574.0_dp,-13362590354.0_dp,&
         651619640.7_dp,-11214424.18_dp,77392.33017_dp,&
         -184.9052456_dp/)
    REAL(dp), DIMENSION(6)        :: s = (/57568490411.0_dp,1029532985.0_dp,&
         9494680.718_dp,59272.64853_dp,267.8532712_dp,1.0_dp/)
    
    mask = (ABS(x) < 8.0)
    WHERE (mask)
       y=x**2
       bessj0_v=poly(y,r,mask)/poly(y,s,mask)
    ELSEWHERE
       ax=ABS(x)
       z=8.0_dp/ax
       y=z**2
       xx=ax-0.785398164_dp
       bessj0_v = SQRT(0.636619772_dp/ax) * (COS(xx)* &
            poly(y,p,.NOT. mask)-z*SIN(xx)*poly(y,q,.NOT.mask))
    END WHERE
    
  END FUNCTION bessj0_v
  
  !-----------------------------------------------------------------------!
  
  ! Returns the Bessel function J1(x) for any real x.
  
  FUNCTION bessj1_s(x)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)   :: x
    REAL(dp)               :: bessj1_s
    REAL(dp)               :: ax,xx,z
    REAL(dp)               :: y
    REAL(dp), DIMENSION(6) :: r = (/72362614232.0_dp,&
         -7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
         15704.48260_dp,-30.16036606_dp/)
    REAL(dp), DIMENSION(6) :: s = (/144725228442.0_dp,2300535178.0_dp,&
         18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
    REAL(dp), DIMENSION(5) :: p = (/1.0_dp,0.183105e-2_dp,&
         -0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
    REAL(dp), DIMENSION(5) :: q = (/0.04687499995_dp,&
         -0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
         0.105787412e-6_dp/)
    
    IF (ABS(x) < 8.0_dp) THEN
       y=x**2
       bessj1_s=x*(poly(y,r)/poly(y,s))
    else
       ax=ABS(x)
       z=8.0_dp/ax
       y=z**2
       xx=ax-2.356194491_dp
       bessj1_s=SQRT(0.636619772_dp/ax)*(COS(xx) * &
            poly(y,p)-z*SIN(xx)*poly(y,q))*SIGN(1.0_dp,x)
    ENDIF
    
  END FUNCTION bessj1_s
  
!----------------------------------------!
  
  FUNCTION bessj1_v(x)
    
    IMPLICIT NONE
    
    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp), DIMENSION(SIZE(x))       :: bessj1_v
    
    REAL(dp), DIMENSION(SIZE(x))       :: ax,xx,z
    REAL(dp), DIMENSION(SIZE(x))       :: y
    LOGICAL, DIMENSION(SIZE(x))        :: mask
    REAL(dp), DIMENSION(6)             :: r = (/72362614232.0_dp,&
         -7895059235.0_dp,242396853.1_dp,-2972611.439_dp,&
         15704.48260_dp,-30.16036606_dp/)
    REAL(dp), DIMENSION(6)             :: s = (/144725228442.0_dp,2300535178.0_dp,&
         18583304.74_dp,99447.43394_dp,376.9991397_dp,1.0_dp/)
    REAL(dp), DIMENSION(5)             :: p = (/1.0_dp,0.183105e-2_dp,&
         -0.3516396496e-4_dp,0.2457520174e-5_dp,-0.240337019e-6_dp/)
    REAL(dp), DIMENSION(5)             :: q = (/0.04687499995_dp,&
         -0.2002690873e-3_dp,0.8449199096e-5_dp,-0.88228987e-6_dp,&
         0.105787412e-6_dp/)
    
    mask = (ABS(x) < 8.0)
    WHERE (mask)
       y=x**2
       bessj1_v = x*(poly(y,r,mask)/poly(y,s,mask))
    ELSEWHERE
       ax=ABS(x)
       z=8.0_dp/ax
       y=z**2
       xx=ax-2.356194491_dp
       bessj1_v = SQRT(0.636619772_dp/ax)*(COS(xx) * &
            poly(y,p,.NOT. mask)-z*SIN(xx)*poly(y,q,.NOT. mask)) * &
            SIGN(1.0_dp,x)
    END WHERE
    
  END FUNCTION bessj1_v
  
  !-----------------------------------------------------------------------!
  
  !Returns the Bessel functions rj = Jν , ry = Yν and their derivatives rjp = Jν′ ,
  !ryp = Yν′, for positive x and for xnu = ν ≥ 0
  
  SUBROUTINE bessjy(x,xnu,rj,ry,rjp,ryp)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: x,xnu
    REAL(dp), INTENT(OUT)      :: rj,ry,rjp,ryp

    INTEGER, PARAMETER         :: MAXIT=10000
    REAL(dp), PARAMETER        :: XMIN=2.0_dp
    REAL(dp), PARAMETER        :: EPS=1.0e-16_dp
    REAL(dp), PARAMETER        :: FPMIN=1.0e-30_dp
    INTEGER                    :: i, isign, l, nl
    REAL(dp)                   :: a,b,c,d,del,del1,e,f
    REAL(dp)                   :: fact,fact2,fact3,ff,gam,gam1,gam2
    REAL(dp)                   :: gammi,gampl,h,p,pimu,pimu2
    REAL(dp)                   :: q,r,rjl,rjl1,rjmu,rjp1,rjpl,rjtemp
    REAL(dp)                   :: ry1,rymu,rymup,rytemp,sum,sum1
    REAL(dp)                   :: w,x2,xi,xi2,xmu,xmu2
    COMPLEX(dp)                :: aa,bb,cc,dd,dl,pq
    
    
    IF(x.LT.0.0_dp .OR. xnu .LT. 0.0_dp) THEN
       WRITE(*,*) 'The argument and the order of bessel must be greater than zero.'
       STOP
    ENDIF
    
    nl=merge(int(xnu+0.5_dp), max(0,int(xnu-x+1.5_dp)), x < XMIN)
    
    xmu=xnu-nl
    xmu2=xmu*xmu
    xi=1.0_dp/x
    xi2=2.0_dp*xi
    
    w=xi2/pi
    isign=1
    h=xnu*xi
    IF (h < FPMIN) h=FPMIN
    b=xi2*xnu
    d=0.0
    c=h
    DO i=1,MAXIT
       b=b+xi2
       d=b-d
       IF (ABS(d) < FPMIN) d=FPMIN
       c=b-1.0_dp/c

       IF (ABS(c) < FPMIN) c=FPMIN
       d=1.0_dp/d
       del=c*d
       h=del*h
       IF (d < 0.0) isign=-isign
       IF (ABS(del-1.0_dp) < EPS) EXIT
    ENDDO
    IF (i > MAXIT) THEN
       WRITE(*,*) 'x too large in bessjy; try asymptotic expansion'
       STOP
    ENDIF
    
    rjl=isign*FPMIN
    rjpl=h*rjl
    rjl1=rjl
    rjp1=rjpl
    fact=xnu*xi

    DO l=nl,1,-1
       rjtemp=fact*rjl+rjpl
       fact=fact-xi
       rjpl=fact*rjtemp-rjl
       rjl=rjtemp
    ENDDO
    
    IF (rjl == 0.0) rjl=EPS
    f=rjpl/rjl
    IF (x < XMIN) THEN
       x2=0.5_dp*x
       pimu= pi * xmu
       IF (ABS(pimu) < EPS) THEN
          fact=1.0
       ELSE
          fact=pimu/SIN(pimu)
       ENDIF

       d=-LOG(x2)
       e = xmu * d
       IF (ABS(e) < EPS) THEN
          fact2=1.0
       ELSE
          fact2=SINH(e)/e
       ENDIF
       
       CALL beschb(xmu,gam1,gam2,gampl,gammi)

       ff=2.0_dp / pi *fact *(gam1*COSH(e)+gam2*fact2*d)
       e=EXP(e)
       p=e/(gampl*pi)
       q=1.0_dp/(e*pi*gammi)
       pimu2=0.5_dp*pimu
       IF (ABS(pimu2) < EPS) THEN
          fact3=1.0
       ELSE
          fact3=SIN(pimu2)/pimu2
       ENDIF
       r=pi*pimu2*fact3*fact3
       c=1.0
       d=-x2*x2
       sum=ff+r*q
       sum1=p
       
       DO i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*(ff+r*q)
          sum=sum+del
          del1=c*p-i*del
          sum1=sum1+del1
          IF (ABS(del) < (1.0_dp+ABS(sum))*EPS) EXIT
       ENDDO
       
       IF (i > MAXIT) THEN
          WRITE(*,*) 'bessy series failed to converge'
          STOP
       ENDIF
       
       rymu=-sum
       ry1=-sum1*xi2
       rymup=xmu*xi*rymu-ry1
       rjmu=w/(rymup-f*rymu)
       
    ELSE
       a=0.25_dp-xmu2
       pq=CMPLX(-0.5_dp*xi,1.0_dp,kind=dp)
       aa=CMPLX(0.0_dp,xi*a,kind=dp)
       bb=CMPLX(2.0_dp*x,2.0_dp,kind=dp)
       cc=bb+aa/pq
       dd=1.0_dp/bb
       pq=cc*dd*pq
       
       DO i=2,MAXIT
          a=a+2*(i-1)
          bb=bb+CMPLX(0.0_dp,2.0_dp,kind=dp)
          dd=a*dd+bb
          IF (ABSC(dd) < FPMIN) dd=FPMIN
          cc=bb+a/cc
          IF (ABSC(cc) < FPMIN) cc=FPMIN
          dd=1.0_dp/dd
          dl=cc*dd
          pq=pq*dl
          IF (ABSC(dl-1.0_dp) < EPS) EXIT
       ENDDO
       
       IF (i > MAXIT) THEN
          WRITE(*,*) 'cf2 failed in bessjy'
          STOP
       ENDIF
       p=REAL(pq)
       q=AIMAG(pq)
       gam=(p-f)/q
       rjmu=SQRT(w/((p-f)*gam+q))
       rjmu=SIGN(rjmu,rjl)
       rymu=rjmu*gam
       rymup=rymu*(p+q/gam)
       ry1=xmu*xi*rymu-rymup
       
    ENDIF
    fact=rjmu/rjl
    rj=rjl1*fact
    rjp=rjp1*fact
    
    DO i=1,nl
       rytemp=(xmu+i)*xi2*ry1-rymu
       rymu=ry1
       ry1=rytemp
    ENDDO
    ry=rymu
    ryp=xnu*xi*rymu-ry1
    
  END SUBROUTINE bessjy
  
  !-------------------------------------------------------!

  !Returns spherical Bessel functions jn(x), yn(x),
  ! and their derivatives jn′ (x), yn′ (x)
  ! for integer n ≥ 0 and x > 0.
  
  SUBROUTINE sphbessjy_s(n,x,sj,sy,sjp,syp)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)   :: n
    REAL(dp), INTENT(IN)  :: x
    REAL(dp), INTENT(OUT) :: sj,sy,sjp,syp
    REAL(dp)              :: rootpio2
    REAL(dp)              :: factor,order,rj,rjp,ry,ryp
    
    IF(n.LT.0 .AND. x .LT.0.0_dp) THEN
       WRITE(*, *) 'Argument and order must be greater than zero.'
       STOP
    ENDIF
    
    rootpio2 = SQRT( pi / 2.0_dp )
    order=n+0.5_dp
    call bessjy(x,order,rj,ry,rjp,ryp)
    factor=rootpio2/sqrt(x)
    sj=factor*rj
    sy=factor*ry
    sjp=factor*rjp-sj/(2.0_dp*x)
    syp=factor*ryp-sy/(2.0_dp*x)
    
  END SUBROUTINE sphbessjy_s
  
  !----------------------------------------!

  ! The parallel version of the routine above.
  
  SUBROUTINE sphbessjy_v(n,x,sj,sy,sjp,syp)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(IN)   :: n
    REAL(dp), INTENT(IN)  :: x(:)
    REAL(dp), INTENT(OUT) :: sj(:), sy(:), sjp(:), syp(:)
    REAL(dp)              :: rootpio2
    REAL(dp)              :: order
    REAL(dp), DIMENSION(size(x)) :: factor,rj,rjp,ry,ryp
    INTEGER               :: i
    
    IF((SIZE(x).NE.SIZE(sj)) .OR. (SIZE(x).NE.SIZE(sy)) &
         .OR. (SIZE(x).NE.SIZE(sjp)) &
         .OR. (SIZE(x).NE.SIZE(syp))) THEN
       WRITE(*, *) 'Input arrays have not the same length.'
       STOP
    ENDIF
    IF(n.LT.0 .AND. ALL(x < 0.0_dp)) THEN
       WRITE(*, *) 'Argument and order must be greater than zero.'
       STOP
    ENDIF
    
    rootpio2 = SQRT(pi / 2.0_dp )
    
    order=n+0.5_dp
    
    DO i = 1, SIZE(x)
       CALL bessjy(x(i),order,rj(i),ry(i),rjp(i),ryp(i))
    ENDDO
    
    factor=rootpio2/SQRT(x)
    sj=factor*rj
    sy=factor*ry
    sjp=factor*rjp-sj/(2.0_dp*x)
    syp=factor*ryp-sy/(2.0_dp*x)
    
  END SUBROUTINE sphbessjy_v
  
  !---------------------------------------------------------!
  
  ! Evaluates Γ1 and Γ2 by Chebyshev expansion for |x| ≤ 1/2.
  ! Also returns 1/Γ(1 + x) and 1 /Γ(1 − x).

  SUBROUTINE beschb_s(x,gam1,gam2,gampl,gammi)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)   :: x
    REAL(dp), INTENT(OUT)  :: gam1,gam2,gampl,gammi
    INTEGER, PARAMETER     :: NUSE1=7,NUSE2=8
    REAL(dp)               :: xx
    REAL(dp), DIMENSION(7) :: c1=(/-1.142022680371168_dp,&
         6.5165112670737e-3_dp,3.087090173086e-4_dp,-3.4706269649e-6_dp,&
         6.9437664e-9_dp,3.67795e-11_dp,-1.356e-13_dp/)
    REAL(dp), DIMENSION(8) :: c2=(/1.843740587300905_dp,&
         -7.68528408447867e-2_dp,1.2719271366546e-3_dp,&
         -4.9717367042e-6_dp, -3.31261198e-8_dp,2.423096e-10_dp,&
         -1.702e-13_dp,-1.49e-15_dp/)
    
    xx=8.0_dp*x*x-1.0_dp
    gam1=chebev(-1.0_dp,1.0_dp,c1(1:NUSE1),xx)
    gam2=chebev(-1.0_dp,1.0_dp,c2(1:NUSE2),xx)
    gampl=gam2-x*gam1
    gammi=gam2+x*gam1
    
  END SUBROUTINE beschb_s

!--------------------------------------------------!
  
  SUBROUTINE beschb_v(x,gam1,gam2,gampl,gammi)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)          :: x(:)
    REAL(dp), INTENT(OUT)         :: gam1(:), gam2(:), gampl(:), gammi(:)
    
    INTEGER, PARAMETER            :: NUSE1=7,NUSE2=8
    REAL(dp), DIMENSION(size(x))  :: xx
    REAL(dp), DIMENSION(7)        :: c1=(/-1.142022680371168_dp,&
         6.5165112670737e-3_dp,3.087090173086e-4_dp,-3.4706269649e-6_dp,&
         6.9437664e-9_dp,3.67795e-11_dp,-1.356e-13_dp/)
    REAL(dp), DIMENSION(8)        :: c2=(/1.843740587300905_dp,&
         -7.68528408447867e-2_dp,1.2719271366546e-3_dp,&
         -4.9717367042e-6_dp, -3.31261198e-8_dp,2.423096e-10_dp,&
         -1.702e-13_dp,-1.49e-15_dp/)
    
    xx=8.0_dp*x*x-1.0_dp
    gam1=chebev(-1.0_dp,1.0_dp,c1(1:NUSE1),xx)
    gam2=chebev(-1.0_dp,1.0_dp,c2(1:NUSE2),xx)
    gampl=gam2-x*gam1
    gammi=gam2+x*gam1
    
  END SUBROUTINE beschb_v

  !---------------------------------------------------------!

  ! Chebyshev evaluation.
  ! x Chebyshev coefficients
  
  FUNCTION chebev_s(a,b,c,x)

    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: a,b,x
    REAL(dp), INTENT(IN) :: c(:)
    REAL(dp)              :: chebev_s
    
    INTEGER :: j,m
    REAL(dp) :: d,dd,sv,y,y2

    IF ((x-a)*(x-b) > 0.0_dp) THEN
       WRITE(*,*) 'x not in range in chebev_s'
       STOP
    ENDIF
    
    m=size(c)
    d=0.0_dp
    dd=0.0_dp
    y=(2.0_dp*x-a-b)/(b-a)
    y2=2.0_dp*y
    
    do j=m,2,-1
       sv=d
       d=y2*d-dd+c(j)
       dd=sv
    end do
    chebev_s=y*d-dd+0.5_dp*c(1)

  END FUNCTION chebev_s

  !---------------------------------------!
  
  FUNCTION chebev_v(a,b,c,x)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)         :: a,b
    REAL(dp), INTENT(IN)         :: c(:), x(:)
    REAL(dp), DIMENSION(size(x)) :: chebev_v
    
    INTEGER                      :: j,m
    REAL(dp), DIMENSION(size(x)) :: d,dd,sv,y,y2
    
    IF (ANY((x-a)*(x-b) > 0.0)) THEN
       WRITE(*,*) 'x not in range in chebev_v'
       STOP
    ENDIF
    
    m=SIZE(c)
    d=0.0_dp
    dd=0.0_dp
    y=(2.0_dp*x-a-b)/(b-a)
    y2=2.0_dp*y
    DO j=m,2,-1
       sv=d
       d=y2*d-dd+c(j)
       dd=sv
    ENDDO
    chebev_v = y*d-dd+0.5_dp*c(1)
    
  END FUNCTION chebev_v
  
  !-------------------------------------------------------!
  
  ! Returns the sum of the real and the imaginary part
  ! of a complex number.

  FUNCTION absc(z)
    
    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN) :: z
    REAL(sp)                :: absc
    
    absc=ABS( REAL(z) ) + ABS( AIMAG(z) )
    
  END FUNCTION absc

  !-----------------------------------------------------!
  
  ! Polynomial evaluation.
  
  FUNCTION poly_dd(x,coeffs)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(IN) :: coeffs(:)
    REAL(dp)             :: poly_dd
    
    REAL(dp)              :: pow
    REAL(dp), ALLOCATABLE :: vec(:)
    INTEGER               :: i,n,nn
    
    n=SIZE(coeffs)
    IF (n <= 0) THEN
       poly_dd=0.0_dp
    ELSE IF (n < NPAR_POLY) THEN
       poly_dd=coeffs(n)
       DO i=n-1,1,-1
          poly_dd=x*poly_dd+coeffs(i)
       ENDDO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF(nn == 1) EXIT
          pow=pow*pow
          n=nn
       ENDDO
       poly_dd=vec(1)
       DEALLOCATE(vec)
    ENDIF
    
  END FUNCTION poly_dd
  
  !----------------------------!
  
  FUNCTION poly_rc(x,coeffs)

    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN)  :: x
    REAL(dp), INTENT(IN)     :: coeffs(:)
    COMPLEX(dp)              :: poly_rc

    COMPLEX(dp)              :: pow
    COMPLEX(dp), ALLOCATABLE :: vec(:)
    INTEGER                  :: i,n,nn
        
    n=SIZE(coeffs)
    IF (n <= 0) THEN
       poly_rc=0.0_dp
    ELSE IF (n < NPAR_POLY) THEN
       poly_rc=coeffs(n)
       
       DO i=n-1,1,-1
          poly_rc=x*poly_rc+coeffs(i)
       ENDDO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF (nn == 1) EXIT
          pow=pow*pow
          n=nn
       ENDDO
       poly_rc=vec(1)
       DEALLOCATE(vec)
    ENDIF
  END FUNCTION poly_rc
  
  !----------------------------!
  
  FUNCTION poly_cc(x,coeffs)

    IMPLICIT NONE
    
    COMPLEX(dp), INTENT(IN) :: x
    COMPLEX(dp), INTENT(IN) :: coeffs(:)
    COMPLEX(dp)             :: poly_cc

    COMPLEX(dp)             :: pow
    COMPLEX(dp), ALLOCATABLE :: vec(:)
    INTEGER                  :: i,n,nn

    n=SIZE(coeffs)
    IF(n <= 0) THEN
       poly_cc=0.0_dp
    ELSEIF (n < NPAR_POLY) THEN
       poly_cc=coeffs(n)
       do i=n-1,1,-1
          poly_cc=x*poly_cc+coeffs(i)
       ENDDO
    ELSE
       ALLOCATE(vec(n+1))
       pow=x
       vec(1:n)=coeffs
       DO
          vec(n+1)=0.0_dp
          nn=ishft(n+1,-1)
          vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
          IF(nn == 1) EXIT
          pow=pow*pow
          n=nn
       end do
       poly_cc = vec(1)
       DEALLOCATE(vec)
    ENDIF
    
  END FUNCTION poly_cc
  
  !----------------------------!
   
  FUNCTION poly_ddv(x,coeffs)

    IMPLICIT NONE
    
    REAL(DP), INTENT(IN) :: coeffs(:), x(:)
    REAL(DP), DIMENSION(SIZE(x)) :: poly_ddv

    INTEGER                      :: i,n,m
    
    m=SIZE(coeffs)
    n=SIZE(x)

    IF (m <= 0) THEN
       poly_ddv=0.0_dp
    ELSE IF (m < n .or. m < NPAR_POLY) THEN
       poly_ddv=coeffs(m)
       DO i=m-1,1,-1
          poly_ddv=x*poly_ddv+coeffs(i)
       ENDDO
    ELSE
       DO i=1,n
          poly_ddv(i)=poly_dd(x(i),coeffs)
       ENDDO
    ENDIF
    
  END FUNCTION poly_ddv
  
  !----------------------------!
    
  FUNCTION poly_msk_ddv(x,coeffs,mask)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)         :: coeffs(:), x(:)
    LOGICAL, INTENT(IN)          :: mask(:)
    REAL(dp), DIMENSION(SIZE(x)) :: poly_msk_ddv
    
    poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
    
  END FUNCTION poly_msk_ddv
  
  
END MODULE bessel
  
  
