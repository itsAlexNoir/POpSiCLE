MODULE bessel
  
  USE constants
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC       :: bessj0
  PUBLIC       :: bessj1
  PUBLIC       :: bessjy
  PUBLIC       :: sphbessjy
  PUBLIC       :: sphj
  PUBLIC       :: csphjy
  
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

  ! Returns the Bessel FUNCTION J0(x) for any REAL x.
  
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
  
  ! Returns the Bessel FUNCTION J1(x) for any REAL x.
  
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
    ELSE
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
  
  !Returns the Bessel FUNCTIONs rj = Jν , ry = Yν and their derivatives rjp = Jν′ ,
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

  !Returns spherical Bessel FUNCTIONs jn(x), yn(x),
  ! and their derivatives jn′ (x), yn′ (x)
  ! for INTEGER n ≥ 0 and x > 0.
  
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
    
    DO j=m,2,-1
       sv=d
       d=y2*d-dd+c(j)
       dd=sv
    END DO
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
  
  ! Returns the sum of the REAL and the imaginary part
  ! of a COMPLEX number.

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
       DO i=n-1,1,-1
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
       END DO
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
  
  !**************************************************************************


  !*****************************************************************************80
  !
  !! SPHJ computes spherical Bessel FUNCTIONs jn(x) and their derivatives.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  ModIFied:
  !
  !    12 January 2016
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin.
  !    ModIFications suggested by Vincent Lagage, 12 January 2016.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special FUNCTIONs,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  PARAMETERs:
  !
  !    Input, INTEGER ( kind = 4 ) N, the order.
  !
  !    Input, REAL ( kind = 8 ) X, the argument.
  !
  !    Output, INTEGER ( kind = 4 ) NM, the highest order computed.
  !
  !    Output, REAL ( kind = 8 ) SJ(0:N), the values of jn(x).
  !
  !    Output, REAL ( kind = 8 ) DJ(0:N), the values of jn'(x).
  !
  
  SUBROUTINE sphj( n, x, nm, sj, dj )

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)   :: n
    REAL(dp), INTENT(IN)  :: x
    INTEGER, INTENT(OUT)  :: nm
    REAL(dp), INTENT(OUT) :: sj(0:n)
    REAL(dp), INTENT(OUT) :: dj(0:n)
    
    REAL(dp)              :: cs
    REAL(dp)              :: f
    REAL(dp)              :: f0
    REAL(dp)              :: f1
    INTEGER               :: k
    INTEGER               :: m
    REAL(dp)              :: sa
    REAL(dp)              :: sb

    !-------------------------------------------!
    
    nm = n
    !
    !  Original code.
    !
    IF ( .TRUE. ) THEN
       
       IF ( abs(x) <= 1.0D-100 ) THEN
          DO k = 0, n
             sj(k) = 0.0D+00
             dj(k) = 0.0D+00
          ENDDO
          sj(0) = 1.0D+00
          dj(1) = 0.3333333333333333D+00
          RETURN
       ENDIF
       !
       !  Updated code.
       !
    ELSE
       
       IF ( abs(x) <= 1.0D-16 ) THEN
          DO k = 0, n
             sj(k) = 0.0D+00
             dj(k) = 0.0D+00
          ENDDO
          sj(0) = 1.0D+00
          IF ( 0 < n ) THEN
             DO k = 1, n
                sj(k) = sj(k-1) * x / REAL ( 2 * k + 1, dp )
             ENDDO
             dj(1) = 1.0D+00 / 3.0D+00
          END IF
          RETURN
       ENDIF
       
    ENDIF
    
    sj(0) = SIN ( x ) / x
    sj(1) = ( sj(0) - COS ( x ) ) / x
    
    IF ( 2 <= n ) THEN
       
       sa = sj(0)
       sb = sj(1)
       m = msta1 ( x, 200 )
       IF ( m < n ) THEN
          nm = m
       ELSE
          m = msta2 ( x, n, 15 )
       ENDIF
       
       f0 = 0.0D+00
       f1 = 1.0D+00-100
       DO k = m, 0, -1
          f = ( 2.0D+00 * k + 3.0D+00 ) * f1 / x - f0
          IF ( k <= nm ) THEN
             sj(k) = f
          ENDIF
          f0 = f1
          f1 = f
       ENDDO
       
       IF( ABS( sa ) <= ABS( sb ) ) THEN
          cs = sb / f0
       ELSE
          cs = sa / f
       ENDIF
       
       DO k = 0, nm
          sj(k) = cs * sj(k)
       ENDDO
       
    ENDIF
    
    dj(0) = ( COS(x) - SIN(x) / x ) / x
    DO k = 1, nm
       dj(k) = sj(k-1) - ( k + 1.0D+00 ) * sj(k) / x
    ENDDO
    
    RETURN
    
  END SUBROUTINE sphj
  
!*****************************************************************************80
!
!! CSPHJY: spherical Bessel FUNCTIONs jn(z) and yn(z) for COMPLEX argument.
!
!  Discussion:
!
!    This procedure computes spherical Bessel FUNCTIONs jn(z) and yn(z)
!    and their derivatives for a COMPLEX argument.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  ModIFied:
!
!    01 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special FUNCTIONs,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
! 
!  PARAMETERs:
!
!    Input, INTEGER ( kind = 4 ) N, the order of jn(z) and yn(z).
!
!    Input, COMPLEX ( kind = 8 ) Z, the argument.
!
!    Output, INTEGER ( kind = 4 ) NM, the highest order computed.
!
!    Output, COMPLEX ( kind = 8 ) CSJ(0:N0, CDJ(0:N), CSY(0:N), CDY(0:N),
!    the values of jn(z), jn'(z), yn(z), yn'(z).
!

  SUBROUTINE csphjy ( n, z, nm, csj, cdj, csy, cdy )

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)       ::  n
    COMPLEX(dp), INTENT(IN)   ::  z
    INTEGER, INTENT(OUT)      ::  nm
    
    COMPLEX(dp), INTENT(OUT)  :: csj(0:n)
    COMPLEX(dp), INTENT(OUT)  :: cdj(0:n)
    COMPLEX(dp), INTENT(OUT)  :: csy(0:n)
    COMPLEX(dp), INTENT(OUT)  :: cdy(0:n)

    REAL(dp)                  :: a0
    COMPLEX(dp)               :: cf
    COMPLEX(dp)               :: cf0
    COMPLEX(dp)               :: cf1
    COMPLEX(dp)               :: cs
    COMPLEX(dp)               :: csa
    COMPLEX(dp)               :: csb
    INTEGER                   :: k
    INTEGER                   :: m
    
    !-------------------------------------------!
    
    a0 = abs(z)
    nm = n
    
    IF ( a0 < 1.0D-60 ) THEN
       DO k = 0, n
          csj(k) = 0.0_dp
          cdj(k) = 0.0_dp
          csy(k) = -1.0D+300
          cdy(k) = 1.0D+300
       ENDDO
       csj(0) = CMPLX( 1.0_dp, 0.0_dp, dp )
       cdj(1) = CMPLX( 0.333333333333333_dp, 0.0_dp, dp )
       RETURN
    ENDIF
    
    csj(0) = SIN ( z ) / z
    csj(1) = ( csj(0) - COS ( z ) ) / z
    
    IF ( 2 <= n ) THEN
       csa = csj(0)
       csb = csj(1)
       m = msta1 ( a0, 200 )
       IF ( m < n ) THEN
          nm = m
       ELSE
          m = msta2 ( a0, n, 15 )
       ENDIF
       cf0 = 0.0_dp
       cf1 = 1.0_dp-100
       DO k = m, 0, -1
          cf = ( 2.0_dp * k + 3.0_dp ) * cf1 / z - cf0
          IF ( k <= nm ) THEN
             csj(k) = cf
          ENDIF
          cf0 = cf1
          cf1 = cf
       ENDDO
       
       IF ( ABS ( csa ) <= ABS ( csb ) ) THEN
          cs = csb / cf0
       ELSE
          cs = csa / cf
       ENDIF
       
       DO k = 0, nm
          csj(k) = cs * csj(k)
       ENDDO
       
    ENDIF
    
    cdj(0) = ( COS ( z ) - SIN ( z ) / z ) / z
    DO k = 1, nm
       cdj(k) = csj(k-1) - ( k + 1.0_dp ) * csj(k) / z
    ENDDO
    csy(0) = - COS ( z ) / z
    csy(1) = ( csy(0) - SIN ( z ) ) / z
    cdy(0) = ( SIN ( z ) + COS ( z ) / z ) / z
    cdy(1) = ( 2.0_dp * cdy(0) - COS ( z ) )  / z
    
    DO k = 2, nm
       IF ( ABS ( csj(k-2) ) < ABS ( csj(k-1) ) ) THEN 
          csy(k) = ( csj(k) * csy(k-1) - 1.0_dp / ( z * z ) ) / csj(k-1)
       ELSE
          csy(k) = ( csj(k) * csy(k-2) &
               - ( 2.0_dp * k - 1.0_dp ) / z ** 3 ) / csj(k-2)
       ENDIF
    ENDDO
    
    DO k = 2, nm
       cdy(k) = csy(k-1) - ( k + 1.0_dp ) * csy(k) / z
    ENDDO
    
    RETURN
  END SUBROUTINE csphjy
  
  !***********************************!
  
  !*****************************************************************************80
  !
  !! R8_GAMMA_LOG evaluates the logarithm of the gamma FUNCTION.
  !
  !  Discussion:
  !
  !    This routine calculates the LOG(GAMMA) FUNCTION for a positive REAL
  !    argument X.  Computation is based on an algorithm outlined in
  !    references 1 and 2.  The program uses rational FUNCTIONs that
  !    theoretically approximate LOG(GAMMA) to at least 18 signIFicant
  !    decimal digits.  The approximation for X > 12 is from reference
  !    3, while approximations for X < 12.0 are similar to those in
  !    reference 1, but are unpublished.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license.
  !
  !  ModIFied:
  !
  !    15 April 2013
  !
  !  Author:
  !
  !    Original FORTRAN77 version by William Cody, Laura Stoltz.
  !    FORTRAN90 version by John Burkardt.
  !
  !  Reference:
  !
  !    William Cody, Kenneth Hillstrom,
  !    Chebyshev Approximations for the Natural Logarithm of the
  !    Gamma FUNCTION,
  !    Mathematics of Computation,
  !    Volume 21, Number 98, April 1967, pages 198-203.
  !
  !    Kenneth Hillstrom,
  !    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
  !    May 1969.
  !
  !    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
  !    Charles Mesztenyi, John Rice, Henry Thatcher,
  !    Christoph Witzgall,
  !    Computer Approximations,
  !    Wiley, 1968,
  !    LC: QA297.C64.
  !
  !  PARAMETERs:
  !
  !    Input, REAL ( kind = 8 ) X, the argument of the FUNCTION.
  !
  !    Output, REAL ( kind = 8 ) R8_GAMMA_LOG, the value of the FUNCTION.
  !

  FUNCTION r8_gamma_log ( x )
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN)       :: x
    REAL(dp)                   :: r8_gamma_log
    
    
    REAL(dp), DIMENSION ( 7 )  :: c = (/ &
         -1.910444077728D-03, &
         8.4171387781295D-04, &
         -5.952379913043012D-04, &
         7.93650793500350248D-04, &
         -2.777777777777681622553D-03, &
         8.333333333333333331554247D-02, &
         5.7083835261D-03 /)
    REAL(dp)               :: corr
    REAL(dp)               :: d1 = -5.772156649015328605195174D-01
    REAL(dp)               :: d2 = 4.227843350984671393993777D-01
    REAL(dp)               :: d4 = 1.791759469228055000094023_dp
    REAL(dp), PARAMETER    :: frtbig = 2.25D+76
    INTEGER                :: i
    REAL(dp), DIMENSION ( 8 ) :: p1 = (/ &
         4.945235359296727046734888D+00, &
         2.018112620856775083915565D+02, &
         2.290838373831346393026739D+03, &
         1.131967205903380828685045D+04, &
         2.855724635671635335736389D+04, &
         3.848496228443793359990269D+04, &
         2.637748787624195437963534D+04, &
         7.225813979700288197698961D+03 /)
    REAL(dp), DIMENSION ( 8 ) :: p2 = (/ &
         4.974607845568932035012064D+00, &
         5.424138599891070494101986D+02, &
         1.550693864978364947665077D+04, &
         1.847932904445632425417223D+05, &
         1.088204769468828767498470D+06, &
         3.338152967987029735917223D+06, &
         5.106661678927352456275255D+06, &
         3.074109054850539556250927D+06 /)
    REAL(dp), DIMENSION ( 8 ) :: p4 = (/ &
         1.474502166059939948905062D+04, &
         2.426813369486704502836312D+06, &
         1.214755574045093227939592D+08, &
         2.663432449630976949898078D+09, &
         2.940378956634553899906876D+10, &
         1.702665737765398868392998D+11, &
         4.926125793377430887588120D+11, &
         5.606251856223951465078242D+11 /)
    REAL(dp), DIMENSION ( 8 ) :: q1 = (/ &
         6.748212550303777196073036D+01, &
         1.113332393857199323513008D+03, &
         7.738757056935398733233834D+03, &
         2.763987074403340708898585D+04, &
         5.499310206226157329794414D+04, &
         6.161122180066002127833352D+04, &
         3.635127591501940507276287D+04, &
         8.785536302431013170870835D+03 /)
    REAL(dp), DIMENSION ( 8 ) :: q2 = (/ &
         1.830328399370592604055942D+02, &
         7.765049321445005871323047D+03, &
         1.331903827966074194402448D+05, &
         1.136705821321969608938755D+06, &
         5.267964117437946917577538D+06, &
         1.346701454311101692290052D+07, &
         1.782736530353274213975932D+07, &
         9.533095591844353613395747D+06 /)
    REAL(dp), DIMENSION ( 8 ) :: q4 = (/ &
         2.690530175870899333379843D+03, &
         6.393885654300092398984238D+05, &
         4.135599930241388052042842D+07, &
         1.120872109616147941376570D+09, &
         1.488613728678813811542398D+10, &
         1.016803586272438228077304D+11, &
         3.417476345507377132798597D+11, &
         4.463158187419713286462081D+11 /)
    REAL(dp)              :: res
    REAL(dp), PARAMETER   :: sqrtpi = 0.9189385332046727417803297_dp
    REAL(dp), PARAMETER   :: xbig = 2.55D+305
    REAL(dp)              :: xden
    REAL(dp), PARAMETER   :: xinf = 1.79D+308
    REAL(dp)              :: xm1
    REAL(dp)              :: xm2
    REAL(dp)              :: xm4
    REAL(dp)              :: xnum
    REAL(dp)              :: y
    REAL(dp)              :: ysq

    !-------------------------------------------------------!
    
    y = x
    
    IF ( 0.0_dp < y .AND. y <= xbig ) THEN
       
       IF ( y <= EPSILON ( y ) ) THEN
          
          res = - LOG ( y )
          !
          !  EPS < X <= 1.5.
          !
       ELSE IF ( y <= 1.5_dp ) THEN
          
          IF ( y < 0.6796875_dp ) THEN
             corr = -LOG ( y )
             xm1 = y
          ELSE
             corr = 0.0D+00
             xm1 = ( y - 0.5_dp ) - 0.5_dp
          ENDIF
          
          IF ( y <= 0.5_dp .or. 0.6796875_dp <= y ) THEN
             
             xden = 1.0_dp
             xnum = 0.0_dp
             DO i = 1, 8
                xnum = xnum * xm1 + p1(i)
                xden = xden * xm1 + q1(i)
             ENDDO
             
             res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )
             
          ELSE
             
             xm2 = ( y - 0.5_dp ) - 0.5_dp
             xden = 1.0_dp
             xnum = 0.0_dp
             DO i = 1, 8
                xnum = xnum * xm2 + p2(i)
                xden = xden * xm2 + q2(i)
             ENDDO
             
             res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )
             
          ENDIF
          !
          !  1.5 < X <= 4.0.
          !
       ELSE IF ( y <= 4.0_dp ) THEN
          
          xm2 = y - 2.0_dp
          xden = 1.0_dp
          xnum = 0.0_dp
          DO i = 1, 8
             xnum = xnum * xm2 + p2(i)
             xden = xden * xm2 + q2(i)
          ENDDO
          
          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
          !
          !  4.0 < X <= 12.0.
          !
       ELSE IF ( y <= 12.0_dp ) THEN
          
          xm4 = y - 4.0_dp
          xden = -1.0_dp
          xnum = 0.0_dp
          DO i = 1, 8
             xnum = xnum * xm4 + p4(i)
             xden = xden * xm4 + q4(i)
          ENDDO
          
          res = d4 + xm4 * ( xnum / xden )
          !
          !  Evaluate for 12 <= argument.
          !
       ELSE
          
          res = 0.0_dp
          
          IF ( y <= frtbig ) THEN
             
             res = c(7)
             ysq = y * y
             
             DO i = 1, 6
                res = res / ysq + c(i)
             ENDDO
             
          ENDIF
          
          res = res / y
          corr = LOG ( y )
          res = res + sqrtpi - 0.5_dp * corr
          res = res + y * ( corr - 1.0_dp )
          
       ENDIF
       !
       !  Return for bad arguments.
       !
    ELSE
       
       res = xinf
       
    ENDIF
    !
    !  Final adjustments and return.
    !
    r8_gamma_log = res
    
    RETURN
  END FUNCTION r8_gamma_log
  
  !*****************************************************************************80
  !
  !! ENVJ is a utility FUNCTION used by MSTA1 and MSTA2.
  !
  !  Discussion:
  !
  !    ENVJ estimates -log(Jn(x)) from the estimate
  !    Jn(x) approx 1/sqrt(2*pi*n) * ( e*x/(2*n))^n
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  ModIFied:
  !
  !    14 January 2016
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !    ModIFications suggested by Vincent Lafage, 11 January 2016.
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special FUNCTIONs,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  PARAMETERs:
  !
  !    Input, INTEGER ( kind = 4 ) N, the order of the Bessel FUNCTION.
  !
  !    Input, REAL ( kind = 8 ) X, the absolute value of the argument.
  !
  !    Output, REAL ( kind = 8 ) ENVJ, the value.
  !
  
  FUNCTION envj( n, x )
    
    IMPLICIT NONE
    
    REAL(dp)             :: envj
    INTEGER, INTENT(IN)  :: n
    REAL(dp), INTENT(IN) :: x
    
    REAL(dp)             :: logten
    REAL(dp)             :: n_r8
    
    !-----------------------------------!

    !
    !  Original code
    !
    IF ( .TRUE. ) THEN
       
       envj = 0.5_dp * LOG10 ( 6.28_dp * n ) &
            - n * LOG10 ( 1.36_dp * x / n )
       !
       !  ModIFication suggested by Vincent Lafage.
       !
    ELSE
       
       n_r8 = REAL ( n, dp )
       logten = LOG ( 10.0_dp )
       envj = r8_gamma_log ( n_r8 + 1.0_dp ) / logten - n_r8 * LOG10( x )
       
    ENDIF
    
    RETURN
    
  END FUNCTION envj
  
  !******************************!
  
  
  !*****************************************************************************80
  !
  !! MSTA1 determines a backward recurrence starting point for Jn(x).
  !
  !  Discussion:
  !
  !    This procedure determines the starting point for backward  
  !    recurrence such that the magnitude of    
  !    Jn(x) at that point is about 10^(-MP).
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  ModIFied:
  !
  !    08 July 2012
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special FUNCTIONs,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  PARAMETERs:
  !
  !    Input, REAL ( kind = 8 ) X, the argument.
  !
  !    Input, INTEGER ( kind = 4 ) MP, the negative logarithm of the 
  !    desired magnitude.
  !
  !    Output, INTEGER ( kind = 4 ) MSTA1, the starting point.
  !
  
  FUNCTION msta1 ( x, mp )
    
    IMPLICIT NONE

    INTEGER               :: msta1
    REAL(dp), INTENT(IN)  :: x
    INTEGER, INTENT(IN)   :: mp
    
    REAL(dp)              :: a0
    REAL(dp)              :: f
    REAL(dp)              :: f0
    REAL(dp)              :: f1
    INTEGER               :: it
    INTEGER               :: n0
    INTEGER               :: n1
    INTEGER               :: nn
    
    !------------------------------!
    
    a0 = ABS ( x )
    n0 = INT ( 1.1_dp * a0 ) + 1
    f0 = envj ( n0, a0 ) - mp
    n1 = n0 + 5
    f1 = envj ( n1, a0 ) - mp
    DO it = 1, 20       
       nn = n1 - ( n1 - n0 ) / ( 1.0_dp - f0 / f1 )                  
       f = envj ( nn, a0 ) - mp
       IF ( ABS ( nn - n1 ) < 1 ) THEN
          EXIT
       ENDIF
       n0 = n1
       f0 = f1
       n1 = nn
       f1 = f
    ENDDO
    
    msta1 = nn
    
    RETURN
  END FUNCTION msta1
  
  !*****************************************************************************80
  !
  !! MSTA2 determines a backward recurrence starting point for Jn(x).
  !
  !  Discussion:
  !
  !    This procedure determines the starting point for a backward
  !    recurrence such that all Jn(x) has MP signIFicant digits.
  !
  !    Jianming Jin supplied a modIFication to this code on 12 January 2016.
  !
  !  Licensing:
  !
  !    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
  !    they give permission to incorporate this routine into a user program 
  !    provided that the copyright is acknowledged.
  !
  !  ModIFied:
  !
  !    14 January 2016
  !
  !  Author:
  !
  !    Shanjie Zhang, Jianming Jin
  !
  !  Reference:
  !
  !    Shanjie Zhang, Jianming Jin,
  !    Computation of Special FUNCTIONs,
  !    Wiley, 1996,
  !    ISBN: 0-471-11963-6,
  !    LC: QA351.C45.
  !
  !  PARAMETERs:
  !
  !    Input, REAL ( kind = 8 ) X, the argument of Jn(x).
  !
  !    Input, INTEGER ( kind = 4 ) N, the order of Jn(x).
  !
  !    Input, INTEGER ( kind = 4 ) MP, the number of signIFicant digits.
  !
  !    Output, INTEGER ( kind = 4 ) MSTA2, the starting point.
  !
  FUNCTION msta2 ( x, n, mp )
    
    IMPLICIT NONE
    
    INTEGER                :: msta2
    REAL(dp), INTENT(IN)   :: x
    INTEGER, INTENT(IN)    :: n
    INTEGER, INTENT(IN)    :: mp
    
    REAL(dp)               :: a0
    REAL(dp)               :: ejn
    REAL(dp)               :: f
    REAL(dp)               :: f0
    REAL(dp)               :: f1
    REAL(dp)               :: hmp
    INTEGER                :: it
    INTEGER                :: n0
    INTEGER                :: n1
    INTEGER                :: nn
    REAL(dp)               :: obj

    !--------------------------------------!
    
    a0 = ABS ( x )
    hmp = 0.5_dp * mp
    ejn = envj ( n, a0 )
    
    IF ( ejn <= hmp ) THEN
       obj = mp
       !
       !  Original code:
       !
       !   n0 = int ( 1.1D+00 * a0 )
       !
       !  Updated code:
       !
       n0 = INT ( 1.1_dp * a0 ) + 1
    ELSE
       obj = hmp + ejn
       n0 = n
    ENDIF
    
    f0 = envj ( n0, a0 ) - obj
    n1 = n0 + 5
    f1 = envj ( n1, a0 ) - obj
    
    DO it = 1, 20
       nn = n1 - ( n1 - n0 ) / ( 1.0_dp - f0 / f1 )
       f = envj ( nn, a0 ) - obj
       IF ( ABS ( nn - n1 ) < 1 ) THEN
          EXIT
       ENDIF
       n0 = n1
       f0 = f1
       n1 = nn
       f1 = f
    ENDDO
    
    msta2 = nn + 10
    
    RETURN
  END FUNCTION msta2
  
  !*********************************!
  
  !*******************************************************************!
  
END MODULE bessel


