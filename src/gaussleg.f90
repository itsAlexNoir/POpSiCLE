!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: gaussleg
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 24/09/2015
!
! DESCRIPTION:
!> \brief Gauss legendre related routines
!> \details This module contains routines to prepare and perform
!> Gauss-Legendre quadratures. Also, it calculates Legendre polynomials,
!> Clebsch-Gordon coeffcients and spherical harmonics transforms.
!
!------------------------------------------------------------------------------
MODULE gaussleg

  USE constants_pop

  IMPLICIT NONE

  PRIVATE

  PUBLIC      :: get_gauss_stuff
  PUBLIC      :: make_gauss_lobatto
  PUBLIC      :: make_rnormal
  PUBLIC      :: make_legendre
  PUBLIC      :: plgndr
  PUBLIC      :: make_factlog
  PUBLIC      :: clebsch

  INTERFACE plgndr
     MODULE PROCEDURE plgndr_s
     MODULE PROCEDURE plgndr_v
  END INTERFACE plgndr

CONTAINS

  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE get_gauss_stuff
  !
  !>  \brief This subroutine returns the Gauss nodes of a particular
  !>  range of values, and the correspondent Gauss weights
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lower_bound Minimum value of the passed array
  !> \param[in] upper_bound Maximum value of the passed array
  !> \param[out] axis Returned axis with gaussian nodes.
  !> \param[out] weights Returned weights array.
  !
  !=======================================================================
  !=======================================================================

  SUBROUTINE get_gauss_stuff(lower_bound, upper_bound, axis, weights)

    IMPLICIT NONE

    REAL(dp), INTENT(IN)         :: lower_bound
    REAL(dp), INTENT(IN)         :: upper_bound
    REAL(dp), INTENT(OUT)        :: axis(:)
    REAL(dp), INTENT(OUT)        :: weights(:)

    INTEGER                      :: i,j
    INTEGER                      :: n, m, its
    REAL(dp), PARAMETER          :: precision = 1.0e-15_dp
    INTEGER, PARAMETER           :: MAXIT=10
    REAL(dp)                     :: xm, xl
    REAL(dp), DIMENSION((size(axis)+1)/2) :: z1, z
    REAL(dp), DIMENSION((size(axis)+1)/2) :: pp, p3, p2, p1
    LOGICAL, DIMENSION((size(axis)+1)/2)  :: unfinished

    ! Check the size of the arrays. They must be equal
    IF (SIZE(axis).NE.SIZE(weights)) THEN
       WRITE(*,*) 'Axis and weights arrays have not equal length.'
       STOP
    ENDIF

    n=SIZE(axis)
    m=(n+1)/2
    xm=0.5_dp*(upper_bound + lower_bound)
    xl=0.5_dp*(upper_bound - lower_bound)
    z=COS( pi * (arth(1,1,m) - 0.25_dp) / ( n + 0.5_dp))

    unfinished=.TRUE.
    DO its=1,MAXIT
       WHERE (unfinished)
          p1=1.0
          p2=0.0
       END WHERE
       DO j=1,n
          WHERE (unfinished)
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
          END WHERE
       ENDDO

       WHERE (unfinished)
          pp=n*(z*p1-p2)/(z*z-1.0_dp)
          z1=z
          z=z1-p1/pp
          unfinished=(ABS(z-z1) > precision)
       END WHERE
       IF (.NOT. ANY(unfinished)) EXIT
    ENDDO

    IF (its == MAXIT+1) THEN
       WRITE(*,*) 'Too many iterations in gaussleg!'
       STOP
    ENDIF

    axis(1:m)=xm-xl*z
    axis(n:n-m+1:-1)=xm+xl*z
    weights(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
    weights(n:n-m+1:-1)=weights(1:m)

  END SUBROUTINE get_gauss_stuff

  !------------------------------------------------------------------------------!
  
  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE make_gauss_lobatto
  !
  !>  \brief This subroutine returns the axis and weights for
  !   a Gauss-Lobatto quadrature. It is similar to the normal Gauss
  !   quadrature, but it includes the two end points as quadrature points.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] lb Minimum value of the passed array
  !> \param[in] ub Maximum value of the passed array
  !> \param[out] x Returned axis with gaussian nodes.
  !> \param[out] w Returned weights array.
  !
  !=======================================================================
  !=======================================================================
  
  SUBROUTINE make_gauss_lobatto(lb, ub, x, w)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) 	:: lb, ub
    REAL(dp), INTENT(OUT)	:: x(:)
    REAL(dp), INTENT(OUT)	:: w(:)
    
    INTEGER			:: n, nm1, iter, ix, iorder
    REAL(dp), ALLOCATABLE	:: leg(:, :)
    REAL(dp), ALLOCATABLE	:: xold(:)
    INTEGER, PARAMETER           :: maxiterations = 100000000
    REAL, PARAMETER              :: eps = 1.0e-15_dp
    
    ! Check the size of the arrays. They must be equal
    IF (SIZE(x).NE.SIZE(w)) THEN
       WRITE(*,*) 'Axis and weights arrays have not equal length.'
       STOP
    ENDIF
    
    n   = SIZE(x)
    nm1 = n - 1
    
    DO ix = 1, n
       x(ix) = COS(pi * REAL(nm1 + ix - 1, dp) / REAL(nm1, dp))
    ENDDO
    
    ALLOCATE(leg(1:n, 1:n))
    ALLOCATE(xold(1:n))
    
    leg  = 0.0
    xold = 2.0
    
    DO iter = 1, maxiterations
       
       IF (MAXVAL(ABS(x - xold)).LE.eps) EXIT
       
       xold = x;
       
       leg(:, 1) = 1.0_dp
       leg(:, 2) = x
       
       DO iorder = 2, nm1
          DO ix = 1, n
             
             leg(ix, iorder + 1) = (REAL(2 * iorder - 1, dp) * x(ix) *          &
                  leg(ix, iorder) -                           &
                  REAL(iorder - 1, dp) *                      &
                  leg(ix, iorder - 1)) / REAL(iorder, dp)
             
          ENDDO
       ENDDO
       
       DO ix = 1, n
          
          x(ix) = xold(ix) - (x(ix) * leg(ix, n) - leg(ix, nm1)) /             &
               (REAL(n, dp) * leg(ix, n))
          
       ENDDO
       
    ENDDO
    
    DO ix = 1, n   
       
       w(ix) = 2.0_dp / (REAL(nm1 * n, dp) * leg(ix, n) * leg(ix, n))
       
       x(ix) = 0.5_dp * x(ix) * (ub - lb) + 0.5_dp * (ub + lb)
       w(ix) = 0.5_dp * w(ix) * (ub - lb)
       
       WRITE(*, '(1X, I3, 2(1X, F20.16))') ix, x(ix), w(ix)
       
    ENDDO
    
    DEALLOCATE(leg)
    DEALLOCATE(xold)
    
  END SUBROUTINE make_gauss_lobatto
  
  !------------------------------------------------------------------------------!
  
  !=======================================================================
  !=======================================================================
  !
  !   FUNCTION arth
  !
  !>  \brief Returns an array of length n containing an arithmetic
  !>  progression whose first value is first and whose increment is increment.
  !>  If first and increment have rank greater than zero, returns an array
  !>  of one larger rank, with the last subscript having size n and
  !>  indexing the progressions.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[out] arth Array where the progression is stored.
  !> \param[in] first The forst element.
  !> \param[in] increment the increment of the progression.
  !> \param[in] n Number of elements in the progression.
  !
  !=======================================================================
  !=======================================================================

  FUNCTION arth(first,increment,n)

    INTEGER, INTENT(IN)    :: first, increment
    INTEGER, INTENT(IN)    :: n
    REAL(dp), DIMENSION(n) :: arth

    INTEGER                :: k

    IF( n > 0 ) arth(1) = first
    DO k= 2, n
       arth(k) = arth(k-1) + increment
    ENDDO

  END FUNCTION arth

  !------------------------------------------------------------------------------!

  !=======================================================================
  !=======================================================================
  !
  !   SUBROUTINE make_rnormal
  !
  !>  \brief This subroutine calculates the normalization co-efficients
  !>  of the Legendre polynomials for all the {il, im} combinations
  !>  which are required.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[out] rnormal The array which holds the normalization
  !>                   co-efficients as norm(im, il).
  !> \param[in] factlog The log of factorials.
  !> \param[in] MAXL The maximum value of the angular momentum quantum
  !>                 number.
  !> \param[in] NUMPTS The number of points at which the factorials are
  !>                  calculated at.
  !
  !=======================================================================
  !=======================================================================

  SUBROUTINE make_rnormal( rnormal, factlog, MAXL, NUMPTS )

    IMPLICIT NONE

    ! Passed variables
    INTEGER      :: MAXL, NUMPTS
    REAL(dp)     :: factlog(0:NUMPTS)
    REAL(dp)     :: rnormal(0:MAXL, 0:MAXL)

    ! Local variables
    INTEGER       :: il, im
    REAL(dp)      :: rintera, rinterb, rinterc, rinterd
    REAL(dp), PARAMETER :: oneovpi = 1.0_dp / pi

    DO il = 0, MAXL
       DO im = 0, il

          rintera = LOG(0.25_dp * oneovpi * DBLE(2 * il + 1))
          rinterb = factlog(il - im)
          rinterc = factlog(il + im)
          rinterd = 0.5_dp * (rintera + rinterb - rinterc)

          rnormal(im, il) = EXP(rinterd)

       ENDDO
    ENDDO

  END SUBROUTINE make_rnormal

  !------------------------------------------------------------------!

  !=======================================================================
  !=======================================================================
  !
  ! SUBROUTINE make_legendre
  !
  !> \brief This subroutine calculates the Legendre polynomials for all
  !> permissible {il,im} combinations which are given by
  !> il={0,MAXL} => im={0,il}. This routine is based on that given in
  !> Numerical Recipies in FORTRAN, 2nd Edition, page 246-248, in which
  !> the Legendre polynomial was calculated as a function value using a
  !> recurrence relation. This routine improves this relation by
  !> storing the intermediate Legendre polynomials which the
  !> recurrrence relation discarded.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[out] legenpl The array storing the Legendre polynomials given by:
  !>                   legenpl(x, im, il)
  !> \param[in] costheta The array storing the cosine of the angle at which
  !>                   the Legendre polynomial is calculated at, i.e.
  !>                   x={-1,1}
  !> \param[in] MAXL The maximum value of the single electron angular
  !>                   momentum quantum number.
  !> \param[in] MAXTHETAPTS The maximum theta point.
  !
  !=======================================================================
  !=======================================================================

  SUBROUTINE make_legendre( legenpl, costheta, MAXL, MAXTHETAPTS )

    IMPLICIT NONE

    ! Passed variables
    INTEGER      :: MAXL, MAXTHETAPTS
    REAL(dp)     :: legenpl(0:MAXTHETAPTS, 0:MAXL, 0:MAXL)
    REAL(dp)     :: costheta(0:MAXTHETAPTS)

    ! Local variables
    INTEGER       :: i, ill, imagqn, itheta
    REAL(dp)      :: somx2, fact, pmm, pmmp1, pll

    !Three nested do-loops. For each point the value of P(l, m) is
    !calculated.
    DO itheta = 0, MAXTHETAPTS

       DO imagqn = 0, MAXL

          !The value of P(im, im) is calculated first.
          pmm = 1.0_dp

          IF (imagqn.GT.0) THEN

             somx2 = SQRT((1.0_dp - costheta(itheta)) * &
                  (1.0_dp + costheta(itheta)))
             fact = 1.0_dp

             DO i = 1, imagqn
                pmm = -pmm * fact * somx2
                fact = fact + 2.0_dp
             ENDDO

          ENDIF

          legenpl(itheta, imagqn, imagqn) = pmm

          !If im + 1 </= MAXL then the value of  P(im, imm + 1) is
          !calculated.
          IF ((imagqn + 1).LE.MAXL) THEN
             pmmp1 = costheta(itheta) * DBLE(2 * imagqn + 1) * pmm
             legenpl(itheta, imagqn, imagqn + 1) = pmmp1
          ENDIF

          !If im + 2</ = MAXL then the value of P(im, {im+1, imaxm})
          !is calculated.
          IF ((imagqn + 2).LE.MAXL) THEN

             DO ill = imagqn + 2, MAXL
                pll = (costheta(itheta) * DBLE(2 * ill - 1) * pmmp1 - &
                     DBLE(ill + imagqn - 1) * pmm) / &
                     DBLE(ill - imagqn)
                legenpl(itheta, imagqn, ill) = pll
                pmm = pmmp1
                pmmp1 = pll
             ENDDO

          ENDIF

       ENDDO
    ENDDO

  END SUBROUTINE make_legendre

  !-----------------------------------------------------------------!

  FUNCTION plgndr_s(l,m,x)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: l,m
    REAL(dp), INTENT(IN) :: x
    REAL(dp)             :: plgndr_s

    INTEGER              :: ll
    REAL(dp)             :: pll,pmm,pmmp1,somx2

    IF(m<0 .OR. m>l .OR. ABS(x)>1.0_dp ) THEN
       WRITE(*,*) 'Input error in plgndr_s!'
       STOP
    ENDIF
    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0_dp-x)*(1.0_dp+x))
       pmm=PRODUCT(DBLE(arth(1,2,m)))*somx2**m
       IF (MOD(m,2) == 1) pmm=-pmm
    ENDIF

    IF (l == m) THEN
       plgndr_s=pmm
    ELSE
       pmmp1=x*(2*m+1)*pmm

       IF (l == m+1) THEN
          plgndr_s=pmmp1
       ELSE

          DO ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          plgndr_s=pll
       ENDIF
    ENDIF

  END FUNCTION plgndr_s


  !-----------------------------------------------------------------!

  FUNCTION plgndr_v(l,m,x)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: l,m
    REAL(dp), DIMENSION(:), INTENT(IN) :: x
    REAL(dp), DIMENSION(size(x)) :: plgndr_v

    INTEGER :: ll
    REAL(dp), DIMENSION(size(x)) :: pll,pmm,pmmp1,somx2

    IF(m<0 .OR. m>l .OR. ALL(ABS(x)>1.0_dp) ) THEN
       WRITE(*,*) 'Input error in plgndr_v!'
       STOP
    ENDIF

    pmm=1.0
    IF (m > 0) THEN
       somx2=SQRT((1.0_dp-x)*(1.0_dp+x))
       pmm=PRODUCT(DBLE(arth(1,2,m)))*somx2**m
       IF (MOD(m,2) == 1) pmm=-pmm
    ENDIF

    IF (l == m) THEN
       plgndr_v=pmm
    ELSE
       pmmp1=x*(2*m+1)*pmm
       IF (l == m+1) THEN
          plgndr_v=pmmp1
       ELSE
          DO ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          ENDDO
          plgndr_v=pll
       ENDIF
    ENDIF

  END FUNCTION plgndr_v

  !=======================================================================
  !=======================================================================
  !
  ! SUBROUTINE make_factlog
  !
  !> \brief This subroutine calculates the log of factorials. The log of the
  !> factorials is chosen since we can then use addition/subtraction
  !> instead of multiplication/division when calculating quantities
  !> involving factorials. This gives more stability.
  !
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[out] factlog The array storing the log of the factorials.
  !> \param[in]  NUMPTS  The number of factorials we wish to calculate.
  !> The log of the factorials for all integers in the range
  !> {0,NUMPTS} are calculated.
  !
  !=======================================================================
  !=======================================================================

  SUBROUTINE make_factlog( factlog, NUMPTS )

    IMPLICIT NONE

    ! Passed variables
    INTEGER     :: NUMPTS
    REAL(dp)    :: factlog(0:NUMPTS)

    ! Local variables
    INTEGER     ::  i, j
    REAL(dp)    ::  rinter

    ! The log of the factorial of zero is one in this routine.
    factlog(0) = 0.0_dp

    ! Calculate the log of the factorials of all other numbers.
    DO i = 1, NUMPTS
       rinter = 1.0_dp

       DO j = 1, i
          rinter = rinter * DBLE(j)
       ENDDO

       factlog(i) = LOG(rinter)

    ENDDO

  END SUBROUTINE make_factlog

  !----------------------------------------------------------!

  !=======================================================================
  !=======================================================================
  !
  ! SUBROUTINE clebsch
  !
  !> \brief This function calculates the Clebsch-Gordon co-efficient for
  !> given angular momentum and magnetic quantum numbers using an
  !> expression due to Racah---See Rose, M.E., `Elementary Theory of
  !> Angular Momentum', [1963], Wiley, pp 40.
  !>
  !======================SUBROUTINE ARGUMENTS=============================
  !
  !> \param[in] factlog The log of factorials.
  !> \param[in] NUMPTS The number of points at which the factorials are
  !>            calculated at.
  !> \param[in] j1 Single electron angular momentum quantum number for
  !>            electron 1.
  !> \param[in] j2 Single electron angular momentum quantum number for
  !>            electron 2.
  !> \param[in] ibigj Total angular momentum quantum number.
  !> \param[in] m1 Magnetic momentum quantum number for electron 1.
  !> \param[in] m2 Magnetic quantum number for electron 2.
  !> \param[in] ibigm Total magnetic quantum number.
  !
  !=======================================================================
  !=======================================================================

  FUNCTION clebsch( j1, j2, ibigj, m1, m2, ibigm, factlog, NUMPTS )

    IMPLICIT NONE

    ! Passed variables
    INTEGER, INTENT(IN)   :: NUMPTS, j1, m1, j2, m2, ibigj, ibigm
    REAL(dp), INTENT(IN)  :: factlog(0:NUMPTS)
    REAL(dp)              :: clebsch

    ! Local variables
    INTEGER                :: ik, kmin, kmax, ja, jb, jc, jd, je
    REAL(dp)               :: a
    REAL(dp), PARAMETER    :: zero = 0.0_dp
    REAL(dp), PARAMETER    :: half = 0.5_dp
    REAL(dp), PARAMETER    :: one = 1.0_dp

    clebsch = zero

    ! Only proceed if various triangularity constraints for the
    ! quantum numbers are fulfilled.
    IF (ABS(m1).GT.j1 .OR. ABS(m2).GT.j2 .OR. ABS(ibigm).GT.ibigj .OR. &
         (2 * MAX(j1, j2, ibigj)).GT.(j1 + j2 + ibigj) .OR.              &
         (m1 + m2).NE.ibigm) THEN

       clebsch = zero
       RETURN
       !WRITE(*,*) 'Error in quantum numbers'
       !STOP

    ELSE

       a = half * (factlog(j1 + j2 - ibigj) +                          &
            factlog(j2 + ibigj - j1) + factlog(ibigj + j1 - j2) -      &
            factlog(j1 + j2 + ibigj + 1) + factlog(j1 + m1) +          &
            factlog( j1 - m1) + factlog(j2 + m2) + factlog( j2 - m2) + &
            factlog(ibigj + ibigm) + factlog(ibigj - ibigm))

       ja = j1 + j2 - ibigj
       jb = j1 - m1
       jc = j2 + m2
       jd = ibigj - j2 + m1
       je = ibigj - j1 - m2
       kmin = MAX(0, -jd, -je)
       kmax = MIN(ja, jb, jc)

       DO ik = kmin, kmax
          clebsch = clebsch + (-1.0_dp) ** ik * EXP(a - factlog(ik) -  &
               factlog(ja - ik) - factlog(jb - ik) -                 &
               factlog(jc - ik) - factlog(jd + ik) -                 &
               factlog(je + ik))
       ENDDO

       clebsch = clebsch * SQRT(DBLE(2 * ibigj) + one)

    ENDIF

  END FUNCTION clebsch

  !----------------------------------------------------------------!

END MODULE gaussleg
