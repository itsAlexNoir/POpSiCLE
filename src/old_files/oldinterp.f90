!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: interp
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 15/04/2015
!
! DESCRIPTION:
!> \brief Interpolating routines.
!> \details This module contains all the subprograms needed for
!> an interpolation of the wavefunction, including routines for
!> multidimensions grid. Many of the routines in this module have
!> been extract form the numerical recipes books.
!
!------------------------------------------------------------------------------

MODULE interp
  
  USE constants
  USE tools
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC      :: locate
  PUBLIC      :: hunt
  PUBLIC      :: set_dy2
  PUBLIC      :: polint
  PUBLIC      :: splint
  PUBLIC      :: polin2
  PUBLIC      :: splie2
  PUBLIC      :: splin2
  PUBLIC      :: bcucof
  !PUBLIC      :: tcucof
  PUBLIC      :: bcuint
  
CONTAINS

  !----------------------------------------------------------------------------
  !
  ! FUNCTION locate
  !
  !> \brief Searching routine based on the bisection method
  !> \details Given an array xx(1:N), and given a value x,
  !> returns a value j such that x is between xx(j) and xx(j + 1).
  !> xx must be monotonic, either increasing or decreasing.
  !> j = 0 or j = N is returned to indicate that x is out of range.
  !
  !> \param[in] xx Array of tabulated values
  !> \param[in] x Given value to be search
  !> \return locate Array index of xx in which x is centered.
  !
  !----------------------------------------------------------------------------
  
  FUNCTION locate(xx,x)
    
    IMPLICIT NONE
    
    REAL(dp), INTENT(IN) :: xx(:)
    REAL(dp), INTENT(IN) :: x
    INTEGER :: locate
    ! Variables
    INTEGER :: n,jl,jm,ju
    ! ascnd true if ascending, false if descending.
    LOGICAL :: ascnd
    
    n=size(xx)
    ! True if ascending order of table, false otherwise.
    ascnd = (xx(n) >= xx(1))
    !  Initialize lower
    jl=0 
    ! and upper limits.
    ju=n+1
    
    do
       !  Repeat until this condition is satisfied.
       if (ju-jl <= 1) exit
       !  Compute a midpoint,
       jm=(ju+jl)/2    
       if (ascnd .eqv. (x >= xx(jm))) then
          ! and replace either the lower limit
          jl=jm                             
       else
          !  or the upper limit, as appropriate.
          ju=jm
       end if
    end do
    !  Then set the output, being careful with the endpoints.
    if (x == xx(1)) then               
       locate=1
    else if (x == xx(n)) then
       locate=n-1
    else
       locate=jl
    end if
    
  END FUNCTION locate
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE hunt
  !
  !> \brief Searching routine for correlated values
  !> \details Given an array xx(1:N), and given a value x, returns a value jlo
  !> such that x is between xx(jlo) and xx(jlo+1). xx must be monotonic,
  !> either increasing or decreasing.
  !> jlo = 0 or jlo = N is returned to indicate that x is out of range.
  !> jlo on input is taken as the initial guess for jlo on output.
  !
  !> \param[in] xx Array of tabulated values
  !> \param[in] x Given value to be search
  !> \param[out] jlo Array index of xx in which x is centered.
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE hunt(xx,x,jlo)
    
    IMPLICIT NONE
    
    INTEGER, INTENT(INOUT) :: jlo
    REAL(dp), INTENT(IN) :: x
    REAL(dp), DIMENSION(:), INTENT(IN) :: xx
    
    INTEGER :: n,inc,jhi,jm
    LOGICAL :: ascnd
    
    n=size(xx)
    ! True if ascending order of table, false otherwise.
    ascnd = (xx(n) >= xx(1))
    ! Input guess not useful. Go immediately to bisection
    if (jlo <= 0 .or. jlo > n) then
       jlo=0
       jhi=n+1
    else
       inc=1 !< Set the hunting increment.
       ! Hunt up:
       if (x >= xx(jlo) .eqv. ascnd) then
          do
             jhi=jlo+inc
             ! Done hunting, since off end of table.
             if (jhi > n) then             
                jhi=n+1
                exit
             else
                if (x < xx(jhi) .eqv. ascnd) exit
                ! Not done hunting, so double the increment
                jlo=jhi
                inc=inc+inc
             end if
          end do ! And try again
          ! Hunt down:
       else
          jhi=jlo
          do
             jlo=jhi-inc
             ! Done Hunting, since off end of table.
             if (jlo < 1) then
                jlo=0
                exit
             else
                if (x >= xx(jlo) .eqv. ascnd) exit
                ! Not done hunting, so double the increment
                jhi=jlo
                inc=inc+inc
             end if
          end do ! And try again
       end if
    end if   ! Done hunting, value bracketed
    
    ! Hunt is done, so begin the final bisection phase:
    do
       if (jhi-jlo <= 1) then
          if (x == xx(n)) jlo=n-1
          if (x == xx(1)) jlo=1
          exit
       else
          jm=(jhi+jlo)/2
          if (x >= xx(jm) .eqv. ascnd) then
             jlo=jm
          else
             jhi=jm
          end if
       end if
    end do
    
  END SUBROUTINE hunt

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE polint
  !
  !> \brief Polynomial interpolation
  !> \details Given arrays xa and ya of length N, and given a value x,
  !> this routine returns a value y, and an error estimate dy.
  !> If \f$P(x)\f$ is the polynomial of degree N − 1 such that \f$P(xa_i)=ya_i\f$,
  !> \f$i=1,...,N\f$, then the returned value \f$y=P(x)\f$.
  !
  !> \param[in] xa Known points
  !> \param[in] ya Tabulated values
  !> \param[in] x Point where the function is interpolated
  !> \param[out] y Interpolated function value
  !> \param[out] dy Error estimation
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE polint(xa,ya,x,y,dy)
    
    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: xa(:), ya(:)
    REAL(dp), INTENT(IN) :: x
    REAL(dp), INTENT(OUT) :: y,dy
    INTEGER :: m,n,ns
    REAL(dp), DIMENSION(size(xa)) :: c,d,den,ho
    
    n=assert_eq(size(xa),size(ya),'polint')
    ! Initiialize the tableau of c’s and d’s.
    c=ya
    d=ya
    ho=xa-x
    ! find index ns of closest table entry.
    ns=iminloc(abs(x-xa))
    ! This is the initial approximation to y.
    y=ya(ns)
    ns=ns-1
    ! For each colum of the table,
    ! we loop over the current c's and d's and update them
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(den(1:n-m) == 0.0_dp)) &
            ! This error can occur only if two input xa’s are (to within roundoff) identical.
            call nrerror('polint: calculation failure')
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       ! Here the c’s and d’s are updated.
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       ! After each column in the tableau is completed, we decide which correction, c or d,
       ! we want to add to our accumulating value of y, i.e., which path to take through
       ! the tableau—forking up or down. We do this in such a way as to take the
       ! most “straight line” route through the tableau to its apex, updating ns accordingly
       ! to keep track of where we are. This route keeps the partial approxima-
       ! tions centered (insofar as possible) on the target x. The last dy added is thus
       ! the error indication.
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do

    
  END SUBROUTINE polint

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE splint
  !
  !> \brief Spline interpolation
  !> \details Given the arrays xa and ya, which tabulate a function
  !> (with the \f$xa_i\f$’s in increasing or decreasing order), and given
  !> the array \f$y2_a\f$, which is the output from spline above, and given
  !> a value of x, this routine returns a cubic-spline interpolated value.
  !> The arrays xa, ya and y2a are all of the same size.
  !
  !> \param[in] xa Known points
  !> \param[in] ya Tabulated values
  !> \param[in] y2a Second derivative of the tabulated function
  !> \param[in] x Point where the function is interpolated 
  !> \return spline Interpolated value
  !
  !
  !-----------------------------------------------------------------------------
  
  
  FUNCTION splint(xa,ya,y2a,x)
    
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: xa(:), ya(:), y2a(:)
    REAL(dp), INTENT(IN) :: x
    REAL(dp) :: splint

    INTEGER :: khi,klo,n
    REAL(dp) :: a,b,h
    n=assert_eq(size(xa),size(ya),size(y2a),'splint')
    klo=max(min(locate(xa,x),n-1),1)

    ! We will find the right place in the table by means of locate’s bisection algorithm.
    ! This is optimal if sequential calls to this routine are at random values of x.
    ! If sequential calls are in order, and closely spaced, one would do better to store
    ! previous values of klo and khi and test if they remain appropriate on the next call.

    khi=klo+1 ! klo and khi now bracket the input value of x.

    h=xa(khi)-xa(klo)
    ! The xa’s must be distinct.
    if (h == 0.0) call nrerror('bad xa input in splint')
    ! Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h
    splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp

  END FUNCTION splint

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE set_dy2
  !
  !> \brief Calculate derivatives for the tabulated function.
  !> \details Given arrays \f$x\f$ and \f$y\f$ of length N containing a tabulated
  !> function, i.e., \f$y_i=f(x_i)\f$, with \f$x_1<x_2<\ldots< x_N\f$, and given
  !> values \f$yp_1\f$ and \f$yp_n\f$ for the first derivative of the interpolating
  !> function at points 1 and N, respectively, this routine returns an array y2 of length N
  !> that contains the second derivatives of the interpolating function at the
  !> tabulated points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger,
  !> the routine is signaled to set the corresponding boundary condition for
  !> a natural spline, with zero second derivative on that boundary.
  !
  !> \param[in] x Known points
  !> \param[in] y Tabulated values
  !> \param[in] yp1 Value of the first derivative at the first point
  !> \param[in] ypn Value of the first derivative at the last point
  !> \param[out] y2 Second derivative values
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE set_dy2(x,y,yp1,ypn,y2)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN)  :: x,y
    REAL(dp), INTENT(IN)                :: yp1,ypn
    REAL(dp), DIMENSION(:), INTENT(OUT) :: y2
    INTEGER :: n
    REAL(dp), DIMENSION(size(x))        :: a,b,c,r
    
    n=assert_eq(size(x),size(y),size(y2),'spline')
    ! Set up the tridiagonal matrix
    c(1:n-1)=x(2:n)-x(1:n-1)
    r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
    r(2:n-1)=r(2:n-1)-r(1:n-2)
    a(2:n-1)=c(1:n-2)
    b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
    b(1)=1.0_dp
    b(n)=1.0_dp
    ! The lower boundary condition is set either to be “natural”
    if (yp1 > 0.99e30_dp) then
       r(1)=0.0_dp
       c(1)=0.0_dp
       ! or else to have a specified first derivative.
    else
       r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
       c(1)=0.5_dp
    end if
    ! The upper boundary condition is set either to be “natural”
    if (ypn > 0.99e30_dp) then
       r(n)=0.0_dp
       a(n)=0.0_dp
       ! or else to have a specified first derivative.
    else
       r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
       a(n)=0.5_dp
    end if
    call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
    
  END SUBROUTINE set_dy2

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE polin2
  !
  !> \brief Polynomial interpolation in 2D
  !> \details Given arrays x1a of length M and x2a of length N of independent variables,
  !> and an \f$M \times N\f$ array of function values ya, tabulated at the grid points
  !> defined by x1a and x2a, and given values x1 and x2 of the independent variables,
  !> this routine returns an interpolated function value y, and an accuracy indication
  !> dy (based only on the interpolation in the x1 direction, however).
  !
  !> \param[in] x1a Grid points in the first dimension
  !> \param[in] x2a Grid points in the second dimension
  !> \param[in] ya Tabulated values
  !> \param[in] x1 Point where the function is interpolated
  !> \param[in] x2 Point where the function is interpolated
  !> \param[out] y Interpolated function value
  !> \param[out] dy Error estimation
  !
  !-----------------------------------------------------------------------------
 
  SUBROUTINE polin2(x1a,x2a,ya,x1,x2,y,dy)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: x1a,x2a
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: ya
    REAL(dp), INTENT(IN) :: x1,x2
    REAL(dp), INTENT(OUT) :: y,dy
    INTEGER :: j,m,ndum
    REAL(dp), DIMENSION(size(x1a)) :: ymtmp
    REAL(dp), DIMENSION(size(x2a)) :: yntmp
    
    m=assert_eq(size(x1a),size(ya,1),'polin2: m')
    ndum=assert_eq(size(x2a),size(ya,2),'polin2: ndum')
    ! Loop over rows.
    do j=1,m
       ! Copy row into temporary storage.
       yntmp=ya(j,:)
       ! Interpolate answer into temporary storage.
       call polint(x2a,yntmp,x2,ymtmp(j),dy)
    end do
    ! o the final interpolation.
    call polint(x1a,ymtmp,x1,y,dy)
  END SUBROUTINE polin2

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE splie2
  !
  !> \brief Gives the second derivative values of a 2D array
  !> \details Given an M × N tabulated function ya, and N tabulated independent
  !> variables x2a, this routine constructs one-dimensional natural cubic splines
  !> of the rows of ya and returns the second derivatives in the M × N array y2a.
  !> (The array x1a is included in the argument list merely for consistency with
  !> routine splin2.)
  !
  !> \param[in] x1a Grid points in the first dimension
  !> \param[in] x2a Grid points in the second dimension
  !> \param[in] ya Tabulated values
  !> \param[out] y2a Second derivative values of ya matrix
  !
  !-----------------------------------------------------------------------------

  SUBROUTINE splie2(x1a,x2a,ya,y2a)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: x1a,x2a
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: ya
    REAL(dp), DIMENSION(:,:), INTENT(OUT) :: y2a
    INTEGER :: j,m,ndum
    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')

    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
    do j=1,m
       call set_dy2(x2a,ya(j,:),1.0e30_dp,1.0e30_dp,y2a(j,:))
       ! Values 10e30 signal a natural spline.
    end do
  END SUBROUTINE splie2


  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE splin2
  !
  !> \brief Bicubic spline interpolation
  !> \details Given x1a, x2a, ya as described in splie2 and y2a as produced by
  !> that routine; and given a desired interpolating point x1,x2; this routine
  !> returns an interpolated function value by bicubic spline interpolation.
  !
  !> \param[in] x1a Grid points in the first dimension
  !> \param[in] x2a Grid points in the second dimension
  !> \param[in] ya Tabulated values
  !> \param[in] y2a Second derivative values of ya matrix
  !> \param[in] x1 Point where the function is interpolated
  !> \param[in] x2 Point where the function is interpolated
  !> \return splin2 Interpolated function value
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION splin2(x1a,x2a,ya,y2a,x1,x2)
    REAL(dp), DIMENSION(:), INTENT(IN) :: x1a,x2a
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: ya,y2a
    REAL(dp), INTENT(IN) :: x1,x2
    REAL(dp) :: splin2
    INTEGER :: j,m,ndum
    REAL(dp), DIMENSION(size(x1a)) :: yytmp,y2tmp2
    
    m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
    ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
    do j=1,m
       ! Perform m evaluations of the row splines constructed by splie2,
       ! using the one-dimensional spline evaluator splint.
       yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
    end do
    ! Construct the one-dimensional column spline and evaluate it.
    call set_dy2(x1a,yytmp,1.0e30_dp,1.0e30_dp,y2tmp2)
    splin2=splint(x1a,yytmp,y2tmp2,x1)
  END FUNCTION splin2

  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE bcucof
  !
  !> \brief Calculate coefficients for bicubic interpolation
  !> \details Given arrays y, y1, y2, and y12, each of length 4, containing
  !> the function, gradients, and cross derivative at the four grid points of a
  !> rectangular grid cell (numbered counterclockwise from the lower left), and
  !> given d1 and d2, the length of the grid cell in the 1- and 2- directions,
  !> this routine returns the 4 × 4 table c that is used by routine bcuint for
  !> bicubic interpolation.
  !
  !> \param[in] y Tabulated function
  !> \param[in] y1 Tabulated first derivative in x1 direction at four points
  !> \param[in] y2 Tabulated first derivative in x2 direction at four points
  !> \param[in] y12 Cross derivative at points
  !> \param[in] d1 Length of the grid cell in x1 direction
  !> \param[in] d2 Length of the grid cell in x2 direction
  !> \param[out] c coefficients table
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE bcucof(y,y1,y2,y12,d1,d2,c)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: d1,d2
    REAL(dp), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
    REAL(dp), DIMENSION(4,4), INTENT(OUT) :: c
    REAL(dp), DIMENSION(16) :: x
    REAL(dp), DIMENSION(16,16) :: wt
    DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
         8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
         2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
         2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
         -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
         -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
         -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
    ! Pack a temporary vector x.
    x(1:4)=y
    x(5:8)=y1*d1
    x(9:12)=y2*d2
    x(13:16)=y12*d1*d2
    
    ! Matrix multiply by the stored table.
    x=matmul(wt,x)
    ! Unpack the result into the output table.
    c=reshape(x,(/4,4/),order=(/2,1/))
  END SUBROUTINE bcucof

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE tcucof
  !
  !> \brief Calculate coefficients for tricubic interpolation
  !> \details Given arrays y, y1, y2, and y12, each of length 4, containing
  !> the function, gradients, and cross derivative at the four grid points of a
  !> rectangular grid cell (numbered counterclockwise from the lower left), and
  !> given d1 and d2, the length of the grid cell in the 1- and 2- directions,
  !> this routine returns the 4 × 4 table c that is used by routine bcuint for
  !> bicubic interpolation.
  !
  !> \param[in] y Tabulated function
  !> \param[in] y1 Tabulated first derivative in x1 direction at four points
  !> \param[in] y2 Tabulated first derivative in x2 direction at four points
  !> \param[in] y3 Tabulated first derivative in x2 direction at four points
  
  !> \param[in] y123 Cross derivative at points
  !> \param[in] d1 Length of the grid cell in x1 direction
  !> \param[in] d2 Length of the grid cell in x2 direction
  !> \param[in] d3 Length of the grid cell in x3 direction
  !> \param[out] c coefficients table
  !
  !-----------------------------------------------------------------------------
  
!!$  SUBROUTINE tcucof(y,y1,y2,y12,d1,d2,c)
!!$    IMPLICIT NONE
!!$    REAL(dp), INTENT(IN) :: d1,d2
!!$    REAL(dp), DIMENSION(8), INTENT(IN) :: y,y1,y2,y12
!!$    REAL(dp), DIMENSION(8,8), INTENT(OUT) :: c
!!$    REAL(dp), DIMENSION(64) :: x
!!$    REAL(dp), DIMENSION(64,64) :: wt
!!$
!!$    DATA wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
!!$         8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
!!$         2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
!!$         2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
!!$         -3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
!!$         -1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
!!$         -1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
!!$    ! Pack a temporary vector x.
!!$    x(1:8)=y
!!$    x(5:8)=y1*d1
!!$    x(9:12)=y2*d2
!!$    x(13:16)=y12*d1*d2
!!$    
!!$    ! Matrix multiply by the stored table.
!!$    x=matmul(wt,x)
!!$    ! Unpack the result into the output table.
!!$    c=reshape(x,(/4,4/),order=(/2,1/))
!!$  END SUBROUTINE bcucof
    
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE bcuint
  !
  !> \brief Bicubic interpolation
  !> \details Bicubic interpolation within a grid square. Input quantities
  !> are y,y1,y2,y12 (as described in bcucof); x1l and x1u, the lower and upper
  !> coordinates of the grid square in the 1- direction; x2l and x2u likewise for
  !> the 2-direction; and x1,x2, the coordinates of the desired point for the
  !> interpolation. The interpolated function value is returned as ansy, and
  !> the interpolated gradient values as ansy1 and ansy2. This routine calls bcucof.
  !
  !> \param[in] y Tabulated function
  !> \param[in] y1 Tabulated first derivative in x1 direction at four points
  !> \param[in] y2 Tabulated first derivative in x2 direction at four points
  !> \param[in] y12 Cross derivative at points
  !> \param[in] x1l Lower coordinate of the grid cell in the first direction
  !> \param[in] x1u Upper coordinate of the grid cell in the first direction
  !> \param[in] x2l Lower coordinate of the grid cell in the second direction
  !> \param[in] x2u Upper coordinate of the grid cell in the second direction
  !> \param[in] x1 Coordinate in the first direction for the interpolated point
  !> \param[in] x2 Coordinate in the second direction for the interpolated point
  !> \param[out] ansy Interpolated function
  !> \param[out] ansy1 Interpolated first derivative in the first direction
  !> \param[out] ansy2 Interpolated first derivative in the first direction
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
    IMPLICIT NONE
    REAL(dp), DIMENSION(4), INTENT(IN) :: y,y1,y2,y12
    REAL(dp), INTENT(IN) :: x1l,x1u,x2l,x2u,x1,x2
    REAL(dp), INTENT(OUT) :: ansy,ansy1,ansy2
    INTEGER :: i
    REAL(dp) :: t,u
    REAL(dp), DIMENSION(4,4) :: c
    
    ! Get the c’s.
    call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
    if (x1u == x1l .or. x2u == x2l) call &
         nrerror('bcuint: problem with input values - boundary pair equal?')
    t=(x1-x1l)/(x1u-x1l)
    u=(x2-x2l)/(x2u-x2l)
    ansy=0.0_dp
    ansy2=0.0_dp
    ansy1=0.0_dp
    do i=4,1,-1
       ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
       ansy2=t*ansy2+(3.0_dp*c(i,4)*u+2.0_dp*c(i,3))*u+c(i,2)
       ansy1=u*ansy1+(3.0_dp*c(4,i)*t+2.0_dp*c(3,i))*t+c(2,i)
    end do
    ansy1=ansy1/(x1u-x1l)
    ansy2=ansy2/(x2u-x2l)
    
  END SUBROUTINE bcuint

END MODULE interp
