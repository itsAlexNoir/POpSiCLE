!-----------------------------------------------------------------------------
! POpSiCLE (PhOtoelectron SpeCtrum library for Laser-matter intEractions)
!-----------------------------------------------------------------------------
!
! MODULE: tools
!> \author Alejandro de la Calle. Queen's University Belfast.
!> \author Daniel Dundas. Queen's University Belfast.
!> \date 16/04/2015
!
! DESCRIPTION:
!> \brief Auxiliary and utility routines.
!> \details This module contains subprograms needed for support
!> to other subprograms from the rest of the modules of the library.
!> Many (if not almost) all the routines in this module have been
!> extract form the numerical recipes book.
!
!------------------------------------------------------------------------------

MODULE tools

  USE constants

  IMPLICIT NONE

  PRIVATE

  !---------------------------------------------------------------
  ! Module subprograms
  !---------------------------------------------------------------

  PUBLIC  :: imaxloc
  PUBLIC  :: iminloc
  PUBLIC  :: assert_eq
  PUBLIC  :: nrerror
  PUBLIC  :: tridag

  !---------------------------------------------------------------
  ! Module interfaces
  !---------------------------------------------------------------
  
  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE imaxloc
  
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE assert_eq
  
  !> On a purely serial machine, for greater efficiency, remove
  !! the generic name tridag from the following interface,
  !! and put it on the next one after that.
  
  INTERFACE tridag
     MODULE PROCEDURE tridag_par
  END INTERFACE tridag
  
!!$  INTERFACE
!!$     SUBROUTINE tridag_ser(a,b,c,r,u)
!!$       REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,c,r
!!$       REAL(dp), DIMENSION(:), INTENT(OUT) :: u
!!$     END SUBROUTINE tridag_ser
!!$  END INTERFACE
  
     
CONTAINS

  !-----------------------------------------------------------------------------
  ! Routines returning a location as an integer value:
  !-----------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE imaxloc_r
  !
  !> \brief Index of maximun value on an array of reals.
  !
  !> \param[in] arr
  !> \return imaxloc_r index
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION imaxloc_r(arr)
    REAL(dp), INTENT(IN) :: arr(:)
    INTEGER :: imaxloc_r
    INTEGER :: imax(1)
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE imaxloc_i
  !
  !> \brief Index of maximun value on an array of integers.
  !
  !> \param[in] arr
  !> \return imaxloc_i index
  !
  !-----------------------------------------------------------------------------

  FUNCTION imaxloc_i(iarr)
    INTEGER, INTENT(IN) :: iarr(:)
    INTEGER :: imax(1)
    INTEGER :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE imaxloc_r
  !
  !> \brief Index of minimun value on an array of reals.
  !
  !> \param[in] arr
  !> \return imaxloc_r index
  !
  !-----------------------------------------------------------------------------

  
  FUNCTION iminloc(arr)
    REAL(dp), INTENT(IN) :: arr(:)
    INTEGER :: imin(1)
    INTEGER :: iminloc
    imin=minloc(arr(:))
    iminloc=imin(1)
  END FUNCTION iminloc

  !-----------------------------------------------------------------------------
  !  Routines for argument checking and error handling:
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE assert_eq2
  !
  !> \brief Report and die if integers not all equal (used for size checking).
  !
  !> \param[in] n1  Integer to check
  !> \param[in] n2 Second integer to check
  !> \param[in] string Error message
  !> \return assert_eq2
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION assert_eq2(n1,n2,string)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN)          :: n1,n2
    INTEGER                      :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
    
  END FUNCTION assert_eq2
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE assert_eq3
  !
  !> \brief Report and die if integers not all equal (used for size checking).
  !
  !> \param[in] n1  Integer to check
  !> \param[in] n2 Second integer to check
  !> \param[in] n3 Third integer to check
  !> \param[in] string Error message
  !> \return assert_eq3
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION assert_eq3(n1,n2,n3,string)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE assert_eq4
  !
  !> \brief Report and die if integers not all equal (used for size checking).
  !
  !> \param[in] n1  Integer to check
  !> \param[in] n2 Second integer to check
  !> \param[in] n3 Third integer to check
  !> \param[in] n4 Fourth integer to compare
  !> \param[in] string Error message
  !> \return assert_eq4
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE assert_eqn
  !
  !> \brief Report and die if integers not all equal (used for size checking).
  !
  !> \param[in] n  Array of integers to check
  !> \param[in] string Error message
  !> \return assert_eqn
  !
  !-----------------------------------------------------------------------------
  
  FUNCTION assert_eqn(nn,string)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: nn(:)
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn
  
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE nerror
  !
  !> \brief  Report a message, then die.
  !
  !> \param[in] string Error message
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE nrerror(string)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE tridiag_ser
  !
  !> \brief Tridiagonal algorithm, serial version
  !> \details Solves for a vector u of size N the tridiagonal linear set given
  !> by equation (2.4.1) using a serial algorithm. Input vectors b
  !> (diagonal elements) and r (right-hand sides) have size N, while a and c
  !> (off-diagonal elements) are size N − 1.
  !
  !> \param[in] a Super off-diagonal elements
  !> \param[in] b Diagonal elements
  !> \param[in] c Sub off-diagonal elements
  !> \param[in] r Right-hand side vector
  !> \param[out] u Vector solution
  !
  !-----------------------------------------------------------------------------
  
  SUBROUTINE tridag_ser(a,b,c,r,u)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(dp), DIMENSION(:), INTENT(OUT) :: u
    REAL(dp), DIMENSION(size(b)) :: gam
    INTEGER :: n,j
    REAL    :: bet
    
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
    bet=b(1)
    if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
    ! If this happens then you should rewrite your equations as a set of order N − 1, with u2
    ! trivially eliminated.
    u(1)=r(1)/bet
    ! Decomposition and forward substitution.
    do j=2,n
       gam(j)=c(j-1)/bet 
       bet=b(j)-a(j-1)*gam(j)
       ! Algorithm fails.
       if (bet == 0.0) &
            call nrerror('tridag_ser: Error at code stage 2')
       u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do
    ! Backsubstitution.
    do j=n-1,1,-1
       u(j)=u(j)-gam(j+1)*u(j+1)
    end do
  END SUBROUTINE tridag_ser
  
  !-----------------------------------------------------------------------------
  !
  ! SUBROUTINE tridiag_par
  !
  !> \brief Tridiagonal algorithm, parallel version
  !> \details Solves for a vector u of size N the tridiagonal linear set given
  !> by equation (2.4.1) using a parallel algorithm. Input vectors b
  !> (diagonal elements) and r (right-hand sides) have size N , while a and c
  !> (off-diagonal elements) are size N − 1.
  !
  !> \param[in] a Super off-diagonal elements
  !> \param[in] b Diagonal elements
  !> \param[in] c Sub off-diagonal elements
  !> \param[in] r Right-hand side vector
  !> \param[out] u Vector solution
  !
  !-----------------------------------------------------------------------------
  
  RECURSIVE SUBROUTINE tridag_par(a,b,c,r,u)
    IMPLICIT NONE
    REAL(dp), DIMENSION(:), INTENT(IN) :: a,b,c,r
    REAL(dp), DIMENSION(:), INTENT(OUT) :: u
    ! Determines when serial algorithm is invoked.
    INTEGER, PARAMETER :: NPAR_TRIDAG=4
    INTEGER :: n,n2,nm,nx
    REAL(dp), DIMENSION(size(b)/2) :: y,q,piva
    REAL(dp), DIMENSION(size(b)/2-1) :: x,z
    REAL(dp), DIMENSION(size(a)/2) :: pivc
    
    n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
    
    if (n < NPAR_TRIDAG) then
       call tridag_ser(a,b,c,r,u) 
    else
       if (maxval(abs(b(1:n))) == 0.0) &
            call nrerror('tridag_par: possible singular matrix')
       n2=size(y)
       nm=size(pivc)
       nx=size(x)
       ! Zero the odd a’s and even c’s, giving x,y,z,q.
       piva = a(1:n-1:2)/b(1:n-1:2)
       pivc = c(2:n-1:2)/b(3:n:2)
       y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
       q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
       if (nm < n2) then
          y(n2) = b(n)-piva(n2)*c(n-1)
          q(n2) = r(n)-piva(n2)*r(n-1)
       endif
       x = -piva(2:n2)*a(2:n-2:2)
       z = -pivc(1:nx)*c(3:n-1:2)
       ! Recurse and get even u's
       call tridag_par(x,y,z,q,u(2:n:2))
       ! Substitute and get odd u's.
       u(1) = (r(1)-c(1)*u(2))/b(1)
       u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
            -c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
       if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
    end if
  END SUBROUTINE tridag_par
    
END MODULE tools
