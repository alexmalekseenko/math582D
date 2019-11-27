module gaussquad
  !
  ! This module is a revision of the Fortran 77 code gaussq.f
  !
  !   Original version 20 Jan 1975 from Stanford
  !   Modified 21 Dec 1983 by Eric Grosse
  !   f95 version Nov 2005 by Bill McLean
  !
  ! The routines are used to compute the nodes t(j) and weights
  ! w(j) for Gaussian-type quadrature rules with pre-assigned
  ! nodes.  These quantities are used when one wishes to approximate
  !    
  !                 / b
  !                 |
  !                 | f(x) w(x) dx
  !                 |
  !                 / a
  !
  ! by             
  !
  !                  n
  !                 ---
  !                 \
  !                  |  f(t(j)) * w(j).
  !                 /
  !                 ---
  !                 j=1            
  !
  ! (Note w(x) and w(j) have no connection with each other.)
  ! Here w(x) is one of six possible non-negative weight
  ! functions (listed below), and f(x) is the
  ! function to be integrated.  Gaussian quadrature is particularly
  ! useful on infinite intervals (with appropriate weight
  ! functions), since then other techniques often fail.
  !
  ! Associated with each weight function w(x) is a set of
  ! orthogonal polynomials.  The nodes t(j) are just the zeroes
  ! of the proper n-th degree polynomial.
  !
  ! References:
  !
  !      1.  Golub, G. H., and Welsch, J. H., Calculation of gaussian
  !          quadrature rules, Mathematics of Computation 23 (april,
  !          1969), pp. 221-230.
  !      2.  Golub, G. H., Some modified matrix eigenvalue problems,
  !          Siam Review 15 (april, 1973), pp. 318-334 (section 7).
  !      3.  Stroud and Secrest, Gaussian Quadrature Formulas, Prentice-
  !          Hall, Englewood Cliffs, N.J., 1966.

  implicit none

  integer, parameter :: WP = selected_real_kind(15)
  character(len=*), parameter :: VERSION = "2.4"

  ! Use at most MAXITS QL iterations in the symmetric, tridiagonal 
  ! eigenproblem routine.
  integer, parameter :: MAXITS = 30

  real(kind=WP), parameter :: ZERO = 0, HALF = 0.5_WP, &
                              ONE = 1, TWO = 2, FOUR = 4

  intrinsic :: abs, epsilon, sqrt, sign

  !
  ! The f77 function dlgama, from the specfun library, returns the 
  ! logarithm of the gamma function.
  !
  interface 
    double precision function dlgama(x) 
      double precision x
    end function dlgama
  end interface 

  !
  ! The following named constants are used to select the type of
  ! quadrature rule.
  !
  integer, parameter ::   &
    LEGENDRE         = 0, & ! w(x) = 1              on (-1,1)
    CHEBYSHEV_FIRST  = 1, & ! w(x) = 1/sqrt(1-x**2) on (-1,1)
    CHEBYSHEV_SECOND = 2, & ! w(x) = sqrt(1-x**2)   on (-1,1)
    JACOBI           = 3, & ! w(x) = (1-x)**alpha * (1+x)**beta on (-1,1)
    LAGUERRE         = 4, & ! w(x) = exp(-x) * x**alpha on (0,Inf)
    HERMITE          = 5    ! w(x) = exp(-x**2)     on (-Inf,Inf)

  character(len=*), parameter :: RULE_NAME(0:5) = (/   &
     "Legendre               ",                        &
     "Chebyshev (first kind) ",                        &
     "Chebyshev (second kind)",                        &
     "Jacobi                 ",                        &
     "Laguerre               ",                        &
     "Hermite                " /)

contains

  subroutine gauss_rule(icode, n, t, w, work, alpha, beta, endpt, info)
    !
    ! Computes an n-point Gauss quadrature rule.
    !
    integer,          intent(in)  :: icode, n
    real(kind=WP),    intent(out) :: t(n), w(n), work(n)
    real(kind=WP),    intent(in)  :: alpha, beta
    character(len=1), intent(in)  :: endpt
    integer,          intent(out) :: info
    !
    ! The arguments have the following meanings:
    !
    !   icode   The integer code specifying the weight function and interval,
    !           e.g.,  HERMITE.
    !
    !   n       The number of quadrature points.
    !
    !   t, w    Arrays holding the points and weights.
    !
    !   work    A workspace array
    !
    !   alpha   Ignored unless icode = JACOBI or LAGUERRE; defaults to 0.
    !   beta    Ignored unless icode = JACOBI; defaults to 0.
    !
    !   endpt   Specifies whether to compute a plain Gauss rule, in
    !           which all points are strictly in the interior of the
    !           interval, or a Gauss-Radau rule, in which one endpoint
    !           is a quadrature point, or a Gauss-Lobatto rule, in which 
    !           both endpoints are quadrature points.  Put 
    !           endpt = 'N' for no endpoints        (plain Gauss rule)
    !           endpt = 'B' for both endpoints      (Lobatto rule)
    !           endpt = 'L' for left endpoint only  (left Radau rule)
    !           endpt = 'R' for right endpoint only (right Radau rule)
    !           (The endpoint in question must be finite.)
    !
    !   info    = 0 on successful return
    !           = -k if kth argument has an illegal value
    !           > 0 eigensystem routine failed
    !
    real(kind=WP) :: muzero, lo, hi

    info = 0
    if ( n < 1 ) then
      info = -2
      return
    end if
    if ( icode == JACOBI .or. icode == LAGUERRE ) then
      if ( alpha <= -ONE ) then
        info = -6
        return
      end if
    end if
    if ( icode == JACOBI .and. beta <= -ONE ) then
      info = -7
      return
    end if
    if (       endpt /= 'N' .and. endpt /= 'B'        &
         .and. endpt /= 'L' .and. endpt /= 'R' ) then
      print *, endpt
      info = -8
      return
    end if

    call recur_coeffs(icode, n, alpha, beta, t, work, muzero, info)
    if ( info /= 0 ) return
    lo = -ONE
    hi = +ONE
    if ( n == 1 .and. endpt == 'B' ) then  ! Need at least 2 points to 
      info = -8                            ! include both ends.
      return
    end if
    select case(icode)
    case(LAGUERRE)
      select case(endpt)
      case('L')
        lo = ZERO
      case('R','B')
        info = -8
        return
      end select
    case(HERMITE)
      select case(endpt)
      case('L','R','B')
        info = -8
        return
      end select
    end select
    call custom_gauss_rule(n, t, work, w, muzero, endpt, lo, hi, info)
      
  end subroutine gauss_rule
    
  subroutine custom_gauss_rule(n, a, b, w, muzero, endpt, lo, hi, info)
    integer,          intent(in)    :: n
    real(kind=WP),    intent(inout) :: a(n), b(n), w(n)
    character(len=1), intent(in)    :: endpt
    real(kind=WP),    intent(in)    :: muzero, lo, hi
    integer,          intent(out)   :: info
    !
    ! On entry:
    !
    !   a, b hold the coefficients in the 3-term recurrence relation for
    !        the orthonormal polynomials (as given by recur_coeffs for
    !        any of the classical weight functions).
    !
    ! On return:
    !
    !   a, w hold the points and weights of the Gauss rule.
    !
    ! The meanings of the other arguments are as for gauss_rule.
    !
    real(kind=WP) :: g, t1
    integer       :: i

    ! The matrix of coefficients is assumed to be symmetric.  The array t 
    ! contains the diagonal elements, the array work the off-diagonal elements.
    !
    ! Make appropriate changes in the lower right 2 by 2 submatrix.
    !
    select case(endpt)
    case('L')
      if ( n == 1 ) then
        a(1) = lo
      else
        a(n) = solve(n, lo, a, b) * b(n-1)**2 + lo
      end if
    case('R')
      if ( n == 1 ) then
        a(1) = hi 
      else
        a(n) = solve(n, hi, a, b) * b(n-1)**2 + hi
      end if
    case('B')
      if ( n == 1 ) then
        info = -6
        return
      end if
      g = solve(n, lo, a, b)
      t1 = ( lo - hi ) / ( solve(n, hi, a, b) - g )
      b(n-1) = sqrt(t1)
      a(n) = lo + g * t1
    end select
    
    call st_eigenproblem(n, a, b, w, info)
    do i = 1, n
      w(i) = muzero * w(i)**2
    end do

  contains

    function solve(n, shift, a, b) result(s)
      !
      ! This procedure performs elimination to solve for the
      ! n-th component of the solution delta to the equation
      !
      !   (Jn - shift*Identity) * delta  = en,
      !
      ! where en is the vector of all zeroes except for 1 in
      ! the n-th position.  The matrix Jn is symmetric tridiagonal, with 
      ! diagonal elements a(i), off-diagonal elements b(i).  This equation
      ! must be solved to obtain the appropriate changes in the lower
      ! 2 by 2 submatrix of coefficients for orthogonal polynomials.
      !
      integer,       intent(in) :: n
      real(kind=WP), intent(in) :: shift, a(n), b(n)
      real(kind=WP)             :: s

      integer :: nm1, i
      real(kind=WP) :: t

      t = a(1) - shift
      nm1 = n - 1
      do i = 2, nm1
        t = a(i) - shift - b(i-1)**2 / t
      end do
      s = ONE / t
      
    end function solve

  end subroutine custom_gauss_rule

  subroutine recur_coeffs(icode, n, alpha, beta, a, b, muzero, info)
    integer,       intent(in)  :: icode, n
    real(kind=WP), intent(in)  :: alpha, beta
    real(kind=WP), intent(out) :: a(n), b(n), muzero
    integer,       intent(out) :: info
    !
    ! This procedure supplies the coefficients a(j), b(j) of the
    ! recurrence relation
    !
    !             b p (x) = (x - a ) p   (x) - b   p   (x)
    !              j j            j   j-1       j-1 j-2
    !
    ! For the various classical (normalized) orthogonal polynomials,
    ! and the zero-th moment
    !
    !                      / 
    !                      |
    !             muzero = | w(x) dx
    !                      |
    !                      /
    !
    ! of the given weight function w(x).  since the
    ! polynomials are orthonormalized, the tridiagonal matrix is
    ! guaranteed to be symmetric.
    !
    ! The input parameter alpha is used only for Laguerre and
    ! Jacobi polynomials, and the parameter beta is used only for
    ! Jacobi polynomials.  The Laguerre and Jacobi polynomials
    ! require the gamma function.
    !
    real(kind=WP) :: pi, ri, ab, abi, a2b2
    integer :: i

    pi = FOUR * atan(ONE)

    select case ( icode )
    case(LEGENDRE)
      muzero = TWO
      do i = 1, n
        a(i) = ZERO
        ri   = real(i, WP)
        b(i) = ri / sqrt(FOUR*ri*ri-ONE)
      end do
    case(CHEBYSHEV_FIRST)
      muzero = pi
      do i = 1, n
        a(i) = ZERO
        b(i) = HALF
      end do
      b(1) = sqrt(HALF)
    case(CHEBYSHEV_SECOND)
      muzero = HALF * pi
      do i = 1, n
        a(i) = ZERO
        b(i) = HALF
      end do
    case(HERMITE)
      muzero = sqrt(pi)
      do i = 1, n
        a(i) = ZERO
        b(i) = sqrt(HALF*i)
      end do
    case(JACOBI)
      ab = alpha + beta
      abi = TWO + ab
      muzero = TWO**(ab + ONE) * exp( &
               dlgama(alpha+ONE) + dlgama(beta+ONE) - dlgama(abi) )
      a(1) = (beta - alpha)/abi
      b(1) = sqrt( FOUR*(ONE+alpha)*(ONE+beta) / ((abi+ONE)*abi*abi) )
      a2b2 = beta*beta - alpha*alpha
      do i = 2, n
        abi = TWO*i + ab
        a(i) = a2b2 / ( (abi-TWO)*abi )
        b(i) = sqrt( FOUR*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-ONE)*abi*abi) )
      end do
      return
    case(LAGUERRE)
      muzero = exp( dlgama(alpha+ONE) )
      do i = 1, n
        a(i) = TWO*i - ONE + alpha
        b(i) = sqrt(i*(i+alpha))
      end do
    case default ! Bad value of icode.
      info = -1
  end select

  end subroutine recur_coeffs

  subroutine st_eigenproblem(n, d, e, z, info)
    ! Finds the eigenvalues and first components of the eigenvectors of a 
    ! symmetric tridiagonal matrix by the implicit QL method.
    integer,       intent(in)    :: n
    real(kind=WP), intent(inout) :: d(n), e(n)
    real(kind=WP), intent(out)   :: z(n)
    integer,       intent(out)   :: info
    !
    ! This subroutine is a translation of an algol procedure,
    ! Num. Math. 12, 377-383(1968) by Martin and Wilkinson,
    ! as modified in Num. Math. 15, 450(1970) by Dubrulle.
    ! Handbook for Auto. Comp., vol.ii-linear algebra, 241-248(1971).
    ! This is a modified version of the 'eispack' routine imtql2.
    !
    !  n     The order of the matrix.
    !
    !  d     On entry, holds the diagonal elements of the matrix.
    !        On exit,  holds the eigenvalues in ascending order.
    !        If info > 0 on exit, then the eigenvalues are correct but
    !        unordered for indices 1, 2, ..., info-1;

    !
    !  e     On enty,  holds the subdiagonal elements of the matrix
    !        in its first n-1 positions;  e(n) is arbitrary.
    !        On exit, e is overwritten.
    !
    !  z     Holds the first components of the orthonormal eigenvectors
    !        of the symmetric tridiagonal matrix.  if info > 0 on exit,
    !        then z contains the eigenvector components associated with 
    !        the stored eigenvalues.
    !
    !  info  Set to 0 on a normal return.  If info = j > 0 then the j-th
    !        eigenvalue has not been determined after MAXITS iterations.
    !
    integer :: i, j, k, l, m
    real(kind=WP) :: b, c, f, g, p, r, s

    info = 0

    z(1)   = ONE
    z(2:n) = ZERO
 
    if (n == 1) return ! nothing to do for 1x1 matrix.

    e(n) = ZERO
    loop_over_l: do l = 1, n
      loop_over_j: do j = 1, MAXITS
        ! Look for small sub-diagonal element 
        do m = l, n
          if (m == n) exit
          if ( abs(e(m)) <= epsilon(e) * (abs(d(m)) + abs(d(m+1)))) exit
        end do

        p = d(l)
        if (m == l) cycle loop_over_l
        if (j == MAXITS) then
          ! Set error -- no convergence to an eigenvalue after 30 iterations.
          info = l
          exit loop_over_l
        end if
        ! Form shift 
        g = (d(l+1) - p) / (TWO * e(l))
        r = sqrt(g*g+ONE)
        g = d(m) - p + e(l) / (g + sign(r, g))
        s = ONE
        c = ONE
        p = ZERO
        loop_over_i: do i = m-1, l, -1
          f = s * e(i)
          b = c * e(i)
          if ( abs(f) <  abs(g) ) then
            s = f / g
            r = dsqrt(s*s+ONE)
            e(i+1) = g * r
            c = ONE / r
            s = s * c
          else
            c = g / f
            r = sqrt(c*c+ONE)
            e(i+1) = f * r
            s = ONE / r
            c = c * s
          end if
          g = d(i+1) - p
          r = (d(i) - g) * s + TWO * c * b
          p = s * r
          d(i+1) = g + p
          g = c * r - b
          ! Form first component of vector.
          f = z(i+1)
          z(i+1) = s * z(i) + c * f
          z(i) = c * z(i) - s * f
        end do loop_over_i

        d(l) = d(l) - p
        e(l) = g
        e(m) = ZERO
      end do loop_over_j
    end do loop_over_l

    ! Order eigenvalues and eigenvectors.
    do i = 1, n-1
      ! Find p = d(k) = min d(j), i <= j <= n
      k = i
      p = d(i)
      do j = i+1, n
        if ( d(j) >= p ) cycle
        k = j
        p = d(j)
      end do 
      if ( k == i ) cycle
      ! Swap d(i) and d(k), z(i) and z(k)
      d(k) = d(i)
      d(i) = p
      p = z(i)
      z(i) = z(k)
      z(k) = p
    end do

  end subroutine st_eigenproblem

  subroutine orthonormal_polynomials(n, m, x, a, b, muzero, p, ldp)
    integer,       intent(in)  :: n, m, ldp
    real(kind=WP), intent(in)  :: x(m), a(n), b(n), muzero
    real(kind=WP), intent(out) :: p(ldp,0:n)
    !
    ! Returns p(i,j) = value of p_j at x(i), where p_j is the
    ! orthonormal polynomial of degree j determined by the recursion
    ! coefficients a(j), b(j) and zeroth moment muzeror, as computed with 
    ! the subroutine recur_coeffs.  Must have m <= ldp.
    !
    integer :: i, j
    real(kind=WP) :: c, rb

    c = ONE / sqrt(muzero)
    rb = ONE / b(1)
    do i = 1, m
      p(i,0) = c
      p(i,1) = rb * ( x(i) - a(1) ) * c
    end do
    do j = 2, n
      rb = ONE / b(j)
      do i = 1, m
        p(i,j) = rb * ( ( x(i) - a(j) ) * p(i,j-1) - b(j-1) * p(i,j-2) )
      end do
    end do
    
  end subroutine orthonormal_polynomials

end module gaussquad
