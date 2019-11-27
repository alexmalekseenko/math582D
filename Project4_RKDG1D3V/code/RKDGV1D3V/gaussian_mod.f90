module gaussian_mod
implicit none
contains
  SUBROUTINE GauLeg(x1,x2,x,w)
    USE nrtype
    USE nrutil, ONLY : arth,assert_eq,nrerror
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x1,x2
    REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
    REAL(DP), PARAMETER :: EPS=3.0e-14_dp
!!$Given the lower and upper limits of integration x1 and x2, this routine returns arrays x and w
!!$of length N containing the abscissas and weights of the Gauss-Legendre N-point quadrature
!!$formula. The parameter EPS is the relative precision. Note that internal computations are
!!$done in double precision.
    INTEGER(I4B) :: its,j,m,n
    INTEGER(I4B), PARAMETER :: MAXIT=40
    REAL(DP) :: xl,xm
    REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
    LOGICAL(LGTP), DIMENSION((size(x)+1)/2) :: unfinished
    n=assert_eq(size(x),size(w),'GauLeg')
    m=(n+1)/2 !The roots are symmetric in the interval,so
    xm=0.5_dp*(x2+x1) !we only have to find half of them.
    xl=0.5_dp*(x2-x1)
    z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))! Initial approximations to the roots.
    unfinished=.true.
    do its=1,MAXIT !Newtonâs method carried out simultaneously on the roots.
       where (unfinished) 
          p1=1.0
          p2=0.0
       end where
       do j=1,n !Loop up the recurrence relation to get the Legendre polynomials evaluated at z.
          where (unfinished)
             p3=p2
             p2=p1
             p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
          end where
       end do
       !p1 now contains the desired Legendre polynomials. We next compute pp, the derivatives,
       !by a standard relation involving also p2, the polynomials of one lower order.
       where (unfinished)
          pp=n*(z*p1-p2)/(z*z-1.0_dp)
          z1=z
          z=z1-p1/pp !Newtonâs method.
          unfinished=(abs(z-z1) > EPS)
       end where
       if (.not. any(unfinished)) exit
    end do
    if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
    x(1:m)=xm-xl*z !Scale the root to the desired interval,
    x(n:n-m+1:-1)=xm+xl*z !and put in its symmetric counterpart.
    w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2) !Compute the weight
    w(n:n-m+1:-1)=w(1:m) !and its symmetric counterpart.
  END SUBROUTINE GauLeg
end module gaussian_mod
