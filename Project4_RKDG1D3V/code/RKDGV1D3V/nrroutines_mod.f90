module nrroutines_mod
  implicit none
contains

  SUBROUTINE jacobi(a,d,v,nrot)
    USE nrtype; USE nrutil, ONLY : assert_eq,get_diag,nrerror,unit_matrix,&
         & upper_triangle
    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: nrot
    REAL(DP), DIMENSION(:), INTENT(OUT) :: d
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: v
!!$    Computes all eigenvalues and eigenvectors of a real symmetric NxN matrix a. On output,
!!$    elements of a above the diagonal are destroyed. d is a vector of length N that returns the
!!$    eigenvalues of a. v is an N x N matrix whose columns contain, on output, the normalized
!!$    eigenvectors of a. nrot returns the number of Jacobi rotations that were required.
    INTEGER(I4B) :: i,ip,iq,n
    REAL(DP) :: c,g,h,s,sm,t,tau,theta,tresh
    REAL(DP), DIMENSION(size(d)) :: b,z
    n=assert_eq((/size(a,1),size(a,2),size(d),size(v,1),size(v,2)/),'jacobi')
    call unit_matrix(v(:,:)) !Initialize v to the identity matrix.
    b(:)=get_diag(a(:,:)) !Initialize b and d to the diagonal of a
    d(:)=b(:) 
    z(:)=0.0 !This vector will accumulate terms of
    nrot=0 !the form tapq as in eq. (11.1.14).
    do i=1,50
       sm=sum(abs(a),mask=upper_triangle(n,n)) !Sum off-diagonal elements.
       if (sm == 0.0) RETURN
       !       The normal return, which relies on quadratic convergence to machine underflow.
       tresh=merge(0.2_dp*sm/n**2,0.0_dp, i < 4 )
       !       On the first three sweeps, we will rotate only if tresh exceeded.
       do ip=1,n-1
          do iq=ip+1,n
             g=100.0_sp*abs(a(ip,iq))
             !             After four sweeps, skip the rotation if the off-diagonal element is small.
             if ((i > 4) .and. (abs(d(ip))+g == abs(d(ip))) &
                  &.and. (abs(d(iq))+g == abs(d(iq)))) then
                a(ip,iq)=0.0
             else if (abs(a(ip,iq)) > tresh) then
                h=d(iq)-d(ip)
                if (abs(h)+g == abs(h)) then
                   t=a(ip,iq)/h !t = 1/(2Î¸)
                else
                   theta=0.5_sp*h/a(ip,iq) !Equation (11.1.10).
                   t=1.0_sp/(abs(theta)+sqrt(1.0_sp+theta**2))
                   if (theta < 0.0) t=-t
                end if
                c=1.0_sp/sqrt(1+t**2)
                s=t*c
                tau=s/(1.0_sp+c)
                h=t*a(ip,iq)
                z(ip)=z(ip)-h
                z(iq)=z(iq)+h
                d(ip)=d(ip)-h
                d(iq)=d(iq)+h
                a(ip,iq)=0.0
                call jrotate(a(1:ip-1,ip),a(1:ip-1,iq))
                !                Case of rotations 1 <= j < p.
                call jrotate(a(ip,ip+1:iq-1),a(ip+1:iq-1,iq))
                !                Case of rotations p < j < q.
                call jrotate(a(ip,iq+1:n),a(iq,iq+1:n))
                !                Case of rotations q < j <= n.
                call jrotate(v(:,ip),v(:,iq))
                nrot=nrot+1
             end if
          end do
       end do
       b(:)=b(:)+z(:)
       d(:)=b(:) !Update d with the sum of tapq,
       z(:)=0.0 !and reinitialize z.
    end do
    call nrerror('too many iterations in jacobi')
  CONTAINS
    SUBROUTINE jrotate(a1,a2)
      REAL(DP), DIMENSION(:), INTENT(INOUT) :: a1,a2
      REAL(DP), DIMENSION(size(a1)) :: wk1
      wk1(:)=a1(:)
      a1(:)=a1(:)-s*(a2(:)+a1(:)*tau)
      a2(:)=a2(:)+s*(wk1(:)-a2(:)*tau)
    END SUBROUTINE jrotate
  END SUBROUTINE jacobi


  SUBROUTINE gaussj(a,b)
    USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerand,outerprod,swap
    IMPLICIT NONE
    REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a,b
!!$    Linear equation solution by Gauss-Jordan elimination, equation (2.1.1). a is an NÃN input
!!$    coefficient matrix. b is an N ÃM input matrix containing M right-hand-side vectors. On
!!$    output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of
!!$    solution vectors.
    INTEGER(I4B), DIMENSION(size(a,1)) :: ipiv,indxr,indxc
!!$    These arrays are used for bookkeeping on the pivoting.
    LOGICAL(LGTP), DIMENSION(size(a,1)) :: lpiv
    REAL(DP) :: pivinv
    REAL(DP), DIMENSION(size(a,1)) :: dumc
    INTEGER(I4B), TARGET :: irc(2)
    INTEGER(I4B) :: i,l,n
    INTEGER(I4B), POINTER :: irow,icol
    n=assert_eq(size(a,1),size(a,2),size(b,1),'gaussj')
    irow => irc(1)
    icol => irc(2)
    ipiv=0
    do i=1,n !Main loop over columns to be reduced.
       lpiv = (ipiv == 0) !Begin search for a pivot element.
       irc=maxloc(abs(a),outerand(lpiv,lpiv))
       ipiv(icol)=ipiv(icol)+1
       if (ipiv(icol) > 1) call nrerror('gaussj: singular matrix (1)')
!!$       We now have the pivot element, so we interchange rows, if needed, to put the pivot
!!$       element on the diagonal. The columns are not physically interchanged, only relabeled:
!!$       indxc(i), the column of the ith pivot element, is the ith column that is reduced, while
!!$       indxr(i) is the row in which that pivot element was originally located. If indxr(i) =
!!$       indxc(i) there is an implied column interchange. With this form of bookkeeping, the
!!$       solution bâs will end up in the correct order, and the inverse matrix will be scrambled
!!$       by columns.
       if (irow /= icol) then
          call swap(a(irow,:),a(icol,:))
          call swap(b(irow,:),b(icol,:))
       end if
       indxr(i)=irow !We are now ready to divide the pivot row by the pivot
       indxc(i)=icol !element, located at irow and icol.
       if (a(icol,icol) == 0.0) &
            & call nrerror('gaussj: singular matrix (2)')
       pivinv=1.0_sp/a(icol,icol)
       a(icol,icol)=1.0
       a(icol,:)=a(icol,:)*pivinv
       b(icol,:)=b(icol,:)*pivinv
       dumc=a(:,icol) !Next, we reduce the rows, except for the pivot one, of
       a(:,icol)=0.0 !course.
       a(icol,icol)=pivinv
       a(1:icol-1,:)=a(1:icol-1,:)-outerprod(dumc(1:icol-1),a(icol,:))
       b(1:icol-1,:)=b(1:icol-1,:)-outerprod(dumc(1:icol-1),b(icol,:))
       a(icol+1:,:)=a(icol+1:,:)-outerprod(dumc(icol+1:),a(icol,:))
       b(icol+1:,:)=b(icol+1:,:)-outerprod(dumc(icol+1:),b(icol,:))
    end do
!!$    It only remains to unscramble the solution in view of the column interchanges. We do this
!!$    by interchanging pairs of columns in the reverse order that the permutation was built up.
    do l=n,1,-1
       call swap(a(:,indxr(l)),a(:,indxc(l)))
    end do
  END SUBROUTINE gaussj

end module nrroutines_mod
