module nrutil
  use nrtype
  implicit none
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8

  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE

  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eqn!,assert_eq4
  END INTERFACE

  INTERFACE get_diag
     MODULE PROCEDURE get_diag_rv, get_diag_dv
  END INTERFACE

  INTERFACE outerdiff
     !Returns the outer difference of two vectors
     MODULE PROCEDURE outerdiff_r,outerdiff_d,outerdiff_i
  END INTERFACE

  INTERFACE swap
     MODULE PROCEDURE swap_rv,swap_dv
  END INTERFACE

  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE


contains
  FUNCTION arth_r(first,increment,n)
    REAL(SP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(SP), DIMENSION(n) :: arth_r
    INTEGER(I4B) :: k,k2
    REAL(SP) :: temp
    if (n > 0) arth_r(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_r(k)=arth_r(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_r(k)=arth_r(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_r
  !BL
  FUNCTION arth_d(first,increment,n)
    REAL(DP), INTENT(IN) :: first,increment
    INTEGER(I4B), INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: arth_d
    INTEGER(I4B) :: k,k2
    REAL(DP) :: temp
    if (n > 0) arth_d(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_d(k)=arth_d(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_d(k)=arth_d(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_d
  !BL
  FUNCTION arth_i(first,increment,n)
    INTEGER(I4B), INTENT(IN) :: first,increment,n
    INTEGER(I4B), DIMENSION(n) :: arth_i
    INTEGER(I4B) :: k,k2,temp
    if (n > 0) arth_i(1)=first
    if (n <= NPAR_ARTH) then
       do k=2,n
          arth_i(k)=arth_i(k-1)+increment
       end do
    else
       do k=2,NPAR2_ARTH
          arth_i(k)=arth_i(k-1)+increment
       end do
       temp=increment*NPAR2_ARTH
       k=NPAR2_ARTH
       do
          if (k >= n) exit
          k2=k+k
          arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
          temp=temp+temp
          k=k2
       end do
    end if
  END FUNCTION arth_i

  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror

  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    integer (I4B), INTENT(IN) :: n1,n2
    integer (I4B) :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2

  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    integer (I4B), INTENT(IN) :: n1,n2,n3
    integer (I4B) :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3

  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    integer (I4B), DIMENSION(:), INTENT(IN) :: nn
    integer (I4B) :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

  SUBROUTINE unit_matrix(mat)
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: mat
    INTEGER(I4B) :: i,n
    n=min(size(mat,1),size(mat,2))
    mat(:,:)=0.0_dp
    do i=1,n
       mat(i,i)=1.0_dp
    end do
  END SUBROUTINE unit_matrix

  FUNCTION get_diag_rv(mat)
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(SP), DIMENSION(size(mat,1)) :: get_diag_rv
    INTEGER(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_rv')
    do j=1,size(mat,1)
       get_diag_rv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_rv
  !BL
  FUNCTION get_diag_dv(mat)
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: mat
    REAL(DP), DIMENSION(size(mat,1)) :: get_diag_dv
    INTEGER(I4B) :: j
    j=assert_eq2(size(mat,1),size(mat,2),'get_diag_dv')
    do j=1,size(mat,1)
       get_diag_dv(j)=mat(j,j)
    end do
  END FUNCTION get_diag_dv

  FUNCTION upper_triangle(j,k,extra)
!!$When the optional argument extra is zero or absent, returnsa logical mask of
!!$shape (j, k) whose values are true above and to the rightof the diagonal, false
!!$elsewhere (including on the diagonal). When extra is present and positive,
!!$a corresponding number of additional (sub-)diagonals are returned as true.
!!$(extra = 1 makes the main diagonal return true.) When extra is present
!!$and negative, it suppresses a corresponding number of superdiagonals
    INTEGER(I4B), INTENT(IN) :: j,k
    INTEGER(I4B), OPTIONAL, INTENT(IN) :: extra
    LOGICAL(LGTP), DIMENSION(j,k) :: upper_triangle
    INTEGER(I4B) :: n
    n=0
    if (present(extra)) n=extra
    upper_triangle=(outerdiff(arth_i(1,1,j),arth_i(1,1,k)) < n)
  END FUNCTION upper_triangle

  FUNCTION outerdiff_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerdiff_r
    outerdiff_r = spread(a,dim=2,ncopies=size(b)) - &
         &spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_r
  !BL
  FUNCTION outerdiff_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerdiff_d
    outerdiff_d = spread(a,dim=2,ncopies=size(b)) - &
         & spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_d
  !BL
  FUNCTION outerdiff_i(a,b)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: a,b
    INTEGER(I4B), DIMENSION(size(a),size(b)) :: outerdiff_i
    outerdiff_i = spread(a,dim=2,ncopies=size(b)) - &
         & spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerdiff_i


  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
  
  SUBROUTINE swap_dv(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_dv

  FUNCTION outerand(a,b)
    LOGICAL(LGTP), DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL(LGTP), DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand

  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
  !BL
  FUNCTION outerprod_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d

end module nrutil
