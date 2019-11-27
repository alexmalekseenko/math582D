module source_f_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none
   intrinsic Sin

   interface sf1D1D
     module procedure sf1D1D_x_scalar_u_scalar, sf1D1D_x_vector_u_scalar, sf1D1D_x_scalar_u_vector, sf1D1D_x_vector_u_vector
   end interface

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function evaluates the 1D1D source function on vector value for both arguments returns 
! the matrix with the corresponding values of a 2D function
!
! It expects two vector variables x and u  
! 
! EXPECTS: 
!  x -- point(s) in the first variable where the function needs to be evaluated
!  u -- point(s) in the second variable where the function needs to be evaluated  
! RETURNS: 
! y -- value(s), usually a matrix y(size(x), size(u))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface sf1D1D
!!    module procedure sf1D1D_x_scalar_u_scalar, sf1D1D_x_vector_u_scalar, sf1D1D_x_scalar_u_vector, sf1D1D_x_vector_u_vector
!!  end interface

! This is a copy of the function for the vector first variable and vector second 
function sf1D1D_x_vector_u_vector (x, u)  result (y)
   real (DP), dimension (:), intent (in)  :: x     ! vector of values of variable (x) where the function needs to be evaluated  
   real (DP), dimension (:), intent (in)  :: u     ! vector of values of variable (u) where the function needs to be evaluated  
   real (DP), dimension (size(x),size(u)) :: y ! values of the function  
   integer (I4B) :: x_count ! local counter 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y=0
do x_count=1,size(x)
 y(x_count,:) = sin(x(x_count) + u)
end do  
end function sf1D1D_x_vector_u_vector

! This is a copy of the function for the vector first variable and scalar second variable
function sf1D1D_x_vector_u_scalar (x, u) result (y)
   real (DP), dimension (:), intent (in)  :: x     ! vector of values in variable (x) where the function needs to be evaluated  
   real (DP), intent (in)  :: u        ! the value of variable (u) where the function needs to be evaluated  
   real (DP), dimension (size(x)) :: y ! values of the function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y = Sin(x + u)
end function sf1D1D_x_vector_u_scalar

! This is a copy of the function for the scalar first variable and vector second 
function sf1D1D_x_scalar_u_vector (x, u) result (y)
   real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
   real (DP), dimension (:), intent (in)  :: u     ! vector of values in variable (u) where the function needs to be evaluated  
   real (DP), dimension (size(u)) :: y ! values of the function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y = Sin(x+ u)
end function sf1D1D_x_scalar_u_vector

! This is a copy of the function for the scalar first and scalar second 
function sf1D1D_x_scalar_u_scalar (x, u) result (y)
   real (DP), intent (in)  :: x        ! the value of variable (x) where the function needs to be evaluated  
   real (DP), intent (in)  :: u        ! the value of variable (u) where the function needs to be evaluated  
   real (DP) :: y ! the value of the function  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
y = Sin(x + u)
end function sf1D1D_x_scalar_u_scalar

end module source_f_mod