!
!  poly_tool_mod.f90
!
! This module contains procedures and functions to work with polynomials
! The operations supported
!  * Horner's method for evauating polynomials (on vectors) 
!  * Evaluating coefficient of the derivative of the polynomial 
!  *
! 06/25/08  Alex Alekseenko
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module poly_tool_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

   interface EvalHorner
     module procedure EvalHorner_x_scalar, EvalHorner_x_vector, EvalHorner_x_vector_coeff_matrix
   end interface

contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function evaluates the coefficients of the derivative of the polynomial 
! It expects the coefficents in the following format 
! For a polynomial 
! P_{n}(x)=a_{n}x^{n}+a_{n-1}x^{n-1}+\ldots+a_{0},
! the coefficients are in the format
! [a_{n} \ a_{n-1}\ \ldots \ a_{0}].
! 
! EXPECTS: 
! coeff -- 1x(l+1) array of the coefficients of the polynomial
!          in format [a_{n}, a_{n-1}, ... , a_{0}]
! RETURNS: 
! n --- degree of the derivative of the polynomial 
! diff_coeff -- coefficients of the derivative of the polynomial 
!            in format [b_{n-1}, b_{n-2}, ... , b_{0}]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DerivPoly (coeff) result (diff_coeff)
   intrinsic max
   real (DP), dimension (:), intent (in)  :: coeff ! coefficients of a polynomial: assumed shape dummy array has lower bound 1!
   integer                        :: n     ! degree of the polynomial
   integer                        :: loc_i ! local counter   !
   integer                        :: loc_alloc_stat ! local variable
   real (DP), dimension (1: max(size(coeff)-1,1)) :: diff_coeff  ! coefficients of the differentiated polynomial  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
n = DegreePoly(coeff);
if (n>0) then  
  diff_coeff = coeff(1:n)*(/ (loc_i, loc_i=n,1,-1) /)
else
  diff_coeff = (/ 0 /)
end if 
end function DerivPoly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function evaluates the degree of a polynomial 
! It expects the coefficents in the following format 
! For a polynomial 
! P_{n}(x)=a_{n}x^{n}+a_{n-1}x^{n-1}+\ldots+a_{0},
! the coefficients are in the format
! [a_{n} \ a_{n-1}\ \ldots \ a_{0}].
! 
! EXPECTS: 
! coeff -- 1x(l+1) array of the coefficients of the polynomial
!          in format [a_{n}, a_{n-1}, ... , a_{0}]
! RETURNS: 
! n --- degree of the derivative of the polynomial 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DegreePoly (coeff) result (n)
   real (DP), dimension (:), intent (in)  :: coeff ! coefficients of a polynomial: assumed shape dummy array has lower bound 1!
   integer  :: n                ! degree of the polynomial after differentiation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
n = size(coeff,1)-1;
end function DegreePoly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This function evaluates the polynomial using the Horner method
! It expects the coefficents in the following format 
! For a polynomial 
! P_{n}(x)=a_{n}x^{n}+a_{n-1}x^{n-1}+\ldots+a_{0},
! the coefficients are in the format
! [a_{n} \ a_{n-1}\ \ldots \ a_{0}].
! 
! EXPECTS: 
! coeff -- 1x(l+1) array of the coefficients of the polynomial
!          in format [a_{n}, a_{n-1}, ... , a_{0}]
! x -- an array of points where the polynomial P(x) needs to be evaluated
! RETURNS: 
! y -- array of values P(x)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function
!!!interface EvalHorner
!!!  module function EvalHorner_x_scalar, EvalHorner_x_vector
!!!end interface

! This is a copy of the Horner's method for the vector  variable x
function EvalHorner_x_vector (coeff, x) result (y)
   real (DP), dimension (:), intent (in) :: coeff ! coefficients of a polynomial
   real (DP), dimension (:), intent (in) :: x     ! array of point where the polynomial needs to be evaluated  
   real (DP), dimension (size(x)) :: y ! values of the polynomial 
   integer   :: loc_n, loc_k ! loc_n--degree of the polynomial, loc_k -- a counter 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
loc_n = DegreePoly(coeff)
y=0
if (loc_n > -1) then  ! in case an non-empty array came -- do the Horner's method 
   do loc_k = 1, loc_n+1 
      y = coeff(loc_k) + y*x  
   end do    
else 
     print *, "EvalHorner_x_vector: Empty array of coefficients"
end if 
end function EvalHorner_x_vector

! This is a copy of the Horner's method for the scalar variable x
function EvalHorner_x_scalar (coeff, x) result (y)
   real (DP), dimension (:), intent (in) :: coeff ! coefficients of a polynomial
   real (DP), intent (in) :: x     ! the point where the polynomial needs to be evaluated  
   real (DP) :: y ! the value of the polynomial 
   integer   :: loc_n, loc_k ! loc_n--degree of the polynomial, loc_k -- a counter 
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
loc_n = DegreePoly(coeff)
y=0
if (loc_n > -1) then  ! in case an non-empty array came -- do the Horner's method 
   do loc_k = 1, loc_n+1 
      y = coeff(loc_k) + y*x  
   end do    
else 
     print *, "EvalHorner_x_scalar: Empty array of coefficients"
end if 
end function EvalHorner_x_scalar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a copy of the horner's method that will evaluate an array of polynomials
! whose coefficients are suplied as (columns? of) coeff(:,i)
! The function will evaluate each of these polynomials on the vector of points: 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function EvalHorner_x_vector_coeff_matrix (coeff, x) result (y)
   real (DP), dimension (:,:), intent (in) :: coeff ! coefficients of a polynomial (i-th polynomial coefficiens are coeff(:,i))
   real (DP), dimension (:), intent (in) :: x     ! array of point where the polynomials needs to be evaluated  
   real (DP), dimension (size(x),size(coeff,2)) :: y ! values of the polynomials 
   integer   :: loc_k  ! -- a local counter 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (size(coeff,2) > 0) then  ! in case an non-empty array came -- do the Horner's method 
   do loc_k = 1, size(coeff,2) 
      y(:,loc_k) = EvalHorner_x_vector(coeff(:,loc_k),x)  
   end do    
else 
     print *, "EvalHorner_x_vector_coeff_matrix: Empty array of coefficients"
end if 
end function EvalHorner_x_vector_coeff_matrix

end module poly_tool_mod

