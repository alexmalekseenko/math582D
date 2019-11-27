! This Module deals with the functions and procedures for 
! basis functions (Legendre Polynomials)
! The operations supported: 
!  *  Generates Coefficients of Legendre's polynomials
!!!!!!!!!!! 06/24/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
module basis_fun_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none
contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LegendrePoly (n) 
! returns coefficents of the Legendre polynomial of degree n
! P_{n}(x)=a_{n}x^{n}+a_{n-1}x^{n-1}+\ldots+a_{0},
! in the format
! [a_{n} \ a_{n-1}\ \ldots \ a_{0}].
! 
! EXPECTS: 
!    n --- degree of the Legendre polynomial 
! RETURNS: 
!    coeff --- array of the coefficents 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 function LegendrePoly (n) result (coeff)
   integer, intent (in)  :: n
   real (DP), dimension (0:n) :: coeff
! for now just choose from known polynomials
   select case (n)
        case (0)
            coeff = (/ 1.0 /)
        case (1)
            coeff = (/ 1.0, 0.0 /)
        case (2)
            coeff = (/ 3.0, 0.0, -1.0 /) / 2.0
        case (3)
            coeff = (/ 5.0, 0.0, -3.0, 0.0 /) / 2.0 
        case (4)
            coeff = (/ 35.0, 0.0, -30.0, 0.0, 3.0 /) / 8.0
        case (5)
            coeff = (/ 63.0, 0.0, -70.0, 0.0, 15.0, 0.0 /) / 8.0
        case (6)
            coeff = (/ 231.0, 0.0, -315.0, 0.0, 105.0, 0.0, -5.0 /) / 16.0
		case (7)
			coeff = (/ 429.0, 0.0, -693.0, 0.0, 315.0, 0.0, -35.0, 0.0 /) / 16.0 
		case (8)
			coeff = (/ 6435.0, 0.0, -12012.0, 0.0, 6930.0, 0.0, -1260.0, 0.0, 35.0 /) / 128.0 
		case default 
		    print *, "Unsupported Legendre Polynomial degree"
		    coeff = 0
		    stop
    end select 
 end function LegendrePoly

end module basis_fun_mod
   