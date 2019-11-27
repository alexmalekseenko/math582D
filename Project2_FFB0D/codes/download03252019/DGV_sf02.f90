!! THis module makes a bimodal distribution using two maxwellians 
!! one maxwellian on the left and another on the right and some sort of linear? switching in between... 
!! IN 3D in velocity space and 1D in physical space...

!DIMENSIONLESS SEE NOTES FOR THE REDUCTION FORMULAS
module DGV_sf02
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   
   implicit none
   intrinsic Exp,Min
   
   interface fade
     module procedure fade, fade_x_vector
   end interface

   interface f_1D3D
     module procedure f_1D3D, f_1D3D_x_vector, f_1D3D_u_vectors, f_1D3D_x_vector_u_vectors
   end interface

   interface maxwelveldist
     module procedure maxwelveldist, maxwelveldist_T_vector, maxwelveldist_u_vectors, &
                                     maxwelveldist_T_vector_u_vectors
   end interface
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This function mimics the shock wave stationary solution. We use it here to test accuracy of 
!  integration in the initial data. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are the paramters of the gas for the shock wave Mach 3.0!
! Gas constants for Ar:  
! mass = mol_mass! molecular mass
real (DP), parameter, private :: mass=6.63d-26 
! conditions before the shock wave,
real (DP), parameter, private :: T1=0.07500_DP  ! dimensionless temperature  -- ratio to T_{\infyt}, 
real (DP), parameter, private :: d1=1.0001_DP*2.0_DP  ! dimensionless density == ration to initial number density
real (DP), parameter, private :: u1=0.7500014236_DP ! dimensionless bulk velocity  -- ratio to C_{\infty}
real (DP), parameter, private :: v1=0.0_DP ! v1=?
real (DP), parameter, private :: w1=0.0_DP ! w1=?

! condtions after the shock wave 
real (DP), parameter, private :: T2=0.27500_DP   !  dimensionless temperature  -- ratio to T_{\infyt}, 
real (DP), parameter, private :: d2=2.9999_DP*2.0_DP   !  dimensionless density == ration to initial number density
real (DP), parameter, private :: u2=0.2499978913_DP ! dimensionless bulk velocity  -- ratio to C_{\infty}
real (DP), parameter, private :: v2=0.0_DP ! v2=?
real (DP), parameter, private :: w2=0.0_DP ! w2=?
! 
real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0
!! Constants for switching ::
real (DP), parameter, private :: x1 =-0.1_DP   ! first function "begin fade" (before fade coefficient of first maxwellian is 1)
real (DP), parameter, private :: x2 = 0.1_DP   ! first function "end fade"   (after fade coefficient of first maxwellian is 0)
!!! 
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! f_1D3D            !!! to use with the main variable f
!
! This function evaluates a bimodal maxwellian distribution 
! 
! This function evaluates the 1D3D source function on vector value for x and vector values for u,v,w
! 
!
! It expects 1+3 vectors: one in x, and three the same size for u,v,w. And a scalar time t
! 
! EXPECTS: 
!  x -- point(s) in the first variable where the function needs to be evaluated
!  u,v,w -- point(s) in the second variable where the function needs to be evaluated  
!  t -- time is the third variable of the function
! RETURNS: 
! y -- value(s), usually a matrix y(size(x), size(u))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface f_1D3D
!!    module procedure f_1D3D, f_1D3D_x_vector, f_1D3D_u_vectors, f_1D3D_x_vector_u_vectors
!!  end interface
!! 


! This is a copy of the function for the vector first variable and vector second 
function f_1D3D_x_vector_u_vectors (x, u, v, w, t)  result (y)
   real (DP), dimension (:), intent (in)  :: x     ! vector of values of variable (x) where the function needs to be evaluated  
   real (DP), dimension (:), intent (in)  :: u,v,w     ! vectors of values of variable (u,v,w) where the function needs to be evaluated  
   real (DP), intent (in)  :: t        ! the value of time parameter,
   real (DP), dimension (size(x),size(u)) :: y ! values of the function  
   integer (I4B) :: i ! local counter 
   real (DP), dimension (size(x)) :: N1 ! the weight function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N1 = fade(x)
do i = 1, size(u)
y(:,i) = maxwelveldist (T1,u1,v1,w1,d1,u(i),v(i),w(i))*N1 + maxwelveldist (T2,u2,v2,w2,d2,u(i),v(i),w(i))*(1-N1)
end do  
end function f_1D3D_x_vector_u_vectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! f_1D3D            !!! to use with the main variable f
!
! This function evaluates a bimodal maxwellian distribution 
! 
! This function evaluates the 1D3D source function on vector value for x and vector values for u,v,w
! 
!
! It expects 1 scalar x and 3 vectors of same size for u,v,w. And a scalar time t
! 
! EXPECTS: 
!  x -- point(s) in the first variable where the function needs to be evaluated
!  u,v,w -- point(s) in the second variable where the function needs to be evaluated  
!  t -- time is the third variable of the function
! RETURNS: 
! y -- value(s), usually a matrix y(size(x), size(u))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface f_1D3D
!!    module procedure f_1D3D, f_1D3D_x_vector, f_1D3D_u_vectors, f_1D3D_x_vector_u_vectors
!!  end interface
!! 


! This is a copy of the function for the scalar first variable and vector second 
function f_1D3D_u_vectors (x, u, v, w, t)  result (y)
   real (DP), intent (in)  :: x     ! value of variable (x) where the function needs to be evaluated  
   real (DP), dimension (:), intent (in)  :: u,v,w     ! vectors of values of variable (u,v,w) where the function needs to be evaluated  
   real (DP), intent (in)  :: t        ! the value of time parameter,
   real (DP), dimension (size(u)) :: y ! values of the function  
   real (DP) :: N1 ! the weight function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N1 = fade(x)
y = maxwelveldist (T1,u1,v1,w1,d1,u,v,w)*N1 + maxwelveldist (T2,u2,v2,w2,d2,u,v,w)*(1-N1)
end function f_1D3D_u_vectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! f_1D3D            !!! to use with the main variable f
!
! This function evaluates a bimodal maxwellian distribution 
! 
! This function evaluates the 1D3D source function on vector value for x and vector values for u,v,w
! 
!
! It expects 1 vector for x and three scalars for u,v,w. And a scalar time t
! 
! EXPECTS: 
!  x -- point(s) in the first variable where the function needs to be evaluated
!  u,v,w -- point(s) in the second variable where the function needs to be evaluated  
!  t -- time is the third variable of the function
! RETURNS: 
! y -- value(s), usually a matrix y(size(x), size(u))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface f_1D3D
!!    module procedure f_1D3D, f_1D3D_x_vector, f_1D3D_u_vectors, f_1D3D_x_vector_u_vectors
!!  end interface
!! 


! This is a copy of the function for the vector first variable and scalar second
function f_1D3D_x_vector (x, u, v, w, t)  result (y)
   real (DP), dimension (:), intent (in)  :: x     ! vector of values of variable (x) where the function needs to be evaluated  
   real (DP), intent (in)  :: u,v,w     ! vectors of values of variable (u,v,w) where the function needs to be evaluated  
   real (DP), intent (in)  :: t        ! the value of time parameter,
   real (DP), dimension (size(x)) :: y ! values of the function  
   real (DP), dimension (size(x)) :: N1 ! the weight function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N1 = fade(x)
y = maxwelveldist (T1,u1,v1,w1,d1,u,v,w)*N1 + maxwelveldist (T2,u2,v2,w2,d2,u,v,w)*(1-N1)
end function f_1D3D_x_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! f_1D3D            !!! to use with the main variable f
!
! This function evaluates a bimodal maxwellian distribution 
! 
! This function evaluates the 1D3D source function on vector value for x and vector values for u,v,w
! 
!
! It expects 1 vector for x and three scalars for u,v,w. And a scalar time t
! 
! EXPECTS: 
!  x -- point(s) in the first variable where the function needs to be evaluated
!  u,v,w -- point(s) in the second variable where the function needs to be evaluated  
!  t -- time is the third variable of the function
! RETURNS: 
! y -- value(s), usually a matrix y(size(x), size(u))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface f_1D3D
!!    module procedure f_1D3D, f_1D3D_x_vector, f_1D3D_u_vectors, f_1D3D_x_vector_u_vectors
!!  end interface
!! 


! This is a copy of the function for the scalar first variable and scalar second
function f_1D3D (x, u, v, w, t)  result (y)
   real (DP), intent (in)  :: x     ! vector of values of variable (x) where the function needs to be evaluated  
   real (DP), intent (in)  :: u,v,w     ! vectors of values of variable (u,v,w) where the function needs to be evaluated  
   real (DP), intent (in)  :: t        ! the value of time parameter,
   real (DP)  :: y ! values of the function   
   real (DP) :: N1 ! the weight function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
N1 = fade(x)
y = maxwelveldist (T1,u1,v1,w1,d1,u,v,w)*N1 + maxwelveldist (T2,u2,v2,w2,d2,u,v,w)*(1-N1)
end function f_1D3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maxwelveldist (T,u_0,n,u) result (y)
! 
! This function evauates the 1D maxwellian equilibrium distribution with given temperature and average velocity
! Temperature and average velocity can be arrays (corresponding to different points in x variable
!
! This function evaluates 
! f2_{M}(t,x,u)=(2\pi RT(t,x))^{-1/2} \exp(-\frac{(u-\bar{u})^2}{2RT}) 
!
! This is a reloadable function
!!!!!!!!!!!!!!!!!!!!!!!
!
! this is the copy when T is vector and u is scalar
function maxwelveldist_T_vector_u_vectors (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), dimension (:), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: n    ! density parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: u_0,v_0,w_0  ! average velocity (may depend on x)
   !!! T,n,u_0 must be the same size !!! 
   real (DP), dimension (:), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP), dimension (size(T),size(u))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP), dimension (size(T))  :: beta ! local variable  to keep temporary results    
   integer (I4B) :: i ! local counter 
!!!     
beta=sqrt(pi25DT*T)*(pi25DT*T)
do i=1,size(u)
y(:,i) = n*exp(-((u(i)-u_0)*(u(i)-u_0)+(v(i)-v_0)*(v(i)-v_0)+&
                  (w(i)-w_0)*(w(i)-w_0))/max(T,0.0000001_DP))/beta
end do 
end function maxwelveldist_T_vector_u_vectors

! this is the copy when T is vector and u is scalar
function maxwelveldist_T_vector (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), dimension (:), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: n    ! density parameter (may depend on x)
   real (DP), dimension (:), intent (in) :: u_0, v_0,w_0  ! average velocity (may depend on x)
   !!! T,n,u_0 must be the same size !!! 
   real (DP), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP), dimension (size(T))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP), dimension (size(T))  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist_T_vector

! this is the copy when T is scalar and u is vector
function maxwelveldist_u_vectors (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), intent (in) :: T    ! temperature parameter (scalar)
   real (DP), intent (in) :: u_0,v_0,w_0  ! average velocity (scalar)
   real (DP), intent (in) :: n    ! density parameter (scalar)
   real (DP), dimension (:), intent (in) :: u,v,w    ! values of the velocity variable 
   real (DP), dimension (size(u))  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP)  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist_u_vectors

! this is the copy when both T and u are scalars
function maxwelveldist (T,u_0,v_0,w_0,n,u,v,w) result (y)
   real (DP), intent (in) :: T    ! temperature parameter (may depend on x)
   real (DP), intent (in) :: u_0,v_0,w_0  ! average velocity (scalar)
   real (DP), intent (in) :: n    ! density parameter (scalar)
   real (DP), intent (in) :: u,v,w    ! value of the velocity variable 
   real (DP)  :: y ! value of the dencity for this values of u and T   
!!!
   real (DP)  :: beta ! local variable  to keep temporary results    
beta=sqrt(pi25DT*T)*(pi25DT*T)
y = n*exp(-((u-u_0)*(u-u_0)+(v-v_0)*(v-v_0)+&
                  (w-w_0)*(w-w_0))/max(T,0.0000001_DP))/beta
end function maxwelveldist 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fade(x)
! 
! this functions make the fading coefficient. 
! if x <  x1 it returns 1
! if x >  x2 it returns 0
! otherwise, it returns 
!
! x1 and x2 are constants of the module.... 
!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function fade (x) result (y)
real (DP), intent (in) :: x  ! the value of the variable 
real (DP) :: y
!!!!!!!!!!!!!!!!!!
y = 1
if (x .ge. x1) then 
 if (x .ge. x2) then
 y=0 
 else 
 y=(x2-x)/(x2-x1)
 end if 
end if  
!!!!!!!!!!!!!!!!!!
end function fade 

function fade_x_vector (x) result (y)
real (DP), dimension (:), intent (in) :: x  ! the value of the variable 
real (DP), dimension (size(x)) :: y
!
integer (I4B) :: i
!!!!!!!!!!!!!!!!!!
y = 1
do i=1,size(x)
if (x(i) .ge. x1) then 
 if (x(i) .ge. x2) then
 y(i)=0 
 else 
 y(i)=(x2-x(i))/(x2-x1)
 end if 
end if  
end do
!!!!!!!!!!!!!!!!!!
end function fade_x_vector

end module DGV_sf02
