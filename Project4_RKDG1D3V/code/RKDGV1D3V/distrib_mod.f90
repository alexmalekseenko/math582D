!
!  This module contains functions for evaluation of equilibrium distribution
!
module distrib_mod
use nrtype ! contains the kind parameters (DP), (DP), (I4B) etc. 
 implicit none
 intrinsic Exp,Min

 interface MaxwellDistr1D1D
   module procedure MaxwellDistr1D1D, MaxwellDistr1D1D_u_vector
 end interface
 
 interface ExponMaxwell1D1D
   module procedure ExponMaxwell1D1D, ExponMaxwell1D1D_u_vector
 end interface
!!!!!!!!!
real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0
!!!!!!!!!
 
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MaxwellDistr1D1D
!
! This function evaluates the Maxwellian Equilibrium Distirbution function at points where the 
! average velocity (avel) and the temperature (temp) are given  
!
! f_{M}(t,x)=n(t,x)(2\pi R T(t,x))^{-1/2}\exp[-(u-\bar{u})^2/(2RT(t,x))]
!
! It expects two scalar, vector, or matrix variables (avel) and (tempr) of the same size always
! and a scalar or vector variable u, and a scalar R (the gas constant)  
! 
! EXPECTS: 
!  dens --- values of density
!  avel --- values of the average velocity at point(s) x where the disctribution function needs to be calculated, 
!                  ! avel(m,j)  m --- the node on the interval I_{j}, j --is the interval in x from xmesh
!  tempr --- values of the temperature at point(s) x where the disctribution function needs to be calculated,
!                  ! tempr(m,j)  m --- the node on the interval I_{j}, j --is the interval in x from xmesh 
!  u -- point(s) in the velocity variable where the function needs to be evaluated (actual, not scaled) 
!  R --  is the ordinary gas constant 
!
! RETURNS: 
! y -- value of the distribution function usually a matrix y(size(u), size(avel,1))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface MaxwellDistr1D1D
!!    module procedure MaxwellDistr1D1D, MaxwellDistr1D1D_u_vector, MaxwellDistr1D1D_u_vector_avel_vector
!!  end interface

! This is a copy of the function for all scalar variables 
function MaxwellDistr1D1D (u, dens, avel, tempr, R)  result (y)
   real (DP), intent (in)  :: u     ! the value of the variable (u) where the function needs to be evaluated  
   real (DP), intent (in)  :: dens,avel, tempr  ! the value of density, average velocity and temperature at the point in x (and time) where 
									! the distribution function is being evelauted
   real (DP), intent (in)  :: R     ! the ordinary gas constant									
   real (DP) :: y ! values of the function  
   real (DP) :: beta ! local variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
beta = 2.0_DP*R*tempr
y=dens*exp(-(u-avel)*(u-avel)/max(beta,0.000000001_DP))/sqrt(pi25DT*beta)
end function MaxwellDistr1D1D
!

! This is a copy of the function for (u) vector and (avel), (temp)  scalar 
function MaxwellDistr1D1D_u_vector (u, dens, avel, tempr, R)  result (y)
   real (DP), dimension (:), intent (in)  :: u     ! the value of the variable (u) where the function needs to be evaluated  
   real (DP), intent (in)  :: dens, avel, tempr  ! the value of average velocity and temperature at the point in x (and time) where 
									! the distribution function is being evelauted
   real (DP), intent (in)  :: R     ! the ordinary gas constant									
   real (DP), dimension (size(u)) :: y ! values of the function  
   real (DP) :: beta ! local variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
beta = 2.0_DP*R*tempr
y=dens*exp(-(u-avel)*(u-avel)/max(beta,0.000000001_DP))/sqrt(pi25DT*beta)
end function MaxwellDistr1D1D_u_vector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExponMaxwell1D1D
!
! This function evaluates the exponent that is part of the 
! Maxwellian Equilibrium Distirbution function at points where the 
! average velocity (avel) and the temperature (temp) are given  
!
!  e_{M}(t,x)=\exp[-(u-\bar{u})^2/(2RT(t,x))]
!
! This function is an intermediate result to be used in integrals of the maxwellian distribution: the integrals 
! will be multiplied by the factor n(t,x)(2\pi R T(t,x))^{-1/2} on the outside
!
! It expects two scalar, vector, or matrix variables (avel) and (tempr) of the same size always
! and a scalar or vector variable u, and a scalar R (the gas constant)  
! 
! EXPECTS: 
!  avel --- values of the average velocity at point(s) x where the disctribution function needs to be calculated, 
!                  ! avel(m,j)  m --- the node on the interval I_{j}, j --is the interval in x from xmesh
!  tempr --- values of the temperature at point(s) x where the disctribution function needs to be calculated,
!                  ! tempr(m,j)  m --- the node on the interval I_{j}, j --is the interval in x from xmesh 
!  u -- point(s) in the velocity variable where the function needs to be evaluated (actual, not scaled) 
!  R --  is the ordinary gas constant 
!
! RETURNS: 
! y -- value of the distribution function usually a matrix y(size(u), size(avel,1))
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! This is a reloadble function: see the interface statement above>>
!!   interface MaxwellDistr1D1D
!!    module procedure MaxwellDistr1D1D, MaxwellDistr1D1D_u_vector, MaxwellDistr1D1D_u_vector_avel_vector
!!  end interface

! This is a copy of the function for all scalar variables 
function ExponMaxwell1D1D (u, avel, tempr, R)  result (y)
   real (DP), intent (in)  :: u     ! the value of the variable (u) where the function needs to be evaluated  
   real (DP), intent (in)  :: avel, tempr  ! the value of average velocity and temperature at the point in x (and time) where 
									! the distribution function is being evelauted
   real (DP), intent (in)  :: R     ! the ordinary gas constant									
   real (DP) :: y ! values of the function  
   real (DP) :: beta ! local variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
beta = 2.0_DP*R*tempr
y=exp(-(u-avel)*(u-avel)/max(beta,0.0000001_DP))
end function ExponMaxwell1D1D
!

! This is a copy of the function for (u) vector and (avel), (temp)  scalar 
function ExponMaxwell1D1D_u_vector (u, avel, tempr, R)  result (y)
   real (DP), dimension (:), intent (in)  :: u     ! the value of the variable (u) where the function needs to be evaluated  
   real (DP), intent (in)  :: avel, tempr  ! the value of average velocity and temperature at the point in x (and time) where 
									! the distribution function is being evelauted
   real (DP), intent (in)  :: R     ! the ordinary gas constant									
   real (DP), dimension (size(u)) :: y ! values of the function  
   real (DP) :: beta ! local variable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
beta = 2.0_DP*R*tempr
y=exp(-(u-avel)*(u-avel)/max(beta,0.0000001_DP))
end function ExponMaxwell1D1D_u_vector
!

end module distrib_mod