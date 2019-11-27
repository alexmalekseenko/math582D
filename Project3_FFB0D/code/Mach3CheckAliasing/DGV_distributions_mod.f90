!! This module makes a bimodal distribution using two maxwellians 
!! one maxwellian on the left and another on the right and some sort of linear? switching in between... 
!! IN 3D in velocity space and 1D in physical space...

!DIMENSIONLESS SEE NOTES FOR THE REDUCTION FORMULAS
module DGV_distributions_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 

implicit none
   interface maxwelveldist
     module procedure maxwelveldist, maxwelveldist_T_vector, maxwelveldist_u_vectors, &
                                     maxwelveldist_T_vector_u_vectors
   end interface

   interface ESBGK
     module procedure ESBGK_f0
   end interface
   
!!!!!!!!!!!!! parameters and global variables 
real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0
!!!!!!!!!!


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! maxwelveldist (T,u_0,n,u) result (y)
! 
! This function evauates the 3D velocity dimensionless maxwellian equilibrium distribution with given temperature and average velocity
! Temperature and average velocity can be arrays (corresponding to different points in x variable
!
! This function evaluates 
! f_{M}(t,x,u)=(\pi T(t,x))^{-1/2} \exp(-\frac{(u-\bar{u})^2}{2RT}) 
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ESBGK_f0
! 
! This function evauates the 3D ES-BGK equilibrium distribution with given temperature and average velocity.
! Note that the parameters are dependent on time and velocity but not in the spacial variable x
!
! This function evaluates 
! f_0(t,u)=\frac{n(t)}{\sqrt{(\pi^3 det(\mathbb{T})}} \exp(-c'\mathbb{T}c) 
!
!!!!!!!!!!!!!!!!!!!!!!!

function ESBGK_f0 (TensorInv,Determinant,n,u_0,v_0,w_0,nodes_u,nodes_v,nodes_w) result (f0)
! Evaluate f0 for ESBGK
real (DP), dimension(3,3), intent (in) :: TensorInv ! The inverse of the tensor matrix
real (DP), intent (in) :: u_0, v_0, w_0 ! These are the bulk velocities
real (DP), intent (in) :: Determinant ! This is the Determinant of the Tensor matrix
real (DP), intent (in) :: n ! this is the number density
real (DP), dimension(:), intent (in) :: nodes_u,nodes_v,nodes_w ! These are the velocity nodes
integer (I4B) :: Unodes ! the number of nodes in velocity space
real (DP), dimension(3) :: Left, c ! c-v vector
real (DP), parameter :: PI25DT = 3.141592653589793238462643d0
integer :: loc_alloc_stat
real (DP), dimension (:), allocatable  :: Prod ! This is the component that goes into the argument of the exponent
real (DP), dimension (1:size(nodes_u,1)) :: f0 ! result -- the ES-BGK distribution function
!real (DP), dimension(1:size(nodes_u,1)) :: nodes_cu, nodes_cv, nodes_cw
integer (I4B) :: i,j,k ! these are local counters

Unodes = size(nodes_U)   !  Get the number of elements of each of the nodes

allocate (Prod(1:size(nodes_u,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "ESBGK_f0: Allocation error for variables (Prod)"
 stop
end if

! from here, the component in the exponent's argument is computed

do i=1,Unodes
   ! c = \vec{u} - \vec{bar{u}}
   c(1) = nodes_u(i) - u_0
   c(2) = nodes_v(i) - v_0
   c(3) = nodes_w(i) - w_0

   ! used as part of the left operation (i.e. c'*TensorInv)
   Left(1) = c(1)*TensorInv(1,1) + c(2)*TensorInv(2,1) + c(3)*TensorInv(3,1)
   Left(2) = c(1)*TensorInv(1,2) + c(2)*TensorInv(2,2) + c(3)*TensorInv(3,2)
   Left(3) = c(1)*TensorInv(1,3) + c(2)*TensorInv(2,3) + c(3)*TensorInv(3,3)

   Prod(i) = Left(1)*c(1) + Left(2)*c(2) + Left(3)*c(3)! final computation of c'*TensorInv*c
end do

f0 = n/sqrt((PI25DT)**3*Determinant) * exp(-Prod)

deallocate (Prod)

end function ESBGK_f0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shakhov_f0
! 
! This function evauates the dimensionless 3D Shakhov equilibrium distribution function 
! with given temperature and average velocity and the 
! give Prandtl number. 
!
!Note that the parameters are dependent on time and velocity but not in the spacial variable x
!
! This function evaluates 
! f_0(t,u)=f_{M}(t,u)(1+(4/5)(1-Pr)S_{a}c^{a}(c_{a}c^{a}-5/2)) where c_a=(u_a-\bar{u}_a)/\sqrt{T}
! S_a=q_{a}/(nT^{3/2}), q_{a} = \int f (u_{a}-\bar(u)_{a})|u-\bar{u}|^2 du
!
!!!!!!!!!!!!!!!!!!!!!!!

function Shakhov_f0 (alpha,n,u_0,v_0,w_0,T,S_u,S_v,S_w,nodes_u,nodes_v,nodes_w) result (f0)
! Evaluate f0 for Shakhov model
real (DP), intent (in) :: alpha ! the Prandtl number 
real (DP), intent (in) :: T ! This is the temperature of the local maxwellian 
real (DP), intent (in) :: n ! this is the number density 
real (DP), intent (in) :: u_0, v_0, w_0 ! These are the bulk velocities of the local maxwellian
real (DP), dimension(:), intent (in) :: S_u,S_v,S_w ! These are the components of the vector S in the shakov model

real (DP), dimension(:), intent (in) :: nodes_u,nodes_v,nodes_w ! These are the velocity nodes
real (DP), dimension (1:size(nodes_u,1)) :: f0 ! the Shakhov distribution function -- the result
!
!
real (DP), parameter :: PI25DT = 3.141592653589793238462643d0
real (DP), dimension (:), allocatable :: ShCorr,c_u, c_v, c_w! This are useful arrays, 
integer (I4B) :: i,j,k,loc_alloc_stat ! these are local counters


allocate (ShCorr(1:size(nodes_u,1)),c_u(1:size(nodes_u,1)),c_v(1:size(nodes_u,1)),&
                  c_w(1:size(nodes_u,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "Shakhov_f0: Allocation error for variables (ShCorr,c_u,c_v,c_w)"
 stop
end if

! first, we compute the components of vector c:
c_u = (nodes_u-u_0)/sqrt(T)
c_v = (nodes_v-v_0)/sqrt(T)
c_w = (nodes_w-w_0)/sqrt(T)
! Next we compute the Shakhov correction term
ShCorr = (S_u*c_u+S_v*c_v+S_w*c_w)*((c_u*c_u+c_v*c_v+c_w*c_w)-5.0_DP/2.0_DP)
ShCorr = (1+4.0_DP/5.0_DP*((-1)*alpha)/(1-alpha)*ShCorr)
!!!!!!!!
f0 = maxwelveldist(T,u_0,v_0,w_0,n,nodes_u,nodes_v,nodes_w)*ShCorr

deallocate (ShCorr,c_u,c_v,c_w)

end function Shakhov_f0

end module DGV_distributions_mod