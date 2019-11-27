! 09/05/08 Alex
!
! This Module deals with the functions and procedures for 
! evaluating and handling the coefficients of spectral representation 
! of a given function 
! The operations supported: 
!  *  Calculates the coefficients of the spectral represenation by the Legendre's polynomials
!  *  Evaluates the function on a grid from spectral coefficients.
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
module spectral_tools_mod
!
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
use basis_fun_mod
use poly_tool_mod
!
   implicit none
!
interface EvCellLeg1D1D
   module procedure EvCellLeg1D1D, EvCellLeg1D1D_u_vector
end interface

interface EvCellLeg1D
   module procedure EvCellLeg1D, EvCellLeg1D_x_vector
end interface

interface EvalLegPol
   module procedure EvalLegPol, EvalLegPol_x_vector
end interface

interface EvalDerLegPol
   module procedure EvalDerLegPol, EvalDerLegPol_x_vector
end interface

interface EvCellLegStdrd1D1D 
   module procedure EvCellLegStdrd1D1D, EvCellLegStdrd1D1D_u_vector
end interface
!
 intrinsic Epsilon   ! will 
 real (DP), parameter :: eps = epsilon(1.0_DP)*4   ! the epsilon constant to use when comparing real numbers 
!
contains 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DecompLeg1D1D () 
! Returns coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. 
!     
! 
! EXPECTS:
!    xmesh -- meshpoints in variable x
!    umesh -- meshpoints in variable u
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_degree -- order of maximum allowed degree of approximation in any variable 
!    f --- is a dummy external function whose spectral coefficients need to be evaluated. 
!      (ATTN: the function DecompLeg1D1D is coded so that (f) accepts first scalar variable and second vector variable)  
! RETURNS: 
!    coeff --- array of the coefficents 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DecompLeg1D1D (xmesh,umesh,k,s,max_deg,f,curr_time) result (coeff)
  use gaussian_mod
  !!!!!!
  real (DP), dimension (:), intent (in)  :: xmesh ! mesh points in varaible x 
  real (DP), dimension (:), intent (in)  :: umesh ! mesh points in varaible u 
  integer (I4B), intent (in) :: k,s            ! k and s are the higher degrees of polynomal basis in   
  integer (I4B), intent (in) :: max_deg     ! max_degree of the 2D polynomial k and s are the higher degrees of polynomal basis in   
  real (DP), intent (in) :: curr_time ! the time variable -- to pass to the exact solution function f
  !
  real (DP), dimension (:), allocatable :: x_gnodes,x_gweights,u_gnodes,u_gweights ! Gaussian nodes and weights
  real (DP), dimension (:,:), allocatable :: x_wleg_pol,u_wleg_pol !  arrays to store values of leg_polynomials on gaussian nodes multiplied by the weights
  real (DP), dimension (:), allocatable :: x,u  ! local variable, use for in integration in x and u variables   
  integer (I4B) :: i,x_count,u_count,p,l ! local counters
  real (DP), dimension (0:k,0:s,size(umesh)-1,size(xmesh)-1) :: coeff 
  integer :: loc_alloc_stat                     ! local variable
  ! interface block for the dummy function f: (x-scalar, u-vector) 
  ! first variable must be scalar and second must be vector! 
  interface 
   function f (x, u, curr_time) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: curr_time ! the time variable 
    real (DP), dimension (size(u)) :: y ! values of the function  
   end function f
  end interface
  !
  !!!!!! Check for problems with the supplied variables 
  if (size(xmesh) < 2) then  
   print *, "DecompLeg1D1D: ERROR: mesh #1 (xmesh) must have at least two points. stop"
   stop
  end if 
  if (size(umesh) < 2) then  
   print *, "DecompLeg1D1D: ERROR: mesh #2 (xmesh) must have at least two points. stop"
   stop
  end if 
  if ((k < 0) .or. (s < 0)) then  
   print *, "DecompLeg1D1D: ERROR: negative order of polynomial approximation is supplied (k) or (s). stop"
   stop
  end if 
  if (max_deg < max(k,s)) then 
   print *, "DecompLeg1D1D: ERROR: maximal degree of 2D polynomial (max_deg) must be >= 1D degrees (k and s). stop"
   stop
  end if 
!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gnodes(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (x_gnodes)"
     end if 
     !
  allocate (x_gweights(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (x_gweights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  call GauLeg(Real(-1,DP),Real(1,DP),x_gnodes,x_gweights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
!!! Prepare the values of the Legendre polynomials for integration in the first variable 
  allocate (x_wleg_pol(0:k,size(x_gnodes)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (x_wleg_pol)"
     end if 
     !
!!! Evaluate the basis functions on the nodes
  do i=0,k
   x_wleg_pol(i,:) = EvalHorner (LegendrePoly (i), Real(x_gnodes,DP))*Real(x_gweights,DP)
  end do 
!!! End evaluation of basis functions for x variable on gauss nodes  
!!!
!!! Prepare the Gaussian nodes and weights in the second variable (u) !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the second variable
  allocate (u_gnodes(1:max(s+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (u_gnodes)"
     end if 
     !
  allocate (u_gweights(1:max(s+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (u_gweights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the second variable    
  call GauLeg(Real(-1,DP),Real(1,DP),u_gnodes,u_gweights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
!!! Prepare the values of the Legendre polynomials for integration in the first variable 
  allocate (u_wleg_pol(0:s,size(u_gnodes)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D1D: Allocation error for variable (u_wleg_pol)"
     end if 
     !
!!! Evaluate the basis functions on the nodes and muliply by the weights
  do i=0,s
   u_wleg_pol(i,:) = EvalHorner (LegendrePoly (i), Real(u_gnodes,DP))*Real(u_gweights,DP)
  end do 
!!! End evaluation of basis functions for u variable on gauss nodes  
!!!
!!! We will now calculate the coefficients of the spectra decomposition on the entire (rectangular) mesh    
!!! 
!!! Integration BIG Loop over all cells in 2D
!!! 
!!! allocate space for the inside cell nodes: 
  allocate (x(size(x_gnodes)), stat=loc_alloc_stat)
   !
   if (loc_alloc_stat >0) then 
   print *, "DecompLeg1D1D: Allocation error for variable (x)"
   end if 
  allocate (u(size(u_gnodes)), stat=loc_alloc_stat)
   !
   if (loc_alloc_stat >0) then 
   print *, "DecompLeg1D1D: Allocation error for variable (u)"
   end if 
!!! end allocate space for the inside cell nodes    
  do x_count=1,size(xmesh)-1
  do u_count=1,size(umesh)-1
  !!! prepare integration with the cell
    x = (xmesh(x_count+1)-xmesh(x_count))*x_gnodes/2+(xmesh(x_count+1)+xmesh(x_count))/2
    u = (umesh(u_count+1)-umesh(u_count))*u_gnodes/2+(umesh(u_count+1)+umesh(u_count))/2
  !!! Integration within the cell
    do p=0,k
    do l=0,min(max_deg-p,s)
      coeff(p,l,u_count,x_count)=0;
      do i=1,size(x)  
      coeff(p,l,u_count,x_count) = coeff(p,l,u_count,x_count) + sum(f(x(i),u,curr_time)*u_wleg_pol(l,:))*x_wleg_pol(p,i)
      end do 
      coeff(p,l,u_count,x_count) = (2*p+1)*(2*l+1)*coeff(p,l,u_count,x_count)/4     
    end do 
    end do 
  !!! end integration within the cell
  end do  
  end do 
!!! End Integration Big Loop  
deallocate (x,u,u_wleg_pol,u_gweights,u_gnodes,x_wleg_pol,x_gweights,x_gnodes)
end function DecompLeg1D1D 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DecompLeg1D3D ()
!
! this is a copy of the above subroutine adjusted to 3D nodal galerkin discretization in the velocity variable:
! 
! Returns coefficents of the spectral representation by scaled Legendre polynomials
! of a 1D+3D dimensional functions treating the 1D first variable and second 3D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation in x is (k). The desgree of the NodalDG discretization is 
! controlled by the DGV library. The subroutine receives the nodal values of the 3D velocity
!     
! 
! EXPECTS:
!    xmesh -- meshpoints in variable x
!    nodes_u, nodes_v, nodes_w -- meshpoints in variable u,v,w
!    k  --- order of polynomial approximation in variable x
!    f --- is a dummy external function whose spectral coefficients need to be evaluated. 
!      (ATTN: the function DecompLeg1D1D is coded so that (f) accepts first scalar variable and second vector variable)  
! RETURNS: 
!    coeff --- array of the coefficents 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DecompLeg1D3D (xmesh,nodes_u,nodes_v,nodes_w,k,f,curr_time) result (coeff)
  use gaussian_mod
  !!!!!!
  real (DP), dimension (:), intent (in)  :: xmesh ! mesh points in varaible x 
  real (DP), dimension (:), intent (in)  :: nodes_u, nodes_v, nodes_w ! mesh points in varaibles u,v,w 
  integer (I4B), intent (in) :: k            ! k and s are the higher degrees of polynomal basis in   
  real (DP), intent (in) :: curr_time ! the time variable -- to pass to the exact solution function f
  real (DP), dimension (0:k,size(nodes_u,1),size(xmesh,1)-1) :: coeff 
    !
  real (DP), dimension (:), allocatable :: x_gnodes,x_gweights ! Gaussian nodes and weights
  real (DP), dimension (:,:), allocatable :: x_wleg_pol !  arrays to store values of leg_polynomials on gaussian nodes multiplied by the weights
  real (DP), dimension (:), allocatable :: x  ! local variable, use for in integration in x and u variables   
  integer (I4B) :: i,x_count,u_count,p ! local counters
  integer :: loc_alloc_stat                     ! local variable
  ! interface block for the dummy function f: (x-scalar, u-vector) 
  ! first variable must be scalar and second must be vector! 
  interface 
   function f (x, u, v, w, curr_time)  result (y)
    ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
    !
    real (DP), intent (in)  :: x     ! value of variable (x) where the function needs to be evaluated  
    real (DP), intent (in)  :: u,v,w     ! vectors of values of variable (u,v,w) where the function needs to be evaluated  
    real (DP), intent (in)  ::  curr_time       ! the value of time parameter,
    real (DP)  :: y ! values of the function  
   end function f
  end interface
  !
  !!!!!! Check for problems with the supplied variables 
  if (size(xmesh) < 2) then  
   print *, "DecompLeg1D3D: ERROR: mesh #1 (xmesh) must have at least two points. stop"
   stop
  end if 
  if (size(nodes_u) < 2) then  
   print *, "DecompLeg1D3D: ERROR: arrays of nodes in the velocity variable must have at least two points. stop"
   stop
  end if 
  if (k < 0) then  
   print *, "DecompLeg1D3D: ERROR: negative order of polynomial approximation is supplied (k). stop"
   stop
  end if 
!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gnodes(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D3D: Allocation error for variable (x_gnodes)"
     end if 
     !
  allocate (x_gweights(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D3D: Allocation error for variable (x_gweights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  call GauLeg(Real(-1,DP),Real(1,DP),x_gnodes,x_gweights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
!!! Prepare the values of the Legendre polynomials for integration in the first variable 
  allocate (x_wleg_pol(0:k,size(x_gnodes)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "DecompLeg1D3D: Allocation error for variable (x_wleg_pol)"
     end if 
     !
!!! Evaluate the basis functions on the nodes
  do i=0,k
   x_wleg_pol(i,:) = EvalHorner (LegendrePoly (i), Real(x_gnodes,DP))*Real(x_gweights,DP)
  end do 
!!! End evaluation of basis functions for x variable on gauss nodes  
!!!
!!!
!!! We will now calculate the coefficients of the spectra decomposition on the entire (rectangular) mesh    
!!! 
!!! Integration BIG Loop over all cells in variable X
!!! 
!!! allocate space for the inside cell nodes: 
  allocate (x(size(x_gnodes)), stat=loc_alloc_stat)
   !
   if (loc_alloc_stat >0) then 
   print *, "DecompLeg1D1D: Allocation error for variable (x)"
   end if 
!!! end allocate space for the inside cell nodes    
 do x_count = 1,size(xmesh,1)-1
  do u_count = 1,size(nodes_u,1)
    !!! prepare integration with the cell
    x = (xmesh(x_count+1)-xmesh(x_count))*x_gnodes/2+(xmesh(x_count+1)+xmesh(x_count))/2
    !!! Integration within the cell
    do p=0,k
      coeff(p,u_count,x_count)=0;
      do i=1,size(x)  
      coeff(p,u_count,x_count) = coeff(p,u_count,x_count) & 
             + f(x(i),nodes_u(u_count),nodes_v(u_count),nodes_w(u_count),curr_time)*x_wleg_pol(p,i)
      end do 
      coeff(p,u_count,x_count) = (2*p+1)*coeff(p,u_count,x_count)/2.0_DP     
    end do 
    !!! end integration within the cell in x 
  end do  
 end do 
!!! End Integration Big Loop  
deallocate (x,x_wleg_pol,x_gweights,x_gnodes)
end function DecompLeg1D3D 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AssemblLeg1D1D (xxmesh,uumesh,xmesh,umesh,k,s,max_deg,coeff) result (y)
! Evaluates the 1D+1D dimensional function from its spectral representation in Legandre's scaled polynomials.
! The varaible (coeff) contains the coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. The variables (xmesh) and (umesh) determine the spectral grid
! The varaibles (xxmesh) and (uumesh) determine the mesh where the function needs to be evaluated. 
!     
! 
! EXPECTS:
!    xxmesh -- meshpoints in var x where the spectral representation needs to be evaluated 
!    uumesh -- meshpoints in var u where the spectral representation needs to be evaluated 
!    xmesh -- meshpoints of spectral decomposition in variable x
!    umesh -- meshpoints of spectral docomposition in variable u
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_degree -- order of maximum allowed degree of approximation in any variable 
!    coeff --- array of the coefficents 
! RETURNS: 
!    y  -- is a matrix of values  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function AssembleLeg1D1D(xxmesh,uumesh,xmesh,umesh,k,s,max_deg,coeff) result (y)
real (DP), dimension (:), intent (in)  :: xxmesh ! mesh in variable (x) where the spectral representation will be evaluated
real (DP), dimension (:), intent (in)  :: uumesh ! mesh in variable (u) where the spectral representation will be evaluated
real (DP), dimension (:), intent (in)  :: xmesh  ! mesh in variable (x) where the spectral decomposition is defined (ordered)
real (DP), dimension (:), intent (in)  :: umesh  ! mesh in variable (u) where the spectral decomposition is defined (ordered)
integer (I4B), intent (in) :: k,s                ! orders of approximations in (x) and (u) respectively
integer (I4B), intent (in) :: max_deg            ! degree of maximum applowed polynomial interpolation 
real (DP), dimension (0:,0:,:,:), intent (in)  :: coeff  ! coefficients of the spectral decomposition that will be evaluated 
real (DP), dimension (size(xxmesh), size (uumesh)) :: y ! the evaluated of the spectral representation for m=1:Nxx
 ! 
integer (I4B) :: Nxx,Nuu,Nx,Nu
integer (I4B) :: m,n,m_i,n_j        ! useful local counters
!! initialize some variables 
Nxx = size(xxmesh) ! useful local integer constant 
Nuu = size(uumesh) ! useful local integer constant 
Nx = size(xmesh)   ! useful local integer constant 
Nu = size(umesh)   ! useful local integer constant
!! 
 ! pick a point on the drawing mesh
do m=1,Nxx
    do n=1,Nuu
        ! check if the point is off the area on which the spectral coefficients are defined 
        ! we assume that the spectral mesh is ordered! 
        if ((xxmesh(m) <= (xmesh(1) - eps)) .or. (xxmesh(m) >= (xmesh(Nx) + eps)) .or. &  
           (uumesh(n) <= (umesh(1) - eps)) .or. (uumesh(n) >= (umesh(Nu) + eps)) ) then 
            print *, "graphing point is outside the discrete mesh" 
        else 
        ! take the selected point (xxmesh(m), uumesh(n)) and find 
        ! the corresponding cell [xmesh(m_i-1),xmesh(m_i)]\times [umesh(n_j-1),umesh(n_j)] 
        m_i=2
        ! find the interval on spectral mesh for m:
        do while ((xxmesh(m) > xmesh(m_i) + eps) .and. (m_i < Nx)) 
        m_i = m_i+1;
        end do 
        ! we found m_i such that xxmesh(m)\in [xmesh(m_i-1), xmesh(m_i)]
        ! Now we find the interval on spectal mesh for n:
        n_j=2
        do while ((uumesh(n) > umesh(n_j) + eps) .and. (n_j < Nu))
        n_j = n_j+1;
        end do 
        ! we found n_j such that uumesh(n)\in [umesh(n_j-1), umesh(n_j)]

        ! We now evaluate the solution at (xxmesh(m),uumesh(n))
        y(m,n)=EvCellLeg1D1D(xxmesh(m),uumesh(n),xmesh(m_i-1),xmesh(m_i),& 
                                    umesh(n_j-1),umesh(n_j),k,s,max_deg, coeff(:,:,n_j-1,m_i-1))
        end if                              
    end do 
end do         
! end lops over meshes (xxmesh) and  (uumesh)
end function AssembleLeg1D1D 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvCellLeg1D1D (x,u,xl,xr,ul,ur,k,s,max_deg, coeff) result (y)
! Evaluates the 1D+1D dimensional function from its spectral representation in Legandre's scaled polynomials.
! Takes care of one cell of the spectral representation 
! 
! The varaible (coeff) contains the coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. The variables (xl),(xr) and (ul), (ur) 
! determine the cell [xl,xr]\times[ul,ur] on the spectal grid where the spectral coefficinets are defined.
! The varaibles (x) and (u) determine the mesh point (inside where the function needs to be evaluated. 
!
! This is an part of an overloaded generic function! 
! 
! EXPECTS:
!    x -- the x coordinate of the point where the spectral representation needs to be evaluated 
!    u -- the u coordinate of the point where the spectral representation needs to be evaluated 
!    xl,xr  -- left and right boundaries of the mesh interval in x
!    ul,ur  -- left and right boundaries of the mesh interval in u
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_deg -- order of maximum allowed degree of approximation in any variable 
! RETURNS: 
!    y  -- is the value of the spectral representation   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function EvCellLeg1D1D (x,u,xl,xr,ul,ur,k,s,max_deg, coeff) result (y)
real (DP), intent (in)  :: x ! mesh value of variable (x) where the spectral representation will be evaluated
real (DP), intent (in)  :: u ! mesh value of variable (u) where the spectral representation will be evaluated
real (DP), intent (in)  :: xl,xr  ! boundaires of the cell in variable (x) where the spectral decomposition is defined (ordered)
real (DP), intent (in)  :: ul,ur  ! boundaries of the cell in variable (u) where the spectral decomposition is defined (ordered)
integer (I4B), intent (in) :: k,s                ! orders of approximations in (x) and (u) respectively
integer (I4B), intent (in) :: max_deg            ! degree of maximum allowed polynomial interpolation 
real (DP), dimension (0:,0:), intent (in)  :: coeff  ! coefficients of the spectral decomposition that will be evaluated 
real (DP) :: y ! the evaluated of the spectral representation for m=1:Nxx
real (DP) :: xx, uu ! useful local variables 
real (DP), dimension (0:k) :: phi    ! local storage for values of the basis polinomials  in x
real (DP), dimension (0:s) :: lambda ! local storage for values of the basis polinomials in u
integer (I4B) :: p,l ! local counters
!! prepare values of x for substitution in the Legendre polynomial :  
  xx = 2*(x-(xr + xl)/2)/(xr-xl);
!! Evaluate the basis functions \varhi(x) on xx 
  do p=0,k
  phi(p) = EvalHorner (LegendrePoly (p), xx)
  end do 
!! prepare values of u for substitution in the Legendre polynomial :  
  uu= 2*(u -(ur+ul)/2)/(ur-ul);
!! Evaluate the basis functions \lambda(u) on uu 
  do l=0,s
  lambda(l) = EvalHorner (LegendrePoly (l), uu)
  end do 
!! now we are ready to evaluate the entire representation 
  y=0 
  do p=0,k
  do l=0,min(max_deg-p,s)
   y=y + coeff(p,l)*phi(p)*lambda(l)
  end do 
  end do  
  end function EvCellLeg1D1D 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvCellLeg1D1D_u_vector(x,u,xl,xr,ul,ur,k,s,max_deg, coeff) result (y)
! 
! This is a light modification of the above function to handle array argument (u)
! 
! Evaluates the 1D+1D dimensional function from its spectral representation in Legandre's scaled polynomials.
! Takes care of one cell of the spectral representation 
! 
! The varaible (coeff) contains the coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. The variables (xl),(xr) and (ul), (ur) 
! determine the cell [xl,xr]\times[ul,ur] on the spectal grid where the spectral coefficinets are defined.
! The varaibles (x) and (u) determine the mesh point (inside where the function needs to be evaluated. 
!
! This is an part of an overloaded generic function! 
! 
! EXPECTS:
!    x -- the x coordinate of the point where the spectral representation needs to be evaluated 
!    u(:) -- array of points in variable u where the spectral representation needs to be evaluated 
!    xl,xr  -- left and right boundaries of the mesh interval in x
!    ul,ur  -- left and right boundaries of the mesh interval in u
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_deg -- order of maximum allowed degree of approximation in any variable 
! RETURNS: 
!    y  -- is the value of the spectral representation   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function EvCellLeg1D1D_u_vector (x,u,xl,xr,ul,ur,k,s,max_deg, coeff) result (y)
real (DP), intent (in)  :: x ! mesh value of variable (x) where the spectral representation will be evaluated
real (DP), dimension (:), intent (in)  :: u ! array of values in variable (u) where the spectral representation will be evaluated
real (DP), intent (in)  :: xl,xr  ! boundaires of the cell in variable (x) where the spectral decomposition is defined (ordered)
real (DP), intent (in)  :: ul,ur  ! boundaries of the cell in variable (u) where the spectral decomposition is defined (ordered)
integer (I4B), intent (in) :: k,s                ! orders of approximations in (x) and (u) respectively
integer (I4B), intent (in) :: max_deg            ! degree of maximum applowed polynomial interpolation 
real (DP), dimension (0:,0:), intent (in)  :: coeff  ! coefficients of the spectral decomposition that will be evaluated 
                                                     ! coeff(p,m) -- p corresponds for basis functions in x variable and
                                                     ! m corresponds to basis functions in u
real (DP), dimension (size(u)) :: y ! the evaluated spectral representation
real (DP) :: xx ! local variables to keep the points were the representation is evaluated  
real (DP), dimension (size(u)) :: uu ! local variables to keep the points were the representation is evaluated  
real (DP), dimension (0:k) :: phi    ! local storage for values of the basis polinomials in x
real (DP), dimension (0:s,size(u)) :: lambda ! local storage for values of the basis polinomials in u
integer (I4B) :: p,l ! local counters
!! prepare values of x for substitution in the Legendre polynomial :  
  xx = 2*(x-(xr + xl)/2)/(xr-xl)
!! Evaluate the basis functions \varhi(x) on xx 
  do p=0,k
  phi(p) = EvalHorner (LegendrePoly (p), xx)
  end do 
!! prepare values of u for substitution in the Legendre polynomial :  
  uu = 2*(u -(ur + ul)/2)/(ur-ul)
!! Evaluate the basis functions \lambda(u) on uu 
  do l=0,s
  lambda(l,:) = EvalHorner (LegendrePoly (l), uu)
  end do 
!! now we are ready to evaluate the entire spectral representation 
  y=0 
  do p=0,k
  do l=0,min(max_deg-p,s)
   y = y + coeff(p,l)*phi(p)*lambda(l,:)
  end do 
  end do  
  end function EvCellLeg1D1D_u_vector 


function EvCellLeg1D(x,xl,xr,coeff) result (f)
!assembles the value of function f at x from its spectral coefficients on [xl,xr]
! x is scalar
  use gaussian_mod
  use basis_fun_mod
  use poly_tool_mod
  real(DP),intent(in)::x,xl,xr
  real(DP),intent(in)::coeff(0:)
  real(DP)::f
  integer :: i
  real(DP)::xx

  xx=(2.*x-xl-xr)/(xr-xl) 
  f=0
  do i=0,size(coeff)-1
     f=f+coeff(i)*EvalHorner (LegendrePoly(i),xx)
  end do   
  
end function EvCellLeg1D

function EvCellLeg1D_x_vector(x,xl,xr,coeff) result (f)
!assembles the value of function f at x from its spectral coefficients on [xl,xr]
! x is vector
  use basis_fun_mod
  use poly_tool_mod
  real(DP),intent(in)::x(:),xl,xr
  real(DP),intent(in)::coeff(0:)
  real(DP)::f(size(x))
  integer :: i
  real(DP)::xx(size(x))

  xx=(2.*x-xl-xr)/(xr-xl) 
  f=0
  do i=0,size(coeff)-1
! for the future -- it make sense to store the values of legendre polynomials at gaussian nodes (xx are presumably
!gaussian quadrature nodes which are not going to change) rather than calculate them every time as this subroutine is
!going to be used multiple times
     f=f+coeff(i)*EvalHorner (LegendrePoly(i),xx)
  end do   
  
end function EvCellLeg1D_x_vector

function EvaluateCellDf(fmi,alm,ul,ur) result (df)
  use gaussian_mod
  real(DP),intent(in)::fmi(:,0:),& ! coefs of spectral decomposition assembled on x, first index corresponds to
                                   ! gaussian quadrature nodes in x, second corresponds to spectral decomposition coefs in u
       &alm(:,:)                       !matrix converting regular function into eigenfunctions
  real(DP)::df(size(fmi,2),size(fmi,1))
  real(DP)::ul,ur !left and right boundaries of velocity interval
  real(DP):: fhat(0:size(fmi,2)-1,size(fmi,1))  !note change of indices order here
  real(DP)::gu(size(fmi,2)),wu(size(fmi,2)) !nodes and weights of gaussian quadrature in u. In future to be stored in module
  integer::i
  fhat=matmul(alm,transpose(fmi))
  call GauLeg(Real(-1,DP),Real(1,DP),gu,wu) ! in future, will be stored rather than calculated
  gu=.5*(gu*(ur-ul)+(ur-ul))
  do i=1,size(fmi,1)
     df(:,i)=EvCellLeg1D_x_vector(Real(gu,DP),ul,ur,fhat(:,i))
  end do
end function EvaluateCellDf

!!!!!!!!!!!!!!!!
! EvalLegPol ()
!
! This one evaluates a useful matrix of values of the polynomials on the 
! point(s) x 
!
! IS A GENERIC FUNCTION! 
! 
! EXPECTS: 
!   k   -- the max degree of the polynomials
!   x  --- the points where the Leg Pol will be evaluated
!
! RETURNS: 
!   y   --- values of the first k leg. polynomials on points x
!!!!!!!!!!!!!!!!!!

function EvalLegPol(k,x) result (y)
use poly_tool_mod
use basis_fun_mod
!!!
integer (I4B), intent (in) :: k   ! max order of legenre's polynomial to be used
real (DP), intent (in) :: x  ! the point where the polynomials need to be evaluated
real (DP), dimension (0:k) :: y  ! values of the polynomials at point x
!!!
integer (I4B) :: i ! a local counter
!!!!!!!!!!!!
do i=0,k 
   y (i) = EvalHorner (LegendrePoly (i), x)
end do 
end function EvalLegPol

! this is a copy of the EvalLegPol for x -- vector

function EvalLegPol_x_vector (k,x) result (y)
use poly_tool_mod
use basis_fun_mod
!!!
integer (I4B), intent (in) :: k   ! max order of legenre's polynomial to be used
real (DP), dimension (:), intent (in) :: x  ! the point where the polynomials need to be evaluated
real (DP), dimension (size(x),0:k) :: y  ! values of the polynomials at point x
!!!
integer (I4B) :: i ! a local counter
!!!!!!!!!!!!
do i=0,k 
   y (:,i) = EvalHorner (LegendrePoly (i), x)
end do 
end function EvalLegPol_x_vector

!!!!!!!!!!!!!!!!
! EvalDerLegPol ()
!
! This one evaluates a useful matrix of values of the derivative of the Legendre's polynomials on the 
! point(s) x 
!
! IS A GENERIC FUNCTION! 
! 
! EXPECTS: 
!   k   -- the max degree of the polynomials
!   x  --- the points where the Der Leg Pol will be evaluated
!
! RETURNS: 
!   y   --- values of the first k leg. polynomials on points x
!!!!!!!!!!!!!!!!!!

function EvalDerLegPol(k,x) result (y)
use poly_tool_mod
use basis_fun_mod
!!!
integer (I4B), intent (in) :: k   ! max order of legenre's polynomial to be used
real (DP), intent (in) :: x  ! the point where the polynomials need to be evaluated
real (DP), dimension (0:k) :: y  ! values of the polynomials at point x
!!!
integer (I4B) :: i ! a local counter
!!!!!!!!!!!!
do i=0,k 
   y (i) = EvalHorner ( DerivPoly (LegendrePoly (i)), x)
end do 
end function EvalDerLegPol

! this is a copy of the EvalLegPol for x -- vector

function EvalDerLegPol_x_vector (k,x) result (y)
use poly_tool_mod
use basis_fun_mod
!!!
integer (I4B), intent (in) :: k   ! max order of legenre's polynomial to be used
real (DP), dimension (:), intent (in) :: x  ! the points where the polynomials need to be evaluated
real (DP), dimension (size(x),0:k) :: y  ! values of the polynomials at point x
!!!
integer (I4B) :: i ! a local counter
!!!!!!!!!!!!
do i=0,k 
   y (:,i) = EvalHorner ( DerivPoly (LegendrePoly (i)), x)
end do 
end function EvalDerLegPol_x_vector

!!!!!!!!!!!!!!!!
! IntLeg_dP_q_P_p ()
! 
! This one evaluates the matrix of scalar producs S_{qp}=\int_{-1}^{1} P'_{q}(y)P_{p}(y) dy
! 
! 
! EXPECTS:
!
! Sqp  --- matrix where the values of the scalar product will be stored, Sqp(q,p)=S_{qp}
!
!!!!!!!!!!!!!!!!

subroutine IntLeg_dP_q_P_p (Sqp)
use gaussian_mod
!
real (DP), dimension (0:,0:), intent (out)  :: Sqp ! the matrix of scalar products Sqp(q,p)=\int_{-1}^{1} P_{p}(y)P'_{q}(y) dy to be evaluated 
!
integer (I4B) :: k ! max degree of the polynomials 
integer (I4B) :: i,j ! local counters
real (DP), dimension (:), allocatable :: x_gnodes,x_gweights ! Gaussian nodes and weights
real (DP), dimension (:,:), allocatable :: Pp, dPq ! values of Legendre's polyn and deriv of Leg Polyn
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!
k=size(Sqp,1)-1

if (k<0) then 
   print *, "IntLeg_dP_q_P_p: zero size of the scalar product matrix."
else 

!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gnodes(1:max(k+1,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLeg_dP_q_P_p: Allocation error for variable (x_gnodes)"
     end if 
     !
  allocate (x_gweights(1:max(k+1,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLeg_dP_q_P_p: Allocation error for variable (x_gweights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  call GauLeg(Real(-1,DP),Real(1,DP),x_gnodes,x_gweights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
!!! prepare arrays for values of Leg Polynomials 
  allocate (Pp(1:max(k+1,2),0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLeg_dP_q_P_p: Allocation error for variable (Pp)"
     end if 
     !
!!! Evaluate Leg Pol and derivatives on the nodes 
  Pp = EvalLegPol_x_vector(k,Real(x_gnodes,DP))
!!! prepare arrays for values of Derivatives of Legendre's Polynomials 
  allocate (dPq(1:max(k+1,2),0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLeg_dP_q_P_p: Allocation error for variable (dPq)"
     end if 
     !
!!! Evaluate Derivative of Legendre's Polynom on the nodes 
  dPq = EvalDerLegPol_x_vector(k,Real(x_gnodes,DP))
!!! now we evaluate the products 
do i=0,k
do j=0,k
Sqp(i,j) = sum(Pp(:,j)*dPq(:,i)*x_gweights)
end do 
end do 
!!!!!!!!!!
deallocate (x_gnodes,x_gweights,Pp,dPq)
end if 

end subroutine IntLeg_dP_q_P_p



!!!!!!!!!!!!!!!!
! IntLegP_q_dP_p ()
! 
! This one evaluates the matrix of scalar producs S_{qp}=\int_{-1}^{1} P_{q}(y)P'_{p}(y) dy
! 
! 
! EXPECTS:
!
! Sqp  --- matrix where the values of the scalar product will be stored Sqp(q,p)=S_{qp}
!
!!!!!!!!!!!!!!!!

subroutine IntLegP_q_dP_p (Sqp)
use gaussian_mod
!
real (DP), dimension (0:,0:), intent (out)  :: Sqp ! the matrix of scalar products Sqp(q,p)=\int_{-1}^{1} P_{q}(y)P'_{p}(y) dy to be evaluated 
!
integer (I4B) :: k ! max degree of the polynomials 
integer (I4B) :: i,j ! local counters
real (DP), dimension (:), allocatable :: x_gnodes,x_gweights ! Gaussian nodes and weights
real (DP), dimension (:,:), allocatable :: Pq, dPp ! values of Legendre's polyn and deriv of Leg Polyn
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!
k=size(Sqp,1)-1

if (k<0) then 
   print *, "IntLegP_q_dP_p: zero size of the scalar product matrix."
else 

!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gnodes(1:max(k+1,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLegP_q_dP_p: Allocation error for variable (x_gnodes)"
     end if 
     !
  allocate (x_gweights(1:max(k+1,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLegP_q_dP_p: Allocation error for variable (x_gweights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  call GauLeg(Real(-1,DP),Real(1,DP),x_gnodes,x_gweights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
!!! prepare arrays for values of Leg Polynomials 
  allocate (Pq(1:max(k+1,2),0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLegP_q_dP_p: Allocation error for variable (Pq)"
     end if 
     !
!!! Evaluate Leg Pol and derivatives on the nodes 
  Pq = EvalLegPol_x_vector(k,Real(x_gnodes,DP))
!!! prepare arrays for values of Derivatives of Legendre's Polynomials 
  allocate (dPp(1:max(k+1,2),0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "IntLegP_q_dP_p: Allocation error for variable (dPp)"
     end if 
     !
!!! Evaluate Derivative of Legendre's Polynom on the nodes 
  dPp = EvalDerLegPol_x_vector(k,Real(x_gnodes,DP))
!!! now we evaluate the products 
do i=0,k
do j=0,k
Sqp(i,j) = sum(Pq(:,i)*dPp(:,j)*x_gweights)
end do 
end do 
!!!!!!!!!!
deallocate (x_gnodes,x_gweights,Pq,dPp)
end if 

end subroutine IntLegP_q_dP_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! added 09/05/08 Alex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DecompLeg1Dvar1Ddecomp()
!
! This function will return coefficients of the spectral decomposition of a two dimensional function 
! f(x,u,t) with respect to the second variable. The first variable will be substituted with a mesh point. 
! 
! The basis functions are scaled legende'r polynomials, the degree of the spectral decomposition is given by 
! variable (s). The mesh for spectral deocmposition is (umesh)
!     
! 
! EXPECTS:
!    x -- points in variable x where spectral decomposition is evaluated
!    umesh -- meshpoints in variable u
!    s --- order of polynomial approximation in variable u
!    f --- is a dummy external function whose spectral coefficients need to be evaluated. 
!      (ATTN: the function DecompLeg1D1D is coded so that (f) accepts first scalar variable and second vector variable)  
!    u_gauss_nodes --- array of guass nodes to be used on each cell for the spectral decomposition
!    u_wleg_pol --- array of values of legendre's polynomials on the (u_gauss_nodes) multiplied by the corresponding (u_gauss_weights)
! 
! RETURNS: 
!    coeff --- array of the coefficents 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function DecompLeg1Dvar1Ddecomp(x,umesh,s,t,f,u_gauss_nodes,u_wleg_pol) result (coeff)
!
real (DP), intent (in) :: x ! the point in varaible x where the spectral representation needs to be evaluated
real (DP), dimension (0:), intent (in) :: umesh ! the mesh in variable u
integer (I4B), intent (in) :: s ! the degree of the polynomial basis in variable u
real (DP), intent (in) :: t ! time (will be passed to the dumy function)
real (DP), dimension (:), intent (in)  :: u_gauss_nodes  ! nodes of the gauss integration formula to be used in evaluation
real (DP), dimension (0:,:), intent (in) :: u_wleg_pol ! leg polynomials evalueated at guassian nodes (u_gauss_nodes)
                                                       ! u_wleg_pol(s,i)
                                                       ! s -- is the degree(number) of the basis function
                                                       ! i is the number fo the gaussian node (u_gauss_node)
real (DP), dimension (0:s,size(umesh)-1) :: coeff ! coefficients of the spectral decomposition   
real (DP), dimension (size(u_gauss_nodes)) :: u ! variable to store points for function evaluation
integer :: l,i ! local counters
 ! interface block for the dummy function f: (x-scalar, u-vector)
 ! first variable must be scalar and second must be vector! 
 interface 
   function f (x, u, t) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: t ! the time variable 
    real (DP), dimension (size(u)) :: y ! values of the function  
   end function f
 end interface                                                       
! now we go through each cell
  do i=1,size(umesh)-1
  !!! prepare integration with the cell
    u = (umesh(i)-umesh(i-1))*Real(u_gauss_nodes,DP)/2.0_DP+(umesh(i) + umesh(i-1))/2.0_DP
  !!! Integration within the cell [umesh(i-1),umesh(i)]
    do l=0,s
      coeff(l,i)= sum(f(x,u,t)*u_wleg_pol(l,:))*(2*l+1)/2.0_DP
    end do 
  !!! end integration within the cell
  end do 
!!! End evaluating spectral decomposition

end function DecompLeg1Dvar1Ddecomp

!!!!!!!! end added 09/05/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! added Alex 12/05/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvCellLegStdrd1D1D (x,u,k,s,max_deg,coeff) result (y)
!
!
! Evaluates the 1D+1D dimensional function from its spectral representation in Legandre's scaled polynomials.
!
! ATTENTION:
! This one evaluates the spectral representation given by (coeff) on [-1,1]\times [-1,1]. 
! If (x) and (y) are not in this interval -- RETURNS GARBAGE!! 
!
! This is an overloaded function  
! 
! The variable (coeff) contains the coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. The variables (xl),(xr) and (ul), (ur) 
! determine the cell [xl,xr]\times[ul,ur] on the spectal grid where the spectral coefficinets are defined.
! The varaibles (x) and (u) determine the mesh point (inside where the function needs to be evaluated. 
!
! This is an part of an overloaded generic function! 
! 
! EXPECTS:
!    x -- the x coordinate of the point where the spectral representation needs to be evaluated 
!    u -- the u coordinate of the point where the spectral representation needs to be evaluated 
!      -- left and right boundaries of the mesh interval in x are assumed to be [-1,1]
!      -- left and right boundaries of the mesh interval in u are assumed to be [-1,1]
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_deg -- order of maximum allowed degree of approximation in any variable 
! RETURNS: 
!    y  -- is the value of the spectral representation   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function EvCellLegStdrd1D1D (x,u,k,s,max_deg, coeff) result (y)
real (DP), intent (in)  :: x ! mesh value of variable (x) where the spectral representation will be evaluated
real (DP), intent (in)  :: u ! mesh value of variable (u) where the spectral representation will be evaluated
integer (I4B), intent (in) :: k,s                ! orders of approximations in (x) and (u) respectively
integer (I4B), intent (in) :: max_deg            ! degree of maximum allowed polynomial interpolation 
real (DP), dimension (0:,0:), intent (in)  :: coeff  ! coefficients of the spectral decomposition that will be evaluated 
real (DP) :: y ! the evaluated of the spectral representation for m=1:Nxx
real (DP), dimension (0:k) :: phi    ! local storage for values of the basis polinomials  in x
real (DP), dimension (0:s) :: lambda ! local storage for values of the basis polinomials in u
integer (I4B) :: p,l ! local counters
!! Evaluate the basis functions \varhi(x) on xx 
  do p=0,k
  phi(p) = EvalHorner (LegendrePoly (p), x)
  end do 
!! Evaluate the basis functions \lambda(u) on uu 
  do l=0,s
  lambda(l) = EvalHorner (LegendrePoly (l), u)
  end do 
!! now we are ready to evaluate the entire representation 
  y=0 
  do p=0,k
  do l=0,min(max_deg-p,s)
   y=y + coeff(p,l)*phi(p)*lambda(l)
  end do 
  end do  
  end function EvCellLegStdrd1D1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvCellLegStdrd1D1D_u_vector(x,u,k,s,max_deg, coeff) result (y)
! 
! This is a light modification of the above function to handle array argument (u)
! 
! Evaluates the 1D+1D dimensional function from its spectral representation in Legandre's scaled polynomials.
! Takes care of one cell of the spectral representation. 
! 
! ATTENTION: this one assumes that the cell is [-1,1] (to eliminate the unnecessry operations).  
! ASSUMES THAT (x,u)\in[-1,1]\times[-1,1]  does not make sense otherwise    !!!  
! 
! The varaible (coeff) contains the coefficents of the spectral fepresentation by scaled Legendre polynomials
! of a 1D+1D dimensional functions treating the 1D first variable and second 1D variable independently. 
! Namely, the grid of the spectral decomposition is rectangular and the basis functions are products of 1D basis 
! functions. The degrees of the highest polynomial approximation are given by the parameters (k) and (s) for 
! the first and second variable independently. The parameter (max_degree) determines which products of basis 
! functions are considered in the spectral deomposition. If (max_degree)>(k)+(s) then all produts of basis 
! functions for first and second variables will be considered. The variables (xl),(xr) and (ul), (ur) 
! determine the cell [xl,xr]\times[ul,ur] on the spectal grid where the spectral coefficinets are defined.
! The varaibles (x) and (u) determine the mesh point (inside where the function needs to be evaluated. 
!
! This is an part of an overloaded generic function! 
! 
! EXPECTS:
!    x -- the x coordinate of the point where the spectral representation needs to be evaluated 
!    u(:) -- array of points in variable u where the spectral representation needs to be evaluated 
!      -- left and right boundaries of the mesh interval in x are assumed to be [-1,1]
!      -- left and right boundaries of the mesh interval in u are assumed to be [-1,1]
!    k  --- order of polynomial approximation in variable x
!    s  --- order of polynomial approximation in variable u
!    max_deg -- order of maximum allowed degree of approximation in any variable 
! RETURNS: 
!    y  -- is the value of the spectral representation   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function EvCellLegStdrd1D1D_u_vector (x,u,k,s,max_deg, coeff) result (y)
real (DP), intent (in)  :: x ! mesh value of variable (x) where the spectral representation will be evaluated
real (DP), dimension (:), intent (in)  :: u ! array of values in variable (u) where the spectral representation will be evaluated
integer (I4B), intent (in) :: k,s                ! orders of approximations in (x) and (u) respectively
integer (I4B), intent (in) :: max_deg            ! degree of maximum applowed polynomial interpolation 
real (DP), dimension (0:,0:), intent (in)  :: coeff  ! coefficients of the spectral decomposition that will be evaluated 
                                                     ! coeff(p,m) -- p corresponds for basis functions in x variable and
                                                     ! m corresponds to basis functions in u
real (DP), dimension (size(u)) :: y ! the evaluated spectral representation
real (DP), dimension (0:k) :: phi    ! local storage for values of the basis polinomials in x
real (DP), dimension (0:s,size(u)) :: lambda ! local storage for values of the basis polinomials in u
integer (I4B) :: p,l ! local counters
!! Evaluate the basis functions \varhi(x) on x 
  do p=0,k
  phi(p) = EvalHorner (LegendrePoly (p), x)
  end do 
!! Evaluate the basis functions \lambda(u) on u 
  do l=0,s
  lambda(l,:) = EvalHorner (LegendrePoly (l), u)
  end do 
!! now we are ready to evaluate the entire spectral representation 
  y=0 
  do p=0,k
  do l=0,min(max_deg-p,s)
   y = y + coeff(p,l)*phi(p)*lambda(l,:)
  end do 
  end do  
  end function EvCellLegStdrd1D1D_u_vector 

!!!!!!!!! end added Alex 12/05/08 !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module spectral_tools_mod

