! 10/12/08 Alex
!
! This Module deals with the functions and procedures for 
! evaluating different integrals of the solution  
! possibly also for calculating the error in the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
module sol_integr_mod
!
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
use basis_fun_mod
use poly_tool_mod
!
implicit none
!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IntMomentsSolSpec1D1D
!
! This subroutine takes the spectral representation of the solution and 
! calculates integral of the first three moments of the solution 
! using the given grids in x and u and the and the given order of the gauss formula and the mesh 
! refinement parameter
!
! EXPECTS 
! int_xmesh -- mesh where the solution needs to be integrated in x
! int_umesh -- mesh where the solution needs to be integrated in u
! moments_x_gauss_nodes, moments_x_gauss_weights  -- the gauss nodes and guass weights in x to be used in x for moments integration
! moments_u_gauss_nodes, moments_u_gauss_weights -- the gauss nodes and weights in u to be used in momenets integration,
! moments_refine_x, moments_refine_u, --- coefficients in mesh refinements will be used to refine cells of (int_xmesh) 
!                          every cell will be divided into refine_int_xmesh\times moment_int_umesh parts
! (moments_refine_x, moments_refine_u) is set in the common_variables_mod !    
! xmesh --- mesh points in x used in the defintition of the spetral solution
! umesh --- mesh points in u used in the defintition of the spetral solution
! k --- order of the spectral represent in x
! s --- order of the spectral represenataion in u.
! max_deg --- max degree of the basis polynomial (in both x and u)
! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
!
! RETURNS
! 
! ans = the first three moments of f 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function IntMomentsSolSpec1D1D(int_xmesh,int_umesh,moments_x_gauss_nodes,moments_x_gauss_weights,&
                         moments_u_gauss_nodes,moments_u_gauss_weights,&
                         moments_refine_x,moments_refine_u,xmesh,umesh,k,s,max_deg,coeff) result (ans)
use spectral_tools_mod  ! this one has the assembling routine
!!
real (DP), dimension (0:), intent (in) :: int_xmesh,int_umesh ! mesh points/cells in x and u used 
                                                              ! for evaluation of integrals (can be 
                                                              ! different from meshes where spetral solution is defined) 
real (DP), dimension (:), intent (in) :: moments_x_gauss_nodes,moments_x_gauss_weights,&
                                         moments_u_gauss_nodes,moments_u_gauss_weights 
                                                              ! gauss nodes and weights in x and u
                                                              ! to be used in evaluation of moments 
integer (I4B), intent (in) :: moments_refine_x,moments_refine_u ! these are coefficients of cell refinements in x and u
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in variable x where the spectal solution is defined
real (DP), dimension (0:), intent (in) :: umesh ! mesh points in variable u where the spectal solution is defined
integer (I4B), intent (in) :: k,s,max_deg ! order in x, order in u and max degree of polynomial approximation 
real (DP), dimension (0:,0:,:,:), intent (in) :: coeff ! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
real (DP), dimension (3) :: ans ! the value of the computed L_2 norm of the first three moments of the solution 
!
integer :: i,j,ci,cj ! local counters
real (DP), dimension (:), allocatable :: xxmesh, uumesh, xxweights, uuweights ! meshes and gauss weights where the spectral representation will be evaluated 
real (DP), dimension (:,:), allocatable :: weighted_sol ! local variable where the solution is stored during the integration 
                                       ! weighted_sol(j,i) j --- node number in x, i -- node number in u 
integer :: loc_alloc_stat ! to keep allocation status 
integer :: x_gaussorder,u_gaussorder ! local variable to keep the order of gauss integration for moments
real (DP) :: du,dx ! to keep temp calculations
real (DP), dimension (3) :: ans_tempi ! to keep temporary calculations
integer (I4B) :: rix,riu ! these are local coefficients of cell refinements in x and u
! 
! just in case, check if garbage came in moments_refine_x, moments_refine_u: ---
if (moments_refine_x < 1) then  
   rix=1
   print *, "IntMomentsSolSpec1D1D: supplied moments_refine_x parameter is zero or negative. Set moments_refine_x=1"
   else 
   rix=moments_refine_x
   end if 
if (moments_refine_u < 1) then  
   riu=1
   print *, "IntMomentsSolSpec1D1D: supplied moments_refine_u parameter is zero or negative. Set moments_refine_u=1"
   else 
   riu=moments_refine_u
   end if 
x_gaussorder=size(moments_x_gauss_nodes)
u_gaussorder=size(moments_u_gauss_nodes)
!! allocate space for inter cell nodes and weights 
 allocate (xxmesh(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsSolSpec1D1D: Allocation error for variable (xxmesh)"
   end if 
allocate (uumesh(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsSolSpec1D1D: Allocation error for variable (uumesh)"
   end if 
allocate (xxweights(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsSolSpec1D1D: Allocation error for variable (xxweights)"
   end if 
allocate (uuweights(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsSolSpec1D1D: Allocation error for variable (uuweights)"
   end if 
!! allocate space for solution evaluated on (xxmesh)\times(uumesh)
allocate (weighted_sol(1:x_gaussorder*rix,1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsSolSpec1D1D: Allocation error for variable (weighted_sol)"
   end if 
!!
! we prepare a larger array of weights (dublicates of gauss weights) to be used on subdivided cells
do cj=1,rix
 xxweights(1+(cj-1)*x_gaussorder:cj*x_gaussorder) = Real (moments_x_gauss_weights, DP)
end do 
do ci=1,riu
 uuweights(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = Real (moments_u_gauss_weights, DP)
end do 
! ready to calculate the moments:
ans=0 
do i=1,size(int_umesh)-1      !loop in integration cells in u
 ! now each cell [int_umesh(i-1),int_umesh(i)] is subdivided into (riu) parts:
 ! and the new inside-cell meshes are recorded. 
 du = ( int_umesh(i) - int_umesh(i-1) ) / Real(riu,DP)
 do ci=1,riu
  uumesh(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = int_umesh(i-1) + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
 end do
 ans_tempi=0 
 do j=1,size(int_xmesh)-1      !loop in integration cells in x
   ! now each cell [int_xmesh(j-1),int_xmesh(j)] is subdivided into (rix) parts:
   ! and the new inside-cell meshes are recorded. 
   dx = ( int_xmesh(j) - int_xmesh(j-1) ) / Real(rix,DP)
   do cj=1,rix
   xxmesh(1+(cj-1)*x_gaussorder:cj*x_gaussorder) = int_xmesh(j-1) + dx*(Real(cj,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
   end do  
   ! 
   ! first we evaluate the solution at nodes (xxmesh,uumesh)
   weighted_sol = AssembleLeg1D1D(xxmesh,uumesh,xmesh,umesh,k,s,max_deg,coeff)
   ! we then multiply the solution by gauss weights in both x and u:
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uuweights*xxweights(cj)
   end do 
   ! Now it is time to calculate the moments 
   ! the first moment:
   ans_tempi(1)=ans_tempi(1) + sum(weighted_sol)*dx
   ! the second moment
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uumesh
   end do 
   ans_tempi(2)=ans_tempi(2) + sum(weighted_sol)*dx
   ! the third moment 
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uumesh
   end do 
   ans_tempi(3)=ans_tempi(3)+sum(weighted_sol)*dx
   ! end calculating moments in cell [int_xmesh(j-1),int_mesh(j)]\times[int_umesh(i-1), int_umesh(i)]
  end do  ! end loop in (j) --- cells in x
 ans = ans + ans_tempi*du
end do  ! end loop in cells in u
ans=ans/Real(4,DP)                   ! division by 4 comes from the integration
deallocate (uumesh,xxmesh,xxweights,uuweights,weighted_sol)
!
end function IntMomentsSolSpec1D1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IntMomementsSolHOV1D
!
! This function computes the moments of the discrete solution. First it calculates the 
! spectral coefficients of the representation from the spectral-characteristic coefficients. Then it calls the
! function that calculates the moments
!
! 
! Do not detouch from the main program
!
! This subrouting is highly dependent on the variables of the main probram and can not be called before 
! the variables k,s,max_deg,xmesh,umesh,Aml,moment_x_gauss_nodes,moment_u_gauss_nodes,
! moment_x_gauss_weights,moments_u_gauss_weights,M,N ... are initialized
!
! Most of the variables are looked up in the common_variables_mod 
!
! USES
! int_xmesh -- mesh where the solution needs to be integrated in x
! int_umesh -- mesh where the solution needs to be integrated in u
! moments_x_gauss_nodes, moments_x_gauss_weights  -- the gauss nodes and guass weights in x to be used in x for moments integration
! moments_u_gauss_nodes, moments_u_gauss_weights -- the gauss nodes and weights in u to be used in momenets integration,
! moments_refine_x, moments_refine_u, --- coefficients in mesh refinements will be used to refine cells of (int_xmesh) 
!                          every cell will be divided into refine_int_xmesh\times moment_int_umesh parts
! (moments_refine_x, moments_refine_u) is set in the common_variables_mod !    
! xmesh --- mesh points in x
! umesh --- mesh points in u 
! k --- order of the spectral represent in x
! s --- order of the spectral represenataion in u.
! max_deg --- max degree of the basis polynomial (in both x and u)
! curr_time --- coordinate time (used in exact solution)
! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
!
! RETURNS
! 
! ans --- is the value of the first three moments of the numerical solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function IntMomentsSolHOV1D(int_xmesh,int_umesh,ftil1,ftil2) result (ans)
use common_variables_mod, only: k,s,max_deg,xmesh,umesh,Aml,moments_refine_x, moments_refine_u,&
    moments_x_gauss_nodes,moments_x_gauss_weights,moments_u_gauss_nodes,moments_u_gauss_weights,&
    M,N
!!!
real (DP), dimension (0:), intent (in) :: int_xmesh, int_umesh ! meshes to use for the evaluation of moments 
real (DP), dimension (0:,0:,:,:), intent (in) :: ftil1,ftil2 ! coeff(p,m,i,j) -- characteristic coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
real (DP), dimension (3,2) :: ans ! the computed values of the moments
real (DP), dimension (0:k,0:s,M,N) :: f  ! auxiliary variable -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  f(p,m,i,j)( ftil(p,m,i,j), fti2(p,m,i,j):)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- i is the cell in u
                                                     ! -- j is the cell in x
!!!
integer  :: i,j,l,q,m_count    ! local counters
!! 
! convert the characteristic coefficients into non-characteristic, that is
! evaluates the coefficients of the spectral representation from the characteristic representation
   f=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f(q,l,i,j)= f(q,l,i,j) + Aml(l,m_count,i)*ftil1(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
   ans(:,1)=IntMomentsSolSpec1D1D(int_xmesh,int_umesh,moments_x_gauss_nodes,moments_x_gauss_weights,&
                         moments_u_gauss_nodes,moments_u_gauss_weights, &
                         moments_refine_x, moments_refine_u,xmesh,umesh,k,s,max_deg,f)
 ! convert the characteristic coefficients into non-characteristic, that is
 ! evaluate the coefficients of the spectral representation from the characteristic representation
   f=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f(q,l,i,j)= f(q,l,i,j) + Aml(l,m_count,i)*ftil2(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
   ans(:,2)=IntMomentsSolSpec1D1D(int_xmesh,int_umesh,moments_x_gauss_nodes,moments_x_gauss_weights,&
                         moments_u_gauss_nodes,moments_u_gauss_weights, &
                         moments_refine_x, moments_refine_u,xmesh,umesh,k,s,max_deg,f)
end function IntMomentsSolHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IntMomentsExactSol1D1D
!
! This subroutine takes the spectral representation of the solution and 
! calculates integral of the first three moments of the solution 
! using the given grids in x and u and the and the given order of the gauss formula and the mesh 
! refinement parameter
!
! EXPECTS 
! int_xmesh -- mesh where the solution needs to be integrated in x
! int_umesh -- mesh where the solution needs to be integrated in u
! moments_x_gauss_nodes, moments_x_gauss_weights  -- the gauss nodes and guass weights in x to be used in x for moments integration
! moments_u_gauss_nodes, moments_u_gauss_weights -- the gauss nodes and weights in u to be used in momenets integration,
! moments_refine_x, moments_refine_u, --- coefficients in mesh refinements will be used to refine cells of (int_xmesh) 
!                          every cell will be divided into refine_int_xmesh\times moment_int_umesh parts
! ex_sol_fun --- name of the external function that gives the exact solution
!             ! it expects ex_sol_fun(x,u,t), 
!             !  x is vector, u is a vector, and t is a scalar!
! RETURNS
! 
! ans = the first three moments of the exact solution 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function IntMomentsExactSol1D1D(int_xmesh,int_umesh, moments_x_gauss_nodes, moments_x_gauss_weights,& 
                                moments_u_gauss_nodes, moments_u_gauss_weights, &
                                moments_refine_x, moments_refine_u, t,ex_sol_fun) result (ans)
use spectral_tools_mod  ! this one has the assembling routine
!!
real (DP), dimension (0:), intent (in) :: int_xmesh,int_umesh ! mesh points/cells in x and u used 
                                                              ! for evaluation of integrals (can be 
                                                              ! different from meshes where spetral solution is defined) 
real (DP), dimension (:), intent (in) :: moments_x_gauss_nodes, moments_x_gauss_weights,&
                                         moments_u_gauss_nodes, moments_u_gauss_weights 
                                                              ! gauss nodes and weights in x and u
                                                              ! to be used in evaluation of moments 
integer (I4B), intent (in) :: moments_refine_x, moments_refine_u! these are coefficients of cell refinements in x and u
real (DP), intent (in) :: t     ! the coordinate time -- will be passed to the exact solution  
real (DP), dimension (3) :: ans ! the value of the computed L_2 norm of the first three moments of the solution 
!
integer :: i,j,ci,cj ! local counters
real (DP), dimension (:), allocatable :: xxmesh, uumesh, xxweights, uuweights ! meshes and gauss weights where the spectral representation will be evaluated 
real (DP), dimension (:,:), allocatable :: weighted_sol ! local variable where the solution is stored during the integration 
                                 ! weighted_sol(j,i) : j-- number of the node in x , i -- number of the node in u
integer :: loc_alloc_stat ! to keep allocation status 
integer :: x_gaussorder,u_gaussorder ! local variable to keep the order of gauss integration for moments
real (DP) :: du,dx ! to keep temp calculations
real (DP), dimension (3) :: ans_tempi ! to keep temporary calculations
integer (I4B) :: rix,riu ! these are local coefficients of cell refinements in x and u
 ! interface block for the dummy function ex_sol_fun: (x-scalar, u-vector)
 ! first variable must be scalar and second must be vector! 
 interface 
   function ex_sol_fun (x, u, t) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), dimension (:), intent (in)  :: x     ! the vector of values of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: t ! the time variable 
    real (DP), dimension (size(x),size(u)) :: y ! values of the function  
   end function ex_sol_fun
 end interface
! 
! just in case, check if garbage came in moments_refine_x, moments_refine_u: ---
if (moments_refine_x < 1) then  
   rix=1
   print *, "IntMomentsExactSol1D1D: supplied moments_refine_x parameter is zero or negative. Set moments_refine_x=1"
   else 
   rix=moments_refine_x
   end if 
if (moments_refine_u < 1) then  
   riu=1
   print *, "IntMomentsExactSol1D1D: supplied moments_refine_u parameter is zero or negative. Set moments_refine_u=1"
   else 
   riu=moments_refine_u
   end if 
x_gaussorder=size(moments_x_gauss_nodes)
u_gaussorder=size(moments_u_gauss_nodes)
!! allocate space for inter cell meshes ... and nodes
 allocate (xxmesh(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsExactSol1D1D: Allocation error for variable (xxmesh)"
   end if 
allocate (uumesh(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsExactSol1D1D: Allocation error for variable (uumesh)"
   end if 
allocate (xxweights(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsExactSol1D1D: Allocation error for variable (xxweights)"
   end if 
allocate (uuweights(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsExactSol1D1D: Allocation error for variable (uuweights)"
   end if 
!! allocate space for solution evaluated on (xxmesh)\times(uumesh)
   allocate (weighted_sol(1:x_gaussorder*rix,1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "IntMomentsExactSol1D1D: Allocation error for variable (weighted_sol)"
   end if    
!! 
! we prepare a larger array of weights (dublicates of gauss weights) to be used on subdivided cells
do cj=1,rix
 xxweights(1+(cj-1)*x_gaussorder:cj*x_gaussorder) = Real (moments_x_gauss_weights, DP)
end do 
do ci=1,riu
 uuweights(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = Real (moments_u_gauss_weights, DP)
end do 
! ready to calculate the moments:
ans=0 
do i=1,size(int_umesh)-1      !loop in integration cells in u
 ! now each cell [int_umesh(i-1),int_umesh(i)] is subdivided into (riu) parts:
 ! and the new inside-cell meshes are recorded. 
 du = ( int_umesh(i) - int_umesh(i-1) ) / Real(riu,DP)
 do ci=1,riu
  uumesh(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = int_umesh(i-1) + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
 end do
 ans_tempi=0 
 do j=1,size(int_xmesh)-1      !loop in integration cells in x
   ! now each cell [int_xmesh(j-1),int_xmesh(j)] is subdivided into (rix) parts:
   ! and the new inside-cell meshes are recorded. 
   dx = ( int_xmesh(j) - int_xmesh(j-1) ) / Real(rix,DP)
   do cj=1,rix
   xxmesh(1+(cj-1)*x_gaussorder:cj*x_gaussorder) = int_xmesh(j-1) + dx*(Real(cj,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
   end do  
   ! 
   ! first we evaluate the EXACT solution at nodes (xxmesh,uumesh)
   weighted_sol = ex_sol_fun (xxmesh, uumesh, t)
   ! we then multiply the solution by gauss weights in both x and u:
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uuweights*xxweights(cj)
   end do 
   ! Now it is time to calculate the moments 
   ! the first moment:
   ans_tempi(1)=ans_tempi(1)+sum(weighted_sol)*dx
   ! the second moment
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uumesh
   end do 
   ans_tempi(2)=ans_tempi(2)+sum(weighted_sol)*dx
   ! the tird moment 
   do cj=1,x_gaussorder*rix
   weighted_sol(cj,:) = weighted_sol(cj,:)*uumesh
   end do 
   ans_tempi(3)=ans_tempi(3)+sum(weighted_sol)*dx
   ! end calculating moments in cell [int_xmesh(j-1),int_mesh(j)]\times[int_umesh(i-1), int_umesh(i)]
  end do  ! end loop in (j) --- cells in x
 ans = ans + ans_tempi*du
end do  ! end loop in cells in u
ans=ans/Real(4,DP)                   ! division by 4 comes from the integration
deallocate (uumesh,xxmesh,xxweights,uuweights,weighted_sol)
!
end function IntMomentsExactSol1D1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalMacroXmeshSolHOV1D
!
!
! This subroutine computes values of macroparameters by callign the macroparamters integrator.
! Essentially, it just puts together a lot of parameters, converts the characteristic variables in 
! regular ones and then calls the EvalMacroSol1D1D. The evaluation is done at nodes 
! on each interval in (xmesh). Nodes in x are based on (moments_x_gauss_nodes) and (moments_x_refine) 
! and are calculated using the parameters (moments_x_gauss_order). For integration in the velocity, 
! (moments_u_gauss_nodes) and (moments_u_refine) are used. 
! 
! Highly dependent on the outside program
!
! Most of the variables are looked up in the Common_variables_mod
!
! do not call untill moments_x_gauss_nodes,moments_x_gauss_weights,&
!                                moments_u_gauss_nodes,moments_u_gauss_weights,& 
!                                moments_refine_x,moments_refine_u, etc are initialized 
!
! EXPECTS: 
! dens(m,j)  --- array to store density
! avel(m,j)  --- array to store average velocity
! tempr(m,j) --- array to store temperature
!                (m,j) m --- the node on the interval I_{j}, j is the interval in x form xmesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine EvalMacroXmeshSolHOV1D(dens,avel,tempr,ftil1,ftil2,x_nodes)
use common_variables_mod, only: k,s,max_deg,xmesh,umesh,Aml,M,N,moments_x_gauss_nodes, &
                               moments_u_gauss_nodes,moments_u_gauss_weights, & 
                               moments_refine_x,moments_refine_u,gasR
!!
real (DP), dimension (:,:), intent (out) :: dens, avel, tempr ! arrays to store density, average velocity and 
                                        ! temperature. E.G., dens(m,j) 
                                        ! m --- the node on the interval I_{j}, j --is the interval in x form xmesh
real (DP), dimension (0:,0:,:,:), intent (in) :: ftil1,ftil2 ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,i,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- i is the cell in u
                                   ! -- j is the cell in x
!!! 
real (DP), dimension (0:k,0:s,M,N) :: f1,f2  ! auxiliary variable -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  f(p,m,i,j)( ftil(p,m,i,j), fti2(p,m,i,j):)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- i is the cell in u
                                                     ! -- j is the cell in x
integer  :: i,j,l,q,m_count    ! local counters
real (DP), dimension (:), intent (in) :: x_nodes ! nodes to be used in the evaluation of macraparameters
                                                 ! obtained by subdiving interval [-1,1] into (moments_refine_x)
                                                 ! parts and scaling using nodes (moments_x_gauss_nodes) on each part
                                                 ! these nodes will be passed to the standard evaluation procedure to 
                                                 ! bypass scaling [x_{j-1/2},x_{j+1/2}] to [-1,1] and back
!!!
! we convert the characteristic coefficients into non-characteristic, that is
! evaluates the coefficients of the spectral representation from the characteristic representation
   f1=0
   f2=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f1(q,l,i,j)= f1(q,l,i,j) + Aml(l,m_count,i)*ftil1(q,m_count,i,j)
     f2(q,l,i,j)= f2(q,l,i,j) + Aml(l,m_count,i)*ftil2(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
  !! now we call the procedure to evaluate the macroparameters. This one evaluates the macroparameters 
  !! based on the meshes (xmesh) and (umesh).
   call EvalMacroSolXmeshXnodes1D(dens,avel,tempr,x_nodes,umesh,moments_u_gauss_nodes,moments_u_gauss_weights,&
                                moments_refine_u,k,s,max_deg,f1,f2,gasR)
end subroutine EvalMacroXmeshSolHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalMacroSolXmeshXnodes1D
!
! This subroutine computes values of macroparameters. The evaluation is done on each I_{j} cell in x 
! (given by [xmesh]) ,  at all gauss nodes (given by [moments_x_gauss_nodes] and [moments_x_refine]).
! 
!  
! Integrals are calculated using the main mesh in the velocity variable (umesh). Nodes in u are based on 
! (moments_u_gauss_nodes) and (moments_u_refine) ((moments_u_gauss_nodes) are calculated using the parameters 
! (moments_u_gauss_order).  
!  
!
! EXPECTS: 
! dens(m,j)  --- array to store density
! avel(m,j)  --- array to store average velocity
! tempr(m,j) --- array to store temperature
!                (m,j) m --- the node on the interval I_{j}, j is the interval in x form xmesh
!                 Attention::
!                 m must have be the dimension (moments_x_gauss_order)x(moments_u_refine) 
!                 j must be size(xmesh)-1
! x_nodes(p) --- nodes on [-1,1] corresponding to points on each [x_{j-1/2},x_{j+1/2}] on which 
!              where the macroparamters need to be evaluated. 1<=p<=moments_x_gauss_nodes*moments_x_refine 
! umesh(i) --- mesh points in u used in the defintition of the spetral solution ,
! moments_u_gauss_nodes,moments_u_gauss_weights --- gauss nodes and weights to be used for integration in u 
! riu -- coefficient of mesh refinement in u
! k,s,max_deg --- k is the degree of spectral approximation in x, s it the degree of polynomial approximation in u, and 
!                 max_deg -- is the maximum allowed degree in the spectral decomposition
! f1,f2 --- spectral coefficients! f12(p,m,i,j) -- p is the index in the basis functions in x 
!                                  ! -- m is the index in the basis functions in u
!                                   ! -- i is the cell in u
!                                   ! -- j is the cell in x
! R --- the gas constant
!
!
! RETURNS
! 
! dens, avel, tempr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine EvalMacroSolXmeshXnodes1D(dens,avel,tempr,x_nodes,umesh,moments_u_gauss_nodes,moments_u_gauss_weights,riu,&
                                k,s,max_deg,f1,f2,R)
use spectral_tools_mod
!
real (DP), dimension (:,:), intent (out) :: dens, avel, tempr ! arrays to store density, average velocity and 
                                        ! temperature. E.G., dens(m,j) 
                                        ! m --- the node on the interval I_{j}, j --is the interval in x form xmesh
real (DP), dimension (:), intent (in) :: x_nodes 
                                        ! nodes in x on interval [-1,1] to be used in macroparameter evaluations. Essentially, 
                                        ! these correpond to the gauss nodes obtained after refining each an interval I_{j} of [xmesh]
                                        ! We assume that the same gauss nodes and the same refining coefficient are used 
                                        ! on all intervals. Since actual x does not show up, we will keep gauss nodes scaled on [-1,1] for speed
real (DP), dimension (0:), intent (in) :: umesh ! main mesh in variable u, where the spectral representation in defined
real (DP), dimension (:), intent (in) :: moments_u_gauss_nodes, moments_u_gauss_weights 
                                                              ! gauss nodes and weights in u
                                                              ! to be used in evaluation of macoparameters 
integer (I4B), intent (in) :: riu ! these are coefficients of cell refinements in x and u
integer (I4B), intent (in) :: k,s,max_deg ! order in x, order in u and max degree of polynomial approximation 
real (DP), dimension (0:,0:,:,:), intent (in) :: f1,f2 ! f12(p,m,i,j) -- coefficients of spectral decomposition (converted from charactristic variables) 
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m is the index in the basis functions in u
                                   ! -- i is the cell in u
                                   ! -- j is the cell in x
real (DP), intent (in) :: R ! the normal gas constant?                                   
!!
integer :: u_gaussorder ! local variable to keep the order of gauss integration for moments
integer :: i,j,ci,p ! local counters 
real (DP) :: du     ! to keep temp calculations
real (DP), dimension (:), allocatable :: u_nodes, uumesh, uuweights ! meshes and gauss weights where the spectral representation will be evaluated 
real (DP), dimension (:), allocatable :: weighted_sol1, weighted_sol2 ! local variable where the solution is stored during the integration 
                                       ! weighted_sol1(m) m -- node number in u form uumesh(m) on U_{i} 
integer :: loc_alloc_stat ! to keep allocation status 

u_gaussorder=size(moments_u_gauss_nodes)

! check if the macroparamer arrays sizes are conforming 
if ((size(dens,1) .ne. size(x_nodes)) .or. (size(avel,1) .ne. size(x_nodes)) &
  .or. (size(tempr,1) .ne. size(x_nodes))) then 
   print *, "EvalMacroSolXmeshXnodes1D: Error! Sizes of the arrays (dens), (avel), (tempr) and (x_nodes) do not conform!"
   end if 
! allocate arrays for integration in u
allocate (uumesh(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodes1D: Allocation error for variable (uumesh)"
   end if 
allocate (u_nodes(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodes1D: Allocation error for variable (u_nodes)"
   end if 
allocate (uuweights(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodes1D: Allocation error for variable (uuweights)"
   end if 
!! allocate space for solution evaluated on (uumesh)
allocate (weighted_sol1(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodes1D: Allocation error for variable (weighted_sol2)"
   end if 
allocate (weighted_sol2(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodes1D: Allocation error for variable (weighted_sol1)"
   end if 
!!
! we prepare a larger array of weights (dublicates of gauss weights) to be used on subdivided cells
du = 2/Real(riu,DP) 
do ci=1,riu
 uuweights(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = Real (moments_u_gauss_weights, DP)
 u_nodes(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = -1 + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
end do 
! ready to calculate the moments:
!!
do j=1,size(f1,4)     !loop in integration cells in x
do p=1,size(x_nodes) ! loop in the nodes the cell in x (actually on [-1,1], not scaled)
dens(p,j)=0
avel(p,j)=0
tempr(p,j)=0
! 
do i=1,size(f1,3)      !loop in integration cells in u
   du = ( umesh(i) - umesh(i-1) )/Real(riu,DP) ! need this for the integration on subintervals
                            !!! IMPORTANT riu --moments_refine_u must be the same for all intervals in U_i
 do ci=1,riu
  uumesh(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = umesh(i-1) + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
 end do
   ! first we evaluate the solution at nodes (xxmesh,uumesh)
   weighted_sol1 = EvCellLegStdrd1D1D(x_nodes(p),u_nodes,k,s,max_deg,f1(:,:,i,j))
   weighted_sol2 = EvCellLegStdrd1D1D(x_nodes(p),u_nodes,k,s,max_deg,f2(:,:,i,j))
   ! we then multiply the solution by gauss weights in u:
   weighted_sol1 = weighted_sol1*uuweights
   weighted_sol2 = weighted_sol2*uuweights
   ! Now it is time to calculate the moments 
   ! the density:
   dens(p,j) = dens(p,j) + sum(weighted_sol1)*du
   ! the average velocity:
   weighted_sol1 = weighted_sol1*uumesh
   avel(p,j) = avel(p,j) + sum(weighted_sol1)*du
   ! the temperature 
   weighted_sol1 = weighted_sol1*uumesh
   tempr(p,j)=tempr(p,j) + sum(weighted_sol1+2*weighted_sol2)*du
   ! end calculating moments in cell [umesh(i-1), umesh(i)]
 end do  ! end loop in (i) --- cells in u
dens(p,j)=dens(p,j)/2  !!! 
avel(p,j)=avel(p,j)/2  
tempr(p,j)=tempr(p,j)/2
!! Now these were just the integrals, the actual macroparameters are obtained next
avel(p,j)=avel(p,j)/dens(p,j)  ! this takes care of the average velocity
tempr(p,j)=(tempr(p,j)/dens(p,j) - avel(p,j)*avel(p,j))/3/R ! this (presumably) takes care of the temperature
end do ! -- end loop in p --- gauss nodes on I_{j}
end do ! -- end loop in j -- the cell I_{j}
!
deallocate (u_nodes,uumesh,uuweights,weighted_sol1,weighted_sol2)
!
end subroutine EvalMacroSolXmeshXnodes1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalMacroSolXmeshXnodesStandard1D
!
! This is an accelerated version of an above subroutine. The acceleration (if any) is achieved by using 
! pre-calculated nodes and weights. Can only be used when the same order of gauss integration is used on all
! cells and when each cell is refined the same amount of times.
!
! This subroutine computes values of macroparameters. The evaluation is done on each I_{j} cell in x 
! (given by [xmesh]) ,  at all gauss nodes (given by [moments_x_gauss_nodes] and [moments_x_refine]).
! 
!  
! Integrals are calculated using the main mesh in the velocity variable (umesh). Nodes in u are based on 
! (moments_u_gauss_nodes) and (moments_u_refine) ((moments_u_gauss_nodes) are calculated using the parameters 
! (moments_u_gauss_order).  
!  
!
! EXPECTS: 
! dens(m,j)  --- array to store density
! avel(m,j)  --- array to store average velocity
! tempr(m,j) --- array to store temperature
!                (m,j) m --- the node on the interval I_{j}, j is the interval in x form xmesh
!                 Attention::
!                 m must have be the dimension (moments_x_gauss_order)x(moments_u_refine) 
!                 j must be size(xmesh)-1
! x_nodes(p) --- nodes on [-1,1] corresponding to points on each [x_{j-1/2},x_{j+1/2}] on which 
!              where the macroparamters need to be evaluated. 1<=p<=moments_x_gauss_nodes*moments_x_refine 
! Unodes(m,i) --- mesh points in the variable u used to evalue the integral in u for moments and the right side
!                 m --- number of the node (on the refined cell), i -- number of the cell in u
! Unodes_unit(m) --- nodes on interval [-1,1] combined for all refined subcells
! Uweights_unit(m) --- gauss weights on [-1,1] combined for all refined subcells   
! k,s,max_deg --- k is the degree of spectral approximation in x, s it the degree of polynomial approximation in u, and 
!                 max_deg -- is the maximum allowed degree in the spectral decomposition
! f1,f2 --- spectral coefficients! f12(p,m,i,j) -- p is the index in the basis functions in x 
!                                  ! -- m is the index in the basis functions in u
!                                   ! -- i is the cell in u
!                                   ! -- j is the cell in x
! R --- the gas constant
!
!
! RETURNS
! 
! dens, avel, tempr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine EvalMacroSolXmeshXnodesStandard1D(dens,avel,tempr,x_nodes,umesh,Unodes,Unodes_unit,Uweights_unit,riu,&
                                k,s,max_deg,f1,f2,R)
use spectral_tools_mod
!
real (DP), dimension (:,:), intent (out) :: dens, avel, tempr ! arrays to store density, average velocity and 
                                        ! temperature. E.G., dens(m,j) 
                                        ! m --- the node on the interval I_{j}, j --is the interval in x form xmesh
real (DP), dimension (:), intent (in) :: x_nodes 
                                        ! nodes in x on interval [-1,1] to be used in macroparameter evaluations. Essentially, 
                                        ! these correpond to the gauss nodes obtained after refining each an interval I_{j} of [xmesh]
                                        ! We assume that the same gauss nodes and the same refining coefficient are used 
                                        ! on all intervals. Since actual x does not show up, we will keep gauss nodes scaled on [-1,1] for speed
real (DP), dimension (0:), intent (in) :: umesh ! main mesh in variable u, where the spectral representation in defined
real (DP), dimension (:,:), intent (in) :: Unodes ! pre-set nodes on velocity cells (after refinement)
real (DP), dimension (:), intent (in) :: Unodes_unit, Uweights_unit 
                                                              ! gauss nodes and weights in u
                                                              ! to be used in evaluation of macoparameters 
integer (I4B), intent (in) :: riu ! these are coefficients of cell refinements in x and u
integer (I4B), intent (in) :: k,s,max_deg ! order in x, order in u and max degree of polynomial approximation 
real (DP), dimension (0:,0:,:,:), intent (in) :: f1,f2 ! f12(p,m,i,j) -- coefficients of spectral decomposition (converted from charactristic variables) 
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m is the index in the basis functions in u
                                   ! -- i is the cell in u
                                   ! -- j is the cell in x
real (DP), intent (in) :: R ! the normal gas constant?                                   
!!
integer :: i,j,ci,p ! local counters 
real (DP) :: du     ! to keep temp calculations
real (DP), dimension (:), allocatable :: weighted_sol1, weighted_sol2 ! local variable where the solution is stored during the integration 
                                       ! weighted_sol1(m) m -- node number in u form uumesh(m) on U_{i} 
integer :: loc_alloc_stat ! to keep allocation status 

! check if the macroparamer arrays sizes are conforming 
if ((size(dens,1) .ne. size(x_nodes)) .or. (size(avel,1) .ne. size(x_nodes)) &
  .or. (size(tempr,1) .ne. size(x_nodes))) then 
   print *, "EvalMacroSolXmeshXnodesStandard1D: Error! Sizes of the arrays (dens), (avel), (tempr) and (x_nodes) do not conform!"
   end if 
!! allocate space for solution evaluated on (uumesh)
allocate (weighted_sol1(1:size(Unodes_unit)), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodesStandard1D: Allocation error for variable (weighted_sol2)"
   end if 
allocate (weighted_sol2(1:size(Unodes_unit)), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "EvalMacroSolXmeshXnodesStandard1D: Allocation error for variable (weighted_sol1)"
   end if 
!!
do j=1,size(f1,4)     !loop in integration cells in x
do p=1,size(x_nodes) ! loop in the nodes the cell in x (actually on [-1,1], not scaled)
dens(p,j)=0
avel(p,j)=0
tempr(p,j)=0
! 
do i=1,size(f1,3)      !loop in integration cells in u
   du = ( umesh(i) - umesh(i-1) )/Real(riu,DP) ! need this for the integration on subintervals
   ! first we evaluate the solution at nodes (xxmesh,uumesh)
   weighted_sol1 = EvCellLegStdrd1D1D(x_nodes(p),Unodes_unit,k,s,max_deg,f1(:,:,i,j))
   weighted_sol2 = EvCellLegStdrd1D1D(x_nodes(p),Unodes_unit,k,s,max_deg,f2(:,:,i,j))
   ! we then multiply the solution by gauss weights in u:
   weighted_sol1 = weighted_sol1*Uweights_unit
   weighted_sol2 = weighted_sol2*Uweights_unit
   ! Now it is time to calculate the moments 
   ! the density:
   dens(p,j) = dens(p,j) + sum(weighted_sol1)*du
   ! the average velocity:
   weighted_sol1 = weighted_sol1*Unodes(:,i)
   avel(p,j) = avel(p,j) + sum(weighted_sol1)*du
   ! the temperature 
   weighted_sol1 = weighted_sol1*Unodes(:,i)
   tempr(p,j)=tempr(p,j) + sum(weighted_sol1+2*weighted_sol2)*du
   ! end calculating moments in cell [umesh(i-1), umesh(i)]
 end do  ! end loop in (i) --- cells in u
dens(p,j)=dens(p,j)/2  !!! 
avel(p,j)=avel(p,j)/2  
tempr(p,j)=tempr(p,j)/2
!! Now these were just the integrals, the actual macroparameters are obtained next
avel(p,j)=avel(p,j)/dens(p,j)  ! this takes care of the average velocity
tempr(p,j)=(tempr(p,j)/dens(p,j) - avel(p,j)*avel(p,j))/3/R ! this (presumably) takes care of the temperature
end do ! -- end loop in p --- gauss nodes on I_{j}
end do ! -- end loop in j -- the cell I_{j}
!
deallocate (weighted_sol1,weighted_sol2)
!
end subroutine EvalMacroSolXmeshXnodesStandard1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalMacroXmeshSolHOV1D
!
! This is an accelerated version of an above subroutine. The acceleration (if any) is achieved by using 
! pre-calculated nodes and weights. Can only be used when the same order of gauss integration is used on all
! cells and when each cell is refined the same amount of times.
!
!
! This subroutine computes values of macroparameters by callign the macroparamters integrator.
! Essentially, it just puts together a lot of parameters, converts the characteristic variables in 
! regular ones and then calls the EvalMacroSol1D1D. The evaluation is done at nodes 
! on each interval in (xmesh). Nodes in x are based on (moments_x_gauss_nodes) and (moments_x_refine) 
! and are calculated using the parameters (moments_x_gauss_order). For integration in the velocity, 
! (moments_u_gauss_nodes) and (moments_u_refine) are used. 
! 
! Highly dependent on the outside program
!
! Most of the variables are looked up in the Common_variables_mod
!
! do not call untill moments_x_gauss_nodes,moments_x_gauss_weights,&
!                                moments_u_gauss_nodes,moments_u_gauss_weights,& 
!                                moments_refine_x,moments_refine_u, etc are initialized 
!
! EXPECTS: 
! dens(m,j)  --- array to store density
! avel(m,j)  --- array to store average velocity
! tempr(m,j) --- array to store temperature
!                (m,j) m --- the node on the interval I_{j}, j is the interval in x form xmesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EvalMacroXmeshSolStandardHOV1D(dens,avel,tempr,ftil1,ftil2,x_nodes)
use common_variables_mod, only: k,s,max_deg,umesh,Aml,M,N,moments_refine_u,gasR,Unodes,Unodes_unit,Uweights_unit
!!
real (DP), dimension (:,:), intent (out) :: dens, avel, tempr ! arrays to store density, average velocity and 
                                        ! temperature. E.G., dens(m,j) 
                                        ! m --- the node on the interval I_{j}, j --is the interval in x form xmesh
real (DP), dimension (0:,0:,:,:), intent (in) :: ftil1,ftil2 ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,i,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- i is the cell in u
                                   ! -- j is the cell in x
!!! 
real (DP), dimension (0:k,0:s,M,N) :: f1,f2  ! auxiliary variable -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  f(p,m,i,j)( ftil(p,m,i,j), fti2(p,m,i,j):)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- i is the cell in u
                                                     ! -- j is the cell in x
integer  :: i,j,l,q,m_count    ! local counters
real (DP), dimension (:), intent (in) :: x_nodes ! nodes to be used in the evaluation of macraparameters
                                                 ! obtained by subdiving interval [-1,1] into (moments_refine_x)
                                                 ! parts and scaling using nodes (moments_x_gauss_nodes) on each part
                                                 ! these nodes will be passed to the standard evaluation procedure to 
                                                 ! bypass scaling [x_{j-1/2},x_{j+1/2}] to [-1,1] and back
!!!
! we convert the characteristic coefficients into non-characteristic, that is
! evaluates the coefficients of the spectral representation from the characteristic representation
   f1=0
   f2=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f1(q,l,i,j)= f1(q,l,i,j) + Aml(l,m_count,i)*ftil1(q,m_count,i,j)
     f2(q,l,i,j)= f2(q,l,i,j) + Aml(l,m_count,i)*ftil2(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
  !! now we call the procedure to evaluate the macroparameters. This one evaluates the macroparameters 
  !! based on the meshes (xmesh) and (umesh).
   call EvalMacroSolXmeshXnodesStandard1D(dens,avel,tempr,x_nodes,umesh,Unodes,Unodes_unit,Uweights_unit,&
                                moments_refine_u,k,s,max_deg,f1,f2,gasR)
end subroutine EvalMacroXmeshSolStandardHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetMomentsXnodesUnodesBasisFunHOV1D
!
! This procedure sets up nodes for integration of moments and for integration in the right side. 
! uses (xmesh) and (umesh) as the basis. Then each cell is divided into subcells using 
! moments_refine_x and moments_refine_u. On each new cell nodes of gaussian quadrature 
! moments_x_gauss_nodes and moments_u_guass_nodes are introduced (scaled to the appropriate interval). 
! Auxiliary variables (rix_moments_x_gauss_weights) (rix_unit_moments_x_gauss_nodes) and 
! (riu_moments_u_gauss_weights), (rix_unit_moments_u_gauss_nodes) introduced to keep nodes on [-1,1] 
! and weights
!
! Dependent on the main program! 
!
! Do not call before moments_x_gauss_nodes and moments_u_gauss_nodes are set up
! 
! Xnodes(p,j) p -- node, j -- cell in x
! Unodes(m,i) m -- node, i -- cell in u 
! XWeightLegPol(p,l) p -- nodes, l -- number of polynomial 
! UWeightLegPol(m,i) m -- nodes, i -- number of polynomial 
! Xnodes_unit(p) p -- nodes
! Unodes_unit(m) m -- nodes
! Xweights_unit(p) p -- nodes
! Uweights_unit(m) m -- nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMomentsXnodesUnodesBasisFunHOV1D
use common_variables_mod, only: Xnodes,Unodes,XWeightLegPol,UWeightLegPol,XLegPol,ULegPol,Xnodes_unit,Unodes_unit,&
                       Xweights_unit,Uweights_unit,moments_x_gauss_nodes,moments_x_gauss_weights,&
                       moments_u_gauss_nodes,moments_u_gauss_weights,moments_refine_x,moments_refine_u,M,N,s,k,&
                       umesh,xmesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: rix,riu ! these are coefficients of cell refinements in x and u
integer :: x_gaussorder, u_gaussorder ! local variable to keep the order of gauss integration for moments
real (DP) :: dx,du ! -- some dummy variables 
integer :: ci,l,i,p,j ! some counters ... 
integer :: loc_alloc_stat ! to keep allocation status
!
riu=moments_refine_u
rix=moments_refine_x
!!!
u_gaussorder=size(moments_u_gauss_nodes)
x_gaussorder=size(moments_x_gauss_nodes)
!
! allocate arrays for integration in u 
allocate (Unodes(1:u_gaussorder*riu,M), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Unodes)"
   end if 
allocate (Unodes_unit(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Unodes_unit)"
   end if 
allocate (Uweights_unit(1:u_gaussorder*riu), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Uweights_unit)"
   end if 
! allocate arrays for integration in x
allocate (Xnodes(1:x_gaussorder*rix,N), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Xnodes)"
   end if 
allocate (Xnodes_unit(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Xnodes_unit)"
   end if 
allocate (Xweights_unit(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (Xweights_unit)"
   end if 
! now space for leg polynomials: 
allocate (XWeightLegPol(1:x_gaussorder*rix,0:k), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (XWeightLegPol)"
   end if 
allocate (UWeightLegPol(1:u_gaussorder*riu,0:s), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (UWeightLegPol)"
   end if 
   ! now space for leg polynomials: 
allocate (XLegPol(1:x_gaussorder*rix,0:k), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (XLegPol)"
   end if 
allocate (ULegPol(1:u_gaussorder*riu,0:s), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesUnodesBasisFunHOV1D: Allocation error for variable (ULegPol)"
   end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
!! Now we will set up these arrays: 
du = 2/Real(riu,DP) 
do ci=1,riu
 Uweights_unit(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = Real (moments_u_gauss_weights, DP)
 Unodes_unit(1+(ci-1)*u_gaussorder:ci*u_gaussorder) = -1 + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
end do 
do l=0,s
 ULegPol(:,l) = EvalHorner(LegendrePoly (l), Unodes_unit)
 UWeightLegPol(:,l) = ULegPol(:,l)*Uweights_unit
end do 
do i=1,M      !loop in integration cells in u
 du = ( umesh(i) - umesh(i-1) )/Real(riu,DP) ! need this for the integration on subintervals
              !!! IMPORTANT riu --moments_refine_u must be the same for all intervals in U_i
 do ci=1,riu
  Unodes(1+(ci-1)*u_gaussorder:ci*u_gaussorder,i) = umesh(i-1) + du*(Real(ci,DP) - .5) + Real(moments_u_gauss_nodes,DP)*du/2 
 end do
end do       !end loop in cells in u
!!
dx = 2/Real(rix,DP) 
do ci=1,rix
 Xweights_unit(1+(ci-1)*x_gaussorder:ci*x_gaussorder) = Real (moments_x_gauss_weights, DP)
 Xnodes_unit(1+(ci-1)*x_gaussorder:ci*x_gaussorder) = -1 + dx*(Real(ci,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
end do 
do p=0,k
 XLegPol(:,p) = EvalHorner(LegendrePoly (p), Xnodes_unit)
 XWeightLegPol(:,p) = XLegPol(:,p)*Xweights_unit
end do 
do j=1,N      !loop in integration cells in x
 dx = ( xmesh(j) - xmesh(j-1) )/Real(rix,DP) ! need this for the integration on subintervals
              !!! IMPORTANT rix --moments_refine_x must be the same for all intervals in X_j
 do ci=1,rix
  Xnodes(1+(ci-1)*x_gaussorder:ci*x_gaussorder,j) = xmesh(j-1) + dx*(Real(ci,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
 end do
end do       !end loop in cells in x
!!! Kinda done...
end subroutine SetMomentsXnodesUnodesBasisFunHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetMomentsXnodesBasisFunHOV1D
!
! This is a copy of the above subroutine with the structures dealing with discretization in the velocity variable deleted.  
! All relevan objects for the velocity discretization can be accessed via DGVlib_commvar  
!
! This procedure sets up nodes for integration of moments and for integration in the right side. 
! uses (xmesh)  as the basis. Then each cell is divided into subcells using 
! moments_refine_x. On each new cell nodes of gaussian quadrature 
! moments_x_gauss_nodes are introduced (scaled to the appropriate interval). 
! Auxiliary variables (rix_moments_x_gauss_weights) (rix_unit_moments_x_gauss_nodes) are introduced to keep nodes on [-1,1] 
! and weights
!
! Dependent on the main program! 
!
! Do not call before moments_x_gauss_nodes and moments_u_gauss_nodes are set up
! 
! Xnodes(p,j) p -- node, j -- cell in x
! XWeightLegPol(p,l) p -- nodes, l -- number of polynomial 
! Xnodes_unit(p) p -- nodes
! Xweights_unit(p) p -- nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMomentsXnodesBasisFunHOV1D
use common_variables_mod, only: Xnodes,XWeightLegPol,XLegPol,Xnodes_unit,&
                       Xweights_unit,moments_x_gauss_nodes,moments_x_gauss_weights, &
                       moments_refine_x,N,k,xmesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: rix ! these are coefficients of cell refinements in x 
integer :: x_gaussorder ! local variable to keep the order of gauss integration for moments
real (DP) :: dx ! -- some dummy variables 
integer :: ci,l,i,p,j ! some counters ... 
integer :: loc_alloc_stat ! to keep allocation status
!
rix=moments_refine_x
!!!
x_gaussorder=size(moments_x_gauss_nodes)
!
! allocate arrays for integration in x
allocate (Xnodes(1:x_gaussorder*rix,N), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesBasisFunHOV1D: Allocation error for variable (Xnodes)"
   end if 
allocate (Xnodes_unit(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesBasisFunHOV1D: Allocation error for variable (Xnodes_unit)"
   end if 
allocate (Xweights_unit(1:x_gaussorder*rix), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesBasisFunHOV1D: Allocation error for variable (Xweights_unit)"
   end if 
! now space for leg polynomials: 
allocate (XWeightLegPol(1:x_gaussorder*rix,0:k), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesBasisFunHOV1D: Allocation error for variable (XWeightLegPol)"
   end if 
   ! now space for leg polynomials: 
allocate (XLegPol(1:x_gaussorder*rix,0:k), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SetMomentsXnodesBasisFunHOV1D: Allocation error for variable (XLegPol)"
   end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                       
!! Now we will set up these arrays: 
dx = 2/Real(rix,DP) 
do ci=1,rix
 Xweights_unit(1+(ci-1)*x_gaussorder:ci*x_gaussorder) = Real (moments_x_gauss_weights, DP)
 Xnodes_unit(1+(ci-1)*x_gaussorder:ci*x_gaussorder) = -1 + dx*(Real(ci,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
end do 
do p=0,k
 XLegPol(:,p) = EvalHorner(LegendrePoly (p), Xnodes_unit)
 XWeightLegPol(:,p) = XLegPol(:,p)*Xweights_unit
end do 
do j=1,N      !loop in integration cells in x
 dx = ( xmesh(j) - xmesh(j-1) )/Real(rix,DP) ! need this for the integration on subintervals
              !!! IMPORTANT rix --moments_refine_x must be the same for all intervals in X_j
 do ci=1,rix
  Xnodes(1+(ci-1)*x_gaussorder:ci*x_gaussorder,j) = xmesh(j-1) + dx*(Real(ci,DP) - .5) + Real(moments_x_gauss_nodes,DP)*dx/2 
 end do
end do       !end loop in cells in x
!!! Kinda done...
end subroutine SetMomentsXnodesBasisFunHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  added Alex 05/24/2010

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set1DHOVMomentsTimeArraysGaussNodesWeights
!
! This subroutine sets up the arrays that will keep gauss nodes, gauss weights for the 
! evaluation of the moments in the solution. The moments are computed by evaluationg the 
! spectral solution at these nodes in both variables.
!
! Also, this subroutine allocates the arrays where the moments of the spectral 
! solution and the exact solution and the times at which these moments were calculate will be stored. 
!! moments_time_sol (i,j,k) 
!! i ---1,2,3 the first three moments, 4 -- time
!! j = 0,1,...  --- # of the moment of saving. 
!! k = 1,2 --- for each component of the solution
!! moments_x_gauss_order, moments_u_gauss_order --- number of nodes in the gaussian formula(not the methods order)!
! varaibles related to generation of non-unifom meshes: 
!
! Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters moments_x_gauss_order, moments_u_gauss_order 
! num_eval_error are selected!! 
! 
! The number of evaluations for the moments is set equal to the number ot evaluations for the error.
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Set1DHOVMomentsTimeArraysGaussNodesWeights
use common_variables_mod, only: moments_x_gauss_order, moments_u_gauss_order, &
                                moments_x_gauss_nodes, moments_u_gauss_nodes, &
                                moments_x_gauss_weights, moments_u_gauss_weights, &
                                num_eval_error, moments_time_sol, moments_time_exact, &
                                moments_refine_x
                                 
                                ! errays to keep evaluated total moments in the sectral solution and in the exact solution 
														  ! moments_sol (i,j,k) 
                                                          ! i --- ## of moment
                                                          ! j = 0,1,...  --- # of the moment of saving. 
                                                          ! k = 1,2 --- for each component of the solution
                                                          
! varaibles related to generation of non-unifom meshes: 
use gaussian_mod
!!!
integer :: loc_alloc_stat ! to keep allocation status
integer :: i ! local counter
!!! 
!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights for the inegration in the x variable
  moments_x_gauss_order = max(moments_x_gauss_order,1) 
  !!!!! Next line is a diagnostic. Comment when not needed
  ! print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Set moments_x_gauss_order=", moments_x_gauss_order
  !!!!!!! End diagnostic
     !
  allocate (moments_x_gauss_nodes(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_x_gauss_nodes)"
     end if 
     !
  allocate (moments_x_gauss_weights(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_x_gauss_weights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (moments_x_gauss_order == 1) then 
   moments_x_gauss_nodes(:)= (/ 0.0_DP /)
   moments_x_gauss_weights(:)= (/ 2.0_DP /)
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_x_gauss_nodes,moments_x_gauss_weights)  ! the number of nodes/precision of the quadrature depends on the size of gnodes,gweights
  end if
!!! Prepare the Gaussian nodes and weights in the second variable (u) !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the second variable
  moments_u_gauss_order = max(moments_u_gauss_order,1) 
  !!!!! Next line is a diagnostic. Comment when not needed
  ! print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Set moments_u_gauss_order=", moments_u_gauss_order
  !!!!! End diagnostic
     !
  allocate (moments_u_gauss_nodes(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_u_gauss_nodes)"
     end if 
     !
  allocate (moments_u_gauss_weights(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_u_gauss_weights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (moments_u_gauss_order == 1) then 
   moments_u_gauss_nodes(:)= (/ 0.0_DP /)
   moments_u_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_u_gauss_nodes,moments_u_gauss_weights)  ! the number of nodes/precision of the quadrature depends on the size of gnodes,gweights
  end if
!!! Now it is time to set up the errays to store the moments: 
allocate (moments_time_sol(20,0:num_eval_error+1,2), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_sol)"
     end if 
     !
allocate (moments_time_exact(4,0:num_eval_error+1,2), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMomentsTimeArraysGaussNodesWeights: Allocation error for variable (moments_exact)"
     end if 
     !
end subroutine Set1DHOVMomentsTimeArraysGaussNodesWeights


!!!!!!!! end added Alex 05/24/2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sol_integr_mod

