!
!  boundary_cond_mod.f90
!
!  This module contains functions implementing diferent boundary conditions for 1D HOV model.
!
!  
!  7/29/2008 5:18:04 PM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module boundary_cond_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PeriodicBCsLeg
!
!  This subroutine implements periodic BCs. Essentially, it takes components of the spectral decomposition, 
!  and assembles it in the boundary points. Then boundary data at the rigth boundary is set to the computed value at the left 
!  boundary and vise versa. The suffix "Leg" in the name means that this subrouting is for the decomposition using Legendre Scaled 
!  polynomials.  
subroutine PeriodicBCsLeg(TilFmiLeft,TilFmiRight,xmesh,ftil) 
use spectral_tools_mod 
! main variables 
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in varaible
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                 ! characteristic variables ftil(p,m,i,j):
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m is the index in the basis functions in u
                                 ! -- j is the cell in x
real (DP), dimension (size(ftil, 2)), intent (out) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{m,i}{t,x} 
											    ! at the right and left boundary points, TilFmiRight(m,i), TilFmiLeft(m,i)
                                                ! -- m is the index in the basis functions in u
                                                ! -- i is the cell in u
! supplementary variables                                  
integer (I4B) :: i,N ! i,j are counters, N is to store the number of cells in x variable 
!!!						
N=size(ftil,3)
!!! begin loops in i and m
do i=1,size(ftil,2)
  TilFmiRight(i) = EvCellLeg1D(xmesh(0),xmesh(0),xmesh(1),ftil(:,i,1))
  TilFmiLeft(i) = EvCellLeg1D(xmesh(N),xmesh(N-1),xmesh(N),ftil(:,i,N))
end do                                  
!!! end loops in i and j                                
end subroutine PeriodicBCsLeg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! added 09/07/08 Alex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ExactBCs1D3DLeg
!
! This subroutine implements the exact boundary conditions. When boundary data is prescribed by the known function 
! AND this function is actually known for all x and all u and all t, so that a 1d+1d+t dimensional function is given 
! this procedure will evaluate the spectral decomposition of the function in the variable u 
! at the left and right boundaries in the variable $x$. The basis functions are scaled Legenre's polynomials. 
! The suffix "Leg" in the name means that this subroutine is for the decomposition using Legendre Scaled 
! polynomials.
!
subroutine ExactBCs1D3DLegOrds (TilFmiLeft,TilFmiRight,xmesh,t,f,nodes_u,nodes_v,nodes_w) 
use spectral_tools_mod 
intrinsic MatMul
! main variables 
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in variable x
real (DP), intent (in) :: t  ! the time varaible -- will go into the function
real (DP), dimension (:), intent (out) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{m}{t,x} 
											    ! at the right and left boundary points, TilFmiRight(m), TilFmiLeft(m)
                                                ! -- m is the index in the basis functions in u
real (DP), dimension (:), intent (in) :: nodes_u, nodes_v, nodes_w ! velocity nodes at which the BC are evaluated supplied to the subroutine

                                                
! supplementary variables                                  
integer (I4B) :: i,N ! i is a counter, N is to store the number of cells in x variable 
 ! interface block for the dummy function f: (x-scalar, u-vector)
 ! first variable must be scalar and second must be vector! 
 interface 
   function f (x, u,v,w, t) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u,v,w     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: t ! the time variable 
    real (DP), dimension (size(u)) :: y ! values of the function  
   end function f
 end interface
!!! 
N=size(xmesh)-1
! first we evaluate coefficients of the spectral decomposition, not in the characteristic variables,					  
TilFmiRight = f(xmesh(N),nodes_u,nodes_v,nodes_w,t)
TilFmiLeft = f(xmesh(0),nodes_u,nodes_v,nodes_w,t) 
end subroutine ExactBCs1D3DLegOrds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ExactBCs1D3DLeg
!
! This subroutine implements the exact boundary conditions. When boundary data is prescribed by the known function 
! AND this function is actually known for all x and all u and all t, so that a 1d+1d+t dimensional function is given 
! this procedure will evaluate the spectral decomposition of the function in the variable u 
! at the left and right boundaries in the variable $x$. The basis functions are scaled Legenre's polynomials. 
! The suffix "Leg" in the name means that this subroutine is for the decomposition using Legendre Scaled 
! polynomials.
!
subroutine ExactBCs1D3DLeg (TilFmiLeft,TilFmiRight,xmesh,t,f) 
use spectral_tools_mod 
use DGV_commvar, only: nodes_u, nodes_v, nodes_w
intrinsic MatMul
! main variables 
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in variable x
real (DP), intent (in) :: t  ! the time varaible -- will go into the function
real (DP), dimension (:), intent (out) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{m}{t,x} 
											    ! at the right and left boundary points, TilFmiRight(m), TilFmiLeft(m)
                                                ! -- m is the index in the basis functions in u
                                                
! supplementary variables                                  
integer (I4B) :: i,N ! i is a counter, N is to store the number of cells in x variable 
 ! interface block for the dummy function f: (x-scalar, u-vector)
 ! first variable must be scalar and second must be vector! 
 interface 
   function f (x, u,v,w, t) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u,v,w     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: t ! the time variable 
    real (DP), dimension (size(u)) :: y ! values of the function  
   end function f
 end interface
!!! 
N=size(xmesh)-1
! first we evaluate coefficients of the spectral decomposition, not in the characteristic variables,					  
TilFmiRight = f(xmesh(N),nodes_u,nodes_v,nodes_w,t)
TilFmiLeft = f(xmesh(0),nodes_u,nodes_v,nodes_w,t) 
end subroutine ExactBCs1D3DLeg


!!!!!!!!!!!!! end added 09/07/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Added Alex 05/17/09  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DiffuiveBCsLeft1D3DLeg (TilFmi1Left,xmesh(0),xmesh(1),max_deg,f1,&
!                            nodes_u,nodes_v,nodes_w,nodes_gwts,u_w, v_w, w_w, T_w,R)  
!
!  This procedure implements the diffusive boundary condtions on the LEFT boundary.
!  Namely, it takes the discrete solution in the most left spatial cell and 
!  calculates the flux of the paticles on the left boundary. To calculate the flux we find cells that 
!  only have incoming velocities and one (possible) cell that has both incoming and outgoing velocities.
!  The part of the flux that corresponds to the incoming only cells is computed using fast procedures, 
!  The cell that has both incoming and outgoing velocities is computed independently (however, there should 
!  not be such cell if non-uniform mesh in u of type =2 is used. This type of mesh will have velocity of the wall u_{w}= 0 
!  as a mesh point). Once the flux is computed, Maxwellian equilibrium distribution with given temperature T_{w} and 
!  velocity u_{w}=0 is used to compute the incoming fields. For the incoming fields, again, the "incoming only" cell 
!  are treated separately and the cell with the mixed velocities is special
! 
!  Attn: works with spectral coefficients, not characteristic variables
! 
! Expects: 
! 
! Fmi1Left, Fmi1Left ! the values of the spectral decomposition of the solution in the variable u
!                    ! 
!  xl,xr     ---  values of the most left interval in x: [xl, xr]
!  umesh     ---  meshpoints in u:
!  max_deg   ---  max_degree of the polynomial interpolation
!  
! 
! f  ---     ! f(p,m,i) is the coefficients of the spectral representation in the most left cell 
!                                 ! -- coefficients of spectral decomposition in 
!                                 ! characteristic variables CORRESPONDING to the most left 
!                                 ! cell --- ftil(p,m,i,1):
!                                 ! -- p is the index in the basis functions in x 
!                                 ! -- m is the index in the basis functions in u
!                                 ! -- i is the cell in u
!                                 ! -- j = 1 is the cell in x
! u_gauss_nodes,u_wleg_pol --- values of the gauss nodes and the values of the legendre's polynomials 
! u_w  --- velocity of the wall (right now only u_w=0 is implemented)
! T_w  --- temperature of the wall
! R is the gas constant
! Unodes_unit,Uweights_unit,moments_refine_u --- these contain nodes for moment integration on [-1,1] after 
!                  ! a refinment of the interval (moments_refine_u) times, dublicated gauss weights corresponding to 
!                  ! a refined interval, and the coefficient for the refinement
!

!subroutine DiffuiveBCsLeft1D3DLeg (TilFmi1Left,xmesh(0),xmesh(1),max_deg,f1,&
!                            nodes_u,nodes_v,nodes_w,nodes_gwts,u_w, v_w, w_w, T_w,R) 
!use spectral_tools_mod
!use distrib_mod 
!! intrinsic MatMul
!real (DP), intent (in) :: xl,xr ! points of the most left cell in x = [xmesh(0),xmesh(1)]
!integer (I4B), intent (in) :: max_deg ! this is the constant defining the maximum possible degree of approximating polynomials
!real (DP), dimension (0:,:), intent (in) :: f1 ! f(p,m) is the main variable 
!                                 ! -- coefficients of spectral decomposition in 
!                                 ! characteristic variables CORRESPONDING to the most left 
!                                 ! cell --- ftil(p,m,1):
!                                 ! -- p is the index in the basis functions in x 
!                                 ! -- m is the index in the basis functions in u
!                                 ! -- j = 1 is the cell in x
!real (DP), dimension (:), intent (in) :: nodes_u, nodes_v, nodes_w, nodes_gwts ! nodes and weights of the gauss integration formula to be used in evaluation
!                                  ! of moments in variable u
!real (DP), intent (in) :: u_w, v_w, w_w, T_w, R ! parameters of the diffusive BCs -- velocity and temperatture of the wall                                                         
!real (DP), dimension (0:,:), intent (out) :: Fmi1Left,Fmi2Left ! values of the f_{m}{t,x} 
!											    ! at the right and left boundary points, FmiLeft(m,i)
!                                                ! -- m is the index in the basis functions in u
!! supplementary variables
!logical :: middle_cell_exists_flag ! This variable is =true if the wall velocity is not one of the mesh points, 
!                                   ! but instead u_w falls in one of the cells in variable u
!integer (I4B) :: k,M ! k -- degree of polynomial approx in x
!                                       ! s --- degree of polynomial approx in u, M -- number of cells in u
!integer :: i,l,q,m_count,m_cell_i,loc_alloc_stat ! Some counters and integer variables                             
!real (DP) :: n_w, h, du ! some dump variables to keep temp calculations; n_w -- will the the density of the 
!                           ! the gas reflected from the wall
!real (DP), dimension (:), allocatable :: sol ! local variable where the solution is stored during the integration 
!                                       ! weighted_sol(m) -- is value on one cell, m -- node number in u form uumesh(m) on U_{i} 
!real (DP), dimension (:), allocatable :: Unodes_mcell ! -- will be used to evaluate integrals in u in the "median" cell in u                                      
!real (DP), dimension (:,:), allocatable :: legpol_mcell ! legpol+mcell(size(x),0:k) k-- degree of the max polynomial. this will be needed to evaluate leg polynomials to deal with the middle cell
!real (DP) :: mol_flux, beta1, beta2 ! temporary variable to keep the value of the molecular flux and 2\pi R T quantities
!
!!!!!!!!!!!!!
!! will need these guys for brevity
!k=size(f1,1)-1 
!m=
! First we scan the mesh in u (umesh) to determine the boundary cell (if any) and the place where the 
! boundary cell is. 
! check if the velocity of the wall u_w is on the velocity mesh
!if ((u_w < umesh(0)) .or. (u_w > umesh(M))) then 
!  print *, "DiffuiveBCsLeft1D1DLeg: Velocity of the wall u_w is outsied the mesh in u"!
!		    stop
!end if   
! now let us find the cell where u_w is: (It is not optimal to search of it on every time step, but don't know how to do it better
!i=1
!middle_cell_exists_flag = .FALSE.
!do while ((i <= M-1) .and. (umesh(i) <= u_w))
!  i=i+1
!end do 
!m_cell_i = i
!du=(umesh(m_cell_i) - umesh(m_cell_i-1))/1024.0_DP
!if (umesh(m_cell_i)-u_w < du ) then
!  middle_cell_exists_flag = .FALSE.
!  !! Let's add a warning if the "reflected cells" are not present:
!    if (m_cell_i >= M) then
!    print *, "DiffuiveBCsLeft1D1DLeg: Warning: no reflected cells on the left boundary found in umesh."
!    end if
!  !! end warning
!else 
!  if ( u_w-umesh(m_cell_i-1) < du ) then 
!     middle_cell_exists_flag = .FALSE.
!     m_cell_i=m_cell_i-1
!  else
!     middle_cell_exists_flag = .TRUE.
!  end if 
!end if 
! Now that we know where the "incident cells in u" are (everything before meshpoint umesh(i)),
! it is time to calculate the flux: 
! first, let us prepare a variable to store intermediate result:
!allocate (sol(1:size(Unodes_unit)), stat=loc_alloc_stat)
!    if (loc_alloc_stat >0) then 
!    print *, "DiffuiveBCsLeft1D1DLeg: Allocation error for variable (sol)"
!    end if
! 
!mol_flux=0.0_DP
! first let us evaluate the portion of the flux that is due to the "strictly incident" cells.
!do i=1,m_cell_i-1      !loop in strictrly "incident" cells in u
!   du = (umesh(i) - umesh(i-1))/Real(moments_refine_u,DP) ! need this for the integration on subintervals
!   ! first we evaluate the solution at nodes (xxmesh,uumesh)
!   sol = EvCellLegStdrd1D1D_u_vector(-1.0_DP,Unodes_unit,k,s,max_deg,f1(:,:,i))
!   ! we then multiply the solution by gauss weights in u:
!   sol = sol*Uweights_unit*Unodes(:,i)
!   ! Now it is time to calculate the flux 
!   mol_flux = mol_flux + sum(sol)*du
!   ! end calculating flux in the cell [umesh(i-1), umesh(i)]
!end do  ! end loop in (i) --- "incident" cells in u
! now we evaluate the portion of the flux for the "middle cell"
!if (middle_cell_exists_flag .eqv. .TRUE.) then 
!   du = ( u_w - umesh(m_cell_i-1) ) / Real(moments_refine_u, DP) 
!   !! allocate space for nodes in u for integration in the "middle cell"
!   allocate (Unodes_mcell(1:size(Unodes_unit)), stat=loc_alloc_stat)
!    if (loc_alloc_stat >0) then 
!    print *, "DiffuiveBCsLeft1D1DLeg: Allocation error for variable (Unodes_mcell)"
!    end if
!   !! end allocate 
!   Unodes_mcell = (u_w + umesh(m_cell_i-1))/2.0_DP + Unodes_unit*du/2.0_DP
!   sol = EvCellLeg1D1D_u_vector(xl,Unodes_mcell,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f1(:,:,m_cell_i))
!   ! we then multiply the solution by gauss weights in u:
!   sol = sol*Uweights_unit*Unodes_mcell
!   ! Now it is time to calculate the flux 
!   mol_flux = mol_flux + sum(sol)*du
!   ! end calculating flux in the middle cell 
!   deallocate (Unodes_mcell)
!else
! treat the cell right before m_cell_i as an incident cell
!   du = (umesh(m_cell_i) - umesh(m_cell_i-1))/Real(moments_refine_u,DP) ! need this for the integration on subintervals
!   ! first we evaluate the solution at nodes (xxmesh,uumesh)
!   sol = EvCellLegStdrd1D1D_u_vector(-1.0_DP,Unodes_unit,k,s,max_deg,f1(:,:,m_cell_i))
!   ! we then multiply the solution by gauss weights in u:
!   sol = sol*Uweights_unit*Unodes(:,m_cell_i)
!   ! Now it is time to calculate the flux 
!   mol_flux = mol_flux + sum(sol)*du
!   ! end calculating flux in the cell [umesh(i-1), umesh(i)]
!end if
!mol_flux = mol_flux /2.0_DP  !!! comes from integration
!deallocate (sol)
! We have calculated the incident flux. Next we need to evaluate the density of the reflected 
! Maxwellian distribution
!n_w = -mol_flux/( sqrt(R*T_w/pi25DT/2.0_DP) + u_w/2.0_DP )
!beta1 = n_w/sqrt(2.0_DP*pi25DT*R*T_w)          ! we will need beta1 and beta2 to calculate distribution of 
!beta2 = n_w*sqrt(R*T_w/2.0_DP/pi25DT)          ! reflected particles
!Fmi1Left=0
!Fmi2Left=0
! first, let us prepare a variable to store intermediate result:
!    allocate (sol(1:size(u_gauss_nodes)), stat=loc_alloc_stat)
!    if (loc_alloc_stat >0) then 
!    print *, "DiffuiveBCsLeft1D1DLeg: Allocation error for variable (sol)"
!    end if
    
! now we can set up the solution values on the left boundary
! We evaluate coefficients of the spectral decomposition, not in the characteristic variables,					  
! now we go through the "strictly reflecfted" cells
!do i=m_cell_i+1,M
!  !!! prepare integration with the cell
!  u = (umesh(i) - umesh(i-1))*u_gauss_nodes/2.0_DP + (umesh(i) + umesh(i-1))/2.0_DP
!  !!! Integration within the cell [umesh(i-1),umesh(i)]
!  sol = ExponMaxwell1D1D(u, u_w, T_w, R)
!  do l=0,s
!      Fmi1Left(l,i) = sum(sol*u_wleg_pol(l,:))*(2.0_DP*Real(l,DP)+1.0_DP)/2.0_DP ! we only need one right now, the second one can be calculated from the first one (see below)
!  end do 
! To get the correct distribution, need to multiply the maxwell's exponent by beta1 and beta2
!Fmi2Left(:,i) = Fmi1Left(:,i)*beta2
!Fmi1Left(:,i) = Fmi1Left(:,i)*beta1
!end do
! now we need to do the middle cell, if it exists... 
!if (middle_cell_exists_flag .eqv. .TRUE.) then 
!  !! allocate space for nodes in values of the legendre's polynomials for integration in the "middle cell"
!    allocate (legpol_mcell(size(u_gauss_nodes),0:s), stat=loc_alloc_stat)
!    if (loc_alloc_stat >0) then 
!    print *, "DiffuiveBCsLeft1D1DLeg: Allocation error for variable (legpol_mcell)"
!    end if
!  !! end allocate 
!  !!! First let us deal with the "incident" part of the cell... 
!  !!! Integration within the cell [umesh(m_cell_i-1),umesh(m_cell_i)] -- only "the incident" portion
!  !!! first, compute "incident" velocity nodes
!    ! Component #1:
!    u = (u_w - umesh(m_cell_i-1))*u_gauss_nodes/Real(2,DP) + (u_w + umesh(m_cell_i-1))/Real(2,DP)
!    sol = EvCellLeg1D1D_u_vector(xl,u,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f1(:,:,m_cell_i))*u_gauss_weights
!  !!! now we need nodes on[-1,1] that correspond to the values of the basis function 
!    u = Real(2,DP)*( u - (umesh(m_cell_i) + umesh(m_cell_i-1))/Real(2,DP) )/(umesh(m_cell_i) - umesh(m_cell_i-1))
!    legpol_mcell = EvalLegPol(s,u)
!    do l=0,s
!      Fmi1Left(l,m_cell_i) = sum(sol*legpol_mcell(:,l))*(2*l+1)/Real(2,DP)*((u_w - umesh(m_cell_i-1))/(umesh(m_cell_i) - umesh(m_cell_i-1)))
!    end do 
!    !  Component #2: 
!    u = (u_w-umesh(m_cell_i-1))*u_gauss_nodes/Real(2,DP) + (u_w + umesh(m_cell_i-1))/Real(2,DP)
!    sol = EvCellLeg1D1D_u_vector(xl,u,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f2(:,:,m_cell_i))*u_gauss_weights
!  !!! now we need nodes on [-1,1] that correspond to the values of the basis function 
!    do l=0,s
!      Fmi2Left(l,m_cell_i) = sum(sol*legpol_mcell(:,l))*(2*l+1)/Real(2,DP)*((u_w - umesh(m_cell_i-1))/(umesh(m_cell_i) - umesh(m_cell_i-1)))
!    end do 
!  !!! Now let us do the "reflected part:
!  !!! Components #1 and #2: 
!  !!! prepare the "reflected nodes" 
!    u = (umesh(m_cell_i)-u_w)*u_gauss_nodes/Real(2,DP) + (umesh(m_cell_i) + u_w)/Real(2,DP)
!  !!! Integration within the cell [umesh(m_cell_i-1),umesh(m_cell_i)] -- only "the reflected" portion
!    sol = ExponMaxwell1D1D(u, u_w, T_w, R)*u_gauss_weights 
!  !!! now the nodes on [-1,1] for the basis function  
!    u = Real(2,DP)*( u - (umesh(m_cell_i) + umesh(m_cell_i-1) )/Real(2,DP) )/(umesh(m_cell_i) - umesh(m_cell_i-1))
!    legpol_mcell = EvalLegPol(s,u)
!    do l=0,s
!      du = sum(sol*legpol_mcell(:,l))*(2*l+1)/Real(2,DP)*( (umesh(m_cell_i)-u_w)/(umesh(m_cell_i) - umesh(m_cell_i-1)) )
!      Fmi1Left(l,m_cell_i) = Fmi1Left(l,m_cell_i) + du*beta1
!      Fmi2Left(l,m_cell_i) = Fmi2Left(l,m_cell_i) + du*beta2
!    end do 
!  !!! end integration within the boundary cell 
!    deallocate (legpol_mcell)
!else 
!   ! treat the cell right before m_cell_i as an incident cell
!   do l = 0,s
!   do q = 0,min(max_deg - l,k)
!   Fmi1Left(l,m_cell_i) = Fmi1Left(l,m_cell_i) + ((-1)**q)*f1(q,l,m_cell_i)
!   Fmi2Left(l,m_cell_i) = Fmi2Left(l,m_cell_i) + ((-1)**q)*f2(q,l,m_cell_i)
!   end do 
!   end do 
!end if  
!deallocate (sol)
!!! Now let us calculate the boundary value at the strictly incident cell   
!!! Actually, we simply assemble the corresponding coefficients of the incident varaibles... 
!do i=1,m_cell_i-1
!   do l = 0,s
!   do q = 0,min(max_deg - l,k)
!   Fmi1Left(l,i) = Fmi1Left(l,i) + ((-1)**q)*f1(q,l,i)
!   Fmi2Left(l,i) = Fmi2Left(l,i) + ((-1)**q)*f2(q,l,i)
!   end do 
!   end do 
!end do    
!
!end subroutine DiffuiveBCsLeft1D1DLeg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DiffuiveBCsLeftHOV1D (TilFmi1Left,TilFmi2Left,ftil1_left,ftil2_left)
!
!  This procedure is a shell to call the function that implements the diffusive boundary condtions 
!  on the LEFT boundary. Specifically, this procedure converts characteristic variables into spectral coefficients. 
!  Then pulls together a lot parameters and call the subroutin that implements the diffusive conditions
!  Then the result is converted to characteristic variables
! 
!  This procedue looks up a lot of variables in the "common_variables_mod"
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine DiffuiveBCsLeftHOV1D (TilFmi1Left,ftil1_left)

!use common_variables_mod, only: xmesh,max_deg,&
!                                   T_w_left,gasR

!use DGV_commvar, only: nodes_u, nodes_v, nodes_w, nodes_gwts
!
!real (DP), dimension (:), intent (out) :: TilFmi1Left ! values of the f_{m}{t,x} 
											    ! at the right and left boundary points, FmiLeft(m)
                                                ! -- m is the index in the basis functions in u
                                               
                                                
!real (DP), dimension (0:,:), intent (in) :: ftil1_left ! f(p,m) -- auxiliary variable -- 
                                       ! coefficients of spectral decomposition 
                                       ! -- p is the index in the basis functions in x 
                                       ! -- m is the index in the basis functions in u
                                       ! -- these coefficients correspond to j=1 (first cell in x) in the ftil(p,m,j) 


                                                
! supplementary variables
!real (DP), dimension (0:k,size(nodes_u,1)) :: f1,f2 ! f(p,m,i) -- auxiliary variable -- 
                                       ! coefficients of spectral decomposition 
                                       ! we will calculate them from the characteristic variables  
                                       !  f(p,m,i) <= ftil(p,m,i)
                                       ! -- p is the index in the basis functions in x 
                                       ! -- m is the index in the basis functions in u
                                       ! -- these coefficients correspond to j=1 (first cell in x) in the ftil(p,m,j) 

! now ftil1_left(q,m) contains the coefficients of spectral decomposition in the first cell in x 
!! call DiffuiveBCsLeft1D3DLeg(TilFmi1Left,xmesh(0),xmesh(1),max_deg,ftil1_left,&
!!                             nodes_u,nodes_v,nodes_w,nodes_gwts,0.0_DP,T_w_left,gasR) 
!!! It seems like we have gotten the spectal coefficeints at the most left point.

!!! 
!end subroutine DiffuiveBCsLeftHOV1D

! DiffuiveBCsRight1D3DLeg (Fmi1Right,Fmi2Right,xl,xr,umesh,max_deg,f1,f2, &
!                                  u_gauss_nodes,u_gauss_weights,u_wleg_pol,u_w,T_w,R,Unodes, & 
!                                  Unodes_unit,Uweights_unit,moments_refine_u)
!
!  This procedure implements the diffusive boundary condtions on the LEFT boundary.
!  Namely, it takes the discrete solution in the most left spatial cell and 
!  calculates the flux of the paticles on the left boundary. To calculate the flux we find cells that 
!  only have incoming velocities and one (possible) cell that has both incoming and outgoing velocities.
!  The part of the flux that corresponds to the incoming only cells is computed using fast procedures, 
!  The cell that has both incoming and outgoing velocities is computed independently (however, there should 
!  not be such cell if non-uniform mesh in u of type =2 is used. This type of mesh will have velocity of the wall u_{w}= 0 
!  as a mesh point). Once the flux is computed, Maxwellian equilibrium distribution with given temperature T_{w} and 
!  velocity u_{w}=0 is used to compute the incoming fields. For the incoming fields, again, the "incoming only" cell 
!  are treated separately and the cell with the mixed velocities is special
! 
!  Attn: works with spectral coefficients, not characteristic variables
! 
! Expects: 
! 
! Fmi1Right, Fmi2Right ! the values of the spectral decomposition of the solution in the variable u
!                    ! 
!  xl,xr     ---  values of the most left interval in x: [xl, xr]
!  umesh     ---  meshpoints in u:
!  max_deg   ---  max_degree of the polynomial interpolation
!  
! 
! f1,f2  ---     ! f1(p,m,i) is the coefficients of the spectral representation in the most left cell 
!                                 ! -- coefficients of spectral decomposition in 
!                                 ! characteristic variables CORRESPONDING to the most left 
!                                 ! cell --- ftil(p,m,i,1):
!                                 ! -- p is the index in the basis functions in x 
!                                 ! -- m is the index in the basis functions in u
!                                 ! -- i is the cell in u
!                                 ! -- j = N is the cell in x
! u_gauss_nodes,u_wleg_pol --- values of the gauss nodes and the values of the legendre's polynomials 
! u_w  --- velocity of the wall (right now only u_w=0 is implemented)
! T_w  --- temperature of the wall
! R is the gas constant
! Unodes_unit,Uweights_unit,moments_refine_u --- these contain nodes for moment integration on [-1,1] after 
!                  ! a refinment of the interval (moments_refine_u) times, dublicated gauss weights corresponding to 
!                  ! a refined interval, and the coefficient for the refinement
!

subroutine DiffuiveBCsRight1D3DLeg (Fmi1Right,Fmi2Right,xl,xr,umesh,max_deg,f1,f2, &
                                   u_gauss_nodes,u_gauss_weights,u_wleg_pol,u_w,T_w,R,Unodes, & 
                                   Unodes_unit,Uweights_unit,moments_refine_u) 
use spectral_tools_mod
use distrib_mod 
! intrinsic MatMul
real (DP), intent (in) :: xl,xr ! points of the most left cell in x = [xmesh(N-1),xmesh(N)]
real (DP), dimension (0:), intent (in) :: umesh ! mesh points in variable u
integer (I4B), intent (in) :: max_deg ! this is the constant defining the maximum possible degree of approximating polynomials
real (DP), dimension (0:,0:,:), intent (in) :: f1,f2 ! f(p,m,i) is the main variable 
                                 ! -- coefficients of spectral decomposition in 
                                 ! characteristic variables CORRESPONDING to the most left 
                                 ! cell --- ftil(p,m,i,1):
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m is the index in the basis functions in u
                                 ! -- i is the cell in u
                                 ! -- j = N is the cell in x
real (DP), dimension (:), intent (in) :: u_gauss_nodes, u_gauss_weights ! nodes and weights of the gauss integration formula to be used in evaluation
real (DP), dimension (0:,:), intent (in) :: u_wleg_pol ! leg polynomials evalueated at guassian nodes (u_gauss_nodes)
                                                       ! u_wleg_pol(s,i)
                                                       ! s -- is the degree(number) of the basis function
                                                       ! i is the number fo the gaussian node (u_gauss_node)
real (DP), intent (in) :: u_w, T_w, R ! parameters of the diffusive BCs -- velocity and temperatture of the wall                                                         
real (DP), dimension (:,:), intent (in) :: Unodes ! pre-computed nodes on the refined cells in u (to be used
                                                  ! in the evaluation of moments)
                                                  ! Unodes(m,i) m -- node, i -- cell in u 
real (DP), dimension (:), intent (in) :: Unodes_unit ! to keep the refined nodes on the "canonical" cell [-1,1]
                                         ! Xnodes_unit(p) p -- nodes
                                         ! Unodes_unit(m) m -- nodes
real (DP), dimension (:), intent (in) :: Uweights_unit ! to keep the refined weights on the "canonical" cell [-1,1]
                                         ! Xweights_unit(p) p -- nodes
                                         ! Uweights_unit(m) m -- nodes
integer (I4B), intent (in) :: moments_refine_u ! the coefficient of the mesh refinement for refining the mesh in variable u
real (DP), dimension (0:,:), intent (out) :: Fmi1Right,Fmi2Right ! values of the f_{m,i}{t,x} 
											    ! at the right and left boundary points, FmiLeft(m,i)
                                                ! -- m is the index in the basis functions in u
                                                ! -- i is the cell in u
! supplementary variables

real (DP), dimension (size(u_gauss_nodes)) :: u ! nodes of the gauss integration formula to be used in evaluation of spectral coefficients
logical :: middle_cell_exists_flag ! This variable is =true if the wall velocity is not one of the mesh points, 
                                   ! but instead u_w falls in one of the cells in variable u
integer (I4B) :: k,s,M  ! k -- degree of polynomial approx in x
                                       ! s --- degree of polynomial approx in u, M -- number of cells in u
integer :: i,l,q,m_count,m_cell_i,loc_alloc_stat ! Some counters and integer variables                             
real (DP) :: n_w, h, du ! some dump variables to keep temp calculations; n_w -- will the the density of the 
                           ! the gas reflected from the wall
real (DP), dimension (:), allocatable :: sol! local variable where the solution is stored during the integration 
                                       ! weighted_sol(m) -- is value on one cell, m -- node number in u form uumesh(m) on U_{i} 
real (DP), dimension (:), allocatable :: Unodes_mcell ! -- will be used to evaluate integrals in u in the "median" cell in u                                      
real (DP), dimension (:,:), allocatable :: legpol_mcell ! legpol+mcell(size(x),0:k) k-- degree of the max polynomial. this will be needed to evaluate leg polynomials to deal with the middle cell
real (DP) :: mol_flux, beta1, beta2 ! temporary variable to keep the value of the molecular flux and the n(t,x)\sqrt{2\pi RT} factors                                    

!!!!!!!!!!!!
! will need these guys for brevity
k=size(f1,1)-1 
s=size(f1,2)-1
M=size(f1,3) 
! First we scan the mesh in u (umesh) to determine the boundary cell (if any) and the place where the 
! boundary cell is. 
! check if the velocity of the wall u_w is on the velocity mesh
if ((u_w < umesh(0)) .or. (u_w > umesh(M))) then 
  print *, "DiffuiveBCsLeft1D1DLeg: Velocity of the wall u_w is outsied the mesh in u"
		    stop
end if   
! now let us find the cell where u_w is: (It is not optimal to search of it on every time step, but don't know how to do it better
i=1
middle_cell_exists_flag = .FALSE.
do while ((i <= M-1) .and. (umesh(i) <= u_w))
  i=i+1
end do 
m_cell_i = i
du=(umesh(m_cell_i) - umesh(m_cell_i-1))/1024.0_DP
if (umesh(m_cell_i)-u_w < du ) then
  middle_cell_exists_flag = .FALSE.
  !! Let's add a warning if the "incident cells" are not present:
    if (m_cell_i >= M) then
    print *, "DiffuiveBCsRight1D1DLeg: Warning: no incident cells on the right boundary found in umesh."
    end if
  !! end warning
else 
  if ( u_w-umesh(m_cell_i-1) < du ) then 
     middle_cell_exists_flag = .FALSE.
     m_cell_i=m_cell_i-1
  else
     middle_cell_exists_flag = .TRUE.
  end if 
end if 
! Now that we know where the "incident cells in u" are (everything after the meshpoint umesh(m_cell_i)),
! if there is a "middle" cell, then also on this cell : [u_{w},umesh(m_cell_i)]
! it is time to calculate the flux: 
! first, let us prepare a variable to store intermediate result:
allocate (sol(1:size(Unodes_unit)), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
    print *, "DiffuiveBCsRight1D1DLeg: Allocation error for variable (sol)"
    end if
! 
mol_flux=0.0_DP
! first let us evaluate the portion of the flux that is due to the "strictly incident" cells.
do i=m_cell_i+1,M      !loop in strictrly "incident" cells in u
   du = (umesh(i) - umesh(i-1))/Real(moments_refine_u,DP) ! need this for the integration on subintervals
   ! first we evaluate the solution at nodes (xxmesh,uumesh)
   sol = EvCellLegStdrd1D1D_u_vector(1.0_DP,Unodes_unit,k,s,max_deg,f1(:,:,i))
   ! we then multiply the solution by gauss weights in u:
   sol = sol*Uweights_unit*Unodes(:,i)
   ! Now it is time to calculate the flux 
   mol_flux = mol_flux + sum(sol)*du
   ! end calculating flux in the cell [umesh(i-1), umesh(i)]
end do  ! end loop in (i) --- "incident" cells in u
! now we evaluate the portion of the flux for the "middle cell"
if (middle_cell_exists_flag .eqv. .TRUE.) then 
   du = ( umesh(m_cell_i) - u_w ) / Real(moments_refine_u, DP) 
   !! allocate space for nodes in u for integration in the "middle cell"
   allocate (Unodes_mcell(1:size(Unodes_unit)), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
    print *, "DiffuiveBCsLeft1D1DLeg: Allocation error for variable (Unodes_mcell)"
    end if
   !! end allocate 
   Unodes_mcell = (umesh(m_cell_i)+ u_w)/Real(2,DP) + Unodes_unit*du/Real(2,DP)
   sol = EvCellLeg1D1D_u_vector(xl,Unodes_mcell,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f1(:,:,m_cell_i))
   ! we then multiply the solution by gauss weights in u and velocity:
   sol = sol*Uweights_unit*Unodes_mcell
   ! Now it is time to calculate the flux 
   mol_flux = mol_flux + sum(sol)*du
   ! end calculating flux in the middle cell 
   deallocate (Unodes_mcell)
end if
mol_flux = mol_flux / Real(2,DP)  !!! comes from integration
deallocate (sol)
! We have calculated the incident flux. Next we need to evaluate the density of the reflected 
! Maxwellian distribution
n_w = mol_flux/( sqrt(R*T_w/pi25DT/Real(2,DP)) - u_w/Real(2,DP))
beta1=n_w/sqrt(Real(2,DP)*pi25DT*R*T_w)          ! we will need beta1 and beta2 to calculate distribution of 
beta2=n_w*sqrt(R*T_w/Real(2,DP)/pi25DT)          ! reflected particles
Fmi1Right=0
Fmi2Right=0
! first, let us prepare a variable to store intermediate result:
    allocate (sol(1:size(u_gauss_nodes)), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
    print *, "DiffuiveBCsRight1D1DLeg: Allocation error for variable (sol)"
    end if
    
! now we can set up the solution values on the left boundary
! We evaluate coefficients of the spectral decomposition, not in the characteristic variables,					  
! now we go through the "strictly reflecfted" cells
do i=1,m_cell_i-1
  !!! prepare integration with the cell
  u = (umesh(i) - umesh(i-1))*u_gauss_nodes/Real(2,DP) + (umesh(i) + umesh(i-1))/Real(2,DP)
  !!! Integration within the cell [umesh(i-1),umesh(i)]
  sol = ExponMaxwell1D1D(u, u_w, T_w, R)
  do l=0,s
      Fmi1Right(l,i) = sum(sol*u_wleg_pol(l,:))*(2.0_DP*Real(l,DP)+1.0_DP)/Real(2,DP) ! we only need one right now, the second one can be calculated from the first one (see below)
  end do 
! To get the correct distribution, need to multiply the maxwell's exponent by beta1 and beta2
Fmi2Right(:,i) = Fmi1Right(:,i)*beta2
Fmi1Right(:,i) = Fmi1Right(:,i)*beta1
end do
! now we need to do the middle cell, if it exists... 
if (middle_cell_exists_flag .eqv. .TRUE.) then 
  !! allocate space for nodes in values of the legendre's polynomials for integration in the "middle cell"
    allocate (legpol_mcell(size(u_gauss_nodes),0:s), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
    print *, "DiffuiveBCsRight1D1DLeg: Allocation error for variable (legpol_mcell)"
    end if
  !! end allocate 
  !!! First let us deal with the "incident" part of the cell... 
  !!! Integration within the cell [umesh(m_cell_i-1),umesh(m_cell_i)] -- only "the incident" portion
  !!! first, compute "incident" velocity coefficient contribution
    ! Component #1:
    u = (umesh(m_cell_i)-u_w)*u_gauss_nodes/Real(2,DP) + (umesh(m_cell_i)+u_w)/Real(2,DP)
    sol = EvCellLeg1D1D_u_vector(xl,u,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f1(:,:,m_cell_i))*u_gauss_weights
  !!! now we need nodes on[-1,1] that correspond to the values of the basis function 
    u = Real(2,DP)*( u - (umesh(m_cell_i) + umesh(m_cell_i-1))/Real(2,DP) )/(umesh(m_cell_i) - umesh(m_cell_i-1))
    legpol_mcell = EvalLegPol(s,u)
    do l=0,s
      Fmi1Right(l,m_cell_i) = sum(sol*legpol_mcell(:,l))*(2.0_DP*Real(l,DP)+1.0_DP)/Real(2,DP) &
                           *((umesh(m_cell_i)-u_w)/(umesh(m_cell_i) - umesh(m_cell_i-1)))
    end do 
    !  Component #2: 
    u = (umesh(m_cell_i)-u_w)*u_gauss_nodes/Real(2,DP) + (umesh(m_cell_i)+u_w)/Real(2,DP)
    sol = EvCellLeg1D1D_u_vector(xl,u,xl,xr,umesh(m_cell_i-1),umesh(m_cell_i),k,s,max_deg,f2(:,:,m_cell_i))*u_gauss_weights
  !!! now we need nodes on [-1,1] that correspond to the values of the basis function 
    do l=0,s
      Fmi2Right(l,m_cell_i) = sum(sol*legpol_mcell(:,l))*(2.0_DP*Real(l,DP)+1.0_DP)/Real(2,DP) &
                            *((umesh(m_cell_i)-u_w)/(umesh(m_cell_i) - umesh(m_cell_i-1)))
    end do 
  !!! Now let us do the "reflected part:
  !!! Components #1 and #2: 
  !!! prepare the "reflected nodes" 
    u = (u_w-umesh(m_cell_i-1))*u_gauss_nodes/Real(2,DP) + (u_w+umesh(m_cell_i-1))/Real(2,DP)
  !!! Integration within the cell [umesh(m_cell_i-1),umesh(m_cell_i)] -- only "the reflected" portion
    sol = ExponMaxwell1D1D(u, u_w, T_w, R)*u_gauss_weights 
  !!! now the nodes on [-1,1] for the basis function  
    u = 2*( u - (umesh(m_cell_i) + umesh(m_cell_i-1) )/Real(2,DP) )/(umesh(m_cell_i) - umesh(m_cell_i-1))
    legpol_mcell = EvalLegPol(s,u)
    do l=0,s
      du = sum(sol*legpol_mcell(:,l))*(2.0_DP*Real(l,DP)+1.0_DP)/Real(2,DP) & 
                   *( (u_w-umesh(m_cell_i-1))/(umesh(m_cell_i) - umesh(m_cell_i-1)) )
      Fmi1Right(l,m_cell_i) = Fmi1Right(l,m_cell_i) + du*beta1
      Fmi2Right(l,m_cell_i) = Fmi2Right(l,m_cell_i) + du*beta2
    end do 
  !!! end integration within the boundary cell 
    deallocate (legpol_mcell)
else 
  ! treat the cell right before m_cell_i as a reflected cell
  !!! prepare integration with the cell
    u = (umesh(m_cell_i) - umesh(m_cell_i-1))*u_gauss_nodes/Real(2,DP) + (umesh(m_cell_i) + umesh(m_cell_i-1))/Real(2,DP)
  !!! Integration within the cell [umesh(m_cell_1-1),umesh(m_cell_i)]
    sol = ExponMaxwell1D1D(u, u_w, T_w, R)
    do l=0,s
      Fmi1Right(l,m_cell_i) = sum(sol*u_wleg_pol(l,:))*(2.0_DP*Real(l,DP)+1.0_DP)/Real(2,DP) ! we only need one right now, the second one can be calculated from the first one (see below)
    end do 
  ! To get the correct distribution, need to multiply the maxwell's exponent by beta1 and beta2
    Fmi2Right(:,m_cell_i) = Fmi1Right(:,m_cell_i)*beta2
    Fmi1Right(:,m_cell_i) = Fmi1Right(:,m_cell_i)*beta1
end if  
deallocate (sol)
!!! Now let us calculate the boundary value at the strictly incident cell   
!!! Actually, we simply assemble the corresponding coefficients of the incident varaibles... 
do i=m_cell_i+1,M
   do l = 0,s
   do q = 0,min(max_deg - l,k)
   Fmi1Right(l,i) = Fmi1Right(l,i) + f1(q,l,i)
   Fmi2Right(l,i) = Fmi2Right(l,i) + f2(q,l,i)
   end do 
   end do 
end do    
!
end subroutine DiffuiveBCsRight1D3DLeg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DiffuiveBCsRightHOV1D (TilFmi1Right,TilFmi2Right,ftil1_right,ftil2_right) 
!
!  This procedure is a shell to call the function that implements the diffusive boundary condtions 
!  on the LEFT boundary. Specifically, this procedure converts characteristic variables into spectral coefficients. 
!  Then pulls together a lot parameters and call the subroutin that implements the diffusive conditions
!  Then the result is converted to characteristic variables
! 
!  This procedue looks up a lot of variables in the "common_variables_mod"
! real (DP), dimension (0:,0:,:) :: Ainvml, Aml 
!             ! and Aml(:,m,i)=A_{ml,i} is the array of eignevectors of of \hat{T}_{pl,i}=\int_{U_{i}}u\lambda_{l,i}(u)\lambda_{p,i}(u)
!             !          i -- is the interval in velocity U_{i}
!             !          m -- is the eigenvalue number
!             ! Ainvml(:,:,i) is the matrix inverse of Aml(:,:,i)
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DiffuiveBCsRightHOV1D (TilFmi1Right,TilFmi2Right,ftil1_right,ftil2_right)

use common_variables_mod, only: xmesh, umesh, Aml, Ainvml,k,s,M,N,xmesh,umesh,max_deg,&
                                   u_gauss_nodes, u_gauss_weights,u_wleg_pol,T_w_right,gasR,Unodes, & 
                                   Unodes_unit,Uweights_unit,moments_refine_u

real (DP), dimension (0:,:), intent (out) :: TilFmi1Right,TilFmi2Right ! values of the f_{m,i}{t,x} 
											    ! at the right and left boundary points, FmiLeft(m,i)
                                                ! -- m is the index in the basis functions in u
                                                ! -- i is the cell in u 
                                                
real (DP), dimension (0:,0:,:), intent (in) :: ftil1_right,ftil2_right ! f(p,m,i) -- auxiliary variable -- 
                                       ! coefficients of spectral decomposition 
                                       ! we will calculate them from the characteristic variables  
                                       !  f(p,m,i) <= ftil(p,m,i)
                                       ! -- p is the index in the basis functions in x 
                                       ! -- m is the index in the basis functions in u
                                       ! -- i is the cell in u
                                       ! -- these coefficients correspond to j=1 (first cell in x) in the ftil(p,m,i,j) 


                                                
! supplementary variables
real (DP), dimension (0:k,0:s,M) :: f1,f2 ! f(p,m,i) -- auxiliary variable -- 
                                       ! coefficients of spectral decomposition 
                                       ! we will calculate them from the characteristic variables  
                                       !  f(p,m,i) <= ftil(p,m,i)
                                       ! -- p is the index in the basis functions in x 
                                       ! -- m is the index in the basis functions in u
                                       ! -- i is the cell in u
                                       ! -- these coefficients correspond to j=1 (first cell in x) in the ftil(p,m,i,j) 
integer :: i,l,q,m_count ! some counters ... 

! First, we need to convert the characteristic coefficients into non-characteristic, that is
! we evaluatethe coefficients of the spectral representation from the characteristic representation
   f1=0
   f2=0
   do i = 1,M
   do l = 0,s
   do q = 0,min(max_deg - l,k)
     do m_count=0,s
     f1(q,l,i)= f1(q,l,i) + Aml(l,m_count,i)*ftil1_right(q,m_count,i)
     f2(q,l,i)= f2(q,l,i) + Aml(l,m_count,i)*ftil2_right(q,m_count,i)
     end do   
   end do
   end do 
   end do 
! now f(q,l,i) contains the coefficients of spectral decomposition in the first cell in x 
!!!call DiffuiveBCsRight1D1DLeg(TilFmi1Right,TilFmi2Right,xmesh(N-1),xmesh(N),umesh,max_deg,f1,f2, &
!!!                                   u_gauss_nodes,u_gauss_weights,u_wleg_pol,0.0_DP,T_w_right,gasR,Unodes, & 
!!!                                   Unodes_unit,Uweights_unit,moments_refine_u) 
!!! It seems like we have gotten the spectal coefficeints at the most left point.
!!! now they must be converted to the characteristic variables:
do i=1,size(umesh)-1
   TilFmi1Right(:,i)= MatMul(Real(Ainvml(:,:,i),DP),TilFmi1Right(:,i))
   TilFmi2Right(:,i)= MatMul(Real(Ainvml(:,:,i),DP),TilFmi2Right(:,i))
end do 
!!! 
end subroutine DiffuiveBCsRightHOV1D
!!!!! End added Alex 05/17/09     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module boundary_cond_mod 