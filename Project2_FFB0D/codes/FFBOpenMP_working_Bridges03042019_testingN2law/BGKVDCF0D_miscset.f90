!
!  miscset.f90
! 
!  Alex 9/13/2014
!
!  Miscellaneous setup for the BGKVDCF0D code. 
!  constains generation of meshes and more. 
!!!!!!!!!!!!!!!!!!!!!!!

module BGKVDCF0D_miscset

use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGblzmGnodes
! 
! This subroutine sets up the notes used in various integrations and 
! also used as ordimnates/nodal values 
!
! This subroutine uses parameters set by SetUWbgkParams
!
! This subroutine sets arrays used in other programs -- must be exectured before 
! subrouines setting meshes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetDGblzmGnodes
use BGKVDCF0D_commvar, only: x_gauss_n, x_gauss_w, t_gauss_n, t_gauss_w, &
             x_lobat_n, x_lobat_w, t_lobat_n, t_lobat_w, k_c, k_b, d_c, d_b, &
             g_nds_all, g_wts_all, g_max_order, &
             moments_x_gauss_nodes,moments_x_gauss_weights,moments_x_gauss_order
use gaussian_mod
use gaussquad, only: gauss_rule

integer (I4B) :: loc_alloc_stat, info  ! local variable to keep the allocation status
real (DP), dimension (:), allocatable :: ugn, ugw ! scratch arrays to create gauss nodes... 
real (DP), dimension (1:k_b) :: xb ! scratch array to use in the subroutine generating gauss-lobatto nodes. 
real (DP), dimension (1:d_b) :: tb ! scratch array to use in the subroutine to generate gauss-Lobatto nodes. 
integer (I4B) :: i ! local counter

!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gauss_n(1:max(k_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_gauss_n)"
     end if 
     !
  allocate (x_gauss_w(1:max(k_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_gauss_w)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (k_c == 1) then 
   x_gauss_n(:)= (/ 0.0_DP /)
   x_gauss_w(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),x_gauss_n,x_gauss_w)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
!!!
!!! Prepare the Gaussian nodes and weights in u variables !!!!!!!!!!!!!!!!!!!!!!
!   first we need to allocate the arrays to store gauss nodes and weights: 
!
  allocate (g_nds_all(1:g_max_order, 1:g_max_order), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (g_nds_all)"
     end if 
     !
  allocate (g_wts_all(1:g_max_order, 1:g_max_order), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (g_wts_all)"
     end if 
     !
!!! reset the arrays to zeros:
g_nds_all=0
g_wts_all=0
!!! 
g_nds_all(1,1)=0.0_DP
g_wts_all(1,1)=2.0_DP
do i=2,g_max_order 
 ! will need two scrap arrays... 
 allocate (ugn(1:i), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (ugn), i=", i
     end if 
 allocate (ugw(1:i), stat=loc_alloc_stat) 
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (ugw), i=", i
     end if 
call GauLeg(Real(-1,DP),Real(1,DP),ugn,ugw)
g_nds_all(1:i,i)=ugn
g_wts_all(1:i,i)=ugw
deallocate(ugn,ugw)
end do 
!!! Prepare the Gaussian nodes and weights in t variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (t_gauss_n(1:max(d_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_gauss_n)"
     end if 
     !
  allocate (t_gauss_w(1:max(d_c,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_gauss_w)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (d_c == 1) then 
   t_gauss_n(:)= (/ 0.0_DP /)
   t_gauss_w(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),t_gauss_n,t_gauss_w)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
!!!!!!!
! 
! Next we will set the Gauss-Lobatto nodes used in the nodal formulation. 
!!! Allocate arrays for gauss-lobatto nodes and weights
  allocate (t_lobat_n(1:max(d_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_lobat_n)"
     end if 
     !
  allocate (t_lobat_w(1:max(d_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (t_lobat_w)"
     end if 
   call gauss_rule(0, d_b, t_lobat_n, t_lobat_w, tb, 0.0_DP, 0.0_DP, 'B', info)  
!!! Allocate arrays for gauss-lobatto nodes and weights
  allocate (x_lobat_n(1:max(k_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_lobat_n)"
     end if 
     !
  allocate (x_lobat_w(1:max(k_b,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (x_lobat_w)"
     end if 

   call gauss_rule(0, k_b, x_lobat_n, x_lobat_w, xb, 0.0_DP, 0.0_DP, 'B', info)
   ! We set the variables of moments_x_gauss_nodes, moments_x_gauss_weights,
   allocate (moments_x_gauss_nodes(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_x_gauss_nodes)"
     end if 
     !
  allocate (moments_x_gauss_weights(1:moments_x_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_x_gauss_nodes)"
     end if 
 if (moments_x_gauss_order == 1) then 
   moments_x_gauss_nodes(:)= (/ 0.0_DP /)
   moments_x_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_x_gauss_nodes,moments_x_gauss_weights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
 end subroutine SetDGblzmGnodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGVblzmmesh 
!
! This subroutine sets up the meshes. Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters 
! N,M, mesh_x_uniform, moments_x_gauss_order, 
! are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles N,M.
!
! DO not call before arrays moments_x_gauss_nodes,moments_x_gauss_weights_x,moments_x_gauss_order
!          are selected! 
! 
! Most of the varaibles are looked up in the BGKDCF0D_commvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetBGK1Dmesh 
use BGKVDCF0D_commvar, only:  xmesh_l, xmesh_r, N, x_L, x_R, &
                   mesh_x_uniform, x_nonuniform_mesh_type, &
                   moments_x_gauss_nodes

intrinsic Min, Max, MinVal
                   
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer (I4B) :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) :: dx ! local temp variables to store 1/3 smallest distance between nodes.
 ! first the mesh in x
 if (mesh_x_uniform) then 
  allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
  stop
  end if
     !
  xmesh_l = (/ (x_L + (x_R - x_L )/Real(N,DP)*i, i=0,N-1) /)
  !
  allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
  stop
  end if
     !
  xmesh_r = (/ (x_L + (x_R - x_L )/Real(N,DP)*i, i=1,N) /)
  !
 else
  select case (x_nonuniform_mesh_type)
  case (1) ! use "the gauss nodes for evaluating moments in x" to set up the "mesh in x"
           ! interval [x_left,x_right] is divided in subintervals as to have gauss nodes at centerpoints 
           ! (and some extra points need to be introduced to make this possible). 
  mxgo = size(moments_x_gauss_nodes) ! temp remember order of the guass method for integration of moments .... 
  dx=min(minval(moments_x_gauss_nodes(2:mxgo) - moments_x_gauss_nodes(1:mxgo-1)),&
                 moments_x_gauss_nodes(1)+Real(1,DP), Real(1,DP)-moments_x_gauss_nodes(mxgo))/Real(3,DP)
  !! allocate memory for the nodes in x
     allocate (xmesh_l(1:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
     stop
     end if
     !
!! allocate memory for the nodes in x
     allocate (xmesh_r(1:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
     stop
     end if
     !
  !! Now we set up the mesh in x:
  xmesh_l(1)=x_L
  do i=1,mxgo
  xmesh_l(2*i) = (x_R+x_L)/Real(2,DP) + (x_R-x_L)/Real(2,DP)*(moments_x_gauss_nodes(i) - dx)
  xmesh_l(2*i+1) = (x_R+x_L)/Real(2,DP) + (x_R-x_L)/Real(2,DP)*(moments_x_gauss_nodes(i) + dx)
  end do 
  xmesh_r(1:2*mxgo)=xmesh_l(2:2*mxgo+1)
  xmesh_r(mxgo*2+1) = x_R
  N=mxgo*2+1 ! the new value of (N) in case somebody uses it directly... 
  !! mesh in x is ready
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/8-1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=6 !!! 
  !! allocate memory for the nodes in x
  if (N<6) then 
    print *, "SetDGVblzmmesh: For the type 2 nonuniform mesh in x the number of mesh points has to be >=6 (xmesh_l)"
    stop  
  end if  
    allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
    !
    allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
    !
    dx=(x_R-x_L)/(Real(N-4,DP)-Real(1,DP)/Real(4,DP))
  !! Now we will set up the mesh in variable x:
    xmesh_l(1)=x_L
    xmesh_l(2)=x_L+dx/Real(8,DP)
    xmesh_l(3)=xmesh_l(2)+dx/Real(4,DP)
    xmesh_l(4)=xmesh_l(3)+dx/Real(2,DP)
     do i=5,N-3
      xmesh_l(i)=xmesh_l(i-1)+dx
     end do 
    xmesh_r(N)=x_R
    xmesh_l(N)=xmesh_r(N)-dx/Real(8,DP)
    xmesh_l(N-1)=xmesh_l(N)-dx/Real(4,DP)
    xmesh_l(N-2)=xmesh_l(N-1)-dx/Real(2,DP)
    xmesh_r(1:N-1)=xmesh_l(2:N)         
  !! end of the non-uniform mesh type 2. 
  case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=4 !!! 
  !! allocate memory for the nodes in x
  if (N<4) then 
    print *, "SetDGVblzmmesh: For the type 3 nonuniform mesh in x the number of mesh points has to be >=4 (xmesh_l)"
    stop  
  end if  
    allocate (xmesh_l(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
    !
    allocate (xmesh_r(1:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
    !
    dx=(x_R-x_L)/(Real(N-2,DP)-Real(1,DP)/Real(2,DP))
  !! Now we will set up the mesh in variable x:
    xmesh_l(1)=x_L
    xmesh_l(2)=x_L+dx/Real(4,DP)
    xmesh_l(3)=xmesh_l(2)+dx/Real(2,DP)
     do i=4,N-2
      xmesh_l(i)=xmesh_l(i-1)+dx
     end do 
    xmesh_r(N)=x_R
    xmesh_l(N)=xmesh_r(N)-dx/Real(4,DP)
    xmesh_l(N-1)=xmesh_l(N)-dx/Real(2,DP)
    xmesh_r(1:N-1)=xmesh_l(2:N)
  !!! end of the non-uniform mesh type 3. 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  case default 
    allocate (xmesh_l(1:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_l)"
    stop
    end if
     !
    allocate (xmesh_r(1:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "SetDGVblzmmesh: Allocation error for variable (xmesh_r)"
    stop
    end if
     ! 
    xmesh_l = (/ (x_L + (x_R-x_L)/Real(3,DP)*i, i=0,2) /)
    xmesh_r(1:2)=xmesh_l(2:3)
    xmesh_r(3) = x_R 
    N=3 ! the new value of (N) in case somebody uses it directly... 
    print *, "SetDGVblzmmesh: unsupported choice of non-uniform mesh. x_nonuniform_mesh_type= ", x_nonuniform_mesh_type
end select
 end if  
 
 end subroutine SetBGK1Dmesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  Subroutine OldCreateUmesh 
!  This subroutine is used to generate 1D meshes in the velocity variables 
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine OldCreateUmesh(ugrid,ugrid_cap,mesh_u_uniform,u_nonuniform_mesh_type,u_L,u_R) 
!
real (DP), dimension (:), intent (out)  :: ugrid
integer (I4B), dimension (:), intent (out)  :: ugrid_cap
logical, intent(in) :: mesh_u_uniform
integer (I4B), intent (in) :: u_nonuniform_mesh_type
real (DP), intent (out) :: u_L, u_R
! local variables: 
integer (I4B) :: M ! number of cells in the grid
integer (I4B) :: i,j ! local counters
!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: smr, sml ! refinement parameters for mesh of type 3
real (DP) :: du 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
M=size(ugrid,1)-1
ugrid_cap(1) = M+1
! 
if (mesh_u_uniform) then 
  ugrid = (/ (u_L + (u_R-u_L)/Real(M,DP)*i, i=0,M) /) 
  else
  select case (u_nonuniform_mesh_type)
  case (1) ! Unsupported  
     print *, "OldCreateUmesh: Unsupported type of mesh in velocity variable (u_nonuniform_mesh_type=1)"
     stop
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! Only Static wall is implemented at this time! 
           ! 
 ! first, let us check that u_wall =0  is in the velocity interval                     
  if ((u_L > 0.0_DP) .or. (u_R < 0.0_DP)) then 
    u_L = min (0.0_DP,u_L)
    u_R = max (0.0_DP,u_R)
  end if                    
  ! now we start to fill in the mesh:
  i=-M
  j=1
  ugrid(1)=u_L
  do while ( u_R - (u_R-u_L)/Real(M,DP)*i > (u_R-u_L)/Real(M,DP)/2 - 0.00000000001_DP )    
  if (i*(u_R-u_L)/Real(M,DP)-u_L > (u_R-u_L)/Real(M,DP)/2) then 
  j=j+1
  ugrid(j)=i*(u_R-u_L)/Real(M,DP)
  end if 
  i=i+1
  end do 
  ugrid(M+1)=u_R
  ! now we have filled in the mesh:
 case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! also it will surround $u=0$ with small cells. The sizes of the sell are obtained 
           ! by refining the size of biggest cells $smr^{sml}$ times. Here $smr$ -- is the 
           ! refinement factor and $sml$ is the total level of reinements. 
           ! Only Static wall is implemented at this time! 
           ! 
 ! First we set values for $smr$ and $sml$ 
 smr=10
 sml=2
 if (2*sml > M-2) then
  print *, "SetDGVblzmmesh: selected number of refinements 'sml' of mesh is too large for the given M"
  stop
 end if
   ! now we need to evaluate the size of the largest cell. 
 du=(u_R-u_L)/((M-2*sml) + 2.0_DP/Real(smr,DP)* &
             (1.0_DP/(Real(smr,DP)**sml)-1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP))
 ! now we start to fill in the mesh:
  j=1
  ugrid(1) = u_L
  do while (ugrid(j) < - du*(1.0_DP/(Real(smr,DP)**sml) - 1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP)/Real(smr,DP) - 0.00000000001_DP )    
  j=j+1
  ugrid(j)=ugrid(j-1)+du
  end do
  ! Now we fill in the small cells
  do i=1,sml
  j=j+1
  ugrid(j) = ugrid(j-1) + du/(Real(smr,DP)**i)
  end do 
  ugrid(j)=0.0_DP ! Just to have it exact --- otherwise it is not zero, but 10^-13
  do i=sml,1,-1
  j=j+1
  ugrid(j) = ugrid(j-1) + du/(Real(smr,DP)**i)
  end do 
  ! now the rest of the cells.
  do while (ugrid(j) < u_R - du - 0.00000000001_DP)
  j = j+1
  ugrid(j) = ugrid(j-1)+du
  end do
  ugrid(M+1) = u_R
 case default 
     print *, "SetDGVblzmmesh: Unknown type of mesh"
     stop
  end select
 end if  
end subroutine OldCreateUmesh


end module BGKVDCF0D_miscset
