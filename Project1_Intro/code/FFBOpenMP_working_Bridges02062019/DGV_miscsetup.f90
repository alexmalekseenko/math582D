!
!  miscset.f90
! 
!  Alex 9/13/2014
!
!  Miscellaneous setup for the Discontinuous Galerkin Velocity Library code. 
!  constains subroutines for generating grids, meshes, etc. in the velocity space 
!  also subroutines for set up various complementary arrays and variables. 
!!!!!!!!!!!!!!!!!!!!!!!

module DGV_miscset

use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGblzmGnodes
! 
! This subroutine sets up the notes used in various integrations and 
! also used as ordinates/nodal values 
!
! This subroutine uses parameters set by SetDGVParams
!
! This subroutine sets arrays used in other programs -- must be run before 
! SetDGVblzmmesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetGDVGnodes
use DGV_commvar, only: g_nds_all, g_wts_all, g_max_order, moments_u_gauss_nodes, moments_u_gauss_weights, &
             moments_u_gauss_order
use gaussian_mod
use gaussquad, only: gauss_rule

integer (I4B) :: loc_alloc_stat, info  ! local variable to keep the allocation status
real (DP), dimension (:), allocatable :: ugn, ugw ! scratch arrays to create gauss nodes... 
integer (I4B) :: i ! local counter

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
   ! We set the variables of moments_u_gauss_nodes, moments_u_gauss_weights
    
 allocate (moments_u_gauss_nodes(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_u_gauss_nodes)"
     end if 
     !
  allocate (moments_u_gauss_weights(1:moments_u_gauss_order), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetDGblzmGnodes: Allocation error for variable (moments_u_gauss_nodes)"
     end if 
 if (moments_u_gauss_order == 1) then 
   moments_u_gauss_nodes(:)= (/ 0.0_DP /)
   moments_u_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),moments_u_gauss_nodes,moments_u_gauss_weights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 

end subroutine SetGDVGnodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGVblzmmesh 
!
! This subroutine sets up the meshes. Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters 
! su,sv,sw,Mu,Mv,Mw, mesh_u_uniform, 
! moments_u_gauss_order are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles su,Mu,svMv,sw,Mw.
!
! DO not call before arrays 
!          moments_u_gauss_order,moments_u_gauss_nodes,moments_u_gauss_weights are selected! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetDGVblzmmesh 
use DGV_commvar, only:  Mu, Mv, Mw, u_L, u_R, &
                   v_L, v_R, w_L, w_R, grids_cap_u, grids_cap_v, grids_cap_w, grids_u, grids_v, grids_w, & 
                   mesh_u_uniform, u_nonuniform_mesh_type,&
                   moments_u_gauss_nodes, mesh_v_uniform, mesh_w_uniform, & 
                   v_nonuniform_mesh_type,w_nonuniform_mesh_type

intrinsic Min, Max, MinVal
                   
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer (I4B) :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) ::  du ! local temp variables to store 1/3 smallest distance between nodes.
 
 !  we will setup the level zero meshes in u ( grids_u/_v/_w and grids_cap_u/_v/_w )  
 ! first the mesh in u
  allocate (grids_u(1:Mu+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_u)"
  stop
  end if
     !
  allocate (grids_cap_u(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_u)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_u,grids_cap_u,mesh_u_uniform,u_nonuniform_mesh_type,u_L,u_R)
  ! then mesh in v 
  allocate (grids_v(1:Mv+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_v)"
  stop
  end if
     !
  allocate (grids_cap_v(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_v)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_v,grids_cap_v,mesh_v_uniform, v_nonuniform_mesh_type,v_L,v_R)
  ! now the mesh in w
  allocate (grids_w(1:Mw+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_w)"
  stop
  end if
     !
  allocate (grids_cap_w(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmesh: Allocation error for variable (grids_cap_w)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_w,grids_cap_w,mesh_w_uniform,w_nonuniform_mesh_type,w_L,w_R)
     !
 end subroutine SetDGVblzmmesh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGVblzmmeshII
!
! This subroutine sets up the secondary meshes. Do not detouch from the main program 
! 
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
! 
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! This subroutine is highly dependent on the common variables and preparatory subrouines. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters 
! suII,svII,swII,MuII,MvII,MwII, mesh_u_uniform, 
! moments_u_gauss_order are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles suII,MuII,svII, MvII,swII,MwII.
!
! DO not call before arrays 
!          moments_u_gauss_order,moments_u_gauss_nodes,moments_u_gauss_weights are selected! 
! 
! Most of the varaibles are looked up in the common varaibles module (DGV_commvar)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetDGVblzmmeshII 
use DGV_commvar, only:  Mu=>MuII, Mv=>MvII, Mw=>MwII, &
                   u_L, u_R, v_L, v_R, w_L, w_R, &
                   grids_cap_u=>grids_cap_uII, grids_cap_v=>grids_cap_vII, grids_cap_w=>grids_cap_wII,  &
                   grids_u=>grids_uII, grids_v=>grids_vII, grids_w=>grids_wII, & 
                   mesh_u_uniform, u_nonuniform_mesh_type,&
                   moments_u_gauss_nodes, mesh_v_uniform, mesh_w_uniform, & 
                   v_nonuniform_mesh_type,w_nonuniform_mesh_type

intrinsic Min, Max, MinVal
                   
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer (I4B) :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) ::  du ! local temp variables to store 1/3 smallest distance between nodes.
 
 !  we will setup the level zero meshes in u ( grids_u/_v/_w and grids_cap_u/_v/_w )  
 ! first the mesh in u
  allocate (grids_u(1:Mu+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_uII)"
  stop
  end if
     !
  allocate (grids_cap_u(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_cap_uII)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_u,grids_cap_u,mesh_u_uniform,u_nonuniform_mesh_type,u_L,u_R)
  ! then mesh in v 
  allocate (grids_v(1:Mv+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_vII)"
  stop
  end if
     !
  allocate (grids_cap_v(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_cap_vII)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_v,grids_cap_v,mesh_v_uniform, v_nonuniform_mesh_type,v_L,v_R)
  ! now the mesh in w
  allocate (grids_w(1:Mw+1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_wII)"
  stop
  end if
     !
  allocate (grids_cap_w(1:1), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "SetDGVblzmmeshII: Allocation error for variable (grids_cap_wII)"
  stop
  end if
     !   
  call OldCreateUmesh(grids_w,grids_cap_w,mesh_w_uniform,w_nonuniform_mesh_type,w_L,w_R)
     !
 end subroutine SetDGVblzmmeshII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  Subroutine OldCreateUmesh 
!  This subroutine is used to generate 1D meshes in the velocity variables 
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set3DCellsR_DGV
!
! This subroutine creates a 3D velocity mesh. 
!
! The subroutine looks up the grids arrays and creates 3D cells from the grids arrays. 
! If more than one grid is defined, grids are treated as non-related. Their ierarchy is ignored -- actually it is 
! not known from the grids array rigth away and needs to be discovered by a painful process. So we just do not do it. 
! instead, we would like to use this subroutine to create zero level mesh and we will hope that we only have one grid at this point.
! 
! Later, another subroutine may be called to create embedded grids. 
!
! Depends on the main program (other subroutines must be exectured, first) 
! and the DGV_commvar. Before calling make sure su,sv,sw have been initialized!
! Also make sure that 1D grids in velocity have been allocated and initialized. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Set3DCellsR_DGV

use DGV_commvar, only: grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,su,sv,sw
!
integer (I4B) :: iu,iv,iw, ig ! some local counters

integer (I4B) :: mm,mx,igu,igv,igw ! local scrap counters,  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creating cells from whaterwer is already in the grids -- refinement will be applied later 
!

! we assume that at least one grid of level zero is defined...
if (size(grids_u,1) < 1) then 
  print *, "Set3DCellsR_DGV: Error, there seem to be no 1D grids_u/_v/_w defined."
  stop 
end if 
mm=(grids_cap_u(1)-1)*(grids_cap_v(1)-1)*(grids_cap_w(1)-1)
if (mm < 1) then 
  print *, "Set3DCellsR_DGV: Error, at least one of the 1D grids_u/_v/_w defined."
  stop 
end if 
! create the cells arrays... 
call AllocCellsArrsDGV(mm)
! now we will need to populate the cells: 
call FillCellsArrsDGV(1,mm, grids_u(1:grids_cap_u(1)),grids_v(1:grids_cap_v(1)),grids_w(1:grids_cap_w(1)),1,su,sv,sw)
! set some counters in case we have more stuff to process...
mx = mm  ! this will accumulate the total number of the recorded cells ...
igu=grids_cap_u(1)
igv=grids_cap_v(1)
igw=grids_cap_w(1)
! now if there is more grids... we needs to continue to create cells... 
do ig=2,size(grids_cap_u,1)   ! all grids_u/_v/_w arrays are of the same length 
 ! for each grid create cells...  
 mm=(grids_cap_u(ig)-1)*(grids_cap_v(ig)-1)*(grids_cap_w(ig)-1) !! this is how many new cells will be on this grid
 ! extend the cells arrays to fit more cells ...
 call ExtendCellsArrsDGV(mx+mm)
 ! now we will need to populate the cells: 
 call FillCellsArrsDGV(mx+1,mx+mm,grids_u(igu+1:igu+grids_cap_u(ig)),grids_v(igv+1:igv+grids_cap_v(ig)), & 
             grids_w(igw+1:igw+grids_cap_w(ig)),ig,su,sv,sw)
 ! set some counters in case we have more stuff to process...
 mx = mx+mm  ! this will accumulate the total number of the recorded cells ...
 igu=igu+grids_cap_u(ig)
 igv=igv+grids_cap_v(ig)
 igw=igw+grids_cap_w(ig)
 !
end do 

end subroutine Set3DCellsR_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AllocCellsArrsDGV
!! This is a "macro" subrouting to save some space in another subroutine. 
!! given the sizes, it allocates the arrays and returns emply arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AllocCellsArrsDGV(M)
use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent (in) :: M ! The size of the arrays...  
!
integer (I4B) :: loc_alloc_stat
!!!!!!!!!!!!!!
!
allocate (cells_pgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_pgrid)"
  stop
  end if
allocate (cells_cgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_cgrid)"
  stop
  end if
allocate (cells_lu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lu)"
  stop
  end if
allocate (cells_lv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lv)"
  stop
  end if
allocate (cells_lw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_lw)"
  stop
  end if
allocate (cells_ru(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_ru)"
  stop
  end if
allocate (cells_rv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_rv)"
  stop
  end if
allocate (cells_rw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_rw)"
  stop
  end if
allocate (cells_refu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refu)"
  stop
  end if
allocate (cells_refv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refv)"
  stop
  end if
allocate (cells_refw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_refw)"
  stop
  end if
allocate (cells_gou(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gou)"
  stop
  end if
allocate (cells_gov(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gov)"
  stop
  end if
allocate (cells_gow(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGV: Allocation error for variable (cells_gow)"
  stop
  end if
end subroutine AllocCellsArrsDGV

subroutine ExtendCellsArrsDGV(M)
use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent (in) :: M ! The new size of the cells arrays...  
!
integer (I4B) :: loc_alloc_stat
integer (I4B) :: Md ! a crap variable to keed the old size 
!!!!!!!!!!!!!!
real (DP), dimension(:), allocatable ::cpg,ccg,clu,clv,clw,cru,crv,crw,&
                                       crefu,crefv,crefw,cgou,cgov,cgow  

Md=size(cells_pgrid, 1)
if (Md>M) then 
 print *, "ExtendCellsArrsDGV: new size for arrays cells is smaller then the old one... "
 stop
end if
!!!!
!!!! 
allocate (cpg(1:Md),ccg(1:Md),clu(1:Md),clv(1:Md),clw(1:Md),cru(1:Md),crv(1:Md),crw(1:Md), stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variables cpg,ccg,clu,clv,clw,cru,crv,crw"
  stop
  end if
allocate (crefu(1:Md),crefv(1:Md),crefw(1:Md),cgou(1:Md),cgov(1:Md),cgow(1:Md),stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variables cells_refw, cells_gow, cells_gou, cells_gov"
  stop
  end if
! save the cells arrays in the temp arrays... 
cpg=cells_pgrid
ccg=cells_cgrid
clu=cells_lu 
clv=cells_lv 
clw=cells_lw
cru=cells_ru
crv=cells_rv
crw=cells_rw
crefu=cells_refu
crefv=cells_refv
crefw=cells_refw 
cgow=cells_gow
cgou=cells_gou 
cgov=cells_gov
! deallocate the old cells arrays... 
deallocate (cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov)
! now we allocate with the new size M   
allocate (cells_pgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_pgrid)"
  stop
  end if
allocate (cells_cgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_cgrid)"
  stop
  end if
allocate (cells_lu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lu)"
  stop
  end if
allocate (cells_lv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lv)"
  stop
  end if
allocate (cells_lw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_lw)"
  stop
  end if
allocate (cells_ru(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_ru)"
  stop
  end if
allocate (cells_rv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_rv)"
  stop
  end if
allocate (cells_rw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_rw)"
  stop
  end if
allocate (cells_refu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_refu)"
  stop
  end if
allocate (cells_refv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV Allocation error for variable (cells_refv)"
  stop
  end if
allocate (cells_refw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_refw)"
  stop
  end if
allocate (cells_gou(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gou)"
  stop
  end if
allocate (cells_gov(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gov)"
  stop
  end if
allocate (cells_gow(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGV: Allocation error for variable (cells_gow)"
  stop
  end if
!! and the inverse assignment:
cells_pgrid(1:Md)=cpg
cells_cgrid(1:Md)=ccg
cells_lu(1:Md)=clu
cells_lv(1:Md)=clv
cells_lw(1:Md)=clw
cells_ru(1:Md)=cru
cells_rv(1:Md)=crv
cells_rw(1:Md)=crw
cells_refu(1:Md)=crefu
cells_refv(1:Md)=crefv
cells_refw(1:Md)=crefw
cells_gow(1:Md)=cgow
cells_gou(1:Md)=cgou
cells_gov(1:Md)=cgov
!! the rest will be filled outside. 
deallocate(cpg,ccg,clu,clv,clw,cru,crv,crw,crefu,crefv,crefw,cgou,cgov,cgow)
  
end subroutine ExtendCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Set3DCellsR_DGVII
!
! This subroutine creates a 3D velocity mesh. 
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! The subroutine looks up the grids arrays and creates 3D cells from the grids arrays. 
! If more than one grid is defined, grids are treated as non-related. Their ierarchy is ignored -- actually it is 
! not known from the grids array rigth away and needs to be discovered by a painful process. So we just do not do it. 
! instead, we would like to use this subroutine to create zero level mesh and we will hope that we only have one grid at this point.
! 
! Later, another subroutine may be called to create embedded grids. 
!
! Depends on the main program (other subroutines must be exectured, first) 
! and the D4GV_commvar. Before calling make sure suII,svII,swII have been initialized!
! Also make sure that secondary 1D grids in velocity have been allocated and initialized. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Set3DCellsR_DGVII

use DGV_commvar, only: grids_cap_u=>grids_cap_uII, grids_cap_v=>grids_cap_vII, &
                       grids_cap_w=>grids_cap_wII, & 
                       grids_u=>grids_uII, grids_v=>grids_vII, grids_w=>grids_wII, &
                       su=>suII, sv=>svII, sw=>swII
!
integer (I4B) :: iu,iv,iw, ig ! some local counters

integer (I4B) :: mm,mx,igu,igv,igw ! local scrap counters,  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creating cells from whaterwer is already in the grids -- refinement will be applied later 
!

! we assume that at least one grid of level zero is defined...
if (size(grids_u,1) < 1) then 
  print *, "Set3DCellsR_DGVII: Error, there seem to be no 1D grids_u/_v/_w defined."
  stop 
end if 
mm=(grids_cap_u(1)-1)*(grids_cap_v(1)-1)*(grids_cap_w(1)-1)
if (mm < 1) then 
  print *, "Set3DCellsR_DGVII: Error, at least one of the 1D grids_u/_v/_w defined."
  stop 
end if 
! create the cells arrays... 
call AllocCellsArrsDGVII(mm)
! now we will need to populate the cells: 
call FillCellsArrsDGVII(1,mm, grids_u(1:grids_cap_u(1)),grids_v(1:grids_cap_v(1)),grids_w(1:grids_cap_w(1)),1,su,sv,sw)
! set some counters in case we have more stuff to process...
mx = mm  ! this will accumulate the total number of the recorded cells ...
igu=grids_cap_u(1)
igv=grids_cap_v(1)
igw=grids_cap_w(1)
! now if there is more grids... we needs to continue to create cells... 
do ig=2,size(grids_cap_u,1)   ! all grids_u/_v/_w arrays are of the same length 
 ! for each grid create cells...  
 mm=(grids_cap_u(ig)-1)*(grids_cap_v(ig)-1)*(grids_cap_w(ig)-1) !! this is how many new cells will be on this grid
 ! extend the cells arrays to fit more cells ...
 call ExtendCellsArrsDGVII(mx+mm)
 ! now we will need to populate the cells: 
 call FillCellsArrsDGVII(mx+1,mx+mm,grids_u(igu+1:igu+grids_cap_u(ig)),grids_v(igv+1:igv+grids_cap_v(ig)), & 
             grids_w(igw+1:igw+grids_cap_w(ig)),ig,su,sv,sw)
 ! set some counters in case we have more stuff to process...
 mx = mx+mm  ! this will accumulate the total number of the recorded cells ...
 igu=igu+grids_cap_u(ig)
 igv=igv+grids_cap_v(ig)
 igw=igw+grids_cap_w(ig)
 !
end do 

end subroutine Set3DCellsR_DGVII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! AllocCellsArrsDGVII
!!
!!  This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
!!  This is accomplished by re-naming variable in the USE commvar, only: statement. 
!!
!!
!! The suffix II in the name suggests that the subroutine works with the secondary mesh
!!
!! This is a "macro" subrouting to save some space in another subroutine. 
!! given the sizes, it allocates the arrays and returns emply arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AllocCellsArrsDGVII(M)

use DGV_commvar, only: cells_pgridII, cells_cgridII, &
                   cells_luII, cells_lvII, cells_lwII, & 
                   cells_ruII, cells_rvII, cells_rwII, &
                   cells_refuII, cells_refvII, cells_refwII, &
                   cells_gowII, cells_gouII, cells_govII

integer (I4B), intent (in) :: M ! The size of the arrays...  
!
integer (I4B) :: loc_alloc_stat
!!!!!!!!!!!!!!
!
allocate (cells_pgridII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_pgridIIII)"
  stop
  end if
allocate (cells_cgridII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_cgridIIII)"
  stop
  end if
allocate (cells_luII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_luIIII)"
  stop
  end if
allocate (cells_lvII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_lvIIII)"
  stop
  end if
allocate (cells_lwII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_lwIIII)"
  stop
  end if
allocate (cells_ruII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_ruIIII)"
  stop
  end if
allocate (cells_rvII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_rvIIII)"
  stop
  end if
allocate (cells_rwII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_rwIIII)"
  stop
  end if
allocate (cells_refuII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_refuIIII)"
  stop
  end if
allocate (cells_refvII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_refvIIII)"
  stop
  end if
allocate (cells_refwII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_refwIIII)"
  stop
  end if
allocate (cells_gouII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_gouIIII)"
  stop
  end if
allocate (cells_govII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_govIIII)"
  stop
  end if
allocate (cells_gowII(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "AllocCellsArrsDGVII: Allocation error for variable (cells_gowIIII)"
  stop
  end if
end subroutine AllocCellsArrsDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ExtendCellsArrsDGVII
!! This is a "macro" subroutine to save some space in another subroutine. 
!! given the sizes, this subroutine extends the cells arrays
!!
!! This is a copy of an above subroutine, with the only change that it works with secondary mehses. 
!! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!!
!!
!! The suffix II in the name suggests that the subroutine works with the secondary mesh
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendCellsArrsDGVII(M)

use DGV_commvar, only: cells_pgrid => cells_pgridII, cells_cgrid => cells_cgridII,  &
    cells_lu => cells_luII, cells_lv => cells_lvII, cells_lw => cells_lwII, &
    cells_ru => cells_ruII, cells_rv => cells_rvII, cells_rw => cells_rwII, &
    cells_refu => cells_refuII, cells_refv => cells_refvII, cells_refw => cells_refwII, &
                   cells_gou=>cells_gouII, cells_gov=>cells_govII, cells_gow=>cells_gowII 

integer (I4B), intent (in) :: M ! The new size of the cells arrays...  
!
integer (I4B) :: loc_alloc_stat
integer (I4B) :: Md ! a crap variable to keed the old size 
!!!!!!!!!!!!!!
real (DP), dimension(:), allocatable ::cpg,ccg,clu,clv,clw,cru,crv,crw,&
                                       crefu,crefv,crefw,cgou,cgov,cgow  

Md=size(cells_pgrid, 1)
if (Md>M) then 
 print *, "ExtendCellsArrsDGVII: new size for arrays cells is smaller then the old one... "
 stop
end if
!!!!
!!!! 
allocate (cpg(1:Md),ccg(1:Md),clu(1:Md),clv(1:Md),clw(1:Md),cru(1:Md),crv(1:Md),crw(1:Md), stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variables cpg,ccg,clu,clv,clw,cru,crv,crw"
  stop
  end if
allocate (crefu(1:Md),crefv(1:Md),crefw(1:Md),cgou(1:Md),cgov(1:Md),cgow(1:Md),stat=loc_alloc_stat)
!
if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variables crefu, crefv, crefw, cgou, cgov, cgow"
  stop
  end if
! save the cells arrays in the temp arrays... 
cpg=cells_pgrid
ccg=cells_cgrid
clu=cells_lu 
clv=cells_lv 
clw=cells_lw
cru=cells_ru
crv=cells_rv
crw=cells_rw
crefu=cells_refu
crefv=cells_refv
crefw=cells_refw 
cgow=cells_gow
cgou=cells_gou 
cgov=cells_gov
! deallocate the old cells arrays... 
deallocate (cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov)
! now we allocate with the new size M   
allocate (cells_pgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_pgridII)"
  stop
  end if
allocate (cells_cgrid(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_cgridII)"
  stop
  end if
allocate (cells_lu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_luII)"
  stop
  end if
allocate (cells_lv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_lvII)"
  stop
  end if
allocate (cells_lw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_lwII)"
  stop
  end if
allocate (cells_ru(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_ruII)"
  stop
  end if
allocate (cells_rv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_rvII)"
  stop
  end if
allocate (cells_rw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_rwII)"
  stop
  end if
allocate (cells_refu(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_refuII)"
  stop
  end if
allocate (cells_refv(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_refvII)"
  stop
  end if
allocate (cells_refw(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_refwII)"
  stop
  end if
allocate (cells_gou(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_gouII)"
  stop
  end if
allocate (cells_gov(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_govII)"
  stop
  end if
allocate (cells_gow(1:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "ExtendCellsArrsDGVII: Allocation error for variable (cells_gowII)"
  stop
  end if
!! and the inverse assignment:
cells_pgrid(1:Md)=cpg
cells_cgrid(1:Md)=ccg
cells_lu(1:Md)=clu
cells_lv(1:Md)=clv
cells_lw(1:Md)=clw
cells_ru(1:Md)=cru
cells_rv(1:Md)=crv
cells_rw(1:Md)=crw
cells_refu(1:Md)=crefu
cells_refv(1:Md)=crefv
cells_refw(1:Md)=crefw
cells_gow(1:Md)=cgow
cells_gou(1:Md)=cgou
cells_gov(1:Md)=cgov
!! the rest will be filled outside. 
deallocate(cpg,ccg,clu,clv,clw,cru,crv,crw,crefu,crefv,crefw,cgou,cgov,cgow)
  
end subroutine ExtendCellsArrsDGVII 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine deallocates the cells arrays 
!

subroutine  DeAllocCellsArrsDGV
use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

deallocate(cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov)

end subroutine DeAllocCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine deallocates the cells arrays 
!

subroutine  DeAllocCellsArrsDGVII
use DGV_commvar, only: cells_pgridII, cells_cgridII, cells_luII, cells_lvII, cells_lwII, & 
                   cells_ruII, cells_rvII, cells_rwII, cells_refuII, cells_refvII, & 
                   cells_refwII, cells_gowII, cells_gouII, cells_govII

deallocate(cells_pgridII, cells_cgridII, cells_luII, cells_lvII, cells_lwII, & 
                   cells_ruII, cells_rvII, cells_rwII, cells_refuII, cells_refvII, & 
                   cells_refwII, cells_gowII, cells_gouII, cells_govII)

end subroutine DeAllocCellsArrsDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FillCellsArrsDGV
! 
! This subroutine populates the certain portion of the cells_arrays.
! 
! Expects to have cells_arrays to exists and to be of a proper size
! 
! cells_cgrid = number of the grid. If the cell is not refined then the value is -1, if the cell is refined than 
!               cells_cgrid gives the number of the grid. 

subroutine FillCellsArrsDGV(ibeg,iend,umesh,vmesh,wmesh,pgrid,gou,gov,gow)

use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

integer (I4B), intent(in) :: ibeg,iend ! beginning and end of the range to fill the cells array
real (DP), dimension (0:), intent (in) :: umesh,vmesh,wmesh ! 1D meshes to be used
integer (I4B), intent (in) :: pgrid ! the number of the parent grid --- grid where these cells belong. 
integer (I4B), intent (in) :: gou,gov,gow ! the max order of gauss lagrange basis to be assigned to the cells -- all cells in this range will be assigned the same order in all 
! 
integer (I4B) :: iu,iv,iw,mm ! local counters
! first some elementary checks
if ((size(cells_pgrid,1)< ibeg) .or. (size(cells_pgrid,1)< iend)) then 
 print *, "FillCellsArrsDGV: supplied range does not work for the cells arrays. "
 stop
end if
if ((ibeg > iend) .or. (iend-ibeg + 1  /= (size(vmesh,1)-1)*(size(wmesh,1)-1)*(size(umesh,1)-1)))   then 
 print *, "FillCellsArrsDGV: supplied range does not work: ibeg>iend or iend-ibeg + 1 /= size(vmesh)*size(wmesh)*size(umesh). "
 stop
end if  
! ok, let us populate the cells
! some properties will be the same for the entire range: 
cells_pgrid(ibeg:iend) = pgrid
cells_cgrid(ibeg:iend) = -1
cells_refu(ibeg:iend) = 1
cells_refv(ibeg:iend) = 1
cells_refw(ibeg:iend) = 1
cells_gou(ibeg:iend) = gou
cells_gov(ibeg:iend) = gov
cells_gow(ibeg:iend) = gow
! some properties will be different: 
mm=ibeg
do iu=0,size(umesh,1)-2
do iv=0,size(vmesh,1)-2
do iw=0,size(wmesh,1)-2
   cells_lu(mm) = umesh(iu)
   cells_lv(mm) = vmesh(iv)
   cells_lw(mm) = wmesh(iw)
   cells_ru(mm) = umesh(iu+1)
   cells_rv(mm) = vmesh(iv+1)
   cells_rw(mm) = wmesh(iw+1)
   mm=mm+1
end do 
end do 
end do 
end subroutine FillCellsArrsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FillCellsArrsDGVII
! 
! This subroutine populates the certain portion of the cells_arrays.
! 
! 
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! Expects to have cells_arraysII to exists and to be of a proper size
! 
! cells_cgridII = number of the grid. If the cell is not refined then the value is -1, if the cell is refined than 
!               cells_cgrid gives the number of the grid. 
!!!!!!!!!!!!!!!!!!!!!

subroutine FillCellsArrsDGVII(ibeg,iend,umesh,vmesh,wmesh,pgrid,gou,gov,gow)

use DGV_commvar, only: cells_pgrid=>cells_pgridII, cells_cgrid=>cells_cgridII, &
                   cells_lu=>cells_luII, cells_lv=>cells_lvII, cells_lw=>cells_lwII, & 
                   cells_ru=>cells_ruII, cells_rv=>cells_rvII, cells_rw=>cells_rwII, &
                   cells_refu=>cells_refuII, cells_refv=>cells_refvII, cells_refw=>cells_refwII, & 
                   cells_gow=>cells_gowII, cells_gou=>cells_gouII, cells_gov=>cells_govII

integer (I4B), intent(in) :: ibeg,iend ! beginning and end of the range to fill the cells array
real (DP), dimension (0:), intent (in) :: umesh,vmesh,wmesh ! 1D meshes to be used
integer (I4B), intent (in) :: pgrid ! the number of the parent grid --- grid where these cells belong. 
integer (I4B), intent (in) :: gou,gov,gow ! the max order of gauss lagrange basis to be assigned to the cells -- all cells in this range will be assigned the same order in all 
! 
integer (I4B) :: iu,iv,iw,mm ! local counters
! first some elementary checks
if ((size(cells_pgrid,1)< ibeg) .or. (size(cells_pgrid,1)< iend)) then 
 print *, "FillCellsArrsDGVII: supplied range does not work for the cells arrays. "
 stop
end if
if ((ibeg > iend) .or. (iend-ibeg + 1  /= (size(vmesh,1)-1)*(size(wmesh,1)-1)*(size(umesh,1)-1)))   then 
 print *, "FillCellsArrsDGVII: supplied range does not work: ibeg>iend or iend-ibeg + 1 .neqv. size(vmesh)*size(wmesh)*size(umesh). "
 stop
end if  
! ok, let us populate the cells
! some properties will be the same for the entire range: 
cells_pgrid(ibeg:iend) = pgrid
cells_cgrid(ibeg:iend) = -1
cells_refu(ibeg:iend) = 1
cells_refv(ibeg:iend) = 1
cells_refw(ibeg:iend) = 1
cells_gou(ibeg:iend) = gou
cells_gov(ibeg:iend) = gov
cells_gow(ibeg:iend) = gow
! some properties will be different: 
mm=ibeg
do iu=0,size(umesh,1)-2
do iv=0,size(vmesh,1)-2
do iw=0,size(wmesh,1)-2
   cells_lu(mm) = umesh(iu)
   cells_lv(mm) = vmesh(iv)
   cells_lw(mm) = wmesh(iw)
   cells_ru(mm) = umesh(iu+1)
   cells_rv(mm) = vmesh(iv+1)
   cells_rw(mm) = wmesh(iw+1)
   mm=mm+1
end do 
end do 
end do 
end subroutine FillCellsArrsDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!! CellsRefineDGV
!! 
!! This subrotine refines cells that are on the supplied list
!!
!! Depends on the main program. The cells arrays must be set up before calling this subroutine
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CellsRefineDGV(reflist,refu,refv,refw)

use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w

integer (I4B), dimension (:), intent (in) :: reflist  !list of cells to refine at this step
integer (I4B), intent(in) :: refu,refv,refw ! refinement factors in dimension u, v, and w
!!!!!!!!
!! IMPORTANT: 
!!
integer (I4B) :: go_inc=0 !! THIS IS A PARAMETER: the refined cells will have orders of the parent + go_inc
!!
!!!!!!!! 
integer (I4B) :: i,ic,mm ! local counters 
integer (I4B) :: isc ! the number of the cell that is being refined
integer (I4B) :: cgrid ! a scrap number to keep the number of the new grid
integer (I4B) :: gou,gov,gow ! the scrap variables to keep the gauss order of the basis
integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

real (DP), dimension (:), allocatable :: dumu,dumv,dumw 
! a quick check for errors: 
if ((refu*refv*refw < 1 ) .or. (size(reflist,1)<1)) then 
  print *,"CellsRefineDGV: error in incoming parameters, refu*refv*refw < 1 or size(reflist,1)<1"
end if 
! also, we need to have a shorthand for how many cells we will add
mm=(refu)*(refv)*(refw)
!
do i=1,size(reflist,1)
  isc = reflist(i)
  ! We need to construnct the new refined grids
  allocate (dumu(1:refu+1),dumv(1:refv+1),dumw(1:refw+1), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "CellsRefineDGV: Allocation error for variable (dumu),(dumv),or (dumw)"
  stop
  end if
     !
  dumu = (/ (cells_lu(isc) + (cells_ru(isc) - cells_lu(isc))/Real(refu,DP)*i, i=0,refu) /)
  dumv = (/ (cells_lv(isc) + (cells_rv(isc) - cells_lv(isc))/Real(refv,DP)*i, i=0,refv) /)
  dumw = (/ (cells_lw(isc) + (cells_rw(isc) - cells_lw(isc))/Real(refw,DP)*i, i=0,refw) /)
  ! Now we need to extend the grids arrays to include the new grids... 
  call ExtendGridsDGV(dumu,dumv,dumw,cgrid)
  ! Now we mark the cell with number (isc) to have a child grid...
  cells_cgrid(isc) = cgrid
  cells_refu(isc)= refu
  cells_refv(isc)= refv
  cells_refw(isc)= refw
  ! We also want to know what is the order of the Lagrange basis on that cell: 
  gou=cells_gou(isc)
  gov=cells_gov(isc)
  gow=cells_gow(isc)
  ! we need to know how long was cells before we added new: 
  ic=size(cells_lu,1)
  ! Now we need to extend the cells arrays to include the new cells... 
  call ExtendCellsArrsDGV(ic+mm)
  ! now we need to populate the new cells
  call FillCellsArrsDGV(ic+1,ic+mm,dumu,dumv,dumw,cgrid,gou+go_inc,gov+go_inc,gow+go_inc)
  ! ready to do the next cell 
end do 


end subroutine CellsRefineDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendGridsDGV
!
! This subroutine extends the grids_cap_u/_v/_w and the grids_u/_v/_w arrays
! Depends on the main program. 
!
! dumu,dumv,dumw -- these are the new grids that need to be added to the old grids... 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendGridsDGV(dumu,dumv,dumw,cgrid)

use DGV_commvar, only: grids_cap_u, grids_cap_v, grids_cap_w, grids_u, grids_v, grids_w

real (DP), dimension (:), intent (in) :: dumu,dumv,dumw ! the new velocity values that will be added to grids
integer (I4B), intent (out) :: cgrid ! the number of the child grid that will be created 
!
integer (I4B) :: iu,iv,iw ! some scrap variables... 
integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

real (DP), dimension (:), allocatable :: tmp_u,tmp_v,tmp_w ! temporary arrays to keep the old grids

! a quick check for errors
if ((size(dumu,1)< 3 ) .or. (size(dumv,1)< 3) .or. (size(dumw,1)< 3 )) then 
  print *,"ExtendGridsDGV: error in incoming parameters, at least one of the arrays does not require a refinement < 1"
end if 
! create the dump arrays to save current grids_cap_u/v/w:
iu=size(grids_cap_u,1)
allocate (tmp_u(1:iu),tmp_v(1:iu),tmp_w(1:iu), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error 1 for variable (tmp_u),(tmp_v),or (tmp_w)"
  stop
  end if
! save the old grids_cap arrays
tmp_u=grids_cap_u
tmp_v=grids_cap_v
tmp_w=grids_cap_w
! deallocate the grids_cap 
deallocate (grids_cap_u,grids_cap_v,grids_cap_w)
! allocate them 1 cell longer... 
allocate (grids_cap_u(1:iu+1),grids_cap_v(1:iu+1),grids_cap_w(1:iu+1), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error for variable (grids_cap_u/_v/_w)"
  stop
  end if
! restore the old grids_cap arrays
grids_cap_u(1:iu)=tmp_u
grids_cap_v(1:iu)=tmp_v
grids_cap_w(1:iu)=tmp_w
! record the new cells
grids_cap_u(iu+1)=size(dumu,1)
grids_cap_v(iu+1)=size(dumv,1)
grids_cap_w(iu+1)=size(dumw,1)
! done with grids_cap arrays
cgrid = iu+1 ! this will be passed out as the number of the new created grid
deallocate(tmp_u,tmp_v,tmp_w)
! now we extend the grids_u/_v/_w arrays... 
iu=size(grids_u,1)
iv=size(grids_v,1)
iw=size(grids_w,1)
allocate (tmp_u(1:iu),tmp_v(1:iv),tmp_w(1:iw), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error 2 for variable  (tmp_u),(tmp_v),or (tmp_w)"
  stop
  end if
! save the old grids_u/_v/_w arrays
tmp_u=grids_u
tmp_v=grids_v
tmp_w=grids_w
! deallocate the grids_cap 
deallocate(grids_u,grids_v,grids_w)
! allocate them again but longer... 
allocate (grids_u(1:iu+size(dumu,1)),grids_v(1:iv+size(dumv,1)),grids_w(1:iw+size(dumw,1)), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendGridsDGV: Allocation error for variable (grids_cap_u/_v/_w)"
  stop
  end if
! restore the old grids_cap arrays
grids_u(1:iu)=tmp_u
grids_v(1:iv)=tmp_v
grids_w(1:iw)=tmp_w
! record the new cells
grids_u(iu+1:iu+size(dumu,1)) = dumu
grids_v(iv+1:iv+size(dumv,1)) = dumv
grids_w(iw+1:iw+size(dumw,1)) = dumw
! done with grids_cap arrays
deallocate(tmp_u,tmp_v,tmp_w)
end subroutine ExtendGridsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine SetNodesDGV 
!
! This subroutine sets all nodes in velocity that will ever be used. 
!
! nodes_cells contains the # of the cell where this node belongs
! nodex_x/_y/_z contains the coordinates of the nodes.
! nodes_weights contains the product of the weights of gauss quadrature for each node so that the 
!               functions can be evaluated on the nodes and multiplied by the weigth and summed. 
!
! grids_cap_u -- stores the grid capasity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetNodesDGV

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w, nodes_gwts,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   g_wts_all, nodes_ui, nodes_vi, nodes_wi, &
				   nodes_ut, nodes_vt, nodes_wt
                   
intrinsic SUM 

integer (I4B) :: ic,im,jm,iu,iv,iw ! some scrap variables... 
real (DP) :: ru,qu,rv,qv,rw,qw ! some scrap variables ...

integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

im=0
do ic=1,size(cells_pgrid,1)
 if (cells_cgrid(ic) == -1) then 
  im=im+cells_gou(ic)*cells_gov(ic)*cells_gow(ic)
 end if
end do 
!
! We now allocating the nodes arrays
allocate (nodes_pcell(1:im),nodes_ui(1:im),nodes_vi(1:im),nodes_wi(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGV: Allocation error for variable (nodes_pcell,nodes_ui/_vi/_wi)"
  stop
  end if

allocate (nodes_ut(1:im),nodes_vt(1:im),nodes_wt(1:im),stat=loc_alloc_stat)
	if (loc_alloc_stat >0) then
	print *, "SetNodesArDGV: Allocation error for variables (nodes_ut, nodes_vt, nodes_wt)"
	stop
	end if

allocate (nodes_u(1:im),nodes_v(1:im),nodes_w(1:im),nodes_gwts(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGV: Allocation error for variables (nodes_u/_v/_w/_gwts)"
  stop
  end if
! Now we will populate them. 
jm=1
do ic=1,size(cells_pgrid,1)
if (cells_cgrid(ic) == -1) then 
 do iu=1,cells_gou(ic)
   ru=g_nds_all(iu,cells_gou(ic))*(cells_ru(ic)-cells_lu(ic))/2.0_DP + (cells_ru(ic)+cells_lu(ic))/2.0_DP
   qu=g_wts_all(iu,cells_gou(ic))
 do iv=1,cells_gov(ic)
   rv=g_nds_all(iv,cells_gov(ic))*(cells_rv(ic)-cells_lv(ic))/2.0_DP + (cells_rv(ic)+cells_lv(ic))/2.0_DP
   qv=g_wts_all(iv,cells_gov(ic))
 do iw=1,cells_gow(ic)
   rw=g_nds_all(iw,cells_gow(ic))*(cells_rw(ic)-cells_lw(ic))/2.0_DP + (cells_rw(ic)+cells_lw(ic))/2.0_DP
   qw=g_wts_all(iw,cells_gow(ic))
! ready to fill the nodes! 
nodes_pcell(jm) = ic   
nodes_u(jm)= ru 
nodes_v(jm)= rv
nodes_w(jm)= rw
nodes_gwts(jm)=qu*qv*qw*(cells_ru(ic)-cells_lu(ic))*(cells_rv(ic)-cells_lv(ic))*(cells_rw(ic)-cells_lw(ic))/8.0_DP
nodes_ui(jm)=iu
nodes_vi(jm)=iv
nodes_wi(jm)=iw
jm=jm+1
 end do 
 end do 
 end do 
end if  
end do 
!
if (jm-1 /= im) then 
  print *, "SetNodesArDGV: Error, jm-1 .neqv. im. Nodes arrays may be populated wrong."
  stop
  end if

end subroutine SetNodesDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine SetNodesDGVII 
!
! This subroutine sets all nodes in velocity in secondary grid that will ever be used. 
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! nodes_cellsII contains the # of the cell where this node belongs
! nodex_x/_y/_zII contains the coordinates of the nodes.
! nodes_weightsII contains the product of the weights of gauss quadrature for each node so that the 
!               functions can be evaluated on the nodes and multiplied by the weigth and summed. 
!
! grids_cap_uII -- stores the grid capasity 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetNodesDGVII

use DGV_commvar, only: nodes_pcell=>nodes_pcellII, &
                   nodes_u=>nodes_uII, nodes_v=>nodes_vII, nodes_w=>nodes_wII, &
                   nodes_gwts=>nodes_gwtsII, cells_pgrid=>cells_pgridII, cells_cgrid=>cells_cgridII, &
                   cells_lu=>cells_luII, cells_lv=>cells_lvII, cells_lw=>cells_lwII, & 
                   cells_ru=>cells_ruII, cells_rv=>cells_rvII, cells_rw=>cells_rwII, & 
                   cells_refu=>cells_refuII, cells_refv=>cells_refvII, cells_refw=>cells_refwII, &
                   cells_gow=>cells_gowII, cells_gou=>cells_gouII, cells_gov=>cells_govII, &
                   grids_cap_u=>grids_cap_uII, grids_cap_v=>grids_cap_vII, grids_cap_w=>grids_cap_wII, &
                   grids_u=>grids_uII, grids_v=>grids_vII, grids_w=>grids_wII,&
                   g_nds_all, g_wts_all, &
                   nodes_ui=>nodes_uiII, nodes_vi=>nodes_viII, nodes_wi=>nodes_wiII, &
				   nodes_ut=>nodes_utII, nodes_vt=>nodes_vtII, nodes_wt=>nodes_wtII
                   
intrinsic SUM 

integer (I4B) :: ic,im,jm,iu,iv,iw ! some scrap variables... 
real (DP) :: ru,qu,rv,qv,rw,qw ! some scrap variables ...

integer (I4B) :: loc_alloc_stat  ! local variable to keep the allocation status

im=0
do ic=1,size(cells_pgrid,1)
 if (cells_cgrid(ic) == -1) then 
  im=im+cells_gou(ic)*cells_gov(ic)*cells_gow(ic)
 end if
end do 
!
! We now allocating the nodes arrays
allocate (nodes_pcell(1:im),nodes_ui(1:im),nodes_vi(1:im),nodes_wi(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGVII: Allocation error for variable (nodes_pcell,nodes_ui/_vi/_wiII)"
  stop
  end if

allocate (nodes_ut(1:im),nodes_vt(1:im),nodes_wt(1:im),stat=loc_alloc_stat)
	if (loc_alloc_stat >0) then
	print *, "SetNodesArDGVII: Allocation error for variables (nodes_utII, nodes_vtII, nodes_wtII)"
	stop
	end if

allocate (nodes_u(1:im),nodes_v(1:im),nodes_w(1:im),nodes_gwts(1:im),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetNodesArDGVII: Allocation error for variables (nodes_u/_v/_w/_gwtsII)"
  stop
  end if
! Now we will populate them. 
jm=1
do ic=1,size(cells_pgrid,1)
if (cells_cgrid(ic) == -1) then 
 do iu=1,cells_gou(ic)
   ru=g_nds_all(iu,cells_gou(ic))*(cells_ru(ic)-cells_lu(ic))/2.0_DP + (cells_ru(ic)+cells_lu(ic))/2.0_DP
   qu=g_wts_all(iu,cells_gou(ic))
 do iv=1,cells_gov(ic)
   rv=g_nds_all(iv,cells_gov(ic))*(cells_rv(ic)-cells_lv(ic))/2.0_DP + (cells_rv(ic)+cells_lv(ic))/2.0_DP
   qv=g_wts_all(iv,cells_gov(ic))
 do iw=1,cells_gow(ic)
   rw=g_nds_all(iw,cells_gow(ic))*(cells_rw(ic)-cells_lw(ic))/2.0_DP + (cells_rw(ic)+cells_lw(ic))/2.0_DP
   qw=g_wts_all(iw,cells_gow(ic))
! ready to fill the nodes! 
nodes_pcell(jm) = ic   
nodes_u(jm)= ru 
nodes_v(jm)= rv
nodes_w(jm)= rw
nodes_gwts(jm)=qu*qv*qw*(cells_ru(ic)-cells_lu(ic))*(cells_rv(ic)-cells_lv(ic))*(cells_rw(ic)-cells_lw(ic))/8.0_DP
nodes_ui(jm)=iu
nodes_vi(jm)=iv
nodes_wi(jm)=iw
jm=jm+1
 end do 
 end do 
 end do 
end if  
end do 
!
if (jm-1 .ne. im) then 
  print *, "SetNodesArDGVII: Error, jm-1 /= im. Nodes arrays may be populated wrong."
  stop
  end if

end subroutine SetNodesDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetCellsUGIAshiftArrays
!
!
! This subroutine is setting up the arrays cells_gui, _gvi, _gwi that 
! store the relative integer adress of the cell on the grid. 
!
! Only wil work if there is only one grid in the mesh.
!
! Also this subroutine sets up an arrays that gives the adress shift for different nodes to 
! read staff from the A-Array.
!
! Note that all these arrays can be created in other places of the algorithm, but we did nto think about them 
! earlier. So we put them here.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetCellsUGINdsSftArrs(ccell)

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi, cells_ugi, cells_vgi, cells_wgi, &
             cells_pgrid, grids_cap_u, grids_cap_v, grids_cap_w, nodes_Ashift, nodes_phican,&
             cells_gou,cells_gov,cells_gow,nodes_dui,nodes_dvi,nodes_dwi,A_capphi
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,Ashift,phishift,cshift ! scrap variable...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer (I4B) :: phicell ! the number of the cell where phi belongs
integer (I4B) :: gou,gov,gow ! scrap variables to keep the number of nodes in a cell and the shift in nodes to the canonical cell
integer (I4B) :: pgcu,pgcv,pgcw ! scrap variables to keep the number of cells on the grid in each dimension
integer (I4B) :: ccell_ugi,ccell_vgi,ccell_wgi ! scrap variables to keep the integer address of the canonical cell on the grid
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"EvalCollisionPeriodicA_DGV: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
if (size(grids_cap_u,1) /= 1) then 
     print *, "SetCellsUGIAshiftArrays: Error. The number of grids must be 1. Stop"
     stop
end if  
!!! First, we need to calculate some useful constants:
gou=cells_gou(ccell)
gov=cells_gov(ccell)
gow=cells_gow(ccell)
cshift=(ccell-1)*gou*gov*gow ! this will be the number of the node right before the first node in canonical cell
! first we find the number of cells in each dimension
pgcu=grids_cap_u(cells_pgrid(ccell))-1
pgcv=grids_cap_v(cells_pgrid(ccell))-1
pgcw=grids_cap_w(cells_pgrid(ccell))-1
! now we will calculate the relative integer address of all cells on the grid.
! first, we create a space where to store this information: 
i=size(cells_pgrid,1)
allocate (cells_ugi(1:i),cells_vgi(1:i),cells_wgi(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArrays: Allocation error for variable (cells_ugi/_vgi/_wgi)"
  stop
  end if
! Nex we will calculate the adresses
! IMPORTANT! Notice that this algorithm uses the convention that when enumerating cells the outside loop is in x, the 
! middle loop is in y and the innermost is in z.  

do phicell = 1,pgcu*pgcv*pgcw ! loop in all cells
!
ccell_ugi=1
ccell_vgi=1
ccell_wgi=1
i=1
do while (i<phicell)  
 i=i+1
 if (ccell_wgi == pgcw) then 
   ccell_wgi = 1 
   if (ccell_vgi == pgcv) then 
     ccell_vgi = 1
     if (ccell_ugi == pgcu) then  
      print *,"EvalCollisionPeriodicA_DGV: Error the index of the canonical cell is not fouind. Value of (gui) is too big."
      stop
     else 
      ccell_ugi=ccell_ugi+1
     end if
   else 
     ccell_vgi=ccell_vgi+1
   end if     
 else
   ccell_wgi = ccell_wgi+1
 end if           
end do 
! finshed finding the cells addresses, now recording:
cells_ugi(phicell) = ccell_ugi
cells_vgi(phicell) = ccell_vgi
cells_wgi(phicell) = ccell_wgi
end do 
! finished finding the local adress of the cells in the case of one grid only 
! Next find the nodes_Ashift and nodes_dui/_dvi/_dwi arrays. 
! First we allocated them... 
i=size(nodes_vi,1)
allocate (nodes_Ashift(1:i),nodes_dui(1:i),nodes_dvi(1:i),nodes_dwi(1:i),nodes_phican(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArrays: Allocation error for variable (nodes_Ashift,nodes_dui/dvi/dwi)"
  stop
  end if
! next we populate the nodes_dui/_dvi/dwi -arrays... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_vi,1)
    phicell = nodes_pcell(iphi)
    if (cells_pgrid(phicell) /= cells_pgrid(ccell)) then 
     print *, "SetCellsUGIAshiftArrays: Error. The parent cell of phi belongs to a different grid than the canonical cell. Stop"
     stop
    end if 
    ! next we record the shifts in integer index in from the cell containting phi to the canonical
    nodes_dui(iphi)=cells_ugi(phicell) - cells_ugi(ccell)
    nodes_dvi(iphi)=cells_vgi(phicell) - cells_vgi(ccell)
    nodes_dwi(iphi)=cells_wgi(phicell) - cells_wgi(ccell)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_Ashift number and the nodes_phican number 
    ! First we calculate the number of the node on the canonical cell that corresponds to the node with number iphi 
    phishift = cshift + (nodes_ui(iphi)-1)*gov*gow + (nodes_vi(iphi)-1)*gow + nodes_wi(iphi) ! this is the number of the node
    nodes_phican(iphi) = phishift ! we will record this number. 
    ! Now we will calculate the Ashift number
    !! ATTENTION: when no A file present, the next line generates an access violation code! Program crashes
    Ashift = sum(A_capphi(1:phishift-1))
    nodes_Ashift(iphi)=Ashift
    !!!!! all done
end do 
end subroutine SetCellsUGINdsSftArrs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetCellsUGIAshiftArraysII
!
!
! This subroutine is setting up the arrays cells_guiII, _gviII, _gwiII that 
! store the relative integer adress of the cell on the grid. 
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
!
! Only wil work if there is only one grid in the mesh.
!
! Also this subroutine sets up an arrays that gives the adress shift for different nodes to 
! read staff from the A-Array.
!
! Note that all these arrays can be created in other places of the algorithm, but we did nto think about them 
! earlier. So we put them here.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetCellsUGINdsSftArrsII(ccell)

use DGV_commvar, only: nodes_pcell=>nodes_pcellII, &
    nodes_ui=>nodes_uiII, nodes_vi=>nodes_viII, nodes_wi=>nodes_wiII, & 
    cells_ugi=>cells_ugiII, cells_vgi=>cells_vgiII, cells_wgi=>cells_wgiII, &
    cells_pgrid=>cells_pgridII, & 
    grids_cap_u=>grids_cap_uII, grids_cap_v=>grids_cap_vII, grids_cap_w=>grids_cap_wII, &
    nodes_Ashift=>nodes_AshiftII, nodes_phican=>nodes_phicanII, A_capphi=>A_capphiII, &
    cells_gou=>cells_gouII, cells_gov=>cells_govII, cells_gow=>cells_gowII, &
    nodes_dui=>nodes_duiII, nodes_dvi=>nodes_dviII, nodes_dwi=>nodes_dwiII
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,Ashift,phishift,cshift ! scrap variable...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer (I4B) :: phicell ! the number of the cell where phi belongs
integer (I4B) :: gou,gov,gow ! scrap variables to keep the number of nodes in a cell and the shift in nodes to the canonical cell
integer (I4B) :: pgcu,pgcv,pgcw ! scrap variables to keep the number of cells on the grid in each dimension
integer (I4B) :: ccell_ugi,ccell_vgi,ccell_wgi ! scrap variables to keep the integer address of the canonical cell on the grid
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"EvalCollisionPeriodicA_DGVII: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
if (size(grids_cap_u,1) /= 1) then 
     print *, "SetCellsUGIAshiftArraysII: Error. The number of grids must be 1. Stop"
     stop
end if  
!!! First, we need to calculate some useful constants:
gou=cells_gou(ccell)
gov=cells_gov(ccell)
gow=cells_gow(ccell)
cshift=(ccell-1)*gou*gov*gow ! this will be the number of the node right before the first node in canonical cell
! first we find the number of cells in each dimension
pgcu=grids_cap_u(cells_pgrid(ccell))-1
pgcv=grids_cap_v(cells_pgrid(ccell))-1
pgcw=grids_cap_w(cells_pgrid(ccell))-1
! now we will calculate the relative integer address of all cells on the grid.
! first, we create a space where to store this information: 
i=size(cells_pgrid,1)
allocate (cells_ugi(1:i),cells_vgi(1:i),cells_wgi(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArraysII: Allocation error for variable (cells_ugi/_vgi/_wgiII)"
  stop
  end if
! Nex we will calculate the adresses
! IMPORTANT! Notice that this algorithm uses the convention that when enumerating cells the outside loop is in x, the 
! middle loop is in y and the innermost is in z.  

do phicell = 1,pgcu*pgcv*pgcw ! loop in all cells
!
ccell_ugi=1
ccell_vgi=1
ccell_wgi=1
i=1
do while (i<phicell)  
 i=i+1
 if (ccell_wgi == pgcw) then 
   ccell_wgi = 1 
   if (ccell_vgi == pgcv) then 
     ccell_vgi = 1
     if (ccell_ugi == pgcu) then  
      print *,"EvalCollisionPeriodicA_DGVII: Error the index of the canonical cell is not fouind. Value of (ugi) is too big."
      stop
     else 
      ccell_ugi=ccell_ugi+1
     end if
   else 
     ccell_vgi=ccell_vgi+1
   end if     
 else
   ccell_wgi = ccell_wgi+1
 end if           
end do 
! finshed finding the cells addresses, now recording:
cells_ugi(phicell) = ccell_ugi
cells_vgi(phicell) = ccell_vgi
cells_wgi(phicell) = ccell_wgi
end do 
! finished finding the local adress of the cells

! Next find the nodes_Ashift and nodes_dui/_dvi/_dwi arrays. 
! First we allocated them... 
i=size(nodes_vi,1)
allocate (nodes_Ashift(1:i),nodes_dui(1:i),nodes_dvi(1:i),nodes_dwi(1:i),nodes_phican(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGIAshiftArraysII: Allocation error for variable (nodes_Ashift,nodes_dui/dvi/dwiII)"
  stop
  end if
! next we populate the nodes_dui/_dvi/dwi -arrays... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_vi,1)
    phicell = nodes_pcell(iphi)
    if (cells_pgrid(phicell) /= cells_pgrid(ccell)) then 
     print *, "SetCellsUGIAshiftArraysII: Error. The parent cell of phi belongs to a different grid than the canonical cell. Stop"
     stop
    end if 
    ! next we record the shifts in integer index in from the cell containting phi to the canonical
    nodes_dui(iphi)=cells_ugi(phicell) - cells_ugi(ccell)
    nodes_dvi(iphi)=cells_vgi(phicell) - cells_vgi(ccell)
    nodes_dwi(iphi)=cells_wgi(phicell) - cells_wgi(ccell)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_Ashift number and the nodes_phican number 
    ! First we calculate the number of the node on the canonical cell that corresponds to the node with number iphi 
    phishift = cshift + (nodes_ui(iphi)-1)*gov*gow + (nodes_vi(iphi)-1)*gow + nodes_wi(iphi) ! this is the number of the node
    nodes_phican(iphi) = phishift ! we will record this number. 
    ! Now we will calculate the Ashift number
    Ashift = sum(A_capphi(1:phishift-1))
    nodes_Ashift(iphi)=Ashift
    !!!!! all done
end do 
end subroutine SetCellsUGINdsSftArrsII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetCellsUGI_DGV
!
!
! This subroutine is setting up the arrays cells_gui, _gvi, _gwi that 
! store the relative integer adress of the cell on the grid. 
!
! Only will work if there is only one grid in the mesh.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetCellsUGI_DGV(ccell)

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi, cells_ugi, cells_vgi, cells_wgi, &
             cells_pgrid, grids_cap_u, grids_cap_v, grids_cap_w, nodes_phican,&
             cells_gou,cells_gov,cells_gow,nodes_dui,nodes_dvi,nodes_dwi
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,phishift,cshift ! scrap variable...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer (I4B) :: phicell ! the number of the cell where phi belongs
integer (I4B) :: gou,gov,gow ! scrap variables to keep the number of nodes in a cell and the shift in nodes to the canonical cell
integer (I4B) :: pgcu,pgcv,pgcw ! scrap variables to keep the number of cells on the grid in each dimension
integer (I4B) :: ccell_ugi,ccell_vgi,ccell_wgi ! scrap variables to keep the integer address of the canonical cell on the grid
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"SetCellsUGI_DGV: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
if (size(grids_cap_u,1) /= 1) then 
     print *, "SetCellsUGI_DGV: Error. The number of grids must be 1. Stop"
     stop
end if  
!!! First, we need to calculate some useful constants:
gou=cells_gou(ccell)
gov=cells_gov(ccell)
gow=cells_gow(ccell)
cshift=(ccell-1)*gou*gov*gow ! this will be the number of the node right before the first node in canonical cell
! first we find the number of cells in each dimension
pgcu=grids_cap_u(cells_pgrid(ccell))-1
pgcv=grids_cap_v(cells_pgrid(ccell))-1
pgcw=grids_cap_w(cells_pgrid(ccell))-1
! now we will calculate the relative integer address of all cells on the grid.
! first, we create a space where to store this information: 
i=size(cells_pgrid,1)
allocate (cells_ugi(1:i),cells_vgi(1:i),cells_wgi(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGI_DGV: Allocation error for variable (cells_ugi/_vgi/_wgi)"
  stop
  end if
! Nex we will calculate the adresses
! IMPORTANT! Notice that this algorithm uses the convention that when enumerating cells the outside loop is in x, the 
! middle loop is in y and the innermost is in z.  

do phicell = 1,pgcu*pgcv*pgcw ! loop in all cells
!
ccell_ugi=1
ccell_vgi=1
ccell_wgi=1
i=1
do while (i<phicell)  
 i=i+1
 if (ccell_wgi == pgcw) then 
   ccell_wgi = 1 
   if (ccell_vgi == pgcv) then 
     ccell_vgi = 1
     if (ccell_ugi == pgcu) then  
      print *,"SetCellsUGI_DGV: Error the index of the canonical cell is not fouind. Value of (gui) is too big."
      stop
     else 
      ccell_ugi=ccell_ugi+1
     end if
   else 
     ccell_vgi=ccell_vgi+1
   end if     
 else
   ccell_wgi = ccell_wgi+1
 end if           
end do 
! finshed finding the cells addresses, now recording:
cells_ugi(phicell) = ccell_ugi
cells_vgi(phicell) = ccell_vgi
cells_wgi(phicell) = ccell_wgi
end do 
! finished finding the local adress of the cells in the case of one grid only 
! Next find the nodes_Ashift and nodes_dui/_dvi/_dwi arrays. 
! First we allocated them... 
i=size(nodes_vi,1)
allocate (nodes_dui(1:i),nodes_dvi(1:i),nodes_dwi(1:i),nodes_phican(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetCellsUGI_DGV: Allocation error for variable (nodes_dui/dvi/dwi)"
  stop
  end if
! next we populate the nodes_dui/_dvi/dwi -arrays... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_vi,1)
    phicell = nodes_pcell(iphi)
    if (cells_pgrid(phicell) /= cells_pgrid(ccell)) then 
     print *, "SetCellsUGI_DGV: Error. The parent cell of phi belongs to a different grid than the canonical cell. Stop"
     stop
    end if 
    ! next we record the shifts in integer index in from the cell containting phi to the canonical
    nodes_dui(iphi)=cells_ugi(phicell) - cells_ugi(ccell)
    nodes_dvi(iphi)=cells_vgi(phicell) - cells_vgi(ccell)
    nodes_dwi(iphi)=cells_wgi(phicell) - cells_wgi(ccell)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_Ashift number and the nodes_phican number 
    ! First we calculate the number of the node on the canonical cell that corresponds to the node with number iphi 
    phishift = cshift + (nodes_ui(iphi)-1)*gov*gow + (nodes_vi(iphi)-1)*gow + nodes_wi(iphi) ! this is the number of the node
    nodes_phican(iphi) = phishift ! we will record this number. 
    !!!!! all done
end do 
end subroutine SetCellsUGI_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetAshift_DGV
!
!
! This subroutine is setting up the arrays nodes_Ashift  
!
! Only will work if there is only one grid in the mesh.
!
! Also this subroutine sets up an arrays that gives the adress shift for different nodes to 
! read staff from the A-Array.
!
! Note that all these arrays can be created in other places of the algorithm, but we did nto think about them 
! earlier. So we put them here.
!
! ATTENTION: This subroutine depends on arrays cells_gui,_gvi,_gwi and the nodes_dui/_dvi/_dwi
!            Make sure SetCellsUGI_DGV is called before this subroutine is called 
! 
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetAshift_DGV(ccell)

use DGV_commvar, only: nodes_pcell, cells_pgrid, nodes_Ashift, nodes_phican, A_capphi
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,Ashift,phishift ! scrap variable...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer (I4B) :: phicell ! the number of the cell where phi belongs
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"SetAshift_DGV: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
! first, we create a space where to store Ashift information: 
i=size(nodes_pcell,1)
allocate (nodes_Ashift(1:i),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetAshift_DGV: Allocation error for variable (nodes_Ashift)"
  stop
  end if
! next we populate the Ashift array 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_pcell,1)
    phicell = nodes_pcell(iphi)
    if (cells_pgrid(phicell) /= cells_pgrid(ccell)) then 
     print *, "SetAshift_DGV: Error. The parent cell of phi belongs to a different grid than the canonical cell. Stop"
     stop
    end if 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_Ashift number 
    !! ATTENTION: when no A file present, the next line generates an access violation code! Program crashes
    phishift = nodes_phican(iphi)
    Ashift = sum(A_capphi(1:phishift-1))
    nodes_Ashift(iphi) = Ashift
    !!!!! all done
end do 
end subroutine SetAshift_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetAkorshiftAllNets_DGV
!
!
! This subroutine is setting up the array nodes_AkorshiftAllNets  
!
! Only will work if there is only one grid in the mesh.
!
! Also this subroutine sets up an arrays that gives the adress shift for different nodes to 
! read staff from the A-Array.
!
! Note that all these arrays can be created in other places of the algorithm, but we did nto think about them 
! earlier. So we put them here.
!
! ATTENTION: This subroutine depends on arrays cells_gui,_gvi,_gwi and the nodes_dui/_dvi/_dwi
!            Make sure SetCellsUGI_DGV is called before this subroutine is called 
! 
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetAkorshiftAllNets_DGV(ccell)

use DGV_commvar, only: nodes_phican,nodes_AkorshiftAllNets,AkorAllNets_capphi,numKornets
!!!!!!!!!!
integer (I4B), intent(in) :: ccell ! the number of the canonical cell -- a cell in the middle of the domain for which A is computed 
!!!!!!!!!!
integer (I4B) :: i,mm,AKorshift,phican  ! scrap variables...
integer (I4B) :: iphi ! the number of the basis function. Recall that each node generates a basis function.  
integer :: loc_alloc_stat ! scrap variable to keep the allocation status
!!!!!!!!!!!!!!!!!!!!
if (ccell < 0) then  ! a quick check for garbage in data: 
 print *,"SetAkorshiftAllNets_DGV: Error. The provided value of (ccell) is negative or invalid. Stop."
end if 
!!!!
mm = size(nodes_phican,1)
!!!!
do i=1,numKornets
 allocate (nodes_AkorshiftAllNets(i)%p(1:mm), stat=loc_alloc_stat)        !
    !
 if (loc_alloc_stat >0) then 
  print *, "SetAkorshiftAllNets_DGV: Allocation error for variable  nodes_AkorshiftAllNets(i)%p(1:mm), netnum=", i
  stop
 end if
 !! Now that the shift array is allocated, we need to populate it
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iphi=1,mm
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Now we need to calculate the nodes_AKorshift number 
    !! ATTENTION: when no Akor file present, the next line generates an access violation code! Program crashes
    phican = nodes_phican(iphi)
    AKorshift = sum(AkorAllNets_capphi(i)%p(1:phican-1))
    nodes_AkorshiftAllNets(i)%p(iphi) = AKorshift
    !!!!! all done
 end do 
end do 
end subroutine SetAkorshiftAllNets_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrouine InitDGV0D
!
! This is a conglomerate subrouine that 
! performes initialization of the 
! parameters, contants and the variables 
! of the DGV Library
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitDGV0D

use DGV_readwrite

use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,&
                  MvII,MuII,MwII,suII,svII,swII, &    
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,&
                  min_sens,Trad,fm,fmA,run_mode,Cco,&
				  Order_nu,Order,SecondaryVelGridInUse,&
				  MoRatesArry0DAllocated,MoRatesArry0D,&
				  NuNextUpdateTime0D,NuLastDens0D,NuLastTemp0D,numchnksII,&
                  A_capphiII,nodes_pcellII


use DGV_dgvtools_mod

!!!!!!!!!!!!!! Variables
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: i,canonical_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the next group of commands will initiate the key variables of the 
! DGV library. Change this commands if you need 
! different functionality of the library -- say if it is 
! desired to have a different initial velocity grid 
! or additional DGV variables need to be intitialized. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First we read parameters from the DGV library parameter file.
! the file name "DGVparametes.dat" is currently reserved for the DGV parameters. 
! However, this name can be changes
! by default, the file should be located in the same directory where the executable is started
!
call SetDGVParams("DGVparameters.dat",0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Mv=Mv_list(1)
Mu=Mu_list(1)
Mw=Mw_list(1)
su=su_list(1)
sv=sv_list(1)
sw=sw_list(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SetGDVGnodes    ! set some supplementary arrays to be used in building meshes
call SetDGVblzmmesh  ! sets one-dimensional grids in u,v, and w
! newt we build the velocity discretization 
call Set3DCellsR_DGV ! We create initial the velocity cells
call SetNodesDGV ! this will populate velocity cells with the 
! nodal - DG velocity nodes (see the description of the nodal-DG approach)

!!!!!!
! If you need to refine the velocity cells the following subroutines can be used: 
! in the next two lines the cell number 14 isdivived in 8 cubcells and in the next line
! the cell number 41 in the new mech is divided into 27 cubcells
!!! call CellsRefineDGV((/ 14 /),2,2,2) ! 
!!! call CellsRefineDGV((/ 27 + 14 /),3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!
! if one wants to record the cells, gridds and nodes information, use the following commands
! the parameters for the files names are taken from the DGV library parameter file
! the grids, cells, and nodes can be used in Matlab graphing subroutines. 
! they are usually not used in restart. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!
call WriteGridsDGV
call WriteCellsDGV
call WriteNodesDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines: 
! call WriteI1DGV(0.0_DP,0.0_DP,0.0_DP) ! this one will print the number of velocity nodes that 
                ! is contained in the cell that has the given velocity point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CREATING NEW OPERATORS A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a call that usually is commented. 
! the next subroutines are an OpenMP subroutines
! that are included in this library and will create a 
! new instance of Operator A or Akor.  Specifically, 
! Operators A and Akor depends on the collision model and 
! the DG discretization. Note that subrouitnes creating the operator 
! A and Akor use the primary mesh. (In contrast, when we use pre-computed Operators A and Akor, 
! we can use eithr primary or secondary meshes) Once the (primary) DG discretization is set, 
! on can compute the Operators A or Akor. To do that, uncomment the corresponding line below
! Be aware that pre-comupting A and Akor is a very slow process and 
! usually should not be attempted on one processor for 
! meshes exceeding 33^3 velocity notes. 
!  
! call SetA_DGV        ! call evaluation of A on Gauss nodes
! call SetAKorobov_DGV ! call evaluation of A on Korobov nodes
!!!!!!!!!!!!!!!!!
! For lagre meshes, one should use MPI versio of the 
! subroutine SetA_DGV . 
!
! TO BE ADDED LATER...  
!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Once the Operator A arrays are computed, 
! we need to save then on the hard drive in the 
! directory specified in the DGV parameter file
! To save Operator A on thehard Drive use the 
! the following subroutine 
!!!!!!!
! call WriteAarraysDGV      ! write A on Gauss nodes on hard drive 
! call WriteAKorArraysDGV   ! write A on Korobov nodes on hard drive 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! uncomment the stop directive here is only computing A or Akor arrays....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! stop 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END CREATING NEW OPERATORS A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ pre-computed operators A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This will read Operator A arrays using the name specified in the DGV parameter file
call ReadAarraysDGV 
!!!! This will read Operator A arrays is they are broken into a number of chunks
!! call ReadAarraysChnksDGV(25) ! use 9 if the A-array is divided into 9 chunks numbered 0-8. So, 
                                 ! this parameter is the number of chunks, not the number of the last 
                                 ! chunk. Example: to read chunks _ch000 to _ch008 use value 9 here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Sometimes, it is beneficial to truncate A based on |u-u_1| 
!!! Uncomment if desired to truncate based |u-u_1|. Any value of A that corresponds to 
!!! |u-u_1|> then the provided number is truncated
call TruncAarrsRadius_DGV (3.0_DP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! For FFB algorithm, it is benefitial to truncate the kernel 
!!!!! to reduce aliasing. Use the following subroutine to truncate based on distance from the origin.
!!!!! Comment: need to be replace with the distance from the center of the canonical cell.
!call TruncAarrsRadiusFromOrigin_DGV(2.0_DP) 

!!!! This will read a collection of Akor arrays that correspond to several Korobov nets.
!!!! Comment if not using Korobov Quadrature evaluation of the collision integral 
! call ReadAKorArraysAllNets_DGV 
!!!!!!! Random number will be used later in the code for Korobov integration 
call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END READ pre-computed operators A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Once the A and Akor arrays are in the memory, they are still can not be used untill a few additional arrays are 
! set up. The next subroutines set up those arrays.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UNCOMMENT IF USING EITHER A or AKOR ON THE PRIMARY MESH
! This subroutine will set arrays nodes_dui, _dvi, dwi, and nodes_phican
call SetCellsUGI_DGV(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
! UNCOMMENT IF USING A ON PRIMARY MESH
! This subroutine sets up arrays Ashift
call SetAshift_DGV(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
! UNCOMMENT IF USING Akor on primary mesh (integrating Boltzmann collision operator using Korobov quadratures)
! This subroutine sets up arrays AkorAllNets_shift
! call SetAkorshiftAllNets_DGV(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following few lines set up necessary variables to work with the Boltzmann collision integral. 
! We note that in this code, the Botlzmann collision integral is using secondary mesh. This means that
! the unknown solution will have to be mapped between the primary and secondary meshes. The secondary mesh, 
! in principle can be the same as primary -- then no interpolation is needed. However, in general, secondary mesh is coarser a
! and is use to evaluate the full Boltzmann Collision integral (it might be that the task is just too expensive on the primary mesh
!
! Because the collisio operator will use secondary mesh we check the flag SecondaryVelGridInUse whether 
! the secondary mesh is in use. If it is, we create the mesh and then may read the operator A
!
!!! DEBUG --- uncomment to run Pure Boltzmann on primary grid
SecondaryVelGridInUse=.false.
!!! END DEBUG
if (SecondaryVelGridInUse) then 
! We create variables of the secondary mesh
! A quick check of we have paparameters of the secondary mesh defined: 
 if ((size(Mu_list,1)<2) .or. (size(Mv_list,1)<2) .or.(size(Mw_list,1)<2) .or. &
           (size(su_list,1)<2) .or. (size(sv_list,1)<2) .or. (size(sw_list,1)<2) ) then 
  ! Error -- secondary mesh requested, but at least one parameter is missing for it 
  print *, "InitDGV0D: Secondary velocity meshhas been requested, but at least one parameter is missing for secondary mesh."  
  stop
  ! 
 end if
 ! set up the parameters of the secondary nodal-DG discretization. Secondary and primary discretizations use the same 
 ! boundaries in the velocity space
 MvII=Mv_list(2)
 MuII=Mu_list(2)
 MwII=Mw_list(2)
 suII=su_list(2)
 svII=sv_list(2)
 swII=sw_list(2)
 ! ready to create secondary meshes 
 ! NOTE: the suffix II in the name indicates that the subroutine works with the second mesh 
 call SetDGVblzmmeshII  ! sets one-dimensional grids in u,v, and w 
 ! newt we build the velocity discretization 
 call Set3DCellsR_DGVII ! We create initial velocity cells
 call SetNodesDGVII ! this will populate velocity cells with the 
 ! nodal - DG velocity nodes (see the description of the nodal-DG approach)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Initializing a PRE-COMPUTED Operator A
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Note that Pre-computed Operator A wil be used for evaluating collision integral and estimating relaxation rates. 
 ! These rates will be used in the velocity dependent BGK model. For this model evaluatin will be done using secondary grid. 
 ! Therefore one needs to indicate that secondary velocoity grid is used.
 ! For that the flag SecondaryVelGridInUse needs to be set true in the DGVparameters file.
 ! 
 ! It may happen that the secondary grid is the same as the primary, but in general we 
 ! will expect that the solution is interpolated between the grids
 !!!!!!!!!!!!!!! 

 ! The evaluation of the Boltzmann collision operator in the method of 
 !  Alekseenko-Josyula requires the 
 ! knowledge of pre-computed Operator A. 
 ! There has been a library of instances of 
 ! Operator A computed for different parameters of DG discretizations
 ! The pre-computed operator is loaded using the next batch of commands:
 !
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!!! This will read Operator A arrays using the name specified in the DGV parameter file
 call ReadAarraysDGVII 
 !!!! This will read Operator A arrays is they are broken into a number of chunks
 !! call ReadAarraysChnksDGVII(numchnksII) ! use 9 if the A-array is divided into 9 chunks numbered 0-8. So, 
                                 ! this parameter is the number of chunks, not the number of the last 
                                 ! chunk. Example: to read chunks _ch000 to _ch008 use value 9 here
 !!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Once the Operator A is loaded into the memory, it can be 
 ! truncated somwhat by using two means of truncation:
 ! Truncation using the entry treshhold and truncation to distance between the 
 ! velocity nodes. In the first method, all entried of operator A that
 ! fall below the truncation threshold by the absolute value (min_sens) are eliminated
 ! The following subroutine is sued to cut values based on the minimum value:
 !!!!!
      !call TruncAarrsTrhls_DGVII ( 200.0_DP ) !DGVblzmPl.f90
 !!!!!!
 ! Another way to truncate is to eliminate all values of operator 
 ! A(v_1,v_2,\phi) that correspond to v_1 and v_2 separated by more than 
 ! certain distance. Specifically if \| v_1 - v_2 \|<Trad then the value of 
 ! A is eliminated from the arrays
 !
 ! call TruncAarrsRadiusFromOrigin_DGVII ( 3.0_DP ) 
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! WARNING: the use of both truncation subroutines has to be careful -- 
 ! it is easy to "cut" too much -- usual symptoms  of "cutting too much" 
 ! are too slow relaxation while maintaining correct conservation laws. 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !
 ! Now we need to prepare some helpful arrays -- these will be used in the evaluation 
 ! of the Boltzmann collision operator, evaluation of the linearized Botlzmann,  and evaluation of the decomposed Botlzmann operator
 do i=1,size(A_capphiII)
   if (A_capphiII(i) > 0) then
     canonical_cell=nodes_pcellII(i)
     exit
   end if
 end do
   
 call SetCellsUGINdsSftArrsII(canonical_cell) ! miscset.f90
 !
 ! Another preparatory call. This one sets up array (nodes_primcellII) used for projecting the primary solution on the secondary mesh
 call prepare_sDGVpMap_DGV ! this subroutine nodes_primcellII
 
 ! the next variable keeps linearized Botlzmann operator in it  -- 
 ! fmAis uses a lot of memory (su*Mu*sv*Mv*sw*Mw)^2*doublereal
 ! unless linearized Boltzmann is used on a single node, do not uncomment it
 !!!!!!!!!!     
 !! allocate (fmA(1:size(nodes_u,1),1:size(nodes_u,1)), stat=loc_alloc_stat)
 !! ! check if the allocation was successful    
 !! if (loc_alloc_stat >0) then 
 !!  print *, "InitDGV0D: Allocation error for variables (fmA)"
 !! end if 
 !!!!!!!!!!! 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
! End of the section dedicated to secondary mesh and end of preparatory work for
! Using full Botlzmann collision operator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
allocate (fm(1:size(nodes_u,1)), stat=loc_alloc_stat) ! local Maxwellian and array to keep the velocity dependent collision frequency
! check if the allocation was successful    
if (loc_alloc_stat >0) then 
 print *, "InitDGV0D: Allocation error for variables (fm)"
end if 
!
!!!! DIAGNOSTIC for Prakash can remove later... !!!! 
! allocate the arrays used for preserving the coefficients of the velocity dependent collision frequency
allocate (Cco(Order_nu), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV0D: Allocation error for (Cco)"
end if
Cco=0 ! reset the coefficients to zero
!!!!!!!!!!!!!!!!!!!!!!! END Diagnostic for prakash. Can remove later... 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we check if the storage for the moments rates that are used in the velocity dependent collision frequency has been created yet
if (.not. MoRatesArry0DAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (MoRatesArry0D(1:Order), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "InitDGV0D: Error allocation arrays to store coefficient of the VDCF (MoRatesArry0D)"
       stop
      end if 
 MoRatesArry0DAllocated = .true.
! Also, we will set up the other arrays to make sure that the coefficients are updated:
MoRatesArry0D = 0 ! reset the array
NuNextUpdateTime0D = - 100000.0   ! these bogus values will set off the criteria for update
NuLastDens0D = -100000.0 
NuLastTemp0D = -100000.0 
! 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!    
!! INITIALIZE THE RUN_MODE to ZERO -- in the beginning the regime is strongly non-linear
  run_mode=0
!!!!!!!!!!!!!!!!!!!!!!!

end subroutine InitDGV0D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrouine InitDGV1D(nspatialcells)
!
! This is copy of the above subroutine, except it creates some arrays to 
! support multiple spatial cells. 
!
! ATTENTION: you will need to knwo the total number of the spatial cells that the discretization will have 
! before you call this subroutine.
!
! This is a conglomerate subrouine that 
! performes initialization of the 
! parameters, contants and the variables 
! of the DGV Library
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitDGV1D(nspatcells)

use DGV_readwrite

use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,&
                  MvII,MuII,MwII,suII,svII,swII, &    
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,&
                  min_sens,Trad,fm,fmA,run_mode,Cco1D,&
				  Order_nu,Order,SecondaryVelGridInUse,&
				  MoRatesArry1DAllocated,MoRatesArry1D,&
				  NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D,Nu1DUpdateArrysAllocated,&
				  TotNumSpatCells,MoRatesReliable1D,run_mode_1D


use DGV_dgvtools_mod


!!!!!!!!!!!!!! Variables
integer (I4B), intent (in) :: nspatcells ! the total number of the spatial cells --- arrays will be created to store information for each spatial cell. 
                                         ! these arrays support velocity depemdent collision frequency model 

integer (I4B) :: loc_alloc_stat ! to keep allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the next group of commands will initiate the key variables of the 
! DGV library. Change this commands if you need 
! different functionality of the library -- say if it is 
! desired to have a different initial velocity grid 
! or additional DGV variables need to be intitialized. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First we read parameters from the DGV library parameter file.
! the file name "DGVparametes.dat" is currently reserved for the DGV parameters. 
! However, this name can be changes
! by default, the file should be located in the same directory where the executable is started
!
call SetDGVParams("DGVparameters.dat",0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Mv=Mv_list(1)
Mu=Mu_list(1)
Mw=Mw_list(1)
su=su_list(1)
sv=sv_list(1)
sw=sw_list(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SetGDVGnodes    ! set some supplementary arrays to be used in building meshes
call SetDGVblzmmesh  ! sets one-dimensional grids in u,v, and w
! newt we build the velocity discretization 
call Set3DCellsR_DGV ! We create initial the velocity cells
call SetNodesDGV ! this will populate velocity cells with the 
! nodal - DG velocity nodes (see the description of the nodal-DG approach)

!!!!!!
! If you need to refine the velocity cells the following subroutines can be used: 
! in the next two lines the cell number 14 isdivived in 8 cubcells and in the next line
! the cell number 41 in the new mech is divided into 27 cubcells
!!! call CellsRefineDGV((/ 14 /),2,2,2) ! 
!!! call CellsRefineDGV((/ 27 + 14 /),3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!
! if one wants to record the cells, gridds and nodes information, use the following commands
! the parameters for the files names are taken from the DGV library parameter file
! the grids, cells, and nodes can be used in Matlab graphing subroutines. 
! they are usually not used in restart. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!
call WriteGridsDGV
call WriteCellsDGV
call WriteNodesDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines: 
! call WriteI1DGV(0.0_DP,0.0_DP,0.0_DP) ! this one will print the number of velocity nodes that 
                ! is contained in the cell that has the given velocity point
! call InitWriteDet    ! velocity dependent collision frequency. These subroutines only create an empty file. 
                ! later calls of related subroutines will save some diagnostic data in the created 
                ! files. These files make problems in many clusters becuase of the 
                ! writing protection. If the file is not deleted, it can not be overwritten by 
                ! a new instance of the running software -- we commented the use of the subroutines for now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CREATING NEW OPERATORS A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a call that usually is commented. 
! the next subroutines are an OpenMP subroutines
! that are included in this library and will create a 
! new instance of Operator A or Akor.  Specifically, 
! Operators A and Akor depends on the collision model and 
! the DG discretization. Note that subrouitnes creating the operator 
! A and Akor use the primary mesh. (In contrast, when we use pre-computed Operators A and Akor, 
! we can use eithr primary or secondary meshes) Once the (primary) DG discretization is set, 
! on can compute the Operators A or Akor. To do that, uncomment the corresponding line below
! Be aware that pre-comupting A and Akor is a very slow process and 
! usually should not be attempted on one processor for 
! meshes exceeding 33^3 velocity notes. 
!  
! call SetA_DGV ! dgvtools_mod.f90
! call SetAKorobov_DGV
!!!!!!!!!!!!!!!!!
! For lagre meshes, one should use MPI versio of the 
! subroutine SetA_DGV . 
!
! TO BE ADDED LATER...  
!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Once the Operator A arrays are computed, 
! we need to save then on the hard drive in the 
! directory specified in the DGV parameter file
! To save Operator A on thehard Drive use the 
! the following subroutine 
!!!!!!!
! call WriteAarraysDGV
! call WriteAKorArraysDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ pre-computed operators A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This will read Operator A arrays using the name specified in the DGV parameter file
!!  call ReadAarraysDGV 
!!!! This will read Operator A arrays is they are broken into a number of chunks
!! call ReadAarraysChnksDGV(25) ! use 9 if the A-array is divided into 9 chunks numbered 0-8. So, 
                                 ! this parameter is the number of chunks, not the number of the last 
                                 ! chunk. Example: to read chunks _ch000 to _ch008 use value 9 here
!!!!!!!!!!!!!!!!!!!!!!!!
!!!! This will read a collection of Akor arrays that correspond to several Korobov nets.
!!!! Comment if not using Korobov Quadrature evaluation of the collision integral 
 !call ReadAKorArraysAllNets_DGV 
!!!!!!! Random number will be used later in the code for Korobov integration 
call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END READ pre-computed operators A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Once the A and Akor arrays are in the memory, they are still can not be used untill a few additional arrays are 
! set up. The next subroutines set up those arrays.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UNCOMMENT IF USING EITHER A or AKOR ON THE PRIMARY MESH
! This subroutine will set arrays nodes_dui, _dvi, dwi, and nodes_phican
call SetCellsUGI_DGV(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
! UNCOMMENT IF USING A ON PRIMARY MESH
! This subroutine sets up arrays Ashift
! call SetAshift_DGV(FindCellContainsPoint_DGVII(0.0_DP,0.0_DP,0.0_DP))
! UNCOMMENT IF USING Akor on primary mesh (integrating Boltzmann collision operator using Korobov quadratures)
! This subroutine sets up arrays AkorAllNets_shift
 !call SetAkorshiftAllNets_DGV(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following few lines set up necessary variables to work with the Boltzmann collision integral. 
! We note that in this code, the Botlzmann collision integral is using secondary mesh. This means that
! the unknown solution will have to be mapped between the primary and secondary meshes. The secondary mesh, 
! in principle can be the same as primary -- then no interpolation is needed. However, in general, secondary mesh is coarser a
! and is use to evaluate the full Boltzmann Collision integral (it might be that the task is just too expensive on the primary mesh
!
! Because the collisio operator will use secondary mesh we check the flag SecondaryVelGridInUse whether 
! the secondary mesh is in use. If it is, we create the mesh and then may read the operator A
!
if (SecondaryVelGridInUse) then 
! We create variables of the secondary mesh
! A quick check of we have paparameters of the secondary mesh defined: 
 if ((size(Mu_list,1)<2) .or. (size(Mv_list,1)<2) .or.(size(Mw_list,1)<2) .or. &
           (size(su_list,1)<2) .or. (size(sv_list,1)<2) .or. (size(sw_list,1)<2) ) then 
  ! Error -- secondary mesh requested, but at least one parameter is missing for it 
  print *, "InitDGV0D: Secondary velocity meshhas been requested, but at least one parameter is missing for secondary mesh."  
  stop
  ! 
 end if
 ! set up the parameters of the secondary nodal-DG discretization. Secondary and primary discretizations use the same 
 ! boundaries in the velocity space
 MvII=Mv_list(2)
 MuII=Mu_list(2)
 MwII=Mw_list(2)
 suII=su_list(2)
 svII=sv_list(2)
 swII=sw_list(2)
 ! ready to create secondary meshes 
 ! NOTE: the suffix II in the name indicates that the subroutine works with the second mesh 
 call SetDGVblzmmeshII  ! sets one-dimensional grids in u,v, and w 
 ! newt we build the velocity discretization 
 call Set3DCellsR_DGVII ! We create initial velocity cells
 call SetNodesDGVII ! this will populate velocity cells with the 
 ! nodal - DG velocity nodes (see the description of the nodal-DG approach)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Initializing a PRE-COMPUTED Operator A
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Note that Pre-computed Operator A wil be used for evaluating collision integral and estimating relaxation rates. 
 ! These rates will be used in the velocity dependent BGK model. For this model evaluatin will be done using secondary grid. 
 ! Therefore one needs to indicate that secondary velocoity grid is used.
 ! For that the flag SecondaryVelGridInUse needs to be set true in the DGVparameters file.
 ! 
 ! It may happen that the secondary grid is the same as the primary, but in general we 
 ! will expect that the solution is interpolated between the grids
 !!!!!!!!!!!!!!! 

 ! The evaluation of the Boltzmann collision operator in the method of 
 !  Alekseenko-Josyula requires the 
 ! knowledge of pre-computed Operator A. 
 ! There has been a library of instances of 
 ! Operator A computed for different parameters of DG discretizations
 ! The pre-computed operator is loaded using the next batch of commands:
 !
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!!! This will read Operator A arrays using the name specified in the DGV parameter file
 call ReadAarraysDGVII 
 !!!! This will read Operator A arrays is they are broken into a number of chunks
 !! call ReadAarraysChnksDGVII(25) ! use 9 if nine chunks are in A-arrays, so it is always the number of chunks
 !!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Once the Operator A is loaded into the memory, it can be 
 ! truncated somwhat by using two means of truncation:
 ! Truncation using the entry treshhold and truncation to distance between the 
 ! velocity nodes. In the first method, all entried of operator A that
 ! fall below the truncation threshold by the absolute value (min_sens) are eliminated
 ! The following subroutine is sued to cut values based on the minimum value:
 !!!!!
 !     call TruncAarrsTrhls_DGVII ( min_sens ) !DGVblzmPl.f90
 !!!!!!
 ! Another way to truncate is to eliminate all values of operator 
 ! A(v_1,v_2,\phi) that correspond to v_1 and v_2 separated by more than 
 ! certain distance. Specifically if \| v_1 - v_2 \|<Trad then the value of 
 ! A is eliminated from the arrays
 !
 !!     call TruncAarrsRadius_DGVII ( Trad ) 
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! WARNING: the use of both truncation subroutines has to be careful -- 
 ! it is easy to "cut" too much -- usual symptoms  of "cutting too much" 
 ! are too slow relaxation while maintaining correct conservation laws. 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 !
 ! Now we need to prepare some helpful arrays -- these will be used in the evaluation 
 ! of the Boltzmann collision operator, evaluation of the linearized Botlzmann,  and evaluation of the decomposed Botlzmann operator
 call SetCellsUGINdsSftArrsII(FindCellContainsPoint_DGVII(0.0_DP,0.0_DP,0.0_DP)) ! miscset.f90
 !
 ! Another preparatory call. This one sets up array (nodes_primcellII) used for projecting the primary solution on the secondary mesh
 call prepare_sDGVpMap_DGV ! this subroutine nodes_primcellII
 
 ! the next variable keeps linearized Botlzmann operator in it  -- 
 ! fmAis uses a lot of memory (su*Mu*sv*Mv*sw*Mw)^2*doublereal
 ! unless linearized Boltzmann is used on a single node, do not uncomment it
 !!!!!!!!!!     
 !! allocate (fmA(1:size(nodes_u,1),1:size(nodes_u,1)), stat=loc_alloc_stat)
 !! ! check if the allocation was successful    
 !! if (loc_alloc_stat >0) then 
 !!  print *, "InitDGV0D: Allocation error for variables (fmA)"
 !! end if 
 !!!!!!!!!!! 

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
! End of the section dedicated to secondary mesh and end of preparatory work for
! Using full Botlzmann collision operator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
allocate (fm(1:size(nodes_u,1)), stat=loc_alloc_stat) ! local Maxwellian and array to keep the velocity dependent collision frequency
! check if the allocation was successful    
if (loc_alloc_stat >0) then 
 print *, "InitDGV1D: Allocation error for variables (fm)"
end if 

TotNumSpatCells=nspatcells ! record the total number of cells for future checks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Atention: this array may be used in the BGK model with velocity dependent collision frequency. 
!! However, currently, the version of the subrouine does not make use if this array: therefore its allocation is commented. 
!
! allocate the arrays used for preserving the coeeficients of the velocity dependent collision frequency
!allocate (Cco1D(Order_nu,nspatcells), stat=loc_alloc_stat)
!if (loc_alloc_stat >0) then
! print *, "InitDGV1D: Allocation error for (Cco1D)"
!end if
!Cco1D=0 ! reset the coefficients to zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! next we check if the storage for the moments rates that are used in the velocity dependent collision frequency has been created yet
if (.not. MoRatesArry1DAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (MoRatesArry1D(1:Order,nspatcells), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "InitDGV1D: Error allocation arrays to store coefficient of the VDCF (MoRatesArry0D)"
       stop
      end if 
 MoRatesArry1DAllocated = .true.
! Also, we will set up the other arrays to make sure that the coefficients are updated:
MoRatesArry1D = 0 ! reset the array
end if 
! Allocate arrays that will keep the flag is the relaxation rates are reliable in the spatial cell: 
allocate (MoRatesReliable1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D: Allocation error for (MoRatesReliable1D)"
end if
MoRatesReliable1D=.false. ! reset the coefficients to false
!
! Allocate arrays that will keep the run_mode for ech spatial cell: 
allocate (run_mode_1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D: Allocation error for (run_mode_1D)"
end if
run_mode_1D = 0 ! reset the coefficients to mode=0 = full Boltzmann!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we check if the storage has been created for information about the last update of the relaxation rates
if (.not. Nu1DUpdateArrysAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (NuNextUpdateTime1D(nspatcells),NuLastDens1D(nspatcells),NuLastTemp1D(nspatcells), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "InitDGV1D: Error allocation arrays to store information about last rates apdate: (NuNextUpdateTime1D)," //  &
               "(NuLastDens1D),(NuLastTemp1D)" 
       stop
      end if 
Nu1DUpdateArrysAllocated = .true.
! Also, we will initialize the arrays to make sure that the coefficients are updated:
NuNextUpdateTime1D = - 100000.0   ! these bogus values will set off the criteria for update
NuLastDens1D = -100000.0 
NuLastTemp1D = -100000.0 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!    THIS VARIABLE SHOULD NOT BE USED --- CHECK AND REMOVE LATER>
!! INITIALIZE THE RUN_MODE to ZERO -- in the beginning the regime is strongly non-linear
  run_mode=0
!!!!!!!!!!!!!!!!!!!!!!!

end subroutine InitDGV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! prepare_sDGVpMap_DGV
!
! This subrouine sets some supplementary arrays that will be used to project the solution from the primary to the secondary meshes
! 
! In particular, this subroutine sets the array nodesII_primcell -- this array identifies the cell on the primary node that  This subroutine sets an inmatrix for converting the solutions from the primary velocity mesh to the secondary velocity mesh 
! In the case when the both meshes are uniform. 
!
! ATTENTION! do not call before both primary and secondary meshes are set up! 
subroutine prepare_sDGVpMap_DGV

use DGV_commvar, only: nodes_uII,nodes_vII,nodes_wII,nodes_primcellII

use DGV_dgvtools_mod

!!!!!
integer (I4B):: i ! 
integer  :: loc_alloc_stat ! to keep allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to allocate the nodes_primcellII array
allocate (nodes_primcellII(size(nodes_uII,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "prepare_sDGVpMap_DGV: Allocation error for (nodes_primcellII)"
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we will now go over the secondary mesh and find out where the secondary nodes falls on the primary mesh.
! the number of the cell is recoded into the array nodes_primcellII
do i=1,size(nodes_uII,1)
  nodes_primcellII(i) = FindCellContainsPoint_DGV(nodes_uII(i),nodes_vII(i),nodes_wII(i))
end do 
! all done 
end subroutine prepare_sDGVpMap_DGV

end module DGV_miscset
