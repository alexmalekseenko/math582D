!
! DGV_dgv_mod.f90
! 
! This module contains subroutines for working with Discontinuous Galerkin velocity discretization. 
! It contains definitions and routines for working with the basis functions 
! 
! Subroutines of this module complement the subroutines of setting up grids and cells and nodes in miscset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!include 'lapack.f90'
    
include 'mkl_dfti.f90' 

module DGV_dgvtools_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none


 interface lagrbasfun
     module procedure lagrbasfun, lagrbasfun_xvec
   end interface

real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0


contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! lagrbasfun
!
! This subroutine evaluates the lagrange basis function given the x and the 
! array of nodes kappa
! 
! \Prod_{j\neq i} \frac{kappa(j)-x}{kappa(j)-kappa(i)}
!
!!!!!!!!!!!!!!!!!!!!!!!

function lagrbasfun(i,x,kappa) result (y)

integer (I4B) :: i ! the number of the node where the lagrange basis function is one. 
real (DP) :: x ! the point where the function needs to be evaluated
real (DP), dimension(:) :: kappa ! the nodes of the lagrange basis functions
real (DP) :: y ! the value of the function 
!
integer (I4B) :: j ! some local counter


!!!!!!!!!!!!!1
! a quick consistancy check
if ((i<1) .or. (i>size(kappa,1))) then 
 print *," lagrbasfun: error. (i<1) .or. (i>size(kappa,1)) no Lagrange basis function with this number."
 stop
endif  
!
y=1.0_DP
do j=1,i-1
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 
do j=i+1,size(kappa,1)
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 

end function lagrbasfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lagrbasfun
!
! This subroutine evaluates the lagrange basis function given the x and the 
! array of nodes kappa
! 
! \Prod_{j\neq i} \frac{kappa(j)-x}{kappa(j)-kappa(i)}
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function lagrbasfun_xvec(i,x,kappa) result (y)

integer (I4B) :: i ! the number of the node where the lagrange basis function is one. 
real (DP), dimension(:) :: x ! the point where the function needs to be evaluated
real (DP), dimension(:) :: kappa ! the nodes of the lagrange basis functions
real (DP), dimension(1:size(x,1)) :: y ! the value of the function 
!
integer (I4B) :: j ! some local counter


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick consistancy check
if ((i<1) .or. (i>size(kappa,1))) then 
 print *," lagrbasfun: error. (i<1) .or. (i>size(kappa,1)) no Lagrange basis function with this number."
 stop
endif  
!
y=1.0_DP
do j=1,i-1
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 
do j=i+1,size(kappa,1)
y = y*(kappa(j)-x)/(kappa(j)-kappa(i))
enddo 

end function lagrbasfun_xvec
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EvalLagrBasisFunByNdsDGblzm
! 
! This subroutine evaluates the basis function identified by a node 
! for a given a value of velocity. To do this, first the cell is found 
! to which this velocity belongs. Then the cell is found to which the basis function belongs and 
! the numbers of the 1D-basis functions are identified. Then if the basis function is defined on this 
! cell, then the basis function is evaluated using the function lagrbasfun. If the velocity and the 
! basis function belong to different cells then the value of the basis function is zero
!
! This function directly uses arrays from commvar.f90. Specifically, grids, cells and nodes.
!
! u,v,w = components of the velocity where the basis function needs to be evaluated. 
! 
! i = the number of the node associated with this basis function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalLagrBasisFunByNdsDGblzm(u,v,w,i) result (y)

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity
integer (I4B) :: i  ! the number of the node in the nodes array. this the number of the node associated with 
                    ! the given basis function. Thus basis functions may be numbered using nodes... 
real (DP) ::  y ! the value of the basis function on the given velocity                    

!!!!!!!!!!!!!!                      
integer (I4B) :: j ! a counter -- usually denotes the number of the grid we are working on. 
integer (I4B) :: gi_v,gi_u,gi_w ! counters to help count nodes in the grids array
integer (I4B) :: celli, cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
integer (I4B) :: ui,vi,wi ! are the local numbers on the particular grids where u belongs. Because we have ierarchicval grids, 
             ! there may be different grids that have ui,vi,wi defined. If one of these numbers is zero -- then the 
             ! velocity u,v,w is not on the grid. Unless a mistake was done in setting the grids either all three of them 
             ! are zero -- that means that the velocity is not on this (level zero grid) or all three are non-zero  -- means that
             ! the velocity has been found. Mized results may indicate a breach in data integrity. 
             !!!!!
integer (I4B) :: jj ! a counter
real (DP) :: unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
! we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
! if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
!!!!!!!!
!
j=1
!! set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0
gi_v=0
gi_w=0
do while (j<=size(grids_cap_u,1))
  call finduiviwi(grids_u(gi_u+1:gi_u+grids_cap_u(j)),grids_v(gi_v+1:gi_v+grids_cap_v(j)),&
                       grids_w(gi_w+1:gi_w+grids_cap_w(j)),u,v,w,ui,vi,wi)
  if (ui*vi*wi /= 0) then 
     ! now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1
     do jj=1,j-1
      celli = celli + (grids_cap_u(jj)-1)*(grids_cap_v(jj)-1)*(grids_cap_w(jj)-1)
     end do
     celli=celli + (ui-2)*(grids_cap_v(j)-1)*(grids_cap_w(j)-1)+(vi-2)*(grids_cap_w(j)-1)+(wi-2)
     !  check if this cell is refined
     if ((cells_refu(celli)>1) .or. (cells_refv(celli)>1) .or. (cells_refw(celli)>1)) then  
        j=cells_cgrid(celli)
        ! set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0
        gi_v=0
        gi_w=0
        do jj=1,j-1
         gi_u = gi_u + grids_cap_u(jj)
         gi_v = gi_v + grids_cap_v(jj)
         gi_w = gi_w + grids_cap_w(jj)
        end do
     else 
        exit
     endif
  else 
     celli=0
     exit
  endif                      
enddo
! now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)
! next step is to find the numbers that will help us compute the value of the basis function
! first we need to know to what cell in u,v,w this basis function belongs. This is simple becasue this 
! information is stored in the nodes arrays -- nodes_pcell
cellp=nodes_pcell(i)
! A quick check, if the velocity is not on the cell than the value of the basis function is zero and we are all done
if ((cellp /= celli) .or. (celli == 0)) then 
y=0.0_DP
                    else 
! if the velocity is on the cell, then we need to compute the value of $y$:
y=1.0_DP    
! next we need to know the three local indices that tell what velocity nodal values correspond to this 
! basis function. this is also simple since this information is also stored in the Nodes Arrays.
unor = ( u - (cells_ru(cellp) + cells_lu(cellp))/2.0_DP )/(cells_ru(cellp) - cells_lu(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_ui(i),unor,g_nds_all(:cells_gou(cellp),cells_gou(cellp)))
!
vnor = ( v - (cells_rv(cellp) + cells_lv(cellp))/2.0_DP )/(cells_rv(cellp) - cells_lv(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_vi(i),vnor,g_nds_all(:cells_gov(cellp),cells_gov(cellp)))
!
wnor = ( w - (cells_rw(cellp) + cells_lw(cellp))/2.0_DP )/(cells_rw(cellp) - cells_lw(cellp))*2.0_DP 
y=y*lagrbasfun(nodes_wi(i),wnor,g_nds_all(:cells_gow(cellp),cells_gow(cellp)))
!   
!!!!!!!!!!!!
end if
 
end function EvalLagrBasisFunByNdsDGblzm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EvalLagrBasisFunByNdsDGblzmEZ
! 
! This subroutine evaluates the basis function identified by a node 
! for a given a value of velocity. To do this, first the cell is found 
! to which this velocity belongs. We check if the velocity falls within the 
! cell where the basis function is defined. The 
! numbers of the 1D-basis functions are identified. Then the basis function is 
! evaluated using the function lagrbasfun. If the velocity and the 
! basis function belong to different cells then the value of the basis function is zero
!
! This function directly uses arrays from commvar.f90. Specifically, grids, cells and nodes.
!
! u,v,w = components of the velocity where the basis function needs to be evaluated. 
! 
! i = the number of the node associated with this basis function 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalLagrBasisFunByNdsDGblzmEZ(u,v,w,i) result (y)

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity
integer (I4B) :: i  ! the number of the node in the nodes array. this the number of the node associated with 
                    ! the given basis function. Thus basis functions may be numbered using nodes... 
real (DP) ::  y ! the value of the basis function on the given velocity                    

!!!!!!!!!!!!!!                      
integer (I4B) :: cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
real (DP) :: lu,ru,lv,rv,lw,rw,unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given the number of the node/basis function, we know the cell where the function belongs:
cellp= nodes_pcell(i) ! this is the number of the nodes parent cell. 
! Now we check is the velocity is within the cell bounds: 
lu=cells_lu(cellp); ru=cells_ru(cellp) 
lv=cells_lv(cellp); rv=cells_rv(cellp) 
lw=cells_lw(cellp); rw=cells_rw(cellp) 
if ((u<lu) .or. (ru<u) .or. (v<lv) .or. (rv<v) .or. (w<lw) .or. (rw<w)) then 
 y=0.0_DP
else 
 ! if we are here then the velocity belongs to the cell where the basis function is defined:
 ! next step is to find the numbers that will help us compute the value of the basis function
 ! first we need to know to what cell in u,v,w this basis function belongs. This is simple becasue this 
 ! information is stored in the nodes arrays -- nodes_pcell (see above)
 y=1.0_DP 
 ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
 ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
 unor = ( u - (ru+lu)/2.0_DP )/(ru - lu)*2.0_DP 
 y=y*lagrbasfun(nodes_ui(i),unor,g_nds_all(:cells_gou(cellp),cells_gou(cellp)))
 !
 vnor = ( v - (rv+lv)/2.0_DP )/(rv - lv)*2.0_DP 
 y=y*lagrbasfun(nodes_vi(i),vnor,g_nds_all(:cells_gov(cellp),cells_gov(cellp)))
 !
 wnor = ( w - (rw+lw)/2.0_DP )/(rw - lw)*2.0_DP 
 y=y*lagrbasfun(nodes_wi(i),wnor,g_nds_all(:cells_gow(cellp),cells_gow(cellp)))
 !!!!!!!!!!!!
end if
 
end function EvalLagrBasisFunByNdsDGblzmEZ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalLagrBasisFunByNdsDGblzmEZAccyCheck(u,v,w,i) result (y)
! This subroutine can be used instead of the above subroutine to estimate the conservative properties of the computed operator A
! use this subroutine intead of the above in SetA_DGV or SetAKorobov_DGV to obtain indicators of 
! integration errors. Analytically, the value of computed operator should be zero for any of the five functions below 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalLagrBasisFunByNdsDGblzmEZAccyCheck(u,v,w,i) result (y)

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity
integer (I4B) :: i  ! the number of the node in the nodes array. this the number of the node associated with 
                    ! the given basis function. Thus basis functions may be numbered using nodes... 
real (DP) ::  y ! the value of the basis function on the given velocity                    

!!!!!!!!!!!!!!                      
integer (I4B) :: cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
real (DP) :: lu,ru,lv,rv,lw,rw,unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! instead of evaluating the true basis function, we will evaluate the collision invariants. 
! the analitical value of the bilinear form is zero. However, numerical value will not be -- 
! this is a way to estimate the accuracy of the integral 

!!! Uncomment this line if want the zeros's moment 
!y=1.0_DP 
!!! Uncomment one of the three lines for the first moment:
!y = u
!y = v
!y = w 
!!! Uncomment the next line to get the second moment
y=u*u + v*v + w*w  

! Using this one of these substitutes for the value of the basis function will shows how conservative
! is the numerical descritization of the bilinear collision integral 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end function EvalLagrBasisFunByNdsDGblzmEZAccyCheck

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalSolVeloPtCell_DGV 
!
! This function will evaluate the solution at a velocity point on a specified velocity cell.
! The value of the velocity node is assumed to belong to the cell. If it is not, large interpolation 
! errors are expected. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalSolVeloPtCell_DGV(f,u,v,w,pcn) result (y)

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi,&
				   cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, & 
                   cells_gow, cells_gou, cells_gov, g_nds_all

real (DP), dimension (:), intent (in) :: f ! the solution 
real (DP), intent (in)     :: u,v,w    ! the components of the given velocity
integer (I4B), intent (in) :: pcn  ! the number of the cell to which this velocity belong
real (DP) ::  y       ! the value of the solution at the given point. 

!!!!!!!!!!!!!!                      
real (DP) :: unor,vnor,wnor,yy  ! scrap variables to keep the normalized velocity 
integer (I4B) :: j,gou,gov,gow
real (DP), dimension (:), allocatable  :: g_nds_u, g_nds_v, g_nds_w
integer :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!!!
unor = (u - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
vnor = (v - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
wnor = (w - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
! next we will create a temp array that will store the values of the Gauss nodes needed for interpolation on the cell (pcn)
! we record the orders of the gauss nodes in each velocity component for the cell #(pcn)
gou=cells_gou(pcn)
gov=cells_gov(pcn)
gow=cells_gow(pcn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  then we need to allocate arrays to keep the nodes: 
allocate (g_nds_u(gou),g_nds_v(gov),g_nds_w(gow), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "EvalSolVeloPtCell_DGV: Allocation error for (g_nds_u),(g_nds_v),(g_nds_w)"
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!
g_nds_u=g_nds_all(:gou,gou)
g_nds_v=g_nds_all(:gov,gov)
g_nds_w=g_nds_all(:gow,gow)
! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
! basis function for that node and add it to the interpolated value
y=0;
do j=1,size(f,1)
 if (nodes_pcell(j) == pcn) then 
  ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
  ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
  yy=1.0_DP ! reset the value of the basis function
  yy=yy*lagrbasfun(nodes_ui(j),unor,g_nds_u)
  yy=yy*lagrbasfun(nodes_vi(j),vnor,g_nds_v)
  yy=yy*lagrbasfun(nodes_wi(j),wnor,g_nds_w)
  ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
  y = y + f(j)*yy   
 endif
enddo 
deallocate(g_nds_u,g_nds_v,g_nds_w)

end function EvalSolVeloPtCell_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalSolVeloPtCellFast_DGV 
!
! This function is analogous to the previous one except it uses the unifrom structure fo the mech and 
! numbering conventions to speedup the evaluation.
!
! This function will evaluate the solution at a velocity point on a specified velocity cell.
! The value of the velocity node is assumed to belong to the cell. If it is not, large interpolation 
! errors are expected. 
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalSolVeloPtCellFast_DGV(f,u,v,w,pcn) result (y)

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi,&
				   cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, & 
                   g_nds_all,&
                   su,sv,sw

real (DP), dimension (:), intent (in) :: f ! the solution 
real (DP), intent (in)     :: u,v,w    ! the components of the given velocity
integer (I4B), intent (in) :: pcn  ! the number of the cell to which this velocity belong
real (DP) ::  y       ! the value of the solution at the given point. 

!!!!!!!!!!!!!!                      
real (DP) :: unor,vnor,wnor,yy  ! scrap variables to keep the normalized velocity 
integer (I4B) :: j
integer :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!!!
unor = (u - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
vnor = (v - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
wnor = (w - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
! basis function for that node and add it to the interpolated value
y=0;
do j=(pcn-1)*su*sv*sw+1,pcn*su*sv*sw
  ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
  ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
  yy=1.0_DP ! reset the value of the basis function
  yy=yy*lagrbasfun(nodes_ui(j),unor,g_nds_all(:su,su))
  yy=yy*lagrbasfun(nodes_vi(j),vnor,g_nds_all(:su,su))
  yy=yy*lagrbasfun(nodes_wi(j),wnor,g_nds_all(:su,su))
  ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
  y = y + f(j)*yy   
enddo 

end function EvalSolVeloPtCellFast_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalSolVeloPtCellFast_DGV 
!
! This function is analogous to the previous one except it uses the unifrom structure fo the mech and 
! numbering conventions to speedup the evaluation.
!
! This function will evaluate the solution at a velocity point on a specified velocity cell.
! The value of the velocity node is assumed to belong to the cell. If it is not, large interpolation 
! errors are expected. 
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function EvalSolVeloPtCellFastAlt_DGV(f,u,v,w,pcn) result (y)

use DGV_commvar, only: nodes_pcell, &
				   cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, & 
                   g_nds_all,&
                   su,sv,sw

real (DP), dimension (:), intent (in) :: f ! the solution 
real (DP), intent (in)     :: u,v,w    ! the components of the given velocity
integer (I4B), intent (in) :: pcn  ! the number of the cell to which this velocity belong
real (DP) ::  y       ! the value of the solution at the given point. 

!!!!!!!!!!!!!!                      
real (DP) :: unor,vnor,wnor, piunor,pjvnor,plwnor  ! scrap variables to keep the normalized velocity 
integer (I4B) :: i,j,l,node_ijl
integer :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!!!
unor = (u - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
vnor = (v - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
wnor = (w - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
! basis function for that node and add it to the interpolated value
y=0;
node_ijl = (pcn-1)*su*sv*sw
do i=1,su 
 piunor = lagrbasfun(i,unor,g_nds_all(:su,su))
 do j=1,sv
  pjvnor = lagrbasfun(j,vnor,g_nds_all(:sv,sv))
  do l=1,sw 
   ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
   ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
   plwnor = lagrbasfun(l,wnor,g_nds_all(:sw,sw))
  ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
  node_ijl = node_ijl+1
  y = y + f(node_ijl)*piunor*pjvnor*plwnor
enddo
enddo 
enddo  

end function EvalSolVeloPtCellFastAlt_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! QuickCellFindUniformGrid_DGV 
!
! This subroutine uses the assumption that the velocity mesh is uniform. 
! Aslo the function assumes the numbering convention: the frequency of numbering goes 
! w - most frequent changes, v second frequent and u last frequent 
! under this assumptions, the function returns the number of the cell in the primary mesh 
! where the provided velocity point resides. 
! If the point is outside the velocity doamin, the function returns 0  and sets outsidedomain=.true.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function QuickCellFindUniformGrid_DGV(u,v,w,outsidedomain) result (cellnum)

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi, &
				   cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, & 
                   u_L,u_R,v_L,v_R,w_L,w_R,Mu,Mv,Mw
                       
real (DP), intent (in) :: u,v,w ! provided components of the velocity 
logical, intent (out) :: outsidedomain ! = true is the point is outside the velocity doamin 
integer (I4B) :: cellnum ! result 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: ui, vi, wi ! scrap numbers on 1D meshes, 
real (DP) :: du,dv,dw,nd_ru,nd_rv,nd_rw
real (DP), parameter :: tolerance = 1.0D-10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
outsidedomain =.true.
cellnum = 0
if ( (u >= u_L - tolerance) .and. (u <= u_R + tolerance) & 
   .and. (v >= v_L - tolerance) .and. (v <= v_R + tolerance) & 
   .and. (w >= w_L - tolerance) .and. (w <= w_R + tolerance) ) then 
 ! if we are here, then the velocity point is inside the mesh
 ! find ui:
 ui = 1
 nd_ru = u_L + du
 do while ((nd_ru < u - tolerance) .and. (ui < Mu)) 
  ui = ui + 1
  nd_ru = nd_ru + du
 enddo
 if (nd_ru >= u - tolerance) then 
  vi = 1
  nd_rv = v_L + dv
  do while ((nd_rv < v - tolerance) .and. (vi < Mv)) 
   vi = vi + 1
   nd_rv = nd_rv + dv
  enddo
  if (nd_rv >= v - tolerance) then
   wi = 1
   nd_rw = w_L + dw
   do while ((nd_rw < w - tolerance) .and. (wi < Mw)) 
    wi = wi + 1
    nd_rw = nd_rw + dv
   enddo
   if (nd_rw >= w - tolerance) then
    outsidedomain =.false.
    cellnum = (ui-1)*Mv*Mw + (vi-1)*Mw + wi
   endif 
  endif 
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! got the cellnum computed
endif 

end function QuickCellFindUniformGrid_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! QuickCellFindUniformGridIJL_DGV 
!
! This subroutine uses the assumption that the velocity mesh is uniform. 
! Aslo the function assumes the numbering convention: the frequency of numbering goes 
! w - most frequent changes, v second frequent and u last frequent 
! under this assumptions, the function returns the number of the cell in the primary mesh 
! where the provided velocity point resides. 
! If the point is outside the velocity doamin, the function returns 0  and sets outsidedomain=.true.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine QuickCellFindUniformGridIJL_DGV(u,v,w,ui,vi,wi,outsidedomain) 

use DGV_commvar, only: nodes_pcell, nodes_ui, nodes_vi, nodes_wi, &
				   cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, & 
                   u_L,u_R,v_L,v_R,w_L,w_R,Mu,Mv,Mw
                       
real (DP), intent (in) :: u,v,w ! provided components of the velocity 
logical, intent (out) :: outsidedomain ! = true is the point is outside the velocity doamin 
integer (I4B), intent(out) :: ui,vi,wi ! result, dimentional addressed on the uniform velocity mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: du,dv,dw,nd_ru,nd_rv,nd_rw
real (DP), parameter :: tolerance = 1.0D-10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
outsidedomain =.true.
if ( (u >= u_L - tolerance) .and. (u <= u_R + tolerance) & 
   .and. (v >= v_L - tolerance) .and. (v <= v_R + tolerance) & 
   .and. (w >= w_L - tolerance) .and. (w <= w_R + tolerance) ) then 
 ! if we are here, then the velocity point is inside the mesh
 ! find ui:
 ui = 1
 nd_ru = u_L + du
 do while ((nd_ru < u - tolerance) .and. (ui < Mu)) 
  ui = ui + 1
  nd_ru = nd_ru + du
 enddo
 if (nd_ru >= u - tolerance) then 
  vi = 1
  nd_rv = v_L + dv
  do while ((nd_rv < v - tolerance) .and. (vi < Mv)) 
   vi = vi + 1
   nd_rv = nd_rv + dv
  enddo
  if (nd_rv >= v - tolerance) then
   wi = 1
   nd_rw = w_L + dw
   do while ((nd_rw < w - tolerance) .and. (wi < Mw)) 
    wi = wi + 1
    nd_rw = nd_rw + dv
   enddo
   if (nd_rw >= w - tolerance) then
    outsidedomain =.false.
   endif 
  endif 
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! got the cellnum computed
endif 

end subroutine QuickCellFindUniformGridIJL_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! finduiviwi(ugrid,vgrid,wgrid,u,v,w,ui,vi,wi)
! 
! This subroutine lookes trhough three one-dimensional meshes and checks whether u is in the mesh ugrid, 
! v is in the mesh vgrid, w is in the mesh wgrid. If it is, it returns the numbe greater than or equal to 1
! if a zero is returned, the number is not on the mesh
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine finduiviwi(ugrid,vgrid,wgrid,u,v,w,ui,vi,wi)

real (DP), dimension (:), intent (in) :: ugrid,vgrid,wgrid ! the grids in the three variables..
real (DP), intent (in) :: u,v,w ! the components of the velocity in question
integer (I4B), intent (out) :: ui,vi,wi  !  the coordinates of the right endpoint in each dimension. if =0  then the point is not there.. 
!!
real (DP) :: epsill=10d-11 ! a small parameter to help curb the effects of the round off error..  
!!

if ((u < ugrid(1)-epsill) .or. (u > ugrid(size(ugrid,1)) + epsill)) then 
  ui=0
else
  ui=1
  do while ((u >=ugrid(ui)-epsill) .and. (ui<size(ugrid,1)))
  ui=ui+1
  enddo
endif
!
if ((v < vgrid(1)-epsill) .or. (v > vgrid(size(vgrid,1)) + epsill)) then 
  vi=0
else
  vi=1
  do while ((v >= vgrid(vi)-epsill) .and. (vi<size(ugrid,1)))
  vi=vi+1
  enddo
endif
!
if ((w < wgrid(1)-epsill) .or. (w > wgrid(size(wgrid,1)) + epsill)) then 
  wi=0
else
  wi=1
  do while ((w >= wgrid(wi)-epsill) .and. (wi<size(wgrid,1)))
  wi=wi+1
  enddo
endif
end subroutine finduiviwi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! SetA_DGV
!
! This subroutine sets up the collision information operator A(\xi,\xi_1,\varphi^{j}_{p})
! This is a symmetric operator in $\xi$ and $\xi_{1}$, therefore we only are interested in the 
! records for unique unordered pairs (\xi,\xi_{1}). Because all velocities are indexed by one index, say i, 
! we only interested in producing records for pairs $(\xi_{i},\xi_{j})$, $i>j$. Notice that the case \xi=\xi_{i} 
! gives value 0.
!   
! This subroutine depends on the main program
! 
! You should consult the notes n101011.tex for detail on the formulas. 
!
! Evaluation of each entry of A involves a two dimensional integral. If the entry is smaller than some given number, 
! the entry is neglected/nullified by force
! 
! The structure of the arrays related to the operator A: 
! A_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p}) for each basis function with number i
! for each non-zero entry of A this one keeps the index of velocities that produced that non-zero entries. 
! A(i) i -- is the index of nonzero entries, we need to use other A-arrays to restore 
! what velocities and what basis function this index corresponds  
! for example, A_sind_xi(i) gives the first velocity, A_sind_xi1 gives the second velcity
! and A_phi gives the index of the used basis function. 
!                                    
! A(\xi,\xi_{1};\varphi^{j}_{p})=\frac{d^2 |g|}{8}  \int_{0}^{\pi} \int_{0}^{2\pi}
!(\varphi^{j}_{p}(\xi')+\varphi^{j}_{p}(\xi'_{1}))
! d\varepsilon\, \sin \chi d \chi - \frac{d^2 |g| \pi }{2} (\varphi^{j}_{p}(\xi)+\varphi^{j}_{p}(\xi_{1})) 
!
! Takes 
! Trad  == This variable determines the cut off radius for $A$. It should be based on the estimated size of 
!                               ! non-trivial support for the distribution function
! ErrChi == the parameter defining how accurate the evaluation of the integral in Chi  should be. 
! ErrEps == the parameter defining how accurate the evaluation of the integral in Epsilon should be. 
! min_sens == minimum accepted value for A, values below will be neglected
! I1_list == List of basis functions for which A-array is evaluated.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetA_DGV
use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w, nodes_ui, nodes_vi, nodes_wi, &
                   nodes_gwts,A_capphi,A_xi,A_xi1,A_phi,A,&
                   cells_lu, cells_lv, cells_lw, cells_ru, cells_rv, cells_rw,&
                   Trad, ErrChi, ErrEps, min_sens, I1_list,Num_OMP_threads

use gaussian_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), parameter :: pi25DT = 3.141592653589793238462643d0
integer (I4B), parameter :: NdsChi = 10000 ! The max number of cells for integration in chi and epsilon
real (DP) :: ires_e=100.0_DP ! this parameter will determine the resolution for integration in \varepsilon
real (DP) :: ires_t=100.0_DP ! this parameter will determine the resolution in t.
                                   ! the number of points is dictated by the radius of the collision sphere, which is |g|/2
                                   ! the bigger is the radius, the more points is needed to integrate over the sphere.
                                   ! these parameters will represent a lengh of the arc that is desired for application of 
                                   ! a gauss quadrature rule. Then the angular intervals will be broken in portions no bigger than 
                                   ! ires_t(_e)/|g| to ensure sufficient angular resolution
!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (2*NdsChi) :: NodesChi1, NodesChi2 
real (DP), dimension (3*NdsChi) :: FuncChi1, FuncChi2 ! Arrays to store cells and functin values on cells Two copies are needed. Integration in Chi
integer (I4B), dimension (NdsChi) :: CellsChiRef1,CellsChiRef2 ! Arrays to keep refinement flags
!!
integer (I4B) :: A_ct, phi_ct, loc_alloc_stat ! scrap variables. A_ct is a counter of how many records have been created for A
							!phi_ct counts how many records have been created for this particular basis function.  
integer (I4B) :: An ! this varable will keep the current size of the array A.
integer (I4B) :: ni ! this one will keep the size of the nodes array
integer (I4B) :: i1,i2,i3,d,e ! local counters  
integer (I4B) :: pci ! local integer
real (DP), dimension (8) :: uu,vv,ww ! coordinates of the 8 vertices of the support of the basis function
real (DP) :: dphi, dsph ! diameters for the circumscribed sphere for basis function and the collision shpere
real (DP) :: xiu,xiv,xiw,xi1u,xi1v,xi1w,xiupr,xivpr,&
             xiwpr,xi1upr,xi1vpr,xi1wpr,ugu,ugv,ugw ! scrap variables to keep the values of the pre and post collision velocities 
real (DP) :: ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! more screap variables
real (DP) :: ku,kv,kw ! coordinates of the post-coll g'
real (DP) :: dist2, g1,g2 ! useful variables 
integer (I4B)  :: Nchi1 ! numbers of subdivisions in epsilon and chi
real (DP) :: Atemp,quad,quad1,quad2,es,int1,varphi_xi,varphi_xi1 ! scrap variables.
logical :: chiadaptflag
!
integer :: iiii ! a test variable to play with OpenMP runtime functions calls
!!!!!!!!!!!!!!!!!!!!!!!!!! Interface for OpenMP runtime libraries !!!!!!!!!!!!!!!!!!!
interface 
 function omp_get_thread_num() result (y)
   integer :: y 
 end function omp_get_thread_num
 function omp_get_num_threads() result (y)
  integer :: y 
 end function omp_get_num_threads 
 function omp_get_num_procs() result (y)
  integer :: y 
 end function omp_get_num_procs
 function omp_get_stack_size () result (y)
  use nrtype
  integer (I2B) :: y
 end function omp_get_stack_size
end interface  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial Allocation of the A-arrays 
ni=size(nodes_pcell,1)
allocate (A_capphi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetA_DGV: Allocation error for variable  A_capphi"
  stop
  end if
An=10*ni
allocate (A_xi(1:An), A_xi1(1:An), A(1:An), A_phi(1:An), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetA_DGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A_xi=0; A_xi1=0; A=0; A_phi=0; A_capphi=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A_ct=0 ! in the beginning there is zero records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!TEMP: do it just for one velocity nodes:  do i1=1,ni  ! loop in basis functions 
i1=I1_list(1)
!!
phi_ct=0 ! in the beginning there is zero records
pci=nodes_pcell(i1)
!
uu(1)=cells_lu(pci); vv(1)=cells_lv(pci); ww(1)=cells_lw(pci) ! the first one holds the lower left left corner
uu(2)=cells_ru(pci); vv(2)=cells_rv(pci); ww(2)=cells_rw(pci) ! the second holds the upper right right corner
uu(3)=cells_lu(pci); vv(3)=cells_lv(pci); ww(3)=cells_rw(pci)
uu(4)=cells_lu(pci); vv(4)=cells_rv(pci); ww(4)=cells_lw(pci)
uu(5)=cells_lu(pci); vv(5)=cells_rv(pci); ww(5)=cells_rw(pci)
uu(6)=cells_ru(pci); vv(6)=cells_lv(pci); ww(6)=cells_lw(pci)
uu(7)=cells_ru(pci); vv(7)=cells_lv(pci); ww(7)=cells_rw(pci)
uu(8)=cells_ru(pci); vv(8)=cells_rv(pci); ww(8)=cells_lw(pci)
dphi = sqrt((uu(2)-uu(1))**2+(vv(2)-vv(1))**2+(ww(2)-ww(1))**2)
!! we set the resolution parameters based on the size of the velocity cell.
ires_t = max(uu(2)-uu(1),vv(2)-vv(1),ww(2)-ww(1))/32
ires_e = ires_t
!!!
! Add openMP parallel directives here....  
!!!
! iiii=omp_get_num_procs
! OpenMP set the number of threads: 

call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(xiu,xiv,xiw,xi1u,xi1v,xi1w,dsph,dist2,ugu,ugv,ugw, & 
!$OMP    ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,Nchi1,NodesChi1,FuncChi1,d, & 
!$OMP    CellsChiRef1,Atemp,ChiAdaptFlag,int1,NodesChi2,FuncChi2,CellsChiRef2,e, &
!$OMP    quad,quad1,quad2,es,i2,i3,iiii,varphi_xi,varphi_xi1) NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 5)  
do i2=1,ni  ! loop in velocity 1
!!!!
   xiu=nodes_u(i2); xiv=nodes_v(i2); xiw=nodes_w(i2)
   varphi_xi = EvalLagrBasisFunByNdsDGblzmEZ(xiu,xiv,xiw,i1) ! this is the value of the basis function on the first velocity -- will be passed to integrator
do i3=i2+1,ni ! loop in velocity 2
   xi1u=nodes_u(i3); xi1v=nodes_v(i3); xi1w=nodes_w(i3)
   varphi_xi1 = EvalLagrBasisFunByNdsDGblzmEZ(xi1u,xi1v,xi1w,i1) ! this is the value of the basis function on the second velocity -- -- will be passed to integrator
!!!!!!!!!!!! evaluation of A !!!!!!!!!!!!!!!
 ! evaluate the diamieter of the collision sphere and the double distance from the center of the basis function support to the center of the collision sphere
   dsph = sqrt((xiu-xi1u)**2+(xiv-xi1v)**2+(xiw-xi1w)**2)
   dist2 = sqrt((xiu+xi1u-uu(1)-uu(2))**2+(xiv+xi1v-vv(1)-vv(2))**2+(xiw+xi1w-ww(1)-ww(2))**2)
 ! Now we calculate some useful quantities:
   ugu = (xiu-xi1u)/dsph  ! need to introduce a unit vector in the direction of g
   ugv = (xiv-xi1v)/dsph  
   ugw = (xiw-xi1w)/dsph
   ! The components of the vector (ugu,ugv,ugw) should be orthogonal to vecotor g. 
   ! there is only one case when these componets are degenerate. is it when g is parallel to (1,1,1)
   ! in this case, we need to cook a non-trivial vecotor. The next if statement takes case of that: 
   if (abs(ugw-ugv)+abs(ugu-ugw)+abs(ugv-ugu) < 1.0d-6) then 
   ! first, we create two orthogonal non-trivial vectors
   ugux2 = ugv 
   ugvx2 = -ugu 
   ugwx2 = 0
   ugux3 = ugu*ugw
   ugvx3 = ugv*ugw
   ugwx3 = -(ugu)**2-(ugv)**2
   ! 
   g1=sqrt(ugux2**2+ugvx2**2)
   g2=sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   !
   else 
   !first we create the components of two non-trivial orthogonal vectors: 
   ugux2 = (ugw-ugv)
   ugvx2 = (ugu-ugw)
   ugwx2 = (ugv-ugu)
   ugux3 = ((ugv)**2-(ugv)*(ugu)-(ugw)*(ugu)+(ugw)**2)
   ugvx3 = ((ugw)**2-(ugw)*(ugv)-(ugu)*(ugv)+(ugu)**2)
   ugwx3 = ((ugu)**2-(ugu)*(ugw)-(ugv)*(ugw)+(ugv)**2)
   !  
   g1 = sqrt(ugux2**2 + ugvx2**2 + ugwx2**2)
   g2 = sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   end if
!  First we estimate if A is zero by cheking the overlap of spheres
   if ((dist2 > dsph+dphi) .or. (dist2+dphi < dsph) .or. (dsph > Trad)) then 
        cycle ! collision shpere does not hit the support of the basis function continue with the next velocity
              ! or the collision shpere is so large that it is not possible for distribution function to be nonzero for both velocities 
   end if 
 ! if we got here then the collision sphere has some possible overlap with the support of basis function
 ! To evaluate the integral, we will implement adaptive quadrature in both directions using Simpson's rule. 
 !
 ! begin integration in \chi
 ! first, we set up the initial mesh
   Nchi1=FLOOR(pi25DT*dsph/ires_t/2.0_DP)+1
   if (Nchi1>NdsChi) then 
     Nchi1=NdsChi
     print *,"SetA_DGV: Warning! Number of nodes (NdsChi) for integration in epsilon gives insufficient resolution"  
   end if 
   ! we set the initial mesh
   do d=1,Nchi1
    NodesChi1(2*(d-1)+1)=(d-1)*pi25DT/Real(Nchi1,DP)
    NodesChi1(2*d)=d*pi25DT/Real(Nchi1,DP)
   end do
   ! now we evaluate the integrand on the initial mesh
   FuncChi1=0
   do d=1,Nchi1
   ! each d gives one cell
   ! First, we evaluate the integrand on the left node of the cell
    if (d>1) then
     if (NodesChi1(2*d-1) == NodesChi1(2*d-2)) then 
      FuncChi1(3*d-2)= FuncChi1(3*d-3) ! if the node is repeating, the value has been computed already
     else 
      FuncChi1(3*d-2) = A_IntEpsilon(NodesChi1(2*d-1),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
     end if
    else 
     FuncChi1(3*d-2) = A_IntEpsilon(NodesChi1(2*d-1),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    end if
   ! This takes care of the left node in the cell... 
   !  
   !Now we evaluate the solution on the right node of the cell and at the center. 
   FuncChi1(3*d) = A_IntEpsilon(NodesChi1(2*d),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   FuncChi1(3*d-1) = A_IntEpsilon((NodesChi1(2*d-1) + NodesChi1(2*d))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,&
              ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   ! repeat for all cells... 
   end do      
   ! finally, we mark all cells for the refinement in the beginning... 
   CellsChiRef1=0 
   CellsChiRef1(1:Nchi1)=1
   ! 
   !Next we adaptively refine and evaluate the integral:    
   Atemp=0
   ChiAdaptFlag =.true.
   int1 = 100*ErrChi ! initialize with something big for the intial iteration. 
   do while (ChiAdaptFlag)
   ChiAdaptFlag =.false.
   ! Now we will go over the array of the integrand values and will cross 
   ! out cells where the integrand is zero, and those that are not marked for the refinment 
   ! on those cells we summ the integrand and keep the sum.
   ! Cells marked for the refinement we divide into halves 
   NodesChi2=0
   FuncChi2=0
   CellsChiRef2=0
   e=0
   if (int1>ErrChi) then  ! this checks that the integral over the mesh NodesChi1 will be small anyway
   int1=0 ! this will keep the estimate of the integral on the mesh NodesChi1 so that we stop resolving when this number is small
   do d=1,Nchi1
   if ((abs(FuncChi1(3*d-2))+abs(FuncChi1(3*d-1))+abs(FuncChi1(3*d)) > 1.0d-15) .and. (CellsChiRef1(d)==1)) then 
    if (e > NdsChi-2) then 
     print *,"SetA_DGV: Number of nodes in chi is too big (e>NdsChi-2)"
     stop
    end if
    ! save the nodes and the integrand values
    NodesChi2(2*e+1)=NodesChi1(2*d-1)
    FuncChi2(3*e+1)=FuncChi1(3*d-2)
    NodesChi2(2*e+4)=NodesChi1(2*d)
    FuncChi2(3*e+6)=FuncChi1(3*d)
    ! Evaluate the midpoint node: 
    NodesChi2(2*e+2)=(NodesChi1(2*d-1)+NodesChi1(2*d))/2
    NodesChi2(2*e+3)=NodesChi2(2*e+2)
    ! Save the midpoint value of the integrand:
    FuncChi2(3*e+3)=FuncChi1(3*d-1)
    FuncChi2(3*e+4)=FuncChi2(3*e+3)
    ! Evaluate the integrand on the new nodes
    FuncChi2(3*e+2) = A_IntEpsilon((NodesChi2(2*e+2)+NodesChi2(2*e+1))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    FuncChi2(3*e+5) = A_IntEpsilon((NodesChi2(2*e+4)+NodesChi2(2*e+3))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    ! Now it is time to compare the quadratures and to mark cells for the refinement 
    quad  = (FuncChi1(3*d-2)+4*FuncChi1(3*d-1)+FuncChi1(3*d))*(NodesChi1(2*d)-NodesChi1(2*d-1))/6
    quad1 = (FuncChi2(3*e+1)+4*FuncChi2(3*e+2)+FuncChi2(3*e+3))*(NodesChi2(2*e+2)-NodesChi2(2*e+1))/6
    quad2 = (FuncChi2(3*e+4)+4*FuncChi2(3*e+5)+FuncChi2(3*e+6))*(NodesChi2(2*e+4)-NodesChi2(2*e+3))/6
    es = abs(quad1+quad2-quad)/15/(NodesChi1(2*d)-NodesChi1(2*d-1)) ! error indicator
    int1=int1+abs(quad1)+abs(quad2) ! this quantity estimates the integral over the entire mesh NodesChi2 that will be NodesChi1 soon.. 
    if (es > ErrChi) then 
     ChiAdaptFlag = .true.
     CellsChiRef2(e+1) = 1
     CellsChiRef2(e+2) = 1
    else
     CellsChiRef2(e+1) = 0
     CellsChiRef2(e+2) = 0
     Atemp = Atemp + quad1 + quad2 
    end if
    e=e+2
   ! end refinement 
   end if 
   end do 
   !
   end if ! end of the check that the integral over the mesh NodesChi1 will be small anyway
   ! Finally, we need to replace the first mesh with the refined mesh 
   FuncChi1=0
   FuncChi1=FuncChi2
   NodesChi1=0
   NodesChi1=NodesChi2
   CellsChiRef1=0
   CellsChiRef1=CellsChiRef2
   Nchi1=e
   end do
   !!!!!!!!!!!!!!!!!! Looks like we are done with the integration in chi !!!!! 
   ! now it is time to check if the value of A is non-zero. If it is bigger than some specified level of 
   ! minimum sensitivity, then is it added to storage, Otherwise we continue to the next velocity #2.
   if ( ABS(Atemp) > min_sens ) then     
!$omp critical    
     A_ct = A_ct+1               ! we count the record for A
     phi_ct=phi_ct + 1           ! we count the record for this basis function. 
     if (An < A_ct) then
        call ExtendAarraysDGV(An,ni)
     end if  
      ! now we need to add the volume elements for xi and xi1 
      ! ATTENTION: molecular diameter mol_diam is removed in the dimensionless formulation)
      ! ATTENTION: for evaluation of A operator using the dimensional code use mol_diam = 1.0:
      ! ATTENTION: OLD (dimenional) code Atemp = Atemp*mod_diam^2*(nodes_gwts(i3))**(nodes_gwts(i2))
      ! Dimensionless code. Molecular diamter is accounted for in the spatial operator. 
      ! xi1:
      Atemp = Atemp*(nodes_gwts(i3)) ! this takes care of the volume elements and the weight for xi1 
      ! xi			
      Atemp = Atemp*(nodes_gwts(i2)) ! this takes care of the volume elements and the weight for xi 
      !		
      dsph=dsph/8.0_DP  ! this takes care of the fraction 1/8 still need to add |g|
      A(A_ct) = Atemp*dsph ! this takes care of |g|/8
      A_xi(A_ct) = i2
      A_xi1(A_ct) = i3 
      A_phi(A_ct) = i1
!!
!iiii = omp_get_thread_num()
!print *, "I2=", i2, "Thread", iiii, "A_ct=", A_ct, i3        
!$omp end critical           
   end if                    
!!!!!!!!!!!! end evaluation of A !!!!!!!!!!!
end do  ! END LOOP in I3
! add this two lines to track progress and Print a hello message to Check if parallel stuff worked... !!!!
!iiii = omp_get_thread_num()
!print *, "I2=", i2, "Thread", iiii, "A_ct=%i8", A_ct    
! 
end do ! END LOOP IN I2
A_capphi(i1)=phi_ct
print *, "Set_A i1=", i1, "A_ct=", A_ct
!! end do !! TEMPORARY do it for just one node...  
!!!!!!!!!!!!
call ShrinkAarraysDGV(A_ct)
!!!!!!!!!!!!
end subroutine SetA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendAarraysDGV(An,ni)
!
! This subroutine extends arrays A, A_phi,A_xi,A_xi1 for additional ni records

subroutine ExtendAarraysDGV(An,ni)

use DGV_commvar, only: A,A_phi,A_xi,A_xi1

integer (I4B), intent (out) :: An ! An is the length of the arrays that will be updated.  
integer (I4B), intent (in) :: ni ! ni is the number of records to add

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B) , dimension (:), allocatable :: A_iscr  ! integer scrap array
integer (I4B) :: nn ! scrap variable
!
integer (I4B) :: loc_alloc_stat
!
nn=size(A,1)
An=nn+ni
! extending A ... 
allocate (A_rscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_rscr, nn=", nn
  stop
  end if
A_rscr=A
deallocate(A)
allocate (A(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A, nn=", nn+ni
  stop
  end if
A(1:nn)=A_rscr
deallocate(A_rscr)
! end extending A

! extending A_xi
allocate(A_iscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_iscr, nn=", nn
  stop
  end if
A_iscr=A_xi
deallocate(A_xi)
allocate (A_xi(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi, nn=", nn+ni
  stop
  end if
A_xi(1:nn)=A_iscr
! extending A_xi1
A_iscr=A_xi1
deallocate(A_xi1)
allocate (A_xi1(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi1, nn=", nn+ni
  stop
  end if
A_xi1(1:nn)=A_iscr
!extending A_phi
A_iscr=A_phi
deallocate(A_phi)
allocate (A_phi(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_phi, nn=", nn+ni
  stop
  end if
A_phi(1:nn)=A_iscr
! end extending A_phi
deallocate(A_iscr)

end subroutine ExtendAarraysDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendAarraysDGVII(An,ni)
!
! This subroutine extends arrays AII, A_phiII,A_xiII,A_xi1II for additional ni records
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendAarraysDGVII(An,ni)

use DGV_commvar, only: AII, A_phiII, A_xiII, A_xi1II

integer (I4B), intent (out) :: An ! An is the length of the arrays that will be updated.  
integer (I4B), intent (in) :: ni ! ni is the number of records to add

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B) , dimension (:), allocatable :: A_iscr  ! integer scrap array
integer (I4B) :: nn ! scrap variable
!
integer (I4B) :: loc_alloc_stat
!
nn=size(AII,1)
An=nn+ni
! extending A ... 
allocate (A_rscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_rscr, nn=", nn
  stop
  end if
A_rscr=AII
deallocate(AII)
allocate (AII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A, nn=", nn+ni
  stop
  end if
AII(1:nn)=A_rscr
deallocate(A_rscr)
! end extending A

! extending A_xi
allocate(A_iscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_iscr, nn=", nn
  stop
  end if
A_iscr=A_xiII
deallocate(A_xiII)
allocate (A_xiII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi, nn=", nn+ni
  stop
  end if
A_xiII(1:nn)=A_iscr
! extending A_xi1
A_iscr=A_xi1II
deallocate(A_xi1II)
allocate (A_xi1II(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_xi1, nn=", nn+ni
  stop
  end if
A_xi1II(1:nn)=A_iscr
!extending A_phi
A_iscr=A_phiII
deallocate(A_phiII)
allocate (A_phiII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAarraysDGV: Allocation error for variable  A_phi, nn=", nn+ni
  stop
  end if
A_phiII(1:nn)=A_iscr
! end extending A_phi
deallocate(A_iscr)

end subroutine ExtendAarraysDGVII

subroutine ShrinkAarraysDGV(ni)

use DGV_commvar, only: A,A_phi,A_xi,A_xi1

integer (I4B), intent (in) :: ni ! ni is the new size of the arrays to be shrunk

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B), dimension (:), allocatable :: A_iscr  ! integer scrap array
!
integer (I4B) :: loc_alloc_stat 

! shrinking A ... 
allocate (A_rscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_rscr, nn=", ni
  stop
  end if
A_rscr=A(1:ni)
deallocate(A)
allocate (A(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A, nn=", ni
  stop
  end if
A=A_rscr
deallocate(A_rscr)
! end shrinking A

! shrinking A_xi
allocate (A_iscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_iscr, nn=", ni
  stop
  end if
A_iscr=A_xi(1:ni)
deallocate(A_xi)
allocate (A_xi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xi, nn=", ni
  stop
  end if
A_xi=A_iscr
! extending A_xi1
A_iscr=A_xi1(1:ni)
deallocate(A_xi1)
allocate (A_xi1(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xi1, nn=", ni
  stop
  end if
A_xi1=A_iscr
!extending A_phi
A_iscr=A_phi(1:ni)
deallocate(A_phi)
allocate (A_phi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_phi, nn=", ni
  stop
  end if
A_phi=A_iscr
! end extending A_phi
deallocate(A_iscr)
end subroutine ShrinkAarraysDGV

subroutine ShrinkAarraysDGVII(ni)

use DGV_commvar, only: AII,A_phiII,A_xiII,A_xi1II

integer (I4B), intent (in) :: ni ! ni is the new size of the arrays to be shrunk

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B), dimension (:), allocatable :: A_iscr  ! integer scrap array
!
integer (I4B) :: loc_alloc_stat 

! shrinking AII ... 
allocate (A_rscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGVII: Allocation error for variable  A_rscr, nn=", ni
  stop
  end if
A_rscr=AII(1:ni)
deallocate(AII)
allocate (AII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGVII: Allocation error for variable  AII, nn=", ni
  stop
  end if
AII=A_rscr
deallocate(A_rscr)
! end shrinking AII

! shrinking A_xiII
allocate (A_iscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_iscr, nn=", ni
  stop
  end if
A_iscr=A_xiII(1:ni)
deallocate(A_xiII)
allocate (A_xiII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xiII, nn=", ni
  stop
  end if
A_xiII=A_iscr
! shrinking A_xi1II
A_iscr=A_xi1II(1:ni)
deallocate(A_xi1II)
allocate (A_xi1II(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_xi1II, nn=", ni
  stop
  end if
A_xi1II=A_iscr
!shrinking A_phiII
A_iscr=A_phiII(1:ni)
deallocate(A_phiII)
allocate (A_phiII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAarraysDGV: Allocation error for variable  A_phi, nn=", ni
  stop
  end if
A_phiII=A_iscr
! end shrinking A_phiII
deallocate(A_iscr)
end subroutine ShrinkAarraysDGVII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! SetAKorobov_DGV
!
! This subroutine sets up the collision information operator A(\xi,\xi_1,\varphi^{j}_{p}) for the case 
! when \xi and \xi_1 run over the korobov nodes.
!   
! This subroutine depends on the main program
! 
! Evaluation of each entry of A involves a two dimensional integral. If the entry is smaller than some given number, 
! the entry is neglected/nullified by force
! 
! The structure of the arrays related to the operator A: 
!
! Akor_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p}) for each basis function with number i
! for each non-zero entry of Akor this one keeps the index of velocities that produced that non-zero entries. 
!
! Akor(i) i -- is the index of nonzero entries, we need to use other A-arrays to restore 
! what velocities and what basis function this index corresponds  
! for example, Akor_k(i) gives the number of the korobov node for which Akor(i) was computed 
! Akor_phi(i) gives the index of the used basis function. 
!                                    
! A(\xi,\xi_{1};\varphi^{j}_{p})=\frac{d^2 |g|}{8}  \int_{0}^{\pi} \int_{0}^{2\pi}
!(\varphi^{j}_{p}(\xi')+\varphi^{j}_{p}(\xi'_{1}))
! d\varepsilon\, \sin \chi d \chi - \frac{d^2 |g| \pi }{2} (\varphi^{j}_{p}(\xi)+\varphi^{j}_{p}(\xi_{1})) 
!
! Takes 
! Trad  == This variable determines the cut off radius for $A$. It should be based on the estimated size of 
!                               ! non-trivial support for the distribution function
! ErrChi == the parameter defining how accurate the evaluation of the integral in Chi  should be. 
! ErrEps == the parameter defining how accurate the evaluation of the integral in Epsilon should be. 
! min_sens == minimum accepted value for A, values below will be neglected
! I1_list == List of basis functions for which A-array is evaluated.
! I2_from and I2_to -- firt and last kobov node to compute
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetAKorobov_DGV

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w, nodes_ui, nodes_vi, nodes_wi, &
                   nodes_gwts,Akor_capphi,Akor,Akor_k, Akor_phi, &
                   cells_lu, cells_lv, cells_lw, cells_ru, cells_rv, cells_rw,&
                   Trad, ErrChi, ErrEps, min_sens, I1_list, I2_from, I2_to, Num_OMP_threads,& 
                   korob_net_param, u_L,u_R,v_L,v_R,w_L,w_R

use gaussian_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), parameter :: pi25DT = 3.141592653589793238462643d0
integer (I4B), parameter :: NdsChi = 10000 ! The max number of cells for integration in chi and epsilon
real (DP) :: ires_e=100.0_DP ! this parameter will determine the resolution for integration in \varepsilon
real (DP) :: ires_t=100.0_DP ! this parameter will determine the resolution in t.
                                   ! the number of points is dictated by the radius of the collision sphere, which is |g|/2
                                   ! the bigger is the radius, the more points is needed to integrate over the sphere.
                                   ! these parameters will represent a lengh of the arc that is desired for application of 
                                   ! a gauss quadrature rule. Then the angular intervals will be broken in portions no bigger than 
                                   ! ires_t(_e)/|g| to ensure sufficient angular resolution
!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (2*NdsChi) :: NodesChi1, NodesChi2 
real (DP), dimension (3*NdsChi) :: FuncChi1, FuncChi2 ! Arrays to store cells and functin values on cells Two copies are needed. Integration in Chi
integer (I4B), dimension (NdsChi) :: CellsChiRef1,CellsChiRef2 ! Arrays to keep refinement flags
!!
integer (I4B) :: A_ct, phi_ct, loc_alloc_stat ! scrap variables. A_ct is a counter of how many records have been created for A
							!phi_ct counts how many records have been created for this particular basis function.  
integer (I4B) :: An ! this varable will keep the current size of the array A.
integer (I4B) :: ni ! this one will keep the size of the nodes array
integer (I4B) :: i1,i2,d,e ! local counters  
integer (I4B) :: pci ! local integer
real (DP), dimension (8) :: uu,vv,ww ! coordinates of the 8 vertices of the support of the basis function
real (DP) :: dphi, dsph ! diameters for the circumscribed sphere for basis function and the collision shpere
real (DP) :: xiu,xiv,xiw,xi1u,xi1v,xi1w,xiupr,xivpr,&
             xiwpr,xi1upr,xi1vpr,xi1wpr,ugu,ugv,ugw ! scrap variables to keep the values of the pre and post collision velocities 
real (DP) :: ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! more screap variables
real (DP) :: ku,kv,kw ! coordinates of the post-coll g'
real (DP) :: dist2, g1,g2 ! useful variables 
integer (I4B)  :: Nchi1 ! numbers of subdivisions in epsilon and chi
real (DP) :: Atemp,quad,quad1,quad2,es,int1,varphi_xi,varphi_xi1,frac_part ! scrap variables.
logical :: chiadaptflag
!
integer :: iiii ! a test variable to play with OpenMP runtime functions calls
!!!!!!!!!!!!!!!!!!!!!!!!!! Interface for OpenMP runtime libraries !!!!!!!!!!!!!!!!!!!
interface 
 function omp_get_thread_num() result (y)
   integer :: y 
 end function omp_get_thread_num
 function omp_get_num_threads() result (y)
  integer :: y 
 end function omp_get_num_threads 
 function omp_get_num_procs() result (y)
  integer :: y 
 end function omp_get_num_procs
 function omp_get_stack_size () result (y)
  use nrtype
  integer (I2B) :: y
 end function omp_get_stack_size
end interface  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial Allocation of the A-arrays 
ni=size(nodes_pcell,1)
allocate (Akor_capphi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetAKorobov_DGV: Allocation error for variable  Akor_capphi"
  stop
  end if
An = I2_to-I2_from+1 !This is the number of the Korobov points that needs to be processed
allocate (Akor(1:An), Akor_k(1:An), Akor_phi(1:An), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetAKorobov_DGV: Allocation error for variable  Akor, Akor_k, Akor_phi"
  stop
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Akor=0; Akor_k=0; Akor_phi=0; Akor_capphi=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
A_ct=0 ! in the beginning there is zero records
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Uncomment the loop if want to compute for several basuis functions in one take. Otherwise 
!! compute recoprds for one basis functions at a time. 
!! uncomment if want to compute components for velocity nodes/basis functions I1_list(1),..., I1_list(last) !!:!!   do ixxx=1:size(I1_list,1)  ! loop in basis functions 
!!  
i1=I1_list(1)   ! if using the loop above, replace with i1=I1_list(ixxx)
!!
phi_ct=0 ! in the beginning there is zero records++
pci=nodes_pcell(i1)
!
uu(1)=cells_lu(pci); vv(1)=cells_lv(pci); ww(1)=cells_lw(pci) ! the first one holds the lower left left corner
uu(2)=cells_ru(pci); vv(2)=cells_rv(pci); ww(2)=cells_rw(pci) ! the second holds the upper right right corner
uu(3)=cells_lu(pci); vv(3)=cells_lv(pci); ww(3)=cells_rw(pci)
uu(4)=cells_lu(pci); vv(4)=cells_rv(pci); ww(4)=cells_lw(pci)
uu(5)=cells_lu(pci); vv(5)=cells_rv(pci); ww(5)=cells_rw(pci)
uu(6)=cells_ru(pci); vv(6)=cells_lv(pci); ww(6)=cells_lw(pci)
uu(7)=cells_ru(pci); vv(7)=cells_lv(pci); ww(7)=cells_rw(pci)
uu(8)=cells_ru(pci); vv(8)=cells_rv(pci); ww(8)=cells_lw(pci)
dphi = sqrt((uu(2)-uu(1))**2+(vv(2)-vv(1))**2+(ww(2)-ww(1))**2)
!! we set the resolution parameters based on the size of the velocity cell.
ires_t = max(uu(2)-uu(1),vv(2)-vv(1),ww(2)-ww(1))/32
ires_e = ires_t
!!!
! Add openMP parallel directives here....  
!!!
! iiii=omp_get_num_procs
! OpenMP set the number of threads: 

 call omp_set_num_threads(Num_OMP_threads) ! overrides enviromental variable. if commented will use enviromental variable instead
!$OMP PARALLEL DO PRIVATE(xiu,xiv,xiw,xi1u,xi1v,xi1w,dsph,dist2,ugu,ugv,ugw, & 
!$OMP    ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,Nchi1,NodesChi1,FuncChi1,d, & 
!$OMP    CellsChiRef1,Atemp,ChiAdaptFlag,int1,NodesChi2,FuncChi2,CellsChiRef2,e, &
!$OMP    quad,quad1,quad2,es,i2,iiii,varphi_xi,varphi_xi1,frac_part) NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 5)  
do i2=I2_from,I2_to  ! loop in specified Korobov nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we will compute a korobov velocity node in (R^3)x(R^3)
!
frac_part = Real( (korob_net_param(2)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(2)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xiu = (frac_part - 0.5_DP)*(u_R-u_L)+(u_R+u_L)/2.0_DP 
frac_part = Real( (korob_net_param(3)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(3)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xiv = (frac_part - 0.5_DP)*(v_R-v_L)+(v_R+v_L)/2.0_DP 
frac_part = Real( (korob_net_param(4)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(4)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xiw = (frac_part - 0.5_DP)*(w_R-w_L)+(w_R+w_L)/2.0_DP 
frac_part = Real( (korob_net_param(5)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(5)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xi1u = (frac_part - 0.5_DP)*(u_R-u_L)+(u_R+u_L)/2.0_DP 
frac_part = Real( (korob_net_param(6)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(6)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xi1v = (frac_part - 0.5_DP)*(v_R-v_L)+(v_R+v_L)/2.0_DP 
frac_part = Real( (korob_net_param(7)*i2) ,DP)/Real(korob_net_param(1),DP)  - & 
         Real(FLOOR(( Real((korob_net_param(7)*i2),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
xi1w = (frac_part - 0.5_DP)*(w_R-w_L)+(w_R+w_L)/2.0_DP 
!        
!!!!
   varphi_xi = EvalLagrBasisFunByNdsDGblzmEZ(xiu,xiv,xiw,i1) ! this is the value of the basis function on the first velocity -- will be passed to integrator
   varphi_xi1 = EvalLagrBasisFunByNdsDGblzmEZ(xi1u,xi1v,xi1w,i1) ! this is the value of the basis function on the second velocity -- -- will be passed to integrator
!!!!!!!!!!!! evaluation of A !!!!!!!!!!!!!!!
 ! evaluate the diamieter of the collision sphere and the double distance from the center of the basis function support to the center of the collision sphere
   dsph = sqrt((xiu-xi1u)**2+(xiv-xi1v)**2+(xiw-xi1w)**2)
   dist2 = sqrt((xiu+xi1u-uu(1)-uu(2))**2+(xiv+xi1v-vv(1)-vv(2))**2+(xiw+xi1w-ww(1)-ww(2))**2)
 ! Now we calculate some useful quantities:
   ugu = (xiu-xi1u)/dsph  ! need to introduce a unit vector in the direction of g
   ugv = (xiv-xi1v)/dsph  
   ugw = (xiw-xi1w)/dsph
   ! The components of the vector (ugu,ugv,ugw) should be orthogonal to vecotor g. 
   ! there is only one case when these componets are degenerate. is it when g is parallel to (1,1,1)
   ! in this case, we need to cook a non-trivial vecotor. The next if statement takes case of that: 
   if (abs(ugw-ugv)+abs(ugu-ugw)+abs(ugv-ugu) < 1.0d-6) then 
   ! first, we create two orthogonal non-trivial vectors
   ugux2 = ugv 
   ugvx2 = -ugu 
   ugwx2 = 0
   ugux3 = ugu*ugw
   ugvx3 = ugv*ugw
   ugwx3 = -(ugu)**2-(ugv)**2
   ! 
   g1=sqrt(ugux2**2+ugvx2**2)
   g2=sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   !
   else 
   !first we create the components of two non-trivial orthogonal vectors: 
   ugux2 = (ugw-ugv)
   ugvx2 = (ugu-ugw)
   ugwx2 = (ugv-ugu)
   ugux3 = ((ugv)**2-(ugv)*(ugu)-(ugw)*(ugu)+(ugw)**2)
   ugvx3 = ((ugw)**2-(ugw)*(ugv)-(ugu)*(ugv)+(ugu)**2)
   ugwx3 = ((ugu)**2-(ugu)*(ugw)-(ugv)*(ugw)+(ugv)**2)
   !  
   g1 = sqrt(ugux2**2 + ugvx2**2 + ugwx2**2)
   g2 = sqrt(ugux3**2+ugvx3**2+ugwx3**2)
   end if
!  First we estimate if A is zero by cheking the overlap of spheres
   if ((dist2 > dsph+dphi) .or. (dist2+dphi < dsph) .or. (dsph > Trad)) then 
        cycle ! collision shpere does not hit the support of the basis function continue with the next velocity
              ! or the collision shpere is so large that it is not possible for distribution function to be nonzero for both velocities 
   end if 
 ! if we got here then the collision sphere has some possible overlap with the support of basis function
 ! To evaluate the integral, we will implement adaptive quadrature in both directions using Simpson's rule. 
 !
 ! begin integration in \chi
 ! first, we set up the initial mesh
   Nchi1=FLOOR(pi25DT*dsph/ires_t/2.0_DP)+1
   if (Nchi1>NdsChi) then 
     Nchi1=NdsChi
     print *,"SetA_DGV: Warning! Number of nodes (NdsChi) for integration in epsilon gives insufficient resolution"  
   end if 
   ! we set the initial mesh
   do d=1,Nchi1
    NodesChi1(2*(d-1)+1)=(d-1)*pi25DT/Real(Nchi1,DP)
    NodesChi1(2*d)=d*pi25DT/Real(Nchi1,DP)
   end do
   ! now we evaluate the integrand on the initial mesh
   FuncChi1=0
   do d=1,Nchi1
   ! each d gives one cell
   ! First, we evaluate the integrand on the left node of the cell
    if ((d>1) .and. (NodesChi1(2*d-1) == NodesChi1(2*d-2))) then 
     FuncChi1(3*d-2)= FuncChi1(3*d-3) ! if the node is repeating, the value has been computed already
    else 
     FuncChi1(3*d-2) = A_IntEpsilon(NodesChi1(2*d-1),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    end if
   ! This takes care of the left node in the cell... 
   !  
   !Now we evaluate the solution on the right node of the cell and at the center. 
   FuncChi1(3*d) = A_IntEpsilon(NodesChi1(2*d),ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   FuncChi1(3*d-1) = A_IntEpsilon((NodesChi1(2*d-1) + NodesChi1(2*d))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,&
              ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
   ! repeat for all cells... 
   end do      
   ! finally, we mark all cells for the refinement in the beginning... 
   CellsChiRef1=0 
   CellsChiRef1(1:Nchi1)=1
   ! 
   !Next we adaptively refine and evaluate the integral:    
   Atemp=0
   ChiAdaptFlag =.true.
   int1 = 100*ErrChi ! initialize with something big for the intial iteration. 
   do while (ChiAdaptFlag)
   ChiAdaptFlag =.false.
   ! Now we will go over the array of the integrand values and will cross 
   ! out cells where the integrand is zero, and those that are not marked for the refinment 
   ! on those cells we summ the integrand and keep the sum.
   ! Cells marked for the refinement we divide into halves 
   NodesChi2=0
   FuncChi2=0
   CellsChiRef2=0
   e=0
   if (int1>ErrChi) then  ! this checks that the integral over the mesh NodesChi1 will be small anyway
   int1=0 ! this will keep the estimate of the integral on the mesh NodesChi1 so that we stop resolving when this number is small
   do d=1,Nchi1
   if ((abs(FuncChi1(3*d-2))+abs(FuncChi1(3*d-1))+abs(FuncChi1(3*d)) > 1.0d-15) .and. (CellsChiRef1(d)==1)) then 
    if (e > NdsChi-2) then 
     print *,"SetAKorobov_DGV: Number of nodes in chi is too big (e>NdsChi-2)"
     stop
    end if
    ! save the nodes and the integrand values
    NodesChi2(2*e+1)=NodesChi1(2*d-1)
    FuncChi2(3*e+1)=FuncChi1(3*d-2)
    NodesChi2(2*e+4)=NodesChi1(2*d)
    FuncChi2(3*e+6)=FuncChi1(3*d)
    ! Evaluate the midpoint node: 
    NodesChi2(2*e+2)=(NodesChi1(2*d-1)+NodesChi1(2*d))/2
    NodesChi2(2*e+3)=NodesChi2(2*e+2)
    ! Save the midpoint value of the integrand:
    FuncChi2(3*e+3)=FuncChi1(3*d-1)
    FuncChi2(3*e+4)=FuncChi2(3*e+3)
    ! Evaluate the integrand on the new nodes
    FuncChi2(3*e+2) = A_IntEpsilon((NodesChi2(2*e+2)+NodesChi2(2*e+1))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    FuncChi2(3*e+5) = A_IntEpsilon((NodesChi2(2*e+4)+NodesChi2(2*e+3))/2.0_DP,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,&
               ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1)
    ! Now it is time to compare the quadratures and to mark cells for the refinement 
    quad  = (FuncChi1(3*d-2)+4*FuncChi1(3*d-1)+FuncChi1(3*d))*(NodesChi1(2*d)-NodesChi1(2*d-1))/6
    quad1 = (FuncChi2(3*e+1)+4*FuncChi2(3*e+2)+FuncChi2(3*e+3))*(NodesChi2(2*e+2)-NodesChi2(2*e+1))/6
    quad2 = (FuncChi2(3*e+4)+4*FuncChi2(3*e+5)+FuncChi2(3*e+6))*(NodesChi2(2*e+4)-NodesChi2(2*e+3))/6
    es = abs(quad1+quad2-quad)/15/(NodesChi1(2*d)-NodesChi1(2*d-1)) ! error indicator
    int1=int1+abs(quad1)+abs(quad2) ! this quantity estimates the integral over the entire mesh NodesChi2 that will be NodesChi1 soon.. 
    if (es > ErrChi) then 
     ChiAdaptFlag = .true.
     CellsChiRef2(e+1) = 1
     CellsChiRef2(e+2) = 1
    else
     CellsChiRef2(e+1) = 0
     CellsChiRef2(e+2) = 0
     Atemp = Atemp + quad1 + quad2 
    end if
    e=e+2
   ! end refinement 
   end if 
   end do 
   !
   end if ! end of the check that the integral over the mesh NodesChi1 will be small anyway
   ! Finally, we need to replace the first mesh with the refined mesh 
   FuncChi1=0
   FuncChi1=FuncChi2
   NodesChi1=0
   NodesChi1=NodesChi2
   CellsChiRef1=0
   CellsChiRef1=CellsChiRef2
   Nchi1=e
   end do
   !!!!!!!!!!!!!!!!!! Looks like we are done with the integration in chi !!!!! 
   ! now it is time to check if the value of A is non-zero. If it is bigger than some specified level of 
   ! minimum sensitivity, then is it added to storage, Otherwise we continue to the next velocity #2.
   if ( ABS(Atemp) > min_sens ) then     
!$omp critical    
     A_ct = A_ct+1               ! we count the record for A
     phi_ct=phi_ct + 1           ! we count the record for this basis function. 
     if (An < A_ct) then
        call ExtendAKorArraysDGV(An,ni)
     end if  
      ! now we need to add the volume elements for xi and xi1 
      ! ATTENTION: molecular diameter mol_diam is removed in the dimensionless formulation)
      ! Dimensionless code. Molecular diameter is accounted for in the spatial operator. 
      Atemp = Atemp*( ((u_R-u_L)*(v_R-v_L)*(w_R-w_L))**2 )/Real(korob_net_param(1),DP) ! this takes care of the volume elements and the weight 1/p for xi,xi1 
      !		
      dsph=dsph/8.0_DP  ! this takes care of the fraction 1/8 still need to add |g|
      Akor(A_ct) = Atemp*dsph ! this takes care of |g|/8
      Akor_k(A_ct) = i2
      Akor_phi(A_ct) = i1
!!
iiii = omp_get_thread_num()
print *, "SetAKorobov_DGV: I2=", i2, "Thread", iiii, "A_ct=", A_ct       
!$omp end critical           
   end if                    
!!!!!!!!!!!! end evaluation of Akor !!!!!!!!!!!
end do ! END LOOP IN I2
Akor_capphi(i1)=phi_ct
print *, "SetAKorobov_DGV: i1=", i1, "A_ct=", A_ct
!! end do !! TEMPORARY do it for just one node...  
!!!!!!!!!!!!
call ShrinkAKorArraysDGV(A_ct)
!!!!!!!!!!!!
end subroutine SetAKorobov_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendAKorArraysDGV(An,ni)
!
! This subroutine extends arrays Akor, Akor_phi,Akor_k for additional ni records

subroutine ExtendAKorArraysDGV(An,ni)

use DGV_commvar, only: Akor,Akor_phi,Akor_k

integer (I4B), intent (out) :: An ! An is the length of the arrays that will be updated.  
integer (I4B), intent (in) :: ni ! ni is the number of records to add

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B) , dimension (:), allocatable :: A_iscr  ! integer scrap array
integer (I4B) :: nn ! scrap variable
!
integer (I4B) :: loc_alloc_stat
!
nn=size(Akor,1)
An=nn+ni
! extending A ... 
allocate (A_rscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGV: Allocation error for variable  A_rscr, nn=", nn
  stop
  end if
A_rscr=Akor
deallocate(Akor)
allocate (Akor(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGV: Allocation error for variable  Akor, nn=", nn+ni
  stop
  end if
Akor(1:nn)=A_rscr
deallocate(A_rscr)
! end extending Akor

! extending Akor_k
allocate(A_iscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGV: Allocation error for variable  A_iscr, nn=", nn
  stop
  end if
A_iscr=Akor_k
deallocate(Akor_k)
allocate (Akor_k(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGV: Allocation error for variable  Akor_k, nn=", nn+ni
  stop
  end if
Akor_k(1:nn)=A_iscr
!extending Akor_phi
A_iscr=Akor_phi
deallocate(Akor_phi)
allocate (Akor_phi(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGV: Allocation error for variable  Akor_phi, nn=", nn+ni
  stop
  end if
Akor_phi(1:nn)=A_iscr
! end extending A_phi
deallocate(A_iscr)

end subroutine ExtendAKorArraysDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ExtendAKorArraysDGVII(An,ni)
!
! This subroutine extends arrays AkorII, Akor_phiII, Akor_kII for additional ni records

subroutine ExtendAKorArraysDGVII(An,ni)

use DGV_commvar, only: AkorII,Akor_phiII,Akor_kII

integer (I4B), intent (out) :: An ! An is the length of the arrays that will be updated.  
integer (I4B), intent (in) :: ni ! ni is the number of records to add

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B) , dimension (:), allocatable :: A_iscr  ! integer scrap array
integer (I4B) :: nn ! scrap variable
!
integer (I4B) :: loc_alloc_stat
!
nn=size(AkorII,1)
An=nn+ni
! extending A ... 
allocate (A_rscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGVII: Allocation error for variable  A_rscr, nn=", nn
  stop
  end if
A_rscr=AkorII
deallocate(AkorII)
allocate (AkorII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGVII: Allocation error for variable  AkorII, nn=", nn+ni
  stop
  end if
AkorII(1:nn)=A_rscr
deallocate(A_rscr)
! end extending Akor

! extending Akor_k
allocate(A_iscr(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGVII: Allocation error for variable  A_iscr, nn=", nn
  stop
  end if
A_iscr=Akor_kII
deallocate(Akor_kII)
allocate (Akor_kII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGVII: Allocation error for variable  Akor_kII, nn=", nn+ni
  stop
  end if
Akor_kII(1:nn)=A_iscr
!extending Akor_phi
A_iscr=Akor_phiII
deallocate(Akor_phiII)
allocate (Akor_phiII(1:nn+ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ExtendAKorArraysDGVII: Allocation error for variable  Akor_phiII, nn=", nn+ni
  stop
  end if
Akor_phiII(1:nn)=A_iscr
! end extending A_phi
deallocate(A_iscr)

end subroutine ExtendAKorArraysDGVII

subroutine ShrinkAKorArraysDGV(ni)

use DGV_commvar, only: Akor,Akor_phi,Akor_k

integer (I4B), intent (in) :: ni ! ni is the new size of the arrays to be shrunk

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B), dimension (:), allocatable :: A_iscr  ! integer scrap array
!
integer (I4B) :: loc_alloc_stat 

! shrinking Akor ... 
allocate (A_rscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGV: Allocation error for variable  A_rscr, nn=", ni
  stop
  end if
A_rscr=Akor(1:ni)
deallocate(Akor)
allocate (Akor(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGV: Allocation error for variable  Akor, nn=", ni
  stop
  end if
Akor=A_rscr
deallocate(A_rscr)
! end shrinking Akor

! shrinking Akor_k
allocate (A_iscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGV: Allocation error for variable  A_iscr, nn=", ni
  stop
  end if
A_iscr=Akor_k(1:ni)
deallocate(Akor_k)
allocate (Akor_k(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGV: Allocation error for variable  Akor_k, nn=", ni
  stop
  end if
Akor_k=A_iscr
!shrinking Akor_phi
A_iscr=Akor_phi(1:ni)
deallocate(Akor_phi)
allocate (Akor_phi(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGV: Allocation error for variable  Akor_phi, nn=", ni
  stop
  end if
Akor_phi=A_iscr
! end shrinking Akor_phi
deallocate(A_iscr)
end subroutine ShrinkAKorArraysDGV


subroutine ShrinkAKorArraysDGVII(ni)

use DGV_commvar, only: AkorII,Akor_phiII,Akor_kII

integer (I4B), intent (in) :: ni ! ni is the new size of the arrays to be shrunk

!!!
real (DP), dimension (:), allocatable :: A_rscr ! scrap array
integer (I4B), dimension (:), allocatable :: A_iscr  ! integer scrap array
!
integer (I4B) :: loc_alloc_stat 

! shrinking Akor ... 
allocate (A_rscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGVII: Allocation error for variable  A_rscr, nn=", ni
  stop
  end if
A_rscr=AkorII(1:ni)
deallocate(AkorII)
allocate (AkorII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGVII: Allocation error for variable  AkorII, nn=", ni
  stop
  end if
AkorII=A_rscr
deallocate(A_rscr)
! end shrinking AkorII

! shrinking Akor_kII
allocate (A_iscr(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGVII: Allocation error for variable  A_iscr, nn=", ni
  stop
  end if
A_iscr=Akor_kII(1:ni)
deallocate(Akor_kII)
allocate (Akor_kII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGVII: Allocation error for variable  Akor_kII, nn=", ni
  stop
  end if
Akor_kII=A_iscr
!shrinking Akor_phiII
A_iscr=Akor_phiII(1:ni)
deallocate(Akor_phiII)
allocate (Akor_phiII(1:ni), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ShrinkAKorArraysDGVII: Allocation error for variable  Akor_phiII, nn=", ni
  stop
  end if
Akor_phiII=A_iscr
! end shrinking Akor_phiII
deallocate(A_iscr)
end subroutine ShrinkAKorArraysDGVII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! IphiDGV
!
! This subroutine evaluates the moment of the collision integral for basis functions. 
! The basis functions are defined for each velocity node. The Collision Information Operator
! is calculated for all (or few basis functions). 
! The subroutine will take calculate the moment for each basis function recorded in A
!
! Depends on the main program. A-arrays must be defined as well as the disftribution function f
! 
!!!!!!!!!!!!!!!!!!!!!!!!

function IphiDGV(f) result (Iphi)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,nodes_gwts,mol_diam,L_inf,N_inf

real (DP), dimension(:) :: f ! the main variable -- distribution function -- one component for each velocity node
real (DP), dimension (size(A_capphi,1)) :: Iphi ! the result of integration
real (DP), dimension(:), allocatable :: A_SP, F_SP
!!! 
integer (I4B) :: i,j, A_ct,count ! scrap variables,
!!!
allocate (A_SP(size(A,1)), F_SP(size(f,1)))
A_SP = A
f_SP = f
!!!
count=0
Iphi=0
A_ct=0
do i=1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   Iphi(i)=0
   do j=1+A_ct,A_ct+A_capphi(i)
   Iphi(i)=Iphi(i)+A_SP(j)*f_SP(A_xi(j))*f_SP(A_xi1(j))
   if (ABS(A_SP(j)*f_SP(A_xi(j))*f_SP(A_xi1(j)))<1.0d-8) then 
   count=count+1
   end if 
   end do 
   Iphi(i)=2*Iphi(i)/nodes_gwts(i)*((mol_diam/L_inf)**2*N_inf)
   A_ct=A_ct+A_capphi(i)
end do 
!
deallocate (A_SP, f_SP)
end function IphiDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! IphiDGV_decomp
!
! This is a diagnistics subroutine 
!
! This subroutine evaluates the moment of the collision integral for basis functions. 
! using the decomposition method. It is to test the advantage of decomposition -- if any 
! the solution is split into a maxwellian adn the rest. The collision operator is evcaluated.
! The basis functions are defined for each velocity node. The Collision Information Operator
! is calculated for all (or few basis functions). 
! The subroutine will take calculate the moment for each basis function recorded in A
!
! Depends on the main program. A-arrays must be defined as well as the disftribution function f
! 
!!!!!!!!!!!!!!!!!!!!!!!!

function IphiDGV_decomp(f) result (Iphi)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,nodes_gwts,mol_diam,L_inf,N_inf,nodes_u, nodes_v, nodes_w

use DGV_distributions_mod



real (DP), dimension(:) :: f ! the main variable -- distribution function -- one component for each velocity node
real (DP), dimension (size(A_capphi,1)) :: Iphi ! the result of integration

!!! 
integer (I4B) :: i,j, A_ct ! scrap variables,
real (DP), dimension(size(f,1)) :: fmaxwellNew,Df
real (DP) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr 
!!!!!!!
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
Df=f-fMaxwellNew ! evaluate the perturbation from the maxwellian. 

Iphi=0
A_ct=0
do i=1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   Iphi(i)=0
   do j=1+A_ct,A_ct+A_capphi(i)
   Iphi(i)=Iphi(i) + A(j)*fMaxwellNew(A_xi(j))*Df(A_xi1(j))+A(j)*Df(A_xi(j))*fMaxwellNew(A_xi1(j)) + A(j)*Df(A_xi(j))*Df(A_xi1(j))
   end do 
   Iphi(i)=2*Iphi(i)/nodes_gwts(i)*((mol_diam/L_inf)**2*N_inf)
   A_ct=A_ct+A_capphi(i)
end do 
!
end function IphiDGV_decomp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!
!EvalPhiPostCollIntegrand
!
! This function helps to calculate the post collision velocities and calls evaluation of the basis function
! on the post collision velocitis. The function is mainly to help coding -- to reduce the number of repeating lines. 
!
!

function EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
             ugux3,ugvx3,ugwx3,epsil,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) result (y)
!
!
real (DP), intent (in) :: sinchi, coschi, g1,g2                           ! some useful numbers 
real (DP), intent (in) :: ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! useful coefficients
real (DP), intent (in) :: epsil ! the angle for which the integral must be evaluated.
real (DP), intent (in) :: dsph,xiu,xiv,xiw,xi1u,xi1v,xi1w ! the pre-collision velocities 
integer (I4B), intent (in) :: i1 ! the number of the basis fucntion to evaluate
real (DP), intent (in) :: varphi_xi,varphi_xi1 ! the values of the basis function on xi and xi1
real (DP) :: y ! the result

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: sinchicoseps, sinchisineps ! useful quantities
real (DP) :: xiupr,xivpr,xiwpr,xi1upr,xi1vpr,xi1wpr ! components of post-collision velocity
real (DP) :: z,ku,kv,kw ! compoents of the unit vector in the direction of post collision relative speed
!!!!!!
     sinchicoseps = cos(epsil)*sinchi/g1
     sinchisineps = sin(epsil)*sinchi/g2
   ! Now it is time to evaluate post-collision velocities... 
   ! First we evaluate the direction of $g'$ : 
     ku = (ugu)*coschi - ugux2*sinchicoseps - ugux3*sinchisineps
     kv = (ugv)*coschi - ugvx2*sinchicoseps - ugvx3*sinchisineps
     kw = (ugw)*coschi - ugwx2*sinchicoseps - ugwx3*sinchisineps
   ! Now we calculate the post collision velocities: 
   !
     xiupr=((xiu + xi1u) + dsph*ku)/2.0_DP
     xi1upr=((xiu + xi1u) - dsph*ku)/2.0_DP
     xivpr=((xiv + xi1v) + dsph*kv)/2.0_DP
     xi1vpr=((xiv + xi1v) - dsph*kv)/2.0_DP
     xiwpr=((xiw + xi1w) + dsph*kw)/2.0_DP
     xi1wpr=((xiw + xi1w) - dsph*kw)/2.0_DP

   ! now we evaluate basis the basis functions on these velocities 
     z = EvalLagrBasisFunByNdsDGblzmEZ(xiupr,xivpr,xiwpr,i1) + &
       EvalLagrBasisFunByNdsDGblzmEZ(xi1upr,xi1vpr,xi1wpr,i1) !- varphi_xi - varphi_xi1 ! comment the last two terms if 
                                                              ! want to compute only positive terms in operator A. (the last two correspond to the values of the basis function on Pre-collision velocities)
     y=z  
end function EvalPhiPostCollIntegrand     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A_IntEpsilon
! This function is to shorten the code. It is pretty much a piece of code that is cut out
! and paste in the form of a function to make the rest of the code read easier. 
! This portion implements the integration in \varepsilon in the evaluation of operator A
!
! This code evaluates \sin\chi \int_{0}^{2\pi} (\varphi(\xi')+\varphi(\xi'_{1}) d\varepsilon 
!
!
! The result is the integral for this particular angle \chi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function A_IntEpsilon(chi,ires_e,xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,g1,g2,dsph,i1,ErrEps,varphi_xi,varphi_xi1) result (y)
             
real (DP), intent (in) :: chi,g1,g2,ires_e                    ! some useful numbers 
real (DP), intent (in) :: ugu,ugv,ugw,ugux2,ugvx2,ugwx2,ugux3,ugvx3,ugwx3 ! useful coefficients
real (DP), intent (in) :: dsph,xiu,xiv,xiw,xi1u,xi1v,xi1w ! the pre-collision velocities 
integer (I4B) :: i1 ! the number of the basis function which is to evaluate. 
real (DP), intent (in) :: ErrEps ! the max set error of integral evauation. 
real (DP), intent (in) :: varphi_xi,varphi_xi1 ! the values of the basis function on xi and xi1
!
real (DP) :: y ! the result
             

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: sinchi,coschi ! useful scrap constants to remember..
real (DP), parameter :: pi25DT = 3.141592653589793238462643d0 
integer (I4B), parameter :: NdsEps=10000 ! The max number of cells for integration in chi and epsilon

real (DP), dimension (2*NdsEps) :: NodesEps1, NodesEps2
real (DP), dimension (3*NdsEps) :: FuncEps1, FuncEps2 ! Arrays to store cells and functin values on cells Two copies are needed. Integration in Eps
integer (I4B), dimension (NdsEps) :: CellsEpsRef1,CellsEpsRef2 ! Arrays to keep refinement flags
real (DP) :: Atemp_temp,quad,quad1,quad2,es,int1 ! scrap variables.
!!!!!
integer (I4B) :: Neps1,g ! variables to keep the number of integration cells
integer (I4B) :: f ! local counter
logical :: EpsAdaptFlag ! a logical variable to tell if any of the refinements need to be done. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   sinchi=sin(chi) ! will be useful to remember this one... 
   coschi=cos(chi) ! will be useful to remember this one too... 
! Begin Integration in epsilon ...  
   Neps1=FLOOR(pi25DT*dsph*sinchi/ires_e)+1
   if (Neps1>NdsEps) then 
     Neps1=NdsEps
     print *,"SetA_DGV: Warning! Number of nodes (NdsEps) for integration in epsilon gives insufficient resolution"  
   end if 
   ! we set the initial mesh
   do f=1,Neps1
    NodesEps1(2*(f-1)+1)=(f-1)*2*pi25DT/Real(Neps1,DP)
    NodesEps1(2*f)=f*2*pi25DT/Real(Neps1,DP)
   end do
   ! now we evaluate the integrand on the initial mesh
   FuncEps1=0
   do f=1,Neps1
   ! each f gives one cell
   ! First, we evaluate the integrand on the left node of the cell
    if (f>1) then 
     if (NodesEps1(2*f-1) == NodesEps1(2*f-2)) then 
      FuncEps1(3*f-2)= FuncEps1(3*f-3) ! if the node is repeating, the value has been computed already
     else 
      FuncEps1(3*f-2) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
                ugux3,ugvx3,ugwx3,NodesEps1(2*f-1),sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1)
     end if 
     FuncEps1(3*f-2) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
                ugux3,ugvx3,ugwx3,NodesEps1(2*f-1),sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1)
    end if
   ! This takes care of the left node in the cell... 
   !  
   !Now we evaluate the solution on the right node of the cell and at the center. 
   FuncEps1(3*f) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,NodesEps1(2*f),sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
   FuncEps1(3*f-1) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
              ugux3,ugvx3,ugwx3,(NodesEps1(2*f)+NodesEps1(2*f-1))/2,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1)
   ! repeat for all cells... 
   end do      
   ! finally, we mark all cells for the refinement in the beginning... 
   CellsEpsRef1=0 
   CellsEpsRef1(1:Neps1)=1
   ! 
   !Next we adaptively refine and evaluate the integral:    
   int1=10*ErrEps !! set it to a large number fo the rist run
   EpsAdaptFlag =.true.
   Atemp_temp=0
   do while (EpsAdaptFlag)
   EpsAdaptFlag =.false.
   ! Now we will go over the array of the integrand values and will cross 
   ! out cells where the integrand is zero, and those that are not marked for the refinment 
   ! on those cells we summ the integrand and keep the sum.
   ! Cells marked for the refinement we divide into halves 
   NodesEps2=0
   FuncEps2=0
   CellsEpsRef2=0
   g=0
   if (int1 > ErrEps) then   ! check is the integral over NodesEps1 is expected to be small. If it is small then the refinement is not performed
   int1=0 ! this will keep the estimate of the integral on the mesh NodesEps1 so that we stop resolving when this number is small
   do f=1,Neps1
   if ((abs(FuncEps1(3*f-2))+abs(FuncEps1(3*f-1))+abs(FuncEps1(3*f)) > 1.0d-15) .and. (CellsEpsRef1(f)==1)) then 
    if (g > NdsEps-2) then 
     print *,"SetA_DGV: Number of nodes in epsilon is too big (g>NdsEps-2)"
     stop
    end if
    ! save the nodes and the integrand values
    NodesEps2(2*g+1)=NodesEps1(2*f-1)
    FuncEps2(3*g+1)=FuncEps1(3*f-2)
    NodesEps2(2*g+4)=NodesEps1(2*f)
    FuncEps2(3*g+6)=FuncEps1(3*f)
    ! Evaluate the midpoint node: 
    NodesEps2(2*g+2)=(NodesEps1(2*f-1)+NodesEps1(2*f))/2
    NodesEps2(2*g+3)=NodesEps2(2*g+2)
    ! Save the midpoint value of the integrand:
    FuncEps2(3*g+3)=FuncEps1(3*f-1)
    FuncEps2(3*g+4)=FuncEps2(3*g+3)
    ! Evaluate the integrand on the new nodes
    FuncEps2(3*g+2) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
    ugux3,ugvx3,ugwx3,(NodesEps2(2*g+2)+NodesEps2(2*g+1))/2.0_DP,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
    FuncEps2(3*g+5) = EvalPhiPostCollIntegrand(xiu,xiv,xiw,xi1u,xi1v,xi1w,ugu,ugv,ugw,ugux2,ugvx2,ugwx2,&
    ugux3,ugvx3,ugwx3,(NodesEps2(2*g+4)+NodesEps2(2*g+3))/2.0_DP,sinchi,coschi,g1,g2,dsph,i1,varphi_xi,varphi_xi1) 
    ! Now it is time to compare the quadratures and to mark cells for the refinement 
    quad  = (FuncEps1(3*f-2)+4*FuncEps1(3*f-1)+FuncEps1(3*f))*(NodesEps1(2*f)-NodesEps1(2*f-1))/6
    quad1 = (FuncEps2(3*g+1)+4*FuncEps2(3*g+2)+FuncEps2(3*g+3))*(NodesEps2(2*g+2)-NodesEps2(2*g+1))/6
    quad2 = (FuncEps2(3*g+4)+4*FuncEps2(3*g+5)+FuncEps2(3*g+6))*(NodesEps2(2*g+4)-NodesEps2(2*g+3))/6
    es = abs(quad1+quad2-quad)/15/(NodesEps1(2*f)-NodesEps1(2*f-1)) ! local error indicator for simpson's rule.
    int1 = int1 + abs(quad1) + abs(quad2) ! This will calculate the estimate from above to the total integral on the current mesh 
    if (es > ErrEps) then 
     EpsAdaptFlag = .true.
     CellsEpsRef2(g+1) = 1
     CellsEpsRef2(g+2) = 1
    else
     CellsEpsRef2(g+1) = 0
     CellsEpsRef2(g+2) = 0
     Atemp_temp = Atemp_temp + quad1 + quad2 
    end if                   
    g=g+2
   ! end refinement 
   end if 
   end do 
   !
   end if ! end check if the integral over NodesEps1 is too small to worry about it. 
   ! Finally, we need to replace the first mesh with the refined mesh 
   FuncEps1=0
   FuncEps1=FuncEps2
   NodesEps1=0
   NodesEps1=NodesEps2
   CellsEpsRef1=0
   CellsEpsRef1=CellsEpsRef2
   Neps1=g
   end do
   !!!!!!!!!!!!!!!!!! Looks like we are done with the integration in epsilon !!!!! 
   y=Atemp_temp*sinchi
end function A_IntEpsilon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheck
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! 
!!!!!!!!!!!!!!!!!
function MassCheck (f) result (y)
use DGV_commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP) :: y ! the result (the mass) 
!!!!!!!!!!!!!!!!!!
real (DP) :: n,ubar,vbar,wbar,temp ! number density, av_v
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!
y=ubar ! return this... 

end function MassCheck


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FindCellContainsPoint_DGV
!
! This subroutine find the number of the active cell in the cell arrays that containes
! the point (u,v,w)
!
! the function returns zero if the velocity point is not on any cells. 
!
!
! the subroutine accesses arrays
! nodes_pcell, nodes_u, nodes_v, nodes_w,&
!				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
!                  cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
!                  cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
!                  grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
!                  nodes_ui, nodes_vi, nodes_wi
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FindCellContainsPoint_DGV(u,v,w) result (y)

use DGV_commvar, only: nodes_pcell, nodes_u, nodes_v, nodes_w,&
				   cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov, grids_cap_u, &
                   grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,g_nds_all, &
                   nodes_ui, nodes_vi, nodes_wi

real (DP) :: u,v,w ! the components of the given velocity point
integer (I4B) ::  y ! the value of the index in the cell arrays that correcponds to the active cell containing point (u,v,w)

!!!!!!!!!!!!!!                      
integer (I4B) :: j ! a counter -- usually denotes the number of the grid we are working on. 
integer (I4B) :: gi_v,gi_u,gi_w ! counters to help count nodes in the grids array
integer (I4B) :: celli, cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
integer (I4B) :: ui,vi,wi ! are the local numbers on the particular grids where u belongs. Because we have ierarchicval grids, 
             ! there may be different grids that have ui,vi,wi defined. If one of these numbers is zero -- then the 
             ! velocity u,v,w is not on the grid. Unless a mistake was done in setting the grids either all three of them 
             ! are zero -- that means that the velocity is not on this (level zero grid) or all three are non-zero  -- means that
             ! the velocity has been found. Mized results may indicate a breach in data integrity. 
             !!!!!
integer (I4B) :: jj ! a counter
real (DP) :: unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
! we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
! if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
!!!!!!!!
!
j=1
!! set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0
gi_v=0
gi_w=0
do while (j<=size(grids_cap_u,1))
  call finduiviwi(grids_u(gi_u+1:gi_u+grids_cap_u(j)),grids_v(gi_v+1:gi_v+grids_cap_v(j)),&
                       grids_w(gi_w+1:gi_w+grids_cap_w(j)),u,v,w,ui,vi,wi)
  if (ui*vi*wi /= 0) then 
     ! now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1
     do jj=1,j-1
      celli = celli + (grids_cap_u(jj)-1)*(grids_cap_v(jj)-1)*(grids_cap_w(jj)-1)
     end do
     celli=celli + (ui-2)*(grids_cap_v(j)-1)*(grids_cap_w(j)-1)+(vi-2)*(grids_cap_w(j)-1)+(wi-2)
     !  check if this cell is refined
     if ((cells_refu(celli)>1) .or. (cells_refv(celli)>1) .or. (cells_refw(celli)>1)) then  
        j=cells_cgrid(celli)
        ! set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0
        gi_v=0
        gi_w=0
        do jj=1,j-1
         gi_u = gi_u + grids_cap_u(jj)
         gi_v = gi_v + grids_cap_v(jj)
         gi_w = gi_w + grids_cap_w(jj)
        end do
     else 
        exit
     endif
  else 
     celli=0
     exit
  endif                      
enddo
! now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)
! we return this number 
y=celli
! 
end function FindCellContainsPoint_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FindCellContainsPoint_DGVII
!
! This is a copy of the above subroutine except it works with the secondary mesh
!
! This subroutine find the number of the active cell in the cell arrays that containes
! the point (u,v,w)
!
! the function returns zero if the velocity point is not on any cells. 
!
!
! the subroutine accesses arrays
! nodes_pcellII, nodes_uII, nodes_vII, nodes_wII,&
!				   cells_pgridII, cells_cgridII, cells_luII, cells_lvII, cells_lwII, & 
!                  cells_ruII, cells_rvII, cells_rwII, cells_refuII, cells_refvII, & 
!                  cells_refwII, cells_gowII, cells_gouII, cells_govII, grids_cap_uII, &
!                  grids_cap_vII,grids_cap_wII,grids_uII,grids_vII,grids_wII,g_nds_all, &
!                  nodes_uiII, nodes_viII, nodes_wiII
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FindCellContainsPoint_DGVII(u,v,w) result (y)

use DGV_commvar, only: nodes_pcellII, nodes_uII, nodes_vII, nodes_wII,&
				   cells_pgridII, cells_cgridII, cells_luII, cells_lvII, cells_lwII, & 
                   cells_ruII, cells_rvII, cells_rwII, cells_refuII, cells_refvII, & 
                   cells_refwII, cells_gowII, cells_gouII, cells_govII, grids_cap_uII, &
                   grids_cap_vII,grids_cap_wII,grids_uII,grids_vII,grids_wII,g_nds_all, &
                   nodes_uiII, nodes_viII, nodes_wiII

real (DP), intent (in) :: u,v,w ! the components of the given velocity point
integer (I4B) ::  y ! the value of the index in the cell arrays that correcponds to the active cell containing point (u,v,w)

!!!!!!!!!!!!!!                      
integer (I4B) :: j ! a counter -- usually denotes the number of the grid we are working on. 
integer (I4B) :: gi_v,gi_u,gi_w ! counters to help count nodes in the grids array
integer (I4B) :: celli, cellp ! these will store the number of the cell where the velocity belongs and the number of the 
            ! parent cell where the basis function belongs. 
integer (I4B) :: ui,vi,wi ! are the local numbers on the particular grids where u belongs. Because we have ierarchicval grids, 
             ! there may be different grids that have ui,vi,wi defined. If one of these numbers is zero -- then the 
             ! velocity u,v,w is not on the grid. Unless a mistake was done in setting the grids either all three of them 
             ! are zero -- that means that the velocity is not on this (level zero grid) or all three are non-zero  -- means that
             ! the velocity has been found. Mized results may indicate a breach in data integrity. 
             !!!!!
integer (I4B) :: jj ! a counter
real (DP) :: unor,vnor,wnor ! scrap variables to store the normalized velus of the components u,v,w             


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
! we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
! if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
!!!!!!!!
!
j=1
!! set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0
gi_v=0
gi_w=0
do while (j<=size(grids_cap_uII,1))
  call finduiviwi(grids_uII(gi_u+1:gi_u+grids_cap_uII(j)),grids_vII(gi_v+1:gi_v+grids_cap_vII(j)),&
                       grids_wII(gi_w+1:gi_w+grids_cap_wII(j)),u,v,w,ui,vi,wi)
  if (ui*vi*wi /= 0) then 
     ! now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1
     do jj=1,j-1
      celli = celli + (grids_cap_uII(jj)-1)*(grids_cap_vII(jj)-1)*(grids_cap_wII(jj)-1)
     end do
     celli=celli + (ui-2)*(grids_cap_vII(j)-1)*(grids_cap_wII(j)-1)+(vi-2)*(grids_cap_wII(j)-1)+(wi-2)
     !  check if this cell is refined
     if ((cells_refuII(celli)>1) .or. (cells_refvII(celli)>1) .or. (cells_refwII(celli)>1)) then  
        j=cells_cgridII(celli)
        ! set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0
        gi_v=0
        gi_w=0
        do jj=1,j-1
         gi_u = gi_u + grids_cap_uII(jj)
         gi_v = gi_v + grids_cap_vII(jj)
         gi_w = gi_w + grids_cap_wII(jj)
        end do
     else 
        exit
     endif
  else 
     celli=0
     exit
  endif                      
enddo
! now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)
! we return this number 
y=celli
! 
end function FindCellContainsPoint_DGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FindI1sByCellNum_DGV
!
! This functions looks at the nodes arrays and finds the numbers of the nodes that
! belong to cell with the number i
! 
! The numbers will be recorded in the array of results with the first number indicating the number of useful records/
! 
! For example: if the result array has a zero in its first place, then the cell has no nodes -- which is not expected, really... 
! 
! if the array has number 5 in its firs place, than elements 2--6 contain the I1s --- these numbers will be used to build the A-array.
! 
!  if there is no cell with the number i -- the program retrurns zero records and prints a worning
!  if the result array is too short to fit all the numbers i1, the program prints the error message and stops.
!
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FindI1sByCellNum_DGV(i,y)

use DGV_commvar, only: nodes_pcell,cells_lu
 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), intent (in) :: i ! the number of the cell where I1s need to be looked at. 
integer (I4B), dimension (:), intent (out) :: y ! the numbers of the nodes (I1 -- in our sleng...) that belong to the cell with number i 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (i4B) :: j,k,pcell_i,sizey ! local counters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sizey=size(y,1)
y=0 ! nullify the thing... Just In Case 
k=0 ! in the beginning there is no records///
if (i>size(cells_lu,1)) then  ! check if i is too big
 y=0 !
 print *,"FindI1sByCellNum_DGV: the cell with the provided number i do not exist -- i is too big. Returned zero records"
else ! if i is not too big, 
 do j=1,size(nodes_pcell,1)
  pcell_i=nodes_pcell(j) 
  if (pcell_i==i) then ! check if the node belongs to cell i 
   k=k+1                ! if it does, check if can still records in y
   if (k+1<=sizey) then  
    y(1+k)=j     !record
   else  ! othersize print the error message and stop
    print *,"FindI1sByCellNum_DGV: The size fo the result array is too small. Can not put all I1s. Stop"
    stop
   end if
  end if  
 end do 
y(1)=k ! the first records is reserved for the number of found I1s
end if 
end subroutine FindI1sByCellNum_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsTrhld_DGV (trhld)
!
! This subroutine trancates the A-arrays. All entries of A that are below the provided 
! treshhold value (trhld) are being cut, arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! ATTENTION: array nodes_Ashift need to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsTrhls_DGV (trhld)

use DGV_commvar, only: A, A_xi, A_xi1, A_phi, A_capphi 

!!!!
real (DP), intent (in) :: trhld ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
!
integer (I4B) :: loc_alloc_stat 

nold = size(A,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=A; A_xiscr=A_xi; A_xi1scr=A_xi1; A_phiscr=A_phi
! Now clear the A-arrays;  
A=0.0_DP; A_xi=0; A_xi1=0; A_phi=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
do i = 1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphi(i)
    if (ABS(Ascr(j))>trhld) then  ! if records in A array are above treshhold, they are saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     A(nnew)=Ascr(j)
     A_xi(nnew)=A_xiscr(j)
     A_xi1(nnew)=A_xi1scr(j)
     A_phi(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphi(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphi(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGV(nnew) !dgvtools_mod.f90

end subroutine TruncAarrsTrhls_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsTrhld_DGVII (trhld)
!
! This subroutine trancates the A-arrays. All entries of A that are below the provided 
! treshhold value (trhld) are being cut, arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! ATTENTION: array nodes_Ashift need to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsTrhls_DGVII (trhld)

use DGV_commvar, only: AII, A_xiII, A_xi1II, A_phiII, A_capphiII 

!!!!
real (DP), intent (in) :: trhld ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
!
integer (I4B) :: loc_alloc_stat 

nold = size(AII,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGVII: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=AII; A_xiscr=A_xIIi; A_xi1scr=A_xi1II; A_phiscr=A_phiII
! Now clear the A-arrays;  
AII=0.0_DP; A_xiII=0; A_xi1II=0; A_phiII=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
do i = 1,size(A_capphiII,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphiII(i)
    if (ABS(Ascr(j))>trhld) then  ! if records in A array are above treshhold, they are saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     AII(nnew)=Ascr(j)
     A_xiII(nnew)=A_xiscr(j)
     A_xi1II(nnew)=A_xi1scr(j)
     A_phiII(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphiII(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphiII(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGVII(nnew) !dgvtools_mod.f90

end subroutine TruncAarrsTrhls_DGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsRadius_DGV (rad)
!
! This subroutine trancates the A-arrays based on the distance between xi and xi1. If the distance is larger than 
! the specified radus (rad), the record is discarded from A and arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! ATTENTION: array nodes_Ashift needs to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsRadius_DGV (rad)

use DGV_commvar, only: A, A_xi, A_xi1, A_phi, A_capphi, nodes_u, nodes_v, nodes_w 

!!!!
real (DP), intent (in) :: rad ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap,ixi,ixi1 ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
real (DP) :: dsq, radsq ! scrap variables to keep distance squared and radius squared
!
integer (I4B) :: loc_alloc_stat 

nold = size(A,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=A; A_xiscr=A_xi; A_xi1scr=A_xi1; A_phiscr=A_phi
! Now clear the A-arrays;  
A=0.0_DP; A_xi=0; A_xi1=0; A_phi=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
radsq=rad**2
do i = 1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphi(i)
    ixi=A_xiscr(j)
    ixi1=A_xi1scr(j)
    dsq = (nodes_u(ixi)-nodes_u(ixi1))**2 + (nodes_v(ixi)-nodes_v(ixi1))**2 + (nodes_w(ixi)-nodes_w(ixi1))**2
    if (dsq <= radsq) then  ! if vectors xi and xi1 are at or closer than distance rad, the record in A is saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     A(nnew)=Ascr(j)
     A_xi(nnew)=A_xiscr(j)
     A_xi1(nnew)=A_xi1scr(j)
     A_phi(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphi(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphi(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGV(nnew)

end subroutine TruncAarrsRadius_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsRadius_DGVII (rad)
!
! This subroutine trancates the A-arrays based on the distance between xi and xi1. If the distance is larger than 
! the specified radus (rad), the record is discarded from A and arrays Axi, Axi1, Aphi, A_capphi are updated correspondingly... 
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! ATTENTION: array nodes_Ashift needs to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsRadius_DGVII (rad)

use DGV_commvar, only: AII, A_xiII, A_xi1II, &
                       A_phiII, A_capphiII, &
                       nodes_uII, nodes_vII, nodes_wII 

!!!!
real (DP), intent (in) :: rad ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap,ixi,ixi1 ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
real (DP) :: dsq, radsq ! scrap variables to keep distance squared and radius squared
!
integer (I4B) :: loc_alloc_stat 

nold = size(AII,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=AII; A_xiscr=A_xiII; A_xi1scr=A_xi1II; A_phiscr=A_phiII
! Now clear the A-arrays;  
AII=0.0_DP; A_xiII=0; A_xi1II=0; A_phiII=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
radsq=rad**2
do i = 1,size(A_capphiII,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphiII(i)
    ixi=A_xiscr(j)
    ixi1=A_xi1scr(j)
    dsq = (nodes_uII(ixi)-nodes_uII(ixi1))**2 + (nodes_vII(ixi)-nodes_vII(ixi1))**2 + (nodes_wII(ixi)-nodes_wII(ixi1))**2
    if (dsq <= radsq) then  ! if vectors xi and xi1 are at or closer than distance rad, the record in A is saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     AII(nnew)=Ascr(j)
     A_xiII(nnew)=A_xiscr(j)
     A_xi1II(nnew)=A_xi1scr(j)
     A_phiII(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphiII(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphiII(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGVII(nnew)

end subroutine TruncAarrsRadius_DGVII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncAarrsRadiusFromOrigin_DGV (rad)
!
! This subroutine truncates the A-arrays based on the distance of xi and xi1 from the origin. If either is too far, then the
! A array is truncated
!
! ATTENTION: array nodes_Ashift needs to be re-created after A has been truncated.
!
!
!!!!!!!!!!!!!!!!!!!
subroutine TruncAarrsRadiusFromOrigin_DGV (rad)

use DGV_commvar, only: A, A_xi, A_xi1, A_phi, A_capphi, nodes_u, nodes_v, nodes_w 

!!!!
real (DP), intent (in) :: rad ! the treshhold at which to cut A:
!!!
integer (I4B) :: nold,nnew ! are the old and new sizes of the A-arrays
integer (I4B) :: i,j,phicap,ixi,ixi1 ! scrap indices
!!!
real (DP), dimension (:), allocatable :: Ascr ! real scrap array
integer (I4B), dimension (:), allocatable :: A_xiscr, A_xi1scr, A_phiscr  ! integer scrap arrays
real (DP) :: dsqi, dsqi1, radsq ! scrap variables to keep distance squared and radius squared
!
integer (I4B) :: loc_alloc_stat 

nold = size(A,1)
! shrinking A ... 
allocate (Ascr(1:nold), A_xiscr(1:nold), A_xi1scr(1:nold), A_phiscr(1:nold), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "TruncAarrsTrhls_DGV: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
  end if
! Save the A arrays .... 
Ascr=A; A_xiscr=A_xi; A_xi1scr=A_xi1; A_phiscr=A_phi
! Now clear the A-arrays;  
A=0.0_DP; A_xi=0; A_xi1=0; A_phi=0
! Now start to fill them in, but omitting all records that are below the treshhold value
nold=0 ! this will count the original records of A  === similar to A_ct index
nnew=0 ! this will count the records after the truncation. 
radsq=rad**2
do i = 1,size(A_capphi,1) ! loop in basis functions --- one node=one basis function
   phicap = 0
   do j = 1+nold,nold+A_capphi(i)
    ixi=A_xiscr(j)
    ixi1=A_xi1scr(j)
    dsqi = nodes_u(ixi)**2 + nodes_v(ixi)**2 + nodes_w(ixi)**2
    dsqi1 = nodes_u(ixi1)**2 + nodes_v(ixi1)**2 + nodes_w(ixi1)**2
    if (dsqi <= radsq .and. dsqi1 <= radsq) then  ! if vectors xi and xi1 are at or closer than distance rad, the record in A is saved, otherwise ignored
     phicap = phicap+1
     nnew = nnew+1 
     A(nnew)=Ascr(j)
     A_xi(nnew)=A_xiscr(j)
     A_xi1(nnew)=A_xi1scr(j)
     A_phi(nnew)=A_phiscr(j)
    end if 
   end do 
   nold = nold+A_capphi(i) ! shift the start index to the place where record sfo rthe next basis function start.
   A_capphi(i) = phicap ! Update the A_capphi --- it now stores the number of records after the truncation
end do 
deallocate(Ascr,A_xiscr,A_xi1scr,A_phiscr) !

call ShrinkAarraysDGV(nnew)

end subroutine TruncAarrsRadiusFromOrigin_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TruncFRadius(f, rad, ubar, vbar, wbar)
!
! This subroutine truncates the f array based on the radius of f from the bulk velocity. rad is the minimum threshold
! required for truncation. 
!
! This function works on the primary mesh
!
!!!!!!!!!!!!!!!!!!!

subroutine TruncFRadius(f, rad, ubar, vbar, wbar)

use DGV_commvar, only: nodes_u, nodes_v, nodes_w, nodes_gwts

real (DP), dimension(:), intent (inout) :: f
real (DP), intent (in) :: rad !minimum radius required to truncate
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: rad2 !square of the rad
real (DP), dimension(:), allocatable :: rads
real (DP) :: n !density
real (DP) :: ubar, vbar, wbar !bulk velocities
integer(I4B) :: i !loop
integer(I4B) :: loc_alloc_stat ! allocation check

allocate(rads(1:size(nodes_u,1)),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "TruncFRadius: Allocation error for variables  Ascr, A_xiscr, _xi1scr, _phiscr" 
  stop
end if



rad2=rad*rad
rads=((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2)
do i=1,size(nodes_u,1)
  if (rads(i) > rad2) then
    if (f(i) > 1.0D-5) then
      print *,"TruncFRadius: Truncating large value at position",i
    end if
    f(i)=0
  end if
end do


end subroutine TruncFRadius


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRec
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! 
!!!!!!!!!!!!!!!!!
subroutine MassCheckRec (f,n,ubar,vbar,wbar,temp)
use DGV_commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp ! number density, av_v
!!!!!!!!!!!!!!!!!!
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!
end subroutine MassCheckRec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRecPlus
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! This one also calculate momentum, temperature and the directional temperature 
!!!!!!!!!!!!!!!!!
subroutine MassCheckRecPlus (f,n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w)
use DGV_commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w ! number density, av_v
!!!!!!!!!!!!!!!!!!
 
integer (I4B) :: i ! scrap index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
! check directional temperatures 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_u = sum(f*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tx
temp_v = sum(f*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Ty
temp_w = sum(f*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tz
!!!!!!!!!!
end subroutine MassCheckRecPlus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MassCheckRecPlusPlus
! 
! this fucntions will try to calculate mass of the disctribution function. This is a debug subroutine
! This one also calculates momentum, temperature and the directional temperature as well as additional moments
! and will update the moments history
!!!!!!!!!!!!!!!!!
subroutine MassCheckRecPlusPlus (f,n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w,mom3_u,mom3_v,mom3_w, & 
                   mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)

use DGV_commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,gasR

real (DP), dimension (:), intent (in) :: f !the vector of nodal values of the distribution function
real (DP), intent (out) :: n,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w,mom3_u,mom3_v,mom3_w,&
                   mom4_u,mom4_v, mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w ! number density, av_v
!!!!!!!!!!!!!!!!!!

integer (I4B) :: i ! scrap index

!!!!!!!!!
! check mass
n=sum(f*nodes_gwts)
!!!!!!!!!
! check momentum 
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
! check directional temperatures 
!!!!!!!!!
temp_u = sum(f*nodes_gwts*(nodes_u-ubar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tx
temp_v = sum(f*nodes_gwts*(nodes_v-vbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Ty
temp_w = sum(f*nodes_gwts*(nodes_w-wbar)**2)/n/3.0_DP*2.0_DP ! dimensionless directional temperature Tz
!!!!!!!!!
! Moments for Q = u^3
mom3_u = sum(f*nodes_gwts*(nodes_u-ubar)**3)/n
mom3_v = sum(f*nodes_gwts*(nodes_v-vbar)**3)/n
mom3_w = sum(f*nodes_gwts*(nodes_w-wbar)**3)/n
!!!!!!!!!
! Moments for Q = u^4
mom4_u = sum(f*nodes_gwts*(nodes_u-ubar)**4)/n
mom4_v = sum(f*nodes_gwts*(nodes_v-vbar)**4)/n
mom4_w = sum(f*nodes_gwts*(nodes_w-wbar)**4)/n
!!!!!!!!!
! Moments for Q = u^5
mom5_u = sum(f*nodes_gwts*(nodes_u-ubar)**5)/n
mom5_v = sum(f*nodes_gwts*(nodes_v-vbar)**5)/n
mom5_w = sum(f*nodes_gwts*(nodes_w-wbar)**5)/n
!!!!!!!!!
! Moments for Q = u^6
mom6_u = sum(f*nodes_gwts*(nodes_u-ubar)**6)/n
mom6_v = sum(f*nodes_gwts*(nodes_v-vbar)**6)/n
mom6_w = sum(f*nodes_gwts*(nodes_w-wbar)**6)/n
!!!!!!!!!

end subroutine MassCheckRecPlusPlus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Next, we compute the Tensor and also the Determinant and inverse of the tensor for use later
! in the ES-BGK distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MakeTensor (f,Tensor)
! Compute the tensor T = (1-alpha)*TRI + 1/n*< c (x) c >
!
!          u   v   w
!     u | T11 T12 T13 |
! T = v | T21 T22 T23 |
!     w | T31 T32 T33 |
!
! This subroutine produces the tensor that we see above associated with the ES-BGK distribution

use DGV_commvar, only: nodes_gwts,nodes_u,nodes_v,nodes_w,alpha,gasR !need to introduce alpha in the read/write parameters.dat file

!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (3,3), intent (out) :: Tensor ! both the tensor and the inverse tensor for the ES-BGK
!!!!!!!!!!!

real (DP) :: ubar, vbar, wbar ! the bulk velocities
real (DP) :: n, Temp ! the density and temperature

!!!!!!!!
! density
n=sum(f*nodes_gwts)
!!!!!!!!
! momentum
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
!!!!!!!!
! temperature
Temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!
Tensor(1,1) = sum(f*nodes_gwts*(nodes_u-ubar)**2)/n*alpha*2.0_DP
Tensor(2,2) = sum(f*nodes_gwts*(nodes_v-vbar)**2)/n*alpha*2.0_DP
Tensor(3,3) = sum(f*nodes_gwts*(nodes_w-wbar)**2)/n*alpha*2.0_DP
Tensor(1,2) = sum(f*nodes_gwts*(nodes_u-ubar)*(nodes_v-vbar))/n*alpha*2.0_DP
Tensor(1,3) = sum(f*nodes_gwts*(nodes_u-ubar)*(nodes_w-wbar))/n*alpha*2.0_DP
Tensor(2,3) = sum(f*nodes_gwts*(nodes_v-vbar)*(nodes_w-wbar))/n*alpha*2.0_DP
Tensor(2,1) = Tensor(1,2)
Tensor(3,1) = Tensor(1,3)
Tensor(3,2) = Tensor(2,3)

Tensor(1,1) = Tensor(1,1) + (1-alpha)*Temp
Tensor(2,2) = Tensor(2,2) + (1-alpha)*Temp
Tensor(3,3) = Tensor(3,3) + (1-alpha)*Temp

end subroutine MakeTensor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Here, the Determinant of the tensor is computed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DetTensor(Tensor) result (Det)
! Evaluate the determinant of the tensor
real (DP), dimension(3,3), intent (in) :: Tensor
real (DP) :: Det
real (DP) :: a, b, c, d, e, f

a = Tensor(1,1)*Tensor(2,2)*Tensor(3,3)
b = Tensor(1,1)*Tensor(2,3)*Tensor(3,2)
c = Tensor(1,2)*Tensor(2,1)*Tensor(3,3)
d = Tensor(1,2)*Tensor(2,3)*Tensor(3,1)
e = Tensor(1,3)*Tensor(2,1)*Tensor(3,2)
f = Tensor(1,3)*Tensor(2,2)*Tensor(3,1)

Det = a-b-c+d+e-f
end function DetTensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
!
!!!!!!!!!!!!!!!
function inv(A) result(Ainv)
  real(dp), dimension(:,:), intent(in) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info=0

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)

  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
end function inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! linleastsqrs(A,b) 
! 
! This subbroutine solves linear least squares problem
! 
! \| Ac-b \|^2 -> min 
! 
! For a rectangular matrix A # rows > # columes 
! 
! Ths solution consists of perfroming SVD decomposition, regularization by eliminating small singular values and 
! solving the least squares problem for the regularized system 
! 
! Returns the solution to the linear least squares problem. Depends on LAPACK and BLAS.
!
!!!!!!!!!!!!!!!
function linleastsqrs(A,b) result (c)

intrinsic INT,Real,LOG,MIN
  
!
  real (DP), dimension(:,:), intent(in) :: A
  real (DP), dimension(:), intent(in) :: b
  real (DP), dimension(size(A,2)) :: c

  real (DP), dimension(1:size(A,1),1:size(A,2)) :: A_temp ! copy of A to pass in LAPACK driver will be overwritten
  real (DP), dimension(1:size(b,1)) :: b_temp ! copy of b to pass to LAPACK driver -- will be owerriden by the solution
  real (DP), dimension(1:size(A,1)) :: S ! this matrix will store singular values of A
  integer :: rank  ! this will store the effective rank of A
    
  real (DP), dimension(:), allocatable :: work  ! work array for LAPACK
  integer, dimension(:), allocatable :: iwork ! work array for LAPACK
   
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: m,n,loc_alloc_stat, info=0
  !!!!!!!!!!!!!!!!!!!!! PARAMETERS OF THE LAPACK DRIVER
  real (DP) :: rcond = 1.0D-5 ! this is the treshlold value for SVD regularization. All singular values 
                  ! s(i)<rcond*s(1) will be treated as 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters that are used by LAPACK. Work array sizes and related 
  integer :: SVDblksz,SVDnlvl,SVDlwork,SVDliwork
  ! End of parametersused by LAPACK
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
! External procedures defined in LAPACK
external DGELSD
  
! Interface for the external subrouinte defined in LAPACK
interface 
 function ILAENV (ispec, name, opts, n1, n2, n3, n4 ) result (y)
  integer :: ispec
  character (len=*) :: name
  character (len=*) :: opts
  integer :: n1, n2, n3, n4
  integer :: y
 end function ILAENV
end interface  
! 
  
m=size(A,1); n=size(A,2) ! record the number of raws and columns
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Before a LAPACK subroutine can be called, a few parameters need 
! to be determined. These parameters determine the size of the 
! workspaces used by the LAPACK driver DGELSD. Subrouine ILAENV is used to determin these 
! parametes.
!!!!!!!!!!!!!
 SVDblksz = ILAENV(9,'dgelsd','',m,n,1,-1) ! see LAPACK manual for details on ILAENV
 SVDnlvl = 10 !not sure how this works. Line suggested in LAPACK (did not work): !INT( LOG(Real(MIN(m,n),SP)/Real((SVDblksz + 1),SP))/LOG(2.0_SP) ) + 1
 SVDlwork = 12*n+2*n*SVDblksz + 8*n*SVDnlvl + n*1 + (SVDblksz+1)**2
 SVDliwork = 3*MIN(m,n)*SVDnlvl + 11*MIN(m,n)
!!!!!!!!!!!!
! in general, after the first call, the values will be updated based on the output of the work arrays
! however, because of the potential conflicts when subroutine is called from an OpenMP parallel region, 
! we choose not to keep the updated values. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate the work arrays: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (work(1:SVDlwork), iwork(1:SVDliwork), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "l2min: Allocation error for variables (work,iwork)"
 stop
end if
!!!!!!!!!!!!!!!!!!!!
! Store A in A_temp and b in b_temp to prevent it from being overwritten by LAPACK
A_temp = A; b_temp=b
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DGELSD is the driver for linear least squares problem with svd decompolision and regularization
! the SVD is computed using the divide and conqer method. See LAPACK manual for details on the subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call DGELSD(m,n,1,A_temp,m,b_temp,m,S,rcond,rank,work,SVDlwork,iwork,info)
if (info /= 0) then
 print *, "l2min: Least squares problem fail. info=", info
 stop 
end if
! if info=0 then recommended values of SVDlwork and SVDliwork should be located in work(1) and iwork(1) respectively
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! UNCOMMEND AND MOVE variables TO COMMON BLOCK/SAVE BLOCK TO RE-USE in consequent calls if not calling subroutine from a parallel region 
!!!! SVDlwork = INT(work(1))  ! these are values recommended by LAPACK for future use
!!!! SVDliwork = iwork(1)  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the solution to the least squares problem is in b_temp(1:n)
c = b_temp(1:n)
! all done
deallocate (work, iwork)
end function linleastsqrs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! getnu(c_0,Cco,psi_u,psi_v,psi_w) 
!
! this subroutine assembles/reconstructs/evaluates the velocity dependent collision frequency
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getnu(c_0,Cco,nodes_u,nodes_v,nodes_w,nuB)

real (DP), intent(in) :: c_0                ! the nu_BGK is in this one 
real (DP), dimension(:), intent(in) :: Cco  ! the coefficients of the velocity dependent collision frequency (VDCF)
real (DP), dimension(:), intent (in) :: nodes_u, nodes_v, nodes_w ! Velocity modes where VDCF needs to be evaluated 
real (DP), dimension(:), intent (out) :: nuB
!!!
integer (I4B) :: i ! scrap variable

!!!!!
nuB(:) = c_0
do i=1,size(Cco,1)
  nuB(:) = nuB(:) + Cco(i)*psi_basis_funcs(i,nodes_u,nodes_v,nodes_w)
end do 
! all done 
end subroutine getnu



function psi_basis_funcs(n,nodes_u,nodes_v,nodes_w) result (y)

integer (I4B), intent (in) :: n ! the number of the basis function. 
real (DP), dimension(:) :: nodes_u,nodes_v,nodes_w ! the nodes at which the basis function is evaluated. 
real (DP), dimension(1:size(nodes_u,1)) :: y ! the value of the function 
!
!!
select case (n)
 case(1)
  y = 1 + nodes_u*0.0d0
 case(2)
  y = nodes_u
 case(3)
  y = nodes_v
 case(4)
  y = nodes_w
 case(5)  
  y = nodes_u**2+nodes_v**2+nodes_w**2
 case(6)
  y = 0.5d0*(3.0d0*nodes_u**2 - 1.0d0)
 case(7)
  y = 0.5d0*(5.0d0*nodes_u**3 - 3.0d0*nodes_u)
 case(8)
  y = 0.5d0*(5.0d0*nodes_v**3 - 3.0d0*nodes_v)
 case(9)
  y = 0.5d0*(3.0d0*nodes_u**2 - 1.0d0)*nodes_v
 case(10)
  y = 0.5d0*(3.0d0*nodes_v**2 - 1.0d0)*nodes_u
 case(11)
  y = 0.125d0*(35.0d0*nodes_u**4 - 30.0d0*nodes_u**2 + 3.0d0)
 case(12)
  y = 0.125d0*(35.0d0*nodes_v**4 - 30.0d0*nodes_v**2 + 3.0d0)
 case(13)
  y = 0.5d0*(3.0d0*nodes_u**2 - 1.0d0)*0.5d0*(3.0d0*nodes_v**2 - 1.0d0) 
 case(14)
  y = 0.5d0*(5.0d0*nodes_u**3 - 3.0d0*nodes_u)*nodes_v
 case(15)
  y = 0.5d0*(5.0d0*nodes_v**3 - 3.0d0*nodes_v)*nodes_u
 case default
  print *, "Can not process: n=", n
  stop
end select 
end function psi_basis_funcs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kernls_enfrsd_moms(n,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar) result (y)
! 
! This function is used to evaluate kernels of the enforced moments. 
! 
! The subroutine is desiged to provide a flexibility in which moments to enforce.
! Change the kernel functions below to determined the enformced moment. 
! 
! Moments are enforced beginning from 1 and finishing with some K < MaxNumEnforcedMoments
! with relaxation of all moments between 
!
!
! Attention! This version is designed to work with EvalColRelES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function kernls_enfrsd_moms(n,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar) result (y)

integer (I4B), intent (in) :: n ! the number of the basis function. 
real (DP), dimension(:) :: nodes_u,nodes_v,nodes_w ! the nodes at which the basis function is evaluated. 
real (DP), intent (in) :: ubar,vbar,wbar ! values of the bulc velocity to be used in the evaluation of the moment kernels 
real (DP), dimension(1:size(nodes_u,1)) :: y ! the value of the function 
!
!integer (I4B) :: j ! some local counter


!!
select case (n)
 case(1)
  y = (nodes_u - ubar)**2 ! moment (u-\bar{u)) (u-\bar{u))  
 case(2) 
  y = (nodes_v - vbar)**2 ! moment (v-\bar{v)) (v-\bar{v)) 
 case(3)
  y = (nodes_w - wbar)**2 ! moment (w-\bar{w)) (w-\bar{w)) 
 case(4)
  y = (nodes_u - ubar)*(nodes_v-vbar)  ! moment (u-\bar{u)) (v-\bar{v))
 case(5)  
  y = (nodes_u - ubar)*(nodes_w-wbar)  ! moment (u-\bar{u)) (w-\bar{w))
 case(6)
  y = (nodes_v - vbar)*(nodes_w-wbar)  ! moment (v-\bar{u)) (v-\bar{v))
 case default
  print *, "kernls_enfrsd_moms: Can not process: n=", n
  stop
end select 
end function kernls_enfrsd_moms

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ComputeRelaxRatesBCI_DGV
!! 
!!
!! This is a subroutine that computes relaxation rates for a selected group of moments from the provided  
!! values of the Boltzmann collision integral (BCI).
!!
!! Values of the collision integral are expected to be evaluated on the secondary DG mesh. Values of the solution will be 
!! expected on the primary mesh.
!!
!!
!! The derivatives of moments are evaluated for the problem of spatially homogeneous relaxation. 
!! The derivatives are thus computed by taking the moment of the collision integral. Also, derivatives of the 
!! local maxwellian are computed for true problem of spatial relaxatiopn this derivative is zero, however, due to the numerical erros, the 
!! maxwellian will change. These errors will tell us how bad is the violation of the conservation laws 
!! 
!! The relaxation speeds will be computed from the following definition: 
!! \nu_{\phi} = -\partial_{t} \ln | f_{\varphi}(t) - f_{\varphi}^M (t)|
!! and using the following formula:L
!! \nu_{\phi} = -\int(Q(f)\phi dv)/ | f_{\varphi}(t) - f_{\varphi}^M (t)|
!!
!! from commvar: fcolII -- the evaluated collision operator on the secondary velocity mesh 
!! from commvar: f -- the value of the solution on the primary mesh
!!
!! 
!! L1_SENS_TRSHLD = 100.0 ! This parameters determines the level of sensitivety of expression for evaluation of the relaxation rates. 
!!          ! the energy moment of the Botlzmann collision integral has to be zero. If it is not zero, this is only due to truncation errors =e_{2}. Thus this moment 
!!         ! can be used to judge about the truncation moments in moments of the collision operator. (The errors expeced to be larger for higher moments). 
!!          ! when evaluating the relaxation rates, the moments $m=\int_{R^3} Q(f,f)\phi$ evaluated with error $e_{m}$. We will hope that $e_{m}$ is comparable to $e_2$. 
!!          ! We will use this parameter to make sure that 
!!          ! m> L1_SENS_TRSHLD * e_{2}. Similarly we will compare $f_{\phi}-f^{M}_{\phi}$  -- the 
!!          ! difference between the moment and the same moment evaluated on a local Maxwellan -- to $L1_SENS_TRSHLD em_{2}$. If this difference is at least (L1_SENS_TRSHLD \cdot e_{2}) 
!!          ! , that is it is at least L1_SENS_TRSHLD times larger than the expected errors then 
!!          ! rate = (m+e_{m})/(f_{\phi}-f^{M}_{\phi})= true rate + (e_m/e_2)1/L1_SENS_TRSHLD
!!         ! 
!!          ! thus, in a sense, L1_SENS_TRSHLD gives the number of correct digits in the computed rate IF e_m\approx e_2 
!!          ! 
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ComputeRelaxRatesBCI_DGV(fcolII,f,momsrates,L1_SENS_TRSHLD)

use DGV_commvar, only:	nodes_u,nodes_v,nodes_w,nodes_gwts,&
						nodes_uII,nodes_vII,nodes_wII,nodes_gwtsII,&
						MaxNumEnforcedMoments, Order, &
						gasR, mol_diam, kBoltzmann,C_inf,gas_viscosity, & ! Parameters of dimensionless reduction and the gas parameters.
						flag_use_VSS_diam,gas_alpha,mol_mass,gas_T_reference
						
use DGV_distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension(:), intent(in) :: fcolII ! the value of the collision integral  -- just the evaluation of the collision integral -- coefficient of the dimensionless reduction to make it into the derivative of the solution have not been added yet.
                          !  frhsII comes defined on the secondary velocity grid
real (DP), dimension(:), intent(in) :: f ! the solution(velocity distribution) on one spatial cell
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent (in) :: L1_SENS_TRSHLD ! sensitivity threshhold parameters for evaluation of the relaxation times
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: i ! scrap indices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: der_mom_tresh, dif_mom_tresh ! this is a scrap variable to use in evaluation of relaxation times
real (DP) :: n,ubar,vbar,wbar,temp ! scrap variables to store the five mandatory moments
real (DP) :: dif_mom, der_mom ! temporary variables to store values of the moments obtained form the solution, from the local 
				! maxwellian and from the collisio noperator  
real (DP), dimension (:), allocatable :: fm ! this is the storage for the local Maxwellian and the derivative of the local maxwellian 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional variables for computing the reference diameter of the VSS model to reflect gas viscosity law??? 
real (DP) :: dimensional_temp,gas_current_viscosity,ref_mol_diam_VSS
integer :: loc_alloc_stat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we compute the moments, we will need them  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!
! check mass
n = sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! check momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
ubar=sum(f*nodes_gwts*nodes_u)/n
vbar=sum(f*nodes_gwts*nodes_v)/n
wbar=sum(f*nodes_gwts*nodes_w)/n
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
temp = sum(f*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we compute the required moments of the solution. Notive that we actualy never need the moments of the 
! function. We only need the differences: the moment of the function minus the moment of the maxwellian. 
! so we will compute the differences betwwen the moment of the fuction and the moment of the local maxwellian.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will calculate the local maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (fm(1:size(f,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "ComputeRelaxRatesBCI_DGV: Allocation error for variables (fm)"
end if 
fm = maxwelveldist(temp,ubar,vbar,wbar,n,nodes_u,nodes_v,nodes_w) ! now we have the maxwellian with the same macroparamters.
!!!!! 
! Now we compute the moments of (f-f_{M}). We notice that the kernes for the moments are
! programmed in the function: kernls_enfrsd_moms. If one desires to change the moments that are enforced -- the change can be 
! accomplished by changing the kernels. One simply has to substitute the kernels that need to be enforced starting with the first kernel. 
! so, let us evaluate the (f_{\phi}-f^{M}_{\phi})
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check for erroneous values 
if ((Order > MaxNumEnforcedMoments) .or. (Order>size(momsrates,1))) then !size(momsrates) should be equal to Order
 print *, "ComputeRelaxRatesBCI_DGV: Value of parameter Order exceeds the maximum allowed", Order, & 
                        MaxNumEnforcedMoments, size(momsrates,1)
 stop
end if 
! evaluate the relaxation rates for momnents programmed in kernls_enfrsd_moms:
fm = f-fm ! now fm keeps the difference
!! compute error indicator for truncation errors in the evaluation of the collision integral
der_mom_tresh = abs(L1_SENS_TRSHLD*sum(fcolII*nodes_gwtsII*((nodes_uII-ubar)**2+(nodes_vII-vbar)**2+(nodes_wII-wbar)**2))/n) ! this is the energy of the collision operator. This will be our error indicator. Theoretically 
								! this quantity much be zero, but roundoff errors will prevet it. errors in the high order moments are expected to be worse. 
								! but at least we have some little idea on the errors. 
dif_mom_tresh = abs(L1_SENS_TRSHLD*sum(fm*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))/n) ! this has to be zero. is the energy of the (f-f_{M}). This will be our error i
								! ndicator for the error of evaluating moments. Theoretically 
								! this quantity much be zero, but truncation errors will prevent this. Errors in the high order moments are expected to be worse,
								! but at least we have some little idea on the errors. 
!!
do i=1,Order
 dif_mom = sum(fm*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar)) ! that is the difference between the moment of the solution and the moment of the local maxwellian (division by density /n is omitted.
 der_mom = sum(fcolII*nodes_gwtsII*kernls_enfrsd_moms(i,nodes_uII,nodes_vII,nodes_wII,ubar,vbar,wbar)) ! that is the moment of the collision operator. Computed on secondary mesh. division by density n is onitted too 
 if ((abs(der_mom) > der_mom_tresh) .and. (abs(dif_mom)>dif_mom_tresh)) then 
   momsrates(i) = -der_mom/dif_mom 
   !!
   if (flag_use_VSS_diam) then 
    !!!!!!!!!!!!!!!! THIS MOLECULAR DIAMETER SHOULD BE USED WHEN HARD SPHERES Diameter need to be adjusted to exhibit correct viscosity...
    dimensional_temp = temp*(C_inf**2)/2.0d0/gasR
    gas_current_viscosity = gas_viscosity*(dimensional_temp/gas_T_reference)**gas_alpha
    ref_mol_diam_VSS = sqrt( 15.0d0*mol_mass*sqrt(gasR*dimensional_temp/pi25DT) &
             /2.0d0/(5.0d0 - 2.0d0*gas_alpha)/(7.0d0-2.0d0*gas_alpha)/gas_current_viscosity ) ! use VSS alpha = 1 for hard spheres
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! we need to add a coeffieint to the computed rates: First, we add the dimensionless constant from the Boltzamnn collision integral: (d^2/L_inf^2)N_inf  and then we factor 
    ! out the dimensionless coefficien of the BGK/ES-BGK and shakhov models: (kBotlzm N_inf C_inf)/(2 gasR L_inf^2 gas_viscosity). After cancelling,we end up with 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    momsrates(i) = momsrates(i)*(ref_mol_diam_VSS**2/kBoltzmann)*(gasR/C_inf)*gas_viscosity*2.0_DP ! this one is using diameter of hard spheres that was adjusted to reflect the gas viscosity law
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   else 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! we need to add a coeffieint to the computed rates: First, we add the dimensionless constant from the Boltzamnn collision integral: (d^2/L_inf^2)N_inf  and then we factor 
    ! out the dimensionless coefficien of the BGK/ES-BGK and shakhov models: (kBotlzm N_inf C_inf)/(2 gasR L_inf^2 gas_viscosity). After cancelling,we end up with 
    momsrates(i) = momsrates(i)*(mol_diam**2/kBoltzmann)*(gasR/C_inf)*gas_viscosity*2.0_DP ! Either this line on the two lines below should be uncommented! 
    !!!!!!!!!! DEBUG : The three lines below do the same as one line above
    ! momsrates(i) = momsrates(i)*((mol_diam/L_inf)**2)*N_inf
    ! print *, "rates", momsrates
    ! momsrates(i) = momsrates(i)*gas_viscosity*2.0_DP*gasR*(L_inf**2)/N_inf/kBoltzmann/C_inf
    !!!!!! END DEBUG
   end if  
 else 
   momsrates(i) = 0.0_DP  ! default and not reliable rate ! 
 end if
end do 

deallocate (fm)
!!
! now the final result can be multiplied by the same coefficient as the ES-BGK model.
end subroutine ComputeRelaxRatesBCI_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ComputeRelaxRatesBCIa_DGV (a is for accelerated, just a tiny saving, really)
!! 
!! This is a copy of the above subroutine, except some quantities are recycled from their previous evaluation. Original subroutine computes them directly:
!!
!! This is a subroutine that computes relaxation rates for a selected group of moments from the provided  
!! values of the Boltzmann collision integral (BCI).
!!
!! Values of the collision integral are expected to be evaluated on the secondary DG mesh. Values of the solution will be 
!! expected on the primary mesh.
!!
!!
!! The derivatives of moments are evaluated for the problem of spatially homogeneous relaxation. 
!! The derivatives are thus computed by taking the moment of the collision integral. Also, derivatives of the 
!! local maxwellian are computed for true problem of spatial relaxatiopn this derivative is zero, however, due to the numerical erros, the 
!! maxwellian will change. These errors will tell us how bad is the violation of the conservation laws 
!! 
!! The relaxation speeds will be computed from the following definition: 
!! \nu_{\phi} = -\partial_{t} \ln | f_{\varphi}(t) - f_{\varphi}^M (t)|
!! and using the following formula:L
!! \nu_{\phi} = -\int(Q(f)\phi dv)/ | f_{\varphi}(t) - f_{\varphi}^M (t)|
!!
!! from commvar: fcolII -- the evaluated collision operator on the secondary velocity mesh 
!! from commvar: f -- the value of the solution on the primary mesh
!!
!! 
!! L1_SENS_TRSHLD = 100.0 ! This parameters determines the level of sensitivety of expression for evaluation of the relaxation rates. 
!!          ! the energy moment of the Botlzmann collision integral has to be zero. If it is not zero, this is only due to truncation errors =e_{2}. Thus this moment 
!!         ! can be used to judge about the truncation moments in moments of the collision operator. (The errors expeced to be larger for higher moments). 
!!          ! when evaluating the relaxation rates, the moments $m=\int_{R^3} Q(f,f)\phi$ evaluated with error $e_{m}$. We will hope that $e_{m}$ is comparable to $e_2$. 
!!          ! We will use this parameter to make sure that 
!!          ! m> L1_SENS_TRSHLD * e_{2}. Similarly we will compare $f_{\phi}-f^{M}_{\phi}$  -- the 
!!          ! difference between the moment and the same moment evaluated on a local Maxwellan -- to $L1_SENS_TRSHLD em_{2}$. If this difference is at least (L1_SENS_TRSHLD \cdot e_{2}) 
!!          ! , that is it is at least L1_SENS_TRSHLD times larger than the expected errors then 
!!          ! rate = (m+e_{m})/(f_{\phi}-f^{M}_{\phi})= true rate + (e_m/e_2)1/L1_SENS_TRSHLD
!!         ! 
!!          ! thus, in a sense, L1_SENS_TRSHLD gives the number of correct digits in the computed rate IF e_m\approx e_2 
!!          ! 
!!          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ComputeRelaxRatesBCIa_DGV(fcolII,Df,momsrates,L1_SENS_TRSHLD,Mom_Close_Maxwl_TRSHLD,&
                                        n,ubar,vbar,wbar,temp,nu,momsrates_reliab)

use DGV_commvar, only:	nodes_u,nodes_v,nodes_w,nodes_gwts,&
						nodes_uII,nodes_vII,nodes_wII,nodes_gwtsII,&
						MaxNumEnforcedMoments, Order, &
						gasR, mol_diam, kBoltzmann,C_inf,gas_viscosity,& ! Parameters of dimensionless reduction and the gas parameters.
						L_inf,N_inf,flag_use_VSS_diam,gas_alpha,mol_mass,gas_T_reference
						
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension(:), intent(in) :: fcolII ! the value of the collision integral  -- just the evaluation of the collision integral -- coefficient of the dimensionless reduction to make it into the derivative of the solution have not been added yet.
                          !  frhsII comes defined on the secondary velocity grid
real (DP), dimension(:), intent(in) :: Df ! the difference between the solution and the local maxwellian (computed elsewhere) on the primary mesh on one spatial cell
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent (in) :: L1_SENS_TRSHLD ! sensitivity threshhold parameters for evaluation of the relaxation times
real (DP), intent (in) :: n,ubar,vbar,wbar,temp ! p variables to store the five mandatory moments
real (DP), intent (in) :: nu ! value of the default relaxation rate -- will be assigned to moments for which the ralaxation rates can not be 
                              ! calculated form the Boltzmann collision operator
logical, intent (out) :: momsrates_reliab ! is true if at least one rate has calculated from the Boltzmann collision operator.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), intent (in) :: Mom_Close_Maxwl_TRSHLD         ! This is the second parameter that determines the sensitivity treshhold for evaluation 
                                                         ! of the relaxation rates for moments. It is possible that the moment is already close to 
                                                         ! its final state. In this case, assuming that there is more noice in the derivative of the moment 
                                                         ! than in the difference of the moment (especially if we use course mesh for evaluating the moment 
                                                         ! derivative, the best strategy is not to compute the relaxation rate for this moment. 
                                                         ! In particular, this will avoid computing rates for conserved moments, when sufficent resolution is used.  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: i ! scrap indices 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: der_mom_tresh, dif_mom_tresh ! this is a scrap variable to use in evaluation of relaxation times
real (DP) :: dif_mom, der_mom ! temporary variables to store values of the moments obtained form the solution, from the local 
				! maxwellian and from the collisio noperator  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional variables for computing the reference diameter of the VSS model to reflect gas viscosity law??? 
real (DP) :: dimensional_temp,gas_current_viscosity,ref_mol_diam_VSS


! Now we compute the moments of Df=(f-f_{M}). We notice that the kernes for the moments are
! programmed in the function: kernls_enfrsd_moms. If one desires to change the moments that are enforced -- the change can be 
! accomplished by changing the kernels. One simply has to substitute the kernels that need to be enforced starting with the first kernel. 
! so, let us evaluate the (f_{\phi}-f^{M}_{\phi})
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check for erroneous values 
if ((Order > MaxNumEnforcedMoments) .or. (Order>size(momsrates,1))) then !size(momsrates) should be equal to Order
 print *, "ComputeRelaxRatesBCI_DGV: Value of parameter Order exceeds the maximum allowed", Order, & 
                        MaxNumEnforcedMoments, size(momsrates,1)
 stop
end if 
! evaluate the relaxation rates for momnents programmed in kernls_enfrsd_moms:
!! compute error indicator for truncation errors in the evaluation of the collision integral
der_mom_tresh = abs(L1_SENS_TRSHLD*sum(fcolII*nodes_gwtsII*((nodes_uII)**2+(nodes_vII)**2+(nodes_wII)**2))*2.0d0/3.0d0) ! this is the energy of the collision operator. This will be our error indicator. Theoretically
!der_mom_tresh = abs(L1_SENS_TRSHLD*sum(fcolII*nodes_gwtsII*((nodes_uII-ubar)**2+(nodes_vII-vbar)**2+(nodes_wII-wbar)**2))) ! this is the energy of the collision operator. This will be our error indicator. Theoretically 
!der_mom_tresh = abs(L1_SENS_TRSHLD*sum(fcolII*nodes_gwtsII*((nodes_uII)**2+(nodes_vII)**2+(nodes_wII)**2))/n) ! this is the energy of the collision operator. This will be our error indicator. Theoretically 

								! this quantity much be zero, but roundoff errors will prevet it. errors in the high order moments are expected to be worse. 
								! but at least we have some little idea on the errors. 
dif_mom_tresh = abs(L1_SENS_TRSHLD*sum(Df*nodes_gwts*((nodes_u)**2+(nodes_v)**2+(nodes_w)**2))) ! this has to be zero. is the energy of the (f-f_{M}). This will be our error i
!dif_mom_tresh = abs(L1_SENS_TRSHLD*sum(Df*nodes_gwts*((nodes_u-ubar)**2+(nodes_v-vbar)**2+(nodes_w-wbar)**2))) ! this has to be zero. is the energy of the (f-f_{M}). This will be our error i
!dif_mom_tresh = abs(L1_SENS_TRSHLD*sum(Df*nodes_gwts*((nodes_u)**2+(nodes_v)**2+(nodes_w)**2))/n) ! this has to be zero. is the energy of the (f-f_{M}). This will be our error i
								! ndicator for the error of evaluating moments. Theoretically 
								! this quantity much be zero, but truncation errors will prevent this. Errors in the high order moments are expected to be worse,
								! but at least we have some little idea on the errors. 
!!
momsrates_reliab = .false.  ! reset the flag to false -- no reliable rates
do i=1,Order
 dif_mom = sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar)) ! that is the difference between the moment of the solution and the moment of the local maxwellian (division by density /n is omitted.
 der_mom = sum(fcolII*nodes_gwtsII*kernls_enfrsd_moms(i,nodes_uII,nodes_vII,nodes_wII,ubar,vbar,wbar)) ! that is the moment of the collision operator. Computed on secondary mesh. division by density n is onitted too 
 if ((abs(der_mom) > der_mom_tresh) .and. (abs(dif_mom) > max(dif_mom_tresh, Mom_Close_Maxwl_TRSHLD))) then 
   momsrates(i) = -der_mom/dif_mom
   momsrates_reliab = .true.  ! if at least one rate was calculated using the Boltzmann integra, this flag is turned into .true. 
   !!
   if (flag_use_VSS_diam) then 
    !!!!!!!!!!!!!!!! THIS MOLECULAR DIAMETER SHOULD BE USED WHEN HARD SPHERES Diameter need to be adjusted to exhibit correct viscosity...
    dimensional_temp = temp*(C_inf**2)/2.0d0/gasR
    gas_current_viscosity = gas_viscosity*(dimensional_temp/gas_T_reference)**gas_alpha
    ref_mol_diam_VSS = sqrt( 15.0d0*mol_mass*sqrt(gasR*dimensional_temp/pi25DT) &
             /2.0d0/(5.0d0 - 2.0d0*gas_alpha)/(7.0d0-2.0d0*gas_alpha)/gas_current_viscosity ) ! use VSS alpha = 1 for hard spheres
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! we need to add a coeffieint to the computed rates: First, we add the dimensionless constant from the Boltzamnn collision integral: (d^2/L_inf^2)N_inf  and then we factor 
    ! out the dimensionless coefficien of the BGK/ES-BGK and shakhov models: (kBotlzm N_inf C_inf)/(2 gasR L_inf^2 gas_viscosity). After cancelling,we end up with 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    momsrates(i) = momsrates(i)*(ref_mol_diam_VSS**2/kBoltzmann)*(gasR/C_inf)*gas_viscosity*2.0_DP ! this one is using diameter of hard spheres that was adjusted to reflect the gas viscosity law
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   else 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! we need to add a coeffieint to the computed rates: First, we add the dimensionless constant from the Boltzamnn collision integral: (d^2/L_inf^2)N_inf  and then we factor 
    ! out the dimensionless coefficien of the BGK/ES-BGK and shakhov models: (kBotlzm N_inf C_inf)/(2 gasR L_inf^2 gas_viscosity). After cancelling,we end up with 
    momsrates(i) = momsrates(i)*(mol_diam**2/kBoltzmann)*(gasR/C_inf)*gas_viscosity*2.0_DP ! Either this line on the two lines below should be uncommented! 
    !!!!!!!!!! DEBUG : The three lines below do the same as one line above
    ! momsrates(i) = momsrates(i)*((mol_diam/L_inf)**2)*N_inf
    ! print *, "rates", momsrates
    ! momsrates(i) = momsrates(i)*gas_viscosity*2.0_DP*gasR*(L_inf**2)/N_inf/kBoltzmann/C_inf
    !!!!!! END DEBUG
   end if  
   ! now the final result can be multiplied by the same coefficient as the ES-BGK model.
 else 
   momsrates(i) = nu ! the rate can not be compute reliably -- use provided default value.
 end if 
end do 
end subroutine ComputeRelaxRatesBCIa_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ProjectDGVtoSecMesh(f,fII) 
! 
! This subroutine takes a solution on the primary mesh and projects it to the secondary mesh (on the same velocity domain)
! This is done by usind nodal-DG iterpolation that is built in into the DG velocity discretizxation. Essentially, given the 
! values of the solution on the primary mesh, we need to produce the values on the the seocnday mesh.
! 
!
subroutine ProjectDGVtoSecMesh(f,fII)

use DGV_commvar, only: nodes_primcellII,cells_gou,cells_gov,cells_gow,&
                 cells_ru,cells_rv,cells_rw,cells_lu,cells_lv,cells_lw,&
                 nodes_uII,nodes_vII,nodes_wII,&
                 nodes_pcell,nodes_ui,nodes_vi,nodes_wi,&
                 g_nds_all



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! values of the solution on the primary mesh 
real (DP), dimension (:), intent (out) :: fII  !values of the solution on the secondary mesh
!!!

integer (I4B) :: i,j,pcn, gou, gov, gow ! scrap indices
real (DP) :: fII_nval,y ! scrap variable to accumulate the interpolated value
real (DP) :: unor,vnor,wnor ! scrap variables to keep the normalized value of the velocity node... 
integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
real (DP), dimension(:), allocatable :: g_nds_u, g_nds_v, g_nds_w ! scrap arrays to keep the DG-Gauss nodes 

!!!!!!!!!!!!!!!!!
! We will go node by node of the secondary mesh and we will interpolate the values of the solution on the nodes secondary mesh. 
!!!!!!!!!!!!!!!! 
fII=0 ! reset the solution.
do i=1,size(fII,1)
 fII_nval=0 ! reset the accumulated interpolated value
 pcn=nodes_primcellII(i) ! first we pull out the number of the primary cell for the selected velocity node on the secondary cell
 ! next, we will compute the normalized values of the velocity on the cell: 
 unor = ( nodes_uII(i) - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
 vnor = ( nodes_vII(i) - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
 wnor = ( nodes_wII(i) - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
 ! next we will create a temp array that will store the values of the Gauss nodes needed for interpolation on the cell (pcn)
 ! we record the orders of the gauss nodes in each velocity component for the cell #(pcn)
 gou=cells_gou(pcn)
 gov=cells_gov(pcn)
 gow=cells_gow(pcn)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!
 !  then we need to allocate arrays to keep the nodes: 
 allocate (g_nds_u(gou),g_nds_v(gov),g_nds_w(gow), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then
  print *, "prepare_sDGVpMap_DGV: Allocation error for (g_nds_u),(g_nds_v),(g_nds_w)"
 end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 g_nds_u=g_nds_all(:gou,gou)
 g_nds_v=g_nds_all(:gov,gov)
 g_nds_w=g_nds_all(:gow,gow)
 ! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
 ! basis function for that node and add it to the interpolated value
 do j=1,size(f,1)
  if (nodes_pcell(j) == pcn) then 
   ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
   ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
   y=1.0_DP ! reset the value of the basis function
   y=y*lagrbasfun(nodes_ui(j),unor,g_nds_u)
   y=y*lagrbasfun(nodes_vi(j),vnor,g_nds_v)
   y=y*lagrbasfun(nodes_wi(j),wnor,g_nds_w)
   ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
   fII_nval=fII_nval + f(j)*y   
  end if
 end do 
 ! all nodes on the cell (pcn) have been added.  Can move to the next node on the secondary mesh: 
 fII(i)=fII_nval
 deallocate(g_nds_u,g_nds_v,g_nds_w)
end do 
end subroutine ProjectDGVtoSecMesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setMom_DGV
! 
! This subroutine computes the matrix that is used to determine the coefficients of the velocity dependent collision frequecy 
! (VDCF) in the BGK model
! and the vector of differences between moments and their local equilibrium values  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setMom_DGV(Df,Mom,DifMom,ubar,vbar,wbar) 

use DGV_commvar, only: nodes_gwts, nodes_u, nodes_v, nodes_w

real (DP), dimension (:), intent (in) :: Df ! the vector of nodal values of the difference between Maxwellian and f = f_M-f
real (DP), dimension (:,:), intent (out) :: Mom
real (DP), dimension (:), intent (out) :: DifMom
real (DP), intent (in) :: ubar,vbar,wbar ! values of the local macroparameters
!!!
integer (I4B) :: i,j ! Scrap indices
real (DP), dimension (1:size(Df,1)) :: Dfngwts, phi ! a scrap array

! the following elements of the moments matrix are integrals of the psi(u,v,w)*phi(u,v,w)*(fM - f) terms
Dfngwts = Df*nodes_gwts ! this is a supplementay multiplication, 
do i=1,size(Mom,1) ! this index goes over all enforced moments
 phi = Dfngwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar) ! kernels of the moments
 do j=1,size(Mom,2) ! this index flips through the basis functions in the representation of the velocity dependent collision frequency 
  Mom(i,j) = sum(phi*psi_basis_funcs(j,nodes_u,nodes_v,nodes_w))            ! psi are the basis functions of the VDCF 
 end do
 DifMom(i) = sum(phi) ! these are the values of the differences of the moments from their local equilibrium values
end do 

end subroutine setMom_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CalcFftA_HO
! 
! High order implementation of the fourier convolutions.
!
! This subroutine calculates all the different Fourier combinations of A. The matrices vary on the cells. There is
! a matrix for every pairwise combination of nodes within a cell and node that the element of A was calculated with.
! This comes out to (su*sv*sw)^3 6 dimensional fourier transform. 
!
! Note that subroutine allocates a full complex matrix for A. 
!
! This routine relies of the intel fortran mkl FFT implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftA_HO
use MKL_DFTI
use DGV_commvar, only:A,A_xi,A_xi1,A_phi,nodes_u,Mw,sw,Mv,sv,Mu,su,FA=>FA_poly, &
                      nodes_ui,nodes_vi,nodes_wi,nodes_pcell

complex, dimension(:), pointer :: FA_temp
integer(I4B) :: i,j !loop counters
integer(I4B) :: first_b_index,cur_b_index
integer(I4B) :: loc_alloc_stat
integer(I4B) :: p,p1,g,g1,cb !scraps to store cell, gauss indices, and canonical basis respectively
integer(I4B) :: fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

allocate(FA(1:Mu*Mv*Mw,1:Mu*Mv*Mw,1:(su*sv*sw),1:(su*sv*sw),1:(su*sv*sw)),stat=loc_alloc_stat)
if (loc_alloc_stat > 0) then
  print *,"CalcFftA_HO: Allocation error for FA"
  stop
endif

!Grab the first basis that was computed.
!Make this global so we can get the basis element for later?
first_b_index=-1
first_b_index=minval(A_phi)
if (first_b_index <= 0) then
  print *,"CalcFftA_HO: first_b_index not set correctly"
  stop
end if

FA=0
do i=1,size(A,1)
  !grab the cell and the gauss index
  p=nodes_pcell(A_xi(i))
  p1=nodes_pcell(A_xi1(i))
  g=(nodes_ui(A_xi(i))-1)*sv*sw+(nodes_vi(A_xi(i))-1)*sw+nodes_wi(A_xi(i))
  g1=(nodes_ui(A_xi1(i))-1)*sv*sw+(nodes_vi(A_xi1(i))-1)*sw+nodes_wi(A_xi1(i))
  cb=A_phi(i)-first_b_index+1
  FA(p,p1,g,g1,cb)=A(i)
  FA(p1,p,g1,g,cb)=A(i) ! store transpose
end do

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ Mw,Mv,Mu,Mw,Mv,Mu /))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_NUMBER_OF_TRANSFORMS, su*sv*sw*su*sv*sw*su*sv*sw)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_DISTANCE, size(FA,dim=1)*size(FA,dim=2))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA(:,1,1,1,1))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiFreeDescriptor",fft_status
  stop
end if

end subroutine CalcFftA_HO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CalcFftA_HO_CE
! 
! High order implementation of the fourier convolutions.
!
! This subroutine calculates all the different Fourier combinations of A. The matrices vary on the cells. There is
! a matrix for every pairwise combination of nodes within a cell and node that the element of A was calculated with.
! This comes out to (su*sv*sw)^3 6 dimensional fourier transform. 
!
! Note that subroutine allocates a full complex matrix for A. 
!
! This routine relies of the intel fortran mkl FFT implementation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftA_HO_CE
use MKL_DFTI
use DGV_commvar, only:A,A_xi,A_xi1,A_phi,nodes_u,Mw,sw,Mv,sv,Mu,su,FA=>FA_ce, &
                      nodes_ui,nodes_vi,nodes_wi,nodes_pcell

real (DP), dimension(:,:,:,:,:), allocatable :: FA_temp
complex (DP), dimension(:,:,:,:,:), allocatable :: FA2
complex (DP) :: e1,e2
integer(I4B) :: i,j !loop counters
integer(I4B) :: first_b_index,cur_b_index
integer(I4B) :: loc_alloc_stat
integer(I4B) :: p,p1,g,g1,cb !scraps to store cell, gauss indices, and canonical basis respectively
integer(I4B) :: fft_status
integer(I4B) :: i1,i2,i3,i4,i5
character(len=DFTI_MAX_MESSAGE_LENGTH) :: err_msg
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

allocate(FA(1:((Mu/2+1)*Mv*Mw),1:Mu*Mv*Mw,1:(su*sv*sw),1:(su*sv*sw),1:(su*sv*sw)),&
   FA_temp(1:(Mu*Mv*Mw),1:Mu*Mv*Mw,1:(su*sv*sw),1:(su*sv*sw),1:(su*sv*sw)),&
   FA2(1:(Mu*Mv*Mw),1:Mu*Mv*Mw,1:(su*sv*sw),1:(su*sv*sw),1:(su*sv*sw)),&
   stat=loc_alloc_stat)
if (loc_alloc_stat > 0) then
  print *,"CalcFftA_HO: Allocation error for FA"
  stop
endif

!Grab the first basis that was computed.
!Make this global so we can get the basis element for later?
first_b_index=-1
first_b_index=minval(A_phi)
if (first_b_index <= 0) then
  print *,"CalcFftA_HO: first_b_index not set correctly"
  stop
end if

FA=0
FA_temp=0
do i=1,size(A,1)
  !grab the cell and the gauss index
  p=nodes_pcell(A_xi(i))
  p1=nodes_pcell(A_xi1(i))
  g=(nodes_ui(A_xi(i))-1)*sv*sw+(nodes_vi(A_xi(i))-1)*sw+nodes_wi(A_xi(i))
  g1=(nodes_ui(A_xi1(i))-1)*sv*sw+(nodes_vi(A_xi1(i))-1)*sw+nodes_wi(A_xi1(i))
  cb=A_phi(i)-first_b_index+1
  FA_temp(p,p1,g,g1,cb)=A(i)
  FA_temp(p1,p,g1,g,cb)=A(i) ! store transpose
end do

FA2=FA_temp

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_REAL, 6, (/ Mw,Mv,Mu,Mw,Mv,Mu /))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_NUMBER_OF_TRANSFORMS, su*sv*sw*su*sv*sw*su*sv*sw)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_DISTANCE, size(FA_temp,dim=1)*size(FA_temp,dim=2))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_DISTANCE, size(FA,dim=1)*size(FA,dim=2))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCommitDescriptor",fft_status
  err_msg = DftiErrorMessage(fft_status)
  print *, err_msg
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA_temp(:,1,1,1,1), FA(:,1,1,1,1))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiFreeDescriptor",fft_status
  stop
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULATE IN PLACE

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ Mw,Mv,Mu,Mw,Mv,Mu /))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw, Mw*Mv, Mw*Mv*Mu, Mw*Mv*Mu*Mw, Mw*Mv*Mu*Mw*Mv/))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_NUMBER_OF_TRANSFORMS, (su*sv*sw)**3)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_DISTANCE, size(FA2,dim=1)*size(FA2,dim=2))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA2(:,1,1,1,1))
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "CalcFftA_HO:fft_status=DftiFreeDescriptor",fft_status
  stop
end if


!!!!!!!!!!! Verify the results



do i5=1,su*sv*sw
do i4=1,su*sv*sw
do i3=1,su*sv*sw
do i2=1,Mu*Mv*Mw
do i1=1,Mu*Mv*Mw
  e2=FA2(i1,i2,i3,i4,i5)
  if (i1>Mv*Mw+Mu/2+1) then
    !reconstruct the index
    j=(i1-1)/(Mv*Mw)+1
    j=mod(Mu-j+1,Mu)+1
    i=mod(i1-1,Mu*Mw)+1
    i=j*Mw*Mv+i
    e1=conjg(FA(i,i2,i3,i4,i5))
  else
    e1=FA(i1,i2,i3,i4,i5)
  end if
  
  if(abs(e2-e1)>1e-10) then
    print *,"Not the same"
    stop
  end if
end do
end do
end do
end do
end do
stop

end subroutine CalcFftA_HO_CE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftA
!
!This subroutine calculates the fourier transform of A treating u,v,w,u',v',w' as independent dimensions. The function uses
!intel MKL definitions that used at the top of this file. The intel fortran FFT is used. 
!
!This function allocates memory for the global variable of FA. Once the function is completed, FA stores the value of
! FFT on A. I do not believe that the intel FFT supports any sort of sparse format, so A has to be converted to the full 
! matrix form before it is passed into the intel FFT functions. This might be a large memory demand since FA will be
! full as well. The full version of A will be deallocated after the fourier transforms finish. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcFftA
use MKL_DFTI
use DGV_commvar, only:A,A_xi, A_xi1,nodes_u,Mw,sw,Mv,sv,Mu,su,FA

integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

!Make sure we can allocate the variables we need first.
allocate(FA(1:size(nodes_u),1:size(nodes_u)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "fftOnA: Allocation error for A_full, FA"
endif

FA=0 

!Overwrite the zeros in the non-zero positions with the correct value
do i=1,size(A,1)
    FA(A_xi(i),A_xi1(i))=A(i)
    fA(A_xi1(i),A_xi(i))=A(i) !Also store the transpose value
end do


fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ Mw*sw,Mv*sv,Mu*su,Mw*sw,Mv*sv,Mu*su /))
if(fft_status>0) then
print *, "fft_status=DftiCreateDescriptor: ", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv, Mw*sw*Mv*sv*Mu*su, Mw*sw*Mv*sv*Mu*su*Mw*sw, Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv, Mw*sw*Mv*sv*Mu*su, Mw*sw*Mv*sv*Mu*su*Mw*sw, Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiCommitDescriptor",fft_status
stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA(:,1))
if(fft_status>0) then
print *, "fft_status=DftiComputeForward",fft_status
stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiFreeDescriptor",fft_status
stop
end if

end subroutine CalcFftA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CalcFftZeroPadA
!
! Calculates the fourier transform of A with zero-padding. 
! 
! This function calculates and stores the fourier of A in a global variable. It uses
! the intel mkl FFT to perform it. Zero padding is first performed which will allocate up
! to 64 times the original size of A, so make sure there is plenty of memory
! before calling this function. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftZeroPadA(pad)
use MKL_DFTI
use DGV_commvar, only:A,A_xi,A_xi1,nodes_u,Mw,sw,Mv,sv,Mu,su,FA

integer, intent(in) :: pad !The amount to pad in each dimension

integer :: loc_alloc_stat,i,fft_status,sd1,sd2,sd3,c1,c2
integer :: M2
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

M2=Mu+2*pad
!Make sure we can allocate the variables we need first.
! The allocation requires 8 in each of the two dimensions of the matrix. Since each row represents
! 3 dimensions, and we need to double each dimension, we need 8 times the original size.
allocate(FA(1:M2**3,1:M2**3), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "fftOnA: Allocation error for A_full, FA"
endif

FA=0 

!We assume that Mu = Mv = Mw, so we just need one variable


!Overwrite the zeros in the non-zero positions with the correct value
!To pad by zero, we need to move around the tensor product. 
do i=1,size(A,1)
    sd1=mod(A_xi(i)-1,Mw)+pad
    sd2=(A_xi(i)-1)/Mw
    sd2=mod(sd2,Mv)+pad
    sd3=(A_xi(i)-1)/Mv/Mw+pad
    c1=((sd1*M2+sd2)*M2+sd3)+1
    sd1=mod(A_xi1(i)-1,Mw)+pad
    sd2=(A_xi1(i)-1)/Mw
    sd2=mod(sd2,Mv)+pad
    sd3=(A_xi1(i)-1)/Mv/Mw+pad
    c2=((sd1*M2+sd2)*M2+sd3)+1
    FA(c1,c2)=A(i)
    FA(c2,c1)=A(i) !Also store the transpose value
end do


fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ M2,M2,M2,M2,M2,M2 /))
if(fft_status>0) then
print *, "fft_status=DftiCreateDescriptor: ", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,M2,M2**2,M2**3,M2**4,M2**5/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,M2,M2**2,M2**3,M2**4,M2**5/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiCommitDescriptor",fft_status
stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA(:,1))
if(fft_status>0) then
print *, "fft_status=DftiComputeForward",fft_status
stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiFreeDescriptor",fft_status
stop
end if

end subroutine CalcFftZeroPadA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftAII
!
!This subroutine calculates the fourier transform of A treating u,v,w,u',v',w' as independent dimensions. The function uses
!intel MKL definitions that used at the top of this file. The intel fortran FFT is used. 
!
!This function allocates memory for the global variable of FAII. Once the function is completed, FAII stores the value of
! FFT on A. I do not believe that the intel FFT supports any sort of sparse format, so AII has to be converted to the full 
! matrix form before it is passed into the intel FFT functions. This might be a large memory demand since FAII will be
! full as well. The full version of AII will be deallocated after the fourier transforms finish. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcFftAII
use MKL_DFTI
use DGV_commvar, only:A=>AII,A_xi=>A_xiII, A_xi1=>A_xi1II,nodes_u=>nodes_uII,Mw=>MwII, &
                      sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII,FA=>FAII

integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

!Make sure we can allocate the variables we need first.
allocate(FA(1:size(nodes_u),1:size(nodes_u)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "fftOnA: Allocation error for A_full, FA"
endif

FA=0 

!Overwrite the zeros in the non-zero positions with the correct value
do i=1,size(A,1)
    FA(A_xi(i),A_xi1(i))=A(i)
    fA(A_xi1(i),A_xi(i))=A(i) !Also store the transpose value
end do


fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ Mw*sw,Mv*sv,Mu*su,Mw*sw,Mv*sv,Mu*su /))
if(fft_status>0) then
print *, "fft_status=DftiCreateDescriptor: ", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv, Mw*sw*Mv*sv*Mu*su, Mw*sw*Mv*sv*Mu*su*Mw*sw, Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv, Mw*sw*Mv*sv*Mu*su, Mw*sw*Mv*sv*Mu*su*Mw*sw, Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiCommitDescriptor",fft_status
stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA(:,1))
if(fft_status>0) then
print *, "fft_status=DftiComputeForward",fft_status
stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiFreeDescriptor",fft_status
stop
end if

end subroutine CalcFftAII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! zeroPadFftOnAII
!
! Calculates the fourier transform of A with zero-padding. 
! 
! This function calculates and stores the fourier of A in a global variable. It uses
! the intel mkl FFT to perform it. Zero padding is first performed which will allocate up
! to 64 times the original size of A, so make sure there is plenty of memory
! before calling this function. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine zeroPadFftOnAII
use MKL_DFTI
use DGV_commvar, only:A=>AII,A_xi=>A_xiII, A_xi1=>A_xi1II,nodes_u=>nodes_uII,Mw=>MwII, &
                      sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII,FA=>FAII

integer :: loc_alloc_stat,i,fft_status,sd1,sd2,sd3,c1,c2
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

!Make sure we can allocate the variables we need first.
! The allocation requires 8 in each of the two dimensions of the matrix. Since each row represents
! 3 dimensions, and we need to double each dimension, we need 8 times the original size.
allocate(FA(1:(8*size(nodes_u)),1:(8*size(nodes_u))), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "fftOnA: Allocation error for A_full, FA"
endif

FA=0 

!Overwrite the zeros in the non-zero positions with the correct value
!To pad by zero, we need to move around the tensor product. When we do
do i=1,size(A,1)
    sd1=mod(A_xi(i)-1,Mw)
    sd2=(A_xi(i)-1)/Mv
    sd2=mod(sd2,Mv)
    sd3=(A_xi(i)-1)/Mv/Mu
    c1=((sd1*2*Mu+sd2)*2*Mv+sd3)
    sd1=mod(A_xi1(i)-1,Mw)
    sd2=(A_xi1(i)-1)/Mv
    sd2=mod(sd2,Mv)
    sd3=(A_xi1(i)-1)/Mv/Mu
    c2=((sd1*2*Mu+sd2)*2*Mv+sd3)
    FA(c1+1,c2+1)=A(i)
    FA(c2+1,c1+1)=A(i) !Also store the transpose value
end do


fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 6, (/ 2*Mw*sw,2*Mv*sv,2*Mu*su,2*Mw*sw,2*Mv*sv,2*Mu*su /))
if(fft_status>0) then
print *, "fft_status=DftiCreateDescriptor: ", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,2*Mw*sw,4*Mw*sw*Mv*sv,8*Mw*sw*Mv*sv*Mu*su,16*Mw*sw*Mv*sv*Mu*su*Mw*sw,32*Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,2*Mw*sw,4*Mw*sw*Mv*sv,8*Mw*sw*Mv*sv*Mu*su,16*Mw*sw*Mv*sv*Mu*su*Mw*sw,32*Mw*sw*Mv*sv*Mu*su*Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiCommitDescriptor",fft_status
stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FA(:,1))
if(fft_status>0) then
print *, "fft_status=DftiComputeForward",fft_status
stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiFreeDescriptor",fft_status
stop
end if

end subroutine zeroPadFftOnAII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EigenDecompFFTOnAII
!
! This function should only be called once FAII has been allocated and computed. This function computes the Eigenvalue Decomposition
! of FAII. This is solving a symmetric complex (not symmetric real/ hermitian complex) problem, so involves a lot of steps from LAPACK to
! compute. The sequence of calls is zgebal, zgehrd, zurghr, zhseqr, ztrevc, zgebak
! Based on : https://software.intel.com/en-us/node/469094#4460A4B6-C2F2-42A8-9C21-7A4EBAD7BABC
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EigenDecompFFTOnAII

use DGV_commvar, only: Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII,nodes_u=>nodes_uII, &
                       FARV=>FAIIRV, FAE=>FAIIE

complex (DP), dimension(:, :), allocatable :: FA
real (DP), dimension(1:size(nodes_u,1)) :: scale,rconde,rcondv !Stores the scale from zgebal
real (DP), dimension(1:size(nodes_u,1)*64) :: work,rwork !If this is too small, the program will suddenly exit with exit code 0, so be careful
                                                !The work and rwork arrays are probably way larger than they need to be, but are small relative
                                                !to the FA array, so its probably not a huge deal. 

integer :: n !Rank of FA
integer :: info,ilo,ihi !Required for the call to zgebal
integer :: lwork !length of the work array above
integer :: i,loc_alloc_stat !misc integers
real :: norm, abnrm
lwork=size(work,1)
n=size(nodes_u,1)

allocate(FARV(1:n,1:n), FAE(1:n), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "fftOnA: Allocation error for A_full, FA"
 stop
endif
FARV=0
FAE=0
! This function does all the work of calculating the eigenvalues and the eigenvectors. I only calculate and store the
! right eigenvectors. See the link below for more details. 
! zgeev:  https://software.intel.com/en-us/node/469230
! zgeevx: https://software.intel.com/en-us/node/469232
info=0
!call zgeevx('B','V','V','B',n,FA,n,FAE,FALV,n,farv,n,ilo,ihi,scale,abnrm,rconde,rcondv,work,lwork,rwork,n,info)
call zgeev('N','V',n,FA,n,FAE,%val(0),1,FARV,n,work,lwork,rwork,info)
if(info /= 0) then
    print *, "SvdFAII: Call to zgebal returned info: ", info
    stop
end if 
    
deallocate(FA)
end subroutine EigenDecompFFTOnAII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CalcFftF_HO(f)
!
! This function performs a FFT on a function that is stored in th primary mesh. It assumes
! that the function is 2 dimensional. The first dimension is the index of the cell it is in. The 
! second dimension is the index number. The DFT is taken along the cell indices for each node.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftF_HO(f)

use DGV_commvar, only:Mw,sw,Mv,sv,Mu,su
use MKL_DFTI

complex(DP), dimension(:,:), intent (inout) :: f

integer :: fft_status 
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw,Mv,Mu/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw, Mw*Mv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw, Mw*Mv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_NUMBER_OF_TRANSFORMS, (/0,1,Mw, Mw*Mv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_DISTANCE, size(f,dim=1))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, f(:,1))
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if

end subroutine CalcFftF_HO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftInvF_HO
!
! This function performs an inverse FFT on a function that is stored in th primary mesh. It assumes
! that the function is 2 dimensional. The first dimension is the index of the cell it is in. The 
! second dimension is the index number. The DFT is taken along the cell indices for each node.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftInvF_HO(f)

use MKL_DFTI
use DGV_commvar, only: Mw,sw,Mv,sv,Mu,su

complex(DP), dimension(:,:), intent (inout) :: f
integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw,Mv,Mu/))
if(fft_status>0) then
 print *, "fft_status=DftiCreateDescriptor: ", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw,Mw*Mv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw,Mw*Mv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_BACKWARD_SCALE, 1.0_DP/(Mw*Mv*Mu))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_NUMBER_OF_TRANSFORMS, su*sv*sw)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_DISTANCE, Mu*Mv*Mw)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
 print *, "fft_status=DftiCommitDescriptor",fft_status
 stop
end if
fft_status=DftiComputeBackward(fft_desc_handler, f(:,1))
if(fft_status>0) then
 print *, "fft_status=DftiComputeBackward",fft_status
 stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if

end subroutine CalcFftInvF_HO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftF
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array.
!
! Note that this function relies off of the intel MKL implementation of the FFT!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcFftF(f, FF)
use MKL_DFTI
use DGV_commvar, only:nodes_u,Mw,sw,Mv,sv,Mu,su

complex(DP), dimension(:), intent (in) :: f
complex(DP), dimension(:), intent (out) :: FF
integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler


do i=1,size(nodes_u,1)
  FF(i)=f(i)
end do

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw*sw,Mv*sv,Mu*su/))
if(fft_status>0) then
print *, "fft_status=DftiCreateDescriptor: ", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv/))
if(fft_status>0) then
print *, "fft_status=DftiSetValue:", fft_status
stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiCommitDescriptor",fft_status
stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FF)
if(fft_status>0) then
print *, "fft_status=DftiComputeForward",fft_status
stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
print *, "fft_status=DftiFreeDescriptor",fft_status
stop
end if
end subroutine CalcFftF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftFII
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array. 
!
! This function works on solutions on the secondary mesh.
!
! Note that this function relies off of the intel MKL implementation of the FFT!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcFftFII(f, FF)
use MKL_DFTI
use DGV_commvar, only:nodes_u=>nodes_uII,Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII

complex(DP), dimension(:), intent (in) :: f
complex(DP), dimension(:), intent (out) :: FF
integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler


do i=1,size(nodes_u,1)
  FF(i)=f(i)
end do

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw*sw,Mv*sv,Mu*su/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw, Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FF)
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if
end subroutine CalcFftFII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftInvF
!
! Performs an inverse FFT on any sequence that is defined on the secondary mesh. The results are stored in f
! which is a double precision complex array.
!
! This function works on solutions on the primary mesh.
!
! Note that this function relies off of the intel MKL implementation of the FFT!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftInvF(FF, f)

use MKL_DFTI
use DGV_commvar, only:Mw,sw,Mv,sv,Mu,su

complex(DP), dimension(:), intent (inout) :: FF
complex(DP), dimension(:), intent (out) :: f
integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

f=FF
fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw*sw,Mv*sv,Mu*su/))
if(fft_status>0) then
 print *, "fft_status=DftiCreateDescriptor: ", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw,Mw*sw*Mv*sv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw,Mw*sw*Mv*sv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_BACKWARD_SCALE, 1.0_DP/(Mw*sw*Mv*sv*Mu*su))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
 print *, "fft_status=DftiCommitDescriptor",fft_status
 stop
end if
fft_status=DftiComputeBackward(fft_desc_handler, f)
if(fft_status>0) then
 print *, "fft_status=DftiComputeBackward",fft_status
 stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if

end subroutine CalcFftInvF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftInvZeroPadF
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array. This relies off the intel MKL FFT. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CalcFftInvZeroPadF(f, FF, pad)
use MKL_DFTI
use DGV_commvar, only:nodes_u,Mw,sw,Mv,sv,Mu,su

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
complex(DP), dimension(:), intent (in) :: f !must be zero padded already
complex(DP), dimension(:), intent (out) :: FF
integer (I4B), intent (in) :: pad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: loc_alloc_stat,i,fft_status,sd1,sd2,sd3,c
integer :: M2
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

M2=Mu+2*pad

FF=f
fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ M2,M2,M2/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,M2,M2**2/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,M2,M2**2/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_BACKWARD_SCALE, 1.0_DP/size(FF))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeBackward(fft_desc_handler, FF)
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if
end subroutine CalcFftInvZeroPadF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!CalcFftInvFII
!
! Performs an inverse FFT on any sequence that is defined on the secondary mesh. The results are stored in f
! which is a double precision complex array.
!
! This function works on solutions on the secondary mesh.
!
! Note that this function relies off of the intel MKL implementation of the FFT!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcFftInvFII(FF, f)

use MKL_DFTI
use DGV_commvar, only:Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII

complex(DP), dimension(:), intent (inout) :: FF
complex(DP), dimension(:), intent (out) :: f
integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

f=FF
fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ Mw*sw,Mv*sv,Mu*su/))
if(fft_status>0) then
 print *, "fft_status=DftiCreateDescriptor: ", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,Mw*sw,Mw*sw*Mv*sv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,Mw*sw,Mw*sw*Mv*sv/))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_BACKWARD_SCALE, 1.0_DP/(Mw*sw*Mv*sv*Mu*su))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
 print *, "fft_status=DftiCommitDescriptor",fft_status
 stop
end if
fft_status=DftiComputeBackward(fft_desc_handler, f)
if(fft_status>0) then
 print *, "fft_status=DftiComputeBackward",fft_status
 stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if

end subroutine CalcFftInvFII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!zeroPadFftII
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array. This relies off the intel MKL FFT. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zeroPadFftII(f, FF)
use MKL_DFTI
use DGV_commvar, only:nodes_u=>nodes_uII,Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII

real(DP), dimension(:), intent (in) :: f
complex(DP), dimension(:), intent (out) :: FF
integer :: loc_alloc_stat,i,fft_status,sd1,sd2,sd3,c
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler


FF=0
!zero pad the tensor product
do i=1,size(nodes_u,1)
  sd1=mod((i-1),Mw)
  sd2=(i-1)/Mv
  sd2=mod(sd2,Mv)
  sd3=(i-1)/Mv/Mu
  c=(sd3*2*Mu+sd2)*2*Mv+sd1
  FF(c+1)=f(i)
end do

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ 2*Mw*sw,2*Mv*sv,2*Mu*su/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FF)
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
fft_status=DftiFreeDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiFreeDescriptor",fft_status
  stop
end if
end subroutine zeroPadFftII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!zeroPadInvFftII
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array. This relies off the intel MKL FFT. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zeroPadInvFftII(f, FF)
use MKL_DFTI
use DGV_commvar, only:nodes_u=>nodes_uII,Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII

real(DP), dimension(:), intent (in) :: f
complex(DP), dimension(:), intent (out) :: FF
integer :: loc_alloc_stat,i,fft_status,sd1,sd2,sd3,c
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

FF=0
!zero pad the tensor product
do i=1,size(nodes_u,1)
  sd1=mod((i-1),Mw)
  sd2=(i-1)/Mv
  sd2=mod(sd2,Mv)
  sd3=(i-1)/Mv/Mu
  c=(sd3*2*Mu+sd2)*2*Mv+sd1
  FF(c+1)=f(i)
end do

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ 2*Mw*sw,2*Mv*sv,2*Mu*su/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_BACKWARD_SCALE, 1.0_DP/size(FF))
if(fft_status>0) then
 print *, "fft_status=DftiSetValue:", fft_status
 stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeBackward(fft_desc_handler, FF)
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
end subroutine zeroPadInvFftII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! InvFftOnZeroPaddedF
!
! Performs a FFT on any sequence that is defined on the secondary mesh. The results are stored in FF
! which is a double precision complex array. This relies off the intel MKL FFT. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine FftOnZeroPaddedInvF(f, FF)
use MKL_DFTI
use DGV_commvar, only:nodes_u=>nodes_uII,Mw=>MwII,sw=>swII,Mv=>MvII,sv=>svII,Mu=>MuII,su=>suII

real(DP), dimension(:), intent (out) :: f
complex(DP), dimension(:), intent (inout) :: FF

integer :: loc_alloc_stat,i,fft_status
type(DFTI_DESCRIPTOR), pointer :: fft_desc_handler

fft_status=DftiCreateDescriptor(fft_desc_handler, DFTI_DOUBLE, DFTI_COMPLEX, 3, (/ 2*Mw*sw,2*Mv*sv,2*Mu*su/))
if(fft_status>0) then
  print *, "fft_status=DftiCreateDescriptor: ", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_PLACEMENT, DFTI_INPLACE)
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_INPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiSetValue(fft_desc_handler, DFTI_OUTPUT_STRIDES, (/0,1,2*Mw*sw, 4*Mw*sw*Mv*sv/))
if(fft_status>0) then
  print *, "fft_status=DftiSetValue:", fft_status
  stop
end if
fft_status=DftiCommitDescriptor(fft_desc_handler)
if(fft_status>0) then
  print *, "fft_status=DftiCommitDescriptor",fft_status
  stop
end if
fft_status=DftiComputeForward(fft_desc_handler, FF)
if(fft_status>0) then
  print *, "fft_status=DftiComputeForward",fft_status
  stop
end if
FF=abs(FF)
f=abs(FF)

!TODO Put back to normal

end subroutine FftOnZeroPaddedInvF


end module DGV_dgvtools_mod
