!
! collision_mod.f90
!
! In this module, the subroutines involved in the evaluation of the collision integral are listed.
!
!
! Authors of Subroutines:
! 
! Alex Alekseenko 
! Craig Euler
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module DGV_collision_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

contains 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicA_DGV (f,fc)
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! The subroutine expects that operator A is read from the hard drive
! To read operator A subroutine "ReadAarraysDGV" must be called in the driver just once
! There are many other variables that 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
! These variables help to use the symmetries of operator A with respect to shift. 
! To set the variables, the subroutine "SetCellsUGINdsSftArrs" must be called just once 
! in the driver after "ReadAarraysDGV" has been exectuted
!  
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicA_DGV(f,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: k,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
real (DP), dimension(:), allocatable :: Aphi ! scrap array to keep the portion of A
real (DP), dimension (:), allocatable :: fxi,fxi1 !scrap array to keep the f(xi) and f(xi1) 
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   allocate(fxi(1:phicap),fxi1(1:phicap),Aphi(1:phicap),stat=loc_alloc_stat)
    !
   if (loc_alloc_stat >0) then 
   print *, "EvalCollisionPeriodicA_DGV: Allocation error for variables (fxi,fxi1,Aphi)"
   stop
   end if
   !!!!!!!!!!!!!!!!!
   fxi=0;fxi1=0 ! nullify the arrays
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   k=0 ! counter -- will be used to compress the zero records in arrays... 
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             k=k+1
             fxi(k)=f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )
           !  zz=((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
           !   (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1)
             fxi1(k)=f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )
             Aphi(k)=A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do   
   !! the arrays are set. Proceed to evaluate the matrix product   
   fc(iphi)=2.0_DP*sum(Aphi(1:k)*fxi(1:k)*fxi1(1:k))/nodes_gwts(iphi)     
   deallocate (fxi,fxi1,Aphi)   
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicA_DGV (f,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGV(f,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicAPlus_DGV_OMP (f,fc) 
!
! This is a copy of the above subroutine except OMP pragmas are added.
!
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGV_OMP(f,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts,Num_OMP_threads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
! OpenMP
call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(dui,dvi,dwi,ixi,xicell,iuxicell,ivxicell,iwxicell, &
!$OMP    ixi1,xi1cell,iuxi1cell,ivxi1cell,iwxi1cell,Ashift,phicap,j) &
!$OMP    NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 64)  
!!!!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2.0_DP*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc + &
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc + &
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGV_OMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicAPlus_DGVII_OMP(f,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This is a copy of the above subroutine, except that it works with secondary meshes and uses OpenMP.
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
! 
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGVII_OMP(f,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII,Num_OMP_threads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
! OpenMP
call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(dui,dvi,dwi,ixi,xicell,iuxicell,ivxicell,iwxicell, &
!$OMP    ixi1,xi1cell,iuxi1cell,ivxi1cell,iwxi1cell,Ashift,phicap,j) &
!$OMP    NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 64)  
!!!!
do iphi=1,size(nodes_AshiftII,1)  ! Loop in the nodal points
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2.0_DP*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGVII_OMP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicAPlus_DGVII(f,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This is a copy of the above subroutine, but without OpenMP.
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGVII(f,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial Allocation of the A-arrays 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 

do iphi=1,size(nodes_AshiftII,1)  ! Loop in the nodal points
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2.0_DP*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGV (f,fm,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGV(f,fm,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0_DP*fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGV_OMP (f,fm,fc) 
!
! A copy of the above subroutine with OpenMP pragmas added
!
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGV_OMP(f,fm,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts,Num_OMP_threads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
! OpenMP
call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(iphi,dui,dvi,dwi,ixi,xicell,iuxicell,ivxicell,iwxicell, &
!$OMP    ixi1,xi1cell,iuxi1cell,ivxi1cell,iwxi1cell,Ashift,phicap,j,xi_j,xi1_j) &
!$OMP    NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 64)  
!!!!
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0_DP*fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGV_OMP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGVII (f,fm,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGVII(f,fm,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do iphi=1,size(nodes_AshiftII,1)  ! Loop in the nodal points
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0*fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGVII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGVII_OMP (f,fm,fc) 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGVII_OMP(f,fm,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII,Num_OMP_threads

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
! OpenMP
call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(iphi,dui,dvi,dwi,ixi,xicell,iuxicell,ivxicell,iwxicell, &
!$OMP    ixi1,xi1cell,iuxi1cell,ivxi1cell,iwxi1cell,Ashift,phicap,j,xi_j,xi1_j) &
!$OMP    NUM_THREADS(Num_OMP_threads) &
!$OMP    SCHEDULE(DYNAMIC, 64)  
!!!!
do iphi=1,size(nodes_AshiftII,1)  ! Loop in the nodal points
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0_DP*fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGVII_OMP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicAPlus_DGVII_MPI(f,fc) 
! 
! This is a modification of the subroutine to run in parallel MPI implementation
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This is a copy of the above subroutine, but without OpenMP.
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_DGVII_MPI(f,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII,procs_nodes_wldII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j! scrap counters 
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initial Allocation of the A-arrays 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
do i=1,size(procs_nodes_wldII,1)  ! Loop in the nodal points assigned to this processor
   iphi=procs_nodes_wldII(i)
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2.0_DP*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_DGVII_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicAPlus_MPI_DGV (f,fc) 
!
! This is an MPI analog of the subroutine. On the master processor, only the collective operation is called. On the slave processor, 
! the evaluation of the portion of the collision operator is performed. 
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated on one cell. The values on other cells are
! obtained by a transposition from operator A on the canonical cell --- which should be selected at the 
! center fot he mesh (or close to the center)
! 
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicAPlus_MPI_DGV(f,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts,procs_nodes_wld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
integer :: irank ! the number of the process on which the software is running
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP), dimension (:), allocatable :: fc_buff ! 
integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status
! 
integer :: ierr ! variables for MPI Calls

!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do i=1,size(procs_nodes_wld,1)  ! Loop in the nodal points
   iphi=procs_nodes_wld(i)
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             fc(iphi)=fc(iphi)+2*f( ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) )*&
                               f( ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+&
              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicAPlus_MPI_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGVII_MPI (f,fm,fc) 
!
! This is a modification of the subroutine EvalCollisionPeriodicMixedTermsA_DGVII 
! to run in MPI parallel code
!
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phicanII
!  nodes_duiII,nodes_dviII,nodes_dwiII,
!  cells_ugiII,cells_vgiII,cells_wgiII
!
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGVII_MPI(f,fm,fc)

use DGV_commvar, only: AII,A_capphiII,A_xiII,A_xi1II,cells_gouII,cells_govII,cells_gowII,&
                   grids_cap_uII,grids_cap_vII,grids_cap_wII,nodes_AshiftII,nodes_phicanII,&
                   nodes_duiII,nodes_dviII,nodes_dwiII,nodes_pcellII,cells_ugiII,cells_vgiII,&
                   cells_wgiII,nodes_uiII,nodes_viII,nodes_wiII,nodes_gwtsII,procs_nodes_wldII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j,i! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gouII(1)
gov=cells_govII(1)
gow=cells_gowII(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_uII(1)-1
pgcv=grids_cap_vII(1)-1
pgcw=grids_cap_wII(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
do i=1,size(procs_nodes_wldII,1)  ! Loop in the nodal points assigned to this processor
   iphi=procs_nodes_wldII(i)
   Ashift=nodes_AshiftII(iphi)
   phicap=A_capphiII(nodes_phicanII(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_duiII(iphi)
   dvi=nodes_dviII(iphi)
   dwi=nodes_dwiII(iphi)
   !!
   do j=1,phicap
      ixi = A_xiII(Ashift+j)
      xicell = nodes_pcellII(ixi) 
      iuxicell = cells_ugiII(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgiII(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgiII(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1II(Ashift + j)
         xi1cell = nodes_pcellII(ixi1)       
         iuxi1cell = cells_ugiII(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgiII(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgiII(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_uiII(ixi)-1)*gow*gov + (nodes_viII(ixi)-1)*gow + nodes_wiII(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_uiII(ixi1)-1)*gow*gov + (nodes_viII(ixi1)-1)*gow + nodes_wiII(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*AII(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0_DP*fc(iphi)/nodes_gwtsII(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGVII_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionPeriodicMixedTermsA_DGV_MPI (f,fm,fc) 
!
! this is a modification of EvalCollisionPeriodicMixedTermsA_DGV to enable MPI parallel algorithm
!
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! This SUBROUTINE IS TO BE USED IN THE DECOMPOSITION MODE when storing derivative is not 
! feasible. It evaluates the cross term \int\int f_{i} fm_{j} A^{ij}_{k} 
!
! In addition
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!  nodes_dui,nodes_dvi,nodes_dwi,
!  cells_ugi,cells_vgi,cells_wgi
!
!!!!!!!!!!!!!

subroutine EvalCollisionPeriodicMixedTermsA_DGV_MPI(f,fm,fc)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,nodes_gwts,procs_nodes_wld

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution-maxwellian at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the components of the maxwellian at the current time step.
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j,i! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the collision integral by contracting with a portion of A.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
do i=1,size(procs_nodes_wld,1)  ! Loop in the nodal points
   iphi=procs_nodes_wld(i)
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fc(iphi)=fc(iphi)+( fm(xi_j)*f(xi1_j) + f(xi_j)*fm(xi1_j) )*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do 
   fc(iphi)=2.0_DP*fc(iphi)/nodes_gwts(iphi)
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!
end subroutine EvalCollisionPeriodicMixedTermsA_DGV_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollisionLinear(f,fc) 
!
! This subroutine evaluates the linearized collision kernel. To call this subroutine, 
! a vector of the linearized kernel must be prepared. the vector is 2\int f_{m}(v)A(v,v_1,phi) dv
! it is calcuated in the subroutine PrepareFMA_DGV
! 
! 
! Here the arrays fmA has the linearizised operator in it
! accessed directly from the commvar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine EvalCollisionLinear(f,fc)        ! This evaluates the linear part

use DGV_commvar, only: fmA

intrinsic MATMUL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fc = 2*MATMUL(f,fmA) ! this is a contraction with the gradient , i.e., the derivative of the collision operator at some maxwellian

end subroutine EvalCollisionLinear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine PrepareFMA_DGV (fm)
!
!
! This subroutine prepares the derivative of the collision integral, 2\int f_{M}(v) A(v,v_1) dv
! 
! fm is the current local maxwellian
! 
! This subroutine assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!
! Before calling this function, make sure that  Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! the evaluation of the linearized collision operator.

subroutine PrepareFMA_DGV
use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,cells_gou,cells_gov,cells_gow,&
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_ui,nodes_vi,nodes_wi,fm,fmA,nodes_gwts

!!!!!!!!!!! ALL IMPORTANT VARIABLES ARE ACCESSED DIRECTLY !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: iphi,ixi,ixi1 ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: xi1_j,xi_j,j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: Ashift,phicap ! scrap variables. Ashift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer :: loc_alloc_stat ! variable to keep the allocation status

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! main loop: prepare the dummy arrays, then evaluate the linearized collision operator: fmA=\int_{R^3} f_{M}(v)A(v,v_1,\phi)dv.
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fmA=0 ! nullify the result before computing anything... 
!
do iphi=1,size(nodes_Ashift,1)  ! Loop in the nodal points
   Ashift=nodes_Ashift(iphi)
   phicap=A_capphi(nodes_phican(iphi))
   ! Need to allocate some scrap arrays:
   !!!!!!!!!!!!!!!!!
   dui=nodes_dui(iphi)
   dvi=nodes_dvi(iphi)
   dwi=nodes_dwi(iphi)
   !!
   do j=1,phicap
      ixi = A_xi(Ashift+j)
      xicell = nodes_pcell(ixi) 
      iuxicell = cells_ugi(xicell) + dui
      if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
       ivxicell = cells_vgi(xicell) + dvi
       if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
        iwxicell = cells_wgi(xicell) + dwi
        if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
         ixi1 = A_xi1(Ashift + j)
         xi1cell = nodes_pcell(ixi1)       
         iuxi1cell = cells_ugi(xi1cell) + dui
         if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
          ivxi1cell = cells_vgi(xi1cell) + dvi
          if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
           iwxi1cell = cells_wgi(xi1cell) + dwi
           if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
             ! Now that we know that both shifts were successfull, we set prepare both fxi and fxi1 and aphi records
             xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc+&
                              (nodes_ui(ixi)-1)*gow*gov + (nodes_vi(ixi)-1)*gow + nodes_wi(ixi) 
             xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc+ &
                              (nodes_ui(ixi1)-1)*gow*gov + (nodes_vi(ixi1)-1)*gow + nodes_wi(ixi1) 
             fmA(xi_j,iphi)  = fmA(xi_j,iphi)  + fm(xi1_j)*A(Ashift+j)
             fmA(xi1_j,iphi) = fmA(xi1_j,iphi) + fm(xi_j)*A(Ashift+j)   
             !
           end if
          end if  
         end if 
        end if 
       end if 
      end if 
   end do   
   !! Need to divide by the volume element for the node with the number iphi 
   !! this division appears in the formula for the collision operator
   fmA(:,iphi)=fmA(:,iphi)/nodes_gwts(iphi)
   !! the arrays are set. Proceed to evaluate the linearized solution  
end do ! End of the main loop in nodal points
!
end subroutine PrepareFMA_DGV 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EvalCollisionPeriodicAKor_DGV (f,fc) 
!!
!!! This subroutine is to evaluate the bilinear collision operator using Korobov quadratures.
!!
!! This subroutine assumes a uniform velocity grid and a periodic structure of the 
!! basis functions. Operator Akor is only evaluated for one basis function. The values 
!! of A for other basis functions are obtained by a transposition from operator A on 
!! the canonical cell. (The canonical cell should be selected at the center of the mesh 
!! or close to the center. The value of Akor_phi can be used to determine which basis function 
!! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!! 
!! This subroutine depends on a number of arrays that need to be prepared before the first time step. 
!! Specifically nodes_phicanII
!!  nodes_duiII,nodes_dviII,nodes_dwiII,nodes_phican
!!  cells_ugiII,cells_vgiII,cells_wgiII
!! and  
!!  AkorAllNets_shift 
!! need to be prepared. 
!! To set up the first bunch, call SetCellsUGI_DGV
!! To set up the second bunch, call ...
!! 
!!
!!!!!!!!!!!!!!!
!
subroutine EvalCollisionPeriodicAKor_DGV(f,fc) 
!
use DGV_commvar, only: nodes_phican,nodes_dui,nodes_dvi,nodes_dwi,nodes_gwts, &
                       AkorAllNets, AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                       korob_net_paramAllNets, nodes_AkorshiftAllNets,numKornets, &
                       u_L,u_R,v_L,v_R,w_L,w_R,&
                       cells_ru,cells_lu,cells_rv,cells_lv,cells_rw,cells_lw
                       
use DGV_dgvtools_mod                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: iphi ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: AKorshift,phicap ! scrap variables. AKorshift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer (I4B) :: k ! index of the korobov node
integer :: loc_alloc_stat ! variable to keep the allocation status
integer :: num_nets, curr_net ! a scrap variable to keep the number of all Korobov nets and the number fo the net that is used
real :: harvest ! a scrap variable
real (DP) :: frac_part,xiu,xiv,xiw,xi1u,xi1v,xi1w,du,dv,dw,fval,fval1 ! scrap variables 
integer (I4B), dimension (7) :: korob_net_param ! scrap variable to keep the korobov parameters.
logical :: outsidedomain, outsidedomain1
integer (I4B) :: pcn, pcn1 ! scrap indices to keep numbers of cells for korobov points 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
num_nets = numKornets ! detemine the number of nets
!!!!!!!!!!!!!!!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!
do iphi=1,size(nodes_gwts,1)  ! Loop in the nodal points. Collision operator needs to be evaluated at each node
  !!!! now we randomly pick the number of the net that will be used for integration  
  call RANDOM_NUMBER(harvest)
  if (harvest >= 1.0D0) then !!! to avoid the result harvest  = 1 for sure
   harvest = .999 
  end if
  curr_net = int(harvest*Real(num_nets))+1 ! randomly select the net from 1...num_nets available nets (all supplied nets must be about the same accuracy...)
  !!!!!!! 
  Akorshift = nodes_AkorshiftAllNets(curr_net)%p(nodes_phican(iphi)) ! shift inside the Akor Array to get to correct basis function
  phicap = AkorAllNets_capphi(curr_net)%p(nodes_phican(iphi))               ! number of nonzero entires for this basis function 
  korob_net_param = korob_net_paramAllNets(curr_net)%p          ! save the parameters of the Korobov net into a temp array
  !!!!!!!!!!!!!!!!!
  dui = nodes_dui(iphi) ! We assume that all Korobov nets are obtained for the same velocity discretization. 
  dvi = nodes_dvi(iphi) ! and use the same canonical cell. Otherwise, different nets will require different arrays    
  dwi = nodes_dwi(iphi) ! nodes_dwi
  !!
  do j = 1,phicap  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  Determine the two velocity vectors in the Korobov Quadrature
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   k = AkorAllNets_k(curr_net)%p(j) ! this is the number of the Korobov node for 
   frac_part = Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP) 
   xiu = frac_part*(u_R-u_L) + u_L 
   frac_part = Real( (korob_net_param(3)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(3)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xiv = frac_part*(v_R-v_L) + v_L
   frac_part = Real( (korob_net_param(4)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(4)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xiw = frac_part*(w_R-w_L) + w_L 
   frac_part = Real( (korob_net_param(5)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(5)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1u = frac_part*(u_R-u_L) + u_L
   frac_part = Real( (korob_net_param(6)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(6)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1v = frac_part*(v_R-v_L) + v_L
   frac_part = Real( (korob_net_param(7)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(7)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1w = frac_part*(w_R-w_L) + w_L 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! now we need to shift velocity point to take into account the use of shift in the canonical function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   xiu = xiu + dui*du;   xiv = xiv + dvi*dv;   xiw = xiw + dwi*dw; 
   xi1u = xi1u + dui*du; xi1v = xi1v + dvi*dv; xi1w = xi1w + dwi*dw; 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! now we interpolate the solution at velocity points (xiu,xiv,xiw) and (xi1u,xi1v,xi1w)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   outsidedomain = .false. 
   pcn  =  QuickCellFindUniformGrid_DGV(xiu,xiv,xiw,outsidedomain) ! we locate the cell that contains (xiu,xiv,xiw)
   outsidedomain1 = .false.                                ! if the velocity value outside domain, return 
   pcn1 =  QuickCellFindUniformGrid_DGV(xi1u,xi1v,xi1w,outsidedomain1) ! outsidedomain = .true.
   if (outsidedomain .or. outsidedomain1) then 
    cycle  ! points that are outsied of the mesh are ignored in the sum. 
   end if 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! interpolate the values of the function at (xiu,xiv,xiw) and (xiu1,xiv1,xiw1)
   !!! we compute the normalized values of the velocity on the cell: 
   fval  = EvalSolVeloPtCellFast_DGV(f,xiu,xiv,xiw,pcn)
   fval1 = EvalSolVeloPtCellFast_DGV(f,xi1u,xi1v,xi1w,pcn1)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Add the new value to the integral
   fc(iphi) = fc(iphi) + fval*fval1*AkorAllNets(curr_net)%p(Akorshift+j)   
  end do 
  fc(iphi)=fc(iphi)/nodes_gwts(iphi) ! Add 8/Delta v^{j} / \omega_{i}  
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!!
end subroutine EvalCollisionPeriodicAKor_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EvalCollisionPeriodicAKorOpt_DGV (f,fc) 
!!
!! This is a version of the above subroutine with summations rearranged tryting to 
!! achieve optimaltiy 
!!
!!! This subroutine is to evaluate the bilinear collision operator using Korobov quadratures.
!!
!! This subroutine assumes a uniform velocity grid and a periodic structure of the 
!! basis functions. Operator Akor is only evaluated for one basis function. The values 
!! of A for other basis functions are obtained by a transposition from operator A on 
!! the canonical cell. (The canonical cell should be selected at the center of the mesh 
!! or close to the center. The value of Akor_phi can be used to determine which basis function 
!! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!! 
!! This subroutine depends on a number of arrays that need to be prepared before the first time step. 
!! Specifically nodes_phicanII
!!  nodes_duiII,nodes_dviII,nodes_dwiII,nodes_phican
!!  cells_ugiII,cells_vgiII,cells_wgiII
!! and  
!!  AkorAllNets_shift 
!! need to be prepared. 
!! To set up the first bunch, call SetCellsUGI_DGV
!! To set up the second bunch, call ...
!! 
!!
!!!!!!!!!!!!!!!
!
subroutine EvalCollisionPeriodicAKorOpt_DGV(f,fc) 
!
use DGV_commvar, only: nodes_phican,nodes_dui,nodes_dvi,nodes_dwi,nodes_gwts, &
                       AkorAllNets, AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                       korob_net_paramAllNets, nodes_AkorshiftAllNets,numKornets, &
                       u_L,u_R,v_L,v_R,w_L,w_R,g_nds_all, &
                       cells_ru,cells_lu,cells_rv,cells_lv,cells_rw,cells_lw,&
                       grids_cap_u,grids_cap_v,grids_cap_w, &
                       cells_gou,cells_gov,cells_gow,&
                       nodes_pcell
                       
use DGV_dgvtools_mod                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: iphi ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: AkorNrecs,n_ijl,xi_j,xi1_j,ni,nj,nl,n1i,n1j,n1l ! scrap variables to keep addresses of cells where the velocity point falls
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ii,jj,ll,ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: k ! index of the korobov node
integer :: loc_alloc_stat ! variable to keep the allocation status
integer :: num_nets, curr_net ! a scrap variable to keep the number of all Korobov nets and the number fo the net that is used
real :: harvest ! a scrap variable
real (DP) :: frac_part,xiu,xiv,xiw,xi1u,xi1v,xi1w,du,dv,dw,fxi,fxi1,Akor_val ! scrap variables
real (DP), dimension(:), allocatable :: lagrarry, lagrarry1 
integer (I4B), dimension (7) :: korob_net_param ! scrap variable to keep the korobov parameters.
logical :: outsidedomain, outsidedomain1
real (DP) :: piunor, pjvnor, plwnor,unor,wnor,vnor !scrap values  
integer (I4B) :: pcn, phican ! scrap indices to keep numbers of cell for korobov points, and sracp to keep the index of the basis function 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
num_nets = numKornets ! detemine the number of nets
!!!!!!!!!!!!!!!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!!!!!!!!!!!!
!!! call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!!!!!!!!!!!!
!call RANDOM_NUMBER(harvest)
!if (harvest >= 1.0D0) then !!! to avoid the result harvest  = 1 for sure
! harvest = .999 
!end if
!curr_net = int(harvest*Real(num_nets))+1 ! randomly select the net from 1...num_nets available nets (all supplied nets must be about the same accuracy...)
!!!  Override the above ...
curr_net = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
korob_net_param = korob_net_paramAllNets(curr_net)%p          ! save the parameters of the Korobov net into a temp array
!!! determine the total number of records in the entire Akor array:
AkorNrecs = size(AkorAllNets(curr_net)%p,1)
if (AkorNrecs /= sum(AkorAllNets_capphi(curr_net)%p)) then 
 print *,"EvalCollisionPeriodicAKorOpt_DGV: possible data error. Length of Akor not equal sum(Akor_capphi)" 
 stop 
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (lagrarry(1:gou*gov*gow),lagrarry1(1:gou*gov*gow), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalCollisionPeriodicAKorOpt_DGV: Error allocation arrays (largarray,lagrarry1)"
 stop
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=1,AkorNrecs  ! select  a record from the pre-computed Akor array
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !  Determine the two velocity vectors that correspond to this record 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Akor_val = AkorAllNets(curr_net)%p(j)
 k = AkorAllNets_k(curr_net)%p(j) ! this is the number of the Korobov node for 
 frac_part = Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP) 
 xiu = frac_part*(u_R-u_L) + u_L 
 frac_part = Real( (korob_net_param(3)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(3)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xiv = frac_part*(v_R-v_L) + v_L
 frac_part = Real( (korob_net_param(4)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(4)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xiw = frac_part*(w_R-w_L) + w_L 
 frac_part = Real( (korob_net_param(5)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(5)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1u = frac_part*(u_R-u_L) + u_L
 frac_part = Real( (korob_net_param(6)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(6)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1v = frac_part*(v_R-v_L) + v_L
 frac_part = Real( (korob_net_param(7)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(7)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1w = frac_part*(w_R-w_L) + w_L  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! For selected Korbobov node find were velocities fall on the uniform velocity grid  (not going to work on non-uniforom grids)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 outsidedomain = .false. 
 call QuickCellFindUniformGridIJL_DGV(xiu,xiv,xiw,ni,nj,nl,outsidedomain) ! we locate the cell that contains (xiu,xiv,xiw)
 outsidedomain1 = .false.                                ! if the velocity value outside domain, return 
 call QuickCellFindUniformGridIJL_DGV(xi1u,xi1v,xi1w,n1i,n1j,n1l,outsidedomain1) ! outsidedomain = .true.
 if (outsidedomain .or. outsidedomain1) then 
  print *,"EvalCollisionPeriodicAKorOpt_DGV: data error. Akor contains recs for nodes outside the vel. dom."
  stop  ! points that are outside of the mesh -- must be an error in data or a wrong velocity domain. 
 end if 
 !!! now we know local coordinates of cells where falls the selected korbov velocity node.    
 !!! Next find local coordinates of the canonical basis function for this records
 phican = AkorAllNets_phi(curr_net)%p(j) ! number of the basis function for which entry of Akor was computed 
 !!!! Now we prepare for interpolation of the solution on the shifted korobov nodes... 
 !!! first triple of the node: 
 !!! compute the number of the cell where the node belongs
 pcn = (ni-1)*pgcw*pgcv + (nj-1)*pgcw + nl
 !!! Next we evaluate nodal DG basis fucntion on the cell where the velocities fall. We evaluate then at every nodal point in 
 !!! velcity on that cell and save in an array. This array will be used to later to interpolate the value of the function 
 !!! on each cell. The algorithm depends on the fact that the velocity grid is uniform.   
 unor = (xiu - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
 vnor = (xiv - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
 wnor = (xiw - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
 ! basis function for that node and add it to the interpolated value
 lagrarry=0
 n_ijl = 1
 do ii=1,gou
  piunor = lagrbasfun(ii,unor,g_nds_all(:gou,gou))
  do jj=1,gov
   pjvnor = lagrbasfun(jj,vnor,g_nds_all(:gov,gov))
   do ll=1,gow 
    ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
    ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
    plwnor = lagrbasfun(ll,wnor,g_nds_all(:gow,gow))
    ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
    lagrarry(n_ijl) = piunor*pjvnor*plwnor
    n_ijl = n_ijl+1
   enddo
  enddo 
 enddo  
 !!! second triple of the Korobov node: 
 !!! compute the number of the cell where the node belongs
 pcn = (n1i-1)*pgcw*pgcv + (n1j-1)*pgcw + n1l
 !!! Next we evaluate nodal DG basis fucntion on the cell where the velocities fall. We evaluate then at every nodal point in 
 !!! velcity on that cell and save in an array. This array will be used to later to interpolate the value of the function 
 !!! on each cell. The algorithm depends on the fact that the velocity grid is uniform.   
 unor = (xiu - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
 vnor = (xiv - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
 wnor = (xiw - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
 ! basis function for that node and add it to the interpolated value
 lagrarry1=0
 n_ijl = 1
 do ii=1,gou
  piunor = lagrbasfun(ii,unor,g_nds_all(:gou,gou))
  do jj=1,gov
   pjvnor = lagrbasfun(jj,vnor,g_nds_all(:gov,gov))
   do ll=1,gow 
    ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
    ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
    plwnor = lagrbasfun(ll,wnor,g_nds_all(:gow,gow))
    ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
    lagrarry1(n_ijl) = piunor*pjvnor*plwnor
    n_ijl = n_ijl+1
   enddo
  enddo 
 enddo  
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! We are ready to add this integration node to the collision integral 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iphi=1,size(nodes_gwts,1)  ! Loop in the nodal points. Collision operator needs to be evaluated at each node
   if (nodes_phican(iphi) == phican) then  !(only proceed with evaluation if this node correspond to this canonical basis function 
    !!!!!!!!!!!!!!!!!
    ! Get grid displacements if cell where iphi belongs from the canoncal cell 
    !!!!!!!!!!!!!!!!!
    dui = nodes_dui(iphi) ! We assume that all Korobov nets are obtained for the same velocity discretization. 
    dvi = nodes_dvi(iphi) ! and use the same canonical cell. Otherwise, different nets will require different arrays    
    dwi = nodes_dwi(iphi) ! nodes_dwi
    !!!!!!!!!!!!!!!!
    ! now we add the shift component by component to the velocity and verify that this shift does not 
    ! through use outside fo the veloicity domain 
    iuxicell = ni + dui
    if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
     ivxicell = nj + dvi
     if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
      iwxicell = nl + dwi
      if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
       iuxi1cell = n1i + dui
       if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
        ivxi1cell = n1j + dvi
        if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
         iwxi1cell = n1l + dwi
         if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
            ! Now that we know that both shifts were successfull, we set prepare both f(xi) and f(xi1) 
            ! and add record to the collision integral
            xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc
            xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc
            fxi = sum(f(xi_j+1 : xi_j+dofc)*lagrarry)
            fxi1 = sum(f(xi1_j+1 : xi1_j+dofc)*lagrarry1)
            fc(iphi)=fc(iphi) + fxi*fxi1*Akor_val
            !
         end if
        end if  
       end if 
      end if 
     end if 
    end if
   !    
  end if
 !
 end do ! end loop in velocity nodes
!!!!
end do ! end loop in Akor records 
!!!!
fc=fc/nodes_gwts ! Add 8/Delta v^{j} / \omega_{i}  
!!
deallocate(lagrarry)
!
end subroutine EvalCollisionPeriodicAKorOpt_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EvalCollisionPeriodicAKorMixedTerms_DGV (f,fm,fc) 
!!
!!! This subroutine is to evaluate the bilinear collision operator using Korobov quadratures.
!!
!! THis version of the subroutine uses two distribution functions.
!! 
!! This subroutine assumes a uniform velocity grid and a periodic structure of the 
!! basis functions. Operator Akor is only evaluated for one basis function. The values 
!! of A for other basis functions are obtained by a transposition from operator A on 
!! the canonical cell. (The canonical cell should be selected at the center of the mesh 
!! or close to the center. The value of Akor_phi can be used to determine which basis function 
!! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!! 
!! This subroutine depends on a number of arrays that need to be prepared before the first time step. 
!! Specifically nodes_phicanII
!!  nodes_duiII,nodes_dviII,nodes_dwiII,nodes_phican
!!  cells_ugiII,cells_vgiII,cells_wgiII
!! and  
!!  AkorAllNets_shift 
!! need to be prepared. 
!! To set up the first bunch, call SetCellsUGI_DGV
!! To set up the second bunch, call ...
!! 
!!
!!!!!!!!!!!!!!!
!
subroutine EvalCollisionPeriodicAKorMixedTerms_DGV(f,fm,fc) 
!
use DGV_commvar, only: nodes_phican,nodes_dui,nodes_dvi,nodes_dwi,nodes_gwts, &
                       AkorAllNets, AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                       korob_net_paramAllNets, nodes_AkorshiftAllNets,numKornets, &
                       u_L,u_R,v_L,v_R,w_L,w_R,&
                       cells_ru,cells_lu,cells_rv,cells_lv,cells_rw,cells_lw
                       
use DGV_dgvtools_mod                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the local maxwellian
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: iphi ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: AKorshift,phicap ! scrap variables. AKorshift is to keep the address of the cell rigth before records in A that correspond to basis function phi. 
                               ! phicap keeps the number of records in A for the basis function phi
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: zz,xicell, xi1cell !indices of the cellscontaining xi and xi1
integer (I4B) :: k ! index of the korobov node
integer :: loc_alloc_stat ! variable to keep the allocation status
integer :: num_nets, curr_net ! a scrap variable to keep the number of all Korobov nets and the number fo the net that is used
real :: harvest ! a scrap variable
real (DP) :: frac_part,xiu,xiv,xiw,xi1u,xi1v,xi1w,du,dv,dw,fval,fval1,fmval,fmval1 ! scrap variables 
integer (I4B), dimension (7) :: korob_net_param ! scrap variable to keep the korobov parameters.
logical :: outsidedomain, outsidedomain1
integer (I4B) :: pcn, pcn1 ! scrap indices to keep numbers of cells for korobov points 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
num_nets = numKornets ! detemine the number of nets
!!!!!!!!!!!!!!!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!
call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!
do iphi=1,size(nodes_gwts,1)  ! Loop in the nodal points. Collision operator needs to be evaluated at each node
  !!!! now we randomly pick the number of the net that will be used for integration  
  call RANDOM_NUMBER(harvest)
  if (harvest >= 1.0D0) then !!! to avoid the result harvest  = 1 for sure
   harvest = .999 
  end if
  curr_net = int(harvest*Real(num_nets))+1 ! randomly select the net from 1...num_nets available nets (all supplied nets must be about the same accuracy...)
  !!!!!!! 
  Akorshift = nodes_AkorshiftAllNets(curr_net)%p(nodes_phican(iphi)) ! shift inside the Akor Array to get to correct basis function
  phicap = AkorAllNets_capphi(curr_net)%p(nodes_phican(iphi))               ! number of nonzero entires for this basis function 
  korob_net_param = korob_net_paramAllNets(curr_net)%p          ! save the parameters of the Korobov net into a temp array
  !!!!!!!!!!!!!!!!!
  dui = nodes_dui(iphi) ! We assume that all Korobov nets are obtained for the same velocity discretization. 
  dvi = nodes_dvi(iphi) ! and use the same canonical cell. Otherwise, different nets will require different arrays    
  dwi = nodes_dwi(iphi) ! nodes_dwi
  !!
  do j = 1,phicap  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  Determine the two velocity vectors in the Korobov Quadrature
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   k = AkorAllNets_k(curr_net)%p(j) ! this is the number of the Korobov node for 
   frac_part = Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP) 
   xiu = (frac_part - 0.5_DP)*(u_R-u_L)+(u_R+u_L)/2.0_DP 
   frac_part = Real( (korob_net_param(3)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(3)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xiv = (frac_part - 0.5_DP)*(v_R-v_L)+(v_R+v_L)/2.0_DP 
   frac_part = Real( (korob_net_param(4)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(4)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xiw = (frac_part - 0.5_DP)*(w_R-w_L)+(w_R+w_L)/2.0_DP 
   frac_part = Real( (korob_net_param(5)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(5)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1u = (frac_part - 0.5_DP)*(u_R-u_L)+(u_R+u_L)/2.0_DP 
   frac_part = Real( (korob_net_param(6)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(6)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1v = (frac_part - 0.5_DP)*(v_R-v_L)+(v_R+v_L)/2.0_DP 
   frac_part = Real( (korob_net_param(7)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
          Real(FLOOR(( Real((korob_net_param(7)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
   xi1w = (frac_part - 0.5_DP)*(w_R-w_L)+(w_R+w_L)/2.0_DP 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! now we need to shift velocity point to take into account the use of shift in the canonical function
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   xiu = xiu + dui*du;   xiv = xiv + dvi*dv;   xiw = xiw + dwi*dw; 
   xi1u = xi1u + dui*du; xi1v = xi1v + dvi*dv; xi1w = xi1w + dwi*dw; 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! now we interpolate the solution at velocity points (xiu,xiv,xiw) and (xi1u,xi1v,xi1w)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   outsidedomain = .false. 
   pcn  =  QuickCellFindUniformGrid_DGV(xiu,xiv,xiw,outsidedomain) ! we locate the cell that contains (xiu,xiv,xiw)
   outsidedomain1 = .false.                                ! if the velocity value outside domain, return 
   pcn1 =  QuickCellFindUniformGrid_DGV(xi1u,xi1v,xi1w,outsidedomain1) ! outsidedomain = .true.
   if (outsidedomain .or. outsidedomain1) then 
    cycle  ! points that are outsied of the mesh are ignored in the sum. 
   end if 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! interpolate the values of the function at (xiu,xiv,xiw) and (xiu1,xiv1,xiw1)
   !!! we compute the normalized values of the velocity on the cell: 
   fval  = EvalSolVeloPtCellFast_DGV(f,xiu,xiv,xiw,pcn)
   fval1 = EvalSolVeloPtCellFast_DGV(f,xi1u,xi1v,xi1w,pcn1)
   fmval = EvalSolVeloPtCellFast_DGV(fm,xiu,xiv,xiw,pcn)
   fmval1= EvalSolVeloPtCellFast_DGV(fm,xi1u,xi1v,xi1w,pcn1)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Add the new value to the integral
   fc(iphi) = fc(iphi)+(fval*fmval1+fval1*fmval)*AkorAllNets(curr_net)%p(Akorshift+j)   
  end do 
  fc(iphi)=fc(iphi)/nodes_gwts(iphi) ! the 1/P weight to put in front of the Korbov quadrature  
   !! The value of the collision integral for velocity node $iphi$ is computed
end do ! End of the main loop in nodal points
!!
end subroutine EvalCollisionPeriodicAKorMixedTerms_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EvalCollisionPeriodicAKorMxdTms_DGV (f,fc) 
!!
!! This is a version of the above subroutine with summations rearranged tryting to 
!! achieve optimaltiy 
!!
!!! This subroutine is to evaluate the bilinear collision operator using Korobov quadratures.
!!
!! This subroutine assumes a uniform velocity grid and a periodic structure of the 
!! basis functions. Operator Akor is only evaluated for one basis function. The values 
!! of A for other basis functions are obtained by a transposition from operator A on 
!! the canonical cell. (The canonical cell should be selected at the center of the mesh 
!! or close to the center. The value of Akor_phi can be used to determine which basis function 
!! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
!! 
!! This subroutine depends on a number of arrays that need to be prepared before the first time step. 
!! Specifically nodes_phicanII
!!  nodes_duiII,nodes_dviII,nodes_dwiII,nodes_phican
!!  cells_ugiII,cells_vgiII,cells_wgiII
!! and  
!!  AkorAllNets_shift 
!! need to be prepared. 
!! To set up the first bunch, call SetCellsUGI_DGV
!! To set up the second bunch, call ...
!! 
!!
!!!!!!!!!!!!!!!
!
subroutine EvalCollisionPeriodicAKorMxdTmsOpt_DGV(f,fm,fc) 
!
use DGV_commvar, only: nodes_phican,nodes_dui,nodes_dvi,nodes_dwi,nodes_gwts, &
                       AkorAllNets, AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                       korob_net_paramAllNets, nodes_AkorshiftAllNets,numKornets, &
                       u_L,u_R,v_L,v_R,w_L,w_R,g_nds_all, &
                       cells_ru,cells_lu,cells_rv,cells_lv,cells_rw,cells_lw,&
                       grids_cap_u,grids_cap_v,grids_cap_w, &
                       cells_gou,cells_gov,cells_gow,&
                       nodes_pcell
                       
use DGV_dgvtools_mod                       

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (in) :: fm ! the local maxwellian
real (DP), dimension (:), intent (out) :: fc ! the value of the collision operator for each component of the solution.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: iphi ! the index of the basis function and nodes
integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: j! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: AkorNrecs,n_ijl,xi_j,xi1_j,ni,nj,nl,n1i,n1j,n1l ! scrap variables to keep addresses of cells where the velocity point falls
integer (I4B) :: dui,dvi,dwi ! integer displacements from the cell where iphi is to the cell of the canonical node
integer (I4B) :: iuxicell,ivxicell,iwxicell,iuxi1cell,ivxi1cell,iwxi1cell ! numbers of the cell for translated xi and xi1
integer (I4B) :: ii,jj,ll,ixicell,ixi1cell ! indices of the nodes xi and xi1                                
integer (I4B) :: k ! index of the korobov node
integer :: loc_alloc_stat ! variable to keep the allocation status
integer :: num_nets, curr_net ! a scrap variable to keep the number of all Korobov nets and the number fo the net that is used
real :: harvest ! a scrap variable
real (DP) :: frac_part,xiu,xiv,xiw,xi1u,xi1v,xi1w,du,dv,dw ! scrap variables
real (DP) :: fxi,fxi1,fmxi,fmxi1,Akor_val ! scrap variables
real (DP), dimension(:), allocatable :: lagrarry, lagrarry1 
integer (I4B), dimension (7) :: korob_net_param ! scrap variable to keep the korobov parameters.
logical :: outsidedomain, outsidedomain1
real (DP) :: piunor, pjvnor, plwnor,unor,wnor,vnor !scrap values  
integer (I4B) :: pcn, phican ! scrap indices to keep numbers of cell for korobov points, and sracp to keep the index of the basis function 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
num_nets = numKornets ! detemine the number of nets
!!!!!!!!!!!!!!!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!!!!!!!!!!!!!!!
du=cells_ru(1)-cells_lu(1) ! the mesh is assumed uniform. Therefore values obtained from the 
dv=cells_rv(1)-cells_lv(1) ! first cell can be used for all cells 
dw=cells_rw(1)-cells_lw(1) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! IMPORTANT: We also assume that there is only one grid! 
pgcu=grids_cap_u(1)-1
pgcv=grids_cap_v(1)-1
pgcw=grids_cap_w(1)-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fc=0 ! nullify the result before computing... 
!!!!!!!!!!!!
!!! call RANDOM_SEED ! prepare the random number generator -- refresh the seed vector
!!!!!!!!!!!!
!call RANDOM_NUMBER(harvest)
!if (harvest >= 1.0D0) then !!! to avoid the result harvest  = 1 for sure
! harvest = .999 
!end if
!curr_net = int(harvest*Real(num_nets))+1 ! randomly select the net from 1...num_nets available nets (all supplied nets must be about the same accuracy...)
!!!  Override the above ...
curr_net = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
korob_net_param = korob_net_paramAllNets(curr_net)%p          ! save the parameters of the Korobov net into a temp array
!!! determine the total number of records in the entire Akor array:
AkorNrecs = size(AkorAllNets(curr_net)%p,1)
if (AkorNrecs /= sum(AkorAllNets_capphi(curr_net)%p)) then 
 print *,"EvalCollisionPeriodicAKorOpt_DGV: possible data error. Length of Akor not equal sum(Akor_capphi)" 
 stop 
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (lagrarry(1:gou*gov*gow),lagrarry1(1:gou*gov*gow), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalCollisionPeriodicAKorOpt_DGV: Error allocation arrays (largarray,lagrarry1)"
 stop
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do j=1,AkorNrecs  ! select  a record from the pre-computed Akor array
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !  Determine the two velocity vectors that correspond to this record 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Akor_val = AkorAllNets(curr_net)%p(j)
 k = AkorAllNets_k(curr_net)%p(j) ! this is the number of the Korobov node for 
 frac_part = Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(2)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP) 
 xiu = frac_part*(u_R-u_L) + u_L 
 frac_part = Real( (korob_net_param(3)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(3)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xiv = frac_part*(v_R-v_L) + v_L
 frac_part = Real( (korob_net_param(4)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(4)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xiw = frac_part*(w_R-w_L) + w_L 
 frac_part = Real( (korob_net_param(5)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(5)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1u = frac_part*(u_R-u_L) + u_L
 frac_part = Real( (korob_net_param(6)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(6)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1v = frac_part*(v_R-v_L) + v_L
 frac_part = Real( (korob_net_param(7)*k) ,DP)/Real(korob_net_param(1),DP)  - & 
        Real(FLOOR(( Real((korob_net_param(7)*k),DP)/Real(korob_net_param(1),DP) ),I4B),DP)
 xi1w = frac_part*(w_R-w_L) + w_L  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! For selected Korbobov node find were velocities fall on the uniform velocity grid  (not going to work on non-uniforom grids)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 outsidedomain = .false. 
 call QuickCellFindUniformGridIJL_DGV(xiu,xiv,xiw,ni,nj,nl,outsidedomain) ! we locate the cell that contains (xiu,xiv,xiw)
 outsidedomain1 = .false.                                ! if the velocity value outside domain, return 
 call QuickCellFindUniformGridIJL_DGV(xi1u,xi1v,xi1w,n1i,n1j,n1l,outsidedomain1) ! outsidedomain = .true.
 if (outsidedomain .or. outsidedomain1) then 
  print *,"EvalCollisionPeriodicAKorOpt_DGV: data error. Akor contains recs for nodes outside the vel. dom."
  stop  ! points that are outside of the mesh -- must be an error in data or a wrong velocity domain. 
 end if 
 !!! now we know local coordinates of cells where falls the selected korbov velocity node.    
 !!! Next find local coordinates of the canonical basis function for this records
 phican = AkorAllNets_phi(curr_net)%p(j) ! number of the basis function for which entry of Akor was computed 
 !!!! Now we prepare for interpolation of the solution on the shifted korobov nodes... 
 !!! first triple of the node: 
 !!! compute the number of the cell where the node belongs
 pcn = (ni-1)*pgcw*pgcv + (nj-1)*pgcw + nl
 !!! Next we evaluate nodal DG basis fucntion on the cell where the velocities fall. We evaluate then at every nodal point in 
 !!! velcity on that cell and save in an array. This array will be used to later to interpolate the value of the function 
 !!! on each cell. The algorithm depends on the fact that the velocity grid is uniform.   
 unor = (xiu - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
 vnor = (xiv - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
 wnor = (xiw - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
 ! basis function for that node and add it to the interpolated value
 lagrarry=0
 n_ijl = 1
 do ii=1,gou
  piunor = lagrbasfun(ii,unor,g_nds_all(:gou,gou))
  do jj=1,gov
   pjvnor = lagrbasfun(jj,vnor,g_nds_all(:gov,gov))
   do ll=1,gow 
    ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
    ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
    plwnor = lagrbasfun(ll,wnor,g_nds_all(:gow,gow))
    ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
    lagrarry(n_ijl) = piunor*pjvnor*plwnor
    n_ijl = n_ijl+1
   enddo
  enddo 
 enddo  
 !!! second triple of the Korobov node: 
 !!! compute the number of the cell where the node belongs
 pcn = (n1i-1)*pgcw*pgcv + (n1j-1)*pgcw + n1l
 !!! Next we evaluate nodal DG basis fucntion on the cell where the velocities fall. We evaluate then at every nodal point in 
 !!! velcity on that cell and save in an array. This array will be used to later to interpolate the value of the function 
 !!! on each cell. The algorithm depends on the fact that the velocity grid is uniform.   
 unor = (xiu - (cells_ru(pcn) + cells_lu(pcn))/2.0_DP )/(cells_ru(pcn) - cells_lu(pcn))*2.0_DP 
 vnor = (xiv - (cells_rv(pcn) + cells_lv(pcn))/2.0_DP )/(cells_rv(pcn) - cells_lv(pcn))*2.0_DP 
 wnor = (xiw - (cells_rw(pcn) + cells_lw(pcn))/2.0_DP )/(cells_rw(pcn) - cells_lw(pcn))*2.0_DP 
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! next we will go over all velocity nodes on the primary cell. If we find a node that belongs to the cell with number (primecellnum) we will assemble the 
 ! basis function for that node and add it to the interpolated value
 lagrarry1=0
 n_ijl = 1
 do ii=1,gou
  piunor = lagrbasfun(ii,unor,g_nds_all(:gou,gou))
  do jj=1,gov
   pjvnor = lagrbasfun(jj,vnor,g_nds_all(:gov,gov))
   do ll=1,gow 
    ! next we need to know the three local indices that tell what velocity nodal values correspond to this 
    ! basis function. this is also simple since this information is also stored in the Nodes Arrays.
    plwnor = lagrbasfun(ll,wnor,g_nds_all(:gow,gow))
    ! now y contains the value of the basis function for the node "j". It is time to add the node J to interpolation: 
    lagrarry1(n_ijl) = piunor*pjvnor*plwnor
    n_ijl = n_ijl+1
   enddo
  enddo 
 enddo  
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 ! We are ready to add this integration node to the collision integral 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!
 do iphi=1,size(nodes_gwts,1)  ! Loop in the nodal points. Collision operator needs to be evaluated at each node
   if (nodes_phican(iphi) == phican) then  !(only proceed with evaluation if this node correspond to this canonical basis function 
    !!!!!!!!!!!!!!!!!
    ! Get grid displacements if cell where iphi belongs from the canoncal cell 
    !!!!!!!!!!!!!!!!!
    dui = nodes_dui(iphi) ! We assume that all Korobov nets are obtained for the same velocity discretization. 
    dvi = nodes_dvi(iphi) ! and use the same canonical cell. Otherwise, different nets will require different arrays    
    dwi = nodes_dwi(iphi) ! nodes_dwi
    !!!!!!!!!!!!!!!!
    ! now we add the shift component by component to the velocity and verify that this shift does not 
    ! through use outside fo the veloicity domain 
    iuxicell = ni + dui
    if ((iuxicell >= 1) .and. (iuxicell <= pgcu)) then ! check if xi outside of bounds in u
     ivxicell = nj + dvi
     if ((ivxicell >= 1) .and. (ivxicell <= pgcv)) then ! check if xi outside of bounds in v 
      iwxicell = nl + dwi
      if ((iwxicell >= 1) .and. (iwxicell <= pgcw)) then ! check if xi outside of bounds in w	
       iuxi1cell = n1i + dui
       if ((iuxi1cell >= 1) .and. (iuxi1cell <= pgcu)) then ! check if xi1 outside of bounds in u
        ivxi1cell = n1j + dvi
        if ((ivxi1cell >= 1) .and. (ivxi1cell <= pgcv)) then ! check if xi1 outside of bounds in v 
         iwxi1cell = n1l + dwi
         if ((iwxi1cell >= 1) .and. (iwxi1cell <= pgcw)) then ! check if xi1 outside of bounds in w	
            ! Now that we know that both shifts were successfull, we set prepare both f(xi) and f(xi1) 
            ! and add record to the collision integral
            xi_j = ((iuxicell-1)*pgcw*pgcv + (ivxicell-1)*pgcw + iwxicell-1)*dofc
            xi1_j = ((iuxi1cell-1)*pgcw*pgcv + (ivxi1cell-1)*pgcw + iwxi1cell-1)*dofc
            fxi = sum(f(xi_j+1 : xi_j+dofc)*lagrarry)
            fmxi = sum(fm(xi_j+1 : xi_j+dofc)*lagrarry)
            fxi1 = sum(f(xi1_j+1 : xi1_j+dofc)*lagrarry1)
            fmxi1 = sum(fm(xi1_j+1 : xi1_j+dofc)*lagrarry1)
            fc(iphi)=fc(iphi) + (fxi*fmxi1+fmxi*fxi1)*Akor_val
            !
         end if
        end if  
       end if 
      end if 
     end if 
    end if
   !    
  end if
 !
 end do ! end loop in velocity nodes
!!!!
end do ! end loop in Akor records 
!!!!
fc=fc/nodes_gwts ! Add 8/Delta v^{j} / \omega_{i}  
!!
deallocate(lagrarry)
!
end subroutine EvalCollisionPeriodicAKorMxdTmsOpt_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right hand side calculation for the ES-BGK distribution and collision operator
! Here, the RHS = nu * (f0 - f)
! where
! nu = collision frequency
! f0 = ESBGK distribution
! f = solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColESBGK(f,RHS)

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_commvar, only: nodes_gwts, nodes_u, nodes_v, nodes_w, &
				   alpha, gas_viscosity, gas_T_reference, gas_alpha, C_inf, gasR

real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: RHS ! the value of the collision operator for each component of the solution.

real (DP) :: u_0, v_0, w_0 ! bulk velocities
real (DP) :: n ! density
real (DP) :: T ! temperature
!real (DP) :: Determinant ! the determinant of the Tensor
real (DP) :: Pressure
real (DP) :: nu ! this is the collision frequency term
real (DP), dimension (:), allocatable :: f0 ! Distribution function

real (DP), parameter :: kBoltzmann = 1.3806503D-23
real (DP), dimension (3,3) :: Tensor, TensorInv ! both the tensor and the inverse tensor for the ES-BGK
real (DP) :: Determinant ! the determinant of the tensor
integer :: loc_alloc_stat

! compute the macroparameters for use in the tensor computation
!!!!!!!!
! density (number density)
n=sum(f*nodes_gwts)
!!!!!!!!
! momentum 
u_0=sum(f*nodes_gwts*nodes_u)/n
v_0=sum(f*nodes_gwts*nodes_v)/n
w_0=sum(f*nodes_gwts*nodes_w)/n
!!!!!!!!
! temperature
T = sum(f*nodes_gwts*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!

! the tensor and its corresponding inverse and determinant is computed here
call MakeTensor(f,Tensor) ! alpha only in common variabls
TensorInv = inv(Tensor)
Determinant = DetTensor(Tensor)

allocate (f0(1:size(f,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalColESBGK: Allocation error for variables (f0)"
 stop
end if

f0 = ESBGK_f0(TensorInv,Determinant,n,u_0,v_0,w_0,nodes_u,nodes_v,nodes_w)

! now to evaluate the collision requency term
Pressure = n*T ! dimensionless Pressure is computed here
nu = Pressure/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/T)**gas_alpha ! final dimensionless nu ! have gas_T_reference be dimensionless?
RHS = nu * (f0 - f)

deallocate (f0)
!
end subroutine EvalColESBGK



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right hand side calculation for the Shakhov model of the collision operator
! Here, the RHS = nu * (f0^s - f)
! where
! nu = collision frequency = P/mu
! f0^s = Shakhov model of the equilibrium distribution
! f = solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColShakov(f,RHS)
use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_commvar, only: nodes_gwts, nodes_u, nodes_v, nodes_w, &
				   alpha, gas_viscosity, gas_T_reference, gas_alpha, C_inf, gasR

real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: RHS ! the value of the collision operator for each component of the solution.

real (DP) :: u_0, v_0, w_0 ! bulk velocities
real (DP) :: n ! density
real (DP) :: T ! temperature
!real (DP) :: Determinant ! the determinant of the Tensor
real (DP) :: Pressure
real (DP) :: nu ! this is the collision frequency term
real (DP), dimension (:), allocatable :: f0 ! Distribution function
real (DP), dimension (:), allocatable :: S_u,S_v,S_w ! components of the vector S in the Shakhov model 
integer :: loc_alloc_stat
real (DP), parameter :: kBoltzmann = 1.3806503D-23

! compute the macroparameters for use in the tensor computation
!!!!!!!!
! density (number density)
n=sum(f*nodes_gwts)
!!!!!!!!
! momentum 
u_0=sum(f*nodes_gwts*nodes_u)/n
v_0=sum(f*nodes_gwts*nodes_v)/n
w_0=sum(f*nodes_gwts*nodes_w)/n
!!!!!!!!
! temperature
T = sum(f*nodes_gwts*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!!!!!!!!!


allocate (f0(1:size(f,1)),S_u(1:size(f,1)),S_v(1:size(f,1)),S_w(1:size(f,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalColShakov: Allocation error for variables (f0,S_u,S_v,S_w)"
 stop
end if

! vector S from the Shakhov model
S_u=sum(f*nodes_gwts*(nodes_u-u_0)*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/(T**(3.0_DP/2.0_DP))
S_v=sum(f*nodes_gwts*(nodes_v-v_0)*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/(T**(3.0_DP/2.0_DP))
S_w=sum(f*nodes_gwts*(nodes_w-w_0)*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/(T**(3.0_DP/2.0_DP))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we compute the Shakhov equilibrium distribution:
f0 = Shakhov_f0(alpha,n,u_0,v_0,w_0,T,S_u,S_v,S_w,nodes_u,nodes_v,nodes_w)
! now to evaluate the collision requency term
Pressure = n*T ! dimensionless Pressure is computed here
nu = Pressure*((gas_T_reference/C_inf**2*2.0d0*gasR)/T)**gas_alpha ! final dimensionless nu ! have gas_T_reference be dimensionless?

RHS = nu * (f0 - f);

deallocate (f0,S_u,S_v,S_w)

end subroutine EvalColShakov

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Original coding C. Euler, heavily modified by Alex Alekseenko
!
! Collision operator implementeing the ES-BGK distribution with velocity dependent collision 
! frequency 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColVelES(f,Df,fcol,time)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w, &
                   !!!!! Diagnostic ADDED THIS TO save the coefficinets on file. Later remove! 
                   Cco 
                   !!!!! END DIAGNOSTIC 
                   
intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar ! scrap variables to keep the values of the local bulk velocity
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
integer :: loc_alloc_stat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab) ! the subroutine returns moment relaxation rates
											! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu
!!!!!!!!!!!!!!!!!!!!!
! We will now attempt to update the cofficients. For that the Matrix Mom has too small of a L_1 matrix norm, or if
! Mom fails to be inverted, we will not use this model. Instead, we will call the default "fall back" model (ES-BGK currently) 
! If Mom is small, the linear system we solve becoms degenerate. Since in this case the coefficients are not accurate anyway, 
! we will just replace our model with a fall back model. Also, we return to fall back model if no reliable relaxation rates could be obtained. 
! A typical situation when mom is small or rates can not be obtained is when the solution is near continuum. So, we recommend to use the model 
! only when a significant deviation from the local Maxwellian happens. 
	! to determine the coefficients of the velocity dependent collision frequency, we solve the system
	! b = Mom*C 
	! C -- is the vector of coefficients, i.e. C = array Cco
	! Mom is the matrix obtained multiplying the BGK-collision term with velocity dependent collision frequency by the kernel of the enforced moment and integrating) -- 
	! Mom is the synonim of the stiffness matrix. 
	! b = Dphi correcpond to the desired derivatives of the enforced moments: 
	! b is obtined by taking the relaxation rate for the moment and mupliplying it by the difference betweenthe moment and its "equilibrium" state. The 
	! equilibrium state is computed by taking the moment of the local maxwellian. Well, actually, the difference in the moment is computed by taking the moment of Df  
	! See additional detail in the method description (Alekseenko Euler, the BGK model with velocity dependent collision frequency). 
	! Let us get started. First we obtain the stifness matrix:
call setMom_DGV (-Df,Mom,DifMom,ubar,vbar,wbar) ! psi (nu), phi (moments)
! Calculate \| Mom \|_{L_1}
mom_L1=MAXVAL(SUM(ABS(Mom),1))
! Next,check if the mom_l1 is large enough. If it is and the enforced moments relaxation rates are reliable, then the Cco vector is updated 
! and the model is used. Otherwise the fallback model is used. 
if ((mom_L1>Mom_trshld) .and. (momsrates_reliab)) then  
 !calculate the b_i terms
 ! ATTENTION: If you desire to enforce conservation of some moments then you should set the correcponsind components of Dphi to zero. 
 Dphi = (momsrates-nu)*DifMom
 ! This subroutine uses the same number of enforced moments and the basis functions in the representation of the 
 ! In the future we will develop more elaborate subroutine. Right now, cases (Order_nu .neq. Order) are not considered
 if (Order_nu == Order) then 
  MomInv = inv(Mom) ! get the inverse of the moments matrix
  ! solve for c_i
  call DGEMM('N', 'N', Order, 1, Order, 1.0d0, MomInv, Order, Dphi, Order, 0.0d0, Cco_temp, Order) ! Cco = MomInv*Dphi
  ! now Cco contain the desired coefficients of the velocity dependent collision frequency (VDCF)
 else if (Order_nu < Order) then 
  ! solve linear least squares to determine coefficients
  Cco_temp = linleastsqrs(Mom,Dphi) 
  ! now Cco_temp contain the desired coefficients of the velocity dependent collision frequency (VDCF)
 else 
  print *, "EvalColVelES: Error: The number of basis functions in VDCF is bigger than the number of enforced moments. Stop"  
  stop
 end if 
 ! the new Cco have beem computed, now let us accemple the velocity dependent collision frequency   
 allocate (nuBolt(1:size(f,1)), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "EvalColVelES: Allocation error for variables (nuBolt)"
  stop
 end if
 !!!!!!!!!!!!! DIAGNOSTICK for PRAKASH !!!!!!!!! save values of Cco into global variable REMOVE later
  Cco=Cco_temp 
 !!!!!!!!!!!!! END DIAGNOSTIC !!!!!!!!!!!!!!!
 call getnu(nu,Cco_temp,nodes_u,nodes_v,nodes_w,nuBolt) ! psi (nu)   
 ! finally, we get the RHS
 fcol = -nuBolt*Df
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! DIAGNOSTICS  
 ! do i=1,size(nuBolt,1)
 !  nuBolt(i) = maxval(nuBolt(i),0.0d0)
 ! end do 
 min_nu = minval(nuBolt) 
 max_nu = maxval(nuBolt)
 if (min_nu < 0.0) then 
  PRINT *, "EvalColVelES: collision frequency is negative at least at one point"
 end if 
 !!! END DIAGNOSTIC
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 deallocate (nuBolt)
 ! The BGK model collision operator with velocity-dependent collision frequency have been evaluated. 
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 PRINT *, "EvalColVelES: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab", mom_L1,Mom_trshld,momsrates_reliab
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
if (temp_dp > 0.0d0) then 
  entr = (-1)*sum(Real(log(f),DP)*fcol)
  if (entr< 0.0) then 
   PRINT *, "EvalColVelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColVelES: entropy:", entr 
  end if 
else 
 PRINT *, "EvalColVelES: velocity distribution is negative at least at one point"
 entr = 0.0d0 
 do i=1,size(f,1)
  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
 end do
 if (entr< 0.0) then 
   PRINT *, "EvalColVelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColVelES: entropy:", entr 
  end if 
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColVelES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColVelES1Donecell
!
! This is a copy of the above subroutine adjusted to be used with solutions on mupliple 
! spatial cells. One cell at a time is evaluates. The number os the spatial cell is passed in the variable cellnum
!
! this subroutine to be used with the multidimensional applications
! 
! Collision operator implementeing the ES-BGK distribution with velocity dependent collision 
! frequency 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColVelES1Donecell(f,Df,fcol,time,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                       nodes_gwts,nodes_u,nodes_v,nodes_w,Cco1D

intrinsic MAXVAL, ABS


real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: nu ! scrap variable to keep the default collision frequency
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
real (DP) :: nden,ubar,vbar,wbar,temp ! scrap variable to keep the macroparameter of pressure,number density and temperature
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable  :: nuBolt ! the values of the velocity dependent collision frequency 
!!!!
integer (I4B) :: i ! scrap index 
integer :: loc_alloc_stat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates in this cell. 
					!If true then moments rates are reliable and can be used 
                    !for the evaluation of velocity dependent collision frequency 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab,cellnum) ! the subroutine returns coefficients nucoefs and also te value of the L1_norm of the difference between the soltuion and the local maxwellian
!!!!!!!!!!!!!!!!!!!!!
! We will now attempt to update the cofficients. For that the Matrix Mom has too small of a L_1 matrix norm, or if
! Mom fails to be inverted, we will not use this model. Instead, we will call the default "fall back" model (ES-BGK currently) 
! If Mom is small, the linear system we solve becoms degenerate. Since in this case the coefficients are not accurate anyway, 
! we will just replace our model with a fall back model. Also, we return to fall back model if no reliable relaxation rates could be obtained. 
! A typical situation when mom is small or rates can not be obtained is when the solution is near continuum. So, we recommend to use the model 
! only when a significant deviation from the local Maxwellian happens. 
	! to determine the coefficients of the velocity dependent collision frequency, we solve the system
	! b = Mom*C 
	! C -- is the vector of coefficients, i.e. C = array Cco
	! Mom is the matrix obtained multiplying the BGK-collision term with velocity dependent collision frequency by the kernel of the enforced moment and integrating) -- 
	! Mom is the synonim of the stiffness matrix. 
	! b = Dphi correcpond to the desired derivatives of the enforced moments: 
	! b is obtined by taking the relaxation rate for the moment and mupliplying it by the difference betweenthe moment and its "equilibrium" state. The 
	! equilibrium state is computed by taking the moment of the local maxwellian. Well, actually, the difference in the moment is computed by taking the moment of Df  
	! See additional detail in the method description (Alekseenko Euler, the BGK model with velocity dependent collision frequency). 
	! Let us get started. First we obtain the stifness matrix:
call setMom_DGV (-Df,Mom,DifMom,ubar,vbar,wbar) ! psi (nu), phi (moments)
! Calculate \| Mom \|_{L_1}
mom_L1=MAXVAL(SUM(ABS(Mom),1))
! Next,check if the mom_l1 is large enough. If it is and the enforced moments relaxation rates are reliable, then the Cco vector is updated 
! and the model is used. Otherwise the fallback model is used. 
if ((mom_L1>Mom_trshld) .and. (momsrates_reliab)) then  
 !calculate the b_i terms
 ! ATTENTION: If you desire to enforce conservation of some moments then you should set the correcponsind components of Dphi to zero. 
 Dphi = (momsrates-nu)*DifMom
 ! This subroutine uses the same number of enforced moments and the basis functions in the representation of the 
 ! In the future we will develop more elaborate subroutine. Right now, cases (Order_nu .neq. Order) are not considered
 if (Order_nu == Order) then 
  MomInv = inv(Mom) ! get the inverse of the moments matrix
  ! solve for c_i
  call DGEMM('N', 'N', Order, 1, Order, 1.0d0, MomInv, Order, Dphi, Order, 0.0d0, Cco_temp, Order) ! Cco_temp = MomInv*Dphi
  ! now Cco_temp contain the desired coefficients of the velocity dependent collision frequency (VDCF)
 else if (Order_nu < Order) then 
  ! solve linear least squares to determine coefficients
  Cco_temp = linleastsqrs(Mom,Dphi) 
  ! now Cco_temp contain the desired coefficients of the velocity dependent collision frequency (VDCF)
 else 
  print *, "EvalColVelES1Donecell: The number of basis functions in VDCF is bigger than the number of enforced moments. Stop"  
  stop
 end if 
 ! the new Cco have beem computed, now let us accemple the velocity dependent collision frequency   
 allocate (nuBolt(1:size(f,1)), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "EvalColVelES1Donecell: Allocation error for variables (nuBolt)"
  stop
 end if
 call getnu(nu,Cco_temp,nodes_u,nodes_v,nodes_w,nuBolt) ! psi (nu)   
 ! finally, we get the RHS
 fcol = -nuBolt*Df
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! DIAGNOSTICS  
 ! do i=1,size(nuBolt,1)
 !  nuBolt(i) = max(nuBolt(i),0.0d0)
 ! end do 
 ! min_nu = minval(nuBolt)
 ! max_nu = maxval(nuBolt)
 ! if (min_nu<0.0) then 
 !  PRINT *, "collision frequency is negative at least at one point"
 ! end if 
 !!! END DIAGNOSTIC
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 deallocate(nuBolt)
 ! The BGK model collision operator with velocity-dependent collision frequency have been evaluated. 
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 !!!!!!!!!!!!!!!!!!!!!
 ! DEBUG: comment the next two lines when no tracking of fall back model is desired 
 PRINT *, "EvalColVelES1Donecell: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab, cellnum", mom_L1,Mom_trshld,&
                                 momsrates_reliab, cellnum
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                 
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: can comment if not wanted - tells if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!temp_dp=minval(f)
!if (temp_dp > 0.0d0) then 
!  entr = (-1)*sum(log(f)*fcol)
!  if (entr< 0.0) then 
!   PRINT *, "entropy is negative!:", entr 
!  else  
!   PRINT *, "entropy:", entr 
!  end if 
!else 
! PRINT *, "velocity distribution is negative at least at one point"
! entr = 0.0d0 
! do i=1,size(f,1)
!  entr = entr - log(max(f(i),0.0000001d0))*fcol(i)
! end do
! if (entr< 0.0) then 
!   PRINT *, "entropy is negative!:", entr 
!  else  
!   PRINT *, "entropy:", entr 
!  end if 
!end if  
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColVelES1Donecell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColRelES(f,Df,fcol,time)
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T
!
! 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColRelES(f,Df,fcol,time)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w, &
                   !!!!! Diagnostic ADDED THIS TO save the coefficinets on file. Later remove! 
                   Cco 
                   !!!!! END DIAGNOSTIC 
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: Theta, Theta_inv ! scrap variable to keep the artificial stress tensor and its inverse 
real (DP), dimension(1:6) :: scrtheta ! 
real (DP), dimension(1:6) :: fij_arr, dfij_arr !Just used to capture the fij and dfij in the bugger to copy over to excel.
real (DP) :: fij, dfij,dens ! scrap variables
! 
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar ! scrap variables to keep the values of the local bulk velocity
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab) ! the subroutine returns moment relaxation rates
											! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 dens=sum(f*nodes_gwts)
 do i=1,6
   dfij = sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij = sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij_arr(i) = fij
   dfij_arr(i) = dfij
   scrtheta(i) = (fij - (momsrates(i)/nu)*dfij)*2.0_DP/dens
 end do 
 ! now the components of artificial stress tensor have been computed. We will transfer the values into the artificial stress tensor
 Theta(1,1) = scrtheta(1)
 Theta(2,2) = scrtheta(2)
 Theta(3,3) = scrtheta(3)
 Theta(2,1) = scrtheta(4)
 Theta(1,2) = scrtheta(4)
 Theta(3,1) = scrtheta(5)
 Theta(1,3) = scrtheta(5)
 Theta(3,2) = scrtheta(6)
 Theta(2,3) = scrtheta(6)
 ! we have got the artificial stress tensor set up 
 Theta_inv = inv(Theta) ! get the inverse of the moments matrix
 !!!!!!!!!!!!!!!!!!
 fij = Theta(1,1)*Theta(2,2)*Theta(3,3)+Theta(2,1)*Theta(3,2)*Theta(1,3)+Theta(3,1)*Theta(1,2)*Theta(2,3) - &
       Theta(1,3)*Theta(2,2)*Theta(3,1)-Theta(2,3)*Theta(3,2)*Theta(1,1)-Theta(3,3)*Theta(1,2)*Theta(2,1) ! this will keep the value of the determinant
 ! finally, we get the RHS
 fcol = nu*(ESBGK_f0 (Theta_inv,fij,dens,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w) - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 PRINT *, "EvalColRelES: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab", mom_L1,Mom_trshld,momsrates_reliab
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
if (temp_dp > 0.0d0) then 
  entr = (-1)*sum(Real(log(f),DP)*fcol)
  if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
else 
 PRINT *, "EvalColRelES: velocity distribution is negative at least at one point"
 entr = 0.0d0 
 do i=1,size(f,1)
  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
 end do
 if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColRelES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColRelES_1Donecell(f,Df,fcol,time)
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T
!
! 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColRelES_1Donecell(f,Df,fcol,time,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: Theta, Theta_inv ! scrap variable to keep the artificial stress tensor and its inverse 
real (DP), dimension(1:6) :: scrtheta ! 
real (DP), dimension(1:6) :: fij_arr, dfij_arr !Just used to capture the fij and dfij in the bugger to copy over to excel.
real (DP) :: fij, dfij,dens ! scrap variables
! 
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar ! scrap variables to keep the values of the local bulk velocity
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab,cellnum) ! the subroutine returns momrates
                                           	! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 dens=sum(f*nodes_gwts)
 do i=1,6
   dfij = sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij = sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij_arr(i)=fij
   dfij_arr(i)=dfij
   scrtheta(i) = (fij - (momsrates(i)/nu)*dfij)*2.0_DP/dens
 end do 
 ! now the components of artificial stress tensor have been computed. We will transfer the values into the artificial stress tensor
 Theta(1,1) = scrtheta(1)
 Theta(2,2) = scrtheta(2)
 Theta(3,3) = scrtheta(3)
 Theta(2,1) = scrtheta(4)
 Theta(1,2) = scrtheta(4)
 Theta(3,1) = scrtheta(5)
 Theta(1,3) = scrtheta(5)
 Theta(3,2) = scrtheta(6)
 Theta(2,3) = scrtheta(6)
 ! we have got the artificial stress tensor set up 
 Theta_inv = inv(Theta) ! get the inverse of the moments matrix
 !!!!!!!!!!!!!!!!!!
 fij = Theta(1,1)*Theta(2,2)*Theta(3,3)+Theta(2,1)*Theta(3,2)*Theta(1,3)+Theta(3,1)*Theta(1,2)*Theta(2,3) - &
       Theta(1,3)*Theta(2,2)*Theta(3,1)-Theta(2,1)*Theta(1,2)*Theta(3,3)-Theta(1,1)*Theta(3,2)*Theta(2,3) ! this will keep the value of the determinant
 ! finally, we get the RHS
 fcol = nu * (ESBGK_f0 (Theta_inv,fij,dens,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w) - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
!$OMP CRITICAL
 PRINT *, "EvalColRelES1Donecell: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab, cellnum", mom_L1,Mom_trshld,&
                                 momsrates_reliab, cellnum
!$OMP END CRITICAL
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
!$OMP CRITICAL
!if (temp_dp > 0.0d0) then 
!  entr = (-1)*sum(Real(log(f),DP)*fcol)
!  if (entr< 0.0) then 
!   PRINT *, "EvalColRelES: entropy is negative!:", entr 
!  else  
!   PRINT *, "EvalColRelES: entropy:", entr 
!  end if 
!else 
! PRINT *, "EvalColRelES: velocity distribution is negative at least at one point"
! entr = 0.0d0 
! do i=1,size(f,1)
!  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
! end do
! if (entr< 0.0) then 
!   PRINT *, "EvalColRelES: entropy is negative!:", entr 
!  else  
!   PRINT *, "EvalColRelES: entropy:", entr 
!  end if 
!end if  
!$OMP END CRITICAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColRelES_1Donecell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColRelES_SC(f,Df,fcol,time)
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T
!
! ATTENTION: MAKE SURE THAT function kernls_enfrsd_moms set up to evaluate 
! components of the stress tensor in the order t_11, t_22, t_33, t_12, t_13, t_23 
!
! This is a modification of the above model with some features introduced to enforce stability 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColRelES_SC(f,Df,fcol,time)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w, &
                   !!!!! Diagnostic ADDED THIS TO save the coefficinets on file. Later remove! 
                   Cco 
                   !!!!! END DIAGNOSTIC 
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: Theta, Theta_inv ! scrap variable to keep the artificial stress tensor and its inverse 
real (DP), dimension(1:6) :: scrtheta ! 
real (DP), dimension(1:6) :: fij_arr, dfij_arr !Just used to capture the fij and dfij in the bugger to copy over to excel.
real (DP) :: fij, dfij,dens ! scrap variables
! 
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar,maxfactor,factor,mvmaxnorm,en_corr ! scrap variables to keep the values of the local bulk velocity
real (DP), parameter :: Stress_TOL=0.5_DP ! The maximum allowed relative perturbation component vise for the corrections of the stress tensor
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab) ! the subroutine returns moment relaxation rates
											! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 dens=sum(f*nodes_gwts)
 mvmaxnorm=0
 ! here we compute the components of the stress tensor and the correction stress tensor
 do i=1,6
   dfij_arr(i) = sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij_arr(i) = sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   !! Compute the max norm of stress tensor, treating matrix as a vector.
   if (abs(fij_arr(i))>mvmaxnorm) then 
    mvmaxnorm=abs(fij_arr(i))
   end if
 end do
 !!! here we compute the correction to enforce conservation of temperature 
 !!! this correction will be subtracted on the diagonal of the RelES stress tensor
 en_corr = (momsrates(1)*dfij_arr(1)+momsrates(2)*dfij_arr(2)+momsrates(3)*dfij_arr(3))/3.0_DP
 !!! here we estimate the component-wise perturbation of the proposed corrections.  
 !!! maxfactor will measure the lagest correction relative to the vector matrix norm
 !!! This factor is 1 by default and does not affect the model. However, if the corrections are large
 !!! they will be scaled down to keep the target distribution locally supported. The scaling factor 
 !!! is added back in the last step. The controlled moments are not affected, but the 
 !!! rest of the momets will be forced to relax faster than the ES_BGK collision frequency
 !!!  
 maxfactor = 1.0_DP
 do i=1,3
   dfij= abs((momsrates(i)*dfij_arr(i) - en_corr)/nu)
   if (dfij > Stress_TOL*mvmaxnorm) then 
    factor = dfij/Stress_TOL*mvmaxnorm
    if (factor > maxfactor) then 
      maxfactor = factor
    end if
   end if    
 end do 
 do i=4,6
   dfij= abs((momsrates(i)/nu)*dfij_arr(i))
   if (dfij > Stress_TOL*mvmaxnorm) then 
    factor = dfij/Stress_TOL*mvmaxnorm
    if (factor > maxfactor) then 
      maxfactor = factor
    end if
   end if    
 end do 
 !!! BEGIN Diagnostics comment for production runs
 if ((maxfactor > 4.0) .and. (maxfactor <=16.0)) then 
  print *, "EvalColRelES_SC_1Donecell: Attention! correction factor exceeded 4" 
 end if 
 if (maxfactor > 16.0) then 
  print *, "EvalColRelES_SC_1Donecell: Attention! correction factor exceeded 16" 
 end if 
 !!! END diagnostics
 do i=1,3
  scrtheta(i) = (fij_arr(i) - (momsrates(i)*dfij_arr(i) - en_corr)/nu/maxfactor)*2.0_DP/dens
 end do 
 do i=4,6
  scrtheta(i) = (fij_arr(i) - (momsrates(i)/nu/maxfactor)*dfij_arr(i))*2.0_DP/dens
 end do 
 ! now the components of artificial stress tensor have been computed. We will transfer the values into the artificial stress tensor
 Theta(1,1) = scrtheta(1)
 Theta(2,2) = scrtheta(2)
 Theta(3,3) = scrtheta(3)
 Theta(2,1) = scrtheta(4)
 Theta(1,2) = scrtheta(4)
 Theta(3,1) = scrtheta(5)
 Theta(1,3) = scrtheta(5)
 Theta(3,2) = scrtheta(6)
 Theta(2,3) = scrtheta(6)
 ! we have got the artificial stress tensor set up 
 Theta_inv = inv(Theta) ! get the inverse of the moments matrix
 !!!!!!!!!!!!!!!!!!
 fij = Theta(1,1)*Theta(2,2)*Theta(3,3)+Theta(2,1)*Theta(3,2)*Theta(1,3)+Theta(3,1)*Theta(1,2)*Theta(2,3) - &
       Theta(1,3)*Theta(2,2)*Theta(3,1)-Theta(2,3)*Theta(3,2)*Theta(1,1)-Theta(3,3)*Theta(1,2)*Theta(2,1) ! this will keep the value of the determinant
 ! finally, we get the RHS
 fcol = nu*maxfactor*(ESBGK_f0 (Theta_inv,fij,dens,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w) - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 PRINT *, "EvalColRelES_SC_1Donecell: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab", mom_L1,Mom_trshld,momsrates_reliab
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
if (temp_dp > 0.0d0) then 
  entr = (-1)*sum(Real(log(f),DP)*fcol)
  if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
else 
 PRINT *, "EvalColRelES: velocity distribution is negative at least at one point"
 entr = 0.0d0 
 do i=1,size(f,1)
  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
 end do
 if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColRelES_SC


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColRelES_SC_1Donecell(f,Df,fcol,time)
!
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T
! 
! This is a modification of the above model with some features introduced to enforce stability 
! 
!
! ATTENTION: MAKE SURE THAT function kernls_enfrsd_moms set up to evaluate 
! components of the stress tensor in the order t_11, t_22, t_33, t_12, t_13, t_23 
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColRelES_SC_1Donecell(f,Df,fcol,time,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: Theta, Theta_inv ! scrap variable to keep the artificial stress tensor and its inverse 
real (DP), dimension(1:6) :: scrtheta ! 
real (DP), dimension(1:6) :: fij_arr, dfij_arr !Just used to capture the fij and dfij in the bugger to copy over to excel.
real (DP) :: fij, dfij,dens ! scrap variables
! 
real (DP) :: mom_L1 ! scrap variable to keep L1-norm of matrix Mom
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
real (DP), parameter :: Mom_trshld = 1.0d-8 ! This constant determines when the coefficents are updated in the velocity dependent collisio nfrequency. 
                     ! if \| Mom \|_{L1}> Mom_trshld then the coefficients are updated. Otherwise, the "fall back" model is used.  
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar,maxfactor,factor,mvmaxnorm,en_corr ! scrap variables to keep the values of the local bulk velocity
real (DP), parameter :: Stress_TOL=1.0_DP ! The maximum allowed relative perturbation component vise for the corrections of the stress tensor
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! Check for NaNs:
 do i=1,size(f)
  if (f(i).ne.f(i)) then
    print *, "EvalColRelES_SC_1Donecell: NaN in the passed solution at point", i  
  end if
 end do 
!!!!!!!!!!!!!! End checking for NaNs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab,cellnum) ! the subroutine returns momrates
                                           	! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!! Begin Diagnostics
 if (nu<1.0e-6) then 
  print *, "EvalColRelES_SC_1Donecell: Attention, zero or small values of the collision freq. nu!"
  stop
 end if  
 !!!!!!!!! End diagnostics
 dens=sum(f*nodes_gwts)
 mvmaxnorm=0
 ! here we compute the components of the stress tensor and the correction stress tensor
 do i=1,6
   dfij_arr(i) = sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   fij_arr(i) = sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
   !! Compute the max norm of stress tensor, treating matrix as a vector.
   if (abs(fij_arr(i)) > mvmaxnorm) then 
    mvmaxnorm = abs(fij_arr(i))
   end if
 end do
 !!!!!!!!!! Begin diagnostics
 if ((fij_arr(1)<=0) .or. (fij_arr(2)<=0) .or. (fij_arr(3)<=0)) then 
  print *, "EvalColRelES_SC_1Donecell: At leas one temperature component is negative"
  stop
 end if 
 !!!!!!!!!! End diagnostics
 !!!!!!!!!! Begin Diagnostics
 if (mvmaxnorm<1.0E-10) then 
  print *, "EvalColRelES_SC_1Donecell: Attention, zero or small values in the stress tensor!"
  stop
 end if 
 !!!!!!!!! End diagnostics 
 
 !!! here we compute the correction to enforce conservation of temperature 
 !!! this correction will be subtracted on the diagonal of the RelES stress tensor
 en_corr = (momsrates(1)*dfij_arr(1) + momsrates(2)*dfij_arr(2) + momsrates(3)*dfij_arr(3))/3.0_DP
 
 !!!!!!!!! Begin diagnostics 
 if (abs(en_corr)> mvmaxnorm*0.1) then 
  print *, "EvalColRelES_SC_1Donecell: Attention, loss of temperature conservation more than 1%!"
  stop
 end if 
 !!!!!!!!!!! End diagnostics
  
 !!! here we estimate the component-wise perturbation of the proposed corrections.  
 !!! here we estimate the component-wise perturbation of the proposed corrections.  
 !!! maxfactor will measure the lagest correction relative to the vector matrix norm
 !!! This factor is 1 by default and does not affect the model. However, if the corrections are large
 !!! they will be scaled down to keep the target distribution locally supported. The scaling factor 
 !!! is added back in the last step. The controlled moments are not affected, but the 
 !!! rest of the momets will be forced to relax faster than the ES_BGK collision frequency
 !!!  
 maxfactor = 1.0_DP
 do i=1,3
   dfij = abs( (momsrates(i)*dfij_arr(i) - en_corr) / nu)
   if (dfij > Stress_TOL*mvmaxnorm) then 
    factor = dfij/(Stress_TOL*mvmaxnorm)
    if (factor > maxfactor) then 
      maxfactor = factor
    end if
   end if    
 end do 
 do i=4,6
   dfij= abs( (momsrates(i)/nu)*dfij_arr(i) )
   if (dfij > Stress_TOL*mvmaxnorm) then 
    factor = dfij/(Stress_TOL*mvmaxnorm)
    if (factor > maxfactor) then 
      maxfactor = factor
    end if
   end if    
 end do 
 !!! BEGIN Diagnostics comment for production runs
 if ((maxfactor > 4.0) .and. (maxfactor <=16.0)) then 
  print *, "EvalColRelES_SC_1Donecell: Attention! correction factor exceeded 4" 
 end if 
 if (maxfactor > 16.0) then 
  print *, "EvalColRelES_SC_1Donecell: Attention! correction factor exceeded 16" 
 end if 
 !!! END diagnostics
 do i=1,3
  scrtheta(i) = (fij_arr(i) - (momsrates(i)*dfij_arr(i) - en_corr)/nu/maxfactor)*2.0_DP/dens
 end do 
 do i=4,6
  scrtheta(i) = (fij_arr(i) - (momsrates(i)/nu/maxfactor)*dfij_arr(i))*2.0_DP/dens
 end do 
 ! now the components of artificial stress tensor have been computed. We will transfer the values into the artificial stress tensor
 Theta(1,1) = scrtheta(1)
 Theta(2,2) = scrtheta(2)
 Theta(3,3) = scrtheta(3)
 Theta(2,1) = scrtheta(4)
 Theta(1,2) = scrtheta(4)
 Theta(3,1) = scrtheta(5)
 Theta(1,3) = scrtheta(5)
 Theta(3,2) = scrtheta(6)
 Theta(2,3) = scrtheta(6)
 ! we have got the artificial stress tensor set up 
 Theta_inv = inv(Theta) ! get the inverse of the moments matrix
 !!!!!!!!!!!!!!!!!!
 fij = Theta(1,1)*Theta(2,2)*Theta(3,3)+Theta(2,1)*Theta(3,2)*Theta(1,3)+Theta(3,1)*Theta(1,2)*Theta(2,3) - &
       Theta(1,3)*Theta(2,2)*Theta(3,1)-Theta(2,1)*Theta(1,2)*Theta(3,3)-Theta(1,1)*Theta(3,2)*Theta(2,3) ! this will keep the value of the determinant
 !!!!!!!!!! Begin diagnostics
 if (fij<=0)  then 
  print *, "EvalColRelES_SC_1Donecell: det of modified stress tensor is zero or negative"
  stop
 end if 
 !!!!!!!!!! End diagnostics
       
 ! finally, we get the RHS
 fcol = nu*maxfactor*(ESBGK_f0 (Theta_inv,fij,dens,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w) - f )
 !!!!!!!!!! Begin diagnostics
 if (maxval(abs(fcol))>1.0E+2) then 
  print *, "EvalColRelES_SC_1Donecell: high values in the collision operator!", maxval(fcol), minval(fcol)
 end if 
 !!!!!!!!!! End diagnostics
 !!!!!!!!!! Begin diagnostics
 if (maxval(abs(fcol))>1.0E+2) then 
  print *, "EvalColRelES_SC_1Donecell: high values in the collision operator!", maxval(abs(fcol))
 end if 
 !!!!!!!!!! End diagnostics
 !!!!!!!!!! Check for NaNs:
 do i=1,size(fcol)
  if (fcol(i).ne.fcol(i)) then
    print *, "EvalColRelES_SC_1Donecell: NaN in the collision operator at point", i  
  end if
 end do 
!!!!!!!!!!!!!! End checking for NaNs
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
!$OMP CRITICAL
 PRINT *, "EvalColRelES1Donecell: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab, cellnum", mom_L1,Mom_trshld,&
                                 momsrates_reliab, cellnum
!$OMP END CRITICAL
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
!$OMP CRITICAL
!if (temp_dp > 0.0d0) then 
!  entr = (-1)*sum(Real(log(f),DP)*fcol)
!  if (entr< 0.0) then 
!   PRINT *, "EvalColRelES: entropy is negative!:", entr 
!  else  
!   PRINT *, "EvalColRelES: entropy:", entr 
!  end if 
!else 
! PRINT *, "EvalColRelES: velocity distribution is negative at least at one point"
! entr = 0.0d0 
! do i=1,size(f,1)
!  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
! end do
! if (entr< 0.0) then 
!   PRINT *, "EvalColRelES: entropy is negative!:", entr 
!  else  
!   PRINT *, "EvalColRelES: entropy:", entr 
!  end if 
!end if  
!$OMP END CRITICAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColRelES_SC_1Donecell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColGrad(f,Df,fcol,time)
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T
!
! 
! Written by J. Limbacher 2017
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColGrad(f,Df,fcol,time)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w, &
                   !!!!! Diagnostic ADDED THIS TO save the coefficinets on file. Later remove! 
                   Cco 
                   !!!!! END DIAGNOSTIC 
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: sigma, sigma2 ! scrap variable to keep the traceless artificial stress tensor and its inverse 
real (DP), dimension(1:3) :: q !heat flux array
real (DP), dimension(Order) :: moments !Store the moments
real (DP) :: fij, dfij,dens,trace_sigma ! scrap variables
! 
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar ! scrap variables to keep the values of the local bulk velocity
real (DP), dimension (size(nodes_u)) :: vel_trace
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab) ! the subroutine returns moment relaxation rates
											! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 dens=sum(f*nodes_gwts)
 do i=1,10
  fij=sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
  dfij=sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
  moments(i)= (fij - (momsrates(i)/(nu))*dfij)
 end do
 !moments(1) is pressure
 sigma(1,1)=moments(2)
 sigma(2,2)=moments(3)
 sigma(3,3)=moments(4)
 sigma(1,2)=moments(5)
 sigma(2,1)=moments(5)
 sigma(1,3)=moments(6)
 sigma(3,1)=moments(6)
 sigma(2,3)=moments(7)
 sigma(3,2)=moments(7)
 q(1)=moments(8)
 q(2)=moments(9)
 q(3)=moments(10)
 
 fcol = nu * (EvalGrad13f0(dens,moments(1),sigma,q,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w)  - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 PRINT *, "EvalColRelES: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab", momsrates_reliab 
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DIAGNOSTIC: pleae comment if not wanted - next few lines tell if the solution is still physical  
!! EVALUATION OF ENTROPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_dp=minval(f)
if (temp_dp > 0.0d0) then 
  entr = (-1)*sum(Real(log(f),DP)*fcol)
  if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
else 
 PRINT *, "EvalColRelES: velocity distribution is negative at least at one point"
 entr = 0.0d0 
 do i=1,size(f,1)
  entr = entr - Real(log(max(f(i),0.0000001d0)),DP)*fcol(i)
 end do
 if (entr< 0.0) then 
   PRINT *, "EvalColRelES: entropy is negative!:", entr 
  else  
   PRINT *, "EvalColRelES: entropy:", entr 
  end if 
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! END DIAGNOSTICS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine EvalColGrad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColGrad_1Donecell(f,Df,fcol,time,cellnum)
!
! This is a copy of the above subroutine adjusted for 1D spatial solvers
! 
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T 
!
! 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColGrad_1Donecell(f,Df,fcol,time,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                   nu, nodes_gwts,nodes_u,nodes_v,nodes_w
                   
use DGV_distributions_mod

intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: sigma, sigma2 ! scrap variable to keep the traceless artificial stress tensor and its inverse 
real (DP), dimension(1:3) :: q !heat flux array
real (DP), dimension(Order) :: moments !Store the moments
real (DP) :: fij, dfij,dens,trace_sigma ! scrap variables
! 
real (DP) :: L1_err  ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
real (DP) :: ubar,vbar,wbar ! scrap variables to keep the values of the local bulk velocity
real (DP), dimension (size(nodes_u)) :: vel_trace
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,ubar,vbar,wbar,nu,momsrates_reliab,cellnum) ! the subroutine returns momrates
                                           	! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu


if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! In the RelES model, the relaxation rates are enforced by providing coefficients of the (artificial) stress tensor that shows up in the target 
 ! distribution function in the form of inhomogeneous guassian n/(pi^3 \det(T))^1/2 exp(-(u-\bar{u})^{T} T^{-1}(u-\bar{u}))
 !
 ! The coefficients of T are detrmined from this equation: 
 !
 ! ADD EQUATION
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 dens=sum(f*nodes_gwts)
 do i=1,10
  fij=sum(f*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
  dfij=sum(Df*nodes_gwts*kernls_enfrsd_moms(i,nodes_u,nodes_v,nodes_w,ubar,vbar,wbar))
  moments(i)= (fij - (momsrates(i)/(nu))*dfij)
 end do
 !moments(1) is pressure
 sigma(1,1)=moments(2)
 sigma(2,2)=moments(3)
 sigma(3,3)=moments(4)
 sigma(1,2)=moments(5)
 sigma(2,1)=moments(5)
 sigma(1,3)=moments(6)
 sigma(3,1)=moments(6)
 sigma(2,3)=moments(7)
 sigma(3,2)=moments(7)
 q(1)=moments(8)
 q(2)=moments(9)
 q(3)=moments(10)
 
 fcol = nu * (EvalGrad13f0(dens,moments(1),sigma,q,ubar,vbar,wbar,nodes_u,nodes_v,nodes_w)  - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 !$OMP CRITICAL
 PRINT *, "EvalColGrad_1Donecell: Invoke fall back model. mom_L1, Mom_trshld, momsrates_reliab, cellnum", &
                                 momsrates_reliab, cellnum
 !$OMP END CRITICAL
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK(f,fcol)
end if
! all done

end subroutine EvalColGrad_1Donecell


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab)
!
! This subroutine returns the relaxation rates to be used in the 
! model with velocity-dependent collision frequency
! 
! The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!
! f - is the solution in a particular spatial cell (of f is the solution to the spatially homogeneous problem
! momsrates  -- are the values of the relaxation rates to be inforced in the model with velocity dependent collisio nfrequency
! L1_err -- also returns the L1 nort of Df 
!
! 03222017 --- evaluation of nu moved from inside of if -- now evaluation is done at every call
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetRelaxRates0D_DGV(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab)

use DGV_commvar, only: MaxNumEnforcedMoments, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       MoRatesArry0D, MoRatesReliable0D,&
                       alpha, gas_alpha, gas_T_reference, C_inf, gasR

use DGV_distributions_mod
use DGV_dgvtools_mod

real (DP), dimension(:), intent(in) :: f ! the solution(velocity distribution) on one spatial cell
real (DP), dimension(:), intent(in) :: Df ! the difference between the solution and the local maxwellian on the primary mesh on one spatial cell (given, computed elsewhere)
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent (in) :: time ! value of the dimensionless time variable
real (DP), intent (out) :: L1_err ! value of the relative L_1 norm of the differnce between the solution and the local Maxwellian
real (DP), intent (out) :: LocUbar,LocVbar,LocWbar ! the values of the local bulk velocity
real (DP), intent (out) :: nu ! value of the default relaxation rate -- will be assigned to moments for which the ralaxation rates can not be 
                              ! calculated form the Boltzmann collision operator
logical, intent (out) :: momsrates_reliab ! is true if at least one rate has calculated from the Boltzmann collision operator.
!!!!!!!!!
!! Atention: these parameters define how relaxatoin rates are determines. 
real (DP), parameter  :: L1_MAX = 0.5 ! This coefficients determines which form of the Boltzmann collision integral is used. See description below
real (DP), parameter  :: L1_SENS_TRSHLD = 0.2 ! This parameters determines the level of sensitivety of expression for evaluation of the relaxation rates. 
          ! the energy moment of the Botlzmann collision integral has to be zero. If it is not zero, this is only due to truncation errors =e_{2}. Thus this moment 
          ! can be used to judge about the truncation moments in moments of the collision operator. (The errors expeced to be larger for higher moments). 
          ! when evaluating the relaxation rates, the moments $m=\int_{R^3} Q(f,f)\phi$ evaluated with error $e_{m}$. We will hope that $e_{m}$ is comparable to $e_2$. 
          ! We will use this parameter to make sure that 
          ! m> L1_SENS_TRSHLD * e_{2}. Similarly we will compare $f_{\phi}-f^{M}_{\phi}$  -- the 
          ! difference between the moment and the same moment evaluated on a local Maxwellan -- to $L1_SENS_TRSHLD em_{2}$. If this difference is at least (L1_SENS_TRSHLD \cdot e_{2}) 
          ! , that is it is at least L1_SENS_TRSHLD times larger than the expected errors then 
          ! rate = (m+e_{m})/(f_{\phi}-f^{M}_{\phi})= true rate + (e_m/e_2)1/L1_SENS_TRSHLD
          ! 
          ! thus, in a sense, L1_SENS_TRSHLD gives the number of correct digits in the computed rate IF e_m\approx e_2 
          ! 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), parameter :: Mom_Close_Maxwl_TRSHLD = 1.0D-4 ! This is the second parameter that determines the sensitivity treshhold for evaluation 
                                                         ! of the relaxation rates for moments. It is possible that the moment is already close to 
                                                         ! its final state. In this case, assuming that there is more noice in the derivative of the moment 
                                                         ! than in the difference of the moment (especially if we use course mesh for evaluating the moment 
                                                         ! derivative, the best strategy is not to compute the relaxation rate for this moment. 
                                                         ! In particular, this will avoid computing rates for conserved moments, when sufficent resolution is used.  

!! End Attention
!!!!!!!!!
real (DP) :: LocDens ! value of the local numerical density is returned 
real (DP) :: LocTempr ! value of the local temperature is returned
real (DP), dimension(:), allocatable :: fII,fcolII,fcol1II,fMII ! temp variables to store the 
 ! value of the collision integral on secondary and the local Maxwellian on the secondary mesh and the 
 ! differnce between the local maxwellian and the local maxwellian on the secondary mesh.
 ! and the difference between the local maxwellian and the solution on the promary mesh
logical :: updateNulcl ! scrap variable to pass the update flag
integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check that the nucoeffs has the right size:
if ((size(momsrates,1) > MaxNumEnforcedMoments) .or. (size(momsrates,1) /= size(MoRatesArry0D,1))) then
 print *, "GetRelaxRates_DGV: Error. Size of the supplied array (momsrates) is incopatible."
 stop ! terminate the program sinse nucoefs has a wrong size
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we check if relaxation speeds in this cell need to be updated.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call CheckNuUpdateNeededF0Da_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err) ! the subroutine accesses all other variables directly from commvar
          ! the subrouine returns a value true if the relaxation speeds for the velocity dependent 
          ! collision frequency need to be updated. 
          ! Also this version of the subroutine returns some useful byproducs of the 
          ! check: Df = value of the difference between f and the local Mawellian, 
          ! relative L_1 norm of the differnce and the Macroparatmers of the local Maxwellian
          !    
          ! This versio of the subroutine works for the solution to 
          ! spatially homogeneous problem (0D). For 1D-2D and 3D problems, use versions of the subroutine that 
          ! work with an arrays in which index corresponding to different spatial cell 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
! Evaluate ES-BGK nu - collision frequency
nu = LocDens*LocTempr/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/LocTempr)**gas_alpha  ! Dimensionless nu -- the classical collision frequency of the ES-BGK model
! classican ES-BGK nu is used as a backup collision frequency 
if (updateNulcl) then
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! to update the relaxation frequencies for the selected group of moments, the 
 ! full Boltzmann collision operator is evaluated using secondary velocity mesh. 
 ! for that the solution needs to be interpolated to the secondary mesh. 
 ! then evaluation of the Boltzmann collision integral is colled. 
 ! once the collision integral is evaluated, the new relxation rates are computed. 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
 allocate (fII(1:size(nodes_uII,1)), fcolII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates0D_DGV: Allocation error for variables (fII),(fcolII)"
      end if 
 ! To project the solution to a (most likely smaller) secondary mesh, we use the following subroutine: 
 call ProjectDGVtoSecMesh(f,fII) ! put the projecteion in fII
 ! Now we will check if the solution is needs to be evaluated in a decomposition mode or a full mode. 
 ! in the full mode $\int_{R^3}\int_{R^3} ffA is evaluated
 ! in the decomposition mode, the solution is split into a sum of the local Maxwellian and another functions
 ! $f=f_{M} + f_{0}$ (note that $f_{0}$ does not have to be small). and the integral 
 ! $\int_{R^3}\int_{R^3} (2f_{M}f_{0}A + f_{0}f_{0}A ) is evaluated instead
 ! this approach reduces errors when the solution is cloase to a Mawellian
 !!! 
 ! because we expect that the solution on the primary mesh is mote accurate, we will check the 
 ! closeness of the solution to the Maxwellian using the primary mesh:
 ! We will use the byproduct of another other check, L1_err -- which is the relative L_1 norm of 
 ! the difference between the distribution f and the local maswellian on the primary mesh
 if (L1_err > L1_MAX) then 
  ! evaluate the Boltzmann collision operator using full mode
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now we will call the subroutine that evaluates the Boltzmann collision operator on the secondary mesh. 
  ! There are two possibilites: an single processor call of the subroutine with possible OpenMP fork and 
  ! call or an MPI parallelization.  
  !!!!!!
  ! to make an single processor or an OpenMP evaluation, uncomment the next lime
  call EvalCollisionPeriodicAPlus_DGVII_OMP(fII,fcolII) ! fcolII contains the value of the collision operator
  ! to make an MPI evaluation, uncomment the next line, but this will involve a lot of preparatory work:
  ! includind setting up MPI universes for the collision integral and such .
  ! this preparatory work should be done in the InitDGV0D/1D subroutine. Make sure all necessary lines are uncommented. 
  ! ADD MPI!!!! ! call EvalCollisionPeriodicAPlusMPI_DGVII(DfII,fcolII)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
 else 
  ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
  allocate (fcol1II(1:size(nodes_uII,1)), fMII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates0D_DGV: Allocation error for variables (fMII),(fcol1II)"
      end if 
  fMII = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_uII,nodes_vII,nodes_wII) ! now we populate the maxwellian with the same macroparamters.
  fII = fII-fMII ! now fII contains the difference between the solution and the local maxwellian
  ! Now we will make to calls of the collision operator
  !!!!!!!!!!
  ! to make an single processor or an OpenMP evaluation, uncomment the next lime
  call EvalCollisionPeriodicAPlus_DGVII_OMP(fII,fcolII) ! fcolII contains the value of the collision operator f_{0}f_{0}A
  call EvalCollisionPeriodicMixedTermsA_DGVII_OMP(fII,fmII,fcol1II) !
  ! to make an MPI evaluation, uncomment the next line, but this will involve a lot of preparatory work:
  ! includind setting up MPI universes for the collision integral and such .
  ! this preparatory work should be done in the InitDGV0D/1D subroutine. Make sure all necessary lines are uncommented. 
  ! ADD MPI!! STILL NOT WORKING!!!! 
  ! call EvalCollisionPeriodicAPlusMPI_DGVII(fII,fcolII)
  ! call EvalCollisionPeriodicMixedTermsAMPI_DGVII(fII,fmII,fcol1II) !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
  fcolII=fcolII+fcol1II
  deallocate (fcol1II,fMII)
 end if     
 deallocate(fII)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! Next we compute the relaxation rates for the selected group of moments from the (Df), and (fcolII).
 ! (Df) on the primary mesh is taken to provide a better computation of the macroparameters.  
 ! fcolII is only available on the secondary mesh  
 !!!!!!!!!!!!!!!!!
 call ComputeRelaxRatesBCIa_DGV(fcolII,Df,momsrates,L1_SENS_TRSHLD,Mom_Close_Maxwl_TRSHLD,&
                     LocDens,LocUbar,LocVbar,LocWbar,LocTempr,nu,momsrates_reliab) ! The first subroutine computes quantities directly
 MoRatesArry0D = momsrates ! save the coefficients for future use.
 MoRatesReliable0D = momsrates_reliab  ! save the reliability flag for the moment relaxation rates. true means that at least one rate is computed based on the Bolzmann collision integral 
 !!!!!!!!!!!!!!!!!
 deallocate(fcolII)
else 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! if the no update is needed, the rates are taken from the storage array
 momsrates = MoRatesArry0D
 momsrates_reliab = MoRatesReliable0D
 !!!! 
end if
end subroutine GetRelaxRates0D_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab,cellnum)
!
! This subroutine returns the relaxation rates to be used in the 
! model with velocity-dependent collision frequency
! 
! The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!
! f - is the solution in a particular spatial cell (of f is the solution to the spatially homogeneous problem
! momsrates  -- are the values of the relaxation rates to be inforced in the model with velocity dependent collisio nfrequency
! L1_err -- also returns the L1 nort of Df 
!
!
! 03222017 --- evaluation of nu moved from inside of if -- now evaluation is done at every call
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetRelaxRates1Donecell_DGV(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab,cellnum)

use DGV_commvar, only: MaxNumEnforcedMoments, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       MoRatesArry1D, MoRatesReliable1D,&
                       alpha, gas_alpha, gas_T_reference, C_inf, gasR

use DGV_distributions_mod
use DGV_dgvtools_mod

real (DP), dimension(:), intent(in) :: f ! the solution(velocity distribution) on one spatial cell
real (DP), dimension(:), intent(in) :: Df ! the difference between the solution and the ;local maxwellian on the primary mesh on one spatial cell (given, computed elsewhere)
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent(in) :: time ! value of the dimensionless time variable
real (DP), intent(out) :: L1_err ! value of the relative L_1 norm of the differnce between the solution and the local Maxwellian
real (DP), intent (out) :: LocUbar,LocVbar,LocWbar ! the values of the local bulk velocity
real (DP), intent (out) :: nu ! value of the default relaxation rate -- will be assigned to moments for which the ralaxation rates can not be 
                              ! calculated form the Boltzmann collision operator (for this cell)
logical, intent (out) :: momsrates_reliab ! is true if at least one rate has calculated from the Boltzmann collision operator. (for this cell)
integer (I4B), intent (in) :: cellnum ! the number of the spatial cell for which the rates are evaluated
!!!!!!!!!
!! Atention: these parameters defines how relaxatoin rates are determines. 
real (DP), parameter  :: L1_MAX = 0.8 ! This coefficients determines which form of the Boltzmann collision integral is used. See description below
real (DP), parameter  :: L1_SENS_TRSHLD = 1.0 ! This parameters determines the level of sensitivety of expression for evaluation of the relaxation rates. 
          ! the energy moment of the Botlzmann collision integral has to be zero. If it is not zero, this is only due to truncation errors =e_{2}. Thus this moment 
          ! can be used to judge about the truncation moments in moments of the collision operator. (The errors expeced to be larger for higher moments). 
          ! when evaluating the relaxation rates, the moments $m=\int_{R^3} Q(f,f)\phi$ evaluated with error $e_{m}$. We will hope that $e_{m}$ is comparable to $e_2$. 
          ! We will use this parameter to make sure that 
          ! m> L1_SENS_TRSHLD * e_{2}. Similarly we will compare $f_{\phi}-f^{M}_{\phi}$  -- the 
          ! difference between the moment and the same moment evaluated on a local Maxwellan -- to $L1_SENS_TRSHLD em_{2}$. If this difference is at least (L1_SENS_TRSHLD \cdot e_{2}) 
          ! , that is it is at least L1_SENS_TRSHLD times larger than the expected errors then 
          ! rate = (m+e_{m})/(f_{\phi}-f^{M}_{\phi})= true rate + (e_m/e_2)1/L1_SENS_TRSHLD
          ! 
          ! thus, in a sense, L1_SENS_TRSHLD gives the number of correct digits in the computed rate IF e_m\approx e_2 
          ! 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), parameter :: Mom_Close_Maxwl_TRSHLD = 1.0D-5 ! This is the second parameter that determines the sensitivity treshhold for evaluation 
                                                         ! of the relaxation rates for moments. It is possible that the moment is already close to 
                                                         ! its final state. In this case, assuming that there is more noice in the derivative of the moment 
                                                         ! than in the difference of the moment (especially if we use course mesh for evaluating the moment 
                                                         ! derivative, the best strategy is not to compute the relaxation rate for this moment. 
                                                         ! In particular, this will avoid computing rates for conserved moments, when sufficent resolution is used.  


!! End Attention
!!!!!!!!!
real (DP) :: LocDens ! value of the local numerical density is returned 
real (DP) :: LocTempr ! value of the local temperature is returned
real (DP), dimension(:), allocatable :: fII,fcolII,fcol1II,fMII ! temp variables to store the 
 ! value of the collision integral on secondary and the local Maxwellian on the secondary mesh and the 
 ! differnce between the local maxwellian and the local maxwellian on the secondary mesh.
 ! and the difference between the local maxwellian and the solution on the promary mesh
logical :: updateNulcl ! scrap variable to pass the update flag
integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check that the nucoeffs has the right size:
if ((size(momsrates,1) > MaxNumEnforcedMoments) .or. (size(momsrates,1) /= size(MoRatesArry1D,1))) then
 print *, "GetRelaxRates_DGV: Error. Size of the supplied array (momsrates) is incopatible."
 stop ! terminate the program sinse nucoefs has a wrong size
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we check if relaxation speeds in this cell need to be updated.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call CheckNuUpdateNeededF1Donecella_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err,cellnum) ! the subroutine accesses all other variables directly from commvar
          ! the subrouine returns a value true if the relaxation speeds for the velocity dependent 
          ! collision frequency need to be updated. 
          ! Also this version of the subroutine returns some useful byproducs of the 
          ! check: Df = value of the difference between f and the local Mawellian, 
          ! relative L_1 norm of the differnce and the Macroparatmers of the local Maxwellian
          !    
          ! This versio of the subroutine works for the solution to 
          ! spatially homogeneous problem (0D). For 1D-2D and 3D problems, use versions of the subroutine that 
          ! work with an arrays in which index corresponding to different spatial cell 
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate ES-BGK nu - collision frequency
nu = LocDens*LocTempr/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/LocTempr)**gas_alpha  ! Dimensionless nu -- the classical collision frequency of the ES-BGK model
! classican ES-BGK nu is used as a backup collision frequency 
if (updateNulcl) then
!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! to update the relaxation frequencies for the selected group of moments, the 
 ! full Boltzmann collision operator is evaluated using secondary velocity mesh. 
 ! for that the solution needs to be interpolated to the secondary mesh. 
 ! then evaluation of the Boltzmann collision integral is colled. 
 ! once the collision integral is evaluated, the new relxation rates are computed. 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
 allocate (fII(1:size(nodes_uII,1)), fcolII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates1Donecell_DGV: Allocation error for variables (fII),(fcolII)"
      end if 
 ! To project the solution to a (most likely smaller) secondary mesh, we use the following subroutine: 
 call ProjectDGVtoSecMesh(f,fII) ! put the projecteion in fII
 ! Now we will check if the solution is needs to be evaluated in a decomposition mode or a full mode. 
 ! in the full mode $\int_{R^3}\int_{R^3} ffA is evaluated
 ! in the decomposition mode, the solution is split into a sum of the local Maxwellian and another functions
 ! $f=f_{M} + f_{0}$ (note that $f_{0}$ does not have to be small). and the integral 
 ! $\int_{R^3}\int_{R^3} (2f_{M}f_{0}A + f_{0}f_{0}A ) is evaluated instead
 ! this approach reduces errors when the solution is cloase to a Mawellian
 !!! 
 ! because we expect that the solution on the primary mesh is mote accurate, we will check the 
 ! closeness of the solution to the Maxwellian using the primary mesh:
 ! We will use the byproduct of another other check, L1_err -- which is the relative L_1 norm of 
 ! the difference between the distribution f and the local maswellian on the primary mesh
 if (L1_err > L1_MAX) then 
  ! evaluate the Boltzmann collision operator using full mode
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now we will call the subroutine that evaluates the Boltzmann collision operator on the secondary mesh. 
  ! There are two possibilites: an single processor call of the subroutine with possible OpenMP fork and 
  ! call or an MPI parallelization.  
  !!!!!!
  ! to make an single processor or an OpenMP evaluation, uncomment the next lime
  call EvalCollisionPeriodicAPlus_DGVII(fII,fcolII) ! fcolII contains the value of the collision operator
  ! to make an MPI evaluation, uncomment the next line, but this will involve a lot of preparatory work:
  ! includind setting up MPI universes for the collision integral and such .
  ! this preparatory work should be done in the InitDGV0D/1D subroutine. Make sure all necessary lines are uncommented. 
  ! ADD MPI!!!! ! call EvalCollisionPeriodicAPlusMPI_DGVII(DfII,fcolII)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
 else 
  ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
  allocate (fcol1II(1:size(nodes_uII,1)), fMII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates1Donecell_DGV: Allocation error for variables (fMII),(fcol1II)"
      end if 
  fMII = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_uII,nodes_vII,nodes_wII) ! now we populate the maxwellian with the same macroparamters.
  fII = fII-fMII ! now fII contains the difference between the solution and the local maxwellian
  ! Now we will make to calls of the collision operator
  !!!!!!!!!!
  ! to make an single processor or an OpenMP evaluation, uncomment the next lime
  call EvalCollisionPeriodicAPlus_DGVII(fII,fcolII) ! fcolII contains the value of the collision operator f_{0}f_{0}A
  call EvalCollisionPeriodicMixedTermsA_DGVII(fII,fMII,fcol1II) !
  ! to make an MPI evaluation, uncomment the next line, but this will involve a lot of preparatory work:
  ! includind setting up MPI universes for the collision integral and such .
  ! this preparatory work should be done in the InitDGV0D/1D subroutine. Make sure all necessary lines are uncommented. 
  ! ADD MPI!! STILL NOT WORKING!!!! 
  ! call EvalCollisionPeriodicAPlusMPI_DGVII(fII,fcolII)
  ! call EvalCollisionPeriodicMixedTermsAMPI_DGVII(fII,fmII,fcol1II) !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
  fcolII=fcolII+fcol1II
  deallocate (fcol1II,fMII)
 end if     
 deallocate(fII)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! Next we compute the relaxation rates for the selected group of moments from the (Df), and (fcolII).
 ! (Df) on the primary mesh is taken to provide a better computation of the macroparameters.  
 ! fcolII is only available on the secondary mesh  
 !!!!!!!!!!!!!!!!!
 ! call ComputeRelaxRatesBCI_DGV(fcolII,f,momsrates,L1_SENS_TRSHLD)    ! Both subroutines do the same, but the second one uses qunatities compluted previously 
 call ComputeRelaxRatesBCIa_DGV(fcolII,Df,momsrates,L1_SENS_TRSHLD,Mom_Close_Maxwl_TRSHLD,LocDens,LocUbar,LocVbar, &
        LocWbar,LocTempr,nu,momsrates_reliab) ! The first subroutine computes quantities directly
 MoRatesArry1D(:,cellnum) = momsrates ! safe the coefficients for future use.
 MoRatesReliable1D(cellnum) = momsrates_reliab ! update the flag is the relaxation rates were computed form the BCI 
 !!!!!!!!!!!!!!!!!
 deallocate(fcolII)
else 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! if the no update is needed, the rates are taken from the storage array
 momsrates = MoRatesArry1D(:,cellnum)
 momsrates_reliab = MoRatesReliable1D(cellnum) ! the flag = true means that rates were computed from the BCI
 !!!! 
end if
end subroutine GetRelaxRates1Donecell_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UniversalCollisionOperator1DonecellDGV(f,fcol,time,dt,cellnum)
! 
! This is a copy of the above subrouine to use in multidimensional applications,
! 
! It is assumed that all spatial cells are numbered with one index. This index is passed 
! in the variable cellnum. Some data about past evaluation is stored for each cell and the index is used to access the data
! Currenlty is only used in the BGK-Velocity-Dependent Collision Model
! 
! 
! This subroutine evaluates the contribution 
! due to particle collisions to the
! evolution of the velocity distribution 
! function in  a kinetic equation.
! 
! The current state of the solution is given in the 
! variable (f), the value of the collision operator is 
! returned in the valriable (fcol). The time step (dt) 
! is multiplied to the value of the collision integral 
! in the process of evaluating the collision operator
!
! The evaluation of the collision operator 
! may use any of the following models:
! full Boltzmann equation; 
! decomposed Boltzmann equation
! linearized Boltzmann equation
! BGK-Type Model with velocity dependent collision frequency
! Classical BGK, ES-BGK and Shakhov models. 
!
! The chocie of what model to use depends on the state of the 
! solution (f) and also on the parameters od evlution
! 
! The subroutine accesses variables of DGV_commvar
! 
!!!!!!!!!!!!!!!!!!!!

subroutine UniversalCollisionOperator1DonecellDGV(f,fcol,time,dt,cellnum)

use DGV_commvar, only: run_mode_1D,mol_diam,L_inf,N_inf,T_inf,C_inf,gas_viscosity,gasR
				   
				   
use DGV_dgvtools_mod

real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (:), intent (out) :: fcol ! value of the right side  
real (DP), intent (in) :: time ! the current time
real (DP), intent (in) :: dt ! The time step
integer (I4B), intent (in) :: cellnum ! the number of the spatial cell/point for which this distribution f function corresponds. 
!!!!!!!!!!!
real (DP), dimension (:), allocatable :: fcol_scr,Df ! scatch variable to keep the right side and the perturbation part
real (DP) :: coef_temp ! Scrap variable
real (DP), parameter :: kBoltzmann = 1.3080648813D-23 ! Boltzmann constant J/K
integer (I4B) :: run_mode ! scrap variable to store run_mode in the celected cell
integer :: loc_alloc_stat

! WARNING: Make sure that Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! evaluation of the collision operator.

!!!!!!!!!!!! set us the allocatable array !
allocate (Df(1:size(fcol,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "UniversalCollisionOperator1DonecellDGV: Allocation error for variables (Df)"
 stop
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call CheckSolutionModeSimple1Donecell_DGV(f,Df,cellnum,.false.) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.

run_mode = run_mode_1D(cellnum)

select case (run_mode) ! run_mode is set outside, when solution is evaluated for closeness to a maxwellian...  
	case (0) ! run_mode=0 means we are very far from a Maxwellian. In this case we just call the collision operator
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! UNCOMMENT if evaluating collision integral using Gauss Quadratures
	  !!!
	  !!! Two choises for the call of collision operator: Uncomment only one of them!  These procedures do the same, but slightly different in implementation
	  !!! call EvalCollisionPeriodicA_DGV(f,fcol1) ! This one uses intermediate arrays and is slower
	  !!! call EvalCollisionPeriodicAPlus_DGV(f,fcol) ! This one is a little faster that the one above...
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! Uncomment if evaluating collision integral using Korobov Quadaratures
	  !!!
!	   call EvalCollisionPeriodicAKorOpt_DGV(f,fcol)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Now we are adding the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
      coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  fcol = fcol*coef_temp
	  !! 
	  !!!!
	case (1) ! this is the non-linear perturbation mode: we need both the linear part and the non-linear 
	  !!!!!!!!!!!! set us the allocatable array !
      allocate (fcol_scr(1:size(fcol,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "UniversalCollisionOperator1DDGV: Allocation error for variables (fcol_scr)"
       stop
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! Uncomment both lines if using Guass integration of the collision integral 
	  !!! call EvalCollisionPeriodicMixedTermsA_DGV(Df,f-Df,fcol_scr)        ! This evaluates the linear part
	  !!! call EvalCollisionPeriodicAPlus_DGV(Df,fcol) ! this evaluates the non-linear part 
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! Uncomment both lines if using Korobov quadratures
!	  call EvalCollisionPeriodicAKorMxdTmsOpt_DGV(Df,f-Df,fcol_scr)
!	  call EvalCollisionPeriodicAKorOpt_DGV(Df,fcol)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  fcol = fcol+fcol_scr
	  deallocate (fcol_scr)
	  !! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2 
	  coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  fcol = fcol*coef_temp
	  !!
	case (2) ! this is the linear mode, we pretty much neglect the quadratic part..
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Guass integration of the collision integral 
	  ! call EvalCollisionPeriodicMixedTermsA_DGV(Df,f-Df,fcol)        ! This evaluates the linear part
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Korobov quadratures
!	  call EvalCollisionPeriodicAKorMixedTerms_DGV(Df,f-Df,fcol)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
	  coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  fcol = fcol*coef_temp
	  !! 
	  !!!!
	case (3) ! the velocity dep. ES-BGK mode
		coef_temp = (kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity
	    !! Now we are adding the dimensionless coefficient
	    !call EvalColESBGK(f,fcol) ! UNCOMMENT IF USING THE classical ES MODEL
	    !call EvalColShakov(f,fcol)  ! UNCOMMENT IF USING THE SHAKOV MODEL
		!call EvalColVelES1Donecell(f,Df,fcol,time,cellnum) ! Uncomment if using velocity dependent model
		call EvalColRelES_SC_1Donecell(f,Df,fcol,time,cellnum)
		fcol = fcol*coef_temp*dt
	case (4) ! ES-BGK mode or Shakhov  
		coef_temp = (kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity
		call EvalColESBGK(f,fcol) ! UNCOMMENT IF USING THE ES MODEL
		!call EvalColShakov(f,fcol)  ! UNCOMMENT IF USING THE SHAKOV MODEL
		fcol = fcol*coef_temp*dt
	case default 
		print *, "UniversalCollisionOperator1DonecellDGV: cannot process value of run_mode", run_mode
		stop
end select

deallocate (Df)

end subroutine UniversalCollisionOperator1DonecellDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UniversalCollisionOperator0DDGV(f,fcol,time,dt)
! 
! This is the version of the subroutine to use in spatially homogeneous problem.
!
! This subroutine evaluates the contribution 
! due to particle collisions to the
! evolution of the velocity distribution 
! function in  a kinetic equation.
! 
! The current state of the solution is given in the 
! variable (f), the value of the collision operator is 
! returned in the valriable (fcol). The time step (dt) 
! is multiplied to the value of the collision integral 
! in the process of evaluating the collision operator
!
! The evaluation of the collision operator 
! may use any of the following models:
! full Boltzmann equation; 
! decomposed Boltzmann equation
! linearized Boltzmann equation
! BGK-Type Model with velocity dependent collision frequency
! Classical BGK, ES-BGK and Shakhov models. 
!
! The chocie of what model to use depends on the state of the 
! solution (f) and also on the parameters od evlution
! 
! The subroutine accesses variables of DGV_commvar
! 
!!!!!!!!!!!!!!!!!!!!

subroutine UniversalCollisionOperator0DDGV(f,fcol,time,dt)

use DGV_commvar, only: run_mode,mol_diam,L_inf,N_inf,T_inf,C_inf,gas_viscosity,gasR,&
				   alpha,gas_T_reference,gas_alpha,&
				   !debug
				   nodes_u,nodes_v,nodes_w
				   
use DGV_dgvtools_mod

real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (:), intent (out) :: fcol ! value of the right side  
real (DP), intent (in) :: time ! the current time
real (DP), intent (in) :: dt ! The time step 
!!!!!!!!!!!
real (DP), dimension (:), allocatable :: fcol_scr,Df,fcol1 ! scatch variable to keep the right side and the perturbation part
real (DP) :: coef_temp,xxx,xxx0,xxx1,xxx2, seconds ! Scrap variable
real (DP), parameter :: kBoltzmann = 1.3080648813D-23 ! Boltzmann constant J/K
integer :: loc_alloc_stat
! WARNING: Make sure that Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! evaluation of the collision operator.

interface 
 function omp_get_wtime() result (y)
   double precision :: y 
 end function omp_get_wtime
end interface  


!!!!!!!!!!!! set us the allocatable array !
allocate (Df(1:size(fcol,1)),fcol1(1:size(fcol,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "UniversalCollisionOperator0DDGV: Allocation error for variables (Df)"
 stop
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call CheckSolutionMode_DGV(f,Df) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.
!
!
select case (run_mode) ! run_mode is set outside, when solution is evaluated for closeness to a maxwellian...  
	case (0) ! run_mode=0 means we are very far from a Maxwellian. In this case we just call the collision operator
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! UNCOMMENT if evaluating collision integral using Gauss Quadratures
	  !!!
	  !!! Two choises for the call of collision operator: Uncomment only one of them!  These procedures do the same, but slightly different in implementation
	  !!! call EvalCollisionPeriodicA_DGV(f,fcol1) ! This one uses intermediate arrays and is slower
	  !!! DEBUG: EVALUATE PROCESSOR TIME NEXT LINE
	  !!seconds = omp_get_wtime ( )
      !!call cpu_time (xxx)
	  !!! END DEBUG
	  call EvalCollisionPeriodicAPlus_DGV_OMP(f,fcol) ! This one is a little faster that the one above...
	  !!! DEBUG: EVALUATE PROCESSOR TIME NEXT LINE
	  !! call cpu_time (xxx0)
	  !! seconds = omp_get_wtime ( ) - seconds;
	  !! print *, "Wall time lapsed in seconds:", seconds
	  !! print *, "Total CPU time for all procs:", xxx0 - xxx
	  !! stop
	  !!! END DEBUG
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !!! Uncomment if evaluating collision integral using Korobov Quadaratures
	  !!!
	  ! call EvalCollisionPeriodicAKorOpt_DGV(f,fcol)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Now we are adding the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
      coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  fcol = fcol*coef_temp
	  !! 
	  !!!!
	case (1) ! this is the non-linear perturbation mode: we need both the linear part and the non-linear 
	  !!!!!!!!!!!! set us the allocatable array !
      allocate (fcol_scr(1:size(fcol,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "UniversalCollisionOperator0DDGV: Allocation error for variables (fcol_scr)"
       stop
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Guass integration of the collision integral 
	  call EvalCollisionPeriodicMixedTermsA_DGV_OMP(Df,f-Df,fcol_scr)        ! This evaluates the linear part
	  call EvalCollisionPeriodicAPlus_DGV_OMP(Df,fcol) ! this evaluates the non-linear part 
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Korobov quadratures
	  !call EvalCollisionPeriodicAKorMxdTmsOpt_DGV(Df,f-Df,fcol_scr)
	  !call EvalCollisionPeriodicAKorOpt_DGV(Df,fcol)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  fcol = fcol+fcol_scr
	  deallocate (fcol_scr)
	  !! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2 
	  coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  !xxx=sum(fcol)
	  fcol = fcol*coef_temp
	  !!
	case (2) ! this is the linear mode, we pretty much neglect the quadratic part..
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Guass integration of the collision integral 
	  call EvalCollisionPeriodicMixedTermsA_DGV_OMP(Df,f-Df,fcol)        ! This evaluates the linear part
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  ! Uncomment if using Korobov quadratures
	  !call EvalCollisionPeriodicAKorMixedTerms_DGV(Df,f-Df,fcol)
	  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  !! Now we are addiing the dimensionless coefficient mol_diam^2*N_inf/(L_inf)^2
	  coef_temp = (mol_diam/L_inf)**2*N_inf*dt
	  fcol = fcol*coef_temp
	case (3) ! the velocity dep. ES-BGK mode
		coef_temp = (kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity
	    !! Now we are adding the dimensionless coefficient
	    !call EvalColESBGK(f,fcol) ! UNCOMMENT IF USING THE classical ES MODEL
	    !call EvalColShakov(f,fcol)  ! UNCOMMENT IF USING THE SHAKOV MODEL
		!call EvalColVelES(f,Df,fcol,time) ! Uncomment if using velocity dependent model
        call EvalColRelES_SC(f,Df,fcol,time) ! Uncomment if using model whith ES-BGK collision frequency and artificial gaussian to enforce second moments
		fcol = fcol*coef_temp*dt
	case (4) ! ES-BGK mode or Shakhov  
		coef_temp = (kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity
		call EvalColESBGK(f,fcol) ! UNCOMMENT IF USING THE ES MODEL
		!call EvalColShakov(f,fcol)  ! UNCOMMENT IF USING THE SHAKOV MODEL
		fcol = fcol*coef_temp*dt
	case default 
		print *, "SpatialOperSH_DGV: cannot process value of run_mode", run_mode
		stop
end select

deallocate (Df)

end subroutine UniversalCollisionOperator0DDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionMode_DGV(f,Df) 
! 
! This subroutine evaluates macroparameters of the solution. Then it subtracts 
! from the solution a Maxwellian with the computed macroparameters. This difference is placed into Df
! then a norm of Df is assessed. If |Df| is large, (rum_mode=0) then we are in a strongly non-equilibrium 
! regime and nothing is done. 
! if the |Df| is moderate, then we switch to mode 1 (run_mode=1) then linear and quadratic contributions are computed separately
! if |Df| is small, we switch to run_mode=2 and only the linear part of collision opertator is computed.   
! 
! This subroutine also checks the solution and calls for the update of the linearized operator fmA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckSolutionMode_DGV(f,Df)

use DGV_commvar, only: run_mode, LocDens, LocUbar, LocVbar, LocWbar, LocTempr, fm, nodes_u, nodes_v, nodes_w, &
                   decomp_lev, linear_lev, vel_lev, ES_lev, nodes_gwts

use DGV_distributions_mod
use DGV_dgvtools_mod

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out):: Df ! the difference between the solution and the local Maxwellian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: n, u_0,v_0,w_0
real (DP), dimension (size(f,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err, L1_err_fm ! a sctratch variables for the errors

 
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
Df = f-fMaxwellNew ! added by Craig
L1_err = SUM(ABS(Df)*nodes_gwts)/LocDens! evaluate the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Print May comment
print*, "L1 error = ", L1_err
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   select case (run_mode)
   case (0) ! if we're in a strongly non-linear regime on the previous time step, then do this: 
      fm=fMaxwellNew ! now fm holds the last computed local maxwellian -- unless we will go back into mode =0 or mode = 1 with some extra condition, 
                     ! then fm will stay the same -- this is something to remember when using the in multiple dimensions this will  
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		print *, "Switching to ES-BGK mode from strongly non-linear"
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		print *, "Switching to vel ES-BGK mode from strongly non-linear regime"
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
        run_mode=2                          ! linear regime
        print *, "Switching to linearized regime from strongly non-linear"
      elseif (L1_err < decomp_lev) then
	    !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_{m}(v)A(v,v_1,phi) dv
        run_mode=1                      ! decomposition regime
        print *, "Switching to decomposition regime from strongly non-linear"
      end if
   case (1) ! if already we're in the decomposition regime
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		print *, "Switching to ES-BGK mode from decomposition regime"
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		print *, "Switching to vel ES-BGK mode from decomposition regime"
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        L1_err_fm=SUM(ABS(fm-fMaxwellNew))/LocDens   ! next we check if the linearized operator needs to be updated
        if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
        end if
        run_mode=2                          ! linear regime
        print *, "Switching to linearized regime from decomposition regime"
      else if (L1_err < decomp_lev) then
        L1_err_fm=SUM(ABS(fm-fMaxwellNew))/LocDens   ! next we check if the linearized operator needs to be updated
        if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
		  ! updated the linear operator
		  print*, "updated the linear operator in mode=1"
        end if
      else
	    run_mode=0
		print*,"switching to strongly non-linear from decomposition regime"
	  end if

   case (2) ! if we're in the linear regime..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		print *, "Switching to ES-BGK mode from linear regime"
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		print *, "Switching to vel ES-BGK mode from linear regime"
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        L1_err_fm=SUM(ABS(fm-fMaxwellNew))/LocDens   ! next we check if the linearized operator needs to be updated
        if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
		  print*, "updated the linear operator in mode=2"
        end if

      else if (L1_err < decomp_lev) then
        L1_err_fm=SUM(ABS(fm-fMaxwellNew))/LocDens   ! next we check if the linearized operator needs to be updated
        if (L1_err_fm > linear_lev/2.0_DP ) then ! we use the linear_level as a measure of two perturbations being close
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
		  ! updated the linear operator
        end if
        run_mode=1                          ! linear regime
        print *, "Switching to decomposition regime from linear regime"
      else
	    run_mode=0
		print*,"switching to strongly non-linear from linear regime"
	  end if

	case (3) ! vel dep ES-BGK mode
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		print *, "Switching to ES-BGK mode from vel-dependent regime"
	  elseif (L1_err < vel_lev) then
	    run_mode=3
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
		  print*, "switching to linear reg from vel-dep"
		  run_mode=2
      else if (L1_err < decomp_lev) then
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
        !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
        ! updated the linear operator
        run_mode=1                          ! linear regime
        print *, "Switching to decomposition regime from vel-dep regime"
      else
	    run_mode=0
		print*,"switching to strongly non-linear from vel-dep regime"
	  end if

	case (4) ! if we're in the ES-BGK mode..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		print *, "Switching to vel-dep mode from ES regime"
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
          fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
          !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
		  print*, "switching to linear reg from ES regime"
		  run_mode=2
      else if (L1_err < decomp_lev) then
        fm=fMaxwellNew ! If the new maxwellian is far away from the old the derivative must be recalculated:
        !!!! call PrepareFMA_DGV ! this will prepare the linearized operator  fmA = 2\int f_m(v)A(v,v_1,phi) dv
        ! updated the linear operator
        run_mode=1                          ! linear regime
        print *, "Switching to decomposition regime from ES regime"
      else
	    run_mode=0
		print*,"switching to strongly non-linear from ES regime"
	  end if

   case default
      print *, "CheckSolutionMode_DGV: Error, can not process the value of run_mode=", run_mode
      stop
   end select
   !! LATER: MAKE SURE THAT fm is always up to date. Then we can use it, instead of fMaxwellNew. For now comment the next line
   !!  Df=f-fm                   ! evaluate the perturbation from maxwellian

end subroutine CheckSolutionMode_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionModeSimple1Donecell_DGV(f,Df,cellnum,slt) 
! 
! This subroutine is a simplified version of the above subroutine. To keep things simple, all preararaory subrouintes are eliminated 
! when swithcing the mode. If slt is .true. then no diagnostic messaging is created. 
!
! This subroutine evaluates macroparameters of the solution. Then it subtracts 
! from the solution a Maxwellian with the computed macroparameters. This difference is placed into Df
! then a norm of Df is assessed. If |Df| is large, (rum_mode=0) then we are in a strongly non-equilibrium 
! regime and nothing is done. 
! if the |Df| is moderate, then we switch to mode 1 (run_mode=1) then linear and quadratic contributions are computed separately
! if |Df| is small, we switch to run_mode=2 and only the linear part of collision opertator is computed.   
! if |Df| is even smaller == run_mode = 3,4 and the model equations will be used. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckSolutionModeSimple1Donecell_DGV(f,Df,cellnum,slt)

use DGV_commvar, only: run_mode_1D, nodes_u, nodes_v, nodes_w, &
                   decomp_lev, linear_lev, vel_lev, ES_lev, nodes_gwts

use DGV_distributions_mod
use DGV_dgvtools_mod

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out):: Df ! the difference between the solution and the local Maxwellian
integer (I4B), intent(in) :: cellnum ! the number of the cell for which run_mode is calculated.
logical, intent (in) :: slt !if =true, no printout is generated

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: n, u_0,v_0,w_0
real (DP), dimension (size(f,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err, L1_err_fm ! a sctratch variables for the errors
integer (I4B) run_mode ! a scratch variable to keep the run mode.
real (DP) :: LocDens, LocUbar, LocVbar, LocWbar, LocTempr ! scratch variables to keep macroparameters, 
 
call MassCheckRec (f,LocDens,LocUbar,LocVbar,LocWbar,LocTempr) ! First, we evaluate the macroparamters of the solution.
fMaxwellNew = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u, nodes_v, nodes_w) ! now we populate the maxwellian with the same macroparamters.
Df = f-fMaxwellNew ! perturbation from the local maxwellian
L1_err = SUM(ABS(Df)*nodes_gwts)/LocDens! evaluate the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Print May comment
if (.not. slt) then 
 print*, "ChSOlModSimple1D: L1 error = ", L1_err, "cellnum=", cellnum
end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
run_mode = run_mode_1D(cellnum) ! restore the old runmode. 
   select case (run_mode)
   case (0) ! if we're in a strongly non-linear regime on the previous time step, then do this: 
      if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then 
		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to ES-BGK mode from strongly non-linear"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		if (.not. slt) then 
 		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to vel ES-BGK mode from strongly non-linear regime"
 		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2                          ! linear regime
        if (.not. slt) then 
 		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to linearized regime from strongly non-linear"
 		end if 
      elseif (L1_err < decomp_lev) then
	    run_mode=1         
	    if (.not. slt) then             ! decomposition regime
         print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to decomposition regime from strongly non-linear"
        end if  
      end if
!
  case (1) ! if already we're in the decomposition regime
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		if (.not. slt) then
		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to ES-BGK mode from decomposition regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to vel ES-BGK mode from decomposition regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then                          ! linear regime
         print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to linearized regime from decomposition regime"
        end if  
      elseif (L1_err < decomp_lev) then
        run_mode=1  
      else
	    run_mode=0
	    if (.not. slt) then 
		 print*, "ChSOlModSimple1D: cell=",cellnum, "Switching to strongly non-linear from decomposition regime"
		end if 
	  end if
!
   case (2) ! if we're in the linear regime..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"ChSOlModSimple1D: cell=",cellnum, "Switching to ES-BGK mode from linear regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *,"ChSOlModSimple1D: cell=",cellnum, "Switching to vel ES-BGK mode from linear regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Determine the new regime:
        run_mode=2
      elseif (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
         print *,"ChSOlModSimple1D: cell=",cellnum, "Switching to decomposition regime from linear regime"
        end if  
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"ChSOlModSimple1D: cell=",cellnum, "Switching to strongly non-linear from linear regime"
		end if  
	  end if
!
	case (3) ! vel dep ES-BGK mode
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"ChSOlModSimple1D: cell=",cellnum,  "Switching to ES-BGK mode from vel-dependent regime"
		end if 
	  elseif (L1_err < vel_lev) then
	    run_mode=3
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then
		 print *,"ChSOlModSimple1D: cell=",cellnum,  "Switching to linear mode from vel-dependent regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *,"ChSOlModSimple1D: cell=",cellnum,  "Switching to decomposition regime from vel-dep regime"
		end if 
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"ChSOlModSimple1D: cell=",cellnum, "Switching to strongly non-linear from vel-dep regime"
		end if 
	  end if
!
	case (4) ! if we're in the ES-BGK mode..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to vel-dep mode from ES regime"
		end if 
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2  
        if (.not. slt) then
		 print*, "ChSOlModSimple1D: cell=",cellnum, "Switching to linear reg from ES regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *, "ChSOlModSimple1D: cell=",cellnum, "Switching to decomposition regime from ES regime"
		end if  
      else
        run_mode=0
		if (.not. slt) then
		 print*, "ChSOlModSimple1D: cell=",cellnum, "Switching to strongly non-linear from ES regime"
		end if 
	  end if

   case default
      print *, "ChSOlModSimple1D: cell=",cellnum, "Error, can not process the value of run_mode=", run_mode
      stop
   end select
run_mode_1D(cellnum) = run_mode  ! save the new runmode.   

end subroutine CheckSolutionModeSimple1Donecell_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CheckNuUpdateNeededF0D_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err)
! 
! This subrouitne is used to check if the coefficients of the velocity dependent collision frequency need to be updated
! 
! The subroutine uses the current solution, current time, some records and some set paramters to decide if the 
! coefficeints neet to be updated. See description of the update criteria in the documnetation, but currently
! three criteria are implemented: first criteria is time. There is an array in which the time of the next update is recorded.
! second criteria is change in density and third criteria is a change in Temperature. Namely when either density or temprature change on more than 10% 
! as compared to the last stored quantity, the coefficients are updates. Additional criterial for update may be introduced in the future. 
! 
! the solution used in the subroutine uses the primary mesh
! 
! the subrouitine returns a few parameters back:
! updateNulcl == logical -- yes if an update is needed. 
! Df - differnce between f and the local Maxwellian
! LocDens,LocUbar,LocVbar,LocWbar,LocTempr, -- values of the local macroparameters
! L1_err - l1 norm of Df
!
! Another important function of this subroutine is that it is also initializing the arrays that store the computed coefficients. 
! If first checks is the array exists and if not, it will create it. It uses flag MoRatesArry0DAllocated to check is the array exits (true=exists) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckNuUpdateNeededF0D_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err)

use DGV_commvar, only: MaxNumBFunctionsCollFreq, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       MoRatesArry0D, MoRatesArry0DAllocated, & 
                       NuNextUpdateTime0D,NuLastDens0D,NuLastTemp0D, &
                       Order_nu,mft_coeff,pi25DT,mol_diam,L_inf,N_inf
                       
use DGV_distributions_mod                       


logical, intent (out) :: updateNulcl  ! the returned flag if =True -- the coefficient need to be updated
real (DP), intent(in) :: time         ! the current time
real (DP), dimension(:), intent (in) :: f ! the current solution on the primary mesh
real (DP), dimension(:), intent (out) :: Df ! the difference between the current solution and the local Maxwellian
real (DP), intent (out) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr ! values of the local macroparameters
real (DP), intent (out) :: L1_err ! L1-norm of Df

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
real (DP), parameter :: dens_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
real (DP), parameter :: tempr_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we begin from checking if the storage for the coefficients of the velocity dependent collision frequency has been created yet
if (.not. MoRatesArry0DAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (MoRatesArry0D(1:Order_nu), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "CheckNuUpdateNeededF0D_DGV: Error allocation arrays to store moment relax. rates of the VDCF (MoRatesArry0D)"
       stop
      end if 
 MoRatesArry0D = 0.0     
 MoRatesArry0DAllocated = .true.
! Also, we will set up the other arrays to make sure that the coefficients are updated:
NuNextUpdateTime0D = time - 100000.0   ! these bogus values will set off the criteria for update
NuLastDens0D = -100000.0 
NuLastTemp0D = -100000.0 
! 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we will need the values of the local macroparameters. also, the subroutine needs to return a byproducts: local macroparameters, Df and l1 norm of Df. 
! we will compute these quantities next: 
! compute mass
LocDens = sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! compute momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocUbar = sum(f*nodes_gwts*nodes_u)/LocDens
LocVbar = sum(f*nodes_gwts*nodes_v)/LocDens
LocWbar = sum(f*nodes_gwts*nodes_w)/LocDens
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocTempr = sum(f*nodes_gwts*((nodes_u-LocUbar)**2+(nodes_v-LocVbar)**2+(nodes_w-LocWbar)**2))/LocDens/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will calculate the local maxwellian, keep it in Df
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Df = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u,nodes_v,nodes_w) ! now we have the maxwellian with the same macroparamters.
!!!!! 
! Now we compute (f-f_{M}): 
Df = f-Df ! now Df keeps the difference
! finally, we need to compute the l1_norm of the Df: 
L1_err = SUM(ABS(Df)*nodes_gwts)/LocDens 
! Done computing auxiliary staff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!
! now we have all variable setup and used at least once, so we will start applying the update criteria:
updateNulcl = .false. ! reset the update flag
! first we check the time of the next update. If we past that time, we need to update: 
if (time > NuNextUpdateTime0D) then 
 updateNulcl = .true.
 ! also, we now need to set the next update time. 
 ! the next update time is obtained by adding a multple of the local mean free time to the 
 ! current time. We use the dimensionless mean free time for this operation:  
 NuNextUpdateTime0D = time + & 
    mft_coeff/LocDens/sqrt(LocTempr)*(L_inf**5)/(N_inf*mol_diam*mol_diam)/4.0/sqrt(pi25DT)*sqrt(2.0_DP)  
end if
! next we check the last saved values of the density 
if (ABS(LocDens - NuLastDens0D) > NuLastDens0D*dens_trshld) then 
 updateNulcl = .true.
 NuLastDens0D = LocDens
end if  
! next we check the last saved values of the temperature 
if (ABS(LocTempr - NuLastTemp0D) > NuLastTemp0D*tempr_trshld) then 
 updateNulcl = .true.
 NuLastTemp0D = LocTempr
end if  
!!!! all done! 
end subroutine CheckNuUpdateNeededF0D_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CheckNuUpdateNeededF0Da_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err)
! 
! A copy of the above subroutine with the exeption that Df - f-f_M is already provided -- so no need to compute it.
! (Df - differnce between f and the local Maxwellian)
!
! This subrouitne is used to check if the coefficients of the velocity dependent collision frequency need to be updated
! 
! The subroutine uses the current solution, current time, some records and some set paramters to decide if the 
! coefficeints neet to be updated. See description of the update criteria in the documnetation, but currently
! three criteria are implemented: first criteria is time. There is an array in which the time of the next update is recorded.
! second criteria is change in density and third criteria is a change in Temperature. Namely when either density or temprature change on more than 10% 
! as compared to the last stored quantity, the coefficients are updates. Additional criterial for update may be introduced in the future. 
! 
! the solution used in the subroutine uses the primary mesh
! 
! the subrouitine returns a few parameters back:
! updateNulcl == logical -- yes if an update is needed. 
! LocDens,LocUbar,LocVbar,LocWbar,LocTempr, -- values of the local macroparameters
! L1_err - l1 norm of Df
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckNuUpdateNeededF0Da_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err)

use DGV_commvar, only: MaxNumBFunctionsCollFreq, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       NuNextUpdateTime0D,NuLastDens0D,NuLastTemp0D, &
                       Order_nu,mft_coeff,pi25DT,mol_diam,L_inf,N_inf
                       
use DGV_distributions_mod                       


logical, intent (out) :: updateNulcl  ! the returned flag if =True -- the coefficient need to be updated
real (DP), intent(in) :: time         ! the current time
real (DP), dimension(:), intent (in) :: f ! the current solution on the primary mesh
real (DP), dimension(:), intent (in) :: Df ! the difference between the current solution and the local Maxwellian
real (DP), intent (out) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr ! values of the local macroparameters
real (DP), intent (out) :: L1_err ! L1-norm of Df

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

real (DP), parameter :: dens_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
real (DP), parameter :: tempr_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we will need the values of the local macroparameters. also, the subroutine needs to return a byproducts: local macroparameters, Df and l1 norm of Df. 
! we will compute these quantities next: 
! compute mass
LocDens = sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! compute momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocUbar = sum(f*nodes_gwts*nodes_u)/LocDens
LocVbar = sum(f*nodes_gwts*nodes_v)/LocDens
LocWbar = sum(f*nodes_gwts*nodes_w)/LocDens
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocTempr = sum(f*nodes_gwts*((nodes_u-LocUbar)**2+(nodes_v-LocVbar)**2+(nodes_w-LocWbar)**2))/LocDens/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finally, we need to compute the l1_norm of the Df: 
L1_err = SUM(ABS(Df)*nodes_gwts)/LocDens 
! Done computing auxiliary staff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!
! now we have all variable setup and used at least once, so we will start applying the update criteria:
updateNulcl = .false. ! reset the update flag
! first we check the time of the next update. If we past that time, we need to update: 
if (time > NuNextUpdateTime0D) then 
 updateNulcl = .true.
 ! also, we now need to set the next update time. 
 ! the next update time is obtained by adding a multple of the local mean free time to the 
 ! current time. We use the dimensionless mean free time for this operation:  
 NuNextUpdateTime0D = time + & 
    mft_coeff/LocDens/sqrt(LocTempr)*(L_inf**5)/(N_inf*mol_diam*mol_diam)/4.0/sqrt(pi25DT)*sqrt(2.0_DP)  
end if
! next we check the last saved values of the density 
if (ABS(LocDens - NuLastDens0D) > NuLastDens0D*dens_trshld) then 
 updateNulcl = .true.
 NuLastDens0D = LocDens
end if  
! next we check the last saved values of the temperature 
if (ABS(LocTempr - NuLastTemp0D) > NuLastTemp0D*tempr_trshld) then 
 updateNulcl = .true.
 NuLastTemp0D = LocTempr
end if  
!!!! all done! 
end subroutine CheckNuUpdateNeededF0Da_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CheckNuUpdateNeededF1Donecella_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err,cellnum)
! 
! A copy of the above subroutine to work with multidimensional solutions. cellnum is the number of the spatial cell
!
! A copy of the above subroutine with the exeption that Df - f-f_M is already provided -- so no need to comput it.
! (Df - differnce between f and the local Maxwellian)
!
! This subrouitne is used to check if the coefficients of the velocity dependent collision frequency need to be updated
! 
! The subroutine uses the current solution, current time, some records and some set paramters to decide if the 
! coefficeints neet to be updated. See description of the update criteria in the documnetation, but currently
! three criteria are implemented: first criteria is time. There is an array in which the time of the next update is recorded.
! second criteria is change in density and third criteria is a change in Temperature. Namely when either density or temprature change on more than 10% 
! as compared to the last stored quantity, the coefficients are updates. Additional criterial for update may be introduced in the future. 
! 
! the solution used in the subroutine uses the primary mesh
! 
! the subrouitine returns a few parameters back:
! updateNulcl == logical -- yes if an update is needed. 
! LocDens,LocUbar,LocVbar,LocWbar,LocTempr, -- values of the local macroparameters
! L1_err - l1 norm of Df
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckNuUpdateNeededF1Donecella_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err,cellnum)

use DGV_commvar, only: MaxNumBFunctionsCollFreq, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D, &
                       Order_nu,mft_coeff,pi25DT,mol_diam,L_inf,N_inf
                       
use DGV_distributions_mod                       


logical, intent (out) :: updateNulcl  ! the returned flag if =True -- the coefficient need to be updated
real (DP), intent(in) :: time         ! the current time
real (DP), dimension(:), intent (in) :: f ! the current solution on the primary mesh
real (DP), dimension(:), intent (in) :: Df ! the difference between the current solution and the local Maxwellian
real (DP), intent (out) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr ! values of the local macroparameters
real (DP), intent (out) :: L1_err ! L1-norm of Df
integer (I4B), intent (in) :: cellnum ! the numbder of the spatial cell for which the need for update is verified

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

real (DP), parameter :: dens_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
real (DP), parameter :: tempr_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we will need the values of the local macroparameters. also, the subroutine needs to return a byproducts: local macroparameters, Df and l1 norm of Df. 
! we will compute these quantities next: 
! compute mass
LocDens = sum(f*nodes_gwts)
!!!!!!!!!!!!!!!!!!!!!!!!!
! compute momentum 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocUbar = sum(f*nodes_gwts*nodes_u)/LocDens
LocVbar = sum(f*nodes_gwts*nodes_v)/LocDens
LocWbar = sum(f*nodes_gwts*nodes_w)/LocDens
! check temperature 
!!!!!!!!!!!!!!!!!!!!!!!!!
LocTempr = sum(f*nodes_gwts*((nodes_u-LocUbar)**2+(nodes_v-LocVbar)**2+(nodes_w-LocWbar)**2))/LocDens/3.0_DP*2.0_DP ! dimensionless temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! finally, we need to compute the l1_norm of the Df: 
L1_err = SUM(ABS(Df)*nodes_gwts)/LocDens 
! Done computing auxiliary staff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!
! now we have all variable setup and used at least once, so we will start applying the update criteria:
updateNulcl = .false. ! reset the update flag
! first we check the time of the next update. If we past that time, we need to update: 
if (time > NuNextUpdateTime1D(cellnum)) then 
 updateNulcl = .true.
 ! also, we now need to set the next update time. 
 ! the next update time is obtained by adding a multple of the local mean free time to the 
 ! current time. We use the dimensionless mean free time for this operation:  
 NuNextUpdateTime1D(cellnum) = time + & 
    mft_coeff/LocDens/sqrt(LocTempr)*(L_inf**5)/(N_inf*mol_diam*mol_diam)/4.0/sqrt(pi25DT)*sqrt(2.0_DP)  
end if
! next we check the last saved values of the density 
if (ABS(LocDens - NuLastDens1D(cellnum)) > NuLastDens1D(cellnum)*dens_trshld) then 
 updateNulcl = .true.
 NuLastDens1D(cellnum) = LocDens
end if  
! next we check the last saved values of the temperature 
if (ABS(LocTempr - NuLastTemp1D(cellnum)) > NuLastTemp1D(cellnum)*tempr_trshld) then 
 updateNulcl = .true.
 NuLastTemp1D(cellnum) = LocTempr
end if  
!!!! all done! 
end subroutine CheckNuUpdateNeededF1Donecella_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CheckNuUpdateNeededF2Donecella_DGV(updateNulcl,time,f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,L1_err,cellnum)
!
! another modification of the above suibroutine, in this version macroparamters are provided to the subroutine
! 
! A copy of the above subroutine to work with multidimensional solutions. cellnum is the number of the spatial cell
!
! A copy of the above subroutine with the exeption that Df - f-f_M is already provided -- so no need to comput it.
! (Df - differnce between f and the local Maxwellian)
!
! This subrouitne is used to check if the coefficients of the velocity dependent collision frequency need to be updated
! 
! The subroutine uses the current solution, current time, some records and some set paramters to decide if the 
! coefficeints neet to be updated. See description of the update criteria in the documnetation, but currently
! three criteria are implemented: first criteria is time. There is an array in which the time of the next update is recorded.
! second criteria is change in density and third criteria is a change in Temperature. Namely when either density or temprature change on more than 10% 
! as compared to the last stored quantity, the coefficients are updates. Additional criterial for update may be introduced in the future. 
! 
! the solution used in the subroutine uses the primary mesh
! 
! the subrouitine returns a few parameters back:
! updateNulcl == logical -- yes if an update is needed. 
! LocDens,LocUbar,LocVbar,LocWbar,LocTempr, -- values of the local macroparameters
! L1_err - l1 norm of Df
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckNuUpdateNeededF2Donecella_DGV(updateNulcl,time,LocDens,LocTempr,cellnum)

use DGV_commvar, only: MaxNumBFunctionsCollFreq, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D, &
                       Order_nu,mft_coeff,pi25DT,mol_diam,L_inf,N_inf
                       
use DGV_distributions_mod                       


logical, intent (out) :: updateNulcl  ! the returned flag if =True -- the coefficient need to be updated
real (DP), intent(in) :: time         ! the current time
real (DP), intent (in) :: LocDens,LocTempr ! values of the local macroparameters
integer (I4B), intent (in) :: cellnum ! the numbder of the spatial cell for which the need for update is verified

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

real (DP), parameter :: dens_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
real (DP), parameter :: tempr_trshld = 0.1 ! the fraction of the last saved density that serves as the threshhold for density triggered update
 
!!!!!!!!!!!!!!!!!
! now we have all variable setup and used at least once, so we will start applying the update criteria:
updateNulcl = .false. ! reset the update flag
! first we check the time of the next update. If we past that time, we need to update: 
if (time > NuNextUpdateTime1D(cellnum)) then 
 updateNulcl = .true.
 ! also, we now need to set the next update time. 
 ! the next update time is obtained by adding a multple of the local mean free time to the 
 ! current time. We use the dimensionless mean free time for this operation:  
 NuNextUpdateTime1D(cellnum) = time + & 
    mft_coeff/LocDens/sqrt(LocTempr)*(L_inf**5)/(N_inf*mol_diam*mol_diam)/4.0/sqrt(pi25DT)*sqrt(2.0_DP)  
end if
! next we check the last saved values of the density 
if (ABS(LocDens - NuLastDens1D(cellnum)) > NuLastDens1D(cellnum)*dens_trshld) then 
 updateNulcl = .true.
 NuLastDens1D(cellnum) = LocDens
end if  
! next we check the last saved values of the temperature 
if (ABS(LocTempr - NuLastTemp1D(cellnum)) > NuLastTemp1D(cellnum)*tempr_trshld) then 
 updateNulcl = .true.
 NuLastTemp1D(cellnum) = LocTempr
end if  
!!!! all done! 
end subroutine CheckNuUpdateNeededF2Donecella_DGV

end module DGV_collision_mod
