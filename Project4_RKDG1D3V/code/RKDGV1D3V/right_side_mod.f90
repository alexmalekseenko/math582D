!
!  This module implements the BGK collision integral 
! 
!
module right_side_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

real (DP), parameter, private :: kBoltzmannconst = 1.380658D-23
real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RightSide_HOV1D (ftil1,tau)
! 
! This function calculates the right side of the 1D BKG model: (see notes)
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RightSide1D3D_HOV1D (r1,ftil1,tau,time)

use common_variables_mod, only: moments_x_gauss_order, moments_refine_x,&
                                XWeightLegPol,XLegPol
                      
use spectral_tools_mod

!!!!!!!!! link to the DGVlib

use DGV_collision_mod


!!!!!!!!!


!!! main variables
real (DP), dimension (0:,:,:), intent (out) :: r1 ! values of the right side multiplied by \tau, 
                                   ! have the same dimensions as the main variable
                                   ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- j is the cell in x
real (DP), dimension (0:,:,:), intent (in) :: ftil1 ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- j is the cell in x

real (DP), intent (in) :: tau      ! time step parameter. Is used in time integration. 
real (DP), intent (in) :: time     ! tie value of the current time variable
!!! auxiliary variables: 
integer (I4B) :: i,j,p,q ! some local counters
! real (DP) :: dx ! some local variables 
integer :: loc_alloc_stat ! to keep allocation status		

real (DP), dimension (:,:), allocatable :: RS1c ! RS1c(l,p) - temp variable to store values of the right side of the kinetic equation on one cell
                                ! p --- number of the node on the cell (uncluding the refinement in x)
						        ! l -- node in velocity

real (DP), dimension (:), allocatable :: fcol,f ! temp variable to store \value of the collision operator at a single point x 
						        
integer (I4B) :: num_x_node ! scrap variable to keep the number of the nodal point x where solution is evaluated to compute galerkin coefficients()

!!!!!!!!!!!!!!!!!
r1=0 ! reset the right side
!!!
! creating the scrap variables -- we assume that all cells in x will have the same amount of Galerkin coefficients.
! and that the same amount of Gauss nodes is used for computing Galerkin coefficients on each cell. 
allocate (RS1c(moments_x_gauss_order*moments_refine_x,size(ftil1,2)), stat=loc_alloc_stat)
 !
 if (loc_alloc_stat >0) then 
  print *, "RightSide1D3D_HOV1D: Allocation error for variable (RS1c)"
 end if 
 !
allocate (f(size(ftil1,2)),fcol(size(ftil1,2)), stat=loc_alloc_stat)
 !
 if (loc_alloc_stat >0) then 
  print *, "RightSide1D3D_HOV1D: Allocation error for variable (f), (fcol)"
 end if 
 !
do j = 1,size(ftil1,3) ! we will have a loop in all cells in x
 ! first,we assemble the solution at all Gauss nodes in x at all velocity nodes in u
 RS1c = 0 ! reset the value of the solution on one cell:
 do i = 1,size(ftil1,2) ! loop in all velocity nodes
  do q = 1,moments_x_gauss_order*moments_refine_x ! loop in points x on the cell 
   do p = 0,size(ftil1,1)-1 ! loop in galerkin coefficients
    !! Evaluate the basis functions \varhi(x) on xx 
    RS1c(q,i) = RS1c(q,i) + XLegPol(q,p)*ftil1(p,i,j)
   end do !end loop in Galerkin coefficients
  end do ! end loop in points x on the cell
 end do ! end loop in all velocity nodes 
 ! the solution is assembled at all points x on one cell in x from the Galerkin coefficients! 
 ! now, we can pass the solutoin to the subroutine evaluating the right hand side.
 do q = 1,moments_x_gauss_order*moments_refine_x ! loop in points x on the cell ! size(Xnodes_unit,1) = moments_x_gauss_order*moments_refine_x
  num_x_node = (j-1)*moments_x_gauss_order*moments_refine_x + q 
  do i = 1,size(ftil1,2) ! loop in all velocity nodes
   f(i) = RS1c(q,i);  ! place solution at one point in x in a temp variable
  end do ! end loop in all velocity nodes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CALL DGVlib subroutine 
  call UniversalCollisionOperator1DonecellDGV(f,fcol,time,tau,num_x_node) ! fcol returns the collision operator multiplied by tau
        ! also, a lot of constants, including the gas constants and the dimensionless reduction constants are used from DGVcommvar 
        ! the constants need to be set up to work properly... 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,size(ftil1,2) ! loop in all velocity nodes
   RS1c(q,i) = fcol(i) ! the value of the collision integral is saved in the temp array 
  end do ! end loop in all velocity nodes
 end do ! end loop in points x on the cell 
 ! lastly, we compute the Galerkin coefficients of the right hand side:
 do i = 1,size(ftil1,2) ! loop in all velocity nodes
  do p = 0,size(ftil1,1)-1 ! loop in galerkin coefficients
   r1(p,i,j) = sum(XWeightLegPol(:,p)*RS1c(:,i))*(2*p+1)/Real(2*moments_refine_x,DP)  ! the right side must have the collision frequency and the time step in it already 
  end do  ! end loop in Galerkin coefficients
 end do ! end loop in all velocity nodes
end do ! end of the loop in all cells in x
!! Kinda done ... 
deallocate (RS1c,f,fcol)

end subroutine RightSide1D3D_HOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RightSide1D3DOMP_HOV1D (ftil1,tau)
!
! This is a version of the above subroutine parallellized using OpenMP
!
! 
! This function calculates the right side of the 1D BKG model: (see notes)
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine RightSide1D3DOMP_HOV1D (r1,ftil1,tau,time)

use common_variables_mod, only: moments_x_gauss_order, moments_refine_x,&
                                XWeightLegPol,XLegPol, Num_OMP_threads
                      
use spectral_tools_mod

!!!!!!!!! link to the DGVlib

use DGV_collision_mod

!!!!!!!!!


!!! main variables
real (DP), dimension (0:,:,:), intent (out) :: r1 ! values of the right side multiplied by \tau, 
                                   ! have the same dimensions as the main variable
                                   ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- j is the cell in x
real (DP), dimension (0:,:,:), intent (in) :: ftil1 ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- m_count is the index in the basis functions in u
                                   ! -- j is the cell in x

real (DP), intent (in) :: tau      ! time step parameter. Is used in time integration. 
real (DP), intent (in) :: time     ! tie value of the current time variable
!!! auxiliary variables: 
integer (I4B) :: i,j,p,q ! some local counters
! real (DP) :: dx ! some local variables 
integer :: loc_alloc_stat ! to keep allocation status		

real (DP), dimension (:,:), allocatable  :: RS1c ! RS1c(l,p) - temp variable to store values of the right side of the kinetic equation on one cell
                                ! p --- number of the node on the cell (uncluding the refinement in x)
						        ! l -- node in velocity

real (DP), dimension (:), allocatable  :: fcol,f ! temp variable to store \value of the collision operator at a single point x 
						        
integer (I4B) :: num_x_node ! scrap variable to keep the number of the nodal point x where solution is evaluated to compute galerkin coefficients()
!!!!!!!!!!!!!!!!!
r1=0 ! reset the right side
!!!

!!!
! Open MP Directives 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OpenMP set the number of threads: 
!!!!! !$OMP call omp_set_num_threads(Num_OMP_threads)

!$OMP PARALLEL DEFAULT (shared) PRIVATE(j,i,q,p,RS1c,num_x_node,f,fcol) NUM_THREADS(Num_OMP_threads) 
allocate (fcol(size(ftil1,2)),f(size(ftil1,2)),&
        RS1c(moments_x_gauss_order*moments_refine_x,size(ftil1,2)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "RightSide1D3DOMP_HOV1D: Allocation error for variable (fcol,f,RS1c)"
 stop
end if

!$OMP DO SCHEDULE(DYNAMIC,1) 
do j = 1,size(ftil1,3) ! we will have a loop in all cells in x
 ! first,we assemble the solution at all Gauss nodes in x at all velocity nodes in u
 RS1c = 0 ! reset the value of the solution on one cell:
 do i = 1,size(ftil1,2) ! loop in all velocity nodes
  do q = 1,moments_x_gauss_order*moments_refine_x ! loop in points x on the cell 
   do p = 0,size(ftil1,1)-1 ! loop in galerkin coefficients
    !! Evaluate the basis functions \varhi(x) on xx 
    RS1c(q,i) = RS1c(q,i) + XLegPol(q,p)*ftil1(p,i,j)
   end do !end loop in Galerkin coefficients
  end do ! end loop in points x on the cell
 end do ! end loop in all velocity nodes 
 ! the solution is assembled at all points x on one cell in x from the Galerkin coefficients! 
 ! now, we can pass the solutoin to the subroutine evaluating the right hand side.
 do q = 1,moments_x_gauss_order*moments_refine_x ! loop in points x on the cell ! size(Xnodes_unit,1) = moments_x_gauss_order*moments_refine_x
  num_x_node = (j-1)*moments_x_gauss_order*moments_refine_x + q 
  do i = 1,size(ftil1,2) ! loop in all velocity nodes
   f(i) = RS1c(q,i);  ! place solution at one point in x in a temp variable
  end do ! end loop in all velocity nodes
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! CALL DGVlib subroutine 
  call UniversalCollisionOperator1DonecellDGV(f,fcol,time,tau,num_x_node) ! fcol returns the collision operator multiplied by tau
        ! and by a number of constants, including the gas constants and the dimensionless reduction constants are used from DGVcommvar 
        ! the constants need to be set up correctly for the method to work properly... 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i = 1,size(ftil1,2) ! loop in all velocity nodes
   RS1c(q,i) = fcol(i) ! the value of the collision integral is saved in the temp array 
  end do ! end loop in all velocity nodes
 end do ! end loop in points x on the cell 
 ! lastly, we compute the Galerkin coefficients of the right hand side:
 do i = 1,size(ftil1,2) ! loop in all velocity nodes
  do p = 0,size(ftil1,1)-1 ! loop in galerkin coefficients
   r1(p,i,j) = sum(XWeightLegPol(:,p)*RS1c(:,i))*(2*p+1)/Real(2*moments_refine_x,DP)  ! the right side must have the collision frequency and the time step in it already 
  end do  ! end loop in Galerkin coefficients
 end do ! end loop in all velocity nodes
end do ! end of the loop in all cells in x
!$OMP END DO
deallocate (RS1c,f,fcol)
!$OMP END PARALLEL 
!! Kinda done ... 

end subroutine RightSide1D3DOMP_HOV1D


end module right_side_mod