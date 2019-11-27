!
!  spat_oper_mod.f90
!
! This module contains routines that are involved in the evaluation of the spatial operator
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module BGKVDCF0D_spat_oper_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SpatialOperSH_DGV (f)
!
! This operator evaluates the spatial operator 
! 
!
! in 0D / the spatially homogeneous case the spatial operator 
! only has the collisoin operator
!
! in 1D and 2D there will be also operators evaluating thr transport part, fluxes and so on 
!
! COMMENT: STILL NEED TO ADD COLLISION FREQUENCY!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SpatialOperSH_DGV (f,fcol,dt)

use BGKVDCF0D_commvar, only: time
use DGV_collision_mod

!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (:), intent (out) :: fcol ! value of the right side  
real (DP), intent (in) :: dt ! The time step 
!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Since there is no transport terms, flux terms and so on, we just call the 
!! right hand side, which in this case is accomplished by calling the 
!! collision operator

call UniversalCollisionOperator0DDGV(f,fcol,time,dt)

!!!

end subroutine SpatialOperSH_DGV



end module BGKVDCF0D_spat_oper_mod 