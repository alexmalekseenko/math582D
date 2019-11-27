!
!  time_integrator_mod.f90
!
!  
!  8/13/2008 3:19:01 PM
!
! This module contains subroutines that implement time integration
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module time_integrator_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TimeIntegratorORK_HOV1D3D
! 
! This subroutine implements Optimal Runge Kutta integration. 
! 
! It is dependent on the main program: uses variables from the 
! common_variables_mod 
!
! the called subroutines look up varaibles in common_variables_mod 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine TimeIntegratorORK_HOV1D3D


use common_variables_mod, only: curr_time,dt,rk,ftil,frhs1
use spatial_operator_mod
!
real (DP), dimension (:,:,:), allocatable :: ft1_1, ft1_2, ft1_3, ft1_4 ! auxiliary variables to keep intermediate time steps
integer :: loc_alloc_stat ! to keep allocation status
!
!!!!!!!!!!!!!!!!

select case (rk)
 case (1) 
  allocate (ft1_1(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_1)"
     end if 
     !
  !Runger Kutta first order (Euler method)
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_1,ftil,curr_time,dt)     
   frhs1 = ft1_1/dt
   ftil = ftil + ft1_1
  ! end Runge Kutta 1st order (Euler method).   
  deallocate(ft1_1)
 case (2)
  allocate (ft1_1(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_1)"
     end if 
     !
  allocate (ft1_2(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_2)"
     end if 
     !
   !Runge Kutta second order (optimal 2nd order method)
  ! t_1:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_1,ftil,curr_time,dt)
   frhs1 = ft1_1/dt
  ! t_2: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil+ft1_1*2/3_DP,&
                      curr_time+2*dt/3_DP,dt)
  !                    
   ftil = ftil + ft1_1/4_DP + ft1_2*3/4_DP
  ! end optimal Runge Kutta second order 
  deallocate (ft1_1, ft1_2)
 case (3)
  allocate (ft1_1(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_1)"
     end if 
     !
  allocate (ft1_2(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_2)"
     end if 
     !
  allocate (ft1_3(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_3)"
     end if 
     !
  ! t_1:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_1,ftil,curr_time,dt)
   frhs1 = ft1_1/dt
  ! t_2: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil+ft1_1/2_DP,&
                      curr_time+dt/2_DP,dt)
  ! t_3:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_3,ftil+ft1_2*3/4_DP,&
                      curr_time+dt*3/4_DP,dt)     
  !
   ftil = ftil + ft1_1*2/9.0_DP + ft1_2*3/9.0_DP + ft1_3*4/9.0_DP
   ! end optimal Runge-Kutta third order
  deallocate(ft1_1, ft1_2, ft1_3)
 case (4)
  allocate (ft1_1(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_1)"
     end if 
     !
  allocate (ft1_2(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_2)"
     end if 
     !
  allocate (ft1_3(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_3)"
     end if 
     !
  ! (Optimal ?) Runge Kutta fourth order (Ralston and Rabinowitz ??)
  ! t_1:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_1,ftil,curr_time,dt)
   frhs1 = ft1_1/dt
  ! t_2: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil+ft1_1/2_DP,&
                      curr_time+dt/2_DP,dt)
  ! t_3:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_3,ftil+ft1_1/4_DP+ft1_2/4_DP,&
                      curr_time+dt/2_DP,dt)     
  ! t_4: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil-ft1_2+ft1_3*2.0_DP,&
                      curr_time+dt,dt)     
  !
   ftil = ftil + ft1_1/6.0_DP + ft1_3*2/3.0_DP + ft1_2/6.0_DP  ! :2nd componenent is re-used for 4th component!
  ! end Runge-Kutta fourth order (Ralston and Rabinowitz ??) 
  deallocate (ft1_1,ft1_2,ft1_3)
 case (5)
 allocate (ft1_1(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_1)"
     end if 
     !
 allocate (ft1_2(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_2)"
     end if 
     !
 allocate (ft1_3(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_3)"
     end if 
     !
 allocate (ft1_4(size(ftil,1),size(ftil,2),size(ftil,3)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "TimeIntegratorORK_HOV1D3D: Allocation error for variable (ft1_4)"
     end if 
     !
  ! Runge Kutta fifth order (Ralston and Rabinowitz ?? compatible with the fourth order)
  ! t_1:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_1,ftil,curr_time,dt)
   frhs1 = ft1_1/dt
  ! t_2: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil+ft1_1/2_DP,&
                      curr_time+dt/2_DP,dt)
  ! t_3:
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_3,ftil+ft1_1/4_DP+ft1_2/4_DP,&
                      curr_time+dt/2_DP,dt)     
  ! t_4: 
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_4,ftil-ft1_2+ft1_3*2,&
                      curr_time+dt,dt)     
  ! we now introduce some memory savings and replace ft1_3 with -ft1_2/5_DP+ft1_3*546/625_DP that is going to be ised in step 6
   ft1_3 = -ft1_2/5_DP+ft1_3*546/625.0_DP
   !t_5:
  ! on this step, we re-use ft1_2 to take place of ft1_5
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_2,ftil+ft1_1*7/27.0_DP+ft1_2*10/27.0_DP+ft1_4/27.0_DP,&
                      curr_time+dt*2/3_DP,dt)
  !t_6: 
  ! on this step, we re-use ft1_3 to take place of ft1_6
   call SpatialOperatorHOV1D3D_UpwindFlux(ft1_3,ftil+ft1_1*28/625.0_DP+ft1_3+ft1_4*54/625.0_DP - ft1_2*378/625.0_DP,&
                      curr_time+dt/5_DP,dt)
  !
   ftil = ftil + ft1_1/24_DP + ft1_4*5/48_DP + ft1_2*27/56_DP + ft1_3*125/336_DP
   deallocate (ft1_1,ft1_2,ft1_3,ft1_4) 
  case default
           print *, "unsupported degree of Runge-Kutta method passed"
 end select 
!
end subroutine TimeIntegratorORK_HOV1D3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TimeIntegratorMTS_SH_DGV
!
! Modified by  Alex 10/21/2011
! 
! This subroutine implements Multiple Step Time integration. 
! 
! It is dependent on the main program: uses variables from the 
! common_variables_mod 
!
! the called subroutines look up varaibles in common_variables_mod 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

SUBROUTINE TimeIntegratorMTS1D3D

use common_variables_mod, only: ftil,ftil1,ftil2,ftil3,ftil4,frhs1,frhs2,frhs3,frhs4, &
                                frhs5,curr_time,dt,rk
use spatial_operator_mod

!!!!!!! DEBUG 
use DGV_commvar, only: run_mode_1D
!!!!!!!! END DEBUG

!!!!!!!!

SELECT CASE (rk)
	CASE (1)
	 !!! MTS FIRST ORDER (EULER METHOD)
	 CALL SpatialOperatorHOV1D3D_UpwindFlux(frhs1,ftil,curr_time,dt)     
     !
     ftil = ftil + frhs1 
       !
	 !!! END MTS 1ST ORDER (EULER METHOD)
	CASE (2)
	 !!! ADAMS-BASHFORTH TWO STEP EXPLICIT METHOD
	 ftil1 = ftil; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperatorHOV1D3D_UpwindFlux(frhs1,ftil,curr_time,dt)
	 ftil = ftil + frhs1*3.0_DP/2.0_DP - dt*frhs2/2.0_DP
	 !!! END ADAMS-BASHFORTH TWO STEP EXPLICIT METHOD
	CASE (3)
	 !!! BEGIN ADAMS-BASHFORTH THREE STEP
	 ftil2 = ftil1; ftil1 = ftil; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperatorHOV1D3D_UpwindFlux(frhs1,ftil,curr_time,dt) 
	 ftil = ftil + frhs1*23.0_DP/12.0_DP - dt*frhs2*4/3.0_DP + dt*frhs3*5.0_DP/12.0_DP 
    !!! END ADAMS-BASHFORTH THREE STEP
	CASE (4)
   	 !!! ADAMS-BASHFORTH FOUR STEP
	 ftil3 = ftil2; ftil2 = ftil1; ftil1 = ftil; frhs4 = frhs3; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperatorHOV1D3D_UpwindFlux(frhs1,ftil,curr_time,dt) 
	 ftil = ftil + frhs1*55.0_DP/24.0_DP - dt*frhs2*59.0_DP/24.0_DP + dt*frhs3*37.0_DP/24.0_DP - dt*frhs4*3.0_DP/8
	 !!! END ADAMS-BASHFORTH FOUR STEP
	CASE (5)
	 !!! ADAMS-BASHFORTH FIVE STEP METHOD
	 ftil4 = ftil3; ftil3 = ftil2; ftil2 = ftil1; ftil1 = ftil 
	 frhs5 = frhs4; frhs4 = frhs3; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperatorHOV1D3D_UpwindFlux(frhs1,ftil,curr_time,dt) 
	 ftil = ftil + frhs1*(1901.0_DP/720.0_DP) - dt*frhs2*(2774.0_DP/720.0_DP) + dt*frhs3*(2616.0_DP/720.0_DP) &
	       - dt*frhs4*(1274.0_DP/720.0_DP) + dt*frhs5*(251.0_DP/720.0_DP)
	 !!! END ADAMS-BASHFORTH FIVE STEP METHOD
    CASE default
	 PRINT *, "TimeIntegratorMTS1D3D: Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT
    frhs1 = frhs1/dt
    curr_time = curr_time + dt ! advance the time for one step
    
    !!!!!! DEBUG 
 !do iii=1,size(run_mode_1D,1)
  print *, "time=", curr_time, "run_mode_1D=", run_mode_1D
 !end do 
 !!!!!!! END DEBUG
    
END SUBROUTINE TimeIntegratorMTS1D3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine PrepareMTS_SH_DGV
!
! Subroutine prepares the intermediate time steps to be used in the multi-step methods.
! f stores the actual solution.
! f1, f2, f3, f4, frhs1, frhs2, frhs3, frhs4 store the value of the 
!	spatial operator and the right hand side on those intermediate time steps 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PrepareMTS1D3D

use common_variables_mod, only: ftil,ftil1,ftil2,ftil3,ftil4,frhs1,frhs2,frhs3,frhs4, &
                                frhs5,curr_time,dt,rk
!!!!!!! DEBUG 
use DGV_commvar, only: run_mode_1D
!!!!!!!! END DEBUG


!!!!!! DECLARE LOCAL VARIABLES !!!!!
INTEGER :: i,s1,s2,s3													! Local counters
INTEGER :: loc_alloc_stat ! to keep the allocation status 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! auxiliary numbers -- sizes of the arrays
s1=size(ftil,1)-1
s2=size(ftil,2)
s3=size(ftil,3)

! First, we need to allocate the storages
SELECT CASE(rk)
	CASE(1)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
	  ALLOCATE(frhs1(0:s1,1:s2,1:s3),ftil1(0:s1,1:s2,1:s3),	stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs1), (ftil1)"
		stop
	  END IF
	CASE(2)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:s1,1:s2,1:s3),frhs2(0:s1,1:s2,1:s3),ftil1(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs1,frsh2,ftil1)"
		stop
	  END IF
	CASE(3)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:s1,1:s2,1:s3),ftil1(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs1,ftil1)"
		stop
	  END IF
	  ALLOCATE(frhs2(0:s1,1:s2,1:s3),frhs3(0:s1,1:s2,1:s3),ftil2(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs2,ftil2,frhs3)"
		stop
	  END IF
	CASE(4) 
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:s1,1:s2,1:s3),ftil1(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs1,ftil1)"
		stop
	  END IF
	  ALLOCATE(frhs2(0:s1,1:s2,1:s3),ftil2(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs2,ftil2)"
		stop
	  END IF
	  ALLOCATE(frhs3(0:s1,1:s2,1:s3),frhs4(0:s1,1:s2,1:s3),ftil3(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs3,ftil3)"
		stop
	  END IF
    CASE(5)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:s1,1:s2,1:s3),ftil1(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs1,ftil1)"
		stop
	  END IF
	  ALLOCATE(frhs2(0:s1,1:s2,1:s3),ftil2(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs2,ftil2)"
		stop
	  END IF
	ALLOCATE(frhs3(0:s1,1:s2,1:s3),ftil3(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs3,ftil3)"
		stop
	  END IF
	  ALLOCATE(frhs4(0:s1,1:s2,1:s3),frhs5(0:s1,1:s2,1:s3),ftil4(0:s1,1:s2,1:s3), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS1D3D: Allocation error for variable (frhs4,ftil4)"
		stop
	  END IF
	CASE default
			PRINT *, "PrepareMTS1D3D: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
			stop
END SELECT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we assume that f contains the intial data. We need to set up the arrays
! If rkmts=1 - no action is needed all data is in place 
! If the rkmts>2 then we need to actually prepare some arrays. 
do i=2,rk
 !! FIRST WE NEED TO PUSH THE SOLUTION AND THE RIGHT SIDE DOWN THE STORAGE TO OPEN SOME SPACE FOR THE NEW DATA 
 SELECT CASE(rk)
 CASE(2)
  ftil1=ftil 
 CASE(3)
  ftil2=ftil1; ftil1=ftil; frhs2=frhs1
 CASE(4)
  ftil3=ftil2; ftil2=ftil1; ftil1=ftil; frhs3=frhs2; frhs2=frhs1
 CASE(5) 
  ftil4=ftil3; ftil3=ftil2; ftil2=ftil1; ftil1=ftil; frhs4=frhs3; frhs3=frhs2; frhs2=frhs1
 CASE default
  PRINT *, "PrepareMTS_SH_DGV: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
  stop
 END SELECT 
! now, the f contains the solution at the current time. We call the time integrator and then f will contain the next time step and 
! frhs1 will contain the solution at the curent time.
 call TimeIntegratorORK_HOV1D3D
 curr_time = curr_time+dt ! advance the time (the ORK does not advance time)
 
 !!!!!! DEBUG 
 !do iii=1,size(run_mode_1D,1)
  print *, "time=", curr_time, "run_mode_1D=", run_mode_1D
 !end do 
 !!!!!!! END DEBUG 
end do 
END SUBROUTINE PrepareMTS1D3D



end module time_integrator_mod

