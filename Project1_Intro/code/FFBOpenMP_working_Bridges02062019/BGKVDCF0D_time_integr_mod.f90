!
!  time_integr_mod.f90
!
! This module contains routines used for time integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BGKVDCF0D_time_integr_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none



contains 

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TimeIntegratorORK_SH_DGV
! 
! This subroutine implements Optimal Runge Kutta integration in the spatially homogeneous case. 
! 
! It is dependent on the main program: uses variables from the 
! commvar 
!
! the called subroutines look up varaibles in commvar
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

SUBROUTINE TimeIntegratorORK_SH_DGV(f,dt)

USE BGKVDCF0D_commvar, only: rkmts, frhs1
USE BGKVDCF0D_spat_oper_mod

REAL (DP), DIMENSION (:), INTENT(OUT) :: f                  					! arrays to store the right hand side of the solution.
REAL (DP), intent(in) :: dt ! the time step

!!! LOCAL ARRAYS USED IN THE RK STEPS - NOT TO BE CONFUSED WITH THE SOLUTION STORAGE ARRAYS USED ELSEWHERE IN time_integrator_mod AND misc_setup_mod 
REAL (DP), DIMENSION (:), ALLOCATABLE :: ft2, ft3, ft4							! auxiliary variables to keep intermediate time steps

!!! MISC. USEFUL LOCAL VARIABLES AND COUNTERS
INTEGER :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SELECT CASE (rkmts)
	!!! RUNGE-KUTTA FIRST ORDER (EULER METHOD)
	CASE (1)
		CALL SpatialOperSH_DGV (f,frhs1,dt)
		f = f + frhs1
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!! RUNGE KUTTA SECOND ORDER (OPTIMAL 2ND ORDER METHOD)
	CASE (2)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE RK STEPS
      ALLOCATE(ft2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft2)"
	  END IF
	  !!!!! BEGIN THE RUNGE-KUTTA STEP !!!!!
	  !!! t_1:
	  CALL SpatialOperSH_DGV (f, frhs1, dt)
      !!! t_2:
	  CALL SpatialOperSH_DGV (f + frhs1*2/3.0_DP, ft2, dt)
	  !!!!
	   f = f + frhs1/4.0_DP + ft2*3.0_DP/4
	  !!!!
	  deallocate(ft2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! RUNGE KUTTA THIRD ORDER (OPTIMAL 3RD ORDER METHOD)
	CASE (3)
		!!! ALLOCATE THE ARRAYS NEEDED FOR THE RK STEPS
	  ALLOCATE(ft2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft2)"
	  END IF
	  ALLOCATE(ft3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft3)"
	  END IF
	  !!!!! BEGIN THE RUNGE-KUTTA STEP !!!!!
	  !!! t_1:
	  CALL SpatialOperSH_DGV (f, frhs1, dt)
	  !!! t_2:
	  CALL SpatialOperSH_DGV (f + frhs1/2.0_DP, ft2, dt)
	  !!! t_3:
	  CALL SpatialOperSH_DGV (f + ft2*3.0_DP/4, ft3, dt)
	 
	  f = f + frhs1*2/9.0_DP + ft2*3.0_DP/9.0_DP + ft3*4/9.0_DP
	  deallocate (ft2,ft3)
	  !!! RUNGE-KUTTA FOURTH ORDER (RALSTON AND RABINOWITZ)
	CASE (4)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE RK STEPS
   	  ALLOCATE(ft2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft2)"
	  END IF
	  ALLOCATE(ft3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft3)"
	  END IF
      !!!!! BEGIN THE RUNGE-KUTTA STEP !!!!!
	  !!! t_1:
	  CALL SpatialOperSH_DGV (f, frhs1, dt)
      !!! t_2: 
	  CALL SpatialOperSH_DGV (f + frhs1/2.0_DP, ft2, dt)
	  !!! t_3:
	  CALL SpatialOperSH_DGV (f + frhs1/4.0_DP + ft2/4.0_DP, ft3, dt)
	  !!! t_4:
	  CALL SpatialOperSH_DGV (f - ft2 + ft3*2, ft2, dt) 
	  f = f + frhs1/6.0_DP + ft3*2/3.0_DP + ft2/6.0_DP  ! :2nd componenent is re-used for 4th component!
	  deallocate (ft2,ft3)
	!!! RUNGE-KUTTA FIFTH ORDER (RALSTON AND RABINOWITZ - COMPATIBLE WITH THE FOURTH ORDER)
	CASE (5)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE RK STEPS
   	  ALLOCATE(ft2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft2)"
	  END IF
	  ALLOCATE(ft3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft3)"
	  END IF
      ALLOCATE(ft4(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"TimeIntegratorORK_SH_DGV: Allocation error for variable (ft4)"
	  END IF
	  !!!!! BEGIN THE RUNGE-KUTTA STEP !!!!!
	  !!! t_1:
	  CALL SpatialOperSH_DGV (f, frhs1, dt)
	  print *, "1"
	  !!! t_2: 
	  CALL SpatialOperSH_DGV (f + frhs1/2.0_DP, ft2, dt)
	  print *, "2"
	  !!! t_3:
	  CALL SpatialOperSH_DGV (f + frhs1/4.0_DP + ft2/4.0_DP, ft3, dt)
	  print *, "3"
	  !!! t_4: 
	  CALL SpatialOperSH_DGV (f - ft2 + ft3*2.0_DP, ft4, dt)
	  print *, "4"
	  !!! WE NOW INTRODUCE SOME MEMORY SAVINGS AND REPLACE FT3 WITH -FT2/5_DP+FT3*546/625_DP THAT IS GOING TO BE USED IN STEP 6
	  ft3 = -ft2/5.0_DP+ft3*546/625.0_DP
      !!! t_5:
	  !!! ON THIS STEP, WE RE-USE ft2 TO TAKE PLACE OF ft5
	  CALL SpatialOperSH_DGV (f + frhs1*7.0_DP/27.0_DP + ft2*10.0_DP/27.0_DP + ft4/27.0_DP, ft2, dt)
	  print *, "5"
	  !!! t_6: 
	  !!! ON THIS STEP, WE RE-USE ft3 TO TAKE PLACE OF ft6
	  CALL SpatialOperSH_DGV (f + frhs1*28.0_DP/625.0_DP + ft3 + ft4*54/625.0_DP - ft2*378/625.0_DP, ft3, dt)
	  print *, "6"
	  f = f + frhs1/24.0_DP + ft4*5.0_DP/48.0_DP + ft2*27.0_DP/56.0_DP + ft3*125/336.0_DP
	  deallocate (ft2,ft3,ft4)
	CASE default
			PRINT *, "Unsupported degree of Runge-Kutta method passed"
			stop

END SELECT 
	   frhs1 = frhs1/dt
END SUBROUTINE TimeIntegratorORK_SH_DGV

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

SUBROUTINE TimeIntegratorMTS_SH_DGV(f,dt)

USE BGKVDCF0D_commvar, only: f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4, &
                    rkmts,frhs5,time
USE BGKVDCF0D_spat_oper_mod                    
                    
REAL (DP), dimension (:), intent (out) :: f ! the solution at the current time step                    
REAL (DP), intent (in) :: dt ! the time step. 

!!!!!!!!

SELECT CASE (rkmts)
	CASE (1)
	 !!! MTS FIRST ORDER (EULER METHOD)
	 CALL SpatialOperSH_DGV (f, frhs1, dt)
	 !
	 f = f + frhs1
     !!! END MTS 1ST ORDER (EULER METHOD)
	CASE (2)
	 !!! ADAMS-BASHFORTH TWO STEP EXPLICIT METHOD
	 f1 = f; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperSH_DGV (f, frhs1, dt)
	 f = f + frhs1*3.0_DP/2.0_DP - dt*frhs2/2.0_DP
	 !!! END ADAMS-BASHFORTH TWO STEP EXPLICIT METHOD
	CASE (3)
	 !!! BEGIN ADAMS-BASHFORTH THREE STEP
	 f2 = f1; f1 = f; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperSH_DGV (f, frhs1, dt)
	 f = f + frhs1*23.0_DP/12.0_DP - dt*frhs2*4/3.0_DP + dt*frhs3*5.0_DP/12.0_DP 
    !!! END ADAMS-BASHFORTH THREE STEP
	CASE (4)
   	 !!! ADAMS-BASHFORTH FOUR STEP
	 f3 = f2; f2 = f1; f1 = f; frhs4 = frhs3; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperSH_DGV (f, frhs1, dt)
	 f = f + frhs1*55.0_DP/24.0_DP - dt*frhs2*59.0_DP/24.0_DP + dt*frhs3*37.0_DP/24.0_DP - dt*frhs4*3.0_DP/8
	 !!! END ADAMS-BASHFORTH FOUR STEP
	CASE (5)
	 !!! ADAMS-BASHFORTH FIVE STEP METHOD
	 f4 = f3; f3 = f2; f2 = f1; f1 = f 
	 frhs5 = frhs4; frhs4 = frhs3; frhs3 = frhs2; frhs2 = frhs1 !! first we need to push the solution into the storage arrays
	 CALL SpatialOperSH_DGV (f, frhs1, dt)
	 f = f + frhs1*(1901.0_DP/720.0_DP) - dt*frhs2*(2774.0_DP/720.0_DP) + dt*frhs3*(2616.0_DP/720.0_DP) &
	       - dt*frhs4*(1274.0_DP/720.0_DP) + dt*frhs5*(251.0_DP/720.0_DP)
	 !!! END ADAMS-BASHFORTH FIVE STEP METHOD
    CASE default
	 PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT
    frhs1 = frhs1/dt
    time = time + dt ! advance the time for one step
END SUBROUTINE TimeIntegratorMTS_SH_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine PrepareMTS_SH_DGV
!
! Subroutine prepares the intermediate time steps to be used in the multi-step methods.
! f stores the actual solution.
! f1, f2, f3, f4, frhs1, frhs2, frhs3, frhs4 store the value of the 
!	spatial operator and the right hand side on those intermediate time steps 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PrepareMTS_SH_DGV (f, dt)

USE BGKVDCF0D_commvar, only:  f1, f2, f3, f4, frhs1, frhs2, frhs3, frhs4, &
							frhs5, rkmts, time
USE BGKVDCF0D_spat_oper_mod

REAL (DP), dimension (:) :: f ! the solution on the initial time step
REAL (DP), intent (in) :: dt ! the time step..

!!!!!! DECLARE LOCAL VARIABLES !!!!!
INTEGER :: i													! Local counters
INTEGER :: loc_alloc_stat ! to keep the allocation status 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to allocate the storages
SELECT CASE(rkmts)
	CASE(1)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1)"
		stop
	  END IF
	CASE(2)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),frhs2(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	CASE(3)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	ALLOCATE(frhs2(1:size(f,1)),frhs3(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
		stop
	  END IF
	CASE(4) 
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
		stop
	  END IF
	  ALLOCATE(frhs3(1:size(f,1)),frhs4(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs3,f3)"
		stop
	  END IF
    CASE(5)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
		stop
	  END IF
	  ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
		stop
	  END IF
	ALLOCATE(frhs3(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs3,f3)"
		stop
	  END IF
	  ALLOCATE(frhs4(1:size(f,1)),frhs5(1:size(f,1)),f4(1:size(f,1)), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs4,f4)"
		stop
	  END IF
	CASE default
			PRINT *, "PrepareMTS_SH_DGV: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
			stop
END SELECT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now we assume that f contains the intial data. We need to set up the arrays
! If rkmts=1 - no action is needed all data is in place 
! If the rkmts>2 then we need to actually prepare some arrays. 
do i=2,rkmts
 !! FIRST WE NEED TO PUSH THE SOLUTION AND THE RIGHT SIDE DOWN THE STORAGE TO OPEN SOME SPACE FOR THE NEW DATA 
 SELECT CASE(rkmts)
 CASE(2)
  f1=f 
 CASE(3)
  f2=f1; f1=f; frhs2=frhs1
 CASE(4)
  f3=f2; f2=f1; f1=f; frhs3=frhs2; frhs2=frhs1
 CASE(5) 
  f4=f3; f3=f2; f2=f1; f1=f; frhs4=frhs3; frhs3=frhs2; frhs2=frhs1
 CASE default
  PRINT *, "PrepareMTS_SH_DGV: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
  stop
 END SELECT 
! now, the f contains the solution at the current time. We call the time integrator and then f will contain the next time step and 
! frhs1 will contain the solution at the curent time.
 call TimeIntegratorORK_SH_DGV(f,dt)
 time = time+dt ! advance the time
end do 
END SUBROUTINE PrepareMTS_SH_DGV

end module BGKVDCF0D_time_integr_mod
