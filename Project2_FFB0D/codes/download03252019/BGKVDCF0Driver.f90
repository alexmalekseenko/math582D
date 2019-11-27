!
!  Collection of subroutines for 3D nodal-DG velocity discretization of the Boltzmann collision integral 
!  Will be integrated with CLAWPACK and with other codes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !Needed for the Intel MKL FFT Interface
      program main
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Envoke modules with type definitions
      use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Envoke Driver Modules
      use BGKVDCF0D_readwrite ! contains useful subroutines for reading/writing. in particular the reading of the problem data file
      use BGKVDCF0D_miscset ! contains various setup subroutines, defines useful arrays, variables etc.
      use BGKVDCF0D_time_integr_mod !
     
      ! Establish access to the common variables 
      use BGKVDCF0D_commvar, only: k_b,k_c,d_c,d_b,N,f,time,dt,t_L,t_R,num_save_solution,&
				num_eval_error,rkmts,need_to_restart,restart_time_txt
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! BEGIN DECLARATION OF DGV LIBRARIES
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Envoke modules of the DGV library
      use DGV_dgvtools_mod
      use DGV_collision_mod
      use DGV_sf02
      use DGV_miscset
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Establish access to the variables in the DGV library
      use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,nodes_v,nodes_w,A_capphi,&
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,&
                  min_sens,Trad,fm,fmA,run_mode,Cco,pad, &
				  Order

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      ! END DECLARATION OF DGV library
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use MKL_DFTI !Needed for the Intel MKL FFT Interface
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Forbid using implicit declarations of variables
      implicit none

      
     ! Variables
 
 real (DP) :: test1 ! scrap variables
 real (DP), dimension(:), allocatable :: fcol ! scrap variable to keep results of the evaluation of collision operator 
 integer (I4B) :: gg,loc_alloc_stat ! to keep them allocation status
 !!!!!!!!!!!
 real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
 !!!!!!!!!!! SAVING SOLUTION AND SUCH !!!!!!!!!!!!!!!
 integer (I4B) :: rec, time_step, err_count
 real (DP) :: next_time_file_record, file_record_period, next_time_error_eval, error_eval_period
 character (len=15) :: suff  ! the string to hold the suffix
 real (DP)  :: ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom4_u,mom4_v,mom4_w,mom6_u,mom6_v,mom6_w ! number density, av_v
 real (DP)  :: mom3_u,mom3_v,mom3_w,mom5_u,mom5_v,mom5_w ! number density, av_v
 real (DP), dimension (1:1000) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a, & 
            mom3_u_a,mom3_v_a,mom3_w_a,mom4_u_a,mom4_v_a,mom4_w_a,mom5_u_a,mom5_v_a,mom5_w_a, & 
            mom6_u_a,mom6_v_a,mom6_w_a  ! scrap arrays to keep the macropameters in.
 !!!!!!!!!!!
Integer :: Status
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialization: we read the parametes to set variables in BGKVDCF0D_commvar
      ! Read the problem parameters from the data file and set the problem variables. 
      call SetUWbgkParams("parameters.dat")
      ! set the order of used gauss-lobatto nodes in x and t
      ! these are not used in the OD driver
      k_b=7 !
      d_b=8
      d_c=2
      k_c=2
      !!!!!!!!!!!
      call SetDGblzmGnodes ! sets some basic arrays of gauss and gauss-lobatto nodes. 
      call SetBGK1Dmesh ! set the meshes in the x and time much of it is not useds in the 0D case
      ! set the dimensions of the level zero meshes... 
      N=3
      !
      ! if (loc_alloc_stat >0) then 
     !print *, "DGVblzm: Allocation error for variables (f), (fcol)"
     !end if 
!!!!!!!! try to check the calculation time for initialization of DGV library 
!   call cpu_time (pr_time_1)
!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of the paramters in the DGV library: 
!!!
      call InitDGV0D
      pad=0 ! This is a padding with zeros parameter -- #pad zeros will be added to the arrays on Both ends  
      call CalcFftZeroPadA(pad)
      !call CalcFftA
!!! 
! End Initialization of the DGV Library:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   call cpu_time (pr_time_2)
!   print *, "Processor time lapsed in seconds on thread 0:", pr_time_2 - pr_time_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
     allocate (f(1:size(nodes_u,1)), fcol(1:size(nodes_u,1)), stat=loc_alloc_stat)
!     allocate (fmA(1:size(nodes_u,1),1:size(nodes_u,1)), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then 
     print *, "DGVblzm: Allocation error for variables (f), (fcol)"
     end if 
     
  !!!!!!!!!!!!!!!!!!!!! TIME EVOLUTION !!!!
     rkmts = 5
	 if (need_to_restart) then
      call RestartSolutionBGK0D(restart_time_txt)
       Read (restart_time_txt, fmt = * ,  iostat = loc_alloc_stat) time ! set the current time to be equal to the time of restore  
	   err_count=0
      else 
       time = t_L !initial time 
       f = f_1D3D(0.0_DP, nodes_u, nodes_v, nodes_w, time) ! set up initial data
       err_count=0
       suff = "i"
         
       
      ! call WriteSol_BGK0D (suff)
     ! compute the macroparameters of the initial data and save it on disk
     !!!!!!!!!!!!!!!!!!!!
     ! Call of a function from the DGV library: 
     call MassCheckRecPlusPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u,& 
           mom3_v,mom3_w,mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
     err_count=err_count+1
     if (err_count<1000) then
        time_a(err_count) = time
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
		mom3_u_a(err_count) = mom3_u
		mom3_v_a(err_count) = mom3_v
		mom3_w_a(err_count) = mom3_w
		mom4_u_a(err_count) = mom4_u
		mom4_v_a(err_count) = mom4_v
		mom4_w_a(err_count) = mom4_w
		mom5_u_a(err_count) = mom5_u
		mom5_v_a(err_count) = mom5_v
		mom5_w_a(err_count) = mom5_w
		mom6_u_a(err_count) = mom6_u
		mom6_v_a(err_count) = mom6_v
		mom6_w_a(err_count) = mom6_w
     else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
     end if
     !!!!!save macroparameters on hard drive
     !call WriteErrPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
     !                         vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
     !                         tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
     call WriteErrPlusPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
							  mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),&
							  mom4_u_a(1:err_count),mom4_v_a(1:err_count),mom4_w_a(1:err_count),&
							  mom5_u_a(1:err_count),mom5_v_a(1:err_count),mom5_w_a(1:err_count),&
							  mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
     !!!!! The initial data  and the initial macroparameters were saved on the hard drive    !!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call cpu_time (pr_time_1) ! Start the account of time
     ! prepare MTS steps.
     time_step=0
     call PrepareMTS_SH_DGV (f, dt)
     call cpu_time (pr_time_2) ! End the account of time
     print *, "Processor time lapsed in seconds for Prepare MTS:", pr_time_2 - pr_time_1 ! Report the account of time
     ! now all the arrays are ready for the time evolution.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
     !!!!!save the solution at the initial stage:                                            !!!!!
       !rec=0
       !write (suff, "(I4)") rec
       write (suff, "(F14.10)") time 
       call WriteSol_BGK0D (suff) ! solution is saved
     !!!!!evaluate macroparameters after the MTS prepare stage
	 call MassCheckRecPlusPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u,&
	               mom3_v,mom3_w,mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
     err_count=err_count+1
     if (err_count<1000) then
        time_a(err_count) = time
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
		mom3_u_a(err_count) = mom3_u
		mom3_v_a(err_count) = mom3_v
		mom3_w_a(err_count) = mom3_w
		mom4_u_a(err_count) = mom4_u
		mom4_v_a(err_count) = mom4_v
		mom4_w_a(err_count) = mom4_w
		mom5_u_a(err_count) = mom5_u
		mom5_v_a(err_count) = mom5_v
		mom5_w_a(err_count) = mom5_w
		mom6_u_a(err_count) = mom6_u
		mom6_v_a(err_count) = mom6_v
		mom6_w_a(err_count) = mom6_w
     else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
     end if
     !!!!!!!!!!!!!! Diagnostics for prakash
     call KeepTrackCco_BGK0D
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!save macroparameters on hard drive
     !call WriteErrPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
     !                         vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
     !                         tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
     call WriteErrPlusPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
							  mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),&
							  mom4_u_a(1:err_count),mom4_v_a(1:err_count),mom4_w_a(1:err_count),&
							  mom5_u_a(1:err_count),mom5_v_a(1:err_count),mom5_w_a(1:err_count),&
							  mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
     !!!!! The initial data  and the initial macroparameters were saved on the hard drive    !!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     end if

      ! measure time:
     call cpu_time (pr_time_1) ! Start the account of time
     ! set up error evaluation and solution recording
     file_record_period = (t_R-t_L)/Real(num_save_solution,DP)
     next_time_file_record = time + file_record_period
     error_eval_period = (t_R-t_L)/Real(num_eval_error,DP)
     next_time_error_eval = time + error_eval_period


!!!!!!! THE TIME LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do while (time <= t_R+0.00000000001_dp)
      time_step=time_step+1
	  
      call TimeIntegratorMTS_SH_DGV(f,dt)
      ! check if it is time to evaluate error or any other quantity (like mass, bulk vel and tempr)
      if (time >= next_time_error_eval) then
       next_time_error_eval = time + error_eval_period
       ! evaluate mass momenum and tempretaure of the solution
       ! end eevaluate mass momenum and tempretaure of the solution
       ! call MassCheckRecPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w) ! put another subroutine for higher moments
	   call MassCheckRecPlusPlus (f,ndens,ubar,vbar,wbar,tempr,tempr_u,tempr_v,tempr_w,mom3_u, &
	         mom3_v,mom3_w,mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w)
       err_count=err_count+1
       if (err_count<1000) then
        time_a(err_count) = time
        ndens_a(err_count) = ndens
        ubar_a(err_count) = ubar
        vbar_a(err_count) = vbar
        wbar_a(err_count) = wbar
        tempr_a(err_count) = tempr
        tempr_u_a(err_count) = tempr_u
        tempr_v_a(err_count) = tempr_v
        tempr_w_a(err_count) = tempr_w
		mom3_u_a(err_count) = mom3_u
		mom3_v_a(err_count) = mom3_v
		mom3_w_a(err_count) = mom3_w
		mom4_u_a(err_count) = mom4_u
		mom4_v_a(err_count) = mom4_v
		mom4_w_a(err_count) = mom4_w
		mom5_u_a(err_count) = mom5_u
		mom5_v_a(err_count) = mom5_v
		mom5_w_a(err_count) = mom5_w
		mom6_u_a(err_count) = mom6_u
		mom6_v_a(err_count) = mom6_v
		mom6_w_a(err_count) = mom6_w
       else
        print *, "the number of error evaluations exceeded the storage, err_count >= 1000"
       end if  
      end if 
      ! Check if it is time to save the solution 
      if (time >= next_time_file_record) then 
		   next_time_file_record = next_time_file_record + file_record_period
		   ! write the solution of disk
		   ! rec=rec+1
		   ! write (suff, "(I4)") rec
		   write (suff, "(F14.10)") time 
		   call WriteSol_BGK0D (suff)
	!       call WriteErrPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
	!                              vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
	!                              tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count))
		   call WriteErrPlusPlus_BGK0D (time_a(1:err_count),ndens_a(1:err_count),ubar_a(1:err_count),&
								  vbar_a(1:err_count),wbar_a(1:err_count),tempr_a(1:err_count),&
								  tempr_u_a(1:err_count),tempr_v_a(1:err_count),tempr_w_a(1:err_count),&
								  mom3_u_a(1:err_count),mom3_v_a(1:err_count),mom3_w_a(1:err_count),&
								  mom4_u_a(1:err_count),mom4_v_a(1:err_count),mom4_w_a(1:err_count),&
								  mom5_u_a(1:err_count),mom5_v_a(1:err_count),mom5_w_a(1:err_count),&
								  mom6_u_a(1:err_count),mom6_v_a(1:err_count),mom6_w_a(1:err_count))
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! will try to track the deviation from the Maxwellian 
      ! For that purpose we will create a global allocatable array that will store the times and the values of L1_err 
      ! It will periodically dump it to disk... all on its own... 
      !!!!!!!!!!!!
      call KeepTrackL1_err_BGK0D(f)
      call KeepTrackCco_BGK0D
      !!!!!!!!!!!!
      
	  print *, "time ", time
     end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     call cpu_time (pr_time_2) ! End the account of time
     print *, "Processor time lapsed in seconds for main run:", pr_time_2 - pr_time_1 ! Report the account of time
! we now test evaluation of the basis functions: 
!    fcol=IphiDGV(f) 
! EVALUATE THE COLLISION INTEGRAL ON A PERIODIC GRID ....
! First we need to prepare some helpuf arrays: 
  !call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP))
! NOW WE call the subrouine that will use periodicity of the grid to evaluate the collision integral      
!!!!!!!!!!! try to check the calculation time
!   call cpu_time (pr_time_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!   call EvalCollisionPeriodicA_DGV(f,fcol)
! end evaluation of the collision integral      !
!   call cpu_time (pr_time_2)
!   print *, "Processor time lapsed in seconds on thread 0:", pr_time_2 - pr_time_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    test1 = MassCheck(fcol) 
     end program main
