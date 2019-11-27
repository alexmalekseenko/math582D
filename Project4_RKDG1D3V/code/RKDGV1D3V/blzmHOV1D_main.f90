! 05/26/2017 Alex
!
!  blzmHOV1D_main
!
!  This file is to put together 1D HOV model
!%%%%%%%%%%%%%%%%
program main
 use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
 use common_variables_mod, only: max_deg,k,N,moments_refine_x,dt,curr_time,moments_x_gauss_order, &          ! contains the main variables and some useful constants 
			k_list,N_list,rk,dt,cfl,restart_time_txt,need_to_restart,initial_time,final_time,&
			num_save_solution, num_eval_error,xmesh,s, M,u_left, u_right, mesh_u_uniform
						
 use read_write_tools_mod ! contains subr. to read from parameters file, 
 use misc_setup_mod       ! contains subrs. to initialize and modify meshes.
 use time_integrator_mod  ! contains Runge Kutta integrators 
 use sol_integr_mod       ! contains integrators and set up for the moments
 use mpi_routines
 
      !use basis_fun_mod
      !use source_f_mod ! the source function sf1D1D is located in sf1D1D_xx_mod
      !use poly_tool_mod
      !use spectral_tools_mod
      !use gaussian_mod
      !use eigenfs_mod
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BEGIN DECLARATION OF DGV LIBRARIES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Envoke modules of the DGV library
      use DGV_miscset
      use DGV_mpiroutines
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Establish access to the variables in the DGV library
      use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,nodes_v,nodes_w,&
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,run_mode,&
                  num_lin_proc,&
                  u_L,u_R,mesh_u_uniformDGV => mesh_u_uniform
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! END DECLARATION OF DGV library
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    
                 

 !
 ! intrinsic functions
 !
intrinsic Real, Max, Min, MinVal, cpu_time
 ! 
 ! variables 
 ! 
 integer (I4B) :: i,time_step ! local counter  
 real (DP) :: file_record_period, next_time_file_record, error_eval_period, next_time_error_eval ! useful variables to manage file records
 integer (I4B) :: error_eval_count ! counter to count evaluations of the error
 integer :: num_cells_index ! to use in selection of number of nodes
 integer :: loc_alloc_stat! to keep allocation status
 character (len=20) :: parchar    ! string to keep char value of the parameter 
 real :: pr_time_1, pr_time_2 ! variable to calculate the processor time 
 !!!!!!!!!!! MPI Stuff goes below... 
 integer :: ierr,irank ! variables for MPI Calls
 !!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! THIS IS AN MPI VERSION OF THE CODE. 
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! INITIALIZE THE MPI environment
 call mpi_init(ierr)
 !!! Now the MPI routines and variables will make sense ..
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on... 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
 ! SET UP PARAMETERS AND INITIALIZE VARIABLES    
 !!!!!!!!!!!!!!!!!!!!
 ! Some variables are obtained from an input file 
 !!! MPI FORK 
 if (irank == 0) then 
     call Set1DHOVParameters("parameters.dat",0) ! Output only produced on the mode with irank=0
 else 
     call Set1DHOVParameters("parameters.dat",1) ! this one does the same but the printout is commented 
 end if 
 !!! END MPI FORK
 !!!! Now the variables are set... 
 ! The obtained varaibles are used to set up other varaibles:
k=k_list(1)            ! order of polynomial approximation in (x)
! absolete - replaced by dgvlib  s=s_list(1)            ! order of polynomial approximation in (u)
max_deg = k ! in this setting, max degree in x and max degree in u,w, are used 
 ! set up the meshes in x. Meshes in u are set up in the DGVlib initialization step: 
 ! first the mesh in x
N=N_list(1)            ! the number of (uniform) mesh cells in x 
! absolete - delete - M=M_list(1)            ! the number of (uniform) mesh cells in u 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!! The rest of setup is done using N:
! Call a piece of code that sets up gaussian nodes for the integration of moments and arrays to store moments
call Set1DHOVMomentsTimeArraysGaussNodesWeights  
! Call a  piece of code that sets up the meshes in variable x initially (may use gauss nodes for the integration of moments) 
call Set1DHOVspatMesh     
 !  Now we set up some useful matrices: 
call Set1DHOVSqp  ! calculate Sqp,j = \int_{-1}^{1} \varphi_{p,j}(x)\partial_{x}\varphi_{q,j}(x) d x
 !!!  Later add here the matrix of values of Legendre's polinomials to be used in integration of the right side...
call SetQ2plus1Altones ! calculate useful vector of coefficinets: (1,3,...,2q+1,...,2k+1), (1,-1,1,-1,...,(-1)^k)
!! this one is absolete -- delete -- these matrices are not needed in the nodal-DG velocity formulation   
! call Set1DHOValphamiAmlAinvml ! calculate \alpha_{m,i}A_{ml,i},A^{ml,i} -- the eigenvectors and eigenvalues 
							  ! of \hat{T}_{pl,i}=\int_{U_{i}}u\lambda_{l,i}(u)\lambda_{p,i}(u)
!!!!!!!!!!!!!!!!!!!!!!!!!!! end - delete !!!!!!!!!!!!!!!!!!!!!

call Set1DHOVLegGaussArraysSol
call SetMomentsXnodesBasisFunHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of the paramters in the DGV library: 
!!!
      call InitDGV1D_MPI(N*(moments_x_gauss_order)*moments_refine_x,irank) ! the evaluation of the collision operator will be 
                              ! conducted in gauss nodes of the spatial cells. We will use moments_x_gauss_order*moments_refine_x nodes in each cell
                              ! so, we will need to store the values for N*moments_x_gauss_order*moments_refine_x points in space in total 
                              
                              ! after this subroutine is called, all variables of the primary and secondary velocity meshes are 
                              ! establised.
!!! 
! End Initialization of the DGV Library:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we initialize the variables of the driver with some relevant values of the dgv library.
! thes evariables of the DGV library are used in names of the files and in some other subroutines, but they are
! not used in the development of the variables.
 s=su
 M=Mu
 u_left=u_L
 u_right=u_R
 mesh_u_uniform = mesh_u_uniformDGV  ! should be "yes" 
!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI FORK 
!!
!!
!! only processors with irank < num_lin_proc will participate in solution of the transport part. 
!! All other processes will be trapped into a loop to help with other tasks or to do nothing.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (irank<num_lin_proc) then 
 !!!!!  
 if (need_to_restart) then
  call RestartSolutionMTS1D3D_MPI(restart_time_txt,irank) ! processors use MPI IO to read the solution
  Read (restart_time_txt, fmt = * ,  iostat = loc_alloc_stat) curr_time ! set the current time to be equal to the time of restore  
 else 
  !!!!! Measure CPU time first take 
  call cpu_time (pr_time_1)
  !!!!!!!!!!!!!  
  curr_time = initial_time       ! initial time  
  call Set1DHOVinitftil_MPI  ! calculates initial data for ftil1 and ftil2: \tilde{f}^{1,2}_{p,j;m,i}(t) 
  call PrepareMTS1D3D_MPI(irank)    ! create additional solution variables and calls RK integrator rk-1 times. After this MTS integration can be called
  ! Record initial solution (record the spectral coefficients 
  write (parchar, "(F11.8)") curr_time
  call WriteSolutionMTS1D3D_MPI(parchar,irank)  ! processors use MPI IO to write the solution
  !!!!!!!! measure CPU time second take and print
  call cpu_time (pr_time_2)
  if (irank == 0) then
   print *, "Processor time lapsed in seconds to run PrepareMTS1D3D_MPI", pr_time_2 - pr_time_1
  end if 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 end if
 ! allocate space for error recording and also gaussian nodes and weights for error evaluation 
 !!! ALEX: need to deside on the arry to store the computed macroparaemters -- may be created in this subroutine! 
 !! call SetErrTimeEvalArraysHOV1D
 !!!!!
 error_eval_count = 1
 !!!!! Calculate the moments of the solution at the initial time 
 !!! STILL NEED TO WRITE THAT SUBROUTINE !!! call Int20MomentsSolHOV1D3D(curr_time, error_eval_count) !  use umesh for HOV1D evaluation 
 !!! STILL NEED TO WRITE THIS SUBROUTINE !!! call WriteMomentsSolHOV1D3D(error_eval_count) 

 !!!!!!!!!!! measure the CPU TIME -- intake initial time: !!!!!
 call cpu_time (pr_time_1)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Main Loop
!
! This block loops the time integration, saves the intermediate solution the given number of times, 
! evaluates the error if exact solution is available and so on.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up error evaluation and solution recording 
file_record_period = (final_time-initial_time)/Real(num_save_solution,DP)
next_time_file_record = curr_time+file_record_period
error_eval_period = (final_time-initial_time)/Real(num_eval_error,DP)
next_time_error_eval = curr_time+error_eval_period
! THE TIME LOOP
time_step = 0
do while (curr_time <= final_time)
   ! time_step=time_step+1 
   ! write (*,*) " time step: ", time_step
   call TimeIntegratorMTS1D3D_MPI
   !!! Check if it is time to calculate the error
   if (curr_time >= next_time_error_eval) then 
    next_time_error_eval = next_time_error_eval + error_eval_period
    error_eval_count = error_eval_count+1
    !!!!! Calculate the moments of the solution at the initial time 
    !!! STILL NEED TO WRITE THAT SUBROUTINE !!! call Int20MomentsSolHOV1D3D(curr_time, error_eval_count) !  use umesh for HOV1D evaluation 
   end if     
   !!! Check if it is time to save the solution 
   if (curr_time >= next_time_file_record) then 
    next_time_file_record = next_time_file_record + file_record_period
    ! write the solution of disk
    write (parchar, "(F11.8)") curr_time
    call WriteSolutionMTS1D3D_MPI(parchar,irank)  ! processors use MPI IO to write the solution
    !!! STILL NEED TO WRITE THIS SUBROUTINE !!! call WriteMomentsSolHOV1D3D(error_eval_count)
   end if 
end do 
! Write one last time
write (parchar, "(F11.8)") curr_time
call WriteSolutionMTS1D3D_MPI(parchar,irank)  ! processors use MPI IO to write the solution
!!
!! This concludes the work on the process with irank < num_lin_procs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else 
!! All processors with irank >= num_lin_proc will land here. 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call FlyTrap1Donecell_DGV_MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This concludes the work on the process with irank > num_lin_procs
end if 

call mpi_finalize(ierr)

if (irank==0) then 
 !!!!!!!!!!! measure the CPU TIME -- intake final time: !!!!!
 call cpu_time (pr_time_2)
 print *, "experiment: k=", k, " max_deg=", max_deg, " rk=", rk, " N=", N
 print *, "Processor time lapsed in seconds", pr_time_2 - pr_time_1
 print *, " |----------------------------------------------------------------| "
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 



end program main
