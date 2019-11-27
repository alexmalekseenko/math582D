!
!  read_write_tools_mod.f90
!
!  7/30/2008 3:58:49 PM
!
! This modle contains useful procedures for reading and writing from the disk
!
!   
! INCLUDES: 
!  Set1DHOVParameters('parameters.dat')
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
module read_write_tools_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

contains 

subroutine Set1DHOVParameters(pfname,slt)
use common_variables_mod
intrinsic Real,Index, Len, Scan, Trim

character (len=*), intent (in) :: pfname ! Name of the file where the parameters are stored 
integer, intent (in) :: slt ! parameter to determine if the output is produced. slt=0 -- output is produced, otherwise no 
character (len=132) :: line              ! string to keep one record 
character (len=10) :: fmtchar            ! string to keep format 
character (len=40) :: line_head          ! string to keep the line header 
integer :: code_file, code_line          ! get the return code from READ
integer :: m_count                            ! dump counter 
integer :: pos                 ! to store position within the string 
integer (I4B), dimension (20) :: i_bulk  ! to store temporarily integers that has been read from the parameter file 
real (DP), dimension (20) :: r_bulk      ! to store temporarily reals that has been  read from the parameter file 
integer :: loc_alloc_stat ! some dump varaible
!
open (15, file=pfname, position ="REWIND", action="READ") ! open the file for reading
code_file=0
do while (code_file == 0)
 read (15, "(A)", iostat=code_file)  line                               ! read one line
 pos = Scan(line, "=")
  if ((pos /= 0) .and. (line(1:1)/="!")) then
  write (fmtchar, "(I3)") pos-1 
  read (line, "(A"//trim(fmtchar)//")", iostat = code_line ) line_head  
  line_head = trim(line_head)
  line(:) = line(pos+1:)   ! remove the heading from the line 
  !
  pos = Scan(line, "!")
  if (pos > 1) then 
  line(:) = trim(line(:pos-1))   ! remove any comments from the line
  end if 
  ! 
  select case (line_head)
   case ("scheme order in x")           ! ready to set up order in variable x
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (k_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (k_list)"
     stop
     end if
     !
     k_list=i_bulk(1:m_count)
    end if  
    if (slt==0) then 
     print *, "k_list=", k_list 
    end if 
   case ("scheme order in u")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (s_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (s_list)"
     stop
     end if
     !
     s_list=i_bulk(1:m_count)
    end if  
    if (slt==0) then 
     print *, "s_list=", s_list 
    end if  
   case ("left endpoint in x")         ! ready to set up left endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_left=r_bulk(1) ! all other inputs values are ignored 
    else 
     x_left=0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "x_left=", x_left 
    end if  
   case ("right endpoint in x")         ! ready to set up right endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_right=r_bulk(1) ! all other input values are ignored 
    else 
     x_right=1_DP ! set up to the default value 
    end if  
    if (slt==0) then
     print *, "x_right=", x_right
    end if  
   case ("left endpoint in u")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_left=r_bulk(1) ! all other inputs values are ignored 
    else 
     u_left=0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "u_left=", u_left
    end if  
   case ("right endpoint in u")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_right=r_bulk(1) ! all other input values are ignored 
    else 
     u_right=1_DP ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "u_right=", u_right
    end if  
   case ("uniform mesh in x")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
   if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_x_uniform = .TRUE.
    else
     mesh_x_uniform = .FALSE.
    end if 
    if (slt==0) then
     print *, "mesh_x_uniform=", mesh_x_uniform
    end if 
   case ("uniform mesh in u")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_u_uniform = .TRUE.
    else
     mesh_u_uniform = .FALSE.
    end if
    if (slt==0) then 
     print *, "mesh_u_uniform=", mesh_u_uniform 
    end if  
   case ("number of cells in x")           ! ready to set up the number of cells in x
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (N_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (N_list)"
     stop
     end if
     !
     N_list=i_bulk(1:m_count)
    end if  
    if (slt==0) then 
     print *, "N_list=", N_list
    end if  
   case ("number of cells in u")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (M_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (M_list)"
     stop
     end if
     !
     M_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then 
     print *, "M_list=", M_list
    end if  
   case ("order of Runge Kutta method")         ! ready to set up selected order of the Runge Kutta time integrator
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    rk=i_bulk(1) ! all other inputs values are ignored 
    else 
    rk=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "rk=", rk   
    end if  
   case ("conditions on the left boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_leftbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_leftbcond=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "selected_leftbcond=", selected_leftbcond  
    end if  
   case ("conditions on the right boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_rightbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_rightbcond=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "selected_rightbcond=", selected_rightbcond  
    end if  
   case ("type of exact solution")         ! ready to set up selected type of boundary data
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_exact_sol=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_exact_sol=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "selected_exact_sol=", selected_exact_sol  
    end if  
   case ("current solution save directory")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_dir = "solution/"
    else
     current_solution_dir = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
      print *, "current_solution_dir=", current_solution_dir 
    end if   
   case ("current solution base name")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_base_name = "bgk1d"
    else
     current_solution_base_name = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
     print *, "current_solution_base_name=", current_solution_base_name 
    end if  
   case ("initial time")       ! ready to set up the initial time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     initial_time=r_bulk(1) ! all other input values are ignored 
    else 
     initial_time=0_DP ! set up to the default value 
    end if
    if (slt==0) then 
     print *, "initial_time=", initial_time
    end if  
   case ("final time")         ! ready to set up the end time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     final_time=r_bulk(1) ! all other input values are ignored 
    else 
     final_time=1_DP ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "final_time=", final_time
    end if  
   case ("instances to evaluate error")  ! ready to set up the number of error evaluations
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1000) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_eval_error=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_eval_error=1000 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "num_eval_error=", num_eval_error 
    end if  
   case ("instances to save solution")  ! ready to set up the number solution recordings
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 10) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_save_solution=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_save_solution=10 ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "num_save_solution=", num_save_solution 
    end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! added Alex 10/15/08
   case ("type of nonuniform mesh in x")  ! ready to set up the type of nonuniform mesh in x
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    x_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    x_nonuniform_mesh_type=1 ! set up to the default value 
    end if
    if (slt==0) then 
     print *, "x_nonuniform_mesh_type=", x_nonuniform_mesh_type 
    end if  
   case ("type of nonuniform mesh in u")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    u_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    u_nonuniform_mesh_type=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "u_nonuniform_mesh_type=", u_nonuniform_mesh_type
    end if  
   case ("gauss order for moments in x")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_x_gauss_order=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_x_gauss_order=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "moments_x_gauss_order=", moments_x_gauss_order 
    end if  
   case ("mesh refinement in x for moments")  ! ready to set up the coefficient of mesh refinement in x for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_refine_x=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_refine_x=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "moments_refine_x=", moments_refine_x
    end if 
   case ("gauss order for moments in u")  ! ready to set up the order of gauss method in u variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_u_gauss_order=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_u_gauss_order=1 ! set up to the default value 
    end if  
    if (slt==0) then
     print *, "moments_u_gauss_order=", moments_u_gauss_order 
    end if  
   case ("mesh refinement in u for moments")  ! ready to set up the coefficient of mesh refinement in u for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_refine_u=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_refine_u=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "moments_refine_u=", moments_refine_u
    end if  
!!! end added Alex 10/15/08
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! added Alex 12/26/08
    ! we are ready to read the ordinary gas constant, reference teperature and viscosity and alpha constant 
    ! for the experiment --- will need in the collision term
  case ("the ordinary gas constant")         ! ready to set up the ordinary gas constant
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasR=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasR=1.0_DP ! set up to the default value (garbage! gasR=1 has no meaning!!)
    end if 
    if (slt==0) then 
      print *, "gasR=", gasR
    end if  
    ! we are ready to read the gas reference temperature for the experiment --- will need in the collision term
  case ("gas reference temperature")         ! ready to set up gas reference temperature
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasTref=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasTref=1.0_DP ! set up to the default value (garbage! gasTref=1 has no meaning!!)
    end if
    if (slt==0) then 
     print *, "gasTref=", gasTref 
    end if  
    ! we are ready to read the gas reference viscosity for the experiment --- will need in the collision term
  case ("gas reference viscosity")         ! ready to set up gas reference viscosity
        !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasmuref=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasmuref=1.0_DP ! set up to the default value (garbage! gasmuref=1 has no meaning!!)
    end if 
    if (slt==0) then 
     print *, "gasmuref=", gasmuref 
    end if  
    ! we are ready to read the gas "alpha" in the gas state law for the experiment --- will need in the collision term
  case ("gas alpha constant")         ! ready to set up gas alpha constant
        !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 1.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     gasalpha=r_bulk(1) ! all other inputs values are ignored 
    else 
     gasalpha=1.0_DP ! set up to the default value (garbage! gasalpha=1 has no meaning!!)
    end if 
    if (slt==0) then 
     print *, "gasalpha=", gasalpha 
    end if  
!!! end added Alex 12/26/08
 case ("time step")         ! ready to set up the time step dt 
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     dt=r_bulk(1) ! all other input values are ignored 
    else 
     dt=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "dt=", dt
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Added Alex 05/22/09 
    ! we are ready to read the temerature of the left wall. This is to be used with diffusive BCs 
  case ("temperature of the left wall")         ! ready to set up the temperature of the left wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     T_w_left=r_bulk(1) ! all other inputs values are ignored 
    else 
     T_w_left=300.0_DP ! set up to the default value (garbage! T_w_left=300K has no particular meaning!!)
    end if 
    if (slt==0) then 
     print *, "T_w_left=", T_w_left
    end if 
    ! we are ready to read the temerature of the right wall. This is to be used with diffusive BCs 
  case ("temperature of the right wall")         ! ready to set up the temperature of the right wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     T_w_right=r_bulk(1) ! all other inputs values are ignored 
    else 
     T_w_right=300.0_DP ! set up to the default value (garbage! T_w_right=300K has no particular meaning!!)
    end if 
    if (slt==0) then 
     print *, "T_w_right=", T_w_right
    end if  
!!!! end added Alex 05/22/09
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! added Alex 08/31/10 
case ("solution restart")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
   if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     need_to_restart = .TRUE.
    else
     need_to_restart = .FALSE.
    end if 
    if (slt==0) then 
     print *, "need_to_restart=", need_to_restart
    end if 
case ("restart time")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     restart_time_txt = "0.00000000"
    else
     restart_time_txt = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
     print *, "restart_time_txt=", restart_time_txt  
    end if  
!!!!! end added Alex 08/31/10 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
case ("number of OMP threads")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     Num_OMP_threads=i_bulk(1) ! all other inputs values are ignored 
    else 
     Num_OMP_threads=1 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "Num_OMP_threads=", Num_OMP_threads    
    end if  
   case default
    if (slt==0) then 
     print *, "Can not process:" // line
    end if  
  end select
 else
  if (slt==0) then 
   print *, line
  end if  
 end if   
 end do 
close(15)
end subroutine Set1DHOVParameters 

subroutine ReadIntegersFromLine (line_head,line,i_bulk, m, defval)
character (len=*), intent (in) :: line_head          ! line with parameter name 
character (len=*), intent (in) :: line               ! line with numbers to be processed    
integer (I4B), intent (in) :: defval 
integer (I4B), dimension (:), intent (out) :: i_bulk ! storage for the varaibles 
integer, intent (out) :: m                           ! counter how many records has been read   
integer :: code_line ! to use with Read statement 
!
character (len=10) :: fmtchar ! to keep format
!
code_line = 0
m=-1
do while (code_line == 0 )
   write (fmtchar,"(I3)") m
   m=m+1
   Read ( line, fmt = * ,  iostat = code_line) i_bulk(1:m+1)
   if ((code_line > 0) .and. (m==0)) then 
     print *, "Set1DHOVParameters: Error read line", code_line, & 
            "Parameter "//line_head//" may have an empty or improper numeric record. I use the default value", defval 
     i_bulk(m+1)=defval                       ! the default value for the parameter
     m=1
   end if 
end do
end subroutine ReadIntegersFromLine

subroutine ReadRealsFromLine (line_head,line,r_bulk, m, defval)
character (len=*), intent (in) :: line_head          ! line with parameter name 
character (len=*), intent (in) :: line               ! line with numbers to be processed    
real (DP), intent (in) :: defval 
real (DP), dimension (:), intent (out) :: r_bulk ! storage for the varaibles 
integer, intent (out) :: m                           ! counter how many records has been read   
integer :: code_line ! to use with Read statement 
!
character (len=10) :: fmtchar ! to keep format
!
code_line = 0
m=-1
do while (code_line == 0 )
   write (fmtchar,"(I3)") m
   m=m+1
   Read ( line, fmt = * ,  iostat = code_line) r_bulk(1:m+1)
   if ((code_line > 0) .and. (m==0)) then 
     print *, "Set1DHOVParameters: Error read line", code_line, & 
            "Parameter "//line_head//" may have an empty or improper numeric record. I use the default value", defval 
     r_bulk(m+1)=defval                       ! the default value for the parameter
     m=1
   end if 
end do
end subroutine ReadRealsFromLine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSolution1D3DHOV_spectral
!
! This procedure dumps the current solution on disk. 
!  
! The coefficients of characteristic variables are transformed into 
! the spectral coefficients and saved. 
! 
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Is mainlyt done to reduce the 
! main program. Do not call before all parameters are set including ftil1 and ftil2, Aml, k,s,max_deg, 
! etc.  
!
! most of the variables are taken directly from 
! common_variables_mod.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteSolution1D3DHOV_dgv
use common_variables_mod, only: xmesh,ftil,k,N,curr_time,rk,dt
!
intrinsic Trim, MatMul, Min, Adjustl
!
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
 ! ftil  -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  ftil(p,j,m)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- j is the cell in x
! first, we prepare the file name to store the solution
call MakeBaseNameHOV1D(file_name)
write (parchar, "(F11.8)") curr_time
file_name = trim(Adjustl(file_name))//trim(Adjustl(parchar))// "time"
!
file_name = trim(file_name)//"_1dx3duspec.dat"
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) k,N,xmesh,size(ftil,2),rk,dt ! the last number is the total number of velocity nodes
! now goes the solution, which we need to convert, first
! we write ftil1 ---> 
write(15) ftil 
!
close (15)
!
end subroutine WriteSolution1D3DHOV_dgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrorHOV1D
!
! This procedure dumps the evaluated error on disk. 
!  
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Is mainlyt done to reduce the 
! main program. Do not call before parameters k, s, max_deg, current_solution_dir, 
! current_solution_base_name, cfl,M,N,rk,curr_time 
! and many many others etc.  
!
! most of the variables are taken directly from 
! common_variables_mod.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteErrorHOV1D(err_sol,m_count)
use common_variables_mod, only: curr_time,initial_time
!
real (DP), dimension (0:,:), intent (in) :: err_sol ! array where the error (or norm of the solution is stored) 
integer (I4B), intent (in) :: m_count  !number of records to be saved          
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
!
! first, we prepare the file name to store the solution
call MakeBaseNameHOV1D(file_name)
file_name = trim(Adjustl(file_name))//"_err.dat"
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) initial_time,curr_time,m_count
! now goes the error itself
write(15) err_sol(0:m_count,:) 
!
close (15)
!
end subroutine WriteErrorHOV1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteMomentsSolHOV1D
!
! This procedure dumps the evaluated moments in discrete solution on disk. 
!  
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Is mainly done to reduce the 
! main program. Do not call before parameters k, s, max_deg, current_solution_dir, 
! current_solution_base_name, cfl,M,N,rk,curr_time 
! and many many others etc.  
!
! most of the variables are taken directly from 
! common_variables_mod.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteMomentsSolHOV1D(moments_sol,m_count)
use common_variables_mod, only: curr_time,initial_time
!
real (DP), dimension (:,0:,:), intent (in) :: moments_sol ! array where the three moments are stored 
                                                          ! moments_sol (i,j,k) 
                                                          ! i --- ## of moment
                                                          ! j = 0,1,...  --- # of the moment of saving. 
                                                          ! k = 1,2 --- for each component of the solution                                                        
integer (I4B), intent (in) :: m_count  !number of records to be saved          
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
! first, we prepare the file name to store the solution
call MakeBaseNameHOV1D(file_name)
file_name = trim(Adjustl(file_name))//"_mosol.dat"  ! this will keep moments of spectal solution
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) initial_time,curr_time,m_count
! now goes the error itself
write(15) moments_sol(:,0:m_count,:) 
!
close (15)
!
end subroutine WriteMomentsSolHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteMomentsExactHOV1D
!
! This procedure dumps the evaluated moments in exact solution on disk. 
!  
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Is mainly done to reduce the 
! main program. Do not call before parameters k, s, max_deg, current_solution_dir, 
! current_solution_base_name, cfl,M,N,rk,curr_time 
! and many many others etc.  
!
! most of the variables are taken directly from 
! common_variables_mod.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteMomentsExactHOV1D(moments_exact,m_count)
use common_variables_mod, only: curr_time,initial_time
!
real (DP), dimension (:,0:,:), intent (in) :: moments_exact ! array where the three moments are stored 
                                                          ! moments_exact (i,j,k) 
                                                          ! i --- ## of moment
                                                          ! j = 0,1,...  --- # of the moment of saving. 
                                                          ! k = 1,2 --- for each component of the solution
                                                          
                                                          
integer (I4B), intent (in) :: m_count  !number of records to be saved          
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
! first, we prepare the file name to store the solution
call MakeBaseNameHOV1D(file_name)
file_name = trim(Adjustl(file_name))//"_moexa.dat"  ! this will keep moments of spectal solution
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) initial_time,curr_time,m_count
! now goes the error itself
write(15) moments_exact(:,0:m_count,:) 
!
close (15)
!
end subroutine WriteMomentsExactHOV1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameHOV1D 
!
! This function creates a recognizable name for the file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k, s, max_deg, current_solution_dir, current_solution_base_name, cfl, &
!          M,N,rk,curr_time,initial_time,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameHOV1D (file_name)
use common_variables_mod, only: k, s, max_deg, current_solution_dir, current_solution_base_name, cfl, &
          M,N,rk,curr_time,initial_time,x_nonuniform_mesh_type,u_nonuniform_mesh_type,mesh_x_uniform,mesh_u_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_solution_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k
file_name = trim(file_name)//trim(Adjustl(parchar))//"k"
write (parchar, "(I3)") s
file_name = trim(file_name)//trim(Adjustl(parchar))//"s"
write (parchar, "(I3)") max_deg
file_name = trim(file_name)//trim(Adjustl(parchar))//"mdgr"
write (parchar, "(I1)") rk
file_name = trim(file_name)//trim(Adjustl(parchar))//"rk"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"U"
else
  select case (x_nonuniform_mesh_type)
  case (1)
  file_name = trim(file_name)//"G"           ! type 1: nonuniform mesh with gauss nodes
  case default
  file_name = trim(file_name)//"U"           ! unsupported type of mesh --- uniform by default
  end select 
end if 
write (parchar, "(I5)") M
file_name = trim(file_name)//trim(Adjustl(parchar))//"M"
if (mesh_u_uniform) then 
file_name = trim(file_name)//"U"
else
  select case (u_nonuniform_mesh_type)
  case (1)
  file_name = trim(file_name)//"G"           ! type 1: nonuniform mesh with gauss nodes
  case default
  file_name = trim(file_name)//"U"           ! unsupported type of mesh --- uniform by default
  end select 
end if 
write (parchar, "(F8.5)") cfl
file_name = trim(file_name)//trim(Adjustl(parchar))//"cfl"
!
end subroutine MakeBaseNameHOV1D
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Added Alex 08/27/2010 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RestartSolution1D3DHOV
!
! This procedure reads a solution from the disk. 
!  
! The spectral coefficients are transformed into the characteristic variables and saved in ftil1, ftil2. 
! 
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Do not call before all parameters are 
! set including ftil1 and ftil2, Aml, k,s,max_deg, 
! etc.  
!
! most of the variables are taken directly from 
! common_variables_mod.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RestartSolution1D3DHOV(restart_time_txt)
use common_variables_mod, only: xmesh,ftil,k,N,rk,dt
!!!
use DGV_commvar, only: nodes_u

!!!
!
intrinsic Trim, MatMul, Min, Adjustl
!
character (len=20) :: restart_time_txt ! the time at which the solution need to be restored
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the parameter 
integer (I4B) :: dump1,dump2     ! string to read some garbage bits from the file
integer (I4B) :: M,M1,k1,N1 ! local counters
integer :: code_line !! local variable to keep the status of reading from file

integer :: loc_alloc_stat ! to keep allocation status


M=size(nodes_u,1) ! this is the total number of nodes in the velocity variable
!  first, we prepare the space for ftil1 and ftil2:   ! ftil(p,m,j):
                                                      ! -- p is the index in the basis functions in x 
                                                      ! -- m is the index in the basis functions in u
                                                      ! -- j is the cell in x
allocate (ftil(0:k,M,N), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "RestartSolution1D3DHOV: Allocation error for variable (ftil)"
     end if 
     !
! Second, we prepare the file name to read the solution
call MakeBaseNameHOV1D(file_name)
file_name = trim(Adjustl(file_name))//trim(Adjustl(restart_time_txt))// "time"
file_name = trim(file_name)//"_1dx3duspec.dat"
!
! now we open the file for record and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) k1,N1,xmesh,M1,rk,dt
! now goes the solution, which we need to convert after reading 
! we read ftil1 and ftil2 ---> 
read (15,  iostat = code_line) ftil ! notice that we need two independent statemenets because this is how the solution was recorded. Each write statement creates a line. Each read statement reads from a new line
close (15)
!!
! A quick check if the saved solution is stored correctly.
if ((k1 /= k) .or. (N1 /= N) .or. (M1 /= M)) then 
 print *, "RestartSolution1D3DHOV: Constants k,M and N in the saved solutions are incompatible."
 stop
end if 
!
end subroutine RestartSolution1D3DHOV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RestartSolutionMTS1D3D
!
! This procedure reads a solution from the disk. 
!
! This subroutine allocates the storage for the solution  
!
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Do not call before all parameters are 
! set especially f -- this has to have been allocated correctly, its dimensions will be used to 
! set dimensions of other arrays. If the dimensions are wrong, the arrays will not read correctly  
! 
!
! most of the variables are taken directly from 
! commvar.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine RestartSolutionMTS1D3D(restart_time_txt)

use common_variables_mod, only: xmesh,k,N,rk,ftil,ftil1,ftil2,ftil3,ftil4,&
                                frhs1,frhs2,frhs3,frhs4,frhs5,dt   
!!!
use DGV_commvar, only: nodes_u
!
!
intrinsic Trim, Adjustl
!
character (len=20) :: restart_time_txt ! the time at which the solution need to be restored
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the parameter 
integer (I4B) :: dump1,dump2     ! string to read some garbage bits from the file
integer (I4B) :: M,M1,k1,N1 ! scrap variables to keep important comstants
integer :: code_line !! local variable to keep the status of reading from file
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! First, we need to allocate the storages
M=size(nodes_u,1) ! number of velocity nodes, also the second dimension of ftil
SELECT CASE(rk)
	CASE(1)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:k,M,N),ftil(0:k,M,N),ftil1(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs1,ftil,ftil1)"
		stop
	  END IF
	CASE(2)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:k,M,N),frhs2(0:k,M,N),ftil(0:k,M,N),ftil1(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs1,frhs2,ftil,ftil1)"
		stop
	  END IF
	CASE(3)
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:k,M,N),ftil(0:k,M,N),ftil1(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs1,ftil,ftil1)"
		stop
	  END IF
	ALLOCATE(frhs2(0:k,M,N),frhs3(0:k,M,N),ftil2(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs2,ffrhs3,ftil2)"
		stop
	  END IF
	CASE(4) 
	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:k,M,N),ftil(0:k,M,N),ftil1(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs1,ftil1,ftil)"
		stop
	  END IF
	ALLOCATE(frhs2(0:k,M,N),ftil2(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs2,ftil2)"
		stop
	  END IF
	  ALLOCATE(frhs3(0:k,M,N),frhs4(0:k,M,N),ftil3(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs3,ftil3)"
		stop
	  END IF
    CASE(5)
      !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
   	  ALLOCATE(frhs1(0:k,M,N),ftil1(0:k,M,N),ftil(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs1,ftil1,ftil)"
		stop
	  END IF
	  ALLOCATE(frhs2(0:k,M,N),ftil2(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs2,ftil2)"
		stop
	  END IF
	ALLOCATE(frhs3(0:k,M,N),ftil3(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs3,ftil3)"
		stop
	  END IF
	  ALLOCATE(frhs4(0:k,M,N),frhs5(0:k,M,N),ftil4(0:k,M,N), stat = loc_alloc_stat)
	  IF(loc_alloc_stat > 0) THEN
		PRINT *,"RestartSolutionMTS1D3D: Allocation error for variable (frhs4,frhs5,ftil4)"
		stop
	  END IF
	CASE default
			PRINT *, "RestartSolutionMTS1D3D: The value of (rk) must be from 1 to 5. No such RK or MTS methods implemented"
			stop
END SELECT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Second, we prepare the file name to read the solution
call MakeBaseNameHOV1D(file_name)
file_name = trim(Adjustl(file_name))//trim(Adjustl(restart_time_txt))// "time"
file_name = trim(file_name)//"_1D3Dsol.dat" ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")

read (15,  iostat = code_line) k1,N1,xmesh,M1,rk,dt
! now goes the solution, which we need to convert after reading 
!!

! A quick check if the saved solution is stored correctly.
if ((k1 /= k) .or. (N1 /= N) .or. (M1 /= M)) then 
 print *, "RestartSolutionMTS1D3D: Constants k,M and N in the saved solutions are incompatible."
 close(15)
 stop
end if 

! Next we save the solution .... 
SELECT CASE(rk)
	CASE(1)
read (15,  iostat = code_line) ftil ! notice that we need two independent statemenets because this is how the solution was recorded. Each write statement creates a line. Each read statement reads from a new line
    CASE(2)
read (15) ftil,ftil1,frhs1
    CASE(3)
read (15) ftil,ftil1,ftil2,frhs1,frhs2
    CASE(4)
read (15) ftil,ftil1,ftil2,ftil3,frhs1,frhs2,frhs3
    CASE(5)
read (15) ftil,ftil1,ftil2,ftil3,ftil4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
!
end subroutine RestartSolutionMTS1D3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSolutionMTS1D3D
!
! This procedure saves a solution to hard drive. 
!
! the directory path is set through parameter current_solution_dir
! the name is composed from the curr_sol_base_name and a values of different parameters 
!
! This procedure is dependent on the main program! Do not call before all parameters are 
! set especially ftil -- this has to have been allocated correctly, its dimensions will be used to 
! set dimensions of other arrays. If the dimensions are wrong, the arrays will not read correctly  
! 
!
! most of the variables are taken directly from 
! commvar.
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteSolutionMTS1D3D

use common_variables_mod, only: xmesh,k,N,rk,ftil,ftil1,ftil2,ftil3,ftil4,&
                                frhs1,frhs2,frhs3,frhs4,frhs5,dt,curr_time,&
                                initial_time   
!!!
!
intrinsic Trim, Adjustl
!
character (len=20) :: restart_time_txt ! the time at which the solution need to be restored
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the parameter 
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! ftil  -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  ftil(p,j,m)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- j is the cell in x
! first, we prepare the file name to store the solution
call MakeBaseNameHOV1D(file_name)
write (parchar, "(F11.8)") curr_time
file_name = trim(Adjustl(file_name))//trim(Adjustl(parchar))// "time"
!
file_name = trim(file_name)// "_1D3Dsol.dat"
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) k,N,xmesh,size(ftil,2),rk,dt ! the last number is the total number of velocity nodes
! now goes the solution, which we need to convert, first
! we write ftil and other arrays ---> 
SELECT CASE(rk)
	CASE(1)
write (15) ftil ! notice that we need two independent statemenets because this is how the solution was recorded. Each write statement creates a line. Each read statement reads from a new line
    CASE(2)
     if (curr_time<=initial_time) then
      write (15) ftil
     else  
      write (15) ftil,ftil1,frhs1
     end if  
    CASE(3)
     if (curr_time<=initial_time) then
      write (15) ftil
     else  
      write (15) ftil,ftil1,ftil2,frhs1,frhs2
     end if  
    CASE(4)
     if (curr_time<=initial_time) then
      write (15) ftil
     else  
      write(15) ftil,ftil1,ftil2,ftil3,frhs1,frhs2,frhs3
     end if  
    CASE(5)
     if (curr_time<=initial_time) then
      write (15) ftil
     else  
      write (15) ftil,ftil1,ftil2,ftil3,ftil4,frhs1,frhs2,frhs3,frhs4
     end if  
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
!
end subroutine  WriteSolutionMTS1D3D

end module read_write_tools_mod

