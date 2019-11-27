!
!  readwrite.f90
!
! Alex 09/12/14
!
! This modle contains useful procedures for reading and writing from the disk
! The subroutines will be used in other functions of the BGKVDCF 0D Driver
! Many of the subroutines are hard linked to the variables of BGKVDCF0D_commvar.f90
! however, some subroutines can be used independently.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

module BGKVDCF0D_readwrite
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetUWbgkParams(pfname)
!
! This subroutine reads the variables from the 
! parameter file on the hard drive. The variables 
! from the common variabl block are accessed directly
!
!!!!!!!!!!!!!!!!!!
subroutine SetUWbgkParams(pfname)
use BGKVDCF0D_commvar
intrinsic Real,Index, Len, Scan, Trim

character (len=*), intent (in) :: pfname ! Name of the file where the parameters are stored 
character (len=132) :: line              ! string to keep one record 
character (len=10) :: fmtchar            ! string to keep format 
character (len=50) :: line_head          ! string to keep the line header 
integer (I4B) :: code_file, code_line          ! get the return code from READ
integer (I4B) :: m_count                            ! dump counter 
integer (I4B) :: pos                 ! to store position within the string 
integer (I4B), dimension (20) :: i_bulk  ! to store temporarily integers that has been read from the parameter file 
real (DP), dimension (20) :: r_bulk      ! to store temporarily reals that has been  read from the parameter file 
integer (I4B) :: loc_alloc_stat ! some dump varaible
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
   case ("degree of local Legendre basis in x")           ! ready to set up local order in variable x
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (k_c_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (k_c_list)"
     stop
     end if
     !
     k_c_list=i_bulk(1:m_count)
    end if  
    print *, "k_c_list=", k_c_list 
   case ("number of G-L boundary nodes in x")           ! ready to set up boundary order in variable x
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count,2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (k_b_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (k_b_list)"
     stop
     end if
     !
     k_b_list=i_bulk(1:m_count)
    end if  
    print *, "k_b_list=", k_b_list 
   case ("degree of local Legendre basis in t")           ! ready to set up local order in variable t
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (d_c_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (d_c_list)"
     stop
     end if
     !
     d_c_list=i_bulk(1:m_count)
    end if  
    print *, "d_c_list=", d_c_list 
   case ("number of G-L boundary nodes in t")           ! ready to set up boundary order in variable t
    !!! First we read the paramters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count,2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the parameters we have just read !!! 
    if (m_count > 0) then 
     allocate (d_b_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (d_b_list)"
     stop
     end if
     !
     d_b_list=i_bulk(1:m_count)
    end if  
    print *, "d_b_list=", d_b_list  
   case ("left endpoint in x")         ! ready to set up left endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     x_L=0.0_DP ! set up to the default value 
    end if  
    print *, "x_L=", x_L 
   case ("right endpoint in x")         ! ready to set up right endpoint in x
    !!! First we read the parameters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     x_R=r_bulk(1) ! all other input values are ignored 
    else 
     x_R=1.0_DP ! set up to the default value 
    end if  
    print *, "x_R=", x_R
   
   case ("uniform mesh in x")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
   if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_x_uniform = .TRUE.
    else
     mesh_x_uniform = .FALSE.
    end if 
     print *, "mesh_x_uniform=", mesh_x_uniform
  
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
    print *, "N_list=", N_list
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   case ("conditions on the left boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_leftbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_leftbcond=1 ! set up to the default value 
    end if  
    print *, "selected_leftbcond=", selected_leftbcond  
   case ("conditions on the right boundary")         ! ready to set up selected type of boundary conditions
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_rightbcond=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_rightbcond=1 ! set up to the default value 
    end if  
    print *, "selected_rightbcond=", selected_rightbcond  
   case ("type of exact solution")         ! ready to set up selected type of boundary data
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    selected_exact_sol=i_bulk(1) ! all other inputs values are ignored 
    else 
    selected_exact_sol=1 ! set up to the default value 
    end if  
    print *, "selected_exact_sol=", selected_exact_sol  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case ("current solution save directory")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_dir = "solution/"
    else
     current_solution_dir = Trim(Adjustl(line))
    end if 
    print *, "current_solution_dir=", current_solution_dir 
   case ("current solution base name")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_solution_base_name = "bgk1d"
    else
     current_solution_base_name = Trim(Adjustl(line))
    end if 
    print *, "current_solution_base_name=", current_solution_base_name 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case ("initial time")       ! ready to set up the initial time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     t_L=r_bulk(1) ! all other input values are ignored 
    else 
     t_L=0.0_DP ! set up to the default value 
    end if  
    print *, "t_L=", t_L
   case ("final time")         ! ready to set up the end time
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     t_R=r_bulk(1) ! all other input values are ignored 
    else 
     t_R=1.0_DP ! set up to the default value 
    end if  
    print *, "t_R=", t_R
   case ("instances to evaluate error")  ! ready to set up the number of error evaluations
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1000) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_eval_error=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_eval_error=1000 ! set up to the default value 
    end if  
    print *, "num_eval_error=", num_eval_error 
   case ("time step")         ! ready to set up the time step dt 
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     dt=r_bulk(1) ! all other input values are ignored 
    else 
     dt=1.0_DP ! set up to the default value 
    end if 
    print *, "dt=", dt
   case ("instances to save solution")  ! ready to set up the number solution recordings
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 10) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_save_solution=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_save_solution=10 ! set up to the default value 
    end if  
    print *, "num_save_solution=", num_save_solution 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   case ("type of nonuniform mesh in x")  ! ready to set up the type of nonuniform mesh in x
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    x_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    x_nonuniform_mesh_type=1 ! set up to the default value 
    end if  
    print *, "x_nonuniform_mesh_type=", x_nonuniform_mesh_type 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
   case ("gauss order for moments in x")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_x_gauss_order=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_x_gauss_order=1 ! set up to the default value 
    end if  
    print *, "moments_x_gauss_order=", moments_x_gauss_order 
   case ("mesh refinement in x for moments")  ! ready to set up the coefficient of mesh refinement in x for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    moments_refine_x=i_bulk(1) ! all other inputs values are ignored 
    else 
    moments_refine_x=1 ! set up to the default value 
    end if  
    print *, "moments_refine_x=", moments_refine_x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case ("temperature of the left wall")         ! ready to set up the temperature of the left wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then 
     T_w_left=r_bulk(1) ! all other inputs values are ignored 
    else
     T_w_left=300.0_DP ! set up to the default value (garbage! T_w_left=300K has no particular meaning!!)
    end if
    print *, "T_w_left=", T_w_left
    ! we are ready to read the temerature of the right wall. This is to be used with diffusive BCs 
  case ("temperature of the right wall")         ! ready to set up the temperature of the right wall
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, 300.0_DP) ! last parameter is the default value for the paramenter to be read (default R=1 is garbage!)
    if (m_count > 0) then
     T_w_right=r_bulk(1) ! all other inputs values are ignored
    else
     T_w_right=300.0_DP ! set up to the default value (garbage! T_w_right=300K has no particular meaning!!)
    end if  
    print *, "T_w_right=", T_w_right
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 case ("solution restart")         ! ready to set up mesh in x is uniform parameter 
    !!! We read the parameter from the input line 
   if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     need_to_restart = .TRUE.
    else
     need_to_restart = .FALSE.
    end if 
     print *, "need_to_restart=", need_to_restart
case ("restart time")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     restart_time_txt = "0.00000000"
    else
     restart_time_txt = Trim(Adjustl(line))
    end if 
    print *, "restart_time_txt=", restart_time_txt 
case ("number of OMP threads")  ! ready to set up the order of gauss method in x variable for the evaluation of moments
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    Num_OMP_threads=i_bulk(1) ! all other inputs values are ignored 
    else 
    Num_OMP_threads=1 ! set up to the default value 
    end if 
    print *, "Num_OMP_threads(Driver)=", Num_OMP_threads    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   case default 
    print *, "Can not process:" // line_head // "=" // line
  end select
 else
 print *, line
 end if   
 end do 
close(15)
end subroutine SetUWbgkParams 

subroutine ReadIntegersFromLine (line_head,line,i_bulk, m, defval)
character (len=*), intent (in) :: line_head          ! line with parameter name 
character (len=*), intent (in) :: line               ! line with numbers to be processed    
integer (I4B), intent (in) :: defval 
integer (I4B), dimension (:), intent (out) :: i_bulk ! storage for the varaibles 
integer (I4B), intent (out) :: m                           ! counter how many records has been read   
integer (I4B) :: code_line ! to use with Read statement 
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
integer (I4B), intent (out) :: m                           ! counter how many records has been read   
integer (I4B) :: code_line ! to use with Read statement 
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameBGK0D 
!
! This function creates a recognizable name for the file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detatch from the main program. 
! 
! Do not call until parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the DGV_commvar and BGKVDCF0D_commv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameBGK0D (file_name)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform

use BGKVDCF0D_commvar, only: k_c, current_solution_dir, current_solution_base_name, &
          N,x_nonuniform_mesh_type,mesh_x_uniform
          !
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_solution_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameBGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameBGK0D2
!
! Modified by Craig 12/31/12, Alex 10/20/2014
! This function creates a recognizable name for the file.
! a lot of parameters is cited in the file name. This function will set up the basic
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detatch from the main program.
!
! Do not call until parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the DGV_commvar and BGKVDCF0D_commvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameBGK0D2 (file_name, directory)

use DGV_commvar, only:  su,sv,sw,Mu,Mv,Mw, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform

use BGKVDCF0D_commvar, only: k_c, current_solution_dir, current_solution_base_name, &
          N,x_nonuniform_mesh_type,mesh_x_uniform 
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name
character (len=8), intent (in) :: directory ! the variable that holds the directory name
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(adjustl(directory))//trim(current_solution_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameBGK0D2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameChnksBGK0D
!
! This function creates a recognizable name for the file. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameChnksBGK0D (file_name,i)

use DGV_commvar, only:  su,sv,sw,Mu,Mv,Mw,&
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform

use BGKVDCF0D_commvar, only: k_c, current_solution_dir, current_solution_base_name,&
         N,x_nonuniform_mesh_type, mesh_x_uniform
!
intrinsic Trim, Min, Adjustl
!
integer (I4B) :: i ! the index of the chunk (starts from zero)
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
write (parchar, "(I3.3)") i
file_name = trim(current_solution_dir)//trim(current_solution_base_name)//"ch"//trim(Adjustl(parchar))//"_"
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameChnksBGK0D

!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperBGK0D 
!
! This function creates a recognizable name for the file to store A-operator in one file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperBGK0D (file_name)

use DGV_commvar, only:  su,sv,sw,Mu,Mv,Mw,current_solution_dir, current_Aoperator_base_name,&
          Mu, Mv, Mw, u_nonuniform_mesh_type, mesh_u_uniform,&
          v_nonuniform_mesh_type, w_nonuniform_mesh_type, mesh_v_uniform, mesh_w_uniform
          
use BGKVDCF0D_commvar, only: k_c, N, x_nonuniform_mesh_type, mesh_x_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameAoperBGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperChnksBGK0D 
!
! This function creates a recognizable name for the files to store chunks of operator A. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters k_c, su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, x_nonuniform_mesh_type, u_nonuniform_mesh_type, mesh_x_uniform, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperChnksBGK0D (file_name,i)

use DGV_commvar, only:  su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_Aoperator_base_name, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform

use BGKVDCF0D_commvar, only:k_c,N,x_nonuniform_mesh_type,mesh_x_uniform
!
intrinsic Trim, Min, Adjustl
!
integer (I4B) :: i ! the index of the chunk (starts from zero)
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
write (parchar, "(I3.3)") i
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)//"ch"//trim(Adjustl(parchar))//"_"
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") k_c
file_name = trim(file_name)//trim(Adjustl(parchar))//"kc"
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
write (parchar, "(I5)") N
file_name = trim(file_name)//trim(Adjustl(parchar))//"N"
if (mesh_x_uniform) then 
file_name = trim(file_name)//"XU"
else
write (parchar, "(I1)") x_nonuniform_mesh_type
file_name = trim(file_name)//"X"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mu
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mu"
if (mesh_u_uniform) then 
 file_name = trim(file_name)//"UU"
else
 write (parchar, "(I1)") u_nonuniform_mesh_type
 file_name = trim(file_name)//"U"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mv
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mv"
if (mesh_v_uniform) then 
 file_name = trim(file_name)//"VU"
else
 write (parchar, "(I1)") v_nonuniform_mesh_type
 file_name = trim(file_name)//"V"//trim(Adjustl(parchar))
end if 
write (parchar, "(I5)") Mw
file_name = trim(file_name)//trim(Adjustl(parchar))//"Mw"
if (mesh_w_uniform) then 
 file_name = trim(file_name)//"WU"
else
 write (parchar, "(I1)") w_nonuniform_mesh_type
 file_name = trim(file_name)//"W"//trim(Adjustl(parchar))
end if 
!
end subroutine MakeBaseNameAoperChnksBGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteSol_BGK0D(suff)
!
! This subroutine will damp on the disk the current state of the solution along with some 
! parameters.
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteSol_BGK0D (suff)

use BGKVDCF0D_commvar, only:  f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5,&
                    time,dt,rkmts

intrinsic Trim
!
character (len=*) :: suff ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D2(file_name, "results/")
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")

write (15) time,dt ! record the current time and the used dt
write (15) rkmts ! record the order of the time integrator.
write (15) size(f,1) ! record the length of the arrays.
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
write (15) f
    CASE(2)
write (15) f,f1,frhs1
    CASE(3)
write (15) f,f1,f2,frhs1,frhs2
    CASE(4)
write (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
write (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
end subroutine WriteSol_BGK0D  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadSol_BGK0D(suff)
!
! This subroutine will read the solution from the disk.
! the solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadSol_BGK0D (suff) 

use BGKVDCF0D_commvar, only:  f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5,&
                    time,dt,rkmts

intrinsic Trim
!
character (len=15) :: suff ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: m ! the length of the arrays
integer (I4B) :: loc_alloc_stat ! some dump varaible


! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")

read (15) time,dt ! record the current time and the used dt
read (15) rkmts ! record the order of the time integrator.
read (15) m ! record the length of the arrays.
! A quick check that the storages are the correct size: 
if (m /= size(f,1)) then 
 print *, "ReadSol_DGV: trying to read a storage of wrong size, m/= size(f,1)"
 close (15)
 stop
end if  
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
read (15) f ! the space for this one needs to be allocated already ... 
    CASE(2)
! We now need to prepare the storage for the data.
allocate (f1(1:m), frhs1(1:m), frhs2(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,frhs1,frhs2"
  close(15)
  stop
  end if
read (15) f,f1,frhs1
    CASE(3)
! We now need to prepare the storage for the data.
allocate (f1(1:m),f2(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,frhs1,frhs2,frhs3"
  close(15)
  stop
  end if
read (15) f,f1,f2,frhs1,frhs2
    CASE(4)
  allocate (f1(1:m),f2(1:m),f3(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,frhs1,frhs2,frhs3,frhs4"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
  allocate (f1(1:m),f2(1:m),f3(1:m),f4(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m),frhs5(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
end subroutine ReadSol_BGK0D  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KeepTrackL1_err_BGK0D(L1_err)
!
! This subroutine saves the values of the relative L1 norm of the deviation of the current solution from the local maxwellian 
!
! Takes the parameters from commvar.mod
!!!!!!!!!!!!!!!!!!

subroutine KeepTrackL1_err_BGK0D(f)

use DGV_distributions_mod
use DGV_dgvtools_mod

use DGV_commvar, only: run_mode,nodes_u, nodes_v, &
       nodes_w,nodes_gwts

use BGKVDCF0D_commvar, only: L1_a, L1_t, time,L1_count,t_L, t_R, &
       num_save_solution,L1_record_period,L1_next_time_record, &
       L1_err_eval_period,L1_next_time_eval,run_mode_array


intrinsic Trim, SUM, ABS
!
real (DP), dimension (:), intent (in) :: f ! the current solution  
!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (size(f,1)) :: locFm  ! a scratch storage for the maxwellian 
real (DP) :: LDens,LUbar,LVbar,LWbar,LTempr ! scratch variables to store macroparameters
real (DP) :: L1_err ! relative L1 error of the deviation of distribution function from the local maxwellian computed elsewhere
integer (I4B), parameter :: c = 1000 ! the total number of records in the error file
integer :: loc_alloc_stat
integer (I4B) :: n ! scrap counter
character (len=132) :: file_name ! the variable to store the file name
!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! EVALUATE L1_error here:
call MassCheckRec (f,LDens,LUbar,LVbar,LWbar,LTempr) ! First, we evaluate the macroparamters of the solution.
locFm = maxwelveldist(LTempr,LUbar,LVbar,LWbar,LDens,nodes_u,nodes_v,nodes_w) ! now we populate the maxwellian with the same macroparamters.
locFm = f-locFm ! 
L1_err = SUM(ABS(locFm)*nodes_gwts)/LDens! evaluate the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!
! check is the error arrays are allocated. if not -- allocate them and reset the counter
if (size(L1_a,1)<2) then
 allocate (L1_a(1:c),L1_t(1:c),run_mode_array(1:c),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then
  print *, "KeepTrackL1_err_DGV: Allocation error for variable (L1_a(1:c),L1_t(1:c))"
  stop
  end if
L1_count = 1
L1_a(1) = L1_err
L1_t(1) = time
run_mode_array(1) = run_mode
!!!!!!!!
! we also set up the time intervals for evaluation and saves
L1_record_period = (t_R-t_L)/Real(num_save_solution,DP)
L1_next_time_record = time + L1_record_period
L1_err_eval_period = (t_R-t_L)/Real(c,DP)
L1_next_time_eval = time + L1_err_eval_period
end if
! once the arrays are allocated, we periodically save the solution
if (time >= L1_next_time_eval) then ! check if it is time to record
if (L1_count < c) then
       L1_count=L1_count+1
       L1_next_time_eval = time + L1_err_eval_period
       L1_a(L1_count) = L1_err
       L1_t(L1_count) = time
	   run_mode_array(L1_count) = run_mode
end if
end if
if (time >= L1_next_time_record) then ! check if it is time to save on disk.
 L1_next_time_record = time + L1_record_period
 ! Writing the error on disk:
 ! first, we prepare the file name to store the solution
 call MakeBaseNameBGK0D2(file_name, "L1_errs/")
 file_name = trim(Adjustl(file_name))//"_L1err.txt"  ! this file will keep the array
 !
 ! now we open the file in formatted write regime for record and save some stuff in it:
 open (15, file=file_name, position ="REWIND", action="WRITE", &
                    form="FORMATTED")
 write(15, *) "This file contains the values of relative L1 error of deviation from local maxwellian and the time stamp:"
 write(15, *) "Time		L1_err		mode"
 ! next we print array using nnn entries in a row
 do n=1,L1_count
   write (15,"(F18.14, A2, F18.15, A2, I2)" ) L1_t(n), "," , L1_a(n), ",", run_mode_array(n)
 end do
 close (15)
end if
end subroutine KeepTrackL1_err_BGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! KeepTrackCco_BGK0D(L1_err)
!
! This subroutine saves the values of the coefficients of the basis used for 
! the velocity dependent collision frequency  
!
! Takes the parameters from commvar.mod
!!!!!!!!!!!!!!!!!!

subroutine KeepTrackCco_BGK0D

use DGV_distributions_mod
use DGV_dgvtools_mod

use DGV_commvar, only: Cco,Order_nu 

use BGKVDCF0D_commvar, only: Cco_rec, Cco_count, time, t_L, t_R, &
       num_save_solution, Cco_record_period,Cco_next_time_record,&
       Cco_eval_period,Cco_next_time_eval

intrinsic Trim, SUM, ABS
!
integer (I4B), parameter :: c = 1000 ! the total number of records in the error file
integer :: loc_alloc_stat
integer (I4B) :: i,n ! scrap counter
character (len=132) :: file_name! the variable to store the file name
character (len=500) :: textline! the variable to write date 
character (len=40) :: char_token! the variable that keeps just one number

!!
!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!
! check is the error arrays are allocated. if not -- allocate them and reset the counter
if (size(Cco_rec,1)<2) then
 allocate (Cco_rec(1:Order_Nu+1,1:c),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then
  print *, "KeepTracCco_BGK0D: Allocation error for variable (Cco_rec)"
  stop
  end if
Cco_count = 1
Cco_rec(1,1) = time 
do i=1,Order_nu
 Cco_rec(i+1,1) = Cco(i)
end do 
!!!!!!!!
! we also set up the time intervals for evaluation and saves
Cco_record_period = (t_R-t_L)/Real(num_save_solution,DP)
Cco_next_time_record = time + Cco_record_period
Cco_eval_period = (t_R-t_L)/Real(c,DP)
Cco_next_time_eval = time + Cco_eval_period
end if
! once the arrays are allocated, we periodically save the solution
if (time >= Cco_next_time_eval) then ! check if it is time to record
if (Cco_count < c) then
       Cco_count=Cco_count+1
       Cco_next_time_eval = time + Cco_eval_period
       Cco_rec(1,Cco_count) = time 
       do i=1,Order_nu
        Cco_rec(i+1,Cco_count) = Cco(i)
       end do 
end if
end if
if (time >= Cco_next_time_record) then ! check if it is time to save on disk.
 Cco_next_time_record = time + Cco_record_period
 ! Writing the error on disk:
 ! first, we prepare the file name to store the solution
 call MakeBaseNameBGK0D2(file_name, "L1_errs/")
 file_name = trim(Adjustl(file_name))//"_Cco.txt"  ! this file will keep the array
 !
 ! now we open the file in formatted write regime for record and save some stuff in it:
 open (15, file=file_name, position ="REWIND", action="WRITE", &
                    form="FORMATTED")
 write(15, *) "This file contains the values of coefficients of velocity dependent collision frequency and the time stamp:"
 write(15, *) "Time		Coefficients 1, ... , Order_nu"
 ! next we print array using nnn entries in a row
 do n=1,Cco_count
  ! for each row of the array we will create a character line that keeps the data. Then we will write it to the file with a default format
  textline = ""
  do i=1,Order_Nu+1
   write (char_token, "(F18.14)") Cco_rec(i,n)
  textline = trim(Adjustl(textline))//trim(Adjustl(char_token))//", "  
  end do  
  write(15, *) trim(Adjustl(textline)) 
 end do
 close (15)
end if
end subroutine KeepTrackCco_BGK0D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RestartSolutionBGK0D
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

subroutine RestartSolutionBGK0D(restart_time_txt)
use BGKVDCF0D_commvar, only: f, f1, f2, f3, f4, frhs1, frhs2, frhs3, frhs4, &
							frhs5, rkmts,dt,time
!
intrinsic Trim, Adjustl
!
character (len=20) :: restart_time_txt ! the time at which the solution need to be restored
character (len=132) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the parameter 
integer (I4B) :: dump1,dump2     ! string to read some garbage bits from the file
integer :: m ! local counters
integer :: code_line !! local variable to keep the status of reading from file
integer :: loc_alloc_stat ! to keep allocation status
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ALEX: Check if this commented area needs to be removed...
! First, we need to allocate the storages
!SELECT CASE(rkmts)
!	CASE(1)
!	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
!   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)!
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1)"
!		stop
!	  END IF
!	CASE(2)
!	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
!  	  ALLOCATE(frhs1(1:size(f,1)),frhs2(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
!		stop
!	  END IF
!	CASE(3)
!	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
!   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
!		stop
!	  END IF
!	ALLOCATE(frhs2(1:size(f,1)),frhs3(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
!		stop
!	  END IF
!	CASE(4) 
!	  !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
!!   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)!
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
!		stop
!	  END IF
!	ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
!		stop
!	  END IF
!	  ALLOCATE(frhs3(1:size(f,1)),frhs4(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs3,f3)"
!		stop
!	  END IF
!    CASE(5)
!      !!! ALLOCATE THE ARRAYS NEEDED FOR THE MTS STEPS
!   	  ALLOCATE(frhs1(1:size(f,1)),f1(1:size(f,1)), stat = loc_alloc_stat)!
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs1,f1)"
!		stop
!	  END IF
!	  ALLOCATE(frhs2(1:size(f,1)),f2(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs2,f2)"
!		stop
!	  END IF
!	ALLOCATE(frhs3(1:size(f,1)),f3(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs3,f3)"
!		stop
!	  END IF
!	  ALLOCATE(frhs4(1:size(f,1)),frhs5(1:size(f,1)),f4(1:size(f,1)), stat = loc_alloc_stat)
!	  IF(loc_alloc_stat > 0) THEN
!		PRINT *,"PrepareMTS_SH_DGV: Allocation error for variable (frhs4,f4)"
!		stop
!	  END IF
!	CASE default
!			PRINT *, "PrepareMTS_SH_DGV: The value of (rkmts) must be from 1 to 5. No such RK or MTS methods implemented"
!			stop
!END SELECT 
!!!!!!!!! END COMMENTED AREA!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Second, we prepare the file name to read the solution
call MakeBaseNameBGK0D(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(restart_time_txt))//"_sol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")

read (15) time,dt ! record the current time and the used dt
read (15) rkmts ! record the order of the time integrator.
read (15) m ! record the length of the arrays.
! A quick check that the storages are the correct size: 
if (m /= size(f,1)) then 
 print *, "ReadSol_DGV: trying to read a storage of wrong size, m/= size(f,1)"
 close (15)
 stop
end if  
! Next we save the solution .... 
SELECT CASE(rkmts)
	CASE(1)
allocate (frhs1(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable frhs1"
  close(15)
  stop
  end if
read (15) f ! the space for this one needs to be allocated already ... 
    CASE(2)
! We now need to prepare the storage for the data.
allocate (f1(1:m), frhs1(1:m), frhs2(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,frhs1,frhs2"
  close(15)
  stop
  end if
read (15) f,f1,frhs1
    CASE(3)
! We now need to prepare the storage for the data.
allocate (f1(1:m),f2(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,frhs1,frhs2,frhs3"
  close(15)
  stop
  end if
read (15) f,f1,f2,frhs1,frhs2
    CASE(4)
  allocate (f1(1:m),f2(1:m),f3(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,frhs1,frhs2,frhs3,frhs4"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,frhs1,frhs2,frhs3
    CASE(5)
  allocate (f1(1:m),f2(1:m),f3(1:m),f4(1:m),frhs1(1:m),frhs2(1:m),frhs3(1:m),frhs4(1:m),frhs5(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadSol_DGV: Allocation error for variable f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5"
  close(15)
  stop
  end if
read (15) f,f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4
    CASE default 
PRINT *, "Unsupported degree of Multiple Time Stepping method passed"
	 stop
	END SELECT    
! end save 
close (15)
!
end subroutine RestartSolutionBGK0D

!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErr_BGK0D (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
!
! This subroutine will damp on the disk some error arrays
! in the 0D solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErr_BGK0D (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a ! these arrays contain time when the errors are recorder and the errors themselves
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: n ! some counters  
character (len=10) :: FMT1 ! the format string 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time		Ndens		Ubar		Vbar		Wbar"
!prepare the format string to print nnn entries in a row
FMT1="(F18.15 "","" F18.5 "","" F18.5 "","" F18.5 "","" F18.5 "","" F18.5)" 
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.7, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15)" ) time_a(n), "," , ndens_a(n),",", &
              ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n)
end do 
close (15)
end subroutine WriteErr_BGK0D 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrPlus_BGK0D (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will dump on the disk some error arrays
! in the 0D solution
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrPlus_BGK0D (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! these arrays contain time when the errors are recorder and the errors themselves
!
character (len=132) :: file_name ! the variable to store the file name
integer (I4B) :: n ! some counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array
!
! now we open the file in formatted write regime for record and save some stuff in it:
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time, Ndens, Ubar, Vbar, Wbar, T, Tx, Ty, Tz"
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.14,A2,F18.15,A2,F18.15,A2,F18.15,A2,F18.15,A2,F18.15,A2,F18.15,A2, F18.15,  A2, F18.15)" ) &
        time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",", tempr_u_a(n),",", tempr_v_a(n),",", tempr_w_a(n)
end do
close (15)
end subroutine WriteErrPlus_BGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrPlusPlus_BGK0D (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will dump on the disk some error arrays
! the solution
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrPlusPlus_BGK0D (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a, &
                                 mom3_u_a,mom3_v_a,mom3_w_a,mom4_u_a,mom4_v_a,mom4_w_a,mom5_u_a,mom5_v_a,&
                                 mom5_w_a,mom6_u_a,mom6_v_a,mom6_w_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! these arrays contain time when the errors are recorder and the errors themselves
real (DP), dimension (:), intent (in) :: mom4_u_a,mom4_v_a,mom4_w_a,mom6_u_a,mom6_v_a,mom6_w_a ! These arrays include the 4th and 6th moments
real (DP), dimension (:), intent (in) :: mom3_u_a,mom3_v_a,mom3_w_a,mom5_u_a,mom5_v_a,mom5_w_a ! These arrays include the 3th and 5th moments
!
character (len=132) :: file_name ! the variable to store the file name
character (len=250) :: format_line 
integer (I4B) :: n ! some counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameBGK0D2(file_name, "moments/")
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array
!
! now we open the file in formatted write regime for record and save some stuff in it:
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time, Ndens, Ubar, Vbar, Wbar, T, Tx, Ty, Tz, ", & 
             "mom3u, mom3v, mom3w, mom4u, mom4v, mom4w, mom5u, mom5v, mom5w, mom6u, mom6v, mom6w"
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
format_line = "(F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2," // &
                "F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2," // &
                "F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10)" 
write(15,format_line)&
        time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",",tempr_u_a(n),",",tempr_v_a(n),",",tempr_w_a(n),",",mom3_u_a(n),",",mom3_v_a(n),",",mom3_w_a(n),",", &
        mom4_u_a(n),",",mom4_v_a(n),",",mom4_w_a(n),",",mom5_u_a(n),",",mom5_v_a(n),",",mom5_w_a(n),",",mom6_u_a(n),&
        ",",mom6_v_a(n),",",mom6_w_a(n)
end do
close (15)
end subroutine WriteErrPlusPlus_BGK0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteDet
!
! writes the determinate values of the Mom matrix for debugging
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteDet(Det)
use BGKVDCF0D_commvar, only: time

real (DP), intent(in) :: Det
intrinsic Trim

character (len=132) :: file_name ! the variable to store the file name
file_name = "Det.txt"

open (7, file=file_name, position ="APPEND", action="WRITE", form="FORMATTED", status="OLD")
	write(7,"(F18.15,A2,F18.10)") time, ",", Det
close(7)
end subroutine WriteDet

end module BGKVDCF0D_readwrite