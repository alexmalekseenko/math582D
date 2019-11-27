!
!  readwrite.f90
!
! Alex 09/12/14
!
! This modle contains useful procedures for reading and writing from the disk
! The subroutines will be used in other functions of the Discontinuous Galerkin Velocity Discretization Library
! Many of the subroutines are hard linked to the variables of DGV_commvar.f90
! however, some subroutines can be used independently.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!1

module DGV_readwrite
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetDGVParams(pfname,slt)
!
! slt == paramters determining if the print out is generated. If s=0 then there is a printout, 
! if slt is any other number, then no
!
! This subroutine reads the variables from the 
! parameter file on the hard drive. The variables 
! from the common variabl block are accessed directly
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetDGVParams(pfname,slt)
use DGV_commvar
intrinsic Real,Index, Len, Scan, Trim

character (len=*), intent (in) :: pfname ! Name of the file where the parameters are stored 
integer, intent (in) :: slt ! parameter detemining if the prinout is generated. If == 0 then printout generated. 

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
   case ("left endpoint in u")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     u_L=0.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "u_L=", u_L
    end if  
   case ("right endpoint in u")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     u_R=r_bulk(1) ! all other input values are ignored 
    else 
     u_R=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "u_R=", u_R
    end if 
    case ("left endpoint in v")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     v_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     v_L=0.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "v_L=", v_L
    end if 
   case ("right endpoint in v")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     v_R=r_bulk(1) ! all other input values are ignored 
    else 
     v_R=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "v_R=", v_R
    end if 
   case ("left endpoint in w")         ! ready to set up left endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(0,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     w_L=r_bulk(1) ! all other inputs values are ignored 
    else 
     w_L=0.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "w_L=", w_L
    end if 
   case ("right endpoint in w")         ! ready to set up right endpoint in u
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line,r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     w_R=r_bulk(1) ! all other input values are ignored 
    else 
     w_R=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "w_R=", w_R
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
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
   case ("uniform mesh in v")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_v_uniform = .TRUE.
    else
     mesh_v_uniform = .FALSE.
    end if
    if (slt==0) then 
     print *, "mesh_v_uniform=", mesh_v_uniform 
    end if 
   case ("uniform mesh in w")         ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     mesh_w_uniform = .TRUE.
    else
     mesh_w_uniform = .FALSE.
    end if
    if (slt==0) then 
     print *, "mesh_w_uniform=", mesh_w_uniform 
    end if 
   case ("use secondary velocity grid")        ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     SecondaryVelGridInUse = .TRUE.
    else
     SecondaryVelGridInUse = .FALSE.
    end if
    if (slt==0) then 
     print *, "SecondaryVelGridInUse=", SecondaryVelGridInUse 
    end if 
   case ("number of cells in u")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mu_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mu_list)"
     stop
     end if
     !
     Mu_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
     print *, "Mu_list=", Mu_list
    end if  
   case ("degree of local Lagrange basis in u")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (su_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (su_list)"
     stop
     end if
     !
     su_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
     print *, "su_list=", su_list 
    end if
   case ("number of cells in v")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mv_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mv_list)"
     stop
     end if
     !
     Mv_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
     print *, "Mv_list=", Mv_list
    end if  
   case ("degree of local Lagrange basis in v")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (sv_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (sv_list)"
     stop
     end if
     !
     sv_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
     print *, "sv_list=", sv_list 
    end if 
   case ("number of cells in w")           ! ready to set up the number of cells in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 2) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (Mw_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (Mw_list)"
     stop
     end if
     !
     Mw_list=i_bulk(1:m_count)
    end if
    if (slt==0) then   
     print *, "Mw_list=", Mw_list
    end if 
   case ("degree of local Lagrange basis in w")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (sw_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (sw_list)"
     stop
     end if
     !
     sw_list=i_bulk(1:m_count)
    end if
    if (slt==0) then 
     print *, "sw_list=", sw_list 
    end if 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
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
   case ("current A operator base name")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_Aoperator_base_name = "noname"
    else
     current_Aoperator_base_name = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
     print *, "current_Aoperator_base_name=", current_Aoperator_base_name
    end if 
   case ("secondary operator A base name")         ! ready to set name of the directory to store soltuion and other files 
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_Aoperator_base_nameII = "noname"
    else
     current_Aoperator_base_nameII = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
     print *, "current_Aoperator_base_nameII=", current_Aoperator_base_nameII  
    end if
case ("current AKor operator base name")         ! ready to set name of the directory to store soltuion and other files  !!modified
    !!! We read the parameter from the input line 
    if ((line == "") .or. (line == " ")) then 
     current_AKoroperator_base_name = "noname"
    else
     current_AKoroperator_base_name = Trim(Adjustl(line))
    end if 
    if (slt==0) then 
     print *, "current_AKoroperator_base_name=", current_AKoroperator_base_name    
    end if 
case ("num of chunks for Aarray")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    numchnks=i_bulk(1) ! all other inputs values are ignored 
    else 
    numchnks=0 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "numchnks=", numchnks 
    end if 
case ("num of chunks for sec Aarray")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    numchnksII=i_bulk(1) ! all other inputs values are ignored 
    else 
    numchnksII=0 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "numchnksII=", numchnksII 
    end if     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
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
   case ("type of nonuniform mesh in v")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    v_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    v_nonuniform_mesh_type=1 ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "v_nonuniform_mesh_type=", v_nonuniform_mesh_type
    end if
   case ("type of nonuniform mesh in w")  ! ready to set up the type of nonuniform mesh in u
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    w_nonuniform_mesh_type=i_bulk(1) ! all other inputs values are ignored 
    else 
    w_nonuniform_mesh_type=1 ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "w_nonuniform_mesh_type=", w_nonuniform_mesh_type
    end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
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
    ! we are ready to read the ordinary gas constant, reference teperature and viscosity and alpha constant 
    ! for the experiment --- will need in the collision term
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case ("cutoff radius of A-array")       ! ready to set up the cutoff radius for A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     Trad=r_bulk(1) ! all other input values are ignored 
    else 
     Trad=100.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "Trad=", Trad 
    end if
case ("error in integral in Chi")       ! ready to set up the error of integral evaluation in Chi in A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     ErrChi=r_bulk(1) ! all other input values are ignored 
    else 
     ErrChi=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "ErrChi=", ErrChi 
    end if 
case ("error in integral in Epsilon")       ! ready to set up the error in integral in Epsilon for evaluation of A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     ErrEps=r_bulk(1) ! all other input values are ignored 
    else 
     ErrEps=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "ErrEps=", ErrEps 
    end if  
case ("cutoff values for A-array")       ! ready to set up the cutoff values of A-array
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     min_sens=r_bulk(1) ! all other input values are ignored 
    else 
     min_sens=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "min_sens=", min_sens
    end if  
case ("list of basis functions for A-array")           ! ready to set up the list of basis functions for evaluation of A-array 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (I1_list(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (I1_list)"
     stop
     end if
     !
     I1_list=i_bulk(1:m_count)
    end if 
    if (slt==0) then  
     print *, "I1_list=", I1_list 
    end if  
 case ("range of node numbers I2")           ! ready to set up the range in the velocity nodes in the second outside loop 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     I2_from =i_bulk(1)
     if (m_count < 2) then
      I2_to =i_bulk(1)
     else 
      I2_to =i_bulk(2) 
     end if
     !
    end if 
    if (slt==0) then  
     print *, "I2_from=", I2_from, "I2_to=", I2_to 
    end if   
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
case ("error maxwell ESBGK")       ! ready to set up the cut off for linearization
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then
     ES_lev=r_bulk(1) ! all other input values are ignored
    else
     ES_lev=1.0_DP ! set up to the default value
    end if
    if (slt==0) then 
     print *, "ES_lev=", ES_lev
    end if  
case ("error maxwell Vel ESBGK")       ! ready to set up the cut off for linearization
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then
     vel_lev=r_bulk(1) ! all other input values are ignored
    else
     vel_lev=1.0_DP ! set up to the default value
    end if
    if (slt==0) then 
     print *, "vel_lev=", vel_lev
    end if  
case ("error maxwell linearization")       ! ready to set up the cut off for linearization
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then
     linear_lev=r_bulk(1) ! all other input values are ignored
    else
     linear_lev=1.0_DP ! set up to the default value
    end if
    if (slt==0) then 
     print *, "linear_lev=", linear_lev
    end if 
case ("error maxwell decomposition")       ! ready to set up the cutoff valiue for non-linear decomposition
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     decomp_lev=r_bulk(1) ! all other input values are ignored 
    else 
     decomp_lev=100.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "decomp_lev=", decomp_lev 
    end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case ("ref termal velocity")       ! ready to set up the reference thermal velocity unsed in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     C_inf=r_bulk(1) ! all other input values are ignored 
    else 
     C_inf=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "C_inf=", C_inf
    end if 
case ("ref characteristic length")       ! ready to set up the reference characteristic length used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     L_inf=r_bulk(1) ! all other input values are ignored 
    else 
     L_inf=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "L_inf=", L_inf
    end if      
    if ((C_inf > 0) .and. (L_inf > 0)) then
     T_inf=L_inf/C_inf 
     if (slt==0) then 
      print *, "T_inf=", T_inf
     end if 
    else 
     T_inf=1.0_DP
    end if
case ("ref number of molecules")       ! ready to set up the reference total number of molecules in the volume used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     N_inf=r_bulk(1) ! all other input values are ignored 
    else 
     N_inf=1.0_DP ! set up to the default value 
    end if
    if (slt==0) then   
     print *, "N_inf=", N_inf     
    end if  
case ("ref molecular diameter")       ! ready to set up the molecular diameter for hard scpere model used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     mol_diam=r_bulk(1) ! all other input values are ignored 
    else 
     mol_diam=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "mol_diam=", mol_diam
    end if  
case ("ref molecular mass")       ! ready to set up the molecular mass used with hard scpere model used in calulating viscosity
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     mol_mass=r_bulk(1) ! all other input values are ignored 
    else 
     mol_mass=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "mol_mass=", mol_mass
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case ("gas alpha")       ! ready to set up the molecular diameter for hard scpere model used in dimensionless reduction
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     gas_alpha=r_bulk(1) ! all other input values are ignored 
    else 
     gas_alpha=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "gas_alpha=", gas_alpha
    end if  
case ("gas viscosity")       ! ready to set up the dimensional viscosity of the gas 
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     gas_viscosity=r_bulk(1) ! all other input values are ignored 
    else 
     gas_viscosity=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "gas_viscosity=", gas_viscosity
    end if  
case ("Temperature reference")       ! ready to set up reference temperature for the viscosity law
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     gas_T_reference=r_bulk(1) ! all other input values are ignored 
    else 
     gas_T_reference=1.0_DP ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "gas_T_reference=", gas_T_reference
    end if 
case ("ESBGK modifier")       ! ready to set up the parameter alpha that sets the Prandtl number in the ES-BGK model. Pr=1/(1-alpha)
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(100,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     alpha=r_bulk(1) ! all other input values are ignored 
    else 
     alpha=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "alpha=", alpha
    end if 
 case ("use VSS diam in relax rates")        ! ready to set up mesh in u is uniform parameter 
    !!! We read the parameter from the input line 
    if ((trim(Adjustl(line)) == "YES") .or. (trim(Adjustl(line)) == "yes") .or. (trim(Adjustl(line)) == "Yes")) then 
     flag_use_VSS_diam = .TRUE.
    else
     flag_use_VSS_diam = .FALSE.
    end if
    if (slt==0) then 
     print *, "flag_use_VSS_diam=", flag_use_VSS_diam 
    end if 
case ("N of enforced moments")  ! ready to set up the number solution recordings
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 10) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    Order=i_bulk(1) ! all other inputs values are ignored 
    else 
    Order=1 ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "Higherst Moment Order Conserved=", Order
    end if  
case ("N of coefficients in VDCF")  ! ready to set up the number solution recordings
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 10) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    Order_nu=i_bulk(1) ! all other inputs values are ignored 
    else 
    Order_nu=1 ! set up to the default value 
    end if 
    if (slt==0) then  
     print *, "Higherst Moment Order Conserved=", Order_nu    
    end if
case ("mean free time rates update")       ! ready to set up the parameter mft_coeff that sets the time interval to update relaxation reates for enforced moments
    !!! First we read the paramters from the input line
    call ReadRealsFromLine(line_head, line, r_bulk, m_count, Real(1.00,DP)) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
     mft_coeff=r_bulk(1) ! all other input values are ignored 
    else 
     mft_coeff=1.0_DP ! set up to the default value 
    end if  
    if (slt==0) then 
     print *, "mft_coeff=", mft_coeff
    end if
case ("parameters of Korobov net")           ! ready to set up order in variable u 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count >=7) then 
     korob_net_param = i_bulk(1:7)
    else 
     korob_net_param = (/ 1536,341,1081,1517,1201,965,361  /)
     print *, "Set1DHOVParameters: Insufficent number of paramters is given for the korobov net (korob_net_param)"
    end if 
    if (slt==0) then  
     print *, "korob_net_param=", korob_net_param 
    end if
case ("number of nets chunks")           ! ready to set up number of nets with corresponding chunks
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line,i_bulk, m_count, 0) ! last parameter is the default value for the paramenter to be read 
    !!! Now let us process the paramters we have just read !!! 
    if (m_count > 0) then 
     allocate (KorNetsChunks(1:m_count), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVParameters: Allocation error for variable (KorNetsChunks)"
     stop
     end if
     !
     KorNetsChunks=i_bulk(1:m_count)
    end if  
    if (slt==0) then
     print *, "KorNetsChunks=", KorNetsChunks   
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! MPI - related parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case ("number of procs for linear problem")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_lin_proc=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_lin_proc=1 ! set up to the default value 
    end if 
    if (num_lin_proc < 1) then
     num_lin_proc = 1 
     if (slt==0) then 
      print *,  "num_lin_proc=", num_lin_proc, "Specified value is too small. Reset to 1."
     end if
    end if 
    if (slt==0) then 
     print *, "num_lin_proc=", num_lin_proc 
    end if 
case ("no update linearization steps wait")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_linstep_wait=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_linstep_wait=1 ! set up to the default value 
    end if 
    if (num_linstep_wait < 1) then
     num_linstep_wait =1 
     if (slt==0) then 
      print *,  "num_linstep_wait=", num_linstep_wait, "Specified value is too small. Reset to 1."
     end if
    end if 
    if (slt==0) then 
     print *, "num_linstep_wait=", num_linstep_wait
    end if 
   case ("number of mpi copies of A")  ! ready to set up the number of chunks in A arrays 
    !!! First we read the parameters from the input line
    call ReadIntegersFromLine(line_head, line, i_bulk, m_count, 1) ! last parameter is the default value for the paramenter to be read 
    if (m_count > 0) then 
    num_Acopies=i_bulk(1) ! all other inputs values are ignored 
    else 
    num_Acopies=0 ! set up to the default value 
    end if 
    if (slt==0) then 
     print *, "num_Acopies=", num_Acopies
    end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
case default
   if (slt==0) then  
    print *, "Can not process:" //line_head // "="// line
   end if   
  end select
 else
  if (slt==0) then  
   print *, line
  end if  
 end if   
 end do 
close(15)
end subroutine SetDGVParams 

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteGridsDGV
! 
! This subroutine writes the grids_arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteGridsDGV

use DGV_commvar, only: grids_u, grids_v, grids_w, grids_cap_u, grids_cap_v, grids_cap_w
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_grids.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(grids_cap_u,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) grids_cap_u, grids_cap_v, grids_cap_w
!
mm=size(grids_u,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) grids_u
!
mm=size(grids_v,1)                   
!
write (15) mm
! now goes the error itself
write(15) grids_v
!
mm=size(grids_w,1)                   
!
write (15) mm
! now goes the error itself
write(15) grids_w
!
close (15)
!
end subroutine WriteGridsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteCellsDGV
! 
! This subroutine writes the cells_arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteCellsDGV

use DGV_commvar, only: cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov
!
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_cells.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(cells_pgrid,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15) cells_pgrid, cells_cgrid, cells_lu, cells_lv, cells_lw, & 
                   cells_ru, cells_rv, cells_rw, cells_refu, cells_refv, & 
                   cells_refw, cells_gow, cells_gou, cells_gov

!
close (15)
!
end subroutine WriteCellsDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteNodesDGV
! 
! This subroutine writes the nodes arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteNodesDGV

use DGV_commvar, only: nodes_u, nodes_v, nodes_w, nodes_gwts, & 
                   nodes_pcell, nodes_ui, nodes_vi, nodes_wi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_nodes.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(nodes_pcell,1)                   
! we record the size of the arrays
write (15) mm
! now goes the error itself
write(15)  nodes_u,nodes_v,nodes_w,nodes_gwts,&
           nodes_pcell,nodes_ui,nodes_vi, nodes_wi
!
close (15)
!
end subroutine WriteNodesDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysDGV
! 
! This subroutine reads A-arrays from the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysDGV

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
real :: zz ! some scrap variable

! A quick check if the Arrays are already allocated
if (size(A,1)>0) then
   print *,"ReadAarraysDGV: A arrays are already allocated. Exit."
   stop
end if    

! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperDGV(file_name)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
end subroutine ReadAarraysDGV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysDGVII
! 
! This subroutine reads A-arrays for secondary velocity mesh from the disk
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysDGVII

use DGV_commvar, only: A=>AII, A_capphi=>A_capphiII,&
             A_xi=>A_xiII, A_xi1=>A_xi1II, A_phi=>A_phiII
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
real :: zz ! some scrap variable

! A quick check if the Arrays are already allocated
if (size(A,1)>0) then
   print *,"ReadAarraysDGVII: AII arrays are already allocated. Exit."
   stop
end if    

! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperDGVII(file_name)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
end subroutine ReadAarraysDGVII
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysChnksDGV
! 
! This subroutine reads A-arrays from the disk when A-arrays are stored on the hard drive in chunks, 
! it assembled them 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysChnksDGV(nchnks)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
use DGV_dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), intent (in) :: nchnks ! the number of chunks of data  
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,An,i ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
!
real (DP), dimension (:), allocatable :: Achnks
integer (I4B), dimension (:), allocatable :: Achnks_capphi, Achnks_xi, Achnks_xi1, Achnks_phi

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysDGV: A arrays are already allocated. Exit."
   stop
end if    
!! Now the chunks of the A-arrays will be read and pieced together in the sequence imbedded in the filennames
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, A_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or A_phi. 
!! 
!! The first chunk will go directly to A-arrays.... 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,0)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
! We have read the first chunk directly into A-arrays. 
!
do i=2,nchnks ! This array will loop until all chunks are read
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,i-1)
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 close (15)

 ! We now need to prepare the storage for the data.
 allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
 close (15)
 !
 ! now we need to record the new chunk into the A-arrays. For that we will need to extend the A-arrays (see below)
 ! first he A_capphi array -- because we expect that I1 do not overlap in chunks and if so, we expect the records be kept 
 ! in the increasing order of I1's --- I1s are stored in A_phi. Because of all that, we can just add Achnks_capphi to A_capphi:
 A_capphi = A_capphi + Achnks_capphi
 ! Now we need to extend the A arrays to add the new records:
 An=size(A,1) ! Let us find out how lond the arrays are now
 call ExtendAarraysDGV(An,nn)  ! we make room for the new records, An is the new size 
 ! Now we add the new reords
 A(An-nn+1:An)=Achnks
 A_xi(An-nn+1:An)=Achnks_xi
 A_xi1(An-nn+1:An)=Achnks_xi1
 A_phi(An-nn+1:An)=Achnks_phi
 ! now we need to destroy the chunks arrays:
 deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)
 ! ready to reead the next chunk
end do 

end subroutine ReadAarraysChnksDGV
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysChnksDGVII
! 
! This is a copy of the above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! This subroutine reads A-arrays from the disk when A-arrays are stored on the hard drive in chunks, 
! it assembled them 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysChnksDGVII(nchnks)

use DGV_commvar, only: A=>AII,A_capphi=>A_capphiII,A_xi=>A_xiII,A_xi1=>A_xi1II,A_phi=>A_phiII
use DGV_dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), intent (in) :: nchnks ! the number of chunks of data  
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,An,i ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
!
real (DP), dimension (:), allocatable :: Achnks
integer (I4B), dimension (:), allocatable :: Achnks_capphi, Achnks_xi, Achnks_xi1, Achnks_phi

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysDGVII: A arrays are already allocated. Exit."
   stop
end if    
!! Now the chunks of the A-arrays will be read and pieced together in the sequence imbedded in the filennames
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, A_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or A_phi. 
!! 
!! The first chunk will go directly to A-arrays.... 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGVII(file_name,0)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGVII: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:nn), A_xi1(1:nn), A(1:nn), A_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGVII: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	A_capphi,A,A_xi,A_xi1,A_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
! We have read the first chunk directly into A-arrays. 
!
do i=2,nchnks ! This array will loop until all chunks are read
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGVII(file_name,i-1)
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 close (15)

 ! We now need to prepare the storage for the data.
 allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGVII: Allocation error for variable  A_capphi"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGVII: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
 close (15)
 !
 ! now we need to record the new chunk into the A-arrays. For that we will need to extend the A-arrays (see below)
 ! first he A_capphi array -- because we expect that I1 do not overlap in chunks and if so, we expect the records be kept 
 ! in the increasing order of I1's --- I1s are stored in A_phi. Because of all that, we can just add Achnks_capphi to A_capphi:
 A_capphi = A_capphi + Achnks_capphi
 ! Now we need to extend the A arrays to add the new records:
 An=size(A,1) ! Let us find out how lond the arrays are now
 call ExtendAarraysDGVII(An,nn)  ! we make room for the new records, An is the new size 
 ! Now we add the new reords
 A(An-nn+1:An)=Achnks
 A_xi(An-nn+1:An)=Achnks_xi
 A_xi1(An-nn+1:An)=Achnks_xi1
 A_phi(An-nn+1:An)=Achnks_phi
 ! now we need to destroy the chunks arrays:
 deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)
 ! ready to reead the next chunk
end do 

end subroutine ReadAarraysChnksDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ScanAarrsChnks4Acapphi(A_capphi_aggreg,Act_aggreg)
! 
! This subroutine reads A-arrays from the disk when A-arrays are stored on the hard drive in chunks, 
! It determines the number of records of A arrays stored in each chunk. This information is recorded in 
! the array chnks_Act
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ScanAarrsChnks4Acapphi(chunks_Act)

use DGV_dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), dimension (:), intent (out) :: chunks_Act
!

character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,i ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nchnks
!
integer (I4B), dimension (:), allocatable :: Achnks_capphi

nchnks=size(chunks_Act,1) ! determine the number of chunks

!! Now the chunks of the A-arrays will be read and their A_capphi pieced together
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,0)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)

! We now need to prepare the storage for the data.
allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ScanAarrsChnks4Acapphi: Allocation error for variables Achnks_capphi"
  stop
  end if

! now we open the file, again, to populate the A_capphi arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	Achnks_capphi ! read the array		A_capphi (other arrays should be still there : A,A_xi,A_xi1,A_phi)
close (15)
!
! We have scanned the first chunk of A-arrays and read the A_capphi for the first chunk. Now we need to move it to A_capphi_aggreg
chunks_Act(1)= sum(Achnks_capphi)
! 
! Now we need to proceed with other chuncks. 

do i=2,nchnks ! This array will loop until all chunks are read
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,i-1)
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi ! read the arrays A_capphi					
 close (15)
 !
 ! now we need to add the number of records to A_capphi . 
 chunks_Act(i) = sum(Achnks_capphi)
 ! ready to reead the next chunk
end do 
deallocate (Achnks_capphi)
!
end subroutine ScanAarrsChnks4Acapphi
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadAarraysChnksOfstNrecDGV
! 
! This subroutine reads a given number of records of A-arrays from the disk when A-arrays are stored on 
! the hard drive in chunks. This subroutine will open chunk with the number fstchnk, it will skip the first 
! (ofst) records, and reads exactly numrec records into A arrays. If the chunk is finished but not eanough records have been read, 
! the software moves to the next chunk untill the desired number of records achieved or untill all chunks were read. 
!
! fstchnk -- the number of the first chunk 
! ofst  -- the number of records to skip in the first chunk (presumably these records are read at another process.
! numrec -- the total number of records to read into A-arrays, 
! nchnks -- the total number of chunks of A.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadAarraysChnksOfstNrecDGV(fstchnk,ofst,numrec,nchnks)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
use DGV_dgvtools_mod
!                   
intrinsic Trim
!
integer (I4B), intent (in) :: fstchnk ! the number of first chunks of data  
integer (I4B), intent (in) :: ofst ! the number of records to skip in the first chunk   
integer (I4B), intent (in) :: numrec ! the number of records that needs to be read  
integer (I4B), intent (in) :: nchnks ! the total number of chunks in A  
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn,i,j,Act, Act_need ! is the scrap variable to keep te size of the A-arrays
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
!
real (DP), dimension (:), allocatable :: Achnks
integer (I4B), dimension (:), allocatable :: Achnks_capphi, Achnks_xi, Achnks_xi1, Achnks_phi

! A quck check if the Arrays are already allocated
if (size(A,1)>0) then 
   print *,"ReadAarraysChnksOfstNrecDGV: A arrays are already allocated. Exit."
   stop
end if    
!! Now the chunks of the A-arrays will be read and pieced together in the sequence imbedded in the filennames
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, A_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or A_phi. 
!! 
!! The first chunk will go directly to A-arrays.... 
!! 
! first, we prepare the file name to read the A-arrays
call MakeBaseNameAoperChnksDGV(file_name,fstchnk-1)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
close (15)
! check if the read was successful:
if (code_line/=0) then 
 print *,  "ReadAarraysChnksOfstNrecDGV: read error for chunk of A-array. Possibly wrong name or damaged file. Stop"
 stop
end if

! We now need to prepare the storage for the data.
allocate (A_capphi(1:mm), Achnks_capphi(1:mm), stat = loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_capphi"
  stop
  end if
allocate (A_xi(1:numrec), A_xi1(1:numrec), A(1:numrec), A_phi(1:numrec), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysChnksOfstNrecDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
! now we open the file, again, to populate the A-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
close (15)
!
! We have read the first chunk directly into Achnks-arrays.
! Now we skip the number of records ofst and put the rest in the permanent A-Arrays.
! first, determine the number of records to read: 
A_capphi = 0
Act = 0
Act_need = min(nn-ofst,numrec-Act)  ! this variable keeps track of how many records needs to be added to the array from this chunk of A 
A(1:Act_need)=Achnks(ofst+1:ofst+Act_need)
A_xi(1:Act_need)=Achnks_xi(ofst+1:ofst+Act_need)
A_xi1(1:Act_need)=Achnks_xi1(ofst+1:ofst+Act_need)
A_phi(1:Act_need)=Achnks_phi(ofst+1:ofst+Act_need)
Act=Act_need ! now we remember how many records we added to the A arrays... 
!
!  Now we need to write an A_capphi that corresponds to what is in A
i=1
do while ((ofst > sum(Achnks_capphi(1:i))) .and. (i<= mm)) ! we need to skip ofst records in this array.
 i=i+1
end do                            ! when this difference is negative, we went over the offset number of recods
if (i>mm) then                    ! a quick check if we had enough records 
 print *, "ReadAarraysChnksOfstNrecDGV: may be rerror in the work array -- no records in the chunk after offset taken"
 stop
end if
                                 
if (sum(Achnks_capphi(1:i)) - ofst < Act_need) then            ! now we need to see if the records for the basis function i is longer than the total number of 
 A_capphi(i) = sum(Achnks_capphi(1:i)) - ofst            !  records that will be read from this chunk. The appropriate number of records is acconted for.
 Act_need = Act_need - sum(Achnks_capphi(1:i)) + ofst
else 
 A_capphi(i) = Act_need  
 Act_need = 0
end if 
 
do while ((Act_need > 0) .and. (i<mm))
 i=i+1
 if (Achnks_capphi(i) >= Act_need) then 
   A_capphi(i) = Act_need
   Act_need = 0
 else 
   A_capphi(i) = Achnks_capphi(i) 
   Act_need = Act_need - Achnks_capphi(i)
 end if   
end do 
if ((i>=mm).and.(Act_need >0)) then                    ! a quick check if we had enough records 
 print *, "ReadAarraysChnksOfstNrecDGV: may be error in the work array -- no records in the chunk after offset taken"
 stop
end if     
deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)                            
! We finished dealing with the first chunk.

! Now if we still have insufficent records in A, we will try to read the next chunk and us records from there..
j = fstchnk+1 ! point the index to the next chunk (remember that chunks are numbered from zero) 

do while ((Act < numrec) .and. (j <= nchnks)) ! This array will loop until enough records is read or until all chunks are read 
 ! first, we prepare the file name to read the A-arrays
 call MakeBaseNameAoperChnksDGV(file_name,j-1) ! the names of chunks are numbered from zero, therefore it is j-1
 file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15,  iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 close (15)
 
 ! We now need to prepare the storage for the data.
 allocate (Achnks_capphi(1:mm), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_capphi"
  stop
  end if
 allocate (Achnks_xi(1:nn), Achnks_xi1(1:nn), Achnks(1:nn), Achnks_phi(1:nn), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadAarraysDGV: Allocation error for variable  A_xi, A_xi1, A_phi or A"
  stop
  end if
 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn ! notice that we read using one read statement because this is 
								! how it was written
 read (15, iostat=code_line)	Achnks_capphi,Achnks,Achnks_xi,Achnks_xi1,Achnks_phi ! read the arrays		A_capphi,A,A_xi,A_xi1,A_phi					
 close (15)
 !
 ! now we need to record the new chunk into the A-arrays. 
 ! We need to put only enough records to meet the total number of records requirement
 Act_need = min(nn, numrec - Act)  ! this calculates the number of records to be saved fron this chunk.  There is no offset aftert the first chunk
 A(Act+1:Act+Act_need) = Achnks(1:Act_need)
 A_xi(Act+1:Act+Act_need) = Achnks_xi(1:Act_need)
 A_xi1(Act+1:Act+Act_need) = Achnks_xi1(1:Act_need)
 A_phi(Act+1:Act+Act_need) = Achnks_phi(1:Act_need)
 Act=Act+Act_need ! now we remember how many records have been added to the A arrays so far
 ! No wwe need to update the A_caphi array. 
 i=0              ! reset the counter 
 do while ((Act_need > 0) .and. (i<mm)) ! will loop untill all records were taken account of or while ran out of records
  i=i+1
  if (Achnks_capphi(i) >= Act_need) then 
   A_capphi(i) = A_capphi(i) + Act_need 
   Act_need = 0
  else 
   A_capphi(i) = A_capphi(i) + Achnks_capphi(i) 
   Act_need = Act_need - Achnks_capphi(i)
  end if   
 end do
 if ((i>=mm).and.(Act_need >0)) then                    ! a quick check if we had enough records 
  print *, "ReadAarraysChnksOfstNrecDGV: may be rerror in the work array -- no records in the chunk after offset taken"
  stop
 end if     
 ! now we need to destroy the chunks arrays:
 deallocate (Achnks,Achnks_xi,Achnks_xi1,Achnks_phi,Achnks_capphi)
 ! ready to reead the next chunk
 j=j+1 ! Advance the number of chunk by one
end do 
! Next we need to set up the supplementary "shift" arrays that will be used for integration of the right side...... 
end subroutine ReadAarraysChnksOfstNrecDGV
!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteAarraysDGV
! 
! This subroutine writes the A arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteAarraysDGV

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_Aarrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(A_capphi,1)
nn=size(A,1)                   
! we record the size of the arrays
write (15) mm,nn
! now goes the error itself
write(15)  A_capphi,A,A_xi,A_xi1,A_phi
!                 
close (15)
!
end subroutine WriteAarraysDGV 
 
                   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameDGV 
!
! This function creates a recognizable name for the file.
! a lot of parameter values is cited in the file name. These are parametes 
! related to the DG velocity discretization. Thus the file name will contain 
! values of important parameters of the DG discretization. It will therefore be easier to 
! recall what numerical simualtion they belong to. Still 
! many parametersa re not represented -- it is advised to keep a copy of the the parameter file 
! with the simulation results. ot all parameters are 
! present in the name. 
!
!The subroutine produced the potion of the file name that 
!will be common for the numericl simulation and will be 
! be used by many output subroutines
!
! DO not detatch from the main program. 
! 
! Do not call until parameters su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          Mu,Mv,Mw,curr_time,u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the DGV_commvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameDGV (file_name)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_solution_base_name, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_solution_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameDGV
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameDGV2
!
! Modified by Craig 12/31/12
! This function creates a recognizable name for the file.
! a lot of parameters is cited in the file name. 
! This function will set up the basic
! file name that will be common for the numerical simulation  
! and will be used by many output procedures.
!
! DO not detatch from the main program.
!
! Do not call until parameters su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          N,Mu,Mv,Mw,curr_time, u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the DGV_commvar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameDGV2 (file_name, directory)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_solution_base_name, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
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
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameDGV2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperDGV 
!
! This function creates a recognizable name for the file to store A-operator in one file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          Mu,Mv,Mw,curr_time, u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperDGV (file_name)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_Aoperator_base_name, &
          u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameAoperDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperDGVII
!
! This function creates a recognizable name for the file to store A-operator in one file.
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! This is a copy of an above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters suII,svII,swII, current_solution_dir, current_solution_base_name,
!          MuII,MvII,MwII, u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperDGVII (file_name)

use DGV_commvar, only: su=>suII,sv=>svII,sw=>swII,Mu=>MuII,Mv=>MvII,Mw=>MwII,&
          current_solution_dir, current_Aoperator_base_name=>current_Aoperator_base_nameII, &
          u_nonuniform_mesh_type, v_nonuniform_mesh_type, w_nonuniform_mesh_type, &
          mesh_u_uniform, mesh_v_uniform, mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
file_name = trim(current_solution_dir)//trim(current_Aoperator_base_name)
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameAoperDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperChnksDGV 
!
! This function creates a recognizable name for the files to store chunks of operator A. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          Mu,Mv,Mw,curr_time,u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperChnksDGV (file_name,i)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_Aoperator_base_name, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
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
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameAoperChnksDGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAoperChnksDGVII
!
! This function creates a recognizable name for the files to store chunks of operator A. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! This is a copy of an above subroutine, with the only change that it works with secondary mehses. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
! DO not detouch from the main program. 
! 
! Do not call untill parameters suII,svII,swII, current_solution_dir, current_solution_base_name,&
!          MuII,MvII,MwII,u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAoperChnksDGVII (file_name,i)

use DGV_commvar, only: su=>suII,sv=>svII,sw=>swII,Mu=>MuII,Mv=>MvII,Mw=>MwII, &
          current_solution_dir, current_Aoperator_base_name=>current_Aoperator_base_nameII, &
          u_nonuniform_mesh_type,v_nonuniform_mesh_type,w_nonuniform_mesh_type, &
          mesh_u_uniform,mesh_v_uniform,mesh_w_uniform
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
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameAoperChnksDGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteDGV2dArray(A,suffix)
!
! This subroutine dumps on the hard drive the 2D array. 
! A is the array of doubles 
! suffix is a the string containing a combination of letters that 
! will be attached to the name fo the file... 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteDGV2dArray(A,suffix)
!                   
intrinsic Trim
!
real (DP), dimension (:,:), intent (in) :: A ! the incoming array of doubles
character (len=10) :: suffix ! the string to hold the suffix
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: m1,m2 ! scrap variables to keep te size of the 2d arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_"//trim(Adjustl(suffix))//".dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
m1=size(A,1)
m2=size(A,2)                   
! we record the sizes of the array
write (15) m1,m2
! now goes the array itself
write(15)  A
!
close (15)
!
end subroutine WriteDGV2dArray

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteI1DGV
!
! This subroutine will write a sequence of node indices to be later used in the evaluation of A operator. 
! 
! I1 is the array of integers that needs to be written on the disk. Notice that the first element in the array 
!  I1(1) contains the total number of records. So then only 2..I1(1)+1 cells in the array are importnat   
!
! u,v,w are the components of the velocity point. Am active cell in velocity space will be found that containes 
! this velocity point and indices of all nodes form this cell will be listed.  
!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteI1DGV(u,v,w)
!
use DGV_dgvtools_mod
!@
intrinsic Trim
!!!!!
real (DP), intent (in) :: u,v,w !  the components of the velocity
!!!!!!!!!!!!!!!!!!!!!!
!
integer (I4B) :: celli ! the number of the cell that contains point (u,v,w)
character (len=132) :: file_name ! the variable to store the file name 
integer, parameter :: nnn=5 !! Number of entries in the row 
integer (I4B), dimension (3000) :: I1 ! an array to store the found nodes' indices.
integer (I4B) :: j,k,n ! some counters  
character (len=20) :: parchar    ! string to keep char value of the paameter 
character (len=10) :: FMT1 ! the format string 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we find the number of the active cell that contains point (u,v,w)
celli = FindCellContainsPoint_DGV(u,v,w)
! we find the indices of all nodes that belong to cell with the number icell
call FindI1sByCellNum_DGV(celli,I1)
! Now it is time to record this infomation 

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_I1"
write (parchar, "(F4.1)") u
file_name = trim(file_name)//trim(Adjustl(parchar))//"u"
write (parchar, "(F4.1)") v
file_name = trim(file_name)//trim(Adjustl(parchar))//"v"
write (parchar, "(F4.1)") w
file_name = trim(file_name)//trim(Adjustl(parchar))//"w.txt"  ! this file will keep the array 
!
! now we open the file in formatted write regime for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the numbers of nodes belonging to the cell containing velocity point:"
write(15, *) "u=", u, "v=", v, "w=", w
write(15, *) "I1= "
!prepare the format string to print nnn entries in a row
write (parchar, "(I9)") nnn
FMT1="("//trim(Adjustl(parchar))//"(I5, A))" 
! next we print array using nnn entries in a row
n=I1(1)
k=0
do while (k + nnn < n)
  write (15, FMT1) (I1(1+j), "," , j=k+1,k+nnn) 
  k=k+nnn
end do 
write (15,"(3(I5, "","" ))") I1(2+k:1+n)
close (15)
end subroutine WriteI1DGV 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErr_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
!
! This subroutine will damp on the disk some error arrays
! in the 0D solution 
!  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErr_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a)
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
call MakeBaseNameDGV(file_name)
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
end subroutine WriteErr_DGV 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrPlus_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will dump on the disk some error arrays
! in the 0D solution
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrPlus_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
intrinsic Trim
!!!!!
real (DP), dimension (:), intent (in) :: time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a ! these arrays contain time when the errors are recorder and the errors themselves
!
character (len=132) :: file_name ! the variable to store the file name
integer (I4B) :: n ! some counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array
!
! now we open the file in formatted write regime for record and save some stuff in it:
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time		Ndens		Ubar		Vbar		Wbar    T     Tx     Ty    Tz"
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
  write (15,"(F18.14, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15, A2, F18.15,  A2, F18.15)" ) &
        time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",", tempr_u_a(n),",", tempr_v_a(n),",", tempr_w_a(n)
end do
close (15)
end subroutine WriteErrPlus_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteErrPlusPlus_DGV (ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a)
!
! This subroutine will dump on the disk some error arrays
! the solution
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteErrPlusPlus_DGV (time_a,ndens_a,ubar_a,vbar_a,wbar_a,tempr_a,tempr_u_a,tempr_v_a,tempr_w_a, &
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
call MakeBaseNameDGV2(file_name, "moments/")
file_name = trim(Adjustl(file_name))//"_macroerr.txt"  ! this file will keep the array
!
! now we open the file in formatted write regime for record and save some stuff in it:
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, *) "This is the file containing the values of density, average velocity and temperature, and the time stamp:"
write(15, *) "Time, Ndens, Ubar, Vbar, Wbar, T, Tx, Ty, Tz, mom3u, mom3v, mom3w, mom4u, mom4v, mom4w, ", &
             "mom5u,  mom5v, mom5w, mom6u, mom6v, mom6w,"
! next we print array using nnn entries in a row
do n=1,size(ndens_a,1)
format_line = "(F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,"// &
               "F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,"// &
               "F18.10,A2,F18.10,A2,F18.10,A2,F18.10,A2,F18.10)"
write(15,trim(Adjustl(format_line))) time_a(n), "," , ndens_a(n),",", ubar_a(n),",", vbar_a(n),",", wbar_a(n),",", tempr_a(n), &
        ",",tempr_u_a(n),",",tempr_v_a(n),",",tempr_w_a(n),",",mom3_u_a(n),",",mom3_v_a(n),",",mom3_w_a(n),",", &
        mom4_u_a(n),",",mom4_v_a(n),",",mom4_w_a(n),",",mom5_u_a(n),",",mom5_v_a(n),",",mom5_w_a(n),",",mom6_u_a(n),",",mom6_v_a(n),",",mom6_w_a(n)
end do
close (15)
end subroutine WriteErrPlusPlus_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  WriteAKorArraysDGV
! 
!  This subroutine writes the Akor arrays on the disk
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WriteAKorArraysDGV

use DGV_commvar, only: Akor,Akor_capphi, Akor_k, Akor_phi, korob_net_param
!                   
intrinsic Trim
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: mm,nn ! is the scrap variable to keep te size of the cells arrays

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_Akor_arrs.dat"  ! this will keep the cells array of the problem
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
mm=size(Akor_capphi,1)
nn=size(Akor,1)                   
! we record the size of the arrays
write (15) mm,nn,korob_net_param
! now goes the error itself
write(15)  Akor_capphi,Akor,Akor_k,Akor_phi
!                 
close (15)
!
end subroutine WriteAKorArraysDGV 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ReadAKorArraysOneNet_DGV (netnum, numchks)
!
!  This subroutine will read Akor arrays for a single net from the hard drive
! 
!  netnum -- this is the number of the net which will be read
!  nuchks -- the total number of chunks to read. 
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadAKorArraysOneNet_DGV (netnum, numchks)


use  DGV_commvar, only:  AkorAllNets, AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                                        korob_net_paramAllNets  

integer (I4B), intent (in) :: netnum ! the number of the net fo which the Akor is read
integer (I4B), intent (in) :: numchks ! the number of chunk in which Akor array is divided for the given Korobov net
!!!!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: nn,mm,length,j,i ! is the scrap variable to keep te size of the Akor-arrays  
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B), dimension (1:7) :: lockorobparams ! scrap arrray to store parameters of the Korobov net
integer (I4B), dimension (:), allocatable :: scrap_capphi ! scrap arrays to read the chunk capphi array 
!
real (DP), dimension (:), pointer ::  preal       
integer (I4B), dimension (:), pointer :: pint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now the chunks of the Akor-arrays will be read and pieced together in the sequence imbedded in the file names
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, Akor_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or Akor_phi. 
!! The first chunk will go directly to Akor-arrays.... 
!! first, we prepare the file name to read the A_kor-arrays
call MakeBaseNameAKorNetChnksDGV(file_name,netnum-1,0)                                                            
file_name = trim(Adjustl(file_name))//"_Akor_arrs.dat"  ! this will keep the cells array of the problem    
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line)  mm,nn ! nn = size of Akor, mm=  num of basis functions/vel. nodes  

close (15)

! We now need to prepare the storage for the data.
allocate (AkorAllNets_capphi(netnum)%p(1:mm), scrap_capphi(1:mm), korob_net_paramAllNets(netnum)%p(1:7), stat=loc_alloc_stat)         !
    !
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable  AkorAllNets_capphi, netnum=", netnum
  stop
 end if
allocate ( AkorAllNets(netnum)%p(1:nn), AkorAllNets_k(netnum)%p(1:nn),&
       AkorAllNets_phi(netnum)%p(1:nn), stat=loc_alloc_stat)    ! loc_alloc_stat  to keep allocation status
    !
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variables  AkorAllNets,AkorAllNets_k,AkorAllNets_phi. netnum=",  netnum
  stop
 end if
! now we open the file, again, to populate the Akor-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn,korob_net_paramAllNets(netnum)%p       

read (15, iostat=code_line)	AkorAllNets_capphi(netnum)%p, AkorAllNets(netnum)%p,&
 AkorAllNets_k(netnum)%p, AkorAllNets_phi(netnum)%p   ! read the arrays  Akor,Akor_k,Akor_phi, Akor_capphi for the net with the number numnet			
close (15)
!
length = nn ! the total length of the Akor array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We have read the first chunk directly into Akor-arrays.                                             
! Now we will read the rest of the chunks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=2,numchks ! This array will loop until all chunks are read in a given net                                           
 ! first, we prepare the file name to read the A-arrays
 call  MakeBaseNameAKorNetChnksDGV(file_name,netnum-1,j-1)																	 
 file_name = trim(Adjustl(file_name))//"_Akor_arrs.dat"  ! this will keep the cells array of the problem			 
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn,lockorobparams
 close (15)
 ! Need to check the consistency of the Akor array data. We verify that the data in the chunk is
 ! obtained for the same Korobov net as in korob_net_paramAllNets(netnum)%p(1:7)
 !
 do i = 1,7
  if (lockorobparams(i) .NE. korob_net_paramAllNets(netnum)%p(i)) then           
     print *, "ReadAKorArraysOneNet_DGV: Inconsistent parameters for Korobov net in chunks of Akor arrays. netnum=", netnum
     stop
  end if
 end do
 !
 ! We now need to prepare the storage for the data. First, we will deal with the real array Akor
 allocate (preal(length+nn), stat=loc_alloc_stat) 
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable preal"
  stop
 end if
 preal(1:length) = AkorAllNets(netnum)%p  ! save the data in AkorAllNets in preal
 deallocate(AkorAllNets(netnum)%p)    ! deallocate the storage for Akor for this net
 AkorAllNets(netnum)%p => preal       ! make the pointed to point to the place where the data Akor is located.  
 nullify (preal)                                        ! free the pointer to the array oif reals.   

 ! Next we will prepare the extended storage for Akor_k and Akor_phi 
 allocate (pint(length+nn), stat=loc_alloc_stat)  
 if (loc_alloc_stat >0) then 
   print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable pint"
   stop
 end if
 pint(1:length)= AkorAllNets_k(netnum)%p     ! copy the content of Akor_k intot he larger storage 
 deallocate(AkorAllNets_k(netnum)%p)     ! deallocate the old storage
 AkorAllNets_k(netnum)%p => pint         ! make the pointer to point to the storage with more space
 nullify (pint)                                              ! free the pointer
 !
 allocate (pint(length+nn), stat=loc_alloc_stat)  
 if (loc_alloc_stat >0) then 
   print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable pint (II)"
   stop
 end if
 pint(1:length) = AkorAllNets_phi(netnum)%p  ! copy the content of Akor_phi to the larger storage
 deallocate(AkorAllNets_phi(netnum)%p)  ! deallocate the old storage
 AkorAllNets_phi(netnum)%p => pint      ! make the pointer to point to the storage with more space
 nullify (pint)                         ! free the pointer
 !	       

 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm, nn, lockorobparams
 read (15, iostat=code_line)  scrap_capphi,AkorAllNets(netnum)%p(length+1:length+nn),&
                              AkorAllNets_k(netnum)%p(length+1:length+nn),&
                              AkorAllNets_phi(netnum)%p(length+1:length+nn)
				
 close (15)                                                                     
 !!!
 ! we will need to extend the A_capphi array and we expect the records be kept 
 ! in the increasing order
 AkorAllNets_capphi(netnum)%p = AkorAllNets_capphi(netnum)%p + scrap_capphi
length = length+nn 
end do
!
deallocate(scrap_capphi) 
!
end subroutine ReadAKorArraysOneNet_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ReadAKorArraysOneNet_DGVII (netnum, numchks)
!
! ATTENTION: the name of the Akor Array is using primary mesh -- need to fix later
!
! This is a copy of the above subroutine, with the only change that it works with secondary meshes. 
! This is accomplished by re-naming variable in the USE commvar, only: statement. 
!
! The suffix II in the name suggests that the subroutine works with the secondary mesh
!
!  This subroutine will read Akor arrays for a single net from the hard drive
! 
!  netnum -- this is the number of the net which will be read
!  nuchks -- the total number of chunks to read. 
! 
! ATTENTION: the name of the Akor Array is using primary mesh -- need to fix later
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadAKorArraysOneNet_DGVII (netnum, numchks)


use  DGV_commvar, only:  AkorAllNets=>AkorAllNetsII, AkorAllNets_k=>AkorAllNets_kII,      &
           AkorAllNets_phi=> AkorAllNets_phiII, AkorAllNets_capphi=>AkorAllNets_capphiII, &
           korob_net_paramAllNets=>korob_net_paramAllNetsII
                                        			 

integer (I4B), intent (in) :: netnum ! the number of the net fo which the Akor is read
integer (I4B), intent (in) :: numchks ! the number of chunk in which Akor array is divided for the given Korobov net
!!!!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: nn,mm,length,j,i ! is the scrap variable to keep te size of the A_kor-arrays   nn=size(A_kor,1)
integer (I4B) :: code_line ! scrap variable
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B), dimension (1:7) :: lockorobparams ! scrap arrray to store parameters of the Korobov net
integer (I4B), dimension (:), allocatable :: scrap_capphi ! scrap arrays to read the chunk capphi array 
!
real (DP), dimension (:), pointer ::  preal       
integer (I4B), dimension (:), pointer :: pint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Now the chunks of the Akor-arrays will be read and pieced together in the sequence imbedded in the file names
!! IT IS IMPORTNAT THAT THE INFORMATION IN STORED IN THE CHUNKS IN THE RIGHT ORDER. 
!! In particular, Akor_capphi array will list the number of records for nodes. It will be expected that records in A 
!! are ordered in the increasing I1, or Akor_phi. 
!! The first chunk will go directly to Akor-arrays.... 
!! first, we prepare the file name to read the A_kor-arrays
call MakeBaseNameAKorNetChnksDGV(file_name,netnum-1,0)                                                            
file_name = trim(Adjustl(file_name))//"_Akor_arrs.dat"  ! this will keep the cells array of the problem    
!
! now we open the file for reading and read some stuff from it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15,  iostat = code_line)  mm,nn ! nn = size of Akor, mm=  num of basis functions/vel. nodes  

close (15)

! We now need to prepare the storage for the data.
allocate (AkorAllNets_capphi(netnum)%p(1:mm), scrap_capphi(1:mm), korob_net_paramAllNets(netnum)%p(1:7), stat=loc_alloc_stat)         !
    !
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable  AkorAllNets_capphi, netnum=", netnum
  stop
 end if
allocate ( AkorAllNets(netnum)%p(1:nn), AkorAllNets_k(netnum)%p(1:nn),&
       AkorAllNets_phi(netnum)%p(1:nn), stat=loc_alloc_stat)    ! loc_alloc_stat  to keep allocation status
    !
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variables  AkorAllNets,AkorAllNets_k,AkorAllNets_phi. netnum=",  netnum
  stop
 end if
! now we open the file, again, to populate the Akor-arrays
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15, iostat = code_line) mm,nn,korob_net_paramAllNets(netnum)%p       

read (15, iostat=code_line)	AkorAllNets_capphi(netnum)%p, AkorAllNets(netnum)%p,&
 AkorAllNets_k(netnum)%p, AkorAllNets_phi(netnum)%p   ! read the arrays  Akor,Akor_k,Akor_phi, Akor_capphi for the net with the number numnet			
close (15)
!
length = nn ! the total length of the Akor array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! We have read the first chunk directly into Akor-arrays.                                             
! Now we will read the rest of the chunks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do j=2,numchks ! This array will loop until all chunks are read in a given net                                           
 ! first, we prepare the file name to read the A-arrays
 call  MakeBaseNameAKorNetChnksDGV(file_name,netnum-1,j-1)																	 
 file_name = trim(Adjustl(file_name))//"_Akor_arrs.dat"  ! this will keep the cells array of the problem			 
 !
 ! now we open the file for reading and read some stuff from it: 
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm,nn,lockorobparams
 close (15)
 ! Need to check the consistency of the Akor array data. We verify that the data in the chunk is
 ! obtained for the same Korobov net as in korob_net_paramAllNets(netnum)%p(1:7)
 !
 do i = 1,7
  if (lockorobparams(i) .NE. korob_net_paramAllNets(netnum)%p(i)) then           
     print *, "ReadAKorArraysOneNet_DGV: Inconsistent parameters for Korobov net in chunks of Akor arrays. netnum=", netnum
     stop
  end if
 end do
 !
 ! We now need to prepare the storage for the data. First, we will deal with the real array Akor
 allocate (preal(length+nn), stat=loc_alloc_stat) 
 if (loc_alloc_stat >0) then 
  print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable preal"
  stop
 end if
 preal(1:length) = AkorAllNets(netnum)%p  ! save the data in AkorAllNets in preal
 deallocate(AkorAllNets(netnum)%p)    ! deallocate the storage for Akor for this net
 AkorAllNets(netnum)%p => preal       ! make the pointed to point to the place where the data Akor is located.  
 nullify (preal)                                        ! free the pointer to the array oif reals.   

 ! Next we will prepare the extended storage for Akor_k and Akor_phi 
 allocate (pint(length+nn), stat=loc_alloc_stat)  
 if (loc_alloc_stat >0) then 
   print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable pint"
   stop
 end if
 pint(1:length)= AkorAllNets_k(netnum)%p     ! copy the content of Akor_k intot he larger storage 
 deallocate(AkorAllNets_k(netnum)%p)     ! deallocate the old storage
 AkorAllNets_k(netnum)%p => pint         ! make the pointer to point to the storage with more space
 nullify (pint)                                              ! free the pointer
 !
 allocate (pint(length+nn), stat=loc_alloc_stat)  
 if (loc_alloc_stat >0) then 
   print *, "ReadAKorArraysOneNet_DGV: Allocation error for variable pint (II)"
   stop
 end if
 pint(1:length) = AkorAllNets_phi(netnum)%p  ! copy the content of Akor_phi to the larger storage
 deallocate(AkorAllNets_phi(netnum)%p)  ! deallocate the old storage
 AkorAllNets_phi(netnum)%p => pint      ! make the pointer to point to the storage with more space
 nullify (pint)                         ! free the pointer
 !	       

 ! now we open the file, again, to populate the A-arrays
 open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
 read (15, iostat = code_line) mm, nn, lockorobparams
 read (15, iostat=code_line)  scrap_capphi,AkorAllNets(netnum)%p(length+1:length+nn),&
                              AkorAllNets_k(netnum)%p(length+1:length+nn),&
                              AkorAllNets_phi(netnum)%p(length+1:length+nn)
				
 close (15)                                                                     
 !!!
 ! we will need to extend the A_capphi array and we expect the records be kept 
 ! in the increasing order
 AkorAllNets_capphi(netnum)%p = AkorAllNets_capphi(netnum)%p + scrap_capphi
length = length+nn 
end do
!
deallocate(scrap_capphi)
!
end subroutine ReadAKorArraysOneNet_DGVII

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine ReadAKorArraysAllNets_DGV 
!  This subroutine will read Akor arrays for all nets from the hard drive
! 
!  KorNetsChunks -- this is the number of nets with corresponding chunks
!  
!  nuchks -- the total number of chunks to read. 
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadAKorArraysAllNets_DGV  

use  DGV_commvar, only:   KorNetsChunks,numKornets  ! the array containing the size of the aray, parameters P, a1,...,a6


integer (I4B) :: numchks ! the number of chunk in which Akor array is divided for the given Korobov net
integer (I4B) :: n,i     !scrap veriable


n=size(KorNetsChunks,1)
numKornets=n ! Remember the number of nets that will be used. 

if (n > 15) then 
   print *, "ReadAKorArraysAllNets_DGV : Number of Korobov nets is too big. Max of 15 nets are currently supported"
   stop
 end if
! Now we start the loop that goes through all nets
do i=1,n
 numchks=KorNetsChunks(i)
 call ReadAKorArraysOneNet_DGV(i, numchks)      ! calling the subroutine that goes through all chunks at a given net number
end do

end subroutine ReadAKorArraysAllNets_DGV  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeBaseNameAKorNetChnksDGV(
!
! This function creates a recognizable name for the files to store chunks of operator A. It should be used when data is writted in several chunks
! The file name will be formed by adding the chunk number 00X to the base name.
!
! a lot of parameters is cited in the file name. This function will set up the basic 
! file name that will be common for the experiment and will be used by many output procedures.
!
! DO not detouch from the main program. 
! 
! Do not call untill parameters su,sv,sw, current_solution_dir, current_solution_base_name, cfl, &
!          Mu,Mv,Mw,curr_time,u_nonuniform_mesh_type, mesh_u_uniform
!           v_nonuniform_mesh_type, mesh_v_uniform, w_nonuniform_mesh_type, mesh_w_uniform
!          are set in the main program
!
! most of the variables are looked up in the common_variables_block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MakeBaseNameAkorNetChnksDGV(file_name,netnumber,i)

use DGV_commvar, only: su,sv,sw,Mu,Mv,Mw, current_solution_dir, current_AKoroperator_base_name, &
          Mu,Mv,Mw,u_nonuniform_mesh_type,mesh_u_uniform, &
          v_nonuniform_mesh_type,w_nonuniform_mesh_type,mesh_v_uniform,mesh_w_uniform
!
intrinsic Trim, Min, Adjustl
!
integer (I4B) :: i ! the index of the chunk (starts from zero)
integer (I4B) :: netnumber !the index of the net (starts from zero)
character (len=132), intent (out) :: file_name ! the variable to store the file name 
character (len=20) :: parchar    ! string to keep char value of the paameter 
!
! first, we prepare the file name to store the solution
write (parchar, "(I2.2)") netnumber
file_name = trim(current_solution_dir)//trim(current_AKoroperator_base_name)//"net"//trim(Adjustl(parchar))//"_"
write (parchar, "(I3.3)") i
file_name = trim(file_name)//"ch"//trim(Adjustl(parchar))//"_"
! add some parameters to the name so it will be easier to recognize the file:
write (parchar, "(I3)") su
file_name = trim(file_name)//trim(Adjustl(parchar))//"su"
write (parchar, "(I3)") sv
file_name = trim(file_name)//trim(Adjustl(parchar))//"sv"
write (parchar, "(I3)") sw
file_name = trim(file_name)//trim(Adjustl(parchar))//"sw"
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
end subroutine MakeBaseNameAKorNetChnksDGV


end module DGV_readwrite
