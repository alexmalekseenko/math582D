!
!  DGV_mpiroutines.f90
! 
!  Alex 05/10/2015
!
!  Miscellaneous subroutines to use with MPI execution of the DGVlib library. 
!  . 
!!!!!!!!!!!!!!!!!!!!!!1

module DGV_mpiroutines

use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 


implicit none

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! PrepMPICollsnDGV_highperf
 !
 ! This subroutines contains an MPI fork 
 ! This subroutine is run on all processors
 !
 ! Thi subroutine prepares MPI environments and data arrays that will be 
 ! used for MPI parallelization for 
 !    (1) evaluation of the Boltzmann collision operator
 !    (2) can be used for parallel evaluation of the linearized Boltzmann collision operator
 !    (3) can be used with a spatial driver to parallilize in velocity space. Specifically, nodes of the 
 !         velocity mesh are divided between MPI communicating processors. Each processor is running 
 !         a copy of spatial operator for a portion of velocity nodes, the processor request values of the collision 
 !         operator for their portion of nodes. 
 ! 
 ! Description: 
 !
 ! This subroutine prepares a number of communicators. 
 !
 !  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
 !
 !  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
 !                                  sharing a copy of A
 ! 
 !  MPI_Acopy_UNIVS(1:num_Acopies), MPI_Acopy_UNIVS_GROUP(1:num_Acopies) - pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors sharing a copy of A
 !  
 !  LinFmA_UNIV_GROUP, LinFmA_UNIV, -- these are scrap variable used to create a universe for a group of processors 
 !                                  that receives information from a single group of processors sharing a copy of A
 !                                  Here is some idea on how these communicators work: 
 !										There is exactly num_Acopies of groups of processors. Let's sasy its number is 0 <= I <=num_Acopies.
 !										Each of the goups has a workload, given by (nodes_proc_groups(1,I), \ldots, nodes_proc_groups(2,I)), 
 !										that is a number of velocity points 
 !										at which this group of processors I evaluates the collision operator.
 !                                      Just in case, this group of processors is included in the communicator MPI_Acopy_UNIVS(I), but that is 
 !										not important at the moment. Values of the collision operator that are computed in the group I   
 !										may be neded on other processors for either evaluation of 
 !										linearized collision operator or for parallel evaluation of the kinetic equations using 
 !										parallelization in the velocity nodes. This is determined by examining array procs_lim, that 
 !										a range of velocity points to each participating processor. If a processor needs a value of the collision 
 !										operator computed in the group I, then the processor 
 !                                      is included in the LinFmA_Univ/_Group. Also, the master node of the MPI_Acopy_UNIVS(I) is included 
 !										included in the LinFmA_Univ(I) (and also made a master node in it?) 
 !										Typically, the master node of MPI_Acopy_UNIVS(I) broadcasts values to all of LinFmA_Univ(I)
 !
 !  MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP -- pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors that receives information 
 !                                  from a single group of processors sharing a copy of A
 !
 !
 !
 !!! MPI FORK
 !!! Now the master process will distribute the chunks of A arrays to other processors. If the  
 !!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
 !!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
 !!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
 !!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
 !!! Therefore, the results will be combined in an integer of (I8B) kind. 
 !!!
 !!! This number will be divided by the number of 
 !!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
 !!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
 !!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
 !!!  
 !!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
 !!! This boradcast will go out with the code 30.
 !!!
 !!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
 !!! 
 !!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
 !!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
 !!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
 !!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
 !!! 
 
 !!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
 !!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
 !!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
 !!!? workload. 
 !!! 
 !!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
 !!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
 !!!
 !!! Next the master processor initiates the work: It simply broadcasts the solution, 
 !!! 
 !!!
 !!! The slave processes will receive the solutions and will start to evaluate their portions of the collision operator. 
 !!! Then they will send their results to the master processor using an aggregate reduce operation.
 !!!
 !!! Finally, let us briefly describe the process of calculating the collisoin opperator by slave 
 !!! processors. once the I1workload is set,  the slave processor will wait for a boradcast of solution. 
 !!! Each slave process then goes over the I1workload and performs the 
 !!! evaluation of the collision integral for each I1 listed there. The results are stored in the array (CollisionInt)
 !!! Once all I1 are processed, this array is aggregately sent to the master node. Operation of the summation is used to combine 
 !!! results for all chunks. 
 !!! 
 !!! Also, the slave process will evaluate the linearized part
 !!! also, sometimes the slave processor will send the linearized part to some other processor for faster assembly...
 !!!
 !!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrepMPICollsnDGV_highperf(irank)

use DGV_commvar, only: numchnks, num_Acopies,nodes_u,A_capphi,procs_nodes_wld,&
                   lin_proc_nodes_wld,MPI_LINEAR_COLL_WORLD, MPI_LINEAR_COLL_WORLD_GROUP,&
                   num_lin_proc,MPI_Acopy_UNIVS,MPI_Acopy_UNIVS_GROUP,Acopy_wrkld,procs_nodes_wld_Acopy_addr,&
                   MPI_LinFmA_UNIVS_GROUP, MPI_LinFmA_UNIVS, My_Acopy_Univ_indx, &
                   Lin_Self_Send_FmA_Flag, linprocndswld_Acopy_univs, linprocndswld_Acopy_addr, &
                   LinFmA_Univs_callorder,linprocndswld_BfmA_addr,Acopy_irank
                  
use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.


!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 
integer (I8B) :: Act_aggreg ! variable to keep the total number of records in all chunks of operator A
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B), dimension (:), allocatable :: chnks_Act ! This is essentially a scrap array that keeps the information about number of records or A operator in each chunk of Aarrays.
integer (I4B), dimension (:,:), allocatable :: mpiCollisionWorkAlloc ! This array will contain info about work allocations for slave processors  
                                                                ! The array is operating like the following: the first index is the number of processor
                                                                ! the second index has 3 fixed fields 
                                                                ! record 1: number of first chunk to read
                                                                ! record 2: offset in the chunks to read
                                                                ! record 3: number of records to read
integer (I4B) ::  mpiAload,  mpiAloadrem, mpiAload_temp  ! Scrap variables to store the averaghe number of records in A array that will be stored  in the 
                                             ! each slave processor and the reminder of the division of the total number of records in A by the number of 
                                             ! processors                                                                
integer (I4B) :: i,nn,size_est,size_est1,max_lin_send ! scrap index
!!!!! MPI variablessortNodesProcsBuffer
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer     ! this buffer will be used in sending and receiving the work load.
integer (I4B), dimension (:), allocatable :: sortNodesProcsBuffer ! this buffer will be used for sorting nodes

integer :: ierr, mpicommworldsize, iprc, sendproc
 ! variables for MPI Calls
integer, dimension (1) :: ranksII_scrap ! This is a scrap variable to check ranks of the processors in the MPI_LINEAR_COLL_WORLD communicator
!!  Acopy_UNIV, Acopy_UNIV_GROUP ! on the master node this is a temporary holder for the pointer to a universe, and communicatio ngroup that contains processors dealing with one copy of A operator
!!  on the slave node this is avariable for the Acopy communicator
integer, dimension (MPI_STATUS_SIZE) :: istatus
integer :: MPI_COMM_WORLD_GROUP ! to keep the handle to the group of MPI_COMM_WORLD
integer, dimension(3) :: ranges ! ranges to create a communicator group
integer, dimension(6) :: ranges2 ! ranges to create a communicator group for the Acopy universes
integer :: Acopy_UNIV_GROUP,Acopy_UNIV,LinFmA_UNIV_GROUP,LinFmA_UNIV ! scrap pointers to Acopy and LinFmA universes
integer (I4B), dimension(:), allocatable :: LinFmA_recv_univs ! This array will hold the list of MPI_LinFmA_UNIVS communicators from which it will expect to receive components of the linearized operator and the velocity nodes
integer (I4B), dimension(:), allocatable :: ranks_FmA_UNIV !
integer (I4B), dimension(:), allocatable :: rashk_temp !
 
!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:), allocatable :: chnks_copies !
integer (I4B), dimension (:), allocatable :: total_A_capphi ! scrap array to distribute processors into groups and to distribute the workload.
integer (I4B) :: jj,j,ii,ofst,leftover,nds,node,fstchnk,numrec,nnds,proc, num_ranks ! scrap counters

integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:,:), allocatable :: procs_lin ! this one will contain the nodes allocation to the group of processors consolidating the linearized operator. 
                          !it is a scrap variable and once the jobs are disctibuted, it will be deallocated...        

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Here we set up some variables that will be used by all processes: 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! we now allocate the buffer for the mpiCollisionWorkAlloc array
 allocate (mpiI4BArryBuffer(1:3*(mpicommworldsize-1)), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf Allocation error for variables mpiI4BArryBuffer"
  stop
  end if
!!!!!!!!!!!!!!
!!  Next we create a communicator that will include the master node (irank =0) and all nodes with rank <=num_lin_proc-1
!!!!!!!!!!!!!!
! first we need to create a communicator group. All processes execute it:
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_WORLD_GROUP, ierr) ! get the hanlde to the group of MPI_COMM_WORLD. 
ranges(1)=0; ranges(2)=num_lin_proc-1; ranges(3)=1; 
call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, MPI_LINEAR_COLL_WORLD_GROUP, ierr) ! MPI_LINEAR_COLL_WORLD_GROUP will be dedicated to communication of 
                           ! processors with numbers 0, ..., num_lin_proc-1
                           ! num_lin_proc is set in DGVparameters.dat 
                           ! In many instances num_lin_proc is set to the number of total MPI processors requested for the run
                           ! But it can be a smalle number too. If num_lin_proc is equal to the total number of MPI processors, then 
                           ! MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_WORLD_GROUP will refer to the same processors                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_highperf: MPI can not create group MPI_LINEAR_COLL_WORLD_GROUP on proc", irank, ". Error=", ierr
         stop
end if
! a quick sanity check if the processor with irank=0 is still with the rank =0 in the new group.
call MPI_GROUP_TRANSLATE_RANKS(MPI_LINEAR_COLL_WORLD_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_highperf: can't translate ranks MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc",irank,".Err=",ierr
 stop
end if
! a quick check that the zeros process is still number zero in the new group: 
if (ranksII_scrap(1) /= 0) then 
 print *,"PrepMPICollsnDGV_highperf: ranks mixed up between  MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc", irank, "stop."
 stop
end if 
! end sanity check 
! next we create a new communicator for the MPI_LINEAR_COLL_WORLD_GROUP: 
call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, MPI_LINEAR_COLL_WORLD, ierr) ! MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for 
                           ! processors with numbers 0, ..., num_lin_proc-1
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_highperf: can't create communicator MPI_LINEAR_COLL_WORLD. Proc=", irank, ". Err=", ierr
 stop
end if
! End create communicator for consolidated handling of linearised collision operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we prepare three important arrays: num_chnks, nodes_proc_groups and procs_lin. 
! These arrays determine the workload of nonlinear evaluation of 
! the collision operator and the workload for evaluation of the linearized collision operator, 
! respectively. 
! First, let us work on the array chnks_copies(2,1:num_Acopies)
allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables chnks_copies. stop"
  stop
end if
chnks_copies=0 ! nullify, just in case.
!! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
!! in particular, a copy of A will be distributed between processors with number chnks_copies(1,1)=1 (master #0 does not get any) 
!! and chnks_copies(2,1)=XX. Another copy will be distributed between chnks_copies(1,2)=XX+1 and chnks_copies(2,2)=YY and so on.
!! 
call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
!! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
!! (the number fo the first and the last processors that contain this particular copy of A. num_Acopies  -- is the number of copies) 
!! 
!! Next we distribute workload between processor groups. It is assumed that since a copy of A is shared 
!! between processors in a group, all of them are needed to evaluate the collision operator for a single 
!! velocity point. This is true for s=1, but for s>1 not all processors in the group are needed to evaluate 
!! collision operator at a velocity node. This can be used later to tune up the work load a bit better later. 
!! At this moment, however, we will take all nodes where the colliision operator needs to be evaluated and we 
!! divide these velocity nodes between the groups of processors, thus creating a workload for each group.
!! Later this array will be used to assign individual workloads for 
!! each processor. 
 nn = size(nodes_u,1) ! Total number of nodes numbered beginning from 1 where the collision integral needs to be evaluated.  
 !! first we divide all velocity nodes between the groups of processors. We have exacly (num_Acopies) processor groups 
 allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables nodes_proc_groups. stop"
  stop
 end if
 call DivideRangeContinuous(nodes_proc_groups,1,nn)
!! we divided velocity nodes between the groups of processors sharing copies of A.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! We make one more preparation for the future: we 
!! setup the array procs_lin that tells the workload for each processor involved in the evaluation of the 
!! linearized collision operator. Similarly to the above, simply divide velocity nodes between the 
!! processors involved in the evaluation of linearied collision operator.  
 allocate (procs_lin(2,1:num_lin_proc), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables total_receive_procs,total_send_procs,procs_lin"
  stop
 end if
!! we populate the linearized group of processors with assgined nodes
 call DivideRangeContinuous(procs_lin,1,nn)
!! end of creating proc_lin: all processors are divided in groups and each group is assgned some number of velocity nodes 
!! Note that if 1D-3D spatial operator is parallelized in velocity variable, proc_lin can be used to determine the work load. 
!! for each processor 
!! arrays procs_lin, nodes_proc_groups and chnks_copies have been created !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! We now create universes to deal with the groups of processors sharing a single copy 
!!! of A and the universes of processors that obtain 
!!! components of the colliion operator from the same group of nodes... 
!!!
!!! These universies will be used for communcation of the components of the linearized operator.
!!! These universes can also be used for evaluation of the collision operator if the spatial operator is parallelized  
!!! 
!!! Begin creating the Acopy_UNIV universes:
!!! first, we allocate storage for the universes... 
 allocate (MPI_Acopy_UNIVS(1:num_Acopies),MPI_Acopy_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
    ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_Acopy_UNIVS=0;MPI_Acopy_UNIVS_GROUP=0; ! Nullify them , just in case...
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Reset the index variable 
 My_Acopy_Univ_indx = 0 ! points to no particular Acopy Universe
 ! now we will create num_Acopies universes that will contain the master node and the nodes indentified in chnks_copies(1:2,j)
 do i=1,num_Acopies
  ranges(1)=chnks_copies(1,i); ranges(2)=chnks_copies(2,i); ranges(3)=1; 
  call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, Acopy_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: MPI can not create group Acopy_UNIVS_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! a quick check if the processor with irank=0 int he new group correspond to the first processor in the group.
  call MPI_GROUP_TRANSLATE_RANKS(Acopy_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: can't translate Acopy_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the zeros process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_highperf: ranks mixed up between Acopy_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, Acopy_UNIV_GROUP, Acopy_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: can't create communicator Acopy_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_Acopy_UNIVS(i)=Acopy_UNIV
  MPI_Acopy_UNIVS_GROUP(i)=Acopy_UNIV_GROUP
  !!!!
  ! We also need to remember to what Acopy universe we belong
  if ((irank >= chnks_copies(1,i)) .and. (irank<= chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! now we know which Acopy universe contains this slave processor. Master processor will still have zero in this variables
  end if 
 end do  
 !! END creating the Acopy_UNIV universes:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will use the costant My_Acopy_Univ_indx to set 
 ! the variable Acopy_irank which defines the rank of the processor in this A copy Universe
 !!!!!!!!!!!!!!!! 
 Acopy_irank = -1
 if (My_Acopy_Univ_indx>0) then 
 !!! MPI FORK
   Acopy_UNIV=MPI_Acopy_UNIVS(My_Acopy_Univ_indx) 
   call mpi_comm_rank(Acopy_UNIV,Acopy_irank,ierr) ! check what processor the program is on... 
   if (ierr /= 0 ) then 
    print *,"PrepMPICollsnDGV_highperf: can't determined rank in Acopy_UNIV.i=",My_Acopy_Univ_indx,"Proc=",irank,". Err=",ierr
   stop
  end if
 end if   
 !!!!!!!!!!!!!!!
 ! end setting Acopy_irank
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !! Begin creating LinFmA_UNIV universes:
 allocate (MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables MPI_LinFmA_UNIVS, MPI_LinFmA_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_LinFmA_UNIVS=0; MPI_LinFmA_UNIVS_GROUP=0; ! Nullify them, just in case ...
!! 
!! Next we need will pick one group of processors sharing a copy of A and find all processors that will 
!! receive components of the linearized operator from this group of processors.  
!! the first processors in the group and all processors receiveing compoents will be added to the 
!! communicator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Lin_Self_Send_FmA_Flag = 0 ! need to reset this flag before start...
 ! need the following temporary arrays
 size_est = 1 + (INT((nn/num_lin_proc),I4B)+2)*2 !this should estimate from above the size of LinFmA_recv_univs
 !  
 allocate (rashk_temp(1:num_lin_proc+1),LinFmA_recv_univs(1:size_est), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables rashk_temp, LinFmA_recv/send_univs. stop. irank", irank 
  stop
 end if
 ! clean up before proceeding
 LinFmA_recv_univs = 0   ! 0 correspond to no communicator containing a single copy of A
 ! begin creating universes MPI_LinFmA_UNIVS
 do i=1,num_Acopies
  ! first we need to set up the range of ranks that will form the communicator group. 
  rashk_temp = -1  ! clean the array
  num_ranks = 1
  rashk_temp(num_ranks) = chnks_copies(1,i) ! this processor will serve as a master processor for this comminucator group.
  ! now we will start screenign the processors participating in the evaluation of the linearized operator. Recall that the 
  ! chnks_copies array contains irank values of the processors that share a copy of A. Also nodes_proc_groups contains 
  ! the first and last of velosity nodes that are assigned to each processor group (aka Acopy group).
  ! We will go over each processor evaluating linearized collision operatot and over each velocity node 
  ! that are assigned to each of them. For each velocity node, we check if that node is assigend to the 
  ! current Acopy group, which is the group containing the processor of rank rashk_temp(1). 
  ! If they do, then this processor will be needing information from that Acopy group. We, therefore, 
  ! add the irank fo the processor to the array of ranks that will comminicate to this Acopy group, the 
  ! rashk_temp() array. 
  do j=1,num_lin_proc ! j=1 points to irank=0
   do jj=procs_lin(1,j),procs_lin(2,j) ! runs over all nodes of linearizaed operator assiged to this processor 
    if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then 
      ! first we check if this this processor (must be a slave proc) is, first, a master proc of an Acopy universe 
      ! and, second, needs to get components of FmA from that universe. In this case it have dual role in FmA -- as the master and the slave -- 
      ! to avoid confusion, we do not add it to the rashk_temp since it is already there. 
      ! also, we set up a spacial flag for that 
      if (j-1 == chnks_copies(1,i)) then 
       if (irank == chnks_copies(1,i)) then 
        Lin_Self_Send_FmA_Flag = 1 !! this means that this processor (must be a slave proc) is a master proc of an Acopy universe and needs to get components of FmA from this universe 
       end if 
      else 
       num_ranks=num_ranks+1
       rashk_temp(num_ranks)=j-1  ! this will correspond to irank
      end if
      exit ! the recieving processor has been added to the ranks. No need to find more nodes that this processor recieves
    end if   
   end do 
   ! next we need to fill some auxiliary array: all nodes that will be recieved will be recorded.
   if (irank == j-1) then
    do jj=procs_lin(1,j),procs_lin(2,j) ! runs over all nodes of linearizaed operator assiged to this processor 
     if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then
       ! Now we fill the temporary receive array # pairs of records # A copy universe # node   # A copy universe # node ... 
       LinFmA_recv_univs(1) = LinFmA_recv_univs(1)+1     ! this will add the communicator (FmA_UNIV) to the list of communicators
       LinFmA_recv_univs(2*LinFmA_recv_univs(1)) = i   ! this records the number of A copy universe from which components of the linearized  operator are coming from
       LinFmA_recv_univs(2*LinFmA_recv_univs(1)+1) = jj ! this records the node to which the componets correspond.     
     end if 
    end do
   end if 
  end do 
  ! we have temporary array of ranks, we now create the final array of ranks:
  allocate (ranks_FmA_UNIV(1:num_ranks),  stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variables ranks_FmA_UNIV. stop. irank", irank 
   stop
  end if
  ranks_FmA_UNIV(:) = rashk_temp(1:num_ranks)
  ! next we create the group for the FmA_UNIV
  call MPI_GROUP_INCL(MPI_COMM_WORLD_GROUP,num_ranks,ranks_FmA_UNIV,LinFmA_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: MPI can not create group LinFmA_UNIV_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! clean up some variables for the next iteration
  deallocate(ranks_FmA_UNIV) 
  num_ranks=0
  rashk_temp=-1
  ! a quick check if the processor with irank=chnks_copies(1,i) has rank =0 in the new group.
  call MPI_GROUP_TRANSLATE_RANKS(LinFmA_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: can't translate LinFmA_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_highperf: ranks mixed up between LinFmA_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, LinFmA_UNIV_GROUP, LinFmA_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: can't create communicator LinFmA_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_LinFmA_UNIVS(i) = LinFmA_UNIV
  MPI_LinFmA_UNIVS_GROUP(i)=LinFmA_UNIV_GROUP
 end do 
 deallocate(rashk_temp)
 !! END setting up MPI communicators MPI_LinFmA_UNIVS and MPI_Acopy_UNIVS
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEBUG 
! print *, "irank:", irank,"LinFMA_univ:", MPI_LinFmA_UNIVS!
!
!!!!!!!!!!!!!!!
!   call mpi_barrier (MPI_COMM_WORLD,ierr)
!!   if (ierr /= 0) then 
 !   print *,"PrepConsolLinColl_estafeta_DGV: MPI_Barrier from proc",irank, "returned error", ierr, "on process", irank
!    stop
!   end if 
!   stop
!!! END DEBUG   
!!!!!!!!!!!!!!!!   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to distribute the workload to the slave processors for 
! evaluating the full collision operator 
!!!!!!!!!!!!!!!!!
 !!  MPI fork
 if (irank==0) then
 !! STEP 1:
 !! Master process code: 
 !! The Master node needs to disctribute operator A. Therefore it will scall all chunks of all arrays and will 
 !! retrieve the total number of the recods in operator A. Then this number will be divided to the number of processors 
 !! to produce the operator A workloads (+/- 1)
 !! begin STEP 1:
  allocate (chnks_Act(1:numchnks), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variables chnks_Act"
   stop
  end if
  !! The next subroutine will scan all chunks of Aarrays and will compile the information about the number of records in 
  !! each chunk in array chnks_Act
  ! call next subroutine to populate chnks_Act
  call ScanAarrsChnks4Acapphi(chnks_Act) ! Act_aggreg will have the total number of records in the Aarrays.
  ! Let us compute the total number of records
  Act_aggreg=0
  do i=1,numchnks
   Act_aggreg=Act_aggreg+chnks_Act(i)
  end do
  !! now we broadcast the array to all processors. 
  !! Now it is time to set up the work allocation array.
  allocate (mpiCollisionWorkAlloc(1:3,1:mpicommworldsize-1), stat=loc_alloc_stat)
     ! 
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperfV: Allocation error for variables mpiCollisionWorkAlloc"
   stop
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mpiCollisionWorkAlloc = 0 ! reset the array
  do jj=1,num_Acopies
   !! Next we divide the total number of records by the total number of available processors to determine the workload per processor
   !! We perform division with reminder
   mpiAload = INT(Act_aggreg/(chnks_copies(2,jj)-chnks_copies(1,jj)+1),I4B)
   mpiAloadrem = INT((Act_aggreg - (INT8(chnks_copies(2,jj)-chnks_copies(1,jj)+1)*mpiAload)), I4B) ! notice that mpiAloarrem < mpicommworldsize-1
   !! A quick check that we have enough memory per processor. If we do not, stop:
   if (mpiAload > 2*10**8) then 
    print *,"PrepMPICollsnDGV_highperf: Error> the number of records of operator A exceeds 2*10^8. Aborting the solver."
    stop
   end if 
   !! Next we populate the work allocation array.
   !! There will be two cases. Case 1 -- number of slave processors= number of chunks 
   !! and Case 2: number of processors neq number of chunks ... 
   if (numchnks == (chnks_copies(2,jj)-chnks_copies(1,jj)+1)) then 
   !! in the case 1 the chunks are assigned to the processor withour rechunking. 
   !! SPECIAL CASE :: NUMBER OF CHUNKS AND NUMBER OF PROCESSORS IS THE SAME! !! IMPORTANT!!
    do i = chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = i - chnks_copies(1,jj) + 1
     mpiCollisionWorkAlloc(2,i) = 0
     mpiCollisionWorkAlloc(3,i) = chnks_Act(i - chnks_copies(1,jj) + 1)
    end do
   else 
   !! in the case 2, the following subroutine determines re-chunking the A-arrays. 
    j=1 ! this is the index that will run over chunks
    ofst=0 ! This is the offset calculator
    leftover=0
    do i=chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = j  
     mpiCollisionWorkAlloc(2,i) = ofst
     if (i <= mpiAloadrem) then
      mpiAload_temp = mpiAload+1
     else 
      mpiAload_temp = mpiAload
     end if   
     do while (mpiCollisionWorkAlloc(3,i) < mpiAload_temp) ! process chunks until i-th processor is loaded. 
      ! a quick check if data is corrupted. if our calculations above arte correct, then there will be exactly as many 
      ! chunks as needed. If we ran oput of chunks -- somewhere there is an error. Then we stop
      if (j > numchnks) then 
       print *, "PrepMPICollsnDGV_highperf: error when trying to populate mpiCollisionWorkAlloc. Ran out of chunks."
       stop
      end if
      ! end error check  
      if ((chnks_Act(j) - ofst) >= (mpiAload_temp - mpiCollisionWorkAlloc(3,i))) then ! if the current chunk has more records than average load, 
       ofst = ofst + mpiAload_temp - mpiCollisionWorkAlloc(3,i)                    ! set up the offset for the next processor (part of the chunk j will go there) 
       mpiCollisionWorkAlloc(3,i) = mpiAload_temp ! assign the load
       ! Need to include a check in case the chunk has exactly as many records as is needed 
       ! in this case we need to move onto the next chunk for other procs...
       if (ofst == chnks_Act(j)) then 
        j=j+1 ! move the pointer to the next chunk
        ofst = 0
       end if  
      else 
       mpiCollisionWorkAlloc(3,i) = mpiCollisionWorkAlloc(3,i) + chnks_Act(j) - ofst ! count in the remaining records from chunk j. However, more records needed for this processor 
       j = j+1                                    ! move on to the next chunk
       ofst = 0                                   ! rest the offset
      end if
     end do   
    end do
   end if 
  end do ! end loop in copies of A 
  ! the array of work allocation has been produced. Now we need to pass it on to the slave processors via a boradcast 
  ! prepare the buffer
  do i=1, mpicommworldsize-1
   do j=1,3
    mpiI4BArryBuffer((i-1)*3+j) = mpiCollisionWorkAlloc(j,i)
   end do 
  end do 
  ! send the buffer via boradcast
  call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_highperf: MPI broadcast from proc 0 returned error", ierr
   stop
  end if  
  !! The workload has been passed onto the slave processors.
  deallocate(mpiI4BArryBuffer,chnks_Act) ! The master processor will not use these anymore... 
 !! the evaluation of collision integral will be performed in an ierarchical fashion: Each group of processors will be 
 !! united in a universe Acopy_UNIV. The first processor in the groupd will recieve rank = 0 and will be assignedas the master for
 !! this univers. The slave processors of a copy of Acopy_UNIV will send data to the master using collective communication.
 !! Also the processor with rank zero can still collect data frorm all toher processor and use broadcasts send the data out to 
 !! slave processors in all Acopy_UNIV universes 
 !! END OF THE STEP 1:
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:  
 !! Now each of the slave processors will read the A-arrays and prepare their local copies of A-arrays. 
 !! The slave processors will produce thier local A-capphi arrays. This information will be important on the next step: 
 !! the determination of each individual processor's job allocation in evaluation of the collision integral. 
  
 !! It is decided to let the slave processors divide the velocity nodes between the groups of processors
 !! essentially each will run acopy of the subroutine and should have the same answer. Then they will pick their portion of 
 !! nodes and will begin process them. In the end they will produce the workloads for themselves. Workload is the list of 
 !! velocity nodes for which the slave processor evaluates the collision operator 
 !!
 !! END OF  STEP 2
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Finally we setup some auxiliary arrays that we will need in the job assignment routin. 
 !! (We really just need the nodes_phican but we want to use the subroutine that sets all the arrays...
 !! first, we need to create the array A_capphi. it is not used on the master node, but we need a dummy array to 
 !! pass it to the subroutine so as to avoid mistakes.
 nn = size(nodes_u,1)
 allocate (A_capphi(1:nn), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable A_capphi. stop"
   stop
  end if
 A_capphi=0 ! This array is not used on the master node
 !
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! End of preparing the array nodes_phican
 !!
 !!!!!!!!!!!!!!!!!!!!! 
else
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!
 !! SLAVE PROCESS CODE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 1: 
 !! 
 !! The slave processor receives the A-array allocation that will be broadcasted from the master node: 
 call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"PrepMPICollsnDGV_highperf: slave processor", irank, "MPI boradcast from proc 0 returned error", ierr
  stop
 end if   
 !! Now we need to retreve three numbers that pertain to this particular processor. 
 !! Specifically, we need to get the fstchnk=mpiCollisionWorkAlloc(irank,1), 
 !! ofst = mpiCollisionWorkAlloc(irank,2) and numrec = mpiCollisionWorkAlloc(irank,3) 
 !! The rest of the array is not important for this process. Therefore, 
 fstchnk=mpiI4BArryBuffer((irank-1)*3+1)
 ofst=mpiI4BArryBuffer((irank-1)*3+2)
 numrec=mpiI4BArryBuffer((irank-1)*3+3)
 !! END OF STEP 1:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:
 !! We proceed processing the received information 
 !! Now each processor needs to read the corresponding part of the A array and get ready for work. 
 call ReadAarraysChnksOfstNrecDGV(fstchnk,ofst,numrec,numchnks)
 !! Now local copies of the A-arrays are set. We are ready to detemine the work load in terms of the velocity nodes 
 !! Need to run this setup subroutine: Presently we really just need the nodes_canphi but we want to use the subroutine that sets all the arrays...
 !! sets up arrays: nodes_Ashift,nodes_phican, nodes_dui,nodes_dvi,nodes_dwi
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! call the subrouine that evaluates the workload: velocity nodes that will be processed on this processor.
 call SetMPILocProcNodesWorkLoad_DGV(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)
 !! the node work allocation is set 
 deallocate(mpiI4BArryBuffer)
 !!
 !! END OF STEP 2.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! looks like we all done on the slave node...
end if 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! End distributing workload for the slave processors for evaluation of full Boltzmann operator.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! STEP 3.
 !!
 !! PARALLELIZATION OF THE LINEARIZATION STAFF 
 !!
 !!
 !! This part of the code replaces the original verison that is using 
 !! send_lin_procs and recv_lin_procs arrays for the construciton of the linearized collision operator
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! ATTENTION: OUTDATED!!!!
 !! The new algorithm for evaluatin the collision operator is as follows:
 !!
 !! A. Slave Processors  compute the components of the linearized operator.
 !! B. Slave processors place these componets in the copy of an array 
 !!    FmA whose rows go over all nodes in the array Acopy_Wrkld that are the same for 
 !!    all processors in the same Acopy universe
 !!    to help them stack the components the array procs_nodes_wkld_Acopy_Addr array us used. 
 !!    procs_nodes_wld_Acopy_addr has the same amount of elements as procs_nodes_wld
 !!    The array procs_nodes_wld_Acopy_addr contains indices of the nodes from 
 !!    procs_nodes_wld in Acopy_Wrkld
 !!    The lead processor of Acopy univ. gathers the FmA to it
 !! 
 !! C. The lead proc of Acopy universe is also the lead processor of the matching LinFmA universe 
 !!    This processor broadcasts FmA to processors involved in the valuation of the collision oeprator. 
 !!    For each Acopy universe there is a LinFmA_Univ universe that combines all processors that need 
 !!    to receive data from processors invloved in Acopy_univ. 
 !! D. The processors pick the data from the boradsact and build the local copy of the linearized operator.
 !!    To do that three arrays are used:
 !!    lin_proc_nodes_wld --- contains the componets of the linearized operator assigned to each processor involved 
 !!    in the evaluation of the linearized operator. 
 !!    linprocndswld_Acopy_univs, --- this array contains informaion about arriving components of the linearized operator:
 !!     # of the universe A copy that will send the components
 !!     # total number of the nodes (each node is one row of the lineaired collision operator
 !!     # Nodes .... 
 !!    then repeat
 !!    linprocndswld_Acopy_addr has the same stucture except instead of nodes, there will be their indices in the corresponding 
 !!    Acopy_Wrkld array
 !!
 !! E. The processors involced in the evaluation of the linearized colliion operator create a copy of BFmA where they store the 
 !!    componets of the consolidated operator.  
 !!    the processors evaluate the parts of the consolidated collision operator and stack the results in the temporary solution 
 !!    array. Finally master processrs gathers all componets of the temp. solution array by a collective communication call. 
 !!!! END OUTDATED  !!
 !!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set up the arrays  
 !! Acopy_Wrkld and procs_nodes_wld_Acopy_addr
 !! These arrays only make sense on slave nodes...
 !!
 !!!!!!!!!!!!!!!!!!!!!!
 if (irank>0) then 
 !! essentially for each Acopy universe, we need to write all the velocity nodes assigned to the Acopy group to 
 !! to Acopy_Wrkld. (Comment: in this context, "assigned velocity node" means the node at which collision integral is evaluated.)
 !! Importance of Acopy_Wrkld is the following. When results for evaluating the collision 
 !! operator are collected across processors sharing one copy of A, a.k.a. an Acopy group, these results are written in an array
 !! with results for different nodes arranged as in Acopy_Wrkld. 
 !! Individual processors in the Acopy group may have different nodes workloads. As a result, the locally computed results, 
 !! will be arranged in an array differently shaped than the results for the Group in total. When data is collected across 
 !! the Acopy group, the results need to be copied to correct places in the results array for the Group. For this, we 
 !! need a translation table array that tells for each entry of a localresult, that is arranged as procs_nodes_wld, where 
 !! the entry should go in the result array for the group which uses the same indexing of work nodes as in Acopy_Wrkld.  
  nn = nodes_proc_groups(2,My_Acopy_Univ_indx) - nodes_proc_groups(1,My_Acopy_Univ_indx) + 1 
  allocate (Acopy_Wrkld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (Acopy_Wrkld) on node", irank, "Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the Acopy_Wrkld
  do i = 1,nn 
   Acopy_Wrkld(i) = nodes_proc_groups(1,My_Acopy_Univ_indx) + i -1 !  Acopy_Wrkld(i) contains all velocity nodes assigned ot this Acopy group beginning 
                   ! from  nodes_proc_groups(1,My_Acopy_Univ_indx) and ending at  nodes_proc_groups(2,My_Acopy_Univ_indx)
  end do 
  !! end populating the Acopy_Wrkld
  !! now let us prepare the array procs_nodes_wld_Acopy_addr 
  nn = size(procs_nodes_wld,1)
  allocate (procs_nodes_wld_Acopy_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (procs_nodes_wld_Acopy_addr) on node",irank,"Stop"
   stop
  end if
  !! we go over the procs_nodes_wld array and find the matcning components in Acopy_Wrkld and record indices
  !! of these components in procs_nodes_wld_Acopy_addr
  procs_nodes_wld_Acopy_addr=0 ! reset, just in case, 
  do j=1,size(procs_nodes_wld,1)
   do ii=1,size(Acopy_Wrkld,1)
    if (Acopy_Wrkld(ii) == procs_nodes_wld(j)) then  
    procs_nodes_wld_Acopy_addr(j) = ii
    exit
    end if 
   end do 
  end do
  !! the procs_nodes_wld_Acopy_addr array is set
 end if   
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set arrays lin_proc_nodes_wld, linprocndswld_Acopy_univs,
 !! and linprocndswld_Acopy_addr
 !!
 !! We recall that these arrays only make sense on processors that are involved in the evaluation of the 
 !! linearized operator
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (irank < num_lin_proc) then 
  nn = procs_lin(2,irank+1) - procs_lin(1,irank+1)+1 ! this is the total number of velocity nodes that will be assigned to this processor 
  allocate (lin_proc_nodes_wld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (lin_proc_nodes_wld) on  node", irank, ". Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the lin_proc_nodes_wld
  do i = procs_lin(1,irank+1),procs_lin(2,irank+1)
   lin_proc_nodes_wld(i-procs_lin(1,irank+1)+1) = i ! lin_proc_nodes_wld contains all velocity nodes (actually, their numbers) beginning from procs_lin(1,irank+1) and ending at procs_lin(2,irank+1)
  end do 
  !! end populating the lin_proc_nodes_wld 
  !! Next we will populate the array linprocndswld_Acopy_univs. The information that will be stored in this array is already in the 
  !! array LinFmA_recv_univs. The data just needs to be re-sorted:
  !! In the buffer the first number is the number of records after that go pairs: (#processor, #node) (coment "node" means velocity node)
  !! What we do next is to rearrange this information in the following format:
  !! (#processor,#number of nodes for this processor,#node1, #node2, and so on, #nodeN) 
  !! then continue onto the next processor
  !! This is done in two steps: First, we prepare the buffer for sorting
  !! we now allocate the buffer for the linprocndswld_Acopy_univs array
  nn = LinFmA_recv_univs(1)
  allocate (sortNodesProcsBuffer(1:3*nn+1), stat=loc_alloc_stat) ! 
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variables sortNodesProcsBuffer on proc",irank," stop"
   stop
  end if
  !! now we are calling the sorting routine. The return variable sortNodesProcsBuffer
  !! will contain sorted records
  call SortNodesProcs(LinFmA_recv_univs,sortNodesProcsBuffer) ! the first element in buffer contains the number of pair. The first element on return array contains the number of non-trivial records.
  !! now we create the array and copy the records into it
  nn = sortNodesProcsBuffer(1) ! the first entry is the number of records in the array except for the itself
  ! now we need to create the local send_lin_procs array. 
  allocate (linprocndswld_Acopy_univs(1:nn), linprocndswld_Acopy_addr(1:nn),linprocndswld_BfmA_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error (linprocndswld_Acopy_univs/_Acopy_addr/_BfmA_addr). Process",irank,"Stop"
   stop
  end if
  ! next we populate the array send_lin_procs ... 
  linprocndswld_Acopy_univs = sortNodesProcsBuffer(2:nn+1)
  ! next we need to populate array linprocndswld_Acopy_addr 
  linprocndswld_Acopy_addr=0 ! reset, just in case
  i=1
  do while (i<size(linprocndswld_Acopy_univs,1))
   linprocndswld_Acopy_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
   linprocndswld_Acopy_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
   do ii=1,linprocndswld_Acopy_univs(i+1)
    j = linprocndswld_Acopy_univs(i)         !!! the index of the Acopy universe
    node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
    linprocndswld_Acopy_addr(i+1+ii) = node - nodes_proc_groups(1,j) + 1  !!!! Attention -- we do not really look into  Acopy_Wrkld, but we look into the array that is used to create Acopy_Wrkld
   end do
   i=i+2+linprocndswld_Acopy_addr(i+1)  
  end do 
 !! the linprocndswld_Acopy_addr arrays are set! 
  ! next we need to populate array linprocndswld_BfmA_addr that contains the addresses on the nodes listed in the 
  ! linprocndswld_Acopy_univs as they appear in the consolidated linearied operator. 
  linprocndswld_BfmA_addr = 0 ! reset, just in case
  i=1
  do while (i<size(linprocndswld_BfmA_addr,1))
   linprocndswld_BfmA_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
   linprocndswld_BfmA_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
   do ii=1,linprocndswld_Acopy_univs(i+1)
    node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
    linprocndswld_BfmA_addr(i+1+ii) = node - procs_lin(1,irank+1) + 1  !!!! Attention -- we do not really look into  lin_proc_nodes_wld, but we look into the array that is used to create lin_proc_nodes_wld
   end do
   i=i+2+linprocndswld_BfmA_addr(i+1)  
  end do 
 !! the linprocndswld_BfmA_addr arrays are set! 
 end if 
 !! clean up 
  deallocate(chnks_copies,procs_lin)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! the last piece is to create the call order array LinFmA_Univs_callorder.
 !! this array contains the sequence in which the LinFmA universes are called
 !! It contains numbers from 1 to num_Acopies that are permuted in some way
 !! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 allocate (LinFmA_Univs_callorder(1:num_Acopies),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error (LinFmA_Univs_callorder).Process",irank,"Stop"
   stop
  end if
 LinFmA_Univs_callorder = 0
 !! Now we will populate the array LinFmA_Univs_callorder
 i=1
 j=1
 do while (i<=num_Acopies) 
    LinFmA_Univs_callorder(i)=j
    i = i + 1
    j = j + 2 
    if (j > num_Acopies) then 
     j = 2
    end if 
 end do 
!! looks like we all done 
end subroutine PrepMPICollsnDGV_highperf

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! PrepMPICollsnDGV_Univs_highperf
 !
 ! The next subroutines are chunks of PrepMPICollsnDGV_highperf. The splitting alows for 
 ! more flexibility in customizing MPI codes developed using the DGv library
 !
 ! This piece creates communication universes for primary velocity mesh. 
 ! 
 !
 ! This subroutines contains an MPI fork 
 ! This subroutine is run on all processors
 !
 ! Thi subroutine prepares MPI environments and data arrays that will be 
 ! used for MPI parallelization for 
 !    (1) evaluation of the Boltzmann collision operator
 !    (2) can be used for parallel evaluation of the linearized Boltzmann collision operator
 !    (3) can be used with a spatial driver to parallilize in velocity space. Specifically, nodes of the 
 !         velocity mesh are divided between MPI communicating processors. Each processor is running 
 !         a copy of spatial operator for a portion of velocity nodes, the processor request values of the collision 
 !         operator for their portion of nodes. 
 ! 
 ! Description: 
 !
 ! This subroutine prepares a number of communicators. 
 !
 !  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
 !
 !  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
 !                                  sharing a copy of A
 ! 
 !  MPI_Acopy_UNIVS(1:num_Acopies), MPI_Acopy_UNIVS_GROUP(1:num_Acopies) - pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors sharing a copy of A
 !  
 !  LinFmA_UNIV_GROUP, LinFmA_UNIV, -- these are scrap variable used to create a universe for a group of processors 
 !                                  that receives information from a single group of processors sharing a copy of A
 !                                  Here is some idea on how these communicators work: 
 !										There is exactly num_Acopies of groups of processors. Let's sasy its number is 0 <= I <=num_Acopies.
 !										Each of the goups has a workload, given by (nodes_proc_groups(1,I), \ldots, nodes_proc_groups(2,I)), 
 !										that is a number of velocity points 
 !										at which this group of processors I evaluates the collision operator.
 !                                      Just in case, this group of processors is included in the communicator MPI_Acopy_UNIVS(I), but that is 
 !										not important at the moment. Values of the collision operator that are computed in the group I   
 !										may be neded on other processors for either evaluation of 
 !										linearized collision operator or for parallel evaluation of the kinetic equations using 
 !										parallelization in the velocity nodes. This is determined by examining array procs_lim, that 
 !										a range of velocity points to each participating processor. If a processor needs a value of the collision 
 !										operator computed in the group I, then the processor 
 !                                      is included in the LinFmA_Univ/_Group. Also, the master node of the MPI_Acopy_UNIVS(I) is included 
 !										included in the LinFmA_Univ(I) (and also made a master node in it?) 
 !										Typically, the master node of MPI_Acopy_UNIVS(I) broadcasts values to all of LinFmA_Univ(I)
 !
 !  MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP -- pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors that receives information 
 !                                  from a single group of processors sharing a copy of A
 !
 !
 !
 !!! MPI FORK
 !!! Now the master process will distribute the chunks of A arrays to other processors. If the  
 !!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
 !!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
 !!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
 !!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
 !!! Therefore, the results will be combined in an integer of (I8B) kind. 
 !!!
 !!! This number will be divided by the number of 
 !!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
 !!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
 !!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
 !!!  
 !!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
 !!! This boradcast will go out with the code 30.
 !!!
 !!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
 !!! 
 !!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
 !!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
 !!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
 !!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
 !!! 
 
 !!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
 !!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
 !!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
 !!!? workload. 
 !!! 
 !!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
 !!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
 !!!
 !!! Next the master processor initiates the work: It simply broadcasts the solution, 
 !!! 
 !!!
 !!! The slave processes will receive the solutions and will start to evaluate their portions of the collision operator. 
 !!! Then they will send their results to the master processor using an aggregate reduce operation.
 !!!
 !!! Finally, let us briefly describe the process of calculating the collisoin opperator by slave 
 !!! processors. once the I1workload is set,  the slave processor will wait for a boradcast of solution. 
 !!! Each slave process then goes over the I1workload and performs the 
 !!! evaluation of the collision integral for each I1 listed there. The results are stored in the array (CollisionInt)
 !!! Once all I1 are processed, this array is aggregately sent to the master node. Operation of the summation is used to combine 
 !!! results for all chunks. 
 !!! 
 !!! Also, the slave process will evaluate the linearized part
 !!! also, sometimes the slave processor will send the linearized part to some other processor for faster assembly...
 !!!
 !!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  PrepMPICollsnDGV_Univs_highperf(irank,chnks_copies,nodes_proc_groups,procs_lin)

use DGV_commvar, only: num_lin_proc, nodes_u, num_Acopies, &
                       MPI_LINEAR_COLL_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, &
                       MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP, My_Acopy_Univ_indx, Acopy_irank, &
                       MPI_LinFmA_UNIVS_GROUP, MPI_LinFmA_UNIVS, Lin_Self_Send_FmA_Flag, &
                       LinFmA_Univs_callorder,&
                       linprocndswld_Acopy_univs,linprocndswld_Acopy_addr,linprocndswld_BfmA_addr
use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:) :: chnks_copies !
integer (I4B), dimension (:,:) :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:,:) :: procs_lin ! this one will contain the nodes allocation to the group of processors consolidating the linearized operator. 
                          !it is a scrap variable and once the jobs are disctibuted, it will be deallocated...   

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nn,i,j,jj,ii,size_est,num_ranks,node ! scrap index
                 
! variables for MPI Calls
integer ::  mpicommworldsize,ierr,iprc,sendproc
integer :: MPI_COMM_WORLD_GROUP ! to keep the handle to the group of MPI_COMM_WORLD
integer, dimension(3) :: ranges ! ranges to create a communicator group
integer :: Acopy_UNIV_GROUP,Acopy_UNIV,LinFmA_UNIV_GROUP,LinFmA_UNIV ! scrap pointers to Acopy and LinFmA universes
!!  Acopy_UNIV, Acopy_UNIV_GROUP ! on the master node this is a temporary holder for the pointer to a universe, and communicatio ngroup that contains processors dealing with one copy of A operator
!!  on the slave node this is avariable for the Acopy communicator
integer, dimension (1) :: ranksII_scrap ! This is a scrap variable to check ranks of the processors in the MPI_LINEAR_COLL_WORLD communicator
integer (I4B), dimension(:), allocatable :: rashk_temp !
integer (I4B), dimension(:), allocatable :: ranks_FmA_UNIV !
integer (I4B), dimension(:), allocatable :: sortNodesProcsBuffer ! this buffer will be used for sorting nodes
integer (I4B), dimension(:), allocatable :: LinFmA_recv_univs ! This array will hold the list of MPI_LinFmA_UNIVS communicators from which it will expect to receive components of the linearized operator and the velocity nodes

!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
if (.NOT. ((size(procs_lin,1)==2) .AND. (size(procs_lin,2)==num_lin_proc))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array procs_lin. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Here we set up some variables that will be used by all processes: 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!!!!!!!!!!!!!!
!!  Next we create a communicator that will include the master node (irank =0) and all nodes with rank <=num_lin_proc-1
!!!!!!!!!!!!!!
! first we need to create a communicator group. All processes execute it:
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_WORLD_GROUP, ierr) ! get the hanlde to the group of MPI_COMM_WORLD. 
ranges(1)=0; ranges(2)=num_lin_proc-1; ranges(3)=1; 
call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, MPI_LINEAR_COLL_WORLD_GROUP, ierr) ! MPI_LINEAR_COLL_WORLD_GROUP will be dedicated to communication of 
                           ! processors with numbers 0, ..., num_lin_proc-1
                           ! num_lin_proc is set in DGVparameters.dat 
                           ! In many instances num_lin_proc is set to the number of total MPI processors requested for the run
                           ! But it can be a smalle number too. If num_lin_proc is equal to the total number of MPI processors, then 
                           ! MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_WORLD_GROUP will refer to the same processors                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group MPI_LINEAR_COLL_WORLD_GROUP on proc", irank, ". Error=", ierr
         stop
end if
! a quick sanity check if the processor with irank=0 is still with the rank =0 in the new group.
call MPI_GROUP_TRANSLATE_RANKS(MPI_LINEAR_COLL_WORLD_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: can't translate ranks MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc",irank,".Err=",ierr
 stop
end if
! a quick check that the zeros process is still number zero in the new group: 
if (ranksII_scrap(1) /= 0) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between  MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc", irank, "stop."
 stop
end if 
! end sanity check 
! next we create a new communicator for the MPI_LINEAR_COLL_WORLD_GROUP: 
call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, MPI_LINEAR_COLL_WORLD, ierr) ! MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for 
                           ! processors with numbers 0, ..., num_lin_proc-1
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator MPI_LINEAR_COLL_WORLD. Proc=", irank, ". Err=", ierr
 stop
end if
! End create communicator for consolidated handling of linearised collision operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! We now create universes to deal with the groups of processors sharing a single copy 
!!! of A and the universes of processors that obtain 
!!! components of the colliion operator from the same group of nodes... 
!!!
!!! These universies will be used for communcation of the components of the linearized operator.
!!! These universes can also be used for evaluation of the collision operator if the spatial operator is parallelized  
!!! 
!!! Begin creating the Acopy_UNIV universes:
!!! first, we allocate storage for the universes... 
 allocate (MPI_Acopy_UNIVS(1:num_Acopies),MPI_Acopy_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
    ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_Acopy_UNIVS=0;MPI_Acopy_UNIVS_GROUP=0; ! Nullify them , just in case...
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Reset the index variable 
 My_Acopy_Univ_indx = 0 ! points to no particular Acopy Universe
 ! now we will create num_Acopies universes that will contain the master node and the nodes indentified in chnks_copies(1:2,j)
 do i=1,num_Acopies
  ranges(1)=chnks_copies(1,i); ranges(2)=chnks_copies(2,i); ranges(3)=1; 
  call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, Acopy_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group Acopy_UNIVS_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! a quick check if the processor with irank=0 int he new group correspond to the first processor in the group.
  call MPI_GROUP_TRANSLATE_RANKS(Acopy_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't translate Acopy_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the zeros process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between Acopy_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, Acopy_UNIV_GROUP, Acopy_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator Acopy_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_Acopy_UNIVS(i)=Acopy_UNIV
  MPI_Acopy_UNIVS_GROUP(i)=Acopy_UNIV_GROUP
  !!!!
  ! We also need to remember to what Acopy universe we belong
  if ((irank >= chnks_copies(1,i)) .and. (irank<= chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! now we know which Acopy universe contains this slave processor. Master processor will still have zero in this variables
  end if 
 end do  
 !! END creating the Acopy_UNIV universes:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will use the costant My_Acopy_Univ_indx to set 
 ! the variable Acopy_irank which defines the rank of the processor in this A copy Universe
 !!!!!!!!!!!!!!!! 
 Acopy_irank = -1
 if (My_Acopy_Univ_indx>0) then 
 !!! MPI FORK
   Acopy_UNIV=MPI_Acopy_UNIVS(My_Acopy_Univ_indx) 
   call mpi_comm_rank(Acopy_UNIV,Acopy_irank,ierr) ! check what processor the program is on... 
   if (ierr /= 0 ) then 
    print *,"PrepMPICollsnDGV_Univs_highperf: can't determined rank in Acopy_UNIV.i=",My_Acopy_Univ_indx,"Proc=",irank,". Err=",ierr
   stop
  end if
 end if   
 !!!!!!!!!!!!!!!
 ! end setting Acopy_irank
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !! Begin creating LinFmA_UNIV universes:
 allocate (MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables MPI_LinFmA_UNIVS, MPI_LinFmA_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_LinFmA_UNIVS=0; MPI_LinFmA_UNIVS_GROUP=0; ! Nullify them, just in case ...
!! 
!! Next we need will pick one group of processors sharing a copy of A and find all processors that will 
!! receive components of the linearized operator from this group of processors.  
!! the first processors in the group and all processors receiveing compoents will be added to the 
!! communicator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Lin_Self_Send_FmA_Flag = 0 ! need to reset this flag before start...
 ! need the following temporary arrays
 allocate (rashk_temp(1:num_lin_proc+1), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables rashk_temp. stop. irank", irank 
  stop
 end if
 ! begin creating universes MPI_LinFmA_UNIVS
 do i=1,num_Acopies
  ! first we need to set up the range of ranks that will form the communicator group. 
  rashk_temp = -1  ! clean the array
  num_ranks = 1
  rashk_temp(num_ranks) = chnks_copies(1,i) ! this processor will serve as a master processor for this comminucator group.
  ! now we will start screenign the processors participating in the evaluation of the linearized operator. Recall that the 
  ! chnks_copies array contains irank values of the processors that share a copy of A. Also nodes_proc_groups contains 
  ! the first and last of velosity nodes that are assigned to each processor group (aka Acopy group).
  ! We will go over each processor evaluating linearized collision operatot and over each velocity node 
  ! that are assigned to each of them. For each velocity node, we check if that node is assigend to the 
  ! current Acopy group, which is the group containing the processor of rank rashk_temp(1). 
  ! If they do, then this processor will be needing information from that Acopy group. We, therefore, 
  ! add the irank fo the processor to the array of ranks that will comminicate to this Acopy group, the 
  ! rashk_temp() array. 
  do j=1,num_lin_proc ! j=1 points to irank=0
   do jj=procs_lin(1,j),procs_lin(2,j) ! runs over all nodes of linearizaed operator assiged to this processor 
    if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then 
      ! first we check if this this processor (must be a slave proc) is, first, a master proc of an Acopy universe 
      ! and, second, needs to get components of FmA from that universe. In this case it have dual role in FmA -- as the master and the slave -- 
      ! to avoid confusion, we do not add it to the rashk_temp since it is already there. 
      ! also, we set up a spacial flag for that 
      if (j-1 == chnks_copies(1,i)) then 
       if (irank == chnks_copies(1,i)) then 
        Lin_Self_Send_FmA_Flag = 1 !! this means that this processor (must be a slave proc) is a master proc of an Acopy universe and needs to get components of FmA from this universe 
       end if 
      else 
       num_ranks=num_ranks+1
       rashk_temp(num_ranks)=j-1  ! this will correspond to irank
      end if
      exit ! the recieving processor has been added to the ranks. No need to find more nodes that this processor recieves
    end if   
   end do 
  end do 
  ! we have temporary array of ranks, we now create the final array of ranks:
  allocate (ranks_FmA_UNIV(1:num_ranks),  stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables ranks_FmA_UNIV. stop. irank", irank 
   stop
  end if
  ranks_FmA_UNIV(:) = rashk_temp(1:num_ranks)
  ! next we create the group for the FmA_UNIV
  call MPI_GROUP_INCL(MPI_COMM_WORLD_GROUP,num_ranks,ranks_FmA_UNIV,LinFmA_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group LinFmA_UNIV_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! clean up some variables for the next iteration
  deallocate(ranks_FmA_UNIV) 
  num_ranks=0
  rashk_temp=-1
  ! a quick check if the processor with irank=chnks_copies(1,i) has rank =0 in the new group.
  call MPI_GROUP_TRANSLATE_RANKS(LinFmA_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't translate LinFmA_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between LinFmA_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, LinFmA_UNIV_GROUP, LinFmA_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator LinFmA_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_LinFmA_UNIVS(i) = LinFmA_UNIV
  MPI_LinFmA_UNIVS_GROUP(i)=LinFmA_UNIV_GROUP
 end do 
 deallocate(rashk_temp)
 !! END setting up MPI communicators MPI_LinFmA_UNIVS and MPI_Acopy_UNIVS
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEBUG 
! print *, "irank:", irank,"LinFMA_univ:", MPI_LinFmA_UNIVS!
!
!!!!!!!!!!!!!!!
!   call mpi_barrier (MPI_COMM_WORLD,ierr)
!!   if (ierr /= 0) then 
 !   print *," PrepMPICollsnDGV_Univs_highperf: MPI_Barrier from proc",irank, "returned error", ierr, "on process", irank
!    stop
!   end if 
!   stop
!!! END DEBUG   
!!!!!!!!!!!!!!!!   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! the last piece is to create the call order array LinFmA_Univs_callorder.
!! this array contains the sequence in which the LinFmA universes are called
!! It contains numbers from 1 to num_Acopies that are permuted in some way
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (LinFmA_Univs_callorder(1:num_Acopies),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "PrepMPICollsnDGV_highperf: Allocation error (LinFmA_Univs_callorder).Process",irank,"Stop"
 stop
end if
LinFmA_Univs_callorder = 0
!! Now we will populate the array LinFmA_Univs_callorder
i=1
j=1
do while (i<=num_Acopies) 
    LinFmA_Univs_callorder(i)=j
    i = i + 1
    j = j + 2 
    if (j > num_Acopies) then 
     j = 2
    end if 
end do 
!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Setting up arrays that facilitate exchange of computed values for parallelization of the 
!! linearized operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI FORK
if (irank < num_lin_proc) then 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set arrays linprocndswld_Acopy_univs,
 !! and linprocndswld_Acopy_addr
 !!
 !! We recall that these arrays only make sense on processors that are involved in the evaluation of the 
 !! linearized operator
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nn=size(nodes_u,1)
 size_est = 1 + (INT((nn/num_lin_proc),I4B)+2)*2 !this should estimate from above the size of LinFmA_recv_univs
 !  
 allocate (LinFmA_recv_univs(1:size_est), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables LinFmA_recv_univs. stop. irank", irank 
  stop
 end if
 ! clean up before proceeding
 LinFmA_recv_univs = 0   ! 0 correspond to no communicator containing a single copy of A
 do i=1,num_Acopies
  do jj=procs_lin(1,irank+1),procs_lin(2,irank+1) ! runs over all nodes of linearizaed operator assiged to this processor 
   if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then
    ! Now we fill the temporary receive array # pairs of records # A copy universe # node   # A copy universe # node ... 
    LinFmA_recv_univs(1) = LinFmA_recv_univs(1)+1     ! this will add the communicator (FmA_UNIV) to the list of communicators
    LinFmA_recv_univs(2*LinFmA_recv_univs(1)) = i   ! this records the number of A copy universe from which components of the linearized  operator are coming from
    LinFmA_recv_univs(2*LinFmA_recv_univs(1)+1) = jj ! this records the node to which the componets correspond.     
   end if 
  end do
 end do   
 !! Next we will populate the array linprocndswld_Acopy_univs. The information that will be stored in this array is already in the 
 !! array LinFmA_recv_univs. The data just needs to be re-sorted:
 !! In the buffer the first number is the number of records after that go pairs: (#processor, #node) (coment "node" means velocity node)
 !! What we do next is to rearrange this information in the following format:
 !! (#processor,#number of nodes for this processor,#node1, #node2, and so on, #nodeN) 
 !! then continue onto the next processor
 !! This is done in two steps: First, we prepare the buffer for sorting
 !! we now allocate the buffer for the linprocndswld_Acopy_univs array
 nn = LinFmA_recv_univs(1)
 allocate (sortNodesProcsBuffer(1:3*nn+1), stat=loc_alloc_stat) ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables sortNodesProcsBuffer on proc",irank," stop"
  stop
 end if
 !! now we are calling the sorting routine. The return variable sortNodesProcsBuffer
 !! will contain sorted records
 call SortNodesProcs(LinFmA_recv_univs,sortNodesProcsBuffer) ! the first element in buffer contains the number of pair. The first element on return array contains the number of non-trivial records.
 !! now we create the array and copy the records into it
 nn = sortNodesProcsBuffer(1) ! the first entry is the number of records in the array except for the itself
 ! now we need to create the local send_lin_procs array. 
 allocate (linprocndswld_Acopy_univs(1:nn), linprocndswld_Acopy_addr(1:nn),linprocndswld_BfmA_addr(1:nn),stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error (linprocndswld_Acopy_univs/_Acopy_addr/_BfmA_addr). Process",irank,"Stop"
  stop
 end if
 ! next we populate the array linprocndswld_Acopy_univs ... 
 linprocndswld_Acopy_univs = sortNodesProcsBuffer(2:nn+1)
 ! next we need to populate array linprocndswld_Acopy_addr 
 linprocndswld_Acopy_addr=0 ! reset, just in case
 i=1
 do while (i<size(linprocndswld_Acopy_univs,1))
  linprocndswld_Acopy_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
  linprocndswld_Acopy_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
  do ii=1,linprocndswld_Acopy_univs(i+1)
   j = linprocndswld_Acopy_univs(i)         !!! the index of the Acopy universe
   node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
   linprocndswld_Acopy_addr(i+1+ii) = node - nodes_proc_groups(1,j) + 1  !!!! Attention -- we do not really look into  Acopy_Wrkld, but we look into the array that is used to create Acopy_Wrkld
  end do
  i=i+2+linprocndswld_Acopy_addr(i+1)  
 end do 
 !! the linprocndswld_Acopy_addr arrays are set! 
 ! next we need to populate array linprocndswld_BfmA_addr that contains the addresses on the nodes listed in the 
 ! linprocndswld_Acopy_univs as they appear in the consolidated linearied operator. 
 linprocndswld_BfmA_addr = 0 ! reset, just in case
 i=1
 do while (i<size(linprocndswld_BfmA_addr,1))
  linprocndswld_BfmA_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
  linprocndswld_BfmA_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
  do ii=1,linprocndswld_Acopy_univs(i+1)
   node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
   linprocndswld_BfmA_addr(i+1+ii) = node - procs_lin(1,irank+1) + 1  !!!! Attention -- we do not really look into  lin_proc_nodes_wld, but we look into the array that is used to create lin_proc_nodes_wld
  end do
  i=i+2+linprocndswld_BfmA_addr(i+1)  
 end do 
 !! the linprocndswld_BfmA_addr arrays are set! 
end if 
!! clean up 
deallocate(LinFmA_recv_univs,sortNodesProcsBuffer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end subroutine  PrepMPICollsnDGV_Univs_highperf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! PrepMPICollsnDGV_Univs_highperf_SelNds
 !
 ! This is a copy of the above subroutine except it is designed to use only a part of the primary velocity nodes. 
 ! Total number of the used nodes is passed as a variable. The nodes listed in the array, say sel_nodes, and will be are divided betwen processors 
 ! the numbers in the arrays nodes_proc_groups,procs_lin correspond to index in array sel_nodes and entries of sel_nodes are the numbers of the 
 ! node in nodes_u/_v/_w etc. arrays.
 ! Note that we do not need to see the arrya of nodes, we lust need to know how many are there, so we do not pass the array sel_nodes itself, just the number of records in it.
 !
 !
 !
 ! The next subroutines are chunks of PrepMPICollsnDGV_highperf. The splitting alows for 
 ! more flexibility in customizing MPI codes developed using the DGv library
 !
 ! This piece creates communication universes for primary velocity mesh. 
 ! 
 !
 ! This subroutines contains an MPI fork 
 ! This subroutine is run on all processors
 !
 ! Thi subroutine prepares MPI environments and data arrays that will be 
 ! used for MPI parallelization for 
 !    (1) evaluation of the Boltzmann collision operator
 !    (2) can be used for parallel evaluation of the linearized Boltzmann collision operator
 !    (3) can be used with a spatial driver to parallilize in velocity space. Specifically, nodes of the 
 !         velocity mesh are divided between MPI communicating processors. Each processor is running 
 !         a copy of spatial operator for a portion of velocity nodes, the processor request values of the collision 
 !         operator for their portion of nodes. 
 ! 
 ! Description: 
 !
 ! This subroutine prepares a number of communicators. 
 !
 !  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
 !
 !  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
 !                                  sharing a copy of A
 ! 
 !  MPI_Acopy_UNIVS(1:num_Acopies), MPI_Acopy_UNIVS_GROUP(1:num_Acopies) - pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors sharing a copy of A
 !  
 !  LinFmA_UNIV_GROUP, LinFmA_UNIV, -- these are scrap variable used to create a universe for a group of processors 
 !                                  that receives information from a single group of processors sharing a copy of A
 !                                  Here is some idea on how these communicators work: 
 !										There is exactly num_Acopies of groups of processors. Let's sasy its number is 0 <= I <=num_Acopies.
 !										Each of the goups has a workload, given by (nodes_proc_groups(1,I), \ldots, nodes_proc_groups(2,I)), 
 !										that is a number of velocity points 
 !										at which this group of processors I evaluates the collision operator.
 !                                      Just in case, this group of processors is included in the communicator MPI_Acopy_UNIVS(I), but that is 
 !										not important at the moment. Values of the collision operator that are computed in the group I   
 !										may be neded on other processors for either evaluation of 
 !										linearized collision operator or for parallel evaluation of the kinetic equations using 
 !										parallelization in the velocity nodes. This is determined by examining array procs_lim, that 
 !										a range of velocity points to each participating processor. If a processor needs a value of the collision 
 !										operator computed in the group I, then the processor 
 !                                      is included in the LinFmA_Univ/_Group. Also, the master node of the MPI_Acopy_UNIVS(I) is included 
 !										included in the LinFmA_Univ(I) (and also made a master node in it?) 
 !										Typically, the master node of MPI_Acopy_UNIVS(I) broadcasts values to all of LinFmA_Univ(I)
 !
 !  MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP -- pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors that receives information 
 !                                  from a single group of processors sharing a copy of A
 !
 !
 !
 !!! MPI FORK
 !!! Now the master process will distribute the chunks of A arrays to other processors. If the  
 !!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
 !!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
 !!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
 !!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
 !!! Therefore, the results will be combined in an integer of (I8B) kind. 
 !!!
 !!! This number will be divided by the number of 
 !!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
 !!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
 !!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
 !!!  
 !!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
 !!! This boradcast will go out with the code 30.
 !!!
 !!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
 !!! 
 !!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
 !!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
 !!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
 !!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
 !!! 
 
 !!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
 !!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
 !!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
 !!!? workload. 
 !!! 
 !!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
 !!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
 !!!
 !!! Next the master processor initiates the work: It simply broadcasts the solution, 
 !!! 
 !!!
 !!! The slave processes will receive the solutions and will start to evaluate their portions of the collision operator. 
 !!! Then they will send their results to the master processor using an aggregate reduce operation.
 !!!
 !!! Finally, let us briefly describe the process of calculating the collisoin opperator by slave 
 !!! processors. once the I1workload is set,  the slave processor will wait for a boradcast of solution. 
 !!! Each slave process then goes over the I1workload and performs the 
 !!! evaluation of the collision integral for each I1 listed there. The results are stored in the array (CollisionInt)
 !!! Once all I1 are processed, this array is aggregately sent to the master node. Operation of the summation is used to combine 
 !!! results for all chunks. 
 !!! 
 !!! Also, the slave process will evaluate the linearized part
 !!! also, sometimes the slave processor will send the linearized part to some other processor for faster assembly...
 !!!
 !!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  PrepMPICollsnDGV_Univs_highperf_SelNds(irank,chnks_copies,nodes_proc_groups,procs_lin,num_sel_nodes)

use DGV_commvar, only: num_lin_proc, nodes_u, num_Acopies, &
                       MPI_LINEAR_COLL_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, &
                       MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP, My_Acopy_Univ_indx, Acopy_irank, &
                       MPI_LinFmA_UNIVS_GROUP, MPI_LinFmA_UNIVS, Lin_Self_Send_FmA_Flag, &
                       LinFmA_Univs_callorder,&
                       linprocndswld_Acopy_univs,linprocndswld_Acopy_addr,linprocndswld_BfmA_addr
use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:) :: chnks_copies !
integer (I4B), dimension (:,:) :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:,:) :: procs_lin ! this one will contain the nodes allocation to the group of processors consolidating the linearized operator. 
                          !it is a scrap variable and once the jobs are disctibuted, it will be deallocated...   
integer (I4B), intent(in)  :: num_sel_nodes ! this is the length of the array (sel_nodes, for instance) that contains indices of velocity nodes that are 
                     ! used in the simulations. Arrays nodes_proc_groups 
                     ! and procs_lin contain values of the indices of the array sel_nodes.                            

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nn,i,j,jj,ii,size_est,num_ranks,node ! scrap index
                 
! variables for MPI Calls
integer ::  mpicommworldsize,ierr,iprc,sendproc
integer :: MPI_COMM_WORLD_GROUP ! to keep the handle to the group of MPI_COMM_WORLD
integer, dimension(3) :: ranges ! ranges to create a communicator group
integer :: Acopy_UNIV_GROUP,Acopy_UNIV,LinFmA_UNIV_GROUP,LinFmA_UNIV ! scrap pointers to Acopy and LinFmA universes
!!  Acopy_UNIV, Acopy_UNIV_GROUP ! on the master node this is a temporary holder for the pointer to a universe, and communicatio ngroup that contains processors dealing with one copy of A operator
!!  on the slave node this is avariable for the Acopy communicator
integer, dimension (1) :: ranksII_scrap ! This is a scrap variable to check ranks of the processors in the MPI_LINEAR_COLL_WORLD communicator
integer (I4B), dimension(:), allocatable :: rashk_temp !
integer (I4B), dimension(:), allocatable :: ranks_FmA_UNIV !
integer (I4B), dimension(:), allocatable :: sortNodesProcsBuffer ! this buffer will be used for sorting nodes
integer (I4B), dimension(:), allocatable :: LinFmA_recv_univs ! This array will hold the list of MPI_LinFmA_UNIVS communicators from which it will expect to receive components of the linearized operator and the velocity nodes

!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
if (.NOT. ((size(procs_lin,1)==2) .AND. (size(procs_lin,2)==num_lin_proc))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array procs_lin. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Here we set up some variables that will be used by all processes: 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!!!!!!!!!!!!!!
!!  Next we create a communicator that will include the master node (irank =0) and all nodes with rank <=num_lin_proc-1
!!!!!!!!!!!!!!
! first we need to create a communicator group. All processes execute it:
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_WORLD_GROUP, ierr) ! get the hanlde to the group of MPI_COMM_WORLD. 
ranges(1)=0; ranges(2)=num_lin_proc-1; ranges(3)=1; 
call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, MPI_LINEAR_COLL_WORLD_GROUP, ierr) ! MPI_LINEAR_COLL_WORLD_GROUP will be dedicated to communication of 
                           ! processors with numbers 0, ..., num_lin_proc-1
                           ! num_lin_proc is set in DGVparameters.dat 
                           ! In many instances num_lin_proc is set to the number of total MPI processors requested for the run
                           ! But it can be a smalle number too. If num_lin_proc is equal to the total number of MPI processors, then 
                           ! MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_WORLD_GROUP will refer to the same processors                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group MPI_LINEAR_COLL_WORLD_GROUP on proc", irank, ". Error=", ierr
         stop
end if
! a quick sanity check if the processor with irank=0 is still with the rank =0 in the new group.
call MPI_GROUP_TRANSLATE_RANKS(MPI_LINEAR_COLL_WORLD_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: can't translate ranks MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc",irank,".Err=",ierr
 stop
end if
! a quick check that the zeros process is still number zero in the new group: 
if (ranksII_scrap(1) /= 0) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between  MPI_LINEAR_COLL_WORLD_GROUP and MPI_COMM_GROUP proc", irank, "stop."
 stop
end if 
! end sanity check 
! next we create a new communicator for the MPI_LINEAR_COLL_WORLD_GROUP: 
call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, MPI_LINEAR_COLL_WORLD, ierr) ! MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for 
                           ! processors with numbers 0, ..., num_lin_proc-1
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator MPI_LINEAR_COLL_WORLD. Proc=", irank, ". Err=", ierr
 stop
end if
! End create communicator for consolidated handling of linearised collision operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! We now create universes to deal with the groups of processors sharing a single copy 
!!! of A and the universes of processors that obtain 
!!! components of the colliion operator from the same group of nodes... 
!!!
!!! These universies will be used for communcation of the components of the linearized operator.
!!! These universes can also be used for evaluation of the collision operator if the spatial operator is parallelized  
!!! 
!!! Begin creating the Acopy_UNIV universes:
!!! first, we allocate storage for the universes... 
 allocate (MPI_Acopy_UNIVS(1:num_Acopies),MPI_Acopy_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
    ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_Acopy_UNIVS=0;MPI_Acopy_UNIVS_GROUP=0; ! Nullify them , just in case...
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Reset the index variable 
 My_Acopy_Univ_indx = 0 ! points to no particular Acopy Universe
 ! now we will create num_Acopies universes that will contain the master node and the nodes indentified in chnks_copies(1:2,j)
 do i=1,num_Acopies
  ranges(1)=chnks_copies(1,i); ranges(2)=chnks_copies(2,i); ranges(3)=1; 
  call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, Acopy_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group Acopy_UNIVS_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! a quick check if the processor with irank=0 int he new group correspond to the first processor in the group.
  call MPI_GROUP_TRANSLATE_RANKS(Acopy_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't translate Acopy_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the zeros process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between Acopy_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, Acopy_UNIV_GROUP, Acopy_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator Acopy_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_Acopy_UNIVS(i)=Acopy_UNIV
  MPI_Acopy_UNIVS_GROUP(i)=Acopy_UNIV_GROUP
  !!!!
  ! We also need to remember to what Acopy universe we belong
  if ((irank >= chnks_copies(1,i)) .and. (irank<= chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! now we know which Acopy universe contains this slave processor. Master processor will still have zero in this variables
  end if 
 end do  
 !! END creating the Acopy_UNIV universes:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will use the costant My_Acopy_Univ_indx to set 
 ! the variable Acopy_irank which defines the rank of the processor in this A copy Universe
 !!!!!!!!!!!!!!!! 
 Acopy_irank = -1
 if (My_Acopy_Univ_indx>0) then 
 !!! MPI FORK
   Acopy_UNIV=MPI_Acopy_UNIVS(My_Acopy_Univ_indx) 
   call mpi_comm_rank(Acopy_UNIV,Acopy_irank,ierr) ! check what processor the program is on... 
   if (ierr /= 0 ) then 
    print *,"PrepMPICollsnDGV_Univs_highperf: can't determined rank in Acopy_UNIV.i=",My_Acopy_Univ_indx,"Proc=",irank,". Err=",ierr
   stop
  end if
 end if   
 !!!!!!!!!!!!!!!
 ! end setting Acopy_irank
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !! Begin creating LinFmA_UNIV universes:
 allocate (MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables MPI_LinFmA_UNIVS, MPI_LinFmA_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_LinFmA_UNIVS=0; MPI_LinFmA_UNIVS_GROUP=0; ! Nullify them, just in case ...
!! 
!! Next we need will pick one group of processors sharing a copy of A and find all processors that will 
!! receive components of the linearized operator from this group of processors.  
!! the first processors in the group and all processors receiveing compoents will be added to the 
!! communicator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 Lin_Self_Send_FmA_Flag = 0 ! need to reset this flag before start...
 ! need the following temporary arrays
 allocate (rashk_temp(1:num_lin_proc+1), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables rashk_temp. stop. irank", irank 
  stop
 end if
 ! begin creating universes MPI_LinFmA_UNIVS
 do i=1,num_Acopies
  ! first we need to set up the range of ranks that will form the communicator group. 
  rashk_temp = -1  ! clean the array
  num_ranks = 1
  rashk_temp(num_ranks) = chnks_copies(1,i) ! this processor will serve as a master processor for this comminucator group.
  ! now we will start screenign the processors participating in the evaluation of the linearized operator. Recall that the 
  ! chnks_copies array contains irank values of the processors that share a copy of A. Also nodes_proc_groups contains 
  ! the first and last of velosity nodes that are assigned to each processor group (aka Acopy group).
  ! We will go over each processor evaluating linearized collision operatot and over each velocity node 
  ! that are assigned to each of them. For each velocity node, we check if that node is assigend to the 
  ! current Acopy group, which is the group containing the processor of rank rashk_temp(1). 
  ! If they do, then this processor will be needing information from that Acopy group. We, therefore, 
  ! add the irank fo the processor to the array of ranks that will comminicate to this Acopy group, the 
  ! rashk_temp() array. 
  do j=1,num_lin_proc ! j=1 points to irank=0
   do jj=procs_lin(1,j),procs_lin(2,j) ! runs over all nodes of linearizaed operator assiged to this processor 
    if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then 
      ! first we check if this this processor (must be a slave proc) is, first, a master proc of an Acopy universe 
      ! and, second, needs to get components of FmA from that universe. In this case it have dual role in FmA -- as the master and the slave -- 
      ! to avoid confusion, we do not add it to the rashk_temp since it is already there. 
      ! also, we set up a spacial flag for that 
      if (j-1 == chnks_copies(1,i)) then 
       if (irank == chnks_copies(1,i)) then 
        Lin_Self_Send_FmA_Flag = 1 !! this means that this processor (must be a slave proc) is a master proc of an Acopy universe and needs to get components of FmA from this universe 
       end if 
      else 
       num_ranks=num_ranks+1
       rashk_temp(num_ranks)=j-1  ! this will correspond to irank
      end if
      exit ! the recieving processor has been added to the ranks. No need to find more nodes that this processor recieves
    end if   
   end do 
  end do 
  ! we have temporary array of ranks, we now create the final array of ranks:
  allocate (ranks_FmA_UNIV(1:num_ranks),  stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables ranks_FmA_UNIV. stop. irank", irank 
   stop
  end if
  ranks_FmA_UNIV(:) = rashk_temp(1:num_ranks)
  ! next we create the group for the FmA_UNIV
  call MPI_GROUP_INCL(MPI_COMM_WORLD_GROUP,num_ranks,ranks_FmA_UNIV,LinFmA_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: MPI can not create group LinFmA_UNIV_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! clean up some variables for the next iteration
  deallocate(ranks_FmA_UNIV) 
  num_ranks=0
  rashk_temp=-1
  ! a quick check if the processor with irank=chnks_copies(1,i) has rank =0 in the new group.
  call MPI_GROUP_TRANSLATE_RANKS(LinFmA_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't translate LinFmA_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: ranks mixed up between LinFmA_UNIV_GROUP and MPI_COMM_GROUP i=",i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, LinFmA_UNIV_GROUP, LinFmA_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperf: can't create communicator LinFmA_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_LinFmA_UNIVS(i) = LinFmA_UNIV
  MPI_LinFmA_UNIVS_GROUP(i)=LinFmA_UNIV_GROUP
 end do 
 deallocate(rashk_temp)
 !! END setting up MPI communicators MPI_LinFmA_UNIVS and MPI_Acopy_UNIVS
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DEBUG 
! print *, "irank:", irank,"LinFMA_univ:", MPI_LinFmA_UNIVS!
!
!!!!!!!!!!!!!!!
!   call mpi_barrier (MPI_COMM_WORLD,ierr)
!!   if (ierr /= 0) then 
 !   print *," PrepMPICollsnDGV_Univs_highperf: MPI_Barrier from proc",irank, "returned error", ierr, "on process", irank
!    stop
!   end if 
!   stop
!!! END DEBUG   
!!!!!!!!!!!!!!!!   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! the last piece is to create the call order array LinFmA_Univs_callorder.
!! this array contains the sequence in which the LinFmA universes are called
!! It contains numbers from 1 to num_Acopies that are permuted in some way
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (LinFmA_Univs_callorder(1:num_Acopies),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "PrepMPICollsnDGV_highperf: Allocation error (LinFmA_Univs_callorder).Process",irank,"Stop"
 stop
end if
LinFmA_Univs_callorder = 0
!! Now we will populate the array LinFmA_Univs_callorder
i=1
j=1
do while (i<=num_Acopies) 
    LinFmA_Univs_callorder(i)=j
    i = i + 1
    j = j + 2 
    if (j > num_Acopies) then 
     j = 2
    end if 
end do 
!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Setting up arrays that facilitate exchange of computed values for parallelization of the 
!! linearized operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI FORK
if (irank < num_lin_proc) then 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set arrays linprocndswld_Acopy_univs,
 !! and linprocndswld_Acopy_addr
 !!
 !! We recall that these arrays only make sense on processors that are involved in the evaluation of the 
 !! linearized operator
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nn=Num_sel_nodes
 size_est = 1 + (INT((nn/num_lin_proc),I4B)+2)*2 !this should estimate from above the size of LinFmA_recv_univs
 !  
 allocate (LinFmA_recv_univs(1:size_est), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperf: Allocation error for variables LinFmA_recv_univs. stop. irank", irank 
  stop
 end if
 ! clean up before proceeding
 LinFmA_recv_univs = 0   ! 0 correspond to no communicator containing a single copy of A
 do i=1,num_Acopies
  do jj=procs_lin(1,irank+1),procs_lin(2,irank+1) ! runs over all nodes of linearizaed operator assiged to this processor 
   if ((jj >= nodes_proc_groups(1,i)) .and. (jj<= nodes_proc_groups(2,i))) then
    ! Now we fill the temporary receive array # pairs of records # A copy universe # node   # A copy universe # node ... 
    LinFmA_recv_univs(1) = LinFmA_recv_univs(1)+1     ! this will add the communicator (FmA_UNIV) to the list of communicators
    LinFmA_recv_univs(2*LinFmA_recv_univs(1)) = i   ! this records the number of A copy universe from which components of the linearized  operator are coming from
    LinFmA_recv_univs(2*LinFmA_recv_univs(1)+1) = jj ! this records the index of the node in sel_nds array to which the componets correspond .     
   end if 
  end do
 end do   
 !! Next we will populate the array linprocndswld_Acopy_univs. The information that will be stored in this array is already in the 
 !! array LinFmA_recv_univs. The data just needs to be re-sorted:
 !! In the buffer the first number is the number of records after that go pairs: (#processor, #node) (coment "node" means velocity node)
 !! What we do next is to rearrange this information in the following format:
 !! (#processor,#number of nodes for this processor,#node1, #node2, and so on, #nodeN) 
 !! then continue onto the next processor
 !! This is done in two steps: First, we prepare the buffer for sorting
 !! we now allocate the buffer for the linprocndswld_Acopy_univs array
 nn = LinFmA_recv_univs(1)
 allocate (sortNodesProcsBuffer(1:3*nn+1), stat=loc_alloc_stat) ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables sortNodesProcsBuffer on proc",irank," stop"
  stop
 end if
 !! now we are calling the sorting routine. The return variable sortNodesProcsBuffer
 !! will contain sorted records
 call SortNodesProcs(LinFmA_recv_univs,sortNodesProcsBuffer) ! the first element in buffer contains the number of pair. The first element on return array contains the number of non-trivial records.
 !! now we create the array and copy the records into it
 nn = sortNodesProcsBuffer(1) ! the first entry is the number of records in the array except for the itself
 ! now we need to create the local send_lin_procs array. 
 allocate (linprocndswld_Acopy_univs(1:nn), linprocndswld_Acopy_addr(1:nn),linprocndswld_BfmA_addr(1:nn),stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error (linprocndswld_Acopy_univs/_Acopy_addr/_BfmA_addr). Process",irank,"Stop"
  stop
 end if
 ! next we populate the array linprocndswld_Acopy_univs ... 
 linprocndswld_Acopy_univs = sortNodesProcsBuffer(2:nn+1)
 ! next we need to populate array linprocndswld_Acopy_addr 
 linprocndswld_Acopy_addr=0 ! reset, just in case
 i=1
 do while (i<size(linprocndswld_Acopy_univs,1))
  linprocndswld_Acopy_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
  linprocndswld_Acopy_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
  do ii=1,linprocndswld_Acopy_univs(i+1)
   j = linprocndswld_Acopy_univs(i)         !!! the index of the Acopy universe
   node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
   linprocndswld_Acopy_addr(i+1+ii) = node - nodes_proc_groups(1,j) + 1  !!!! Attention -- this numbering is relative the sel_nodes array
                                             !!!! attention: we do not really use Acopy_Wrkld, but expect this to be consistent with it. 
  end do                                            
  i=i+2+linprocndswld_Acopy_addr(i+1)  
 end do 
 !! the linprocndswld_Acopy_addr arrays are set! 
 ! next we need to populate array linprocndswld_BfmA_addr that contains the addresses on the nodes listed in the 
 ! linprocndswld_Acopy_univs as they appear in the consolidated linearied operator. 
 linprocndswld_BfmA_addr = 0 ! reset, just in case
 i=1
 do while (i<size(linprocndswld_BfmA_addr,1))
  linprocndswld_BfmA_addr(i) = linprocndswld_Acopy_univs(i) ! the first Acopy universe from where components will arrive
  linprocndswld_BfmA_addr(i+1) = linprocndswld_Acopy_univs(i+1) ! the number of nodes for which the components will arrive
  do ii=1,linprocndswld_Acopy_univs(i+1)
   node = linprocndswld_Acopy_univs(i+1+ii) !!!!!  node for which we determine the address in corresponding copy of Acopy_Wrkld ...
   linprocndswld_BfmA_addr(i+1+ii) = node - procs_lin(1,irank+1) + 1  !!!! Attention -- we do not really look into  lin_proc_nodes_wld, but we look into the array that is used to create lin_proc_nodes_wld
  end do
  i=i+2+linprocndswld_BfmA_addr(i+1)  
 end do 
 !! the linprocndswld_BfmA_addr arrays are set! 
end if 
!! clean up 
deallocate(LinFmA_recv_univs,sortNodesProcsBuffer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end subroutine  PrepMPICollsnDGV_Univs_highperf_SelNds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! PrepMPICollsnDGV_ReadOpA_highperf(irank)
!
! This subroutines is a chunk of PrepMPICollsnDGV_highperf. The splitting alows for 
! more flexibility in customizing MPI codes developed using the DGv library
!
! This piece distributes the A operator between the processors that evaluate the collision operator 
!
!!!!!!!!
! ATTENTION! Must be called after PrepMPICollsnDGV_Univs_highperf or will not work correctly
!!!!!!!!
! This subroutines contains an MPI fork 
! This subroutine is run on all processors
!
! Thi subroutine prepares MPI environments and data arrays that will be 
! used for MPI parallelization for 
!    (1) evaluation of the Boltzmann collision operator
!    (2) can be used for parallel evaluation of the linearized Boltzmann collision operator
!    (3) can be used with a spatial driver to parallilize in velocity space. Specifically, nodes of the 
!         velocity mesh are divided between MPI communicating processors. Each processor is running 
!         a copy of spatial operator for a portion of velocity nodes, the processor request values of the collision 
!         operator for their portion of nodes. 
! 
! Description: 
!
! This subroutine prepares a number of communicators. 
!
!  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
!
!  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
!                                  sharing a copy of A
! 
!  MPI_Acopy_UNIVS(1:num_Acopies), MPI_Acopy_UNIVS_GROUP(1:num_Acopies) - pointer to universes are stored in these arrays.
!                                  each pointer is pointing to a group of proecessors sharing a copy of A
!  
!  LinFmA_UNIV_GROUP, LinFmA_UNIV, -- these are scrap variable used to create a universe for a group of processors 
!                                  that receives information from a single group of processors sharing a copy of A
!                                  Here is some idea on how these communicators work: 
!										There is exactly num_Acopies of groups of processors. Let's sasy its number is 0 <= I <=num_Acopies.
!										Each of the goups has a workload, given by (nodes_proc_groups(1,I), \ldots, nodes_proc_groups(2,I)), 
!										that is a number of velocity points 
!										at which this group of processors I evaluates the collision operator.
!                                      Just in case, this group of processors is included in the communicator MPI_Acopy_UNIVS(I), but that is 
!										not important at the moment. Values of the collision operator that are computed in the group I   
!										may be neded on other processors for either evaluation of 
!										linearized collision operator or for parallel evaluation of the kinetic equations using 
!										parallelization in the velocity nodes. This is determined by examining array procs_lim, that 
!										a range of velocity points to each participating processor. If a processor needs a value of the collision 
!										operator computed in the group I, then the processor 
!                                      is included in the LinFmA_Univ/_Group. Also, the master node of the MPI_Acopy_UNIVS(I) is included 
!										included in the LinFmA_Univ(I) (and also made a master node in it?) 
!										Typically, the master node of MPI_Acopy_UNIVS(I) broadcasts values to all of LinFmA_Univ(I)
!
!  MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP -- pointer to universes are stored in these arrays.
!                                  each pointer is pointing to a group of proecessors that receives information 
!                                  from a single group of processors sharing a copy of A
!
!
!
!!! MPI FORK
!!! Now the master process will distribute the chunks of A arrays to other processors. If the  
!!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
!!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
!!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
!!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
!!! Therefore, the results will be combined in an integer of (I8B) kind. 
!!!
!!! This number will be divided by the number of 
!!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
!!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
!!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
!!!  
!!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
!!! This boradcast will go out with the code 30.
!!!
!!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
!!! 
!!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
!!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
!!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
!!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
!!! 
 
!!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
!!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
!!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
!!!? workload. 
!!! 
!!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
!!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
!!!
!!! Next the master processor initiates the work: It simply broadcasts the solution, 
!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! STEP 3.
 !!
 !! PARALLELIZATION OF THE LINEARIZATION STAFF 
 !!
 !!
 !! This part of the code replaces the original verison that is using 
 !! send_lin_procs and recv_lin_procs arrays for the construciton of the linearized collision operator
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! ATTENTION: OUTDATED!!!!
 !! The new algorithm for evaluatin the collision operator is as follows:
 !!
 !! A. Slave Processors  compute the components of the linearized operator.
 !! B. Slave processors place these componets in the copy of an array 
 !!    FmA whose rows go over all nodes in the array Acopy_Wrkld that are the same for 
 !!    all processors in the same Acopy universe
 !!    to help them stack the components the array procs_nodes_wkld_Acopy_Addr array us used. 
 !!    procs_nodes_wld_Acopy_addr has the same amount of elements as procs_nodes_wld
 !!    The array procs_nodes_wld_Acopy_addr contains indices of the nodes from 
 !!    procs_nodes_wld in Acopy_Wrkld
 !!    The lead processor of Acopy univ. gathers the FmA to it
 !! 
 !! C. The lead proc of Acopy universe is also the lead processor of the matching LinFmA universe 
 !!    This processor broadcasts FmA to processors involved in the valuation of the collision oeprator. 
 !!    For each Acopy universe there is a LinFmA_Univ universe that combines all processors that need 
 !!    to receive data from processors  invloved in Acopy_univ. 
 !! D. The processors pick the data from the boradsact and build the local copy of the linearized operator.
 !!    To do that three arrays are used:
 !!    lin_proc_nodes_wld --- contains the componets of the linearized operator assigned to each processor involved 
 !!    in the evaluation of the linearized operator. 
 !!    linprocndswld_Acopy_univs, --- this array contains informaion about arriving components of the linearized operator:
 !!     # of the universe A copy that will send the components
 !!     # total number of the nodes (each node is one row of the lineaired collision operator
 !!     # Nodes .... 
 !!    then repeat
 !!    linprocndswld_Acopy_addr has the same stucture except instead of nodes, there will be their indices in the corresponding 
 !!    Acopy_Wrkld array
 !!
 !! E. The processors involced in the evaluation of the linearized colliion operator create a copy of BFmA where they store the 
 !!    componets of the consolidated operator.  
 !!    the processors evaluate the parts of the consolidated collision operator and stack the results in the temporary solution 
 !!    array. Finally master processrs gathers all componets of the temp. solution array by a collective communication call. 
 !!!! END OUTDATED  !!
 !!!!!!!!!!!!!!!!!!!!!


!!!!!!!!
! ATTENTION! The subroutine assumes that PrepMPICollsnDGV_Univs_highperf has been called already and important 
!            variable have been set up. Will not work correctly is called before PrepMPICollsnDGV_Univs_highperf
!!!!!!!!

subroutine PrepMPICollsnDGV_ReadOpA_highperf(irank,chnks_copies,nodes_proc_groups)

use DGV_commvar, only: numchnks,num_Acopies,nodes_u,A_capphi

use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:) :: chnks_copies !
integer (I4B), dimension (:,:) :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 

integer (I4B), dimension (:), allocatable :: chnks_Act ! This is essentially a scrap array that keeps the information about number of records or A operator in each chunk of Aarrays.
integer (I4B) :: i,jj,j,ofst,leftover,nn,fstchnk,numrec   ! Scrap index
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I8B) :: Act_aggreg ! variable to keep the total number of records in all chunks of operator A

!!!!! MPI variablessortNodesProcsBuffer
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer     ! this buffer will be used in sending and receiving the work load.
integer (I4B), dimension (:,:), allocatable :: mpiCollisionWorkAlloc ! This array will contain info about work allocations for slave processors  
                                                                ! The array is operating like the following: the first index is the number of processor
                                                                ! the second index has 3 fixed fields 
                                                                ! record 1: number of first chunk to read
                                                                ! record 2: offset in the chunks to read
                                                                ! record 3: number of records to read
integer (I4B) ::  mpiAload, mpiAloadrem, mpiAload_temp  ! Scrap variables to store the averaghe number of records in A array that will be stored  in the 
                                             ! each slave processor and the reminder of the division of the total number of records in A by the number of 
                                             ! processors     
                                             
! variables for MPI Calls
integer :: ierr, mpicommworldsize

!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! we now allocate the buffer for the mpiCollisionWorkAlloc array
allocate (mpiI4BArryBuffer(1:3*(mpicommworldsize-1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables mpiI4BArryBuffer"
 stop
end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Next we need to distribute the workload to the slave processors for 
 ! evaluating the full collision operator. The starting poiint is to divide pre=-computed operator A between the processors.
 !!!!!!!!!!!!!!!!!
 !!  MPI fork
 if (irank==0) then
  !! STEP 1:
  !! Master process code: 
  !! The Master node needs to disctribute operator A. Therefore it will scall all chunks of all arrays and will 
  !! retrieve the total number of the recods in operator A. Then this number will be divided to the number of processors 
  !! to produce the operator A workloads (+/- 1)
  !! begin STEP 1:
  allocate (chnks_Act(1:numchnks), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables chnks_Act"
   stop
  end if
  !! The next subroutine will scan all chunks of Aarrays and will compile the information about the number of records in 
  !! each chunk in array chnks_Act
  ! call next subroutine to populate chnks_Act
  call ScanAarrsChnks4Acapphi(chnks_Act) ! Act_aggreg will have the total number of records in the Aarrays.
  ! Let us compute the total number of records
  Act_aggreg=0
  do i=1,numchnks
   Act_aggreg=Act_aggreg+chnks_Act(i)
  end do
  !! now we broadcast the array to all processors. 
  !! Now it is time to set up the work allocation array.
  allocate (mpiCollisionWorkAlloc(1:3,1:mpicommworldsize-1), stat=loc_alloc_stat)
     ! 
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables mpiCollisionWorkAlloc"
   stop
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mpiCollisionWorkAlloc = 0 ! reset the array
  do jj=1,num_Acopies
   !! Next we divide the total number of records by the total number of available processors to determine the workload per processor
   !! We perform division with reminder
   mpiAload = INT(Act_aggreg/(chnks_copies(2,jj)-chnks_copies(1,jj)+1),I4B)
   mpiAloadrem = INT((Act_aggreg - (INT8(chnks_copies(2,jj)-chnks_copies(1,jj)+1)*mpiAload)), I4B) ! notice that mpiAloarrem < mpicommworldsize-1
   !! A quick check that we have enough memory per processor. If we do not, stop:
   if (mpiAload > 2*10**8) then 
    print *,"PrepMPICollsnDGV_ReadOpA_highperf: Error> the number of records of operator A exceeds 2*10^8. Aborting the solver."
    stop
   end if 
   !! Next we populate the work allocation array.
   !! There will be two cases. Case 1 -- number of slave processors= number of chunks 
   !! and Case 2: number of processors neq number of chunks ... 
   if (numchnks == (chnks_copies(2,jj)-chnks_copies(1,jj)+1)) then 
   !! in the case 1 the chunks are assigned to the processor withour rechunking. 
   !! SPECIAL CASE :: NUMBER OF CHUNKS AND NUMBER OF PROCESSORS IS THE SAME! !! IMPORTANT!!
    do i = chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = i - chnks_copies(1,jj) + 1
     mpiCollisionWorkAlloc(2,i) = 0
     mpiCollisionWorkAlloc(3,i) = chnks_Act(i - chnks_copies(1,jj) + 1)
    end do
   else 
   !! in the case 2, the following subroutine determines re-chunking the A-arrays. 
    j=1 ! this is the index that will run over chunks
    ofst=0 ! This is the offset calculator
    leftover=0
    do i=chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = j  
     mpiCollisionWorkAlloc(2,i) = ofst
     if (i <= mpiAloadrem) then
      mpiAload_temp = mpiAload+1
     else 
      mpiAload_temp = mpiAload
     end if   
     do while (mpiCollisionWorkAlloc(3,i) < mpiAload_temp) ! process chunks until i-th processor is loaded. 
      ! a quick check if data is corrupted. if our calculations above arte correct, then there will be exactly as many 
      ! chunks as needed. If we ran oput of chunks -- somewhere there is an error. Then we stop
      if (j > numchnks) then 
       print *, "PrepMPICollsnDGV_ReadOpA_highperf: error when trying to populate mpiCollisionWorkAlloc. Ran out of chunks."
       stop
      end if
      ! end error check  
      if ((chnks_Act(j) - ofst) >= (mpiAload_temp - mpiCollisionWorkAlloc(3,i))) then ! if the current chunk has more records than average load, 
       ofst = ofst + mpiAload_temp - mpiCollisionWorkAlloc(3,i)                    ! set up the offset for the next processor (part of the chunk j will go there) 
       mpiCollisionWorkAlloc(3,i) = mpiAload_temp ! assign the load
       ! Need to include a check in case the chunk has exactly as many records as is needed 
       ! in this case we need to move onto the next chunk for other procs...
       if (ofst == chnks_Act(j)) then 
        j=j+1 ! move the pointer to the next chunk
        ofst = 0
       end if  
      else 
       mpiCollisionWorkAlloc(3,i) = mpiCollisionWorkAlloc(3,i) + chnks_Act(j) - ofst ! count in the remaining records from chunk j. However, more records needed for this processor 
       j = j+1                                    ! move on to the next chunk
       ofst = 0                                   ! rest the offset
      end if
     end do   
    end do
   end if 
  end do ! end loop in copies of A 
  ! the array of work allocation has been produced. Now we need to pass it on to the slave processors via a boradcast 
  ! prepare the buffer
  do i=1, mpicommworldsize-1
   do j=1,3
    mpiI4BArryBuffer((i-1)*3+j) = mpiCollisionWorkAlloc(j,i)
   end do 
  end do 
  ! send the buffer via boradcast
  call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_ReadOpA_highperf: MPI broadcast from proc 0 returned error", ierr
   stop
  end if  
  !! The workload has been passed onto the slave processors.
  deallocate(chnks_Act,mpiCollisionWorkAlloc) ! The master processor will not use these anymore... 
 !! the evaluation of collision integral will be performed in an ierarchical fashion: Each group of processors will be 
 !! united in a universe Acopy_UNIV. The first processor in the groupd will recieve rank = 0 and will be assignedas the master for
 !! this univers. The slave processors of a copy of Acopy_UNIV will send data to the master using collective communication.
 !! Also the processor with rank zero can still collect data frorm all toher processor and use broadcasts send the data out to 
 !! slave processors in all Acopy_UNIV universes 
 !! END OF THE STEP 1:
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:  
 !! Now each of the slave processors will read the A-arrays and prepare their local copies of A-arrays. 
 !! The slave processors will produce thier local A-capphi arrays. This information will be important on the next step: 
 !! the determination of each individual processor's job allocation in evaluation of the collision integral. 
  
 !! It is decided to let the slave processors divide the velocity nodes between the groups of processors
 !! essentially each will run acopy of the subroutine and should have the same answer. Then they will pick their portion of 
 !! nodes and will begin process them. In the end they will produce the workloads for themselves. Workload is the list of 
 !! velocity nodes for which the slave processor evaluates the collision operator 
 !!
 !! END OF  STEP 2
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Finally we setup some auxiliary arrays that we will need in the job assignment routin. 
 !! (We really just need the nodes_phican but we want to use the subroutine that sets all the arrays...
 !! first, we need to create the array A_capphi. it is not used on the master node, but we need a dummy array to 
 !! pass it to the subroutine so as to avoid mistakes.
 nn = size(nodes_u,1)
 allocate (A_capphi(1:nn), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variable A_capphi. stop"
   stop
  end if
 A_capphi=0 ! This array is not used on the master node
 !
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! End of preparing the array nodes_phican
 !!
 !!!!!!!!!!!!!!!!!!!!! 
else
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!
 !! SLAVE PROCESS CODE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 1: 
 !! 
 !! The slave processor receives the A-array allocation that will be broadcasted from the master node: 
 call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"PrepMPICollsnDGV_ReadOpA_highperf: slave processor", irank, "MPI boradcast from proc 0 returned error", ierr
  stop
 end if   
 !! Now we need to retreve three numbers that pertain to this particular processor. 
 !! Specifically, we need to get the fstchnk=mpiCollisionWorkAlloc(irank,1), 
 !! ofst = mpiCollisionWorkAlloc(irank,2) and numrec = mpiCollisionWorkAlloc(irank,3) 
 !! The rest of the array is not important for this process. Therefore, 
 fstchnk=mpiI4BArryBuffer((irank-1)*3+1)
 ofst=mpiI4BArryBuffer((irank-1)*3+2)
 numrec=mpiI4BArryBuffer((irank-1)*3+3)
 !! END OF STEP 1:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:
 !! We proceed processing the received information 
 !! Now each processor needs to read the corresponding part of the A array and get ready for work. 
 call ReadAarraysChnksOfstNrecDGV(fstchnk,ofst,numrec,numchnks)
 !! Now local copies of the A-arrays are set. We are ready to detemine the work load in terms of the velocity nodes 
 !! Need to run this setup subroutine: Presently we really just need the nodes_canphi but we want to use the subroutine that sets all the arrays...
 !! sets up arrays: nodes_Ashift,nodes_phican, nodes_dui,nodes_dvi,nodes_dwi
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! call the subrouine that evaluates the workload: velocity nodes that will be processed on this processor.
 call SetMPILocProcNodesWorkLoad_DGV(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)
 !! the node work allocation is set 
 !!
 !! END OF STEP 2.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End distributing workload for the slave processors for evaluation of full Boltzmann operator.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(mpiI4BArryBuffer)
end subroutine PrepMPICollsnDGV_ReadOpA_highperf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! PrepMPICollsnDGV_ReadOpA_highperf_SelNds(irank)
!
! This is a copy of the above subroutine except it is designed to use only a part of the primary velocity nodes. 
! The array of nodes used in simulation as a variable sel_nodes. The nodes listed in the array sel_nodes will be divided betwen processors 
! the numbers in the arrays nodes_proc_groups,procs_lin correspond to index in array sel_nodes and entries of sel_nodes are the numbers of the 
! node in nodes_u/_v/_w etc. arrays.
! Note that we do not need to see the arrya of nodes, we lust need to know how many are there, so we do not pass the array sel_nodes itself, just the number of records in it.
!
!
!
! This subroutines is a chunk of PrepMPICollsnDGV_highperf. The splitting alows for 
! more flexibility in customizing MPI codes developed using the DGv library
!
! This piece distributes the A operator between the processors that evaluate the collision operator 
!
!!!!!!!!
! ATTENTION! Must be called after PrepMPICollsnDGV_Univs_highperf or will not work correctly
!!!!!!!!
! This subroutines contains an MPI fork 
! This subroutine is run on all processors
!
! Thi subroutine prepares MPI environments and data arrays that will be 
! used for MPI parallelization for 
!    (1) evaluation of the Boltzmann collision operator
!    (2) can be used for parallel evaluation of the linearized Boltzmann collision operator
!    (3) can be used with a spatial driver to parallilize in velocity space. Specifically, nodes of the 
!         velocity mesh are divided between MPI communicating processors. Each processor is running 
!         a copy of spatial operator for a portion of velocity nodes, the processor request values of the collision 
!         operator for their portion of nodes. 
! 
! Description: 
!
! This subroutine prepares a number of communicators. 
!
!  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
!
!  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
!                                  sharing a copy of A
! 
!  MPI_Acopy_UNIVS(1:num_Acopies), MPI_Acopy_UNIVS_GROUP(1:num_Acopies) - pointer to universes are stored in these arrays.
!                                  each pointer is pointing to a group of proecessors sharing a copy of A
!  
!  LinFmA_UNIV_GROUP, LinFmA_UNIV, -- these are scrap variable used to create a universe for a group of processors 
!                                  that receives information from a single group of processors sharing a copy of A
!                                  Here is some idea on how these communicators work: 
!										There is exactly num_Acopies of groups of processors. Let's sasy its number is 0 <= I <=num_Acopies.
!										Each of the goups has a workload, given by (nodes_proc_groups(1,I), \ldots, nodes_proc_groups(2,I)), 
!										that is a number of velocity points 
!										at which this group of processors I evaluates the collision operator.
!                                      Just in case, this group of processors is included in the communicator MPI_Acopy_UNIVS(I), but that is 
!										not important at the moment. Values of the collision operator that are computed in the group I   
!										may be neded on other processors for either evaluation of 
!										linearized collision operator or for parallel evaluation of the kinetic equations using 
!										parallelization in the velocity nodes. This is determined by examining array procs_lim, that 
!										a range of velocity points to each participating processor. If a processor needs a value of the collision 
!										operator computed in the group I, then the processor 
!                                      is included in the LinFmA_Univ/_Group. Also, the master node of the MPI_Acopy_UNIVS(I) is included 
!										included in the LinFmA_Univ(I) (and also made a master node in it?) 
!										Typically, the master node of MPI_Acopy_UNIVS(I) broadcasts values to all of LinFmA_Univ(I)
!
!  MPI_LinFmA_UNIVS(1:num_Acopies), MPI_LinFmA_UNIVS_GROUP -- pointer to universes are stored in these arrays.
!                                  each pointer is pointing to a group of proecessors that receives information 
!                                  from a single group of processors sharing a copy of A
!
!
!
!!! MPI FORK
!!! Now the master process will distribute the chunks of A arrays to other processors. If the  
!!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
!!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
!!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
!!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
!!! Therefore, the results will be combined in an integer of (I8B) kind. 
!!!
!!! This number will be divided by the number of 
!!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
!!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
!!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
!!!  
!!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
!!! This boradcast will go out with the code 30.
!!!
!!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
!!! 
!!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
!!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
!!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
!!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
!!! 
 
!!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
!!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
!!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
!!!? workload. 
!!! 
!!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
!!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
!!!
!!! Next the master processor initiates the work: It simply broadcasts the solution, 
!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! STEP 3.
 !!
 !! PARALLELIZATION OF THE LINEARIZATION STAFF 
 !!
 !!
 !! This part of the code replaces the original verison that is using 
 !! send_lin_procs and recv_lin_procs arrays for the construciton of the linearized collision operator
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! ATTENTION: OUTDATED!!!!
 !! The new algorithm for evaluatin the collision operator is as follows:
 !!
 !! A. Slave Processors  compute the components of the linearized operator.
 !! B. Slave processors place these componets in the copy of an array 
 !!    FmA whose rows go over all nodes in the array Acopy_Wrkld that are the same for 
 !!    all processors in the same Acopy universe
 !!    to help them stack the components the array procs_nodes_wkld_Acopy_Addr array us used. 
 !!    procs_nodes_wld_Acopy_addr has the same amount of elements as procs_nodes_wld
 !!    The array procs_nodes_wld_Acopy_addr contains indices of the nodes from 
 !!    procs_nodes_wld in Acopy_Wrkld
 !!    The lead processor of Acopy univ. gathers the FmA to it
 !! 
 !! C. The lead proc of Acopy universe is also the lead processor of the matching LinFmA universe 
 !!    This processor broadcasts FmA to processors involved in the valuation of the collision oeprator. 
 !!    For each Acopy universe there is a LinFmA_Univ universe that combines all processors that need 
 !!    to receive data from processors  invloved in Acopy_univ. 
 !! D. The processors pick the data from the boradsact and build the local copy of the linearized operator.
 !!    To do that three arrays are used:
 !!    lin_proc_nodes_wld --- contains the componets of the linearized operator assigned to each processor involved 
 !!    in the evaluation of the linearized operator. 
 !!    linprocndswld_Acopy_univs, --- this array contains informaion about arriving components of the linearized operator:
 !!     # of the universe A copy that will send the components
 !!     # total number of the nodes (each node is one row of the lineaired collision operator
 !!     # Nodes .... 
 !!    then repeat
 !!    linprocndswld_Acopy_addr has the same stucture except instead of nodes, there will be their indices in the corresponding 
 !!    Acopy_Wrkld array
 !!
 !! E. The processors involced in the evaluation of the linearized colliion operator create a copy of BFmA where they store the 
 !!    componets of the consolidated operator.  
 !!    the processors evaluate the parts of the consolidated collision operator and stack the results in the temporary solution 
 !!    array. Finally master processrs gathers all componets of the temp. solution array by a collective communication call. 
 !!!! END OUTDATED  !!
 !!!!!!!!!!!!!!!!!!!!!


!!!!!!!!
! ATTENTION! The subroutine assumes that PrepMPICollsnDGV_Univs_highperf has been called already and important 
!            variable have been set up. Will not work correctly is called before PrepMPICollsnDGV_Univs_highperf
!!!!!!!!

subroutine PrepMPICollsnDGV_ReadOpA_highperf_SelNds(irank,chnks_copies,nodes_proc_groups,sel_nodes)

use DGV_commvar, only: numchnks,num_Acopies,nodes_u,A_capphi

use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:) :: chnks_copies !
integer (I4B), dimension (:,:) :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:), intent (in) :: sel_nodes ! this is the array that contains numbers of the selected velocity nodes on the entire grid. Arrays 
                    ! chnks_copies and nodes_proc_groups contain indices of entries in sel_nodes

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 

integer (I4B), dimension (:), allocatable :: chnks_Act ! This is essentially a scrap array that keeps the information about number of records or A operator in each chunk of Aarrays.
integer (I4B) :: i,jj,j,ofst,leftover,nn,fstchnk,numrec   ! Scrap index
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I8B) :: Act_aggreg ! variable to keep the total number of records in all chunks of operator A

!!!!! MPI variablessortNodesProcsBuffer
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer     ! this buffer will be used in sending and receiving the work load.
integer (I4B), dimension (:,:), allocatable :: mpiCollisionWorkAlloc ! This array will contain info about work allocations for slave processors  
                                                                ! The array is operating like the following: the first index is the number of processor
                                                                ! the second index has 3 fixed fields 
                                                                ! record 1: number of first chunk to read
                                                                ! record 2: offset in the chunks to read
                                                                ! record 3: number of records to read
integer (I4B) ::  mpiAload, mpiAloadrem, mpiAload_temp  ! Scrap variables to store the averaghe number of records in A array that will be stored  in the 
                                             ! each slave processor and the reminder of the division of the total number of records in A by the number of 
                                             ! processors     
                                             
! variables for MPI Calls
integer :: ierr, mpicommworldsize

!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! we now allocate the buffer for the mpiCollisionWorkAlloc array
allocate (mpiI4BArryBuffer(1:3*(mpicommworldsize-1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables mpiI4BArryBuffer"
 stop
end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Next we need to distribute the workload to the slave processors for 
 ! evaluating the full collision operator. The starting poiint is to divide pre=-computed operator A between the processors.
 !!!!!!!!!!!!!!!!!
 !!  MPI fork
 if (irank==0) then
  !! STEP 1:
  !! Master process code: 
  !! The Master node needs to disctribute operator A. Therefore it will scall all chunks of all arrays and will 
  !! retrieve the total number of the recods in operator A. Then this number will be divided to the number of processors 
  !! to produce the operator A workloads (+/- 1)
  !! begin STEP 1:
  allocate (chnks_Act(1:numchnks), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables chnks_Act"
   stop
  end if
  !! The next subroutine will scan all chunks of Aarrays and will compile the information about the number of records in 
  !! each chunk in array chnks_Act
  ! call next subroutine to populate chnks_Act
  call ScanAarrsChnks4Acapphi(chnks_Act) ! Act_aggreg will have the total number of records in the Aarrays.
  ! Let us compute the total number of records
  Act_aggreg=0
  do i=1,numchnks
   Act_aggreg=Act_aggreg+chnks_Act(i)
  end do
  !! now we broadcast the array to all processors. 
  !! Now it is time to set up the work allocation array.
  allocate (mpiCollisionWorkAlloc(1:3,1:mpicommworldsize-1), stat=loc_alloc_stat)
     ! 
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variables mpiCollisionWorkAlloc"
   stop
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mpiCollisionWorkAlloc = 0 ! reset the array
  do jj=1,num_Acopies
   !! Next we divide the total number of records by the total number of available processors to determine the workload per processor
   !! We perform division with reminder
   mpiAload = INT(Act_aggreg/(chnks_copies(2,jj)-chnks_copies(1,jj)+1),I4B)
   mpiAloadrem = INT((Act_aggreg - (INT8(chnks_copies(2,jj)-chnks_copies(1,jj)+1)*mpiAload)), I4B) ! notice that mpiAloarrem < mpicommworldsize-1
   !! A quick check that we have enough memory per processor. If we do not, stop:
   if (mpiAload > 2*10**8) then 
    print *,"PrepMPICollsnDGV_ReadOpA_highperf: Error> the number of records of operator A exceeds 2*10^8. Aborting the solver."
    stop
   end if 
   !! Next we populate the work allocation array.
   !! There will be two cases. Case 1 -- number of slave processors= number of chunks 
   !! and Case 2: number of processors neq number of chunks ... 
   if (numchnks == (chnks_copies(2,jj)-chnks_copies(1,jj)+1)) then 
   !! in the case 1 the chunks are assigned to the processor withour rechunking. 
   !! SPECIAL CASE :: NUMBER OF CHUNKS AND NUMBER OF PROCESSORS IS THE SAME! !! IMPORTANT!!
    do i = chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = i - chnks_copies(1,jj) + 1
     mpiCollisionWorkAlloc(2,i) = 0
     mpiCollisionWorkAlloc(3,i) = chnks_Act(i - chnks_copies(1,jj) + 1)
    end do
   else 
   !! in the case 2, the following subroutine determines re-chunking the A-arrays. 
    j=1 ! this is the index that will run over chunks
    ofst=0 ! This is the offset calculator
    leftover=0
    do i=chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = j  
     mpiCollisionWorkAlloc(2,i) = ofst
     if (i <= mpiAloadrem) then
      mpiAload_temp = mpiAload+1
     else 
      mpiAload_temp = mpiAload
     end if   
     do while (mpiCollisionWorkAlloc(3,i) < mpiAload_temp) ! process chunks until i-th processor is loaded. 
      ! a quick check if data is corrupted. if our calculations above arte correct, then there will be exactly as many 
      ! chunks as needed. If we ran oput of chunks -- somewhere there is an error. Then we stop
      if (j > numchnks) then 
       print *, "PrepMPICollsnDGV_ReadOpA_highperf: error when trying to populate mpiCollisionWorkAlloc. Ran out of chunks."
       stop
      end if
      ! end error check  
      if ((chnks_Act(j) - ofst) >= (mpiAload_temp - mpiCollisionWorkAlloc(3,i))) then ! if the current chunk has more records than average load, 
       ofst = ofst + mpiAload_temp - mpiCollisionWorkAlloc(3,i)                    ! set up the offset for the next processor (part of the chunk j will go there) 
       mpiCollisionWorkAlloc(3,i) = mpiAload_temp ! assign the load
       ! Need to include a check in case the chunk has exactly as many records as is needed 
       ! in this case we need to move onto the next chunk for other procs...
       if (ofst == chnks_Act(j)) then 
        j=j+1 ! move the pointer to the next chunk
        ofst = 0
       end if  
      else 
       mpiCollisionWorkAlloc(3,i) = mpiCollisionWorkAlloc(3,i) + chnks_Act(j) - ofst ! count in the remaining records from chunk j. However, more records needed for this processor 
       j = j+1                                    ! move on to the next chunk
       ofst = 0                                   ! rest the offset
      end if
     end do   
    end do
   end if 
  end do ! end loop in copies of A 
  ! the array of work allocation has been produced. Now we need to pass it on to the slave processors via a boradcast 
  ! prepare the buffer
  do i=1, mpicommworldsize-1
   do j=1,3
    mpiI4BArryBuffer((i-1)*3+j) = mpiCollisionWorkAlloc(j,i)
   end do 
  end do 
  ! send the buffer via boradcast
  call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_ReadOpA_highperf: MPI broadcast from proc 0 returned error", ierr
   stop
  end if  
  !! The workload has been passed onto the slave processors.
  deallocate(chnks_Act,mpiCollisionWorkAlloc) ! The master processor will not use these anymore... 
 !! the evaluation of collision integral will be performed in an ierarchical fashion: Each group of processors will be 
 !! united in a universe Acopy_UNIV. The first processor in the groupd will recieve rank = 0 and will be assignedas the master for
 !! this univers. The slave processors of a copy of Acopy_UNIV will send data to the master using collective communication.
 !! Also the processor with rank zero can still collect data frorm all toher processor and use broadcasts send the data out to 
 !! slave processors in all Acopy_UNIV universes 
 !! END OF THE STEP 1:
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:  
 !! Now each of the slave processors will read the A-arrays and prepare their local copies of A-arrays. 
 !! The slave processors will produce thier local A-capphi arrays. This information will be important on the next step: 
 !! the determination of each individual processor's job allocation in evaluation of the collision integral. 
  
 !! It is decided to let the slave processors divide the velocity nodes between the groups of processors
 !! essentially each will run acopy of the subroutine and should have the same answer. Then they will pick their portion of 
 !! nodes and will begin process them. In the end they will produce the workloads for themselves. Workload is the list of 
 !! velocity nodes for which the slave processor evaluates the collision operator 
 !!
 !! END OF  STEP 2
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Finally we setup some auxiliary arrays that we will need in the job assignment routin. 
 !! (We really just need the nodes_phican but we want to use the subroutine that sets all the arrays...
 !! first, we need to create the array A_capphi. it is not used on the master node, but we need a dummy array to 
 !! pass it to the subroutine so as to avoid mistakes.
 nn = size(nodes_u,1)
 allocate (A_capphi(1:nn), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperf: Allocation error for variable A_capphi. stop"
   stop
  end if
 A_capphi=0 ! This array is not used on the master node
 !
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) ! Need to fix this! replace FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP) with cell number
 !! End of preparing the array nodes_phican
 !!
 !!!!!!!!!!!!!!!!!!!!! 
else
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!
 !! SLAVE PROCESS CODE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 1: 
 !! 
 !! The slave processor receives the A-array allocation that will be broadcasted from the master node: 
 call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"PrepMPICollsnDGV_ReadOpA_highperf: slave processor", irank, "MPI boradcast from proc 0 returned error", ierr
  stop
 end if   
 !! Now we need to retreve three numbers that pertain to this particular processor. 
 !! Specifically, we need to get the fstchnk=mpiCollisionWorkAlloc(irank,1), 
 !! ofst = mpiCollisionWorkAlloc(irank,2) and numrec = mpiCollisionWorkAlloc(irank,3) 
 !! The rest of the array is not important for this process. Therefore, 
 fstchnk=mpiI4BArryBuffer((irank-1)*3+1)
 ofst=mpiI4BArryBuffer((irank-1)*3+2)
 numrec=mpiI4BArryBuffer((irank-1)*3+3)
 !! END OF STEP 1:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:
 !! We proceed processing the received information 
 !! Now each processor needs to read the corresponding part of the A array and get ready for work. 
 call ReadAarraysChnksOfstNrecDGV(fstchnk,ofst,numrec,numchnks)
 !! Now local copies of the A-arrays are set. We are ready to detemine the work load in terms of the velocity nodes 
 !! Need to run this setup subroutine: Presently we really just need the nodes_canphi but we want to use the subroutine that sets all the arrays...
 !! sets up arrays: nodes_Ashift,nodes_phican, nodes_dui,nodes_dvi,nodes_dwi
 call SetCellsUGINdsSftArrs(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! call the subrouine that evaluates the workload: velocity nodes that will be processed on this processor.
 call SetMPILocProcNodesWorkLoad_DGV_SelNds(chnks_copies,nodes_proc_groups,mpicommworldsize,irank,sel_nodes)
 !! the node work allocation is set 
 !!
 !! END OF STEP 2.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End distributing workload for the slave processors for evaluation of full Boltzmann operator.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(mpiI4BArryBuffer)
end subroutine PrepMPICollsnDGV_ReadOpA_highperf_SelNds

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! PrepMPICollsnDGV_Univs_highperfII
 !
 ! This is a modification of the subroutine PrepMPICollsnDGV_Univs_highperf to 
 ! enable MPI parallel evaluation of the collision integral on SECONDARY MESH 
 ! 
 ! NOTE THAT LINEARIZATION IS NOT EXPECTED TO BE USED ON SECONDARY MESH. Therefore, 
 ! the corresponding parts of the subroutine PrepMPICollsnDGV_Univs_highperf were not extended. 
 !
 ! The next subroutines are chunks of PrepMPICollsnDGV_highperf. The splitting alows for 
 ! more flexibility in customizing MPI codes developed using the DGv library
 !
 ! This piece creates communication universes for primary velocity mesh. 
 ! 
 !
 ! This subroutines contains an MPI fork 
 ! This subroutine is run on all processors
 !
 ! Thi subroutine prepares MPI environments and data arrays that will be 
 ! used for MPI parallelization for 
 !    (1) evaluation of the Boltzmann collision operator
 ! 
 ! Description: 
 !
 ! This subroutine prepares a number of communicators. 
 !
 !  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
 !                                  sharing a copy of A
 ! 
 !  MPI_Acopy_UNIVSII(1:num_Acopies), MPI_Acopy_UNIVS_GROUPII(1:num_Acopies) - pointer to universes are stored in these arrays.
 !                                  each pointer is pointing to a group of proecessors sharing a copy of A
 !  
 !
 !
 !!! MPI FORK
 !!! Now the master process will distribute the chunks of A arrays to other processors. If the  
 !!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
 !!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
 !!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
 !!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
 !!! Therefore, the results will be combined in an integer of (I8B) kind. 
 !!!
 !!! This number will be divided by the number of 
 !!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
 !!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
 !!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
 !!!  
 !!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
 !!! This boradcast will go out with the code 30.
 !!!
 !!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
 !!! 
 !!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
 !!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
 !!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
 !!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
 !!! 
 
 !!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
 !!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
 !!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
 !!!? workload. 
 !!! 
 !!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
 !!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
 !!!
 !!! Next the master processor initiates the work: It simply broadcasts the solution, 
 !!! 
 !!!
 !!! The slave processes will receive the solutions and will start to evaluate their portions of the collision operator. 
 !!! Then they will send their results to the master processor using an aggregate reduce operation.
 !!!
 !!! Finally, let us briefly describe the process of calculating the collisoin opperator by slave 
 !!! processors. once the I1workload is set,  the slave processor will wait for a boradcast of solution. 
 !!! Each slave process then goes over the I1workload and performs the 
 !!! evaluation of the collision integral for each I1 listed there. The results are stored in the array (CollisionInt)
 !!! Once all I1 are processed, this array is aggregately sent to the master node. Operation of the summation is used to combine 
 !!! results for all chunks. 
 !!! 
 !!! Also, the slave process will evaluate the linearized part
 !!! also, sometimes the slave processor will send the linearized part to some other processor for faster assembly...
 !!!
 !!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  PrepMPICollsnDGV_Univs_highperfII(irank,chnks_copies,nodes_proc_groups)

use DGV_commvar, only: num_lin_proc, nodes_uII, num_Acopies, &
                       MPI_Acopy_UNIVSII, MPI_Acopy_UNIVS_GROUPII, &
                       My_Acopy_Univ_indxII, Acopy_irankII
use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:), intent (in) :: chnks_copies !
integer (I4B), dimension (:,:), intent (in) :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nn,i,j,jj,ii,size_est,num_ranks,node ! scrap index
                 
! variables for MPI Calls
integer ::  mpicommworldsize,ierr,iprc,sendproc
integer :: MPI_COMM_WORLD_GROUP ! to keep the handle to the group of MPI_COMM_WORLD
integer, dimension(3) :: ranges ! ranges to create a communicator group
integer :: Acopy_UNIV_GROUP,Acopy_UNIV ! scrap pointers to Acopy and LinFmA universes
!!  Acopy_UNIV, Acopy_UNIV_GROUP ! on the master node this is a temporary holder for the pointer to a universe, and communicatio ngroup that contains processors dealing with one copy of A operator
!!  on the slave node this is avariable for the Acopy communicator
integer, dimension (1) :: ranksII_scrap ! This is a scrap variable to check ranks of the processors in the MPI_LINEAR_COLL_WORLD communicator

!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_Univs_highperf: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Here we set up some variables that will be used by all processes: 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_WORLD_GROUP, ierr) ! get the hanlde to the group of MPI_COMM_WORLD.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! We now create universes to deal with the groups of processors sharing a single copy 
!!! of A and the universes of processors that obtain 
!!! components of the colliion operator from the same group of nodes... 
!!!
!!! These universies will be used for communcation of the components of the linearized operator.
!!! These universes can also be used for evaluation of the collision operator if the spatial operator is parallelized  
!!! 
!!! Begin creating the Acopy_UNIV universes:
!!! first, we allocate storage for the universes... 
 allocate (MPI_Acopy_UNIVSII(1:num_Acopies),MPI_Acopy_UNIVS_GROUPII(1:num_Acopies), stat=loc_alloc_stat)
    ! 
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_Univs_highperfII: Allocation error for MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP. irank", irank 
  stop
 end if
 MPI_Acopy_UNIVSII=0;MPI_Acopy_UNIVS_GROUPII=0; ! Nullify them , just in case...
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Reset the index variable 
 My_Acopy_Univ_indxII = 0 ! points to no particular Acopy Universe
 ! now we will create num_Acopies universes that will contain the master node and the nodes indentified in chnks_copies(1:2,j)
 do i=1,num_Acopies
  ranges(1)=chnks_copies(1,i); ranges(2)=chnks_copies(2,i); ranges(3)=1; 
  call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges,Acopy_UNIV_GROUP, ierr) 
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperfII: MPI can not create group Acopy_UNIVS_GROUP on proc", irank, ". Error=", ierr
   stop
  end if
  ! a quick check if the processor with irank=0 int he new group correspond to the first processor in the group.
  call MPI_GROUP_TRANSLATE_RANKS(Acopy_UNIV_GROUP, 1, (/ 0 /), MPI_COMM_WORLD_GROUP, ranksII_scrap, ierr)  
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperfII: can't translate Acopy_UNIV_GROUP=>MPI_COMM_WORLD_GROUP. i=",& 
             i," proc",irank,"Err=",ierr
   stop
  end if
  ! a quick check that the zeros process is still number zero in the new group: 
  if (ranksII_scrap(1) /= chnks_copies(1,i)) then 
   print *,"PrepMPICollsnDGV_Univs_highperfII: ranks mixed up between Acopy_UNIV_GROUP and MPI_COMM_GROUP i=",& 
              i,"proc.",irank,"stop."
   stop
  end if 
  ! end check 
  ! next we create a new communicator: 
  call MPI_COMM_CREATE (MPI_COMM_WORLD, Acopy_UNIV_GROUP, Acopy_UNIV, ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_Univs_highperfII: can't create communicator Acopy_UNIV.i=",i,"Proc=",irank,". Err=",ierr
   stop
  end if
  !! Record the pointers to the Communicator and the Group:
  MPI_Acopy_UNIVSII(i)=Acopy_UNIV
  MPI_Acopy_UNIVS_GROUPII(i)=Acopy_UNIV_GROUP
  !!!!
  ! We also need to remember to what Acopy universe we belong
  if ((irank >= chnks_copies(1,i)) .and. (irank<= chnks_copies(2,i))) then 
   My_Acopy_Univ_indxII = i ! now we know which Acopy universe contains this slave processor. Master processor will still have zero in this variables
  end if 
 end do  
 !! END creating the Acopy_UNIV universes:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will use the costant My_Acopy_Univ_indx to set 
 ! the variable Acopy_irank which defines the rank of the processor in this A copy Universe
 !!!!!!!!!!!!!!!! 
 Acopy_irankII = -1
 if (My_Acopy_Univ_indxII>0) then 
 !!! MPI FORK
   Acopy_UNIV=MPI_Acopy_UNIVSII(My_Acopy_Univ_indxII) 
   call mpi_comm_rank(Acopy_UNIV,Acopy_irankII,ierr) ! check what processor the program is on... 
   if (ierr /= 0 ) then 
    print *,"PrepMPICollsnDGV_Univs_highperfII: can't determined rank in Acopy_UNIV.i=", &
             My_Acopy_Univ_indxII,"Proc=",irank,". Err=",ierr
   stop
  end if
 end if   
 !!!!!!!!!!!!!!!
 ! end setting Acopy_irank
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end subroutine  PrepMPICollsnDGV_Univs_highperfII

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! PrepMPICollsnDGV_ReadOpA_highperfII(irank)
!
! This is a modification of the subroutine PrepMPICollsnDGV_ReadOpA_highperf to 
! enable MPI parallel evaluation of the collision integral on SECONDARY MESH 
!
! NOTE THAT LINEARIZATION IS NOT EXPECTED TO BE USED ON SECONDARY MESH. Therefore, 
! the corresponding parts of the subroutine PrepMPICollsnDGV_Univs_highperf were not extended. 
!
!
! This subroutines is a chunk of PrepMPICollsnDGV_highperf. The splitting alows for 
! more flexibility in customizing MPI codes developed using the DGv library
!
! This piece distributes the A operator between the processors that evaluate the collision operator 
!
!!!!!!!!
! ATTENTION! Must be called after PrepMPICollsnDGV_Univs_highperf or will not work correctly
!!!!!!!!
! This subroutines contains an MPI fork 
! This subroutine is run on all processors
!
! Thi subroutine prepares MPI environments and data arrays that will be 
! used for MPI parallelization for 
!    (1) evaluation of the Boltzmann collision operator
! 
! Description: 
!
! This subroutine prepares a number of communicators. 
!
!  MPI_LINEAR_COLL_WORLD is the  dedicated to communictor for processors with numbers 0, ..., num_lin_proc-1
!
!  Acopy_UNIV, Acopy_UNIV_GROUP -- these are scrap variable used to create a universe for each group of processors 
!                                  sharing a copy of A
! 
!  MPI_Acopy_UNIVSII(1:num_Acopies), MPI_Acopy_UNIVS_GROUPII(1:num_Acopies) - pointer to universes are stored in these arrays.
!                                  each pointer is pointing to a group of proecessors sharing a copy of A
!  
!
!
!!! MPI FORK
!!! Now the master process will distribute the chunks of A arrays to other processors. If the  
!!! if number of chunks is not equal to the number of processors, then A is re-chunked. Specifically, 
!!! The array chnks_Act will be assembles from all chunks. This array will contain the numbers of records of A array 
!!! that is contained in each idnividua chunk. Then the total number records in A will be calculated =A_ct_aggregate. 
!!! Keep in mind that while individual entries of A_capphi will arrive in integer (I4B) kind. Their sum may will exceed this kind. 
!!! Therefore, the results will be combined in an integer of (I8B) kind. 
!!!
!!! This number will be divided by the number of 
!!! processors. Suppose we had P slave processes and the result of this deivision was A_ct=A_ct_per_proc*P+r 
!!! Then the first r processes (r<P) will recieve A_ct_per_proc+1 records and the last P-r will recive A_ct_per_proc records. 
!!! The master process will prepare an array mpiCollisionWorkAlloc with the following information: 
!!!  
!!! Number of processor, number of the chunk to start, offset in the first chunk, total number of records. 
!!! This boradcast will go out with the code 30.
!!!
!!! Notice that with this way of distributing records no two processors that share a copy of A will have use same record of A.   
!!! 
!!! The slave process will read the chunk(s) and sets the displacements arrays based on the information about the grid
!!! NOTE: There is one array, called nodes_Ashift that uses continuous numbering in the array A. This array 
!!! is built based on A and on the A_capphi. This will be different for all slave processes. However, arrays
!!! nodes_phican, nodes_dui, nodes_dvi, nodes_dwi will be the same for all nodes (because they only depend on the grid.
!!! 
 
!!!? Next the slave processor will send out its local copy of A_capphi. This copy will be gathered on the master node
!!!? and used to assign the nodes work load to all the processors. The nodes workload will be stored in the 
!!!? array procs_nodes_wld. Then the pieces of these arrays are sent to individual processors so that they can have their 
!!!? workload. 
!!! 
!!! Next the slave processor recieves the workload: for each I1 in its chunk of A, the processor will have a bunch of 
!!! nodes with the same local address (ui,vi,wi).  This information will be kept in I1workload array. 
!!!
!!! Next the master processor initiates the work: It simply broadcasts the solution, 
!!!


!!!!!!!!
! ATTENTION! The subroutine assumes that PrepMPICollsnDGV_Univs_highperf has been called already and important 
!            variable have been set up. Will not work correctly is called before PrepMPICollsnDGV_Univs_highperf
!!!!!!!!

subroutine PrepMPICollsnDGV_ReadOpA_highperfII(irank,chnks_copies,nodes_proc_groups)

use DGV_commvar, only: numchnksII,num_Acopies,nodes_uII,A_capphiII

use DGV_dgvtools_mod 
use DGV_readwrite 
use DGV_miscset
                  
!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    

integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.
integer (I4B), dimension (:,:), intent (in) :: chnks_copies !
integer (I4B), dimension (:,:), intent (in) :: nodes_proc_groups ! this array will contain the breaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           

!!!!!!!!!!!!!!!!!!!!!!!!!! LOCAL VALRIABLES 

integer (I4B), dimension (:), allocatable :: chnks_Act ! This is essentially a scrap array that keeps the information about number of records or A operator in each chunk of Aarrays.
integer (I4B) :: i,jj,j,ofst,leftover,nn,fstchnk,numrec   ! Scrap index
integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I8B) :: Act_aggreg ! variable to keep the total number of records in all chunks of operator A

!!!!! MPI variablessortNodesProcsBuffer
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer     ! this buffer will be used in sending and receiving the work load.
integer (I4B), dimension (:,:), allocatable :: mpiCollisionWorkAlloc ! This array will contain info about work allocations for slave processors  
                                                                ! The array is operating like the following: the first index is the number of processor
                                                                ! the second index has 3 fixed fields 
                                                                ! record 1: number of first chunk to read
                                                                ! record 2: offset in the chunks to read
                                                                ! record 3: number of records to read
integer (I4B) ::  mpiAload, mpiAloadrem, mpiAload_temp  ! Scrap variables to store the averaghe number of records in A array that will be stored  in the 
                                             ! each slave processor and the reminder of the division of the total number of records in A by the number of 
                                             ! processors     
                                             
! variables for MPI Calls
integer :: ierr, mpicommworldsize



 
!!!!!!!!!!!!!!!!!!!!!!!!!
! A quick array sizes consistency check. For program to operate correctly, this check should pass:
if (.NOT. ((size(chnks_copies,1)==2) .AND. (size(chnks_copies,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperfII: incorrect size of passed array chnks_copies. Stop"
 stop 
end if 
if (.NOT. ((size(nodes_proc_groups,1)==2) .AND. (size(nodes_proc_groups,2)==num_Acopies))) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperfII: incorrect size of passed array nodes_proc_groups. Stop"
 stop 
end if 
! End of array sizes consistency check 
!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! we now allocate the buffer for the mpiCollisionWorkAlloc array
allocate (mpiI4BArryBuffer(1:3*(mpicommworldsize-1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "PrepMPICollsnDGV_ReadOpA_highperfII: Allocation error for variables mpiI4BArryBuffer"
 stop
end if

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Next we need to distribute the workload to the slave processors for 
 ! evaluating the full collision operator. The starting poiint is to divide pre=-computed operator A between the processors.
 !!!!!!!!!!!!!!!!!
 !!  MPI fork
 if (irank==0) then
  !! STEP 1:
  !! Master process code: 
  !! The Master node needs to disctribute operator A. Therefore it will scall all chunks of all arrays and will 
  !! retrieve the total number of the recods in operator A. Then this number will be divided to the number of processors 
  !! to produce the operator A workloads (+/- 1)
  !! begin STEP 1:
  allocate (chnks_Act(1:numchnksII), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperfII: Allocation error for variables chnks_Act"
   stop
  end if
  !! The next subroutine will scan all chunks of Aarrays and will compile the information about the number of records in 
  !! each chunk in array chnks_Act
  ! call next subroutine to populate chnks_Act
  call ScanAarrsChnks4AcapphiII(chnks_Act) ! Act_aggreg will have the total number of records in the Aarrays.
  ! Let us compute the total number of records
  Act_aggreg=0
  do i=1,numchnksII
   Act_aggreg=Act_aggreg+chnks_Act(i)
  end do
  !! now we broadcast the array to all processors. 
  !! Now it is time to set up the work allocation array.
  allocate (mpiCollisionWorkAlloc(1:3,1:mpicommworldsize-1), stat=loc_alloc_stat)
     ! 
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperfII: Allocation error for variables mpiCollisionWorkAlloc"
   stop
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mpiCollisionWorkAlloc = 0 ! reset the array
  do jj=1,num_Acopies
   !! Next we divide the total number of records by the total number of available processors to determine the workload per processor
   !! We perform division with reminder
   mpiAload = INT(Act_aggreg/(chnks_copies(2,jj)-chnks_copies(1,jj)+1),I4B)
   mpiAloadrem = INT((Act_aggreg - (INT8(chnks_copies(2,jj)-chnks_copies(1,jj)+1)*mpiAload)), I4B) ! notice that mpiAloarrem < mpicommworldsize-1
   !! A quick check that we have enough memory per processor. If we do not, stop:
   if (mpiAload > 2*10**8) then 
    print *,"PrepMPICollsnDGV_ReadOpA_highperfII: Error> the number of records of operator A exceeds 2*10^8. Aborting the solver."
    stop
   end if 
   !! Next we populate the work allocation array.
   !! There will be two cases. Case 1 -- number of slave processors= number of chunks 
   !! and Case 2: number of processors neq number of chunks ... 
   if (numchnksII == (chnks_copies(2,jj)-chnks_copies(1,jj)+1)) then 
   !! in the case 1 the chunks are assigned to the processor withour rechunking. 
   !! SPECIAL CASE :: NUMBER OF CHUNKS AND NUMBER OF PROCESSORS IS THE SAME! !! IMPORTANT!!
    do i = chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = i - chnks_copies(1,jj) + 1
     mpiCollisionWorkAlloc(2,i) = 0
     mpiCollisionWorkAlloc(3,i) = chnks_Act(i - chnks_copies(1,jj) + 1)
    end do
   else 
   !! in the case 2, the following subroutine determines re-chunking the A-arrays. 
    j=1 ! this is the index that will run over chunks
    ofst=0 ! This is the offset calculator
    leftover=0
    do i=chnks_copies(1,jj),chnks_copies(2,jj)
     mpiCollisionWorkAlloc(1,i) = j  
     mpiCollisionWorkAlloc(2,i) = ofst
     if (i <= mpiAloadrem) then
      mpiAload_temp = mpiAload+1
     else 
      mpiAload_temp = mpiAload
     end if   
     do while (mpiCollisionWorkAlloc(3,i) < mpiAload_temp) ! process chunks until i-th processor is loaded. 
      ! a quick check if data is corrupted. if our calculations above arte correct, then there will be exactly as many 
      ! chunks as needed. If we ran oput of chunks -- somewhere there is an error. Then we stop
      if (j > numchnksII) then 
       print *, "PrepMPICollsnDGV_ReadOpA_highperfII: error when trying to populate mpiCollisionWorkAlloc. Ran out of chunks."
       stop
      end if
      ! end error check  
      if ((chnks_Act(j) - ofst) >= (mpiAload_temp - mpiCollisionWorkAlloc(3,i))) then ! if the current chunk has more records than average load, 
       ofst = ofst + mpiAload_temp - mpiCollisionWorkAlloc(3,i)                    ! set up the offset for the next processor (part of the chunk j will go there) 
       mpiCollisionWorkAlloc(3,i) = mpiAload_temp ! assign the load
       ! Need to include a check in case the chunk has exactly as many records as is needed 
       ! in this case we need to move onto the next chunk for other procs...
       if (ofst == chnks_Act(j)) then 
        j=j+1 ! move the pointer to the next chunk
        ofst = 0
       end if  
      else 
       mpiCollisionWorkAlloc(3,i) = mpiCollisionWorkAlloc(3,i) + chnks_Act(j) - ofst ! count in the remaining records from chunk j. However, more records needed for this processor 
       j = j+1                                    ! move on to the next chunk
       ofst = 0                                   ! rest the offset
      end if
     end do   
    end do
   end if 
  end do ! end loop in copies of A 
  ! the array of work allocation has been produced. Now we need to pass it on to the slave processors via a boradcast 
  ! prepare the buffer
  do i=1, mpicommworldsize-1
   do j=1,3
    mpiI4BArryBuffer((i-1)*3+j) = mpiCollisionWorkAlloc(j,i)
   end do 
  end do 
  ! send the buffer via boradcast
  call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"PrepMPICollsnDGV_ReadOpA_highperfII: MPI broadcast from proc 0 returned error", ierr
   stop
  end if  
  !! The workload has been passed onto the slave processors.
  deallocate(chnks_Act,mpiCollisionWorkAlloc) ! The master processor will not use these anymore... 
 !! the evaluation of collision integral will be performed in an ierarchical fashion: Each group of processors will be 
 !! united in a universe Acopy_UNIV. The first processor in the groupd will recieve rank = 0 and will be assignedas the master for
 !! this univers. The slave processors of a copy of Acopy_UNIV will send data to the master using collective communication.
 !! Also the processor with rank zero can still collect data frorm all toher processor and use broadcasts send the data out to 
 !! slave processors in all Acopy_UNIV universes 
 !! END OF THE STEP 1:
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:  
 !! Now each of the slave processors will read the A-arrays and prepare their local copies of A-arrays. 
 !! The slave processors will produce thier local A-capphi arrays. This information will be important on the next step: 
 !! the determination of each individual processor's job allocation in evaluation of the collision integral. 
  
 !! It is decided to let the slave processors divide the velocity nodes between the groups of processors
 !! essentially each will run acopy of the subroutine and should have the same answer. Then they will pick their portion of 
 !! nodes and will begin process them. In the end they will produce the workloads for themselves. Workload is the list of 
 !! velocity nodes for which the slave processor evaluates the collision operator 
 !!
 !! END OF  STEP 2
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !! Finally we setup some auxiliary arrays that we will need in the job assignment routin. 
 !! (We really just need the nodes_phican but we want to use the subroutine that sets all the arrays...
 !! first, we need to create the array A_capphi. it is not used on the master node, but we need a dummy array to 
 !! pass it to the subroutine so as to avoid mistakes.
 nn = size(nodes_uII,1)
 allocate (A_capphiII(1:nn), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_ReadOpA_highperfII: Allocation error for variable A_capphi. stop"
   stop
  end if
 A_capphiII=0 ! This array is not used on the master node
 !
 call SetCellsUGINdsSftArrsII(FindCellContainsPoint_DGV(0.0_DP,0.0_DP,0.0_DP)) 
 !! End of preparing the array nodes_phican
 !!
 !!!!!!!!!!!!!!!!!!!!! 
else
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!
 !! SLAVE PROCESS CODE
 !!
 !!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 1: 
 !! 
 !! The slave processor receives the A-array allocation that will be broadcasted from the master node: 
 call mpi_bcast (mpiI4BArryBuffer,(mpicommworldsize-1)*3,MPI_INTEGER4,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"PrepMPICollsnDGV_ReadOpA_highperfII: slave processor", irank, "MPI boradcast from proc 0 returned error", ierr
  stop
 end if   
 !! Now we need to retreve three numbers that pertain to this particular processor. 
 !! Specifically, we need to get the fstchnk=mpiCollisionWorkAlloc(irank,1), 
 !! ofst = mpiCollisionWorkAlloc(irank,2) and numrec = mpiCollisionWorkAlloc(irank,3) 
 !! The rest of the array is not important for this process. Therefore, 
 fstchnk=mpiI4BArryBuffer((irank-1)*3+1)
 ofst=mpiI4BArryBuffer((irank-1)*3+2)
 numrec=mpiI4BArryBuffer((irank-1)*3+3)
 !! END OF STEP 1:
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! STEP 2:
 !! We proceed processing the received information 
 !! Now each processor needs to read the corresponding part of the A array and get ready for work. 
 call ReadAarraysChnksOfstNrecDGVII(fstchnk,ofst,numrec,numchnksII)
 !! Now local copies of the A-arrays are set. We are ready to detemine the work load in terms of the velocity nodes 
 !! Need to run this setup subroutine: Presently we really just need the nodes_canphi but we want to use the subroutine that sets all the arrays...
 !! sets up arrays: nodes_Ashift,nodes_phican, nodes_dui,nodes_dvi,nodes_dwi
 call SetCellsUGINdsSftArrsII(FindCellContainsPoint_DGVII(0.0_DP,0.0_DP,0.0_DP)) 
 !! call the subrouine that evaluates the workload: velocity nodes that will be processed on this processor.
 call SetMPILocProcNodesWorkLoad_DGVII(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)
 !! the node work allocation is set ls
 
 !!
 !! END OF STEP 2.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! End distributing workload for the slave processors for evaluation of full Boltzmann operator.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
deallocate(mpiI4BArryBuffer)
end subroutine PrepMPICollsnDGV_ReadOpA_highperfII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine DivideRangeContinuous
!
! this subroutine will divide the range of integers from (nfrom) to (nto) into N peieces.
! the N is detemined by the length of the array divarr in the first dimension 
!
!

subroutine DivideRangeContinuous (divarr,nfrom,nto)

integer (I4B), dimension (:,:), intent (out) :: divarr ! the arrays that will store the pieces 
					! divarr(:,1) is the first number in the piece  
					! divarr(:,2) is the last number in the piece
integer (I4B), intent (in) :: nfrom  ! this is the first number in the range
integer (I4B), intent (in) :: nto    ! this is the last number in the range.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: numchunks, my_range, rem, i,avlength  ! various  scrap variables

! first, we check that the number of chunck is less than or equal to the entire range of numbers
divarr = 0
numchunks=size(divarr,2)
my_range = nto - nfrom+1
if (numchunks > my_range) then 
 print *,"DivideRangeContinuous: number of chunks is too big for the range, stop."
  stop
end if 
! if the number of chunks is less then or equal to the length of the range, we 
! start the division process: 
avlength = INT(my_range/REAL(numchunks,DP),I4B) ! this should be truncated
if (avlength == 0) then 
 print *, "DivideRangeContinuous: avlength=0 division is not successful. Stop"
 stop 
end if 
rem = my_range - numchunks*avlength
if ((rem < 0) .or. (rem >= numchunks)) then 
 print *, "DivideRangeContinuous: (rem < 0) .or. (rem >= numchunks). stop" 
 stop 
end if 
divarr(1,1)=nfrom
do i=1,rem
 divarr(2,i)=divarr(1,i)+avlength
 divarr(1,i+1) = divarr(2,i)+1
end do 
do i=rem+1,numchunks-1
 divarr(2,i)=divarr(1,i)+avlength-1
 divarr(1,i+1) = divarr(2,i)+1
end do
divarr(2,numchunks)=nto
end subroutine DivideRangeContinuous


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SetMPILocProcNodesWorkLoad_DGV(mpicommworldsize,irank)
!
! This subroutine sets up the array -- proc_nodes_wld
! this array lists all velocity nodes for which the 
! collision integral will be evaluated ON THIS SLAVE PROCESSOR. 
! 
! every row of chnks_copies correspond to a groupd of processors that shares an entire copy of 
! operator A among them. The subroutine will divide the total number of nodes into 
! a number equal chunks correponding to the number of such groups
! 
! within the chunk, the nodes will be distributed among the group so as 
! to provide the evaluation of the collision operator at those nodes. Essentially, 
! if a processor stores a portion of the A that corresponds to a node, this node is assgned to the processor.
! 
! Essentially, the code goes through all nodes and check if 
! the local A array has records that pertain to them. If it does the node 
! goes on the work list. To do that the code goes through all nodes. For each node it finds is canonical node
! from nodes_phican. Then it looks up the local A_capphi arrays and see if the corresponding canonical node is there. 
! IF so, the node goes on the list.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMPILocProcNodesWorkLoad_DGV(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)
use DGV_commvar, only: procs_nodes_wld,nodes_phican,A_capphi,num_Acopies,Acopy_Wrkld,&
                       procs_nodes_wld_Acopy_addr  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:) :: chnks_copies ! array that contains information about division of processors bewtenn different copies of operator A
integer (I4B), dimension (:,:) :: nodes_proc_groups ! array that contain information about division of velocity nodes between groups  of processors
integer, intent (in) :: mpicommworldsize ! the total number of processors in the common world
integer, intent (in) :: irank ! rank of this slave processr in the common world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer ! scrap array to keep found velocity nodes ...
integer (I4B) :: i,mm,nn,j,ii, loc_alloc_stat, nds_ct, My_Acopy_Univ_indx ! scrap index..

!!!!!!!!!!!!!!!!!!!!!!!!!
!! This slave processor belongs to a group that shares a copy of operator A
!! thos group of processors is assigned a bunch of velocity nodes to process. 
!! This subroutine will determine which velocity nodes from this bunch
!! shoul dbe processed on this slave processor. Essentially, it will go over the entire bunch
!! and if it has any components of operator A that are relevant to a velocity node, this node will be 
!! assigned to this processor. THus it is possible that all nodes in the bunch are assigned to each 
!! processor in the gorup. (this is true in the case of 1xN operator A
!! 
!! The first step is to detemrmine what nodes are assigned to this slave processor's group of processors
!! Note that this information could be obtained from the master node. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! for this subroutine we need to know the velocity nodes that are assigned for the processor's bunch.
!!! it is possible to receive these nodes from the master node. However, to avoid communication delays, we 
!!! re-compute the nodes locally. As long as all processors runds the same codes on same data, the results 
!!! should be identical

!!! We use arrays chnks_copies and nodes_proc_groups that are set in the the subroutine PrepMPICollsnDGV. 
 My_Acopy_Univ_indx=0
 do i=1,num_Acopies
  if ((irank >= chnks_copies(1,i)) .and. (irank<=chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! this contains the number of the Acopy group to which this processor belongs
  end if 
 end do 
 !!!! check if we found such number
 if (My_Acopy_Univ_indx==0) then 
  print *, "SetMPILocProcNodesWorkLoad: can not locate chunk to which this slave proc belongs. stop"
  stop
 end if !!!! end check
 !! Now we woill look at the nodes that are assigned to this processor's group.
 !! these nodes are start at nodes_proc_groups(1,My_Acopy_Univ_indx) and end at nodes_proc_groups(2,My_Acopy_Univ_indx)
 !! next we will go over these processors and see if operator A has any records for it. 
 !! we now allocate the buffer to hold the velocity nodes workload 
  mm = nodes_proc_groups(2,My_Acopy_Univ_indx)-nodes_proc_groups(1,My_Acopy_Univ_indx) + 1  ! the total number of nodes? assigned to the group.
  allocate (mpiI4BArryBuffer(1:mm), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
    print *, "SetMPIProcsNodesWorkLoad: Allocation error for variables mpiI4BArryBuffer on proc", irank," stop"
    stop
   end if
 !! end allocate 
 mpiI4BArryBuffer=0 ! clear the buffer
 !! Now we will go over each velcity node in the bunch and check if this processor stores any part of A for this processor.
 !! if it does, we will assign that node to the workload ...  
 nds_ct = 0 ! this is a counter for the found nodes... 
 do i = nodes_proc_groups(1,My_Acopy_Univ_indx),nodes_proc_groups(2,My_Acopy_Univ_indx) ! loop in the processors in the group j
  if ( A_capphi(nodes_phican(i))>0 ) then 
   nds_ct = nds_ct+1 !! coind the total number of found nodes...
   mpiI4BArryBuffer(nds_ct) = i
  end if 
 end do 
 !! Now we populate the local proc work artray:  
 !! First we check that this array exists already. If so, something went wrong and we stop
 if (size(procs_nodes_wld,1) > 0) then 
  print *,"SetMPIProcsNodesWorkLoad: Error. The array procs_nodes_wld already exists on proc", irank, "stop"
  stop 
 end if 
 ! Now if the array does not exist, we will create it and move the selected nodes to it:
 allocate (procs_nodes_wld(1:nds_ct),stat=loc_alloc_stat)
    !
 if (loc_alloc_stat >0) then 
  print *, "SetMPIProcsNodesWorkLoad:  Allocation error for variable (procs_nodes_wld). Proc", irank, "stop"
  stop
 end if
 ! next we populate the array procs_nodes_wld ... 
 procs_nodes_wld = mpiI4BArryBuffer(1:nds_ct) 
 ! the local work load array is set .... 
 deallocate (mpiI4BArryBuffer)
 
 !!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set up the arrays  
 !! Acopy_Wrkld and procs_nodes_wld_Acopy_addr
 !! These arrays only make sense on slave nodes...
 !!
 !!!!!!!!!!!!!!!!!!!!!!
 !! essentially for each Acopy universe, we need to write all the velocity nodes assigned to the Acopy group to 
 !! to Acopy_Wrkld. (Comment: in this context, "assigned velocity node" means the node at which collision integral is evaluated.)
 !! Importance of Acopy_Wrkld is the following. When results for evaluating the collision 
 !! operator are collected across processors sharing one copy of A, a.k.a. an Acopy group, these results are written in an array
 !! with results for different nodes arranged as in Acopy_Wrkld. 
 !! Individual processors in the Acopy group may have different nodes workloads. As a result, the locally computed results, 
 !! will be arranged in an array differently shaped than the results for the Group in total. When data is collected across 
 !! the Acopy group, the results need to be copied to correct places in the results array for the Group. For this, we 
 !! need a translation table array that tells for each entry of a localresult, that is arranged as procs_nodes_wld, where 
 !! the entry should go in the result array for the group which uses the same indexing of work nodes as in Acopy_Wrkld.  
  nn = nodes_proc_groups(2,My_Acopy_Univ_indx) - nodes_proc_groups(1,My_Acopy_Univ_indx) + 1 
  allocate (Acopy_Wrkld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (Acopy_Wrkld) on node", irank, "Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the Acopy_Wrkld
  do i = 1,nn 
   Acopy_Wrkld(i) = nodes_proc_groups(1,My_Acopy_Univ_indx) + i -1 !  Acopy_Wrkld(i) contains all velocity nodes assigned ot this Acopy group beginning 
                   ! from  nodes_proc_groups(1,My_Acopy_Univ_indx) and ending at  nodes_proc_groups(2,My_Acopy_Univ_indx)
  end do 
  !! end populating the Acopy_Wrkld
  !! now let us prepare the array procs_nodes_wld_Acopy_addr 
  nn = size(procs_nodes_wld,1)
  allocate (procs_nodes_wld_Acopy_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (procs_nodes_wld_Acopy_addr) on node",irank,"Stop"
   stop
  end if
  !! we go over the procs_nodes_wld array and find the matcning components in Acopy_Wrkld and record indices
  !! of these components in procs_nodes_wld_Acopy_addr
  procs_nodes_wld_Acopy_addr=0 ! reset, just in case, 
  do j=1,size(procs_nodes_wld,1)
   do ii=1,size(Acopy_Wrkld,1)
    if (Acopy_Wrkld(ii) == procs_nodes_wld(j)) then  
    procs_nodes_wld_Acopy_addr(j) = ii
    exit
    end if 
   end do 
  end do
  !! the procs_nodes_wld_Acopy_addr array is set
!
end subroutine SetMPILocProcNodesWorkLoad_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SetMPILocProcNodesWorkLoad_DGV_SelNds()
!
!
! This is a copy of the above subroutine except it is designed to use only a part of the primary velocity nodes. 
! The array of nodes used in simulation as a variable sel_nodes. The nodes listed in the array sel_nodes will be divided betwen processors 
! the numbers in the arrays nodes_proc_groups,procs_lin correspond to index in array sel_nodes and entries of sel_nodes are the numbers of the 
! node in nodes_u/_v/_w etc. arrays.
! Note that we do not need to see the arrya of nodes, we lust need to know how many are there, so we do not pass the array sel_nodes itself, just the number of records in it.
!
!
! This subroutine sets up the array -- proc_nodes_wld
! this array lists all velocity nodes for which the 
! collision integral will be evaluated ON THIS SLAVE PROCESSOR. 
! 
! every row of chnks_copies correspond to a groupd of processors that shares an entire copy of 
! operator A among them. The subroutine will divide the total number of nodes into 
! a number equal chunks correponding to the number of such groups
! 
! within the chunk, the nodes will be distributed among the group so as 
! to provide the evaluation of the collision operator at those nodes. Essentially, 
! if a processor stores a portion of the A that corresponds to a node, this node is assgned to the processor.
! 
! Essentially, the code goes through all nodes and check if 
! the local A array has records that pertain to them. If it does the node 
! goes on the work list. To do that the code goes through all nodes. For each node it finds is canonical node
! from nodes_phican. Then it looks up the local A_capphi arrays and see if the corresponding canonical node is there. 
! IF so, the node goes on the list.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMPILocProcNodesWorkLoad_DGV_SelNds(chnks_copies,nodes_proc_groups,mpicommworldsize,irank,sel_nodes)
use DGV_commvar, only: procs_nodes_wld,nodes_phican,A_capphi,num_Acopies,Acopy_Wrkld,&
                       procs_nodes_wld_Acopy_addr  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:) :: chnks_copies ! array that contains information about division of processors bewtenn different copies of operator A
integer (I4B), dimension (:,:) :: nodes_proc_groups ! array that contain information about division of velocity nodes between groups  of processors
integer, intent (in) :: mpicommworldsize ! the total number of processors in the common world
integer, intent (in) :: irank ! rank of this slave processr in the common world
integer (I4B), dimension (:), intent (in) :: sel_nodes ! this is the array that contains numbers of the selected velocity nodes on the entire grid. Arrays 
                    ! chnks_copies and nodes_proc_groups contain indices of entries in sel_nodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer ! scrap array to keep found velocity nodes ...
integer (I4B) :: i,mm,nn,j,ii, loc_alloc_stat, nds_ct, My_Acopy_Univ_indx ! scrap index..

!!!!!!!!!!!!!!!!!!!!!!!!!
!! This slave processor belongs to a group that shares a copy of operator A
!! thos group of processors is assigned a bunch of velocity nodes to process. 
!! This subroutine will determine which velocity nodes from this bunch
!! shoul dbe processed on this slave processor. Essentially, it will go over the entire bunch
!! and if it has any components of operator A that are relevant to a velocity node, this node will be 
!! assigned to this processor. THus it is possible that all nodes in the bunch are assigned to each 
!! processor in the gorup. (this is true in the case of 1xN operator A
!! 
!! The first step is to detemrmine what nodes are assigned to this slave processor's group of processors
!! Note that this information could be obtained from the master node. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! for this subroutine we need to know the velocity nodes that are assigned for the processor's bunch.
!!! it is possible to receive these nodes from the master node. However, to avoid communication delays, we 
!!! re-compute the nodes locally. As long as all processors runds the same codes on same data, the results 
!!! should be identical

!!! We use arrays chnks_copies and nodes_proc_groups that are set in the the subroutine PrepMPICollsnDGV. 
 My_Acopy_Univ_indx=0
 do i=1,num_Acopies
  if ((irank >= chnks_copies(1,i)) .and. (irank<=chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! this contains the number of the Acopy group to which this processor belongs
  end if 
 end do 
 !!!! check if we found such number
 if (My_Acopy_Univ_indx==0) then 
  print *, "SetMPILocProcNodesWorkLoad: can not locate chunk to which this slave proc belongs. stop"
  stop
 end if !!!! end check
 !! Now we woill look at the nodes that are assigned to this processor's group.
 !! these nodes are start at nodes_proc_groups(1,My_Acopy_Univ_indx) and end at nodes_proc_groups(2,My_Acopy_Univ_indx)
 !! next we will go over these processors and see if operator A has any records for it. 
 !! we now allocate the buffer to hold the velocity nodes workload 
  mm = nodes_proc_groups(2,My_Acopy_Univ_indx)-nodes_proc_groups(1,My_Acopy_Univ_indx) + 1  ! the total number of nodes? assigned to the group.
  allocate (mpiI4BArryBuffer(1:mm), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
    print *, "SetMPIProcsNodesWorkLoad: Allocation error for variables mpiI4BArryBuffer on proc", irank," stop"
    stop
   end if
 !! end allocate 
 mpiI4BArryBuffer=0 ! clear the buffer
 !! Now we will go over each velcity node in the bunch and check if this processor stores any part of A for this processor.
 !! if it does, we will assign that node to the workload ...  
 nds_ct = 0 ! this is a counter for the found nodes... 
 do i = nodes_proc_groups(1,My_Acopy_Univ_indx),nodes_proc_groups(2,My_Acopy_Univ_indx) ! loop in the processors in the group j
  if ( A_capphi(nodes_phican(i))>0 ) then 
   nds_ct = nds_ct+1 !! coind the total number of found nodes...
   mpiI4BArryBuffer(nds_ct) = sel_nodes(i)
  end if 
 end do 
 !! Now we populate the local proc work artray:  
 !! First we check that this array exists already. If so, something went wrong and we stop
 if (size(procs_nodes_wld,1) > 0) then 
  print *,"SetMPIProcsNodesWorkLoad: Error. The array procs_nodes_wld already exists on proc", irank, "stop"
  stop 
 end if 
 ! Now if the array does not exist, we will create it and move the selected nodes to it:
 allocate (procs_nodes_wld(1:nds_ct),stat=loc_alloc_stat)
    !
 if (loc_alloc_stat >0) then 
  print *, "SetMPIProcsNodesWorkLoad:  Allocation error for variable (procs_nodes_wld). Proc", irank, "stop"
  stop
 end if
 ! next we populate the array procs_nodes_wld ... 
 procs_nodes_wld = mpiI4BArryBuffer(1:nds_ct) 
 ! the local work load array is set .... 
 deallocate (mpiI4BArryBuffer)
 
 !!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set up the arrays  
 !! Acopy_Wrkld and procs_nodes_wld_Acopy_addr
 !! These arrays only make sense on slave nodes...
 !!
 !!!!!!!!!!!!!!!!!!!!!!
 !! essentially for each Acopy universe, we need to write all the velocity nodes assigned to the Acopy group to 
 !! to Acopy_Wrkld. (Comment: in this context, "assigned velocity node" means the node at which collision integral is evaluated.)
 !! Importance of Acopy_Wrkld is the following. When results for evaluating the collision 
 !! operator are collected across processors sharing one copy of A, a.k.a. an Acopy group, these results are written in an array
 !! with results for different nodes arranged as in Acopy_Wrkld. 
 !! Individual processors in the Acopy group may have different nodes workloads. As a result, the locally computed results, 
 !! will be arranged in an array differently shaped than the results for the Group in total. When data is collected across 
 !! the Acopy group, the results need to be copied to correct places in the results array for the Group. For this, we 
 !! need a translation table array that tells for each entry of a localresult, that is arranged as procs_nodes_wld, where 
 !! the entry should go in the result array for the group which uses the same indexing of work nodes as in Acopy_Wrkld.  
  nn = nodes_proc_groups(2,My_Acopy_Univ_indx) - nodes_proc_groups(1,My_Acopy_Univ_indx) + 1 
  allocate (Acopy_Wrkld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (Acopy_Wrkld) on node", irank, "Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the Acopy_Wrkld
  do i = 1,nn 
   Acopy_Wrkld(i) = sel_nodes(nodes_proc_groups(1,My_Acopy_Univ_indx) + i-1) !  Acopy_Wrkld(i) contains all velocity nodes assigned ot this Acopy group beginning 
                   ! from  nodes_proc_groups(1,My_Acopy_Univ_indx) and ending at  nodes_proc_groups(2,My_Acopy_Univ_indx)
  end do 
  !! end populating the Acopy_Wrkld
  !! now let us prepare the array procs_nodes_wld_Acopy_addr 
  nn = size(procs_nodes_wld,1)
  allocate (procs_nodes_wld_Acopy_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (procs_nodes_wld_Acopy_addr) on node",irank,"Stop"
   stop
  end if
  !! we go over the procs_nodes_wld array and find the matcning components in Acopy_Wrkld and record indices
  !! of these components in procs_nodes_wld_Acopy_addr
  procs_nodes_wld_Acopy_addr=0 ! reset, just in case, 
  do j=1,size(procs_nodes_wld,1)
   do ii=1,size(Acopy_Wrkld,1)
    if (Acopy_Wrkld(ii) == procs_nodes_wld(j)) then  
    procs_nodes_wld_Acopy_addr(j) = ii
    exit
    end if 
   end do 
  end do
  !! the procs_nodes_wld_Acopy_addr array is set
!
end subroutine SetMPILocProcNodesWorkLoad_DGV_SelNds


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SetMPILocProcNodesWorkLoad_DGVII(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)
!
! This is a copy of the above subroutine. Except it it modified to work 
! on the secondary mesh 
!
!
! This subroutine sets up the array -- proc_nodes_wld
! this array lists all velocity nodes for which the 
! collision integral will be evaluated ON THIS SLAVE PROCESSOR. 
! 
! every row of chnks_copies correspond to a groupd of processors that shares an entire copy of 
! operator A among them. The subroutine will divide the total number of nodes into 
! a number equal chunks correponding to the number of such groups
! 
! within the chunk, the nodes will be distributed among the group so as 
! to provide the evaluation of the collision operator at those nodes. Essentially, 
! if a processor stores a portion of the A that corresponds to a node, this node is assgned to the processor.
! 
! Essentially, the code goes through all nodes and check if 
! the local A array has records that pertain to them. If it does the node 
! goes on the work list. To do that the code goes through all nodes. For each node it finds is canonical node
! from nodes_phican. Then it looks up the local A_capphi arrays and see if the corresponding canonical node is there. 
! IF so, the node goes on the list.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMPILocProcNodesWorkLoad_DGVII(chnks_copies,nodes_proc_groups,mpicommworldsize,irank)

use DGV_commvar, only: procs_nodes_wld=>procs_nodes_wldII,nodes_phican=>nodes_phicanII,&
                       A_capphi=>A_capphiII,num_Acopies,Acopy_Wrkld=>Acopy_WrkldII,&
                       procs_nodes_wld_Acopy_addr=>procs_nodes_wld_Acopy_addrII  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:) :: chnks_copies ! array that contains information about division of processors bewtenn different copies of operator A
integer (I4B), dimension (:,:) :: nodes_proc_groups ! array that contain information about division of velocity nodes between groups  of processors
integer, intent (in) :: mpicommworldsize ! the total number of processors in the common world
integer, intent (in) :: irank ! rank of this slave processr in the common world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer ! scrap array to keep found velocity nodes ...
integer (I4B) :: i,mm,nn,j,ii, loc_alloc_stat, nds_ct, My_Acopy_Univ_indx ! scrap index..

!!!!!!!!!!!!!!!!!!!!!!!!!
!! This slave processor belongs to a group that shares a copy of operator A
!! thos group of processors is assigned a bunch of velocity nodes to process. 
!! This subroutine will determine which velocity nodes from this bunch
!! shoul dbe processed on this slave processor. Essentially, it will go over the entire bunch
!! and if it has any components of operator A that are relevant to a velocity node, this node will be 
!! assigned to this processor. THus it is possible that all nodes in the bunch are assigned to each 
!! processor in the gorup. (this is true in the case of 1xN operator A
!! 
!! The first step is to detemrmine what nodes are assigned to this slave processor's group of processors
!! Note that this information could be obtained from the master node. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! for this subroutine we need to know the velocity nodes that are assigned for the processor's bunch.
!!! it is possible to receive these nodes from the master node. However, to avoid communication delays, we 
!!! re-compute the nodes locally. As long as all processors runds the same codes on same data, the results 
!!! should be identical

!!! We use arrays chnks_copies and nodes_proc_groups that are set in the the subroutine PrepMPICollsnDGV. 
 My_Acopy_Univ_indx=0
 do i=1,num_Acopies
  if ((irank >= chnks_copies(1,i)) .and. (irank<=chnks_copies(2,i))) then 
   My_Acopy_Univ_indx = i ! this contains the number of the Acopy group to which this processor belongs
  end if 
 end do 
 !!!! check if we found such number
 if (My_Acopy_Univ_indx==0) then 
  print *, "SetMPILocProcNodesWorkLoadII: can not locate chunk to which this slave proc belongs. stop"
  stop
 end if !!!! end check
 !! Now we woill look at the nodes that are assigned to this processor's group.
 !! these nodes are start at nodes_proc_groups(1,My_Acopy_Univ_indx) and end at nodes_proc_groups(2,My_Acopy_Univ_indx)
 !! next we will go over these processors and see if operator A has any records for it. 
 !! we now allocate the buffer to hold the velocity nodes workload 
  mm = nodes_proc_groups(2,My_Acopy_Univ_indx)-nodes_proc_groups(1,My_Acopy_Univ_indx) + 1  ! the total number of nodes? assigned to the group.
  allocate (mpiI4BArryBuffer(1:mm), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
    print *, "SetMPIProcsNodesWorkLoadII: Allocation error for variables mpiI4BArryBuffer on proc", irank," stop"
    stop
   end if
 !! end allocate 
 mpiI4BArryBuffer=0 ! clear the buffer
 !! Now we will go over each velcity node in the bunch and check if this processor stores any part of A for this processor.
 !! if it does, we will assign that node to the workload ...  
 nds_ct = 0 ! this is a counter for the found nodes... 
 do i = nodes_proc_groups(1,My_Acopy_Univ_indx),nodes_proc_groups(2,My_Acopy_Univ_indx) ! loop in the processors in the group j
  if ( A_capphi(nodes_phican(i))>0 ) then 
   nds_ct = nds_ct+1 !! coind the total number of found nodes...
   mpiI4BArryBuffer(nds_ct) = i
  end if 
 end do 
 !! Now we populate the local proc work artray:  
 !! First we check that this array exists already. If so, something went wrong and we stop
 if (size(procs_nodes_wld,1) > 0) then 
  print *,"SetMPIProcsNodesWorkLoadII: Error. The array procs_nodes_wld already exists on proc", irank, "stop"
  stop 
 end if 
 ! Now if the array does not exist, we will create it and move the selected nodes to it:
 allocate (procs_nodes_wld(1:nds_ct),stat=loc_alloc_stat)
    !
 if (loc_alloc_stat >0) then 
  print *, "SetMPIProcsNodesWorkLoadII:  Allocation error for variable (procs_nodes_wld). Proc", irank, "stop"
  stop
 end if
 ! next we populate the array procs_nodes_wld ... 
 procs_nodes_wld = mpiI4BArryBuffer(1:nds_ct) 
 ! the local work load array is set .... 
 deallocate (mpiI4BArryBuffer)
 
 !!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set up the arrays  
 !! Acopy_Wrkld and procs_nodes_wld_Acopy_addr
 !! These arrays only make sense on slave nodes...
 !!
 !!!!!!!!!!!!!!!!!!!!!!
 !! essentially for each Acopy universe, we need to write all the velocity nodes assigned to the Acopy group to 
 !! to Acopy_Wrkld. (Comment: in this context, "assigned velocity node" means the node at which collision integral is evaluated.)
 !! Importance of Acopy_Wrkld is the following. When results for evaluating the collision 
 !! operator are collected across processors sharing one copy of A, a.k.a. an Acopy group, these results are written in an array
 !! with results for different nodes arranged as in Acopy_Wrkld. 
 !! Individual processors in the Acopy group may have different nodes workloads. As a result, the locally computed results, 
 !! will be arranged in an array differently shaped than the results for the Group in total. When data is collected across 
 !! the Acopy group, the results need to be copied to correct places in the results array for the Group. For this, we 
 !! need a translation table array that tells for each entry of a localresult, that is arranged as procs_nodes_wld, where 
 !! the entry should go in the result array for the group which uses the same indexing of work nodes as in Acopy_Wrkld.  
  nn = nodes_proc_groups(2,My_Acopy_Univ_indx) - nodes_proc_groups(1,My_Acopy_Univ_indx) + 1 
  allocate (Acopy_Wrkld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (Acopy_Wrkld) on node", irank, "Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the Acopy_Wrkld
  do i = 1,nn 
   Acopy_Wrkld(i) = nodes_proc_groups(1,My_Acopy_Univ_indx) + i -1 !  Acopy_Wrkld(i) contains all velocity nodes assigned ot this Acopy group beginning 
                   ! from  nodes_proc_groups(1,My_Acopy_Univ_indx) and ending at  nodes_proc_groups(2,My_Acopy_Univ_indx)
  end do 
  !! end populating the Acopy_Wrkld
  !! now let us prepare the array procs_nodes_wld_Acopy_addr 
  nn = size(procs_nodes_wld,1)
  allocate (procs_nodes_wld_Acopy_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (procs_nodes_wld_Acopy_addr) on node",irank,"Stop"
   stop
  end if
  !! we go over the procs_nodes_wld array and find the matcning components in Acopy_Wrkld and record indices
  !! of these components in procs_nodes_wld_Acopy_addr
  procs_nodes_wld_Acopy_addr=0 ! reset, just in case, 
  do j=1,size(procs_nodes_wld,1)
   do ii=1,size(Acopy_Wrkld,1)
    if (Acopy_Wrkld(ii) == procs_nodes_wld(j)) then  
    procs_nodes_wld_Acopy_addr(j) = ii
    exit
    end if 
   end do 
  end do
  !! the procs_nodes_wld_Acopy_addr array is set
!
end subroutine SetMPILocProcNodesWorkLoad_DGVII


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine SortNodesProcs 
! 
! this is a subroutine that sorts arrays of a special format. this soubroutine is useful to 
! create lin_send_array and lin_receive_array
!
! the (pairs) array has the following format
! # first element -- number of pairs  
! the first element is followed by zero or more pairs of the following format
! (#processor where the node will be send or from where it will be recieved, #node that will be sent)
! 
! it may be that the # processor wil appear more than once. Therefore, we will sort this array
! and place the resul in  array (grouped)
!
! The new format is 
! # first element is the toal number of records in this array (not counting itself)
! this is followed by zero or more chunks of the following format
! # processor where the node will be send or from where it will be recieved,
! # the toal number of nodes to be sent to this processor or to be received from this processor
! ##  numbers of the nodes that will be sent IN INCREASING ORDER! 
!
!!!!!!!!!!!!!!!!!!

subroutine SortNodesProcs(pairs,grouped) 

integer (I4B), dimension(:), intent(in)  :: pairs ! array of unsorted pairs
integer (I4B), dimension(:), intent(out) :: grouped ! array of unsorted pairs
!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension(1:size(pairs,1)) :: pairs_scr ! scrap array of unsorted pairs
integer (I4B) :: i,j,jj,jjj,nn,curr_proc,last_rec ! scrap indices
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pairs_scr=pairs ! make a scrap copy
grouped=0 ! clean the result
nn=0 ! reset the record counter - in the beginning there are no records 
last_rec=2*pairs_scr(1)+1 ! this is the total number of pairs in the array
! We now start to read the array
! a check if the array of pairs is too short: 
if ((size(pairs_scr,1)<last_rec) .or. (size(pairs_scr,1)<3)) then 
 print *,"SortNodesProcs: incoming array of pairs is too short for sorting. probably mistake. stop"
 stop
end if   
i=2 ! this points to the first number of the first pair
do while (i < last_rec) 
 !! we find first nonzero record. 
 curr_proc = pairs_scr(i)
 if (curr_proc /= 0) then ! 0 means no processror is assigned, 1=master node 
   !! now we save the new record for this processor in the sorted array.   
   nn = nn + 3   ! there will be two new records at the least 
   grouped(1) = nn   !  the first place holds the total number of non-empty records in this sorted array,
   pairs_scr(i) = 0  ! nullify the records to avoid dublication
   grouped(nn-1) = curr_proc ! this will remember the processor
   grouped(nn) = 1 ! the cell after the process holds the number of found nodes... 
   grouped(nn+1) = pairs_scr(i+1) ! this is the first node that we found for process curr_proc 
   pairs_scr(i+1) = 0 ! nullify the record toi avoid dublication
   !! next we will look if we have more records for the same processor.
   j=i+2 ! now j points to the next record....
   jj=0 ! this is the additional nodes for the same processor counter 
   do while  (j < last_rec)
     !! we find first nonzero record for the same processor. 
     if (pairs_scr(j) == curr_proc) then
      pairs_scr(j) = 0 ! nullify the first component of the pair to avaid dublication
      jj=jj+1 ! a new node found for this processor (curr_proc). Let us count this node in
      !! now we save the nonzero record in the sorted array.   
      grouped(nn) = grouped(nn)+1 ! count the node in in the prouped array:
      ! we also need to be sure that the nodes go in increasing order for each processor
      if (grouped(nn+jj) <= pairs_scr(j+1)) then ! we need to make sure that nodes are added in increasing order... 
       grouped(nn+1+jj) = pairs_scr(j+1)
      else ! if the new node is not bigger than the existing previous, then we need to do a sorting 
       jjj=0
       do while ((jjj < jj) .and. (grouped(nn+jj-jjj) > pairs_scr(j+1)))
        grouped(nn+1+jj-jjj) = grouped(nn+jj-jjj) ! push higher value to the right
        jjj=jjj+1
       end do 
       grouped(nn+1+jj-jjj) = pairs_scr(j+1)
      end if  
      pairs_scr(j+1) = 0 ! nullify the record to avoid duplication 
     end if  
     j=j+2 ! now the local index points to the next pair (it may be zero...)
   end do   
   nn=nn+jj ! if the internal loop found more records, let us shif the pointer so it points to the last node 
   grouped(1) = nn
 end if    
 i=i+2 ! shift the index to the next pair (the pair may have been zeroed by the previous code ...)
end do 
!
end subroutine SortNodesProcs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrouine InitDGV1D_MPI(nspatialcells)
!
! This is copy of the subroutine InitDGV1D from DGV-miscset
! adjusted for an MPI implementation 
!
!
! ATTENTION: you will need to knwo the total number of the spatial cells that the discretization will have 
! before you call this subroutine.
!
! This is a conglomerate subrouine that 
! performes initialization of the 
! parameters, contants and the variables 
! of the DGV Library
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitDGV1D_MPI(nspatcells,irank)

use DGV_readwrite

use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,nodes_v,&
                  nodes_w,nodes_gwts,nodes_uII,&
                  MvII,MuII,MwII,suII,svII,swII, &    
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,&
                  min_sens,Trad,run_mode,Cco1D,&
				  Order_nu,Order,SecondaryVelGridInUse,&
				  MoRatesArry1DAllocated,MoRatesArry1D,&
				  NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D,Nu1DUpdateArrysAllocated,&
				  TotNumSpatCells,MoRatesReliable1D,run_mode_1D,&
				  num_Acopies,num_lin_proc,lin_proc_nodes_wld,&
				  nodes_u_loc,nodes_v_loc,nodes_w_loc,nodes_gwts_loc
				  


use DGV_dgvtools_mod

use DGV_miscset

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    
 
!!!!!!!!!!!!!! Variables
integer (I4B), intent (in) :: nspatcells ! the total number of the spatial cells --- arrays will be created to store information for each spatial cell. 
                                         ! these arrays support velocity depemdent collision frequency model 
integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.

integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nn,i,j ! scrap variables

!!!!!!!!!!!!!! Variables to set up MPI parallelization algorithm
integer (I4B), dimension (:,:), allocatable :: chnks_copies !this arrays gives first and last processor in a group that share a copy of A
integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! this array will contain the breaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:,:), allocatable :: procs_lin ! this one will contain the nodes allocation to the group of processors consolidating the linearized operator. 
                          !it is a scrap variable and once the jobs are disctibuted, it will be deallocated...        


! variables for MPI Calls
integer ::  mpicommworldsize,ierr 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the next group of commands will initiate the key variables of the 
! DGV library. Change this commands if you need 
! different functionality of the library -- say if it is 
! desired to have a different initial velocity grid 
! or additional DGV variables need to be intitialized. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First we read parameters from the DGV library parameter file.
! the file name "DGVparametes.dat" is currently reserved for the DGV parameters. 
! However, this name can be changes
! by default, the file should be located in the same directory where the executable is started
!
if (irank==0) then 
 call SetDGVParams("DGVparameters.dat",0)
else
 call SetDGVParams("DGVparameters.dat",1)
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Mv=Mv_list(1)
Mu=Mu_list(1)
Mw=Mw_list(1)
su=su_list(1)
sv=sv_list(1)
sw=sw_list(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SetGDVGnodes    ! set some supplementary arrays to be used in building meshes
call SetDGVblzmmesh  ! sets one-dimensional grids in u,v, and w
! newt we build the velocity discretization 
call Set3DCellsR_DGV ! We create initial the velocity cells
call SetNodesDGV ! this will populate velocity cells with the 
! nodal - DG velocity nodes (see the description of the nodal-DG approach)

!!!!!!
! If you need to refine the velocity cells the following subroutines can be used: 
! in the next two lines the cell number 14 isdivived in 8 cubcells and in the next line
! the cell number 41 in the new mech is divided into 27 cubcells
!!! call CellsRefineDGV((/ 14 /),2,2,2) ! 
!!! call CellsRefineDGV((/ 27 + 14 /),3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!
! if one wants to record the cells, gridds and nodes information, use the following commands
! the parameters for the files names are taken from the DGV library parameter file
! the grids, cells, and nodes can be used in Matlab graphing subroutines. 
! they are usually not used in restart. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!
if (irank==0) then 
 call WriteGridsDGV
 call WriteCellsDGV
 call WriteNodesDGV
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines: 
! call WriteI1DGV(0.0_DP,0.0_DP,0.0_DP) ! this one will print the number of velocity nodes that 
                ! is contained in the cell that has the given velocity point
! call InitWriteDet    ! velocity dependent collision frequency. These subroutines only create an empty file. 
                ! later calls of related subroutines will save some diagnostic data in the created 
                ! files. These files make problems in many clusters becuase of the 
                ! writing protection. If the file is not deleted, it can not be overwritten by 
                ! a new instance of the running software -- we commented the use of the subroutines for now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CREATING NEW OPERATORS A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a call that usually is commented. 
! the next subroutines are an OpenMP subroutines
! that are included in this library and will create a 
! new instance of Operator A or Akor.  Specifically, 
! Operators A and Akor depends on the collision model and 
! the DG discretization. Note that subrouitnes creating the operator 
! A and Akor use the primary mesh. (In contrast, when we use pre-computed Operators A and Akor, 
! we can use eithr primary or secondary meshes) Once the (primary) DG discretization is set, 
! on can compute the Operators A or Akor. To do that, uncomment the corresponding line below
! Be aware that pre-comupting A and Akor is a very slow process and 
! usually should not be attempted on one processor for 
! meshes exceeding 33^3 velocity notes. 
!  
!!!!!!!!!!!!!!!!!
! For lagre meshes, one should use MPI versio of the 
! subroutine SetA_DGV . 
!
! TO BE ADDED LATER...  
!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Once the Operator A arrays are computed, 
! we need to save then on the hard drive in the 
! directory specified in the DGV parameter file
! To save Operator A on thehard Drive use the 
! the following subroutine 
!!!!!!!
! call WriteAarraysDGV
! call WriteAKorArraysDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI parallelization environment variables and setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! BEGIN quick fix
!! ATTENTION: the 1d3v code uses a lot of collective communications. It is required that all processors in the universe make 
!! the same communication call. Otherswise the code will get stuck. Portion of the libraries were programmed that only processors 
!! with ranks<num_lin_proc participate in the solution. This generate a conflict with making collective communications calls.
!! So, we either replace collective operators with individual (MPI_FILE_write_all with MPI_FILE_WRITE) or make sure that all 
!! processors are calling them at the same time.
!! we will go with the second solution, so we hardcoding a reset of num_lin_proc to be equal to the size of the MPI Common World
!!
!! CHANGED Communicator to MPI_LINEAR_COLL_WORLD in all of these == should be good for now, but will require calling PrepMPICollsnDGV_Univs_highperf
if (num_lin_proc/=mpicommworldsize) then
 num_lin_proc=mpicommworldsize
 if (irank==0) then 
  print*,"InitDGV1D_MPI: Re-setting num_lin_proc=mpicommworldsize:",mpicommworldsize
 end if 
end if 
!! END quick fix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Next we prepare three important arrays: num_chnks, nodes_proc_groups and procs_lin. 
!! These arrays determine the workload of nonlinear evaluation of 
!! the collision operator and the workload for evaluation of the linearized collision operator, 
!! respectively. 
!! First, let us work on the array chnks_copies(2,1:num_Acopies)
allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "InitDGV1D_MPI: Allocation error for variables chnks_copies. stop"
  stop
end if
chnks_copies=0 ! nullify, just in case.
!! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
!! in particular, a copy of A will be distributed between processors with number chnks_copies(1,1)=1 (master #0 does not get any) 
!! and chnks_copies(2,1)=XX. Another copy will be distributed between chnks_copies(1,2)=XX+1 and chnks_copies(2,2)=YY and so on.
!! 
call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
!! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
!! (the number fo the first and the last processors that contain this particular copy of A. num_Acopies  -- is the number of copies) 
!! 
!! Next we distribute workload between processor groups. It is assumed that since a copy of A is shared 
!! between processors in a group, all of them are needed to evaluate the collision operator for a single 
!! velocity point. This is true for s=1, but for s>1 not all processors in the group are needed to evaluate 
!! collision operator at a velocity node. This can be used later to tune up the work load a bit better later. 
!! At this moment, however, we will take all nodes where the colliision operator needs to be evaluated and we 
!! divide these velocity nodes between the groups of processors, thus creating a workload for each group.
!! Later this array will be used to assign individual workloads for 
!! each processor. 
nn = size(nodes_u,1) ! Total number of nodes numbered beginning from 1 where the collision integral needs to be evaluated.  
!! first we divide all velocity nodes between the groups of processors. We have exacly (num_Acopies) processor groups 
allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "InitDGV1D_MPI: Allocation error for variables nodes_proc_groups. stop"
 stop
end if
call DivideRangeContinuous(nodes_proc_groups,1,nn)
!! we divided velocity nodes between the groups of processors sharing copies of A.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! We make one more preparation for the future: we 
!! setup the array procs_lin that tells the workload for each processor involved in the evaluation of the 
!! linearized collision operator. Similarly to the above, simply divide velocity nodes between the 
!! processors involved in the evaluation of linearied collision operator.  
allocate (procs_lin(2,1:num_lin_proc), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "InitDGV1D_MPI: Allocation error for variable procs_lin. Stop"
 stop
end if
!! we populate the linearized group of processors with assgined nodes
call DivideRangeContinuous(procs_lin,1,nn)
!! end of creating proc_lin: all processors are divided in groups and each group is assgned some number of velocity nodes 
!! Note that if 1D-3D spatial operator is parallelized in velocity variable, proc_lin can be used to determine the work load. 
!! for each processor 
!! arrays procs_lin, nodes_proc_groups and chnks_copies have been created !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI FORK -- creating Linearization Workload Array -- is used to parallelize the transport part.
!!! The workload arrays are created only on processors participating in evaluation of linearied collision oparator
!!! Note that the same processor will be used to parallelize transport part in velocity space. The same breaking is uses
!!! for parallelization of the transport part as for evaluation linearized operator
if (irank < num_lin_proc) then 
  nn = procs_lin(2,irank+1) - procs_lin(1,irank+1)+1 ! this is the total number of velocity nodes that will be assigned to this processor 
  allocate (lin_proc_nodes_wld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variable (lin_proc_nodes_wld) on  node", irank, ". Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the lin_proc_nodes_wld
  do i = procs_lin(1,irank+1),procs_lin(2,irank+1)
   lin_proc_nodes_wld(i-procs_lin(1,irank+1)+1) = i ! lin_proc_nodes_wld contains all velocity nodes (actually, their numbers) beginning from procs_lin(1,irank+1) and ending at procs_lin(2,irank+1)
  end do 
  !! end populating the lin_proc_nodes_wld 
  allocate (nodes_u_loc(1:size(lin_proc_nodes_wld,1)),nodes_v_loc(1:size(lin_proc_nodes_wld,1)),&
          nodes_w_loc(1:size(lin_proc_nodes_wld,1)),nodes_gwts_loc(1:size(lin_proc_nodes_wld,1)), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variable (nodes_u_loc,nodes_v_loc,nodes_w_loc,nodes_gwts)"
  end if 
  !!!!!!!!!!!!!!!!
  ! prepare copies of the velocity nodes to be used on this processor
  do i=1,size(lin_proc_nodes_wld,1)
   j=lin_proc_nodes_wld(i)
   nodes_u_loc(i)=nodes_u(j)
   nodes_v_loc(i)=nodes_v(j)
   nodes_w_loc(i)=nodes_w(j)
   nodes_gwts_loc(i)=nodes_gwts(j)
  end do  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! UNCOMMENT BELOW IF USING Boltzmann Collision on primary mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment the following subroutine if desire to have communicators between groups of processors sharing a copy of A
!!! and groups of processors using information from processors sharing as copy of A.
!!! iT ALSO SETS setting us workload/send/recieve/ arrays to make processor work together. 
call PrepMPICollsnDGV_Univs_highperf(irank,chnks_copies,nodes_proc_groups,procs_lin)
!!! Uncomment the following subroutine if desire to use evaluation of the collision operator using primary mesh.   
!!! The subrouine implements reading appropiate portions of operator A into memory of processors 
! call PrepMPICollsnDGV_ReadOpA_highperf(irank,chnks_copies,nodes_proc_groups)
deallocate(chnks_copies,nodes_proc_groups,procs_lin)
!!!
!!!
!! END READ pre-computed operators A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following few lines set up necessary variables to work with the Boltzmann collision integral. 
! We note that in this code, the Botlzmann collision integral is using secondary mesh. This means that
! the unknown solution will have to be mapped between the primary and secondary meshes. The secondary mesh, 
! in principle can be the same as primary -- then no interpolation is needed. However, in general, secondary mesh is coarser a
! and is use to evaluate the full Boltzmann Collision integral (it might be that the task is just too expensive on the primary mesh
!
! Because the collisio operator will use secondary mesh we check the flag SecondaryVelGridInUse whether 
! the secondary mesh is in use. If it is, we create the mesh and then may read the operator A
!
if (SecondaryVelGridInUse) then 
! We create variables of the secondary mesh
! A quick check of we have paparameters of the secondary mesh defined: 
 if ((size(Mu_list,1)<2) .or. (size(Mv_list,1)<2) .or.(size(Mw_list,1)<2) .or. &
           (size(su_list,1)<2) .or. (size(sv_list,1)<2) .or. (size(sw_list,1)<2) ) then 
  ! Error -- secondary mesh requested, but at least one parameter is missing for it 
  print *, "InitDGV1D_MPI: Secondary velocity mesh has been requested, but at least one parameter is missing for secondary mesh."  
  stop
  ! 
 end if
 ! set up the parameters of the secondary nodal-DG discretization. Secondary and primary discretizations use the same 
 ! boundaries in the velocity space
 MvII=Mv_list(2)
 MuII=Mu_list(2)
 MwII=Mw_list(2)
 suII=su_list(2)
 svII=sv_list(2)
 swII=sw_list(2)
 ! ready to create secondary meshes 
 ! NOTE: the suffix II in the name indicates that the subroutine works with the second mesh 
 call SetDGVblzmmeshII  ! sets one-dimensional grids in u,v, and w 
 ! newt we build the velocity discretization 
 call Set3DCellsR_DGVII ! We create initial velocity cells
 call SetNodesDGVII ! this will populate velocity cells with the 
 ! nodal - DG velocity nodes (see the description of the nodal-DG approach)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Initializing a PRE-COMPUTED Operator A
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Note that Pre-computed Operator A wil be used for evaluating collision integral and estimating relaxation rates. 
 ! These rates will be used in the velocity dependent BGK model. For this model evaluatin will be done using secondary grid. 
 ! Therefore one needs to indicate that secondary velocoity grid is used.
 ! For that the flag SecondaryVelGridInUse needs to be set true in the DGVparameters file.
 ! 
 ! It may happen that the secondary grid is the same as the primary, but in general we 
 ! will expect that the solution is interpolated between the grids
 !!!!!!!!!!!!!!! 

 ! The evaluation of the Boltzmann collision operator in the method of 
 !  Alekseenko-Josyula requires the 
 ! knowledge of pre-computed Operator A. 
 ! There has been a library of instances of 
 ! Operator A computed for different parameters of DG discretizations
 ! The pre-computed operator is loaded using the next batch of commands:
 !
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Next we prepare three important arrays: num_chnks, nodes_proc_groups and procs_lin. 
 ! These arrays determine the workload of nonlinear evaluation of 
 ! the collision operator and the workload for evaluation of the linearized collision operator, 
 ! respectively. 
 ! First, let us work on the array2 chnks_copies(2,1:num_Acopies)
 allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variables chnks_copies. stop"
   stop
 end if
 chnks_copies=0 ! nullify, just in case.
 !! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
 !! in particular, a copy of A will be distributed between processors with number chnks_copies(1,1)=1 (master #0 does not get any) 
 !! and chnks_copies(2,1)=XX. Another copy will be distributed between chnks_copies(1,2)=XX+1 and chnks_copies(2,2)=YY and so on.
 !! 
 call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
 !! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
 !! (the number fo the first and the last processors that contain this particular copy of A. num_Acopies  -- is the number of copies) 
 !! 
 !! Next we distribute workload between processor groups. It is assumed that since a copy of A is shared 
 !! between processors in a group, all of them are needed to evaluate the collision operator for a single 
 !! velocity point. This is true for s=1, but for s>1 not all processors in the group are needed to evaluate 
 !! collision operator at a velocity node. This can be used later to tune up the work load a bit better later. 
 !! At this moment, however, we will take all nodes where the colliision operator needs to be evaluated and we 
 !! divide these velocity nodes between the groups of processors, thus creating a workload for each group.
 !! Later this array will be used to assign individual workloads for 
 !! each processor. 
 nn = size(nodes_uII,1) ! Total number of nodes numbered beginning from 1 where the collision integral needs to be evaluated.  
 !! first we divide all velocity nodes between the groups of processors. We have exacly (num_Acopies) processor groups 
 allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "InitDGV1D_MPI: Allocation error for variables nodes_proc_groups. stop"
  stop
 end if
 call DivideRangeContinuous(nodes_proc_groups,1,nn)
 !! we divided velocity nodes between the groups of processors sharing copies of A.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! UNCOMMENT BELOW IF USING Boltzmann Collision on Secondary Mesh
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! Uncomment the following subroutine if desire to have communicators between groups of processors sharing a copy of A
 !!! and groups of processors using information from processors sharing as copy of A.
 !!! iT ALSO SETS setting us workload/send/recieve/ arrays to make processor work together. 
 call PrepMPICollsnDGV_Univs_highperfII(irank,chnks_copies,nodes_proc_groups)
 !!! Uncomment the following subroutine if desire to use evaluation of the collision operator using primary mesh.   
 !!! The subrouine implements reading appropiate portions of operator A into memory of processors 
 call PrepMPICollsnDGV_ReadOpA_highperfII(irank,chnks_copies,nodes_proc_groups)
 deallocate(chnks_copies,nodes_proc_groups)
 !
 ! Another preparatory call. This one sets up array (nodes_primcellII) used for projecting the primary solution on the secondary mesh
 call prepare_sDGVpMap_DGV ! this subroutine nodes_primcellII
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
! End of the section dedicated to secondary mesh and end of preparatory work for
! Using full Botlzmann collision operator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


TotNumSpatCells=nspatcells ! record the total number of cells for future checks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Atention: this array may be used in the BGK model with velocity dependent collision frequency. 
!! However, currently, the version of the subrouine does not make use if this array: therefore its allocation is commented. 
!
! allocate the arrays used for preserving the coeeficients of the velocity dependent collision frequency
!allocate (Cco1D(Order_nu,nspatcells), stat=loc_alloc_stat)
!if (loc_alloc_stat >0) then
! print *, "InitDGV1D: Allocation error for (Cco1D)"
!end if
!Cco1D=0 ! reset the coefficients to zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! next we check if the storage for the moments rates that are used in the velocity dependent collision frequency has been created yet
if (.not. MoRatesArry1DAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (MoRatesArry1D(1:Order,nspatcells), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Error allocation arrays to store coefficient of the VDCF (MoRatesArry0D)"
   stop
 end if 
 MoRatesArry1DAllocated = .true.
 ! Also, we will set up the other arrays to make sure that the coefficients are updated:
 MoRatesArry1D = 0 ! reset the array
end if 
! Allocate arrays that will keep the flag is the relaxation rates are reliable in the spatial cell: 
allocate (MoRatesReliable1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D_MPI: Allocation error for (MoRatesReliable1D)"
end if
MoRatesReliable1D = .false. ! reset the coefficients to false
!
! Allocate arrays that will keep the run_mode for ech spatial cell: 
allocate (run_mode_1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D_MPI: Allocation error for (run_mode_1D)"
end if
run_mode_1D = 0 ! reset the coefficients to mode=0 = full Boltzmann!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we check if the storage has been created for information about the last update of the relaxation rates
if (.not. Nu1DUpdateArrysAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (NuNextUpdateTime1D(nspatcells),NuLastDens1D(nspatcells),NuLastTemp1D(nspatcells), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "InitDGV1D_MPI: Error allocation arrays to store information about last rates apdate: (NuNextUpdateTime1D)," //  &
               "(NuLastDens1D),(NuLastTemp1D)" 
       stop
      end if 
Nu1DUpdateArrysAllocated = .true.
! Also, we will initialize the arrays to make sure that the coefficients are updated:
NuNextUpdateTime1D = - 100000.0   ! these bogus values will set off the criteria for update
NuLastDens1D = -100000.0 
NuLastTemp1D = -100000.0 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!    THIS VARIABLE SHOULD NOT BE USED --- CHECK AND REMOVE LATER>
!! INITIALIZE THE RUN_MODE to ZERO -- in the beginning the regime is strongly non-linear
  run_mode=0
!!!!!!!!!!!!!!!!!!!!!!!

end subroutine InitDGV1D_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrouine InitDGV1DsymW_MPI(nspatialcells)
!
! This is a copy of the above subroutine adjusted for the case when the solution is 
! symmetric in W direction in the velocity variable. In this case, nodes with W<0 can be removed and the full solution can be 
! restored from the nodes with W>0
!
! In this subroutine, the only nodes with non-negative w-component are used in MPI parallelization. 
! Specifically, the lin_proc_nodes_wld only keeps nodes with non-negative w-components and arrays 
! nodes_u_loc, nodes_v_loc, nodes_w_loc and nodes_gwts_loc are populated using only noves 
! with non-negative w-components. 
! 
!
! This is copy of the subroutine InitDGV1D from DGV-miscset
! adjusted for an MPI implementation 
!
!
! ATTENTION: you will need to knwo the total number of the spatial cells that the discretization will have 
! before you call this subroutine.
!
! This is a conglomerate subrouine that 
! performes initialization of the 
! parameters, contants and the variables 
! of the DGV Library
!
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine InitDGV1DsymW_MPI(nspatcells,irank)

use DGV_readwrite

use DGV_commvar, only: Mv,Mu,Mw,su,sv,sw,nodes_u,nodes_v,&
                  nodes_w,nodes_gwts,nodes_uII,&
                  MvII,MuII,MwII,suII,svII,swII, &    
                  Mu_list,Mv_list,Mw_list,su_list,sv_list,sw_list,&
                  min_sens,Trad,run_mode,Cco1D,&
				  Order_nu,Order,SecondaryVelGridInUse,&
				  MoRatesArry1DAllocated,MoRatesArry1D,&
				  NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D,Nu1DUpdateArrysAllocated,&
				  TotNumSpatCells,MoRatesReliable1D,run_mode_1D,&
				  num_Acopies,num_lin_proc,lin_proc_nodes_wld,&
				  nodes_u_loc,nodes_v_loc,nodes_w_loc,nodes_gwts_loc,nodes_symfctr_loc
				  


use DGV_dgvtools_mod

use DGV_miscset

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!                    
 
!!!!!!!!!!!!!! Variables
integer (I4B), intent (in) :: nspatcells ! the total number of the spatial cells --- arrays will be created to store information for each spatial cell. 
                                         ! these arrays support velocity depemdent collision frequency model 
integer, intent (in) :: irank !  MPI variable -- the number of the processor that runs this code.

integer (I4B) :: loc_alloc_stat ! to keep allocation status
integer (I4B) :: nn,i,j,num_sel_nds ! scrap variables
integer (I4B) :: iu,iv,iw,iiu,iiv,iiw !scrap variables

!!!!!!!!!!!!!! Variables to set up MPI parallelization algorithm
integer (I4B), dimension (:,:), allocatable :: chnks_copies !this arrays gives first and last processor in a group that share a copy of A
integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! this array will contain the breaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b           
integer (I4B), dimension (:,:), allocatable :: procs_lin ! this one will contain the nodes allocation to the group of processors consolidating the linearized operator. 
                          !it is a scrap variable and once the jobs are disctibuted, it will be deallocated...        
integer (I4B), dimension (:), allocatable :: nodes_nnegw_num ! scratch array to keep the numbers of nodes with non-negative w-component

! variables for MPI Calls
integer ::  mpicommworldsize,ierr 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the next group of commands will initiate the key variables of the 
! DGV library. Change this commands if you need 
! different functionality of the library -- say if it is 
! desired to have a different initial velocity grid 
! or additional DGV variables need to be intitialized. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First we read parameters from the DGV library parameter file.
! the file name "DGVparametes.dat" is currently reserved for the DGV parameters. 
! However, this name can be changes
! by default, the file should be located in the same directory where the executable is started
!
if (irank==0) then 
 call SetDGVParams("DGVparameters.dat",0)
else
 call SetDGVParams("DGVparameters.dat",1)
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Mv=Mv_list(1)
Mu=Mu_list(1)
Mw=Mw_list(1)
su=su_list(1)
sv=sv_list(1)
sw=sw_list(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call SetGDVGnodes    ! set some supplementary arrays to be used in building meshes
call SetDGVblzmmesh  ! sets one-dimensional grids in u,v, and w
! newt we build the velocity discretization 
call Set3DCellsR_DGV ! We create initial the velocity cells
call SetNodesDGV ! this will populate velocity cells with the 
! nodal - DG velocity nodes (see the description of the nodal-DG approach)

!!!!!!
! If you need to refine the velocity cells the following subroutines can be used: 
! in the next two lines the cell number 14 isdivived in 8 cubcells and in the next line
! the cell number 41 in the new mech is divided into 27 cubcells
!!! call CellsRefineDGV((/ 14 /),2,2,2) ! 
!!! call CellsRefineDGV((/ 27 + 14 /),3,3,3)

!!!!!!!!!!!!!!!!!!!!!!!!
! if one wants to record the cells, gridds and nodes information, use the following commands
! the parameters for the files names are taken from the DGV library parameter file
! the grids, cells, and nodes can be used in Matlab graphing subroutines. 
! they are usually not used in restart. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!
if (irank==0) then 
 call WriteGridsDGV
 call WriteCellsDGV
 call WriteNodesDGV
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Additional subroutines: 
! call WriteI1DGV(0.0_DP,0.0_DP,0.0_DP) ! this one will print the number of velocity nodes that 
                ! is contained in the cell that has the given velocity point
! call InitWriteDet    ! velocity dependent collision frequency. These subroutines only create an empty file. 
                ! later calls of related subroutines will save some diagnostic data in the created 
                ! files. These files make problems in many clusters becuase of the 
                ! writing protection. If the file is not deleted, it can not be overwritten by 
                ! a new instance of the running software -- we commented the use of the subroutines for now
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CREATING NEW OPERATORS A and Akor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is a call that usually is commented. 
! the next subroutines are an OpenMP subroutines
! that are included in this library and will create a 
! new instance of Operator A or Akor.  Specifically, 
! Operators A and Akor depends on the collision model and 
! the DG discretization. Note that subrouitnes creating the operator 
! A and Akor use the primary mesh. (In contrast, when we use pre-computed Operators A and Akor, 
! we can use eithr primary or secondary meshes) Once the (primary) DG discretization is set, 
! on can compute the Operators A or Akor. To do that, uncomment the corresponding line below
! Be aware that pre-comupting A and Akor is a very slow process and 
! usually should not be attempted on one processor for 
! meshes exceeding 33^3 velocity notes. 
!  
!!!!!!!!!!!!!!!!!
! For lagre meshes, one should use MPI versio of the 
! subroutine SetA_DGV . 
!
! TO BE ADDED LATER...  
!
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Once the Operator A arrays are computed, 
! we need to save then on the hard drive in the 
! directory specified in the DGV parameter file
! To save Operator A on thehard Drive use the 
! the following subroutine 
!!!!!!!
! call WriteAarraysDGV
! call WriteAKorArraysDGV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI parallelization environment variables and setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Check the size of the MPI Universs = the total number of processors in the Universe.
call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
!! BEGIN quick fix
!! ATTENTION: the 1d3v code uses a lot of collective communications. It is required that all processors in the universe make 
!! the same communication call. Otherswise the code will get stuck. Portion of the libraries were programmed that only processors 
!! with ranks<num_lin_proc participate in the solution. This generate a conflict with making collective communications calls.
!! So, we either replace collective operators with individual (MPI_FILE_write_all with MPI_FILE_WRITE) or make sure that all 
!! processors are calling them at the same time.
!! we will go with the second solution, so we hardcoding a reset of num_lin_proc to be equal to the size of the MPI Common World
!!
!! CHANGED Communicator to MPI_LINEAR_COLL_WORLD in all of these == should be good for now, but will require calling PrepMPICollsnDGV_Univs_highperf
if (num_lin_proc/=mpicommworldsize) then
 num_lin_proc=mpicommworldsize
 if (irank==0) then 
  print*,"InitDGV1D_MPI: Re-setting num_lin_proc=mpicommworldsize:",mpicommworldsize
 end if 
end if 
!! END quick fix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Next we prepare three important arrays: num_chnks, nodes_proc_groups and procs_lin. 
!! These arrays determine the workload of nonlinear evaluation of 
!! the collision operator and the workload for evaluation of the linearized collision operator, 
!! respectively. 
!! First, let us work on the array chnks_copies(2,1:num_Acopies)
allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "InitDGV1D_MPI: Allocation error for variables chnks_copies. stop"
  stop
end if
chnks_copies=0 ! nullify, just in case.
!! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
!! in particular, a copy of A will be distributed between processors with number chnks_copies(1,1)=1 (master #0 does not get any) 
!! and chnks_copies(2,1)=XX. Another copy will be distributed between chnks_copies(1,2)=XX+1 and chnks_copies(2,2)=YY and so on.
!! 
call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
!! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
!! (the number fo the first and the last processors that contain this particular copy of A. num_Acopies  -- is the number of copies) 
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!
!! next we single out the nodes with non-negative w-component. 
allocate (nodes_nnegw_num(1:size(nodes_w,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "InitDGV1DsymW_MPI: Allocation error for variables nodes_nnegw_num. stop"
  stop
end if
!!!!!!!!!!!!!!!!! uncomment the next 7 lines if it is desired to have nodes enumerated using native velocity mesh.
!!num_sel_nds=0
!!do i=1,size(nodes_w,1)
!! if (nodes_w(i)> -1.0D-11) then 
!!  num_sel_nds=num_sel_nds+1
!!  nodes_nnegw_num(num_sel_nds)=i ! now nodes_nnegw_num contains numbers of the nodes with non-negative w-component
!! end if  
!!end do
!!!!!!!!!!!!!!!!!! end of the select nodes fragments
!!!!!!!!!!!!!!!!!! uncomment the next 7 lines if desire to enumerate nodes consistent with smoke
num_sel_nds = 0
do iu=1,Mu
 do iiu=1,su
  do iv=1,Mv
   do iiv=1,sv
    do iw=1,Mw
     do iiw=1,sw
      i= ((iu-1)*Mv*Mw+(iv-1)*Mw+iw-1)*su*sv*sw + ((iiu-1)*sv*sw+(iiv-1)*sw+iiw)
      if (nodes_w(i)> -1.0D-11) then 
       num_sel_nds=num_sel_nds+1
       nodes_nnegw_num(num_sel_nds)=i ! now nodes_nnegw_num contains numbers of the nodes with non-negative w-component
      end if 
     end do
    end do 
   end do     
  end do  
 end do  
end do
!!!!!!!!!!!!!!!!!! end of the select nodes fragments
!! note that currently nn is equal to the total number of nodeswith non-negative w-component. The variable will be overwrtiiten later.
!!
!! Next we distribute workload between processor groups. It is assumed that since a copy of A is shared 
!! between processors in a group, all of them are needed to evaluate the collision operator for a single 
!! velocity point. This is true for s=1, but for s>1 not all processors in the group are needed to evaluate 
!! collision operator at a velocity node. This can be used later to tune up the work load a bit better later. 
!! At this moment, however, we will take all nodes where the colliision operator needs to be evaluated and we 
!! divide these velocity nodes between the groups of processors, thus creating a workload for each group.
!! Later this array will be used to assign individual workloads for 
!! each processor. 
!! first we divide all velocity nodes between the groups of processors. We have exacly (num_Acopies) processor groups 
allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "InitDGV1D_MPI: Allocation error for variables nodes_proc_groups. stop"
 stop
end if
call DivideRangeContinuous(nodes_proc_groups,1,num_sel_nds)
!! we divided velocity nodes between the groups of processors sharing copies of A.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! We make one more preparation for the future: we 
!! setup the array procs_lin that tells the workload for each processor involved in the evaluation of the 
!! linearized collision operator. Similarly to the above, simply divide velocity nodes between the 
!! processors involved in the evaluation of linearied collision operator.  
allocate (procs_lin(2,1:num_lin_proc), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "InitDGV1D_MPI: Allocation error for variable procs_lin. Stop"
 stop
end if
!! we populate the linearized group of processors with assgined nodes
call DivideRangeContinuous(procs_lin,1,num_sel_nds)
!! end of creating proc_lin: all processors are divided in groups and each group is assgned some number of velocity nodes 
!! Note that if 1D-3D spatial operator is parallelized in velocity variable, proc_lin can be used to determine the work load. 
!! for each processor 
!! arrays procs_lin, nodes_proc_groups and chnks_copies have been created !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI FORK -- creating Linearization Workload Array -- is used to parallelize the transport part.
!!! The workload arrays are created only on processors participating in evaluation of linearied collision oparator
!!! Note that the same processor will be used to parallelize transport part in velocity space. The same breaking is uses
!!! for parallelization of the transport part as for evaluation linearized operator
if (irank < num_lin_proc) then 
  nn = procs_lin(2,irank+1) - procs_lin(1,irank+1)+1 ! this is the total number of velocity nodes that will be assigned to this processor 
  allocate (lin_proc_nodes_wld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variable (lin_proc_nodes_wld) on  node", irank, ". Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the lin_proc_nodes_wld
  do i = procs_lin(1,irank+1),procs_lin(2,irank+1)
   lin_proc_nodes_wld(i-procs_lin(1,irank+1)+1) = nodes_nnegw_num(i) ! lin_proc_nodes_wld contains all velocity nodes (actually, their numbers) beginning from procs_lin(1,irank+1) and ending at procs_lin(2,irank+1)
  end do 
  !! end populating the lin_proc_nodes_wld 
  allocate (nodes_u_loc(1:size(lin_proc_nodes_wld,1)),nodes_v_loc(1:size(lin_proc_nodes_wld,1)),&
          nodes_w_loc(1:size(lin_proc_nodes_wld,1)),nodes_gwts_loc(1:size(lin_proc_nodes_wld,1)),&
          nodes_symfctr_loc(1:size(lin_proc_nodes_wld,1)), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variable (nodes_u_loc,nodes_v_loc,nodes_w_loc,nodes_gwts)"
  end if 
  !!!!!!!!!!!!!!!!
  ! prepare copies of the velocity nodes to be used on this processor
  do i=1,size(lin_proc_nodes_wld,1)
   j=lin_proc_nodes_wld(i)
   nodes_u_loc(i)=nodes_u(j)
   nodes_v_loc(i)=nodes_v(j)
   nodes_w_loc(i)=nodes_w(j)
   nodes_gwts_loc(i)=nodes_gwts(j)
   if (nodes_w(j)>1.0d-11) then 
     nodes_symfctr_loc(i) = 2
   else 
     nodes_symfctr_loc(i) = 1
   end if     
  end do  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! UNCOMMENT BELOW IF USING Boltzmann Collision on primary mesh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Uncomment the following subroutine if desire to have communicators between groups of processors sharing a copy of A
!!! and groups of processors using information from processors sharing as copy of A.
!!! iT ALSO SETS setting us workload/send/recieve/ arrays to make processor work together. 
call PrepMPICollsnDGV_Univs_highperf_SelNds(irank,chnks_copies,nodes_proc_groups,procs_lin,num_sel_nds)
!!! Uncomment the following subroutine if desire to use evaluation of the collision operator using primary mesh.   
!!! The subrouine implements reading appropiate portions of operator A into memory of processors 
! call PrepMPICollsnDGV_ReadOpA_highperf_SelNds(irank,chnks_copies,nodes_proc_groups,nodes_nnegw_num(1:nodes_nnegw_num))
deallocate(chnks_copies,nodes_proc_groups,procs_lin)
!!!
!!!
!! END READ pre-computed operators A
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The following few lines set up necessary variables to work with the Boltzmann collision integral. 
! We note that in this code, the Botlzmann collision integral is using secondary mesh. This means that
! the unknown solution will have to be mapped between the primary and secondary meshes. The secondary mesh, 
! in principle can be the same as primary -- then no interpolation is needed. However, in general, secondary mesh is coarser a
! and is use to evaluate the full Boltzmann Collision integral (it might be that the task is just too expensive on the primary mesh
!
! Because the collisio operator will use secondary mesh we check the flag SecondaryVelGridInUse whether 
! the secondary mesh is in use. If it is, we create the mesh and then may read the operator A
!
if (SecondaryVelGridInUse) then 
! We create variables of the secondary mesh
! A quick check of we have paparameters of the secondary mesh defined: 
 if ((size(Mu_list,1)<2) .or. (size(Mv_list,1)<2) .or.(size(Mw_list,1)<2) .or. &
           (size(su_list,1)<2) .or. (size(sv_list,1)<2) .or. (size(sw_list,1)<2) ) then 
  ! Error -- secondary mesh requested, but at least one parameter is missing for it 
  print *, "InitDGV1D_MPI: Secondary velocity mesh has been requested, but at least one parameter is missing for secondary mesh."  
  stop
  ! 
 end if
 ! set up the parameters of the secondary nodal-DG discretization. Secondary and primary discretizations use the same 
 ! boundaries in the velocity space
 MvII=Mv_list(2)
 MuII=Mu_list(2)
 MwII=Mw_list(2)
 suII=su_list(2)
 svII=sv_list(2)
 swII=sw_list(2)
 ! ready to create secondary meshes 
 ! NOTE: the suffix II in the name indicates that the subroutine works with the second mesh 
 call SetDGVblzmmeshII  ! sets one-dimensional grids in u,v, and w 
 ! newt we build the velocity discretization 
 call Set3DCellsR_DGVII ! We create initial velocity cells
 call SetNodesDGVII ! this will populate velocity cells with the 
 ! nodal - DG velocity nodes (see the description of the nodal-DG approach)

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Initializing a PRE-COMPUTED Operator A
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Note that Pre-computed Operator A wil be used for evaluating collision integral and estimating relaxation rates. 
 ! These rates will be used in the velocity dependent BGK model. For this model evaluatin will be done using secondary grid. 
 ! Therefore one needs to indicate that secondary velocoity grid is used.
 ! For that the flag SecondaryVelGridInUse needs to be set true in the DGVparameters file.
 ! 
 ! It may happen that the secondary grid is the same as the primary, but in general we 
 ! will expect that the solution is interpolated between the grids
 !!!!!!!!!!!!!!! 

 ! The evaluation of the Boltzmann collision operator in the method of 
 !  Alekseenko-Josyula requires the 
 ! knowledge of pre-computed Operator A. 
 ! There has been a library of instances of 
 ! Operator A computed for different parameters of DG discretizations
 ! The pre-computed operator is loaded using the next batch of commands:
 !
 !!!!!!!!!!!!!!!!!!!!!!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Next we prepare three important arrays: num_chnks, nodes_proc_groups and procs_lin. 
 ! These arrays determine the workload of nonlinear evaluation of 
 ! the collision operator and the workload for evaluation of the linearized collision operator, 
 ! respectively. 
 ! First, let us work on the array2 chnks_copies(2,1:num_Acopies)
 allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Allocation error for variables chnks_copies. stop"
   stop
 end if
 chnks_copies=0 ! nullify, just in case.
 !! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
 !! in particular, a copy of A will be distributed between processors with number chnks_copies(1,1)=1 (master #0 does not get any) 
 !! and chnks_copies(2,1)=XX. Another copy will be distributed between chnks_copies(1,2)=XX+1 and chnks_copies(2,2)=YY and so on.
 !! 
 call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
 !! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
 !! (the number fo the first and the last processors that contain this particular copy of A. num_Acopies  -- is the number of copies) 
 !! 
 !! Next we distribute workload between processor groups. It is assumed that since a copy of A is shared 
 !! between processors in a group, all of them are needed to evaluate the collision operator for a single 
 !! velocity point. This is true for s=1, but for s>1 not all processors in the group are needed to evaluate 
 !! collision operator at a velocity node. This can be used later to tune up the work load a bit better later. 
 !! At this moment, however, we will take all nodes where the colliision operator needs to be evaluated and we 
 !! divide these velocity nodes between the groups of processors, thus creating a workload for each group.
 !! Later this array will be used to assign individual workloads for 
 !! each processor. 
 nn = size(nodes_uII,1) ! Total number of nodes numbered beginning from 1 where the collision integral needs to be evaluated.  
 !! first we divide all velocity nodes between the groups of processors. We have exacly (num_Acopies) processor groups 
 allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "InitDGV1D_MPI: Allocation error for variables nodes_proc_groups. stop"
  stop
 end if
 call DivideRangeContinuous(nodes_proc_groups,1,nn)
 !! we divided velocity nodes between the groups of processors sharing copies of A.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! UNCOMMENT BELOW IF USING Boltzmann Collision on Secondary Mesh
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! Uncomment the following subroutine if desire to have communicators between groups of processors sharing a copy of A
 !!! and groups of processors using information from processors sharing as copy of A.
 !!! iT ALSO SETS setting us workload/send/recieve/ arrays to make processor work together. 
 call PrepMPICollsnDGV_Univs_highperfII(irank,chnks_copies,nodes_proc_groups)
 !!! Uncomment the following subroutine if desire to use evaluation of the collision operator using primary mesh.   
 !!! The subrouine implements reading appropiate portions of operator A into memory of processors 
 call PrepMPICollsnDGV_ReadOpA_highperfII(irank,chnks_copies,nodes_proc_groups)
 deallocate(chnks_copies,nodes_proc_groups)
 !
 ! Another preparatory call. This one sets up array (nodes_primcellII) used for projecting the primary solution on the secondary mesh
 call prepare_sDGVpMap_DGV ! this subroutine nodes_primcellII
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end if 
! End of the section dedicated to secondary mesh and end of preparatory work for
! Using full Botlzmann collision operator. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


TotNumSpatCells=nspatcells ! record the total number of cells for future checks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Atention: this array may be used in the BGK model with velocity dependent collision frequency. 
!! However, currently, the version of the subrouine does not make use if this array: therefore its allocation is commented. 
!
! allocate the arrays used for preserving the coeeficients of the velocity dependent collision frequency
!allocate (Cco1D(Order_nu,nspatcells), stat=loc_alloc_stat)
!if (loc_alloc_stat >0) then
! print *, "InitDGV1D: Allocation error for (Cco1D)"
!end if
!Cco1D=0 ! reset the coefficients to zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! next we check if the storage for the moments rates that are used in the velocity dependent collision frequency has been created yet
if (.not. MoRatesArry1DAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (MoRatesArry1D(1:Order,nspatcells), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
   print *, "InitDGV1D_MPI: Error allocation arrays to store coefficient of the VDCF (MoRatesArry0D)"
   stop
 end if 
 MoRatesArry1DAllocated = .true.
 ! Also, we will set up the other arrays to make sure that the coefficients are updated:
 MoRatesArry1D = 0 ! reset the array
end if 
! Allocate arrays that will keep the flag is the relaxation rates are reliable in the spatial cell: 
allocate (MoRatesReliable1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D_MPI: Allocation error for (MoRatesReliable1D)"
end if
MoRatesReliable1D = .false. ! reset the coefficients to false
!
! Allocate arrays that will keep the run_mode for ech spatial cell: 
allocate (run_mode_1D(nspatcells), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "InitDGV1D_MPI: Allocation error for (run_mode_1D)"
end if
run_mode_1D = 0 ! reset the coefficients to mode=0 = full Boltzmann!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! next we check if the storage has been created for information about the last update of the relaxation rates
if (.not. Nu1DUpdateArrysAllocated) then 
 ! we need to create the array for storing the coefficients of the velocity dependent collision frequency: 
 allocate (NuNextUpdateTime1D(nspatcells),NuLastDens1D(nspatcells),NuLastTemp1D(nspatcells), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "InitDGV1D_MPI: Error allocation arrays to store information about last rates apdate: (NuNextUpdateTime1D)," //  &
               "(NuLastDens1D),(NuLastTemp1D)" 
       stop
      end if 
Nu1DUpdateArrysAllocated = .true.
! Also, we will initialize the arrays to make sure that the coefficients are updated:
NuNextUpdateTime1D = - 100000.0   ! these bogus values will set off the criteria for update
NuLastDens1D = -100000.0 
NuLastTemp1D = -100000.0 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!    THIS VARIABLE SHOULD NOT BE USED --- CHECK AND REMOVE LATER>
!! INITIALIZE THE RUN_MODE to ZERO -- in the beginning the regime is strongly non-linear
  run_mode=0
!!!!!!!!!!!!!!!!!!!!!!!
end subroutine InitDGV1DsymW_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UniversalCollisionOperator1DonecellDGV(f,fcol,time,dt,cellnum)
! 
! This is a copy of UniversalCollisionOperator1DonecellDGV adjusted for MPI parallelization
! in this set up the work of solving transport part is split between processors by dividing 
! the components of the transport part between them. Since each component correspond to one velocity node in the nodal-DG setup, 
! division is accomplished by distributing evenly the velocity nodes between processors with rank< num_lin_proc
! 
! Note that in this MPI set up the right side only is evaluated at the velocity nodes that are assigned to that processor. 
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

subroutine UniversalCollisionOperator1DonecellDGV_MPI(f,fcol,time,dt,cellnum)

use DGV_commvar, only: run_mode_1D,mol_diam,L_inf,N_inf,T_inf,C_inf,gas_viscosity,gasR
use DGV_dgvtools_mod
use DGV_collision_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


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
real (DP) :: LocDens, LocUbar, LocVbar, LocWbar, LocTemp  ! scrap variables to keep the local macroparamters
real (DP) :: L1_err ! L1 norm of the difference between the solution and the local maxwellian
logical :: silent
! WARNING: Make sure that Nodes_Ashift, Nodes_ccan and other supplementary arrays are set before calling the 
! evaluation of the collision operator.
!!!!!!!! MPI related
integer :: irank,ierr 
!!!!!!!!!!!! set us the allocatable array !
allocate (Df(1:size(fcol,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "UniversalCollisionOperator1DonecellDGV_MPI: Allocation error for variables (Df)"
 stop
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on... 
if (irank==0) then 
 silent = .false.
else
 silent = .true.
end if
!!!!!!!!!!!!!!! Uncomment next line if using full velocity mesh 
call CheckSolutionModeSimple1Donecell_DGV_MPI(f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTemp,L1_err,cellnum,silent) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.
!!!!!!!!!!!!!!! Uncomment next line if using mesh with symmetry in W component
!call CheckSolutionModeSimple1DonecellSymW_DGV_MPI(f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTemp,L1_err,cellnum,silent) ! This will check the local distribution function against a maxwellian with the same 5 macroparameters.
! debug
if ((irank==0) .and. (cellnum==1)) then 
 print *,"time=",time,"run_mode_1D=", run_mode_1D
end if
! end debug 
run_mode = run_mode_1D(cellnum)
!
select case (run_mode) ! run_mode is set outside, when solution is evaluated for closeness to a maxwellian...  
	case (0) ! run_mode=0 means we are very far from a Maxwellian. In this case we just call the collision operator
	  !!!!!!!!!!!!!!!!!!!!! Stopper ...
	  if (irank==0) then
       print *," mode 0 is not implemented. Stop"
      end if
      stop
      !!!!!!!!!!!!!!!!!!!!!
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
	  !!!!!!!!!!!!!!!!!!!!! Stopper ...
	  if (irank==0) then
       print *," mode 1 is not implemented. Stop"
      end if
      stop
      !!!!!!!!!!!!!!!!!!!!!
	  
	  !!!!!!!!!!!! set us the allocatable array !
      allocate (fcol_scr(1:size(fcol,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "UniversalCollisionOperator1DonecellDGV_MPI: Allocation error for variables (fcol_scr)"
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
	  !!!!!!!!!!!!!!!!!!!!! Stopper ...
	  if (irank==0) then
       print *," mode 2 is not implemented. Stop"
      end if
      stop
      !!!!!!!!!!!!!!!!!!!!!
	  
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
	    ! call EvalColESBGK_1DMPI(f,fcol,LocDens,LocUbar,LocVbar,LocWbar,LocTemp) ! UNCOMMENT IF USING THE ES MODEL
		call EvalColGrad_1Donecell_MPI(f,Df,fcol,time,LocDens,LocUbar,LocVbar,LocWbar,LocTemp,L1_err,cellnum)
		fcol = fcol*coef_temp*dt
	case (4) ! ES-BGK mode or Shakhov  
		coef_temp = (kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity
		! Uncomment one of these operators if using full mesh 
		call EvalColESBGK_1DMPI(f,fcol,LocDens,LocUbar,LocVbar,LocWbar,LocTemp) ! UNCOMMENT IF USING THE ES MODEL
		! Uncomment one of these operators if using Symmetry in W-component (half-mesh)
		! call EvalColESBGKSymW_1DMPI(f,fcol,LocDens,LocUbar,LocVbar,LocWbar,LocTemp) ! UNCOMMENT IF USING THE ES MODEL
		fcol = fcol*coef_temp*dt
	case default 
		print *, "UniversalCollisionOperator1DonecellDGV_MPI: cannot process value of run_mode", run_mode
		stop
end select

deallocate (Df)

end subroutine UniversalCollisionOperator1DonecellDGV_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionModeSimple1Donecell_DGV_MPI(f,Df,cellnum,slt) 
! 
! This is a copy of the subroutine CheckSolutionModeSimple1Donecell_DGV_MPI adjusted for MPI parallelization...
!
! Subroutine returns values of local macroparameters as byproducts LocDens,LocUbar,LocVbar,LocWbar,LocTemp,
!
! Subroutine return value of the L1 difference between solution and the local maxwellian
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

subroutine CheckSolutionModeSimple1Donecell_DGV_MPI(f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTemp,L1_err,cellnum,slt)

use DGV_commvar, only: run_mode_1D, nodes_u_loc, nodes_v_loc, nodes_w_loc, &
                   decomp_lev, linear_lev, vel_lev, ES_lev, nodes_gwts_loc, &
                   MPI_LINEAR_COLL_WORLD 

use DGV_distributions_mod
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out):: Df ! the difference between the solution and the local Maxwellian
integer (I4B), intent(in) :: cellnum ! the number of the cell for which run_mode is calculated.
logical, intent (in) :: slt !if =true, no printout is generated
real (DP), intent (out) :: LocDens, LocUbar, LocVbar, LocWbar, LocTemp ! scratch variables to keep macroparameters, 
real (DP), intent (out) :: L1_err
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: n, u_0,v_0,w_0
real (DP), dimension (size(f,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err_fm ! a sctratch variables for the errors
integer (I4B) :: run_mode ! a scratch variable to keep the run mode.
!!!!!!!!!!!!!!!!!!!!!!!
! MPI-related 
real (DP), dimension (:), allocatable :: send_buff, recv_buff ! send and recveive buffers
integer :: loc_alloc_stat, ierr
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! first we compute the density integral of the portion of solution that is stored on this processor
LocDens = sum(f*nodes_gwts_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! also check momentum integrals
LocUbar=sum(f*nodes_gwts_loc*nodes_u_loc)
LocVbar=sum(f*nodes_gwts_loc*nodes_v_loc)
LocWbar=sum(f*nodes_gwts_loc*nodes_w_loc)
!!! and the local temperature integrals, that is a part of it 
LocTemp = sum(f*nodes_gwts_loc*(nodes_u_loc**2+nodes_v_loc**2+nodes_w_loc**2))/3.0_DP*2.0_DP ! part of dimensionless temperature
!!! next we need to combine these over all processors to have the total value of the density 
allocate(send_buff(5),recv_buff(5), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "CSMSimple1Donecell_DGV_MPI: Error allocation arrays for buffers (send_buiff,recv_buff). Stop" 
 stop
end if 
send_buff(1)=LocDens
send_buff(2)=LocUbar
send_buff(3)=LocVbar
send_buff(4)=LocWbar
send_buff(5)=LocTemp
recv_buff=0 ! clean
! Sum the macroparameters across all processors
call MPI_ALLREDUCE (send_buff,recv_buff,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"CSMSimple1Donecell_DGV_MPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
!! compute the values of the macroparameters based on collective results
LocDens=recv_buff(1)
LocUbar=recv_buff(2)/LocDens
LocVbar=recv_buff(3)/LocDens
LocWbar=recv_buff(4)/LocDens
LocTemp=recv_buff(5)/LocDens-(LocUbar**2+LocVbar**2+LocWbar**2)/3.0_DP*2.0_DP
!!! Now that we know macroparameters, we need to compute the L1-norm of the deviation of the local solution from the local maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
fMaxwellNew = maxwelveldist(LocTemp,LocUbar,LocVbar,LocWbar,LocDens,nodes_u_loc, nodes_v_loc, nodes_w_loc) ! now we populate the maxwellian with the same macroparamters.
Df = f-fMaxwellNew ! perturbation from the local maxwellian (evaluated only on the velocity nodes assigned to this processor)
L1_err = SUM(ABS(Df)*nodes_gwts_loc)/LocDens! evaluate the portion of the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!! 
! Gather the L1_err on all processors.
!!!!
send_buff=0
recv_buff=0
send_buff(1)=L1_err
call mpi_allreduce (send_buff,recv_buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"CSMSimple1Donecell_DGV_MPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
L1_err=recv_buff(1)
!!!!!!!!!!!!!!!!!!!!!!!
deallocate(send_buff,recv_buff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Print May comment
!if (.not. slt) then 
! print*, "CSMSimple1Donecell_DGV_MPI: L1 error = ", L1_err, "cellnum=", cellnum
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
run_mode = run_mode_1D(cellnum) ! restore the old runmode. 
   select case (run_mode)
   case (0) ! if we're in a strongly non-linear regime on the previous time step, then do this: 
      if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then 
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from strongly non-linear"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		if (.not. slt) then 
 		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from strongly non-linear regime"
 		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2                          ! linear regime
        if (.not. slt) then 
 		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linearized regime from strongly non-linear"
 		end if 
      elseif (L1_err < decomp_lev) then
	    run_mode=1         
	    if (.not. slt) then             ! decomposition regime
         print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from strongly non-linear"
        end if  
      end if
!
  case (1) ! if already we're in the decomposition regime
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from decomposition regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from decomposition regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then                          ! linear regime
         print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linearized regime from decomposition regime"
        end if  
      elseif (L1_err < decomp_lev) then
        run_mode=1  
      else
	    run_mode=0
	    if (.not. slt) then 
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from decomposition regime"
		end if 
	  end if
!
   case (2) ! if we're in the linear regime..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from linear regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from linear regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Determine the new regime:
        run_mode=2
      elseif (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
         print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from linear regime"
        end if  
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from linear regime"
		end if  
	  end if
!
	case (3) ! vel dep ES-BGK mode
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to ES-BGK mode from vel-dependent regime"
		end if 
	  elseif (L1_err < vel_lev) then
	    run_mode=3
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to linear mode from vel-dependent regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to decomposition regime from vel-dep regime"
		end if 
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from vel-dep regime"
		end if 
	  end if
!
	case (4) ! if we're in the ES-BGK mode..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel-dep mode from ES regime"
		end if 
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2  
        if (.not. slt) then
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linear reg from ES regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from ES regime"
		end if  
      else
        run_mode=0
		if (.not. slt) then
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from ES regime"
		end if 
	  end if

   case default
      print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Error, can not process the value of run_mode=", run_mode
      stop
   end select
run_mode_1D(cellnum) = run_mode  ! save the new runmode.   

end subroutine CheckSolutionModeSimple1Donecell_DGV_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine CheckSolutionModeSimple1DonecellSymW_DGV_MPI(f,Df,cellnum,slt) 
!
! This is a copy of the subroutine CheckSolutionModeSimple1Donecell_DGV_MPI adjusted to symmetry in W component of velocity 
! (even in w). It is expected that only non-negative nodes in W were kept. 
!
! 
! This is a copy of the subroutine CheckSolutionModeSimple1Donecell_DGV_MPI adjusted for MPI parallelization...
!
! Subroutine returns values of local macroparameters as byproducts LocDens,LocUbar,LocVbar,LocWbar,LocTemp,
!
! Subroutine return value of the L1 difference between solution and the local maxwellian
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

subroutine CheckSolutionModeSimple1DonecellSymW_DGV_MPI(f,Df,LocDens,LocUbar,LocVbar,LocWbar,LocTemp,L1_err,cellnum,slt)

use DGV_commvar, only: run_mode_1D, nodes_u_loc, nodes_v_loc, nodes_w_loc, &
                   decomp_lev, linear_lev, vel_lev, ES_lev, nodes_gwts_loc, &
                   nodes_symfctr_loc, MPI_LINEAR_COLL_WORLD 

use DGV_distributions_mod
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

Intrinsic SUM, ABS

real (DP), dimension (:), intent(in) :: f ! the solution on the current state
real (DP), dimension (:), intent(out):: Df ! the difference between the solution and the local Maxwellian
integer (I4B), intent(in) :: cellnum ! the number of the cell for which run_mode is calculated.
logical, intent (in) :: slt !if =true, no printout is generated
real (DP), intent (out) :: LocDens, LocUbar, LocVbar, LocWbar, LocTemp ! scratch variables to keep macroparameters, 
real (DP), intent (out) :: L1_err
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real (DP) :: n, u_0,v_0,w_0
real (DP), dimension (size(f,1)) :: fMaxwellNew  ! a scratch storage for the maxwellian 
real (DP) :: L1_err_fm ! a sctratch variables for the errors
integer (I4B) :: run_mode,i ! a scratch variable to keep the run mode.
!!!!!!!!!!!!!!!!!!!!!!!
! MPI-related 
real (DP), dimension (:), allocatable :: send_buff, recv_buff ! send and recveive buffers
integer :: loc_alloc_stat, ierr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! first we compute the density integral of the portion of solution that is stored on this processor
LocDens = sum(f*nodes_gwts_loc*nodes_symfctr_loc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! also check momentum integrals
LocUbar=sum(f*nodes_gwts_loc*nodes_u_loc*nodes_symfctr_loc)
LocVbar=sum(f*nodes_gwts_loc*nodes_v_loc*nodes_symfctr_loc)
LocWbar=0
!!! and the local temperature integrals, that is a part of it 
LocTemp = sum(f*nodes_gwts_loc*(nodes_u_loc**2+nodes_v_loc**2+nodes_w_loc**2)*nodes_symfctr_loc)/3.0_DP*2.0_DP ! part of dimensionless temperature
!!! next we need to combine these over all processors to have the total value of the density 
allocate(send_buff(5),recv_buff(5), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "CSMSimple1Donecell_DGV_MPI: Error allocation arrays for buffers (send_buiff,recv_buff). Stop" 
 stop
end if 
send_buff(1)=LocDens
send_buff(2)=LocUbar
send_buff(3)=LocVbar
send_buff(4)=LocWbar
send_buff(5)=LocTemp
recv_buff=0 ! clean
! Sum the macroparameters across all processors
call MPI_ALLREDUCE (send_buff,recv_buff,5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"CSMSimple1Donecell_DGV_MPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
!! compute the values of the macroparameters based on collective results
LocDens=recv_buff(1)
LocUbar=recv_buff(2)/LocDens
LocVbar=recv_buff(3)/LocDens
LocWbar=0
LocTemp=recv_buff(5)/LocDens-(LocUbar**2+LocVbar**2+LocWbar**2)/3.0_DP*2.0_DP
!!! Now that we know macroparameters, we need to compute the L1-norm of the deviation of the local solution from the local maxwellian
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
fMaxwellNew = maxwelveldist(LocTemp,LocUbar,LocVbar,LocWbar,LocDens,nodes_u_loc, nodes_v_loc, nodes_w_loc) ! now we populate the maxwellian with the same macroparamters.
Df = f-fMaxwellNew ! perturbation from the local maxwellian (evaluated only on the velocity nodes assigned to this processor)
L1_err = SUM(ABS(Df)*nodes_gwts_loc*nodes_symfctr_loc)/LocDens! evaluate the portion of the relative L1-norm of the difference (L1 norm of the maxwellian should be density)
!!!! 
! Gather the L1_err on all processors.
!!!!
send_buff=0
recv_buff=0
send_buff(1)=L1_err
call mpi_allreduce (send_buff,recv_buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"CSMSimple1Donecell_DGV_MPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
L1_err=recv_buff(1)
!!!!!!!!!!!!!!!!!!!!!!!
deallocate(send_buff,recv_buff)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostic Print May comment
!if (.not. slt) then 
! print*, "CSMSimple1Donecell_DGV_MPI: L1 error = ", L1_err, "cellnum=", cellnum
!end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
run_mode = run_mode_1D(cellnum) ! restore the old runmode. 
   select case (run_mode)
   case (0) ! if we're in a strongly non-linear regime on the previous time step, then do this: 
      if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then 
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from strongly non-linear"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
		if (.not. slt) then 
 		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from strongly non-linear regime"
 		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2                          ! linear regime
        if (.not. slt) then 
 		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linearized regime from strongly non-linear"
 		end if 
      elseif (L1_err < decomp_lev) then
	    run_mode=1         
	    if (.not. slt) then             ! decomposition regime
         print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from strongly non-linear"
        end if  
      end if
!
  case (1) ! if already we're in the decomposition regime
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
		if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from decomposition regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from decomposition regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then                          ! linear regime
         print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linearized regime from decomposition regime"
        end if  
      elseif (L1_err < decomp_lev) then
        run_mode=1  
      else
	    run_mode=0
	    if (.not. slt) then 
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from decomposition regime"
		end if 
	  end if
!
   case (2) ! if we're in the linear regime..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to ES-BGK mode from linear regime"
		end if  
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel ES-BGK mode from linear regime"
		end if  
      elseif (L1_err < linear_lev) then     ! Determine the new regime:
        run_mode=2
      elseif (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
         print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from linear regime"
        end if  
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from linear regime"
		end if  
	  end if
!
	case (3) ! vel dep ES-BGK mode
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	    if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to ES-BGK mode from vel-dependent regime"
		end if 
	  elseif (L1_err < vel_lev) then
	    run_mode=3
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2
        if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to linear mode from vel-dependent regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum,  "Switching to decomposition regime from vel-dep regime"
		end if 
      else
	    run_mode=0
	    if (.not. slt) then
		 print*,"CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from vel-dep regime"
		end if 
	  end if
!
	case (4) ! if we're in the ES-BGK mode..
	  if (L1_err < ES_lev) then             ! at the lowest L1-error, we switch to ES-BGK mode
	    run_mode=4
	  elseif (L1_err < vel_lev) then
	    run_mode=3
	    if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to vel-dep mode from ES regime"
		end if 
      elseif (L1_err < linear_lev) then     ! Detemine the new regime:
        run_mode=2  
        if (.not. slt) then
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to linear reg from ES regime"
		end if 
      else if (L1_err < decomp_lev) then
        run_mode=1                          ! linear regime
        if (.not. slt) then
		 print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to decomposition regime from ES regime"
		end if  
      else
        run_mode=0
		if (.not. slt) then
		 print*, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Switching to strongly non-linear from ES regime"
		end if 
	  end if

   case default
      print *, "CSMSimple1Donecell_DGV_MPI: cell=",cellnum, "Error, can not process the value of run_mode=", run_mode
      stop
   end select
run_mode_1D(cellnum) = run_mode  ! save the new runmode.   

end subroutine CheckSolutionModeSimple1DonecellSymW_DGV_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! EvalColESBGK_1DMPI
!
! This is a copy of the above subroutine adjusted for MPI implementation of 1D kineitc code. 
! The  difference is that this subroutine only evaluated the collision operator on nodes assigned to this particular processor
!
!
! Right hand side calculation for the ES-BGK distribution and collision operator
! Here, the RHS = nu * (f0 - f)
! where
! nu = collision frequency
! f0 = ESBGK distribution
! f = solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColESBGK_1DMPI(f,RHS,n,u_0,v_0,w_0,T)

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_commvar, only: nodes_u_loc, nodes_v_loc, nodes_w_loc, &
				   alpha, gas_viscosity, gas_T_reference, gas_alpha, C_inf, gasR

real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: RHS ! the value of the collision operator for each component of the solution.

real (DP), intent (in) :: u_0, v_0, w_0 ! bulk velocities
real (DP), intent (in) :: n ! density
real (DP), intent (in) :: T ! temperature
!real (DP) :: Determinant ! the determinant of the Tensor
real (DP) :: Pressure
real (DP) :: nu ! this is the collision frequency term
real (DP), dimension (:), allocatable :: f0 ! Distribution function

real (DP), parameter :: kBoltzmann = 1.3806503D-23
real (DP), dimension (3,3) :: Tensor, TensorInv ! both the tensor and the inverse tensor for the ES-BGK
real (DP) :: Determinant ! the determinant of the tensor
integer :: loc_alloc_stat


! the tensor and its corresponding inverse and determinant is computed here
call MakeTensor_1DMPI (f,Tensor,n,u_0,v_0,w_0,T) ! alpha only in common variabls
TensorInv = inv(Tensor)
Determinant = DetTensor(Tensor)

allocate (f0(1:size(f,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalColESBGK_1dmpi: Allocation error for variables (f0)"
 stop
end if

f0 = ESBGK_f0(TensorInv,Determinant,n,u_0,v_0,w_0,nodes_u_loc,nodes_v_loc,nodes_w_loc)

! now to evaluate the collision requency term
Pressure = n*T ! dimensionless Pressure is computed here
nu = Pressure/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/T)**gas_alpha ! final dimensionless nu ! have gas_T_reference be dimensionless?
RHS = nu * (f0 - f)

deallocate (f0)
!
end subroutine EvalColESBGK_1DMPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! EvalColESBGKSymW_1DMPI
!
!
! This is a copy of EvalColESBGK_1DMPI adjusted for the case of symmetry in W component of velocity  
!
! This is a copy of EvalColESBGK_1D adjusted for MPI implementation of 1D kineitc code. 
! The  difference is that this subroutine only evaluated the collision operator on nodes assigned to this particular processor
!
!
! Right hand side calculation for the ES-BGK distribution and collision operator
! Here, the RHS = nu * (f0 - f)
! where
! nu = collision frequency
! f0 = ESBGK distribution
! f = solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColESBGKSymW_1DMPI(f,RHS,n,u_0,v_0,w_0,T)

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_commvar, only: nodes_u_loc, nodes_v_loc, nodes_w_loc, &
				   alpha, gas_viscosity, gas_T_reference, gas_alpha, C_inf, gasR

real (DP), dimension (:), intent (in) :: f ! the components of the solution at the current time step. 
real (DP), dimension (:), intent (out) :: RHS ! the value of the collision operator for each component of the solution.

real (DP), intent (in) :: u_0, v_0, w_0 ! bulk velocities
real (DP), intent (in) :: n ! density
real (DP), intent (in) :: T ! temperature
!real (DP) :: Determinant ! the determinant of the Tensor
real (DP) :: Pressure
real (DP) :: nu ! this is the collision frequency term
real (DP), dimension (:), allocatable :: f0 ! Distribution function

real (DP), parameter :: kBoltzmann = 1.3806503D-23
real (DP), dimension (3,3) :: Tensor, TensorInv ! both the tensor and the inverse tensor for the ES-BGK
real (DP) :: Determinant ! the determinant of the tensor
integer :: loc_alloc_stat


! the tensor and its corresponding inverse and determinant is computed here
call MakeTensorSymW_1DMPI (f,Tensor,n,u_0,v_0,w_0,T) ! alpha only in common variabls
TensorInv = inv(Tensor)
Determinant = DetTensor(Tensor)

allocate (f0(1:size(f,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "EvalColESBGKSymW_1DMPI: Allocation error for variables (f0)"
 stop
end if

f0 = ESBGK_f0(TensorInv,Determinant,n,u_0,v_0,w_0,nodes_u_loc,nodes_v_loc,nodes_w_loc)

! now to evaluate the collision requency term
Pressure = n*T ! dimensionless Pressure is computed here
nu = Pressure/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/T)**gas_alpha ! final dimensionless nu ! have gas_T_reference be dimensionless?
RHS = nu * (f0 - f)

deallocate (f0)
!
end subroutine EvalColESBGKSymW_1DMPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is a modification of the  MakeTensor adapted to MPI parallelization in velocity 
! vvariable in 1D spatial problem.
!
! 
! Next, we compute the Tensor and also the Determinant and inverse of the tensor for use later
! in the ES-BGK distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MakeTensor_1DMPI (f,Tensor,n,ubar,vbar,wbar,Temp)
! Compute the tensor T = (1-alpha)*TRI + 1/n*< c (x) c >
!
!          u   v   w
!     u | T11 T12 T13 |
! T = v | T21 T22 T23 |
!     w | T31 T32 T33 |
!
! This subroutine produces the tensor that we see above associated with the ES-BGK distribution

use DGV_commvar, only: nodes_gwts_loc,nodes_u_loc,nodes_v_loc,nodes_w_loc,alpha,gasR,&
                       MPI_LINEAR_COLL_WORLD !need to introduce alpha in the read/write parameters.dat file

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (3,3), intent (out) :: Tensor ! both the tensor and the inverse tensor for the ES-BGK
real (DP), intent (in) :: ubar, vbar, wbar ! the bulk velocities
real (DP), intent (in) :: n, Temp ! the density and temperature 
!!!!!!!
! Note that n,ubar,vbar,wbar,Temp can be computed form f, but will require an additional MPI communication. Therefore, we 
! recycle the values from determining the mode of solution. 
!!!!!!!
!!! MPI related
real (DP), dimension (1:6) :: send_buff, recv_buff ! send and recveive buffers
integer :: loc_alloc_stat, ierr
!!!!! Evaluate the local portions of the tensor
send_buff(1) = sum(f*nodes_gwts_loc*(nodes_u_loc-ubar)**2)/n*alpha*2.0_DP
send_buff(2) = sum(f*nodes_gwts_loc*(nodes_v_loc-vbar)**2)/n*alpha*2.0_DP
send_buff(3) = sum(f*nodes_gwts_loc*(nodes_w_loc-wbar)**2)/n*alpha*2.0_DP
send_buff(4) = sum(f*nodes_gwts_loc*(nodes_u_loc-ubar)*(nodes_v_loc-vbar))/n*alpha*2.0_DP
send_buff(5) = sum(f*nodes_gwts_loc*(nodes_u_loc-ubar)*(nodes_w_loc-wbar))/n*alpha*2.0_DP
send_buff(6) = sum(f*nodes_gwts_loc*(nodes_v_loc-vbar)*(nodes_w_loc-wbar))/n*alpha*2.0_DP
recv_buff=0
!!!!!!!
call mpi_allreduce(send_buff,recv_buff,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"MakeTensorSymW_1DMPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
!!!!!!!
Tensor(1,1) = recv_buff(1)
Tensor(2,2) = recv_buff(2)
Tensor(3,3) = recv_buff(3)
Tensor(1,2) = recv_buff(4)
Tensor(1,3) = recv_buff(5)
Tensor(2,3) = recv_buff(6)
Tensor(2,1) = Tensor(1,2)
Tensor(3,1) = Tensor(1,3)
Tensor(3,2) = Tensor(2,3)
!!!!!!!
Tensor(1,1) = Tensor(1,1) + (1-alpha)*Temp
Tensor(2,2) = Tensor(2,2) + (1-alpha)*Temp
Tensor(3,3) = Tensor(3,3) + (1-alpha)*Temp
!!!!!!!
end subroutine MakeTensor_1DMPI




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MakeTensorSymW_1DMPI
!
! This subroutine is a modification of MakeTensor_1DMPI adapted for the case of symmetry in 
! W-component of velocity
!
! This subroutine is a modification of the  MakeTensor adapted to MPI parallelization in velocity 
! vvariable in 1D spatial problem.
!
! 
! Next, we compute the Tensor and also the Determinant and inverse of the tensor for use later
! in the ES-BGK distribution
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MakeTensorSymW_1DMPI (f,Tensor,n,ubar,vbar,wbar,Temp)
! Compute the tensor T = (1-alpha)*TRI + 1/n*< c (x) c >
!
!          u   v   w
!     u | T11 T12 T13 |
! T = v | T21 T22 T23 |
!     w | T31 T32 T33 |
!
! This subroutine produces the tensor that we see above associated with the ES-BGK distribution

use DGV_commvar, only: nodes_gwts_loc,nodes_u_loc,nodes_v_loc,nodes_w_loc,alpha,gasR,&
                       nodes_symfctr_loc,MPI_LINEAR_COLL_WORLD !need to introduce alpha in the read/write parameters.dat file

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


!!!!!!!
real (DP), dimension (:), intent (in) :: f ! the solution at the current  
real (DP), dimension (3,3), intent (out) :: Tensor ! both the tensor and the inverse tensor for the ES-BGK
real (DP), intent (in) :: ubar, vbar, wbar ! the bulk velocities
real (DP), intent (in) :: n, Temp ! the density and temperature 
!!!!!!!
! Note that n,ubar,vbar,wbar,Temp can be computed form f, but will require an additional MPI communication. Therefore, we 
! recycle the values from determining the mode of solution. 
!!!!!!!
!!! MPI related
real (DP), dimension (1:6) :: send_buff, recv_buff ! send and recveive buffers
integer :: loc_alloc_stat, ierr
!!!!! Evaluate the local portions of the tensor
send_buff(1) = sum(f*nodes_gwts_loc*nodes_symfctr_loc*(nodes_u_loc-ubar)**2)/n*alpha*2.0_DP
send_buff(2) = sum(f*nodes_gwts_loc*nodes_symfctr_loc*(nodes_v_loc-vbar)**2)/n*alpha*2.0_DP
send_buff(3) = sum(f*nodes_gwts_loc*nodes_symfctr_loc*(nodes_w_loc-wbar)**2)/n*alpha*2.0_DP
send_buff(4) = sum(f*nodes_gwts_loc*nodes_symfctr_loc*(nodes_u_loc-ubar)*(nodes_v_loc-vbar))/n*alpha*2.0_DP
send_buff(5) = 0 
send_buff(6) = 0 
recv_buff=0
!!!!!!!
call mpi_allreduce(send_buff,recv_buff,6,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
if (ierr /= 0 ) then 
 print *,"MakeTensor_1DMPI: Error computing macroparamters. MPI_REDUCE returned error", ierr
 stop
end if
!!!!!!!
Tensor(1,1) = recv_buff(1)
Tensor(2,2) = recv_buff(2)
Tensor(3,3) = recv_buff(3)
Tensor(1,2) = recv_buff(4)
Tensor(1,3) = 0
Tensor(2,3) = 0
Tensor(2,1) = Tensor(1,2)
Tensor(3,1) = Tensor(1,3)
Tensor(3,2) = Tensor(2,3)
!!!!!!!
Tensor(1,1) = Tensor(1,1) + (1-alpha)*Temp
Tensor(2,2) = Tensor(2,2) + (1-alpha)*Temp
Tensor(3,3) = Tensor(3,3) + (1-alpha)*Temp
!!!!!!!
end subroutine MakeTensorSymW_1DMPI




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColGrad_1Donecell_MPI(f,Df,fcol,time,cellnum)
!
! This is a modification of the subroutine EvalColGrad_1Donecel to prallelize it using MPI.
!
! This is a copy of the above subroutine adjusted for 1D spatial solvers
! 
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColGrad_1Donecell_MPI(f,Df,fcol,time,dens,ubar,vbar,wbar,temp,L1_err,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                     nodes_gwts_loc,nodes_u_loc,nodes_v_loc,nodes_w_loc,&
                     alpha, gas_alpha, gas_T_reference, C_inf, gasR,&
                     MPI_LINEAR_COLL_WORLD
                                        
use DGV_distributions_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated
real (DP), intent (in) :: dens,ubar,vbar,wbar,temp ! macroparamters of the solution computed elsewhere passed to this routine
real (DP), intent (in) :: L1_err ! L1 norm of the difference between the solution and the local maxwellian computed elsewhere passed to this routine
                    ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: sigma, sigma2 ! scrap variable to keep the traceless artificial stress tensor and its inverse 
real (DP), dimension(1:3) :: q !heat flux array
real (DP), dimension(Order) :: moments !Store the moments
real (DP) :: fij, dfij,trace_sigma,nu ! scrap variables
! 
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr, nu_fix ! scrap variable to check positivity of entropy.
real (DP), dimension(1:20) :: dp_buff_send,dp_buff_recv ! buiffer for MPI communications
integer :: ierr,irank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Evaluate ES-BGK nu - collision frequency
nu = dens*temp/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/Temp)**gas_alpha  ! Dimensionless nu -- the classical collision frequency of the ES-BGK model
! classican ES-BGK nu is used as a backup collision frequency 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRates1Donecell_DGV_MPI(momsrates,f,time,L1_err,dens,ubar,vbar,wbar,temp,nu,momsrates_reliab,cellnum) ! the subroutine returns momrates
                                           	! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! The moments are enforced by adjusting the parameters of Grad 13 moment distribution to enforce the 13 moments.
 ! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! BEGIN LATEST FIX
 ! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nu_fix = max(moments(2),moments(3),moments(4))
  if (nu_fix > 0) then 
   nu = nu_fix 
  end if
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! END LATEST FIX
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i = 1,10
  !!!!! compute the portion of the moments that corresponds to the portion of the solution on this processor
  fij=sum(f*nodes_gwts_loc*kernls_enfrsd_moms(i,nodes_u_loc,nodes_v_loc,nodes_w_loc,ubar,vbar,wbar))
  dfij=sum(Df*nodes_gwts_loc*kernls_enfrsd_moms(i,nodes_u_loc,nodes_v_loc,nodes_w_loc,ubar,vbar,wbar))
  !!! now we need to combine this over all processors sharing the solution: 
  dp_buff_send(2*i-1)=fij
  dp_buff_send(2*i)=dfij
 end do  
 dp_buff_recv=0
 !!!!!!!!!!!!!!
 call mpi_allreduce(dp_buff_send,dp_buff_recv,20,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
 if (ierr /= 0 ) then 
  print *," EvalColGrad_1Donecell_MPI: Error gathering moment on prim. mesh. MPI_ALLREDUCE returned error", ierr
  stop
 end if
  moments(1) = dp_buff_recv(1)
 do i = 2,10
 !!!!!!!!!!!!!!
  fij=dp_buff_recv(2*i-1)
  dfij=dp_buff_recv(2*i)
  moments(i) = (fij - (momsrates(i)/(nu))*dfij)
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
 
 fcol = nu * (EvalGrad13f0(dens,moments(1),sigma,q,ubar,vbar,wbar,nodes_u_loc,nodes_v_loc,nodes_w_loc)  - f )
else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
 if (irank==0) then 
  PRINT *, "EvalColGrad_1Donecell_MPI:  Invoke fall back model. momsrates_reliab, cellnum", &
                                 momsrates_reliab, cellnum
 end if 
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGK_1DMPI(f,fcol,dens,ubar,vbar,wbar,temp)
end if
! all done
end subroutine EvalColGrad_1Donecell_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalColGradSymW_1Donecell_MPI(f,Df,fcol,time,cellnum)
!
! This is a modification of the subroutine EvalColGrad_1Donecell_MPI adjusted for symmetry in 
! w-component in velocity 
!
! This is a modification of the subroutine EvalColGrad_1Donecel to prallelize it using MPI.
!
! This is a copy of the above subroutine adjusted for 1D spatial solvers
! 
! This subroutine evaluates the collision operator using the model in which 
! relaxation rates are enforced for moments (u-\bar{u})(u-\bar{u})^T 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine EvalColGradSymW_1Donecell_MPI(f,Df,fcol,time,dens,ubar,vbar,wbar,temp,L1_err,cellnum)

use DGV_dgvtools_mod

use DGV_commvar, only: MaxNumEnforcedMoments,MaxNumBFunctionsCollFreq,Order_nu,Order,&
                     nodes_gwts_loc,nodes_u_loc,nodes_v_loc,nodes_w_loc,nodes_symfctr_loc,&
                     alpha, gas_alpha, gas_T_reference, C_inf, gasR,&
                     MPI_LINEAR_COLL_WORLD
                                        
use DGV_distributions_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


intrinsic MAXVAL, ABS

real (DP), dimension (:), intent (in) :: f,Df ! the components of the solution at the current time step. 
						!Df = f - fM, Df is the difference between the solution and the local Maxwellian
real (DP), dimension (:), intent (out) :: fcol ! the value of the collision operator for each component of the solution.
real (DP), intent (in) :: time ! the current time
integer (I4B) :: cellnum ! number of the spatial cell for which the collision operator is evaluated
real (DP), intent (in) :: dens,ubar,vbar,wbar,temp ! macroparamters of the solution computed elsewhere passed to this routine
real (DP), intent (in) :: L1_err ! L1 norm of the difference between the solution and the local maxwellian computed elsewhere passed to this routine
                    ! Relative L1 norm of the deviation of the solution from the local Maxwellian.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
real (DP), dimension(1:3,1:3) :: sigma, sigma2 ! scrap variable to keep the traceless artificial stress tensor and its inverse 
real (DP), dimension(1:3) :: q !heat flux array
real (DP), dimension(Order) :: moments !Store the moments
real (DP) :: fij, dfij,trace_sigma,nu ! scrap variables
! 
!!!!
logical :: updateNulcl ! a scrap logical variable to use in updating relaxation rates.  
real (DP), dimension (Order) :: momsrates ! a local array to store the relaxation rates for the moments 
real (DP), dimension (Order_nu) :: Cco_temp ! local temporaty variable to store oefficeints for the vel. dep. collision requency
real (DP), dimension (Order,Order_nu) :: Mom, MomInv ! the matrix of the system that is solved to determine the coefficients
real (DP), dimension (Order) :: DifMom, Dphi ! scrap array to contain the values of the differenced between the enforced moments and their eq1uilibrium values. 
real (DP), dimension (:), allocatable :: nuBolt ! the values of the velocity dependent collision frequency 
!!!!
integer (I4B) :: i ! scrap index 
logical :: momsrates_reliab ! scrap variable that keeps the flag of reliability of   moment relaxation rates. If true then moments rates are reliable and can be used 
                   !for the evaluation of velocity dependent collision frequency 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostics variables... 
real (DP) :: min_nu, max_nu ! scrap variable to check if collision frequency goes below zero.
real (DP) :: temp_dp, entr ! scrap variable to check positivity of entropy.
real (DP), dimension(1:20) :: dp_buff_send,dp_buff_recv ! buiffer for MPI communications
integer :: ierr,irank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Evaluate ES-BGK nu - collision frequency
nu = dens*temp/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/Temp)**gas_alpha  ! Dimensionless nu -- the classical collision frequency of the ES-BGK model
! classican ES-BGK nu is used as a backup collision frequency 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Then, we call a conglomerate subroutine that will determine the values of the enforced
! relaxation rates. The subroutine will check if the rates need to be updated. If the rates 
! are not updated, the stored rates are returned.
! if the rates need to be updated, the Boltzmann collision integral is evaluated and 
! the rates are determined from the Boltzmann collision integral
!!!!!!!!!!!!!!!!!!!!!!!!
call GetRelaxRatesSymW1Donecell_DGV_MPI(momsrates,f,time,L1_err,dens,ubar,vbar,wbar,temp,nu,momsrates_reliab,cellnum) ! the subroutine returns momrates
                                           	! and also the value of the L1_norm of the difference between the soltuion and the local maxwellian
											! and also returns a flag momsrates_reliab. If this flag is true then at least one rate was computed 
											! from the Boltzmann collision operator. Otherwise, coded default relaxation rate was used. This rate is 
											! returned in the variable nu

if (momsrates_reliab) then 
 !!!!!!!!!!!!!!!!!!!!!!!
 ! The moments are enforced by adjusting the parameters of Grad 13 moment distribution to enforce the 13 moments.
 ! 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 do i = 1,10
  !!!!! compute the portion of the moments that corresponds to the portion of the solution on this processor
  fij=sum(f*nodes_gwts_loc*nodes_symfctr_loc*kernls_enfrsd_moms(i,nodes_u_loc,nodes_v_loc,nodes_w_loc,ubar,vbar,wbar))
  dfij=sum(Df*nodes_gwts_loc*nodes_symfctr_loc*kernls_enfrsd_moms(i,nodes_u_loc,nodes_v_loc,nodes_w_loc,ubar,vbar,wbar))
  !!! now we need to combine this over all processors sharing the solution: 
  dp_buff_send(2*i-1)=fij
  dp_buff_send(2*i)=dfij
 end do  
 dp_buff_recv=0
 !!!!!!!!!!!!!!
 call mpi_allreduce(dp_buff_send,dp_buff_recv,20,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
 if (ierr /= 0 ) then 
  print *," EvalColGrad_1Donecell_MPI: Error gathering moment on prim. mesh. MPI_ALLREDUCE returned error", ierr
  stop
 end if
 moments(1) = dp_buff_recv(1)
 do i = 2,10
 !!!!!!!!!!!!!!
  fij=dp_buff_recv(2*i-1)
  dfij=dp_buff_recv(2*i)
  moments(i) = (fij - (momsrates(i)/(nu))*dfij)
 end do
 !!!!!!!!!!!!!!
 !moments(1) is pressure
 sigma(1,1)=moments(2)
 sigma(2,2)=moments(3)
 sigma(3,3)=moments(4)
 sigma(1,2)=moments(5)
 sigma(2,1)=moments(5)
 sigma(1,3)=0  ! value is zero because of the symmetry in W-component
 sigma(3,1)=0  ! value is zero because of the symmetry in W-component
 sigma(2,3)=0  ! value is zero because of the symmetry in W-component
 sigma(3,2)=0  ! value is zero because of the symmetry in W-component
 q(1)=moments(8)
 q(2)=moments(9)
 q(3)=0 ! value is zero because of the symmetry in W-component
 
 fcol = nu * (EvalGrad13f0(dens,moments(1),sigma,q,ubar,vbar,wbar,nodes_u_loc,nodes_v_loc,nodes_w_loc)  - f )

else 
 ! If the moments rates are not reliable or if the Mom matrix is very small by L1-norm, we use the fallback model instead:    
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
 if (irank==0) then 
  PRINT *, "EvalColGradSymW_1Donecell_MPI:  Invoke fall back model. momsrates_reliab, cellnum", &
                                 momsrates_reliab, cellnum
 end if 
 ! call ES-BGK model since the model with velocity-dependent collision frequency is expected to fail. 
 call EvalColESBGKSymW_1DMPI(f,fcol,dens,ubar,vbar,wbar,temp)
end if
! all done
end subroutine EvalColGradSymW_1Donecell_MPI



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRelaxRates1Donecell_DGV_MPI(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab,cellnum)
!
! This is a modification of GetRelaxRates1Donecell_DGV adjusted to run in MPI parallel version
!
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

subroutine GetRelaxRates1Donecell_DGV_MPI(momsrates,f,time,L1_err,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,nu,&
                                 momsrates_reliab,cellnum)

use DGV_commvar, only: MaxNumEnforcedMoments, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       MoRatesArry1D, MoRatesReliable1D,&
                       lin_proc_nodes_wld, MPI_LINEAR_COLL_WORLD,num_lin_proc,&
                       NuNextUpdateTime1D

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_collision_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension(:), intent(in) :: f ! the solution(velocity distribution) on one spatial cell
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent (in) :: time ! value of the dimensionless time variable
real (DP), intent (in) :: L1_err ! value of the relative L_1 norm of the differnce between the solution and the local Maxwellian
real (DP), intent (in) :: LocUbar,LocVbar,LocWbar ! the values of the local bulk velocity
real (DP), intent (in) :: LocDens ! value of the local numerical density is returned 
real (DP), intent (in) :: LocTempr ! value of the local temperature is returned                              
real (DP), intent (in) :: nu ! value of the default relaxation rate -- will be assigned to moments for which the ralaxation rates can not be 
                              ! calculated form the Boltzmann collision operator (for this cell)
logical, intent (out) :: momsrates_reliab ! is true if at least one rate has calculated from the Boltzmann collision operator. (for this cell)
integer (I4B), intent (in) :: cellnum ! the number of the spatial cell for which the rates are evaluated
integer :: irank,mpicommworldsize ! rank of the processor that is calling this subroutine,size of the MPI universe 
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
real (DP), dimension(:), allocatable :: f_full,f_full_send ! temp variables to store the entire solution on the primary velocity mesh
real (DP), dimension(:), allocatable :: fII,fcolII,fcol1II,fMII ! temp variables to store the 
 ! value of the collision integral on secondary and the local Maxwellian on the secondary mesh and the 
 ! differnce between the local maxwellian and the local maxwellian on the secondary mesh.
 ! and the difference between the local maxwellian and the solution on the promary mesh
logical :: updateNulcl ! scrap variable to pass the update flag
integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
integer (I4B) :: i,j,mm !scrap counters
integer :: ierr ! scrap to keep result of MPI operation
integer :: mpi_dgv_action ! Message that tells MPI algorithms what to do
integer, dimension(1:1) :: ibuff 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check that the nucoeffs has the right size:
if ((size(momsrates,1) > MaxNumEnforcedMoments) .or. (size(momsrates,1) /= size(MoRatesArry1D,1))) then
 print *, "GetRelaxRates1Donecell_DGV_MPI: Error. Size of the supplied array (momsrates) is incopatible."
 stop ! terminate the program sinse nucoefs has a wrong size
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we check if relaxation speeds in this cell need to be updated.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call CheckNuUpdateNeededF2Donecella_DGV(updateNulcl,time,LocDens,LocTempr,cellnum) ! the subroutine accesses all other variables directly from commvar
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
 allocate (fII(1:size(nodes_uII,1)),fcolII(1:size(nodes_uII,1)),f_full(1:size(nodes_u,1)), &
           f_full_send(1:size(nodes_u,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (fII),(fcolII),(f_full)"
      end if 
 ! First, we need a local copy of the entire solution on the Primary velocity mesh (we only have portion)
 f_full_send=0
 f_full=0
 ! copy the local solution to the send buffer:
 ! prepare copies of the velocity nodes to be used on this processor
 do i=1,size(lin_proc_nodes_wld,1)
  j=lin_proc_nodes_wld(i)
  f_full_send(j)=f(i)
 end do 
 mm = size(nodes_u,1)
 call mpi_allreduce(f_full_send,f_full,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
 if (ierr /= 0 ) then 
  print *,"GetRelaxRates1Donecell_DGV_MPI: Error gathering solution on prim. mesh. MPI_REDUCE returned error", ierr
  stop
 end if
  ! To project the solution to a (most likely smaller) secondary mesh, we use the following subroutine: 
 call ProjectDGVtoSecMesh(f_full,fII) ! put the projecteion in fII
 deallocate(f_full_send)
 !!!!!!! 
 ! Here we address a problem that algorythmical. If this version of the code, only processors with irank<nul_lin_proc 
 ! will be calling the collision operator. The original idea was that this number may be less then the entire 
 ! size of the MPI world comminicator. This leads to some unwanted concequnces, in particular evaluation of the 
 ! is designed so that collision operator i evaluated using all processors, but the one with the rank 0. So,
 ! we need to make sure that processors with irank >= num_lin_proc receive the solution on the secondary mesh, 
 ! that they need to evaluate the collision operator. So, we will have to send the solution over to them. This is taken care of next.
 ! Now, even a bigger problem is that the processors with irank >= num_lin_procs need to be waiting for a command to evaluate their portion 
 ! of the collision operator. This is done by using a special "fly trap" subroutine that makes those processors wait 
 ! for a message to arrive that will tell them what to do.  
 !!!!!!!
 !! Check the size of the MPI Universs = the total number of processors in the Universe.
 call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
 if (mpicommworldsize > num_lin_proc) then  !only tdo this send if there exist processors with irank>num_lin_procs-1 
  ! it is required that processors with irank>num_lin_procs-1 wait for this broadcast
  if (L1_err > L1_MAX) then 
   mpi_dgv_action = 401 ! integer2, request evaluation of the collision operator on the secondary mesh using Boltzmann 
   ibuff(1)=401
  else 
   mpi_dgv_action = 402 ! integer2, request evaluation of the collision operator on the secondary mesh using full mixed model
   ibuff(1)=402
  end if 
  call mpi_bcast (ibuff,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
   stop
  end if 
  ! upon receiving action code 401 or 402 all processor with irank > num_lin_procs will proceed with the following steps: 
  ! 1 -- receive the solution via broadcast from node 0, 2 --- call evaluation of the collision operator. 3-- gather the 
  ! collision operator on processor with rank 0
  select case (mpi_dgv_action)
   case (401) ! need to send solution
    mm = size(nodes_uII,1)
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of solution fII from proc 0 returned error", ierr
     stop
    end if 
   case (402) ! need to send solution and the macroparameters 
    mm = size(nodes_uII,1)
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "(2) MPI broadcast of solution fII from proc 0 returned error", ierr
     stop
    end if
    allocate(f_full_send(1:5), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
     print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (f_full_send)"
    end if 
    f_full_send(1)=LocDens
    f_full_send(2)=LocUbar
    f_full_send(3)=LocVbar
    f_full_send(4)=LocWbar
    f_full_send(5)=LocTempr
    call mpi_bcast (f_full_send,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
     stop
    end if
    deallocate(f_full_send)
   case default 
    print *, "GetRelaxRates1Donecell_DGV_MPI: unrecognized action code. Stop. irank= ", irank
    stop
  end select
  ! all info is send to the idling processors (irank > num_lin_proc) 
 end if  
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will check if the solution is needs to be evaluated in a decomposition mode or a full mode. 
 ! in the full mode $\int_{R^3}\int_{R^3} ffA is evaluated
 ! in the decomposition mode, the solution is split into a sum of the local Maxwellian and another functions
 ! $f=f_{M} + f_{0}$ (note that $f_{0}$ does not have to be small). and the integral 
 ! $\int_{R^3}\int_{R^3} (2f_{M}f_{0}A + f_{0}f_{0}A ) is evaluated instead
 ! this approach reduces errors when the solution is close to a Mawellian
 !!!
 !!! 
 ! because we expect that the solution on the primary mesh is mote accurate, we will check the 
 ! closeness of the solution to the Maxwellian using the primary mesh:
 ! We will use the byproduct of another other check, L1_err -- which is the relative L_1 norm of 
 ! the difference between the distribution f and the local maswellian on the primary mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! The following evaluation only make sense on processors with irank>0
 !!!
 fcolII=0 !reset the resuls, just in case
 if (irank>0) then ! there are no components of A stored on irank=0 processor 
  if (L1_err > L1_MAX) then 
   ! evaluate the Boltzmann collision operator using full mode
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Now we will call the subroutine that evaluates the Boltzmann collision operator on the secondary mesh. 
   ! There are two possibilites: an single processor call of the subroutine with possible OpenMP fork and 
   ! call or an MPI parallelization.  
   !!!!!!
   ! to make an MPI evaluation, uncomment the next lime
   call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
  else 
   ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
   allocate (fcol1II(1:size(nodes_uII,1)), fMII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
       if (loc_alloc_stat >0) then 
        print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (fMII),(fcol1II)"
       end if 
   fMII = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_uII,nodes_vII,nodes_wII) ! now we populate the maxwellian with the same macroparamters.
   fII = fII-fMII ! now fII contains the difference between the solution and the local maxwellian
   ! Now we will make to calls of the collision operator
   !!!!!!!!!!
   ! to make an MPI evaluation, uncomment the next two lines
   call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator f_{0}f_{0}A
   call EvalCollisionPeriodicMixedTermsA_DGVII_MPI(fII,fMII,fcol1II) !
   fcolII=fcolII+fcol1II
   deallocate (fcol1II,fMII)
  end if 
 end if 
 !!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!! The MPI processors have computed portions of the collision operator assigned to them. The portions of the collision operator that are not assigned to them
 !!! should be zero -- this is done insive EvalCollisionPeriodicAPlus_DGVII_MPI and EvalCollisionPeriodicMixedTermsA_DGVII_MPI
 !!! next we need to combine the components. We will do it using MPI collective communication 
 fII=0 ! nullify the array -- we will use it as a recieve buffer
 mm = size(nodes_uII,1)
 if (mpicommworldsize > num_lin_proc) then  !only do this send if there exist processors with irank>num_lin_procs-1 
  call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
   stop
  end if 
 else 
  call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) ! All processors will participate. However, no need to match in the "fly trap" subroutine
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI: (1) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
   stop
  end if
 end if
 ! Note that now, fII has the value of the collision operator!  
 !!!!!!!!!!!!!!!!! done combining  collision operator
 deallocate(fcolII)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! Next we compute the relaxation rates for the selected group of moments from the (Df), and (fcolII).
 ! (Df) on the primary mesh is taken to provide a better computation of the macroparameters.  
 ! fcolII is only available on the secondary mesh  
 !!!!!!!!!!!!!!!!!
 ! call ComputeRelaxRatesBCI_DGV(fcolII,f,momsrates,L1_SENS_TRSHLD)    ! Both subroutines do the same, but the second one uses qunatities compluted previously 
 f_full = f_full - maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u,nodes_v,nodes_w) ! f_full now is Df=f-fm
 call ComputeRelaxRatesBCIa_DGV(fII,f_full,momsrates,L1_SENS_TRSHLD,Mom_Close_Maxwl_TRSHLD,LocDens,LocUbar,LocVbar, &
        LocWbar,LocTempr,nu,momsrates_reliab) ! The first subroutine computes quantities directly
 MoRatesArry1D(:,cellnum) = momsrates ! safe the coefficients for future use.
 MoRatesReliable1D(cellnum) = momsrates_reliab ! update the flag is the relaxation rates were computed form the BCI 
 !!!!!!!!!!!!!!!!!
 deallocate(fII,f_full)
 !!!!!!!!!!!! DEBUG
 if (irank==0) then  
  print *, "cell=",cellnum,"momsrates",momsrates,"reliab=",momsrates_reliab, "next time=", NuNextUpdateTime1D(cellnum)
 end if 
 !!!!!!!!!!!! END Debug 
else 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! if the no update is needed, the rates are taken from the storage array
 momsrates = MoRatesArry1D(:,cellnum)
 momsrates_reliab = MoRatesReliable1D(cellnum) ! the flag = true means that rates were computed from the BCI
 !!!! 
end if
end subroutine GetRelaxRates1Donecell_DGV_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GetRelaxRatesSymW1Donecell_DGV_MPI(momsrates,f,Df,time,L1_err,LocUbar,LocVbar,LocWbar,nu,momsrates_reliab,cellnum)
!
! This is a modification of GetRelaxRates1Donecell_DGV_MPI to adjust for symmetry in W-component of velocity
!
!
! This is a modification of GetRelaxRates1Donecell_DGV adjusted to run in MPI parallel version
!
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

subroutine GetRelaxRatesSymW1Donecell_DGV_MPI(momsrates,f,time,L1_err,LocDens,LocUbar,LocVbar,LocWbar,LocTempr,nu,&
                                 momsrates_reliab,cellnum)

use DGV_commvar, only: MaxNumEnforcedMoments, nodes_uII, nodes_vII, nodes_wII, &
                       nodes_u, nodes_v, nodes_w, nodes_gwts,&
                       Mu,Mv,Mw,su,sv,sw,&
                       MoRatesArry1D, MoRatesReliable1D,&
                       lin_proc_nodes_wld, MPI_LINEAR_COLL_WORLD,num_lin_proc

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_collision_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension(:), intent(in) :: f ! the solution(velocity distribution) on one spatial cell
real (DP), dimension(:), intent(out):: momsrates ! this is the array which will contain the enforced relaxation rates
real (DP), intent (in) :: time ! value of the dimensionless time variable
real (DP), intent (in) :: L1_err ! value of the relative L_1 norm of the differnce between the solution and the local Maxwellian
real (DP), intent (in) :: LocUbar,LocVbar,LocWbar ! the values of the local bulk velocity
real (DP), intent (in) :: LocDens ! value of the local numerical density is returned 
real (DP), intent (in) :: LocTempr ! value of the local temperature is returned                              
real (DP), intent (in) :: nu ! value of the default relaxation rate -- will be assigned to moments for which the ralaxation rates can not be 
                              ! calculated form the Boltzmann collision operator (for this cell)
logical, intent (out) :: momsrates_reliab ! is true if at least one rate has calculated from the Boltzmann collision operator. (for this cell)
integer (I4B), intent (in) :: cellnum ! the number of the spatial cell for which the rates are evaluated
integer :: irank,mpicommworldsize ! rank of the processor that is calling this subroutine,size of the MPI universe 
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
real (DP), dimension(:), allocatable :: f_full,f_full_send ! temp variables to store the entire solution on the primary velocity mesh
real (DP), dimension(:), allocatable :: fII,fcolII,fcol1II,fMII ! temp variables to store the 
 ! value of the collision integral on secondary and the local Maxwellian on the secondary mesh and the 
 ! differnce between the local maxwellian and the local maxwellian on the secondary mesh.
 ! and the difference between the local maxwellian and the solution on the promary mesh
logical :: updateNulcl ! scrap variable to pass the update flag
integer :: loc_alloc_stat ! scrap variable to kee the allocation status. 
integer (I4B) :: i,j,mm !scrap counters
integer :: ierr ! scrap to keep result of MPI operation
integer :: mpi_dgv_action ! Message that tells MPI algorithms what to do
integer, dimension(1:1) :: ibuff 
integer (I4B) :: iu,iv,iw,iiu,iiv,iiw,iwsym,iiwsym,isym ! scrap indices/
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! a quick check that the nucoeffs has the right size:
if ((size(momsrates,1) > MaxNumEnforcedMoments) .or. (size(momsrates,1) /= size(MoRatesArry1D,1))) then
 print *, "GetRelaxRates1Donecell_DGV_MPI: Error. Size of the supplied array (momsrates) is incopatible."
 stop ! terminate the program sinse nucoefs has a wrong size
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we check if relaxation speeds in this cell need to be updated.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call CheckNuUpdateNeededF2Donecella_DGV(updateNulcl,time,LocDens,LocTempr,cellnum) ! the subroutine accesses all other variables directly from commvar
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
 allocate (fII(1:size(nodes_uII,1)),fcolII(1:size(nodes_uII,1)),f_full(1:size(nodes_u,1)), &
           f_full_send(1:size(nodes_u,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (fII),(fcolII),(f_full)"
      end if 
 ! First, we need a local copy of the entire solution on the Primary velocity mesh (we only have portion)
 f_full_send=0
 f_full=0
 ! copy the local solution to the send buffer:
 ! prepare copies of the velocity nodes to be used on this processor
 do i=1,size(lin_proc_nodes_wld,1)
  j=lin_proc_nodes_wld(i)
  f_full_send(j)=f(i)
 end do 
 mm = size(nodes_u,1)
 call mpi_allreduce(f_full_send,f_full,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) !
 if (ierr /= 0 ) then 
  print *,"GetRelaxRates1Donecell_DGV_MPI: Error gathering solution on prim. mesh. MPI_REDUCE returned error", ierr
  stop
 end if
 !!!! Because of thed symmetry, only half of the solution is non-zero. We now will use symmetry to repopulate the 
 !!!! Missing values. 
 !!!! Attention -- the next block assumes that velocity nodes are symmetric around the origin in W-direction 
 !!!! Algorithm will not work if this assumption fails.
 do iu=1,Mu
  do iiu=1,su
   do iv=1,Mv
    do iiv=1,sv
     do iw=1,Mw
      do iiw=1,sw
       i= ((iu-1)*Mv*Mw+(iv-1)*Mw+iw-1)*su*sv*sw + ((iiu-1)*sv*sw+(iiv-1)*sw+iiw)
       if (nodes_w(i) < 0.0_DP) then 
        iwsym = 1 + Mw - iw
        iiwsym = 1 + sw - iiw
        isym = ((iu-1)*Mv*Mw+(iv-1)*Mw+iwsym-1)*su*sv*sw + ((iiu-1)*sv*sw+(iiv-1)*sw+iiwsym)
        f_full(isym)=f_full(i)  ! the solution array was symmetrized
       end if 
      end do 
     end do 
    end do     
   end do  
  end do  
 end do
 !!!! End of the block 
 ! To project the solution to a (most likely smaller) secondary mesh, we use the following subroutine: 
 call ProjectDGVtoSecMesh(f_full,fII) ! put the projecteion in fII
 deallocate(f_full_send)
 !!!!!!! 
 ! Here we address a problem that algorythmical. If this version of the code, only processors with irank<nul_lin_proc 
 ! will be calling the collision operator. The original idea was that this number may be less then the entire 
 ! size of the MPI world comminicator. This leads to some unwanted concequnces, in particular evaluation of the 
 ! is designed so that collision operator i evaluated using all processors, but the one with the rank 0. So,
 ! we need to make sure that processors with irank >= num_lin_proc receive the solution on the secondary mesh, 
 ! that they need to evaluate the collision operator. So, we will have to send the solution over to them. This is taken care of next.
 ! Now, even a bigger problem is that the processors with irank >= num_lin_procs need to be waiting for a command to evaluate their portion 
 ! of the collision operator. This is done by using a special "fly trap" subroutine that makes those processors wait 
 ! for a message to arrive that will tell them what to do.  
 !!!!!!!
 !! Check the size of the MPI Universs = the total number of processors in the Universe.
 call MPI_COMM_SIZE(MPI_COMM_WORLD,mpicommworldsize,ierr) ! check how many processors is in the main communicator
 call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
 if (mpicommworldsize > num_lin_proc) then  !only tdo this send if there exist processors with irank>num_lin_procs-1 
  ! it is required that processors with irank>num_lin_procs-1 wait for this broadcast
  if (L1_err > L1_MAX) then 
   mpi_dgv_action = 401 ! integer2, request evaluation of the collision operator on the secondary mesh using Boltzmann 
   ibuff(1)=401
  else 
   mpi_dgv_action = 402 ! integer2, request evaluation of the collision operator on the secondary mesh using full mixed model
   ibuff(1)=402
  end if 
  call mpi_bcast (ibuff,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
   stop
  end if 
  ! upon receiving action code 401 or 402 all processor with irank > num_lin_procs will proceed with the following steps: 
  ! 1 -- receive the solution via broadcast from node 0, 2 --- call evaluation of the collision operator. 3-- gather the 
  ! collision operator on processor with rank 0
  select case (mpi_dgv_action)
   case (401) ! need to send solution
    mm = size(nodes_uII,1)
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of solution fII from proc 0 returned error", ierr
     stop
    end if 
   case (402) ! need to send solution and the macroparameters 
    mm = size(nodes_uII,1)
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "(2) MPI broadcast of solution fII from proc 0 returned error", ierr
     stop
    end if
    allocate(f_full_send(1:5), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
     print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (f_full_send)"
    end if 
    f_full_send(1)=LocDens
    f_full_send(2)=LocUbar
    f_full_send(3)=LocVbar
    f_full_send(4)=LocWbar
    f_full_send(5)=LocTempr
    call mpi_bcast (f_full_send,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"GetRelaxRates1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
     stop
    end if
    deallocate(f_full_send)
   case default 
    print *, "GetRelaxRates1Donecell_DGV_MPI: unrecognized action code. Stop. irank= ", irank
    stop
  end select
  ! all info is send to the idling processors (irank > num_lin_proc) 
 end if  
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now we will check if the solution is needs to be evaluated in a decomposition mode or a full mode. 
 ! in the full mode $\int_{R^3}\int_{R^3} ffA is evaluated
 ! in the decomposition mode, the solution is split into a sum of the local Maxwellian and another functions
 ! $f=f_{M} + f_{0}$ (note that $f_{0}$ does not have to be small). and the integral 
 ! $\int_{R^3}\int_{R^3} (2f_{M}f_{0}A + f_{0}f_{0}A ) is evaluated instead
 ! this approach reduces errors when the solution is close to a Mawellian
 !!!
 !!! 
 ! because we expect that the solution on the primary mesh is mote accurate, we will check the 
 ! closeness of the solution to the Maxwellian using the primary mesh:
 ! We will use the byproduct of another other check, L1_err -- which is the relative L_1 norm of 
 ! the difference between the distribution f and the local maswellian on the primary mesh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!! The following evaluation only make sense on processors with irank>0
 !!!
 fcolII=0 !reset the resuls, just in case
 if (irank>0) then ! there are no components of A stored on irank=0 processor 
  if (L1_err > L1_MAX) then 
   ! evaluate the Boltzmann collision operator using full mode
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Now we will call the subroutine that evaluates the Boltzmann collision operator on the secondary mesh. 
   ! There are two possibilites: an single processor call of the subroutine with possible OpenMP fork and 
   ! call or an MPI parallelization.  
   !!!!!!
   ! to make an MPI evaluation, uncomment the next lime
   call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!
  else 
   ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
   allocate (fcol1II(1:size(nodes_uII,1)), fMII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
       if (loc_alloc_stat >0) then 
        print *, "GetRelaxRates1Donecell_DGV_MPI: Allocation error for variables (fMII),(fcol1II)"
       end if 
   fMII = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_uII,nodes_vII,nodes_wII) ! now we populate the maxwellian with the same macroparamters.
   fII = fII-fMII ! now fII contains the difference between the solution and the local maxwellian
   ! Now we will make to calls of the collision operator
   !!!!!!!!!!
   ! to make an MPI evaluation, uncomment the next two lines
   call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator f_{0}f_{0}A
   call EvalCollisionPeriodicMixedTermsA_DGVII_MPI(fII,fMII,fcol1II) !
   fcolII=fcolII+fcol1II
   deallocate (fcol1II,fMII)
  end if 
 end if 
 !!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!! The MPI processors have computed portions of the collision operator assigned to them. The portions of the collision operator that are not assigned to them
 !!! should be zero -- this is done insive EvalCollisionPeriodicAPlus_DGVII_MPI and EvalCollisionPeriodicMixedTermsA_DGVII_MPI
 !!! next we need to combine the components. We will do it using MPI collective communication 
 fII=0 ! nullify the array -- we will use it as a recieve buffer
 mm = size(nodes_uII,1)
 if (mpicommworldsize > num_lin_proc) then  !only do this send if there exist processors with irank>num_lin_procs-1 
  call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
   stop
  end if 
 else 
  call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_LINEAR_COLL_WORLD,ierr) ! All processors will participate. However, no need to match in the "fly trap" subroutine
  if (ierr /= 0 ) then 
   print *,"GetRelaxRates1Donecell_DGV_MPI: (1) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
   stop
  end if
 end if
 ! Note that now, fII has the value of the collision operator!  
 !!!!!!!!!!!!!!!!! done combining  collision operator
 deallocate(fcolII)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 ! Next we compute the relaxation rates for the selected group of moments from the (Df), and (fcolII).
 ! (Df) on the primary mesh is taken to provide a better computation of the macroparameters.  
 ! fcolII is only available on the secondary mesh  
 !!!!!!!!!!!!!!!!!
 ! call ComputeRelaxRatesBCI_DGV(fcolII,f,momsrates,L1_SENS_TRSHLD)    ! Both subroutines do the same, but the second one uses qunatities compluted previously 
 f_full = f_full - maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_u,nodes_v,nodes_w) ! f_full now is Df=f-fm
 call ComputeRelaxRatesBCIa_DGV(fII,f_full,momsrates,L1_SENS_TRSHLD,Mom_Close_Maxwl_TRSHLD,LocDens,LocUbar,LocVbar, &
        LocWbar,LocTempr,nu,momsrates_reliab) ! The first subroutine computes quantities directly
 MoRatesArry1D(:,cellnum) = momsrates ! safe the coefficients for future use.
 MoRatesReliable1D(cellnum) = momsrates_reliab ! update the flag is the relaxation rates were computed form the BCI 
 !!!!!!!!!!!!!!!!!
 deallocate(fII,f_full)
 !!!!!!!!!!!! DEBUG
 if (irank==0) then  
  print *, "cell=",cellnum,"momsrates",momsrates,"reliab=",momsrates_reliab
 end if 
 !!!!!!!!!!!! END Debug 
else 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! if the no update is needed, the rates are taken from the storage array
 momsrates = MoRatesArry1D(:,cellnum)
 momsrates_reliab = MoRatesReliable1D(cellnum) ! the flag = true means that rates were computed from the BCI
 !!!! 
end if
end subroutine GetRelaxRatesSymW1Donecell_DGV_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FlyTrap1Donecell_DGV_MPI 
!
! This subroutine will be executed on idling processors. It contains pieces of MPI parallel algorithms that need to be 
! executed on processors with irank>num_lin_proc
!
! The work of this subroutine is determined by action codes. The subroutine is an infinite loop, 
! The first operation in the loop is a broadcast from irank=0. This broadcast will send an action code that 
! determines what the subroutine does next. The subroutine will perform action that are specified by the give action code and will return to the 
! expectation of broadcast.
!
! action codes that are currently implemented:
! 
! -777 exit code
! 401 -- evaluation of collision operator on the secondary mesh using full Boltzmann
! 402 -- evaluation of the collision operator on the seocondary mesh using decomposition f=fm+df
!
!!!!!!!!!!!!!
subroutine FlyTrap1Donecell_DGV_MPI

use DGV_commvar, only: nodes_uII,nodes_vII,nodes_wII
use DGV_distributions_mod
use DGV_collision_mod


!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: mm ! scrap variables
real (DP), dimension(:),allocatable :: fII,fcolII,dpbuff,fcol1II,fmII
integer :: loc_alloc_stat
real (DP) :: LocDens,LocUbar,LocVbar,LocWbar,LocTempr !variables to keep macroparameters

!!!!!!!!!!!!!!!!!!!!!!!
integer :: irank,ierr          ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
!!!!!!!!!!!!!!!!!!!!!!
call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...

mpi_dgv_action=0
do while (mpi_dgv_action /= -777)
 ibuff=0  
 call mpi_bcast (ibuff,1,MPI_INTEGER2,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"FlyTrap1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
  stop
 end if
 mpi_dgv_action=ibuff(1) 
 !!!!!!!!!!!!!!!!! MAIN FORK
 select case (mpi_dgv_action)
  case (401) ! this action code corresponds to evaluation of the collision operator on secondary mesh using full boltzmann
    mm = size(nodes_uII,1)
    allocate (fII(1:size(nodes_uII,1)),fcolII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "FlyTrap1Donecell_DGV_MPI: Allocation error for variables (fII,fcolII)"
      end if
    fII=0 ! clean the array   
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"FlyTrap1Donecell_DGV_MPI:   slave processor", irank,&
         "MPI boradcast of solution fII from proc 0 returned error", ierr
     stop
    end if 
    ! we have the solution. Now we will evaluate the portion of the collision operator
    ! assigned to this processor by calling EvalCollisionPeriodicAPlus_DGVII_MPI
    call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator
    ! only portion of the fcolII is filled in the rest should be set = 0 in the subroutine.
    ! next we call a reduce operation to combine all components of the collision operator 
    fII=0 ! clean, just in case, 
    call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "GetRelaxRates1Donecell_DGV_MPI" subroutine
    if (ierr /= 0 ) then 
     print *,"FlyTrap1Donecell_DGV_MPI:  (0) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
     stop
    end if 
    deallocate (fII,fcolII) ! we do not really need these 
    ! done here
  case (402) ! this action code corresponds to evaluation of the collision operator on secondary mesh using decomposition
    mm = size(nodes_uII,1)
    allocate (fII(1:size(nodes_uII,1)),fcolII(1:size(nodes_uII,1)),dpbuff(1:5), stat=loc_alloc_stat)
      if (loc_alloc_stat >0) then 
       print *, "FlyTrap1Donecell_DGV_MPI: Allocation error for variables (fII,fcolII,dpbuff)"
      end if
    fII=0 ! clean the array   
    call mpi_bcast (fII,mm,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"FlyTrap1Donecell_DGV_MPI:   slave processor", irank,&
         "MPI broadcast (2) of solution fII from proc 0 returned error", ierr
     stop
    end if
    dpbuff=0 ! clean, just in case
    call mpi_bcast (dpbuff,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    if (ierr /= 0 ) then 
     print *,"FlyTrap1Donecell_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
     stop
    end if
    LocDens=dpbuff(1)
    LocUbar=dpbuff(2)
    LocVbar=dpbuff(3)
    LocWbar=dpbuff(4)
    LocTempr=dpbuff(5)
    deallocate(dpbuff) 
    ! we have the solution. Now we will evaluate the portion of the collision operator
    ! assigned to this processor by using decomposition approach 
    ! We will need these arrays to store the solution prohjected to secondary velocity mesh.
    allocate (fcol1II(1:size(nodes_uII,1)), fmII(1:size(nodes_uII,1)), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then 
     print *, "FlyTrap1Donecell_DGV_MPI: Allocation error for variables (fmII),(fcol1II)"
    end if 
    fmII = maxwelveldist(LocTempr,LocUbar,LocVbar,LocWbar,LocDens,nodes_uII,nodes_vII,nodes_wII) ! now we populate the maxwellian with the same macroparamters.
    fII = fII-fMII ! now fII contains the difference between the solution and the local maxwellian
    ! Now we will make to calls of the collision operator
    !!!!!!!!!!
    ! to make an MPI evaluation, uncomment the next two lines
    call EvalCollisionPeriodicAPlus_DGVII_MPI(fII,fcolII) ! fcolII contains the value of the collision operator dfdfA
    call EvalCollisionPeriodicMixedTermsA_DGVII_MPI(fII,fMII,fcol1II) !
    fcolII=fcolII+fcol1II
    deallocate (fcol1II,fMII)
    ! only portion of the fcolII is filled in the rest should be set = 0 in the subroutine.
    ! next we call a reduce operation to combine all components of the collision operator 
    fII=0 ! clean, just in case, 
    call mpi_allreduce(fcolII,fII,mm,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "GetRelaxRates1Donecell_DGV_MPI" subroutine
    if (ierr /= 0 ) then 
     print *,"FlyTrap1Donecell_DGV_MPI:  (1) slave processor", irank,&
          "MPI allreduce of fcolII returned error", ierr
     stop
    end if 
    deallocate (fII,fcolII) ! we do not really need these 
    ! done here
  case (-777)
    print *,"FlyTrap1Donecell_DGV_MPI: slave processor", irank,&
              "end of work signal recieved. Exit now" 
  case default
    print *,"FlyTrap1Donecell_DGV_MPI: slave processor", irank,&
         "Exit Error: unknown action case, mpi_dgv_action=", mpi_dgv_action
           stop
 end select   
end do 
end subroutine FlyTrap1Donecell_DGV_MPI

end module DGV_mpiroutines
