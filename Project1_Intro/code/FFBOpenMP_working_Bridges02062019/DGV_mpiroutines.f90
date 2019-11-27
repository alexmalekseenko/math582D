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
 ! this subroutines makes the MPI fork for some preparatory staff for the evaluation of the collision integral 
 ! 
 ! Description: 
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

use DGV_commvar, only: numchnks,num_Acopies,nodes_u,A_capphi,procs_nodes_wld, &
                   lin_proc_nodes_wld,MPI_LINEAR_COLL_WORLD, MPI_LINEAR_COLL_WORLD_GROUP,&
                   num_lin_proc,MPI_Acopy_UNIVS,MPI_Acopy_UNIVS_GROUP,Acopy_wrkld,procs_nodes_wld_Acopy_addr,&
                   MPI_LinFmA_UNIVS_GROUP, MPI_LinFmA_UNIVS, My_Acopy_Univ_indx, nodes_proc_groups, &
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
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer ! this buffer will be used in sending and receiving the work load.
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
call MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_WORLD_GROUP, ierr) ! get the hanlde to the group of MPI_COMM_WORLD
ranges(1)=0; ranges(2)=num_lin_proc-1; ranges(3)=1; 
call MPI_GROUP_RANGE_INCL(MPI_COMM_WORLD_GROUP,1,ranges, MPI_LINEAR_COLL_WORLD_GROUP, ierr) 
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_highperf: MPI can not create group MPI_LINEAR_COLL_WORLD_GROUP on proc", irank, ". Error=", ierr
         stop
end if
! a quick check if the processor with irank=0 is still with the rank =0 in the new group.
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
! end check 
! next we create a new communicator: 
call MPI_COMM_CREATE (MPI_COMM_WORLD, MPI_LINEAR_COLL_WORLD_GROUP, MPI_LINEAR_COLL_WORLD, ierr)
if (ierr /= 0 ) then 
 print *,"PrepMPICollsnDGV_highperf: can't create communicator MPI_LINEAR_COLL_WORLD. Proc=", irank, ". Err=", ierr
 stop
end if
! End create communicator for consolidated handling of linearised collision operator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we prepare three important arrays: num_chnks,nodes_proc_groups  and procs_lin. These arrays determine the workload of nonlinear evaluation of 
! the collision operator and the workload for evaluation of the linearized collision operator, respectively. 
! first, let us work on the array chnks_copies(2,1:num_Acopies)
allocate (chnks_copies(2,1:num_Acopies), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables chnks_copies. stop"
  stop
end if
chnks_copies=0 ! nullify, just in case.
!! Next we will call the subroutine DivideRangeContinuous to divide the processors between copies of A:
!! for that we need to create the result array 
call DivideRangeContinuous(chnks_copies,1,mpicommworldsize-1)
!! now chnks_copies will contain the information on groups of processors that share one copy of A between them. 
!! 
!! next we distribute A between processors of one such group. Then we will continue with another groupd and so on 
!! untill all groups of processors receive their copy of operator A. 
!! Aslo, we need to assign a workload of velocity nodes for each group of processors. Later this will be used to assign individual workloads for 
!! each processor. We will proceed chunk by chunk... 
!! now we will find the nodes that assigned to the each processor's group...
 nn = size(nodes_u,1)
 !! first we divide velocity nodes between the groups of processors. We have num_Acopies groups altogether
 allocate (nodes_proc_groups(2,1:num_Acopies), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables nodes_proc_groups. stop"
  stop
 end if
 call DivideRangeContinuous(nodes_proc_groups,1,nn)
!! we nodes between the groups of processors. 
!! We make one more prearation for the future: we 
!! setup the array procs_lin that tells the workload for each processor involved in the evaluation of the 
!! linearized collision operator. We simply divide velocity nodes between the processors involved in the evaluation of linearied collision operator. 
 allocate (procs_lin(2,1:num_lin_proc), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "PrepMPICollsnDGV_highperf: Allocation error for variables total_receive_procs,total_send_procs,procs_lin"
  stop
 end if
!! we populate the linearized group of processors with assgined nodes
 call DivideRangeContinuous(procs_lin,1,nn)
!! end of creating proc_lin: all processors are divided in groups and each group is assgned some number of velocity nodes 
!! arrays procs_lin, nodes_proc_groups and chnks_copies have been created, 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! We now create universes to deal with the groups of processors sharing a single copy of A and the universes of processors that obtain 
!!! components of the linearized operator from the same group of nodes... 
!!! these universies will be used in communcation of the components of the linearized operator. 
!!!
!! Begin creating the Acopy_UNIV universes:
!! first, we allocate storage for the universes... 
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
  ! a quick check if the processor with irank=0 is still with the rank =0 in the new group.
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
 allocate (MPI_LinFmA_UNIVS(1:num_Acopies),MPI_LinFmA_UNIVS_GROUP(1:num_Acopies), stat=loc_alloc_stat)
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
  ! now we will start screenign the processors participating in the evaluation of the linearized operetor. 
  ! we will look at the nodes that are assigned to them. Then we look at the nodes that are assigned to the chunk 
  ! containing the processor of rank rashk_temp(1). If they do, then we add them to the ranks array... 
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
 !! the determination of each individual processor's job allocation in evaluation of the collisio integral. 
  
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
 !! 
 !!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!!!!
 !! Next we will set up the arrays  
 !! Acopy_Wrkld and procs_nodes_wld_Acopy_addr
 !! These arrays only make sense on slave nodes...
 !!
 !!!!!!!!!!!!!!!!!!!!!!
 if (irank>0) then 
 !! essentially for each Acopy universe, we need to write all the nodes from  to the Alin_proc_nodes_wld
  nn = nodes_proc_groups(2,My_Acopy_Univ_indx) - nodes_proc_groups(1,My_Acopy_Univ_indx) + 1 
  allocate (Acopy_Wrkld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (Acopy_Wrkld) on node",irank,"Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the lin_proc_nodes_wld
  do i = 1,nn 
   Acopy_Wrkld(i) = nodes_proc_groups(1,My_Acopy_Univ_indx) + i -1 !  Acopy_Wrkld(i) contains all velocity nodes beginning 
                   ! from  nodes_proc_groups(1,My_Acopy_Univ_indx) and ending at  nodes_proc_groups(2,My_Acopy_Univ_indx)
  end do 
  !! end populating the lin_proc_nodes_wld 
  !! now let us preparethe array procs_nodes_wld_Acopy_addr
  nn = size(procs_nodes_wld,1)
  allocate (procs_nodes_wld_Acopy_addr(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (procs_nodes_wld_Acopy_addr) on node",irank,"Stop"
   stop
  end if
  !! we will essentially go over the  procs_nodes_wld array and will find the matcning components in Acopy_Wrkld and record indices o
  !! of these components in 
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
  nn = procs_lin(2,irank+1) - procs_lin(1,irank+1)+1 ! this is the total number of velocity nodes that will be 
  allocate (lin_proc_nodes_wld(1:nn),stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then 
   print *, "PrepMPICollsnDGV_highperf: Allocation error for variable (lin_proc_nodes_wld) on  node", irank, ". Stop"
   stop
  end if
  !! finaly we need to write all the nodes to the lin_proc_nodes_wld
  do i = procs_lin(1,irank+1),procs_lin(2,irank+1)
   lin_proc_nodes_wld(i-procs_lin(1,irank+1)+1) = i ! lin_proc_nodes_wld contains all velocity nodes beginning from procs_lin(1,irank+1) and ending at procs_lin(2,irank+1)
  end do 
  !! end populating the lin_proc_nodes_wld 
  !! Next we will populate the array linprocndswld_Acopy_univs. The information that will be stored in this array is already in the 
  !! Array LinFmA_recv_univs. The data just needs to be re-sorted:
  !! In the buffer the first number is the number of records after that go pairs: (#processor, #node)
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
! NOT DEBUGGED!!! 
!
! SetMPIWorkLoadNodes_DGV
!
! This subroutine sets up the array -- procs_nodes_wld
! this array contains the numbers of all nodes for which the 
! collision integral will be evaluated. 
!
! Essentially, the code goes through all nodes and check if 
! the local A array has records that pertain to them. If it does the node 
! goes on the work list. To do that the code goes through all nodes. For each node it finds is canonical node
! from nodes_phican. Then it looks up the local A_capphi arrays and see if the corresponding canonical node is there. 
! IF so, the node goes on the list.  
!
!
! NOT DEBUGGED!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SetMPIWorkLoadNodes_DGV

use DGV_commvar, only: procs_nodes_wld,nodes_phican,A_capphi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: nodes_array_scratch ! this array will keep the selected nodes temporarily.
integer (I4B) :: i,mm,nds_ct ! scrap index..
integer :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!!!
mm=size(nodes_phican,1)
allocate (nodes_array_scratch(1:mm),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetMPIWorkLoadNodes_DGV: Allocation error for variable (nodes_array_scratch)"
  stop
  end if
! next we populate the array nodes_array_scratch ... 
nds_ct=0
do i=1,mm
 if (A_capphi(nodes_phican(i))>0) then 
  nds_ct=nds_ct+1
  nodes_array_scratch(nds_ct)=i
 end if 
end do 
! now we record the accumulated information in the array procs_nodes_wld
! First we check that this array exists already. If so, something went wrong and we stop
if (size(procs_nodes_wld,1) > 0) then 
 print *,"SetMPIWorkLoadNodes_DGV: Error. The array procs_nodes_wld already exists. Stop."
 stop
end if 
! Now if the array does not exist, we will create it and move the selected nodes to it:
allocate (procs_nodes_wld(1:nds_ct),stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetMPIWorkLoadNodes_DGV: Allocation error for variable (procs_nodes_wld. Stop)"
  stop
  end if
! next we populate the array procs_nodes_wld ... 
procs_nodes_wld = nodes_array_scratch(1:nds_ct) 
!
deallocate(nodes_array_scratch)
!
end subroutine SetMPIWorkLoadNodes_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SetMPIProcsNodesWorkLoad(procs_nodes_wld,chnks_copies,nodes_phican)
!
! This subroutine sets up the array -- proc_nodes_wld
! this array contains the rows numbers of all nodes for which the 
! collision integral will be evaluated. 
! the fist number in the row is the total number of nodes in that row
! zero corresponds to no record
! 
! the number of rows corresponds to the number of processors among which the workload need to be disctributed. 
! chnks_copies contains information about how workload need to be distributed. 
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

subroutine SetMPIProcsNodesWorkLoad(mpi_procs_nodes_wld,chnks_copies,total_A_capphi,nodes_phican)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:) :: mpi_procs_nodes_wld ! this is the array that will hold the work assignment for each processor
								! mpi_procs_nodes_wld(i,j) j is the processor 1=master
								! (1,j) contains the total number of nodes
								! (2:??,j) contains the numbers of nodes.
integer (I4B), dimension (:,:), intent (in) :: chnks_copies ! this stores numbers of groups of processors that have copies of A
integer (I4B), dimension (:), intent (in) :: total_A_capphi ! this stores the information about local copy of A on each process 1 bloc of nn records is on master process and is not real. the rest should be real
integer (I4B), dimension (:), intent (in) :: nodes_phican   ! this contaisn information about the number of the corresponding canonical node for each node 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: nodes_array_scratch ! this array will keep the selected nodes temporarily.
integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! scrap to divide procs between groups 
integer (I4B) :: i,jj,j,mm, loc_alloc_stat,addrs, nds_ct ! scrap index..
!!!!!!!!!!!!!!!!!!!!!!!!!
!! first we divide nodes between the groups of processors. 
allocate (nodes_proc_groups(2,1:size(chnks_copies,2)), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "SetMPIProcsNodesWorkLoad: Allocation error for variables nodes_proc_groups. stop"
  stop
  end if
!!
mm = size(nodes_phican,1)
call DivideRangeContinuous(nodes_proc_groups,1,mm)
!! we divided the nodes between the groups of processors. 
!! now we will go group by group and distribute the processors between nodes. 
mpi_procs_nodes_wld=0 ! wipe the result before beginning
! ready to poulate the workload array 
do j=1,size(chnks_copies,2) ! this look goes trhough each group... 
! 
do i = nodes_proc_groups(1,j),nodes_proc_groups(2,j) ! loop in the nodes in the group j
do jj = chnks_copies(1,j),chnks_copies(2,j) ! loop in the (velocity nodes?) processors? in the group j
 addrs = jj*mm + nodes_phican(i) !! recall that the total_A_capphi has its first block for the master process  and contains zeros there...
 if ( total_A_capphi(addrs)>0 ) then 
  nds_ct = mpi_procs_nodes_wld(1,jj+1)+1 !! jj+1 because the processors are numbered from  0 (master) and so on and the array numbered from 1
  mpi_procs_nodes_wld(1,jj+1) = nds_ct 
  mpi_procs_nodes_wld(nds_ct+1,jj+1) = i
 end if 
end do 
end do 
!
end do 
deallocate (nodes_proc_groups)
! we put here a quick consistencey check -- it mey be useful to catch problem with A arrays:
! we will add the number of nodes in the work array: if this number is less than the total number of 
! nodes -- data for some nodes is missing. 
!
nds_ct=0
do jj=2,size(mpi_procs_nodes_wld,2)
 nds_ct=nds_ct+mpi_procs_nodes_wld(1,jj)
end do 
if (nds_ct < mm) then 
 print *,"SetMPIProcsNodesWorkLoad: Not all nodes were allocated for work. Possible incorrect A-arrys. Stop"
 stop
end if 
! end check
end subroutine SetMPIProcsNodesWorkLoad


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
use DGV_commvar, only: procs_nodes_wld,nodes_phican,A_capphi,num_Acopies

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:,:) :: chnks_copies ! array that contains information about division of processors bewtenn different copies of operator A
integer (I4B), dimension (:,:) :: nodes_proc_groups ! array that contain information about division of velocity nodes between groups  of processors
integer, intent (in) :: mpicommworldsize ! the total number of processors in the common world
integer, intent (in) :: irank ! rank of this slave processr in the common world
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: mpiI4BArryBuffer ! scrap array to keep found velocity nodes ...
integer (I4B) :: i,mm, loc_alloc_stat, nds_ct, my_chnk ! scrap index..

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
 my_chnk=0
 do i=1,num_Acopies
  if ((irank >= chnks_copies(1,i)) .and. (irank<=chnks_copies(2,i))) then 
   my_chnk = i ! this contains the number of the chunk where this processor belongs
  end if 
 end do 
 !!!! check if we found such number
 if (my_chnk==0) then 
  print *, "SetMPILocProcNodesWorkLoad: can not locate chunk to which this slave proc belongs. stop"
  stop
 end if !!!! end check
 !! Now we woill look at the nodes that are assigned to this processor's group.
 !! these nodes are start at nodes_proc_groups(1,my_chnk) and end at nodes_proc_groups(2,my_chnk)
 !! next we will go over these processors and see if operator A has any records for it. 
 !! we now allocate the buffer to hold the velocity nodes workload 
  mm = nodes_proc_groups(2,my_chnk)-nodes_proc_groups(1,my_chnk) + 1  ! the total number of processors assigned to the group.
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
 do i = nodes_proc_groups(1,my_chnk),nodes_proc_groups(2,my_chnk) ! loop in the processors in the group j
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
!
end subroutine SetMPILocProcNodesWorkLoad_DGV

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

end module DGV_mpiroutines