! 09/05/08 Alex 
!
!  misc_setup_mod.f90
!
!  This file deals with initial set up. Essentially these are pieces of codes that were 
!  taken from the main program and formatted as subroutines to make the main program  
!  more managable 
!
!  8/1/2008 10:42:16 AM
!
module misc_setup_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine sets up the meshes. Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main prog ram. 
! 
! The subroutine can not be called before the parameters 
! N,M, mesh_x_uniform, mesh_u_uniform, moments_x_gauss_order, 
! moments_u_gauss_order are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles N,M.
!
! DO not call before arrays moments_x_gauss_nodes,moments_x_gauss_weights_x,moments_x_gauss_order
!          moments_u_gauss_order,moments_u_gauss_nodes,moments_u_gauss_weights are selected! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Set1DHOVMesh 
use common_variables_mod, only: xmesh, umesh, N, M, x_left, x_right, u_left, u_right, &
                   mesh_x_uniform, mesh_u_uniform,x_nonuniform_mesh_type,u_nonuniform_mesh_type,&
                   moments_x_gauss_nodes,moments_u_gauss_nodes

intrinsic Min, Max, MinVal
                   
integer :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) :: dx, du ! local temp variables to store 1/3 smallest distance between nodes.
 ! first the mesh in x
 if (mesh_x_uniform) then 
  allocate (xmesh(0:N), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
  stop
  end if
     !
  xmesh = (/ (x_left + (x_right - x_left )/Real(N,DP)*i, i=0,N) /)
 else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  select case (x_nonuniform_mesh_type)
  case (1) ! use "the gauss nodes for evaluating moments in x" to set up the "mesh in x"
           ! interval [x_left,x_right] is divided in subintervals as to have gauss nodes at centerpoints 
           ! (and some extra points need to be introduced to make this possible). 
  mxgo = size(moments_x_gauss_nodes) ! temp remember order of the guass method for integration of moments .... 
  dx=min(minval(moments_x_gauss_nodes(2:mxgo) - moments_x_gauss_nodes(1:mxgo-1)),&
                 moments_x_gauss_nodes(1)+Real(1,DP), Real(1,DP)-moments_x_gauss_nodes(mxgo))/Real(3,DP)
  !! allocate memory for the nodes in x
     allocate (xmesh(0:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
     stop
     end if
     !
  !! Now we set up the mesh in x:
  xmesh(0)=x_left
  do i=1,mxgo
  xmesh(2*i-1) =(x_right+x_left)/Real(2,DP) + (x_right-x_left)/Real(2,DP)*(moments_x_gauss_nodes(i) - dx)
  xmesh(2*i) = (x_right+x_left)/Real(2,DP) + (x_right-x_left)/Real(2,DP)*(moments_x_gauss_nodes(i) + dx)
  end do 
  xmesh(mxgo*2+1) = x_right
  N=mxgo*2+1 ! the new value of (N) in case somebody uses it directly... 
  !! mesh in x is ready
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Modified Alex 07/07/09 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/8-1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=6 !!! 
  !! allocate memory for the nodes in x
  if (N<6) then 
    print *, "Set1DHOVMesh: For the type 2 nonuniform mesh in x the number of mesh points has to be >=6 (xmesh)"
    stop  
  end if  
    allocate (xmesh(0:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
    !
    dx=(x_right-x_left)/(Real(N-4,DP)-Real(1,DP)/Real(4,DP))
  !! Now we will set up the mesh in variable x:
    xmesh(0)=x_left
    xmesh(1)=x_left+dx/Real(8,DP)
    xmesh(2)=xmesh(1)+dx/Real(4,DP)
    xmesh(3)=xmesh(2)+dx/Real(2,DP)
     do i=4,N-4
      xmesh(i)=xmesh(i-1)+dx
     end do 
    xmesh(N)=x_right
    xmesh(N-1)=xmesh(N)-dx/Real(8,DP)
    xmesh(N-2)=xmesh(N-1)-dx/Real(4,DP)
    xmesh(N-3)=xmesh(N-2)-dx/Real(2,DP)         
  !! end of the non-uniform mesh type 2. 
  case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=4 !!! 
  !! allocate memory for the nodes in x
  if (N<4) then 
    print *, "Set1DHOVMesh: For the type 3 nonuniform mesh in x the number of mesh points has to be >=4 (xmesh)"
    stop  
  end if  
    allocate (xmesh(0:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
    !
    dx=(x_right-x_left)/(Real(N-2,DP)-Real(1,DP)/Real(2,DP))
  !! Now we will set up the mesh in variable x:
    xmesh(0)=x_left
    xmesh(1)=x_left+dx/Real(4,DP)
    xmesh(2)=xmesh(1)+dx/Real(2,DP)
     do i=3,N-3
      xmesh(i)=xmesh(i-1)+dx
     end do 
    xmesh(N)=x_right
    xmesh(N-1)=xmesh(N)-dx/Real(4,DP)
    xmesh(N-2)=xmesh(N-1)-dx/Real(2,DP)
  !!! end of the non-uniform mesh type 3. 
  !! End modified Alex 07/07/09
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  case default 
    allocate (xmesh(0:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
     !
    xmesh = (/ (x_left + (x_right-x_left)/Real(3,DP)*i, i=0,3) /) 
    N=3 ! the new value of (N) in case somebody uses it directly... 
    print *, "Set1DHOVMesh: unsupported choice of non-uniform mesh. x_nonuniform_mesh_type= ", x_nonuniform_mesh_type
end select
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 end if  
 ! now the mesh in u  
 if (mesh_u_uniform) then 
  allocate (umesh(0:M), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "Set1DHOVMesh: Allocation error for variable (umesh)"
  stop
  end if
     !
  umesh = (/ (u_left + (u_right-u_left)/Real(M,DP)*i, i=0,M) /)
 else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  select case (u_nonuniform_mesh_type)
  case (1) ! use "the gauss nodes for evaluating moments in u" to set up the "mesh in u"
           ! interval [u_left,u_right] is divided in subintervals as to have gauss nodes at centerpoints 
           ! (and some extra points need to be introduced to make this possible). 
  mugo = size(moments_u_gauss_nodes) ! temp remember order of the guass method for integration of moments .... 
  du=min(minval(moments_u_gauss_nodes(2:mugo) - moments_u_gauss_nodes(1:mugo-1)),&
                 moments_u_gauss_nodes(1) + Real(1,DP), Real(1,DP) - moments_u_gauss_nodes(mugo))/Real(3,DP)
  !! allocate memory for the nodes in u
     allocate (umesh(0:mugo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMesh: Allocation error for variable (umesh)"
     stop
     end if
     !
  !! Now we set up the mesh in u:
  umesh(0)=u_left
  do i=1,mugo
  umesh(2*i-1) = (u_right+u_left)/Real(2,DP) + (u_right-u_left)/Real(2,DP)*(moments_u_gauss_nodes(i) - du)
  umesh(2*i) =  (u_right+u_left)/Real(2,DP) + (u_right-u_left)/Real(2,DP)*(moments_u_gauss_nodes(i) + du)
  end do 
  umesh(mugo*2+1) = u_right
  M=mugo*2+1 ! the new value of (M) in case somebody uses it directly... 
  !! mesh in u is ready
  case default 
     allocate (umesh(0:3), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMesh: Allocation error for variable (umesh)"
     stop
     end if
     !
  umesh = (/ (u_left + (u_right-u_left)/Real(3,DP)*i, i=0,3) /)
  M=3! the new value of (M) in case somebody uses it directly... 
  print *, "Set1DHOVMesh: unsupported choice of non-uniform mesh. u_nonuniform_mesh_type= ", u_nonuniform_mesh_type
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 05/12/09 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! Only Static wall is implemented at this time! 
           ! 
 ! first, let us check that u_wall =0  is in the velocity interval                     
  if ((u_left > 0.0_DP) .or. (u_right < 0.0_DP)) then 
    u_left = min (0.0_DP,u_left)
    u_right = max (0.0_DP,u_right)
  end if                    
 ! now we reserve the space for the mesh points: 
   allocate (umesh(0:M), stat=loc_alloc_stat)
   !
   if (loc_alloc_stat >0) then 
   print *, "Set1DHOVMesh: Allocation error for variable (umesh)"
   stop
   end if
   !
  ! now we start to fill in the mesh:
  i=-M
  j=0
  umesh(0)=u_left
  do while ( u_right - (u_right-u_left)/Real(M,DP)*i > (u_right-u_left)/Real(M,DP)/2 - 0.00000000001_DP )    
  if (i*(u_right-u_left)/Real(M,DP)-u_left > (u_right-u_left)/Real(M,DP)/2) then 
  j=j+1
  umesh(j)=i*(u_right-u_left)/Real(M,DP)
  end if 
  i=i+1
  end do 
  umesh(M)=u_right
  ! now we have filled in the mesh:
!!! End Modified Alex 05/12/09 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 07/11/11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will include the velocity of the wall, u_{w}=0 as a mesh point
           ! also it will surround $u=0$ with small cells. The sizes of the sell are obtained 
           ! by refining the size of biggest cells $smr^{sml}$ times. Here $smr$ -- is the 
           ! refinement factor and $sml$ is the total level of reinements. 
           ! Only Static wall is implemented at this time! 
           ! 
 ! First we set values for $smr$ and $sml$ 
 smr=10
 sml=2
 if (2*sml > M-2) then
  print *, "Set1DHOVMesh: selected number of refinements 'sml' of mesh is too large for the given M"
  stop
 end if
   ! now we need to evaluate the size of the largest cell. 
 du=(u_right-u_left)/((M-2*sml) + 2.0_DP/Real(smr,DP)* &
             (1.0_DP/(Real(smr,DP)**sml)-1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP))
   ! now we reserve the space for the mesh points: 
   allocate (umesh(0:M), stat=loc_alloc_stat)
   !
   if (loc_alloc_stat >0) then 
   print *, "Set1DHOVMesh: Allocation error for variable (umesh)"
   stop
   end if
   !
   ! now we start to fill in the mesh:
  j=0
  umesh(0)= u_left
  do while (umesh(j) < - du*(1.0_DP/(Real(smr,DP)**sml) - 1.0_DP)/(1.0_DP/Real(smr,DP)-1.0_DP)/Real(smr,DP) - 0.00000000001_DP )    
  j=j+1
  umesh(j)=umesh(j-1)+du
  end do
  ! Now we fill in the small cells
  do i=1,sml
  j=j+1
  umesh(j) = umesh(j-1) + du/(Real(smr,DP)**i)
  end do 
  umesh(j)=0.0_DP ! Just toi have it exact --- otherwise it is not zero, but 10^-13
  do i=sml,1,-1
  j=j+1
  umesh(j) = umesh(j-1) + du/(Real(smr,DP)**i)
  end do 
  ! now the rest of the cells.
  do while (umesh(j) < u_right - 0.00000000001_DP)
  j=j+1
  umesh(j)=umesh(j-1)+du
  end do
  ! now we have filled in the mesh:
!!! End Modified Alex 07/11/11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end select
 end if  
 end subroutine Set1DHOVMesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine Set1DHOVspatMesh
!
! This is a copy of the above subroutine in which the construction of the meshes in velocity space was removed. 
!
! This subroutine sets up the meshes. Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters 
! N,M, mesh_x_uniform, mesh_u_uniform, moments_x_gauss_order, 
! moments_u_gauss_order are selected                               !!!!! 
! CALL this SUBROUTINE before other subroutines use varaibles N,M.
!
! DO not call before arrays moments_x_gauss_nodes,moments_x_gauss_weights_x,moments_x_gauss_order
!          moments_u_gauss_order,moments_u_gauss_nodes,moments_u_gauss_weights are selected! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Set1DHOVspatMesh 
use common_variables_mod, only: xmesh, N,x_left, x_right, &
                   mesh_x_uniform, x_nonuniform_mesh_type,&
                   moments_x_gauss_nodes
   
use DGV_sf02, only: x1,x2               

intrinsic Min, Max, MinVal
                   
integer :: loc_alloc_stat ! to keep allocation status
integer (I4B) ::  i,j       ! local counter 
integer (I4B) :: sml,smr,x1_i,x2_i   ! parameters for umesh with small cells around zero sml -- levels of refinment and smr -- ratio of refinement
integer :: mxgo,mugo      ! local variables to keep the order of the gauss integration for the moments
real (DP) :: dx, du,xscrp ! local temp variables to store 1/3 smallest distance between nodes.
 ! first the mesh in x
 if (mesh_x_uniform) then 
  allocate (xmesh(0:N), stat=loc_alloc_stat)
     !
  if (loc_alloc_stat >0) then 
  print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
  stop
  end if
     !
  xmesh = (/ (x_left + (x_right - x_left )/Real(N,DP)*i, i=0,N) /)
 else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  select case (x_nonuniform_mesh_type)
  case (1) ! use "the gauss nodes for evaluating moments in x" to set up the "mesh in x"
           ! interval [x_left,x_right] is divided in subintervals as to have gauss nodes at centerpoints 
           ! (and some extra points need to be introduced to make this possible). 
  mxgo = size(moments_x_gauss_nodes) ! temp remember order of the guass method for integration of moments .... 
  dx=min(minval(moments_x_gauss_nodes(2:mxgo) - moments_x_gauss_nodes(1:mxgo-1)),&
                 moments_x_gauss_nodes(1)+Real(1,DP), Real(1,DP)-moments_x_gauss_nodes(mxgo))/Real(3,DP)
  !! allocate memory for the nodes in x
     allocate (xmesh(0:mxgo*2+1), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
     stop
     end if
     !
  !! Now we set up the mesh in x:
  xmesh(0)=x_left
  do i=1,mxgo
  xmesh(2*i-1) =(x_right+x_left)/Real(2,DP) + (x_right-x_left)/Real(2,DP)*(moments_x_gauss_nodes(i) - dx)
  xmesh(2*i) = (x_right+x_left)/Real(2,DP) + (x_right-x_left)/Real(2,DP)*(moments_x_gauss_nodes(i) + dx)
  end do 
  xmesh(mxgo*2+1) = x_right
  N=mxgo*2+1 ! the new value of (N) in case somebody uses it directly... 
  !! mesh in x is ready
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Modified Alex 07/07/09 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  case (2) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/8-1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=6 !!! 
  !! allocate memory for the nodes in x
  if (N<6) then 
    print *, "Set1DHOVMesh: For the type 2 nonuniform mesh in x the number of mesh points has to be >=6 (xmesh)"
    stop  
  end if  
    allocate (xmesh(0:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
    !
    dx=(x_right-x_left)/(Real(N-4,DP)-Real(1,DP)/Real(4,DP))
  !! Now we will set up the mesh in variable x:
    xmesh(0)=x_left
    xmesh(1)=x_left+dx/Real(8,DP)
    xmesh(2)=xmesh(1)+dx/Real(4,DP)
    xmesh(3)=xmesh(2)+dx/Real(2,DP)
     do i=4,N-4
      xmesh(i)=xmesh(i-1)+dx
     end do 
    xmesh(N)=x_right
    xmesh(N-1)=xmesh(N)-dx/Real(8,DP)
    xmesh(N-2)=xmesh(N-1)-dx/Real(4,DP)
    xmesh(N-3)=xmesh(N-2)-dx/Real(2,DP)         
  !! end of the non-uniform mesh type 2. 
  case (3) ! This mesh is to be used with diffusive boundary conditions
           ! and will make 1/4-1/2 refinement near the walls
           ! Only Static wall is implemented at this time! 
           ! Requires N>=4 !!! 
  !! allocate memory for the nodes in x
  if (N<4) then 
    print *, "Set1DHOVMesh: For the type 3 nonuniform mesh in x the number of mesh points has to be >=4 (xmesh)"
    stop  
  end if  
    allocate (xmesh(0:N), stat=loc_alloc_stat)
    !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
    !
    dx=(x_right-x_left)/(Real(N-2,DP)-Real(1,DP)/Real(2,DP))
  !! Now we will set up the mesh in variable x:
    xmesh(0)=x_left
    xmesh(1)=x_left+dx/Real(4,DP)
    xmesh(2)=xmesh(1)+dx/Real(2,DP)
     do i=3,N-3
      xmesh(i)=xmesh(i-1)+dx
     end do 
    xmesh(N)=x_right
    xmesh(N-1)=xmesh(N)-dx/Real(4,DP)
    xmesh(N-2)=xmesh(N-1)-dx/Real(2,DP)
  !!! end of the non-uniform mesh type 3. 
  !! End modified Alex 07/07/09
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    case (4) ! This mesh is to be in normal shoc wave problem. It will use values x_1 and x_2 from DGV_sf02
    ! the cells in near x_1 and x_2 will be divided by half
    dx=(x_right-x_left)/(Real(N,DP)) ! this is the basis value of N- this determines the size of the coarsest cells. 
    ! next we will find out what is the number of interval that contains x_1 and  
    ! check error in data
    if (x_left > x1) then 
     print *,"Set1DHOVMesh: x1 in DGV_sf02 seems incompatible with x_left. Stop"
     stop
    end if  
    ! now let us fins on what interval x1 falls
    xscrp=x_left+dx
    i=1
    do while ((xscrp < x1) .and. (i<N)) 
     i = i + 1
     xscrp = xscrp + dx    
    end do
    ! check error in data 
    if ((i==N) .and. (xscrp < x1)) then 
     print *,"Set1DHOVMesh: x1 in DGV_sf02 seems incompatible with x_right. Stop"
     stop
    end if
    x1_i=i ! this is the number of interval containing x1
    !!!!! continue to find the number of the interval that contains x2
    do while ((xscrp < x2) .and. (i<N)) 
     i = i + 1
     xscrp = xscrp + dx    
    end do
    ! check error in data 
    if ((i==N) .and. (xscrp < x2)) then 
     print *,"Set1DHOVMesh: x2 in DGV_sf02 seems incompatible with x_right. Stop"
     stop
    end if
    x2_i=i ! this is the number of interval containing x1
    !!!!! continue to find the number of the interval that contains x2
    !!!!! we found the intervals, now we need to construct the binary tree meshes:
    if (x1_i > x2_i)  then  
     print *, "Set1DHOVMesh: problems constructing the xmesh. stop"
     stop
    end if
    ! Next case is when the cells containing x1 and x2 are the same or next to each other or separated by one or two cells. In this case we will divide each affected cell into two 
    if ((x1_i <= x2_i) .and. (x1_i > x2_i-4)) then  
     allocate (xmesh(0:N+(x2_i-x1_i+1)), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then 
      print *, "Set1DHOVMesh: Allocation error for variable (xmesh) 40"
      stop
     end if
     xmesh(0) = x_left
     do i=1,x1_i-1
      xmesh(i) = xmesh(i-1)+dx
     end do 
     do i=0,x2_i-x1_i
      xmesh(x1_i+2*i)=xmesh(x1_i-1+2*i) + dx/2.0_DP
      xmesh(x1_i+2*i+1)=xmesh(x1_i+2*i) + dx/2.0_DP
     end do 
     do i=x2_i+2+(x2_i-x1_i),N+(x2_i-x1_i+1)
      xmesh(i) = xmesh(i-1)+dx
     end do
     N = N+(x2_i-x1_i+1)
    end if 
    ! Next case is when the cells containing x1 and x2 are separated by at least three cells. In this case we only subdivide four cells 
    if ((x1_i <= x2_i) .and. (x1_i < x2_i-3)) then  
     allocate (xmesh(0:N+4), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then 
      print *, "Set1DHOVMesh: Allocation error for variable (xmesh) 41"
      stop
     end if
     xmesh(0) = x_left
     do i=1,x1_i-1
      xmesh(i) = xmesh(i-1)+dx
     end do 
     xmesh(x1_i)=xmesh(x1_i-1) + dx/2.0_DP
     xmesh(x1_i+1)=xmesh(x1_i) + dx/2.0_DP
     xmesh(x1_i+2)=xmesh(x1_i+1) + dx/2.0_DP
     xmesh(x1_i+3)=xmesh(x1_i+2) + dx/2.0_DP
     do i=x1_i+4,x2_i
      xmesh(i) = xmesh(i-1)+dx
     end do  
     xmesh(x2_i+1)=xmesh(x2_i) + dx/2.0_DP
     xmesh(x2_i+2)=xmesh(x2_i+1) + dx/2.0_DP
     xmesh(x2_i+3)=xmesh(x2_i+2) + dx/2.0_DP
     xmesh(x2_i+4)=xmesh(x2_i+3) + dx/2.0_DP
     do i=x2_i+5,N+4
      xmesh(i) = xmesh(i-1)+dx
     end do  
     N=N+4
    end if 
    
  case (5) ! This mesh is to be in normal shoc wave problem. It will use values x_1 and x_2 from DGV_sf02
    ! the cells in near x_1 and x_2 will be divided by half and possibly by another half
    dx=(x_right-x_left)/(Real(N,DP)) ! this is the basis value of N- this determines the size of the coarsest cells. 
    ! next we will find out what is the number of interval that contains x_1 and  
    ! check error in data
    if (x_left > x1) then 
     print *,"Set1DHOVMesh: x1 in DGV_sf02 seems incompatible with x_left. Stop"
     stop
    end if  
    ! now let us fins on what interval x1 falls
    xscrp=x_left+dx
    i=1
    do while ((xscrp < x1) .and. (i<N)) 
     i = i + 1
     xscrp = xscrp + dx    
    end do
    ! check error in data 
    if ((i==N) .and. (xscrp < x1)) then 
     print *,"Set1DHOVMesh: x1 in DGV_sf02 seems incompatible with x_right. Stop"
     stop
    end if
    x1_i=i ! this is the number of interval containing x1
    !!!!! continue to find the number of the interval that contains x2
    do while ((xscrp < x2) .and. (i<N)) 
     i = i + 1
     xscrp = xscrp + dx    
    end do
    ! check error in data 
    if ((i==N) .and. (xscrp < x2)) then 
     print *,"Set1DHOVMesh: x2 in DGV_sf02 seems incompatible with x_right. Stop"
     stop
    end if
    x2_i=i ! this is the number of interval containing x1
    !!!!! continue to find the number of the interval that contains x2
    !!!!! we found the intervals, now we need to construct the binary tree meshes:
    if (x1_i > x2_i)  then  
     print *, "Set1DHOVMesh: problems constructing the xmesh. stop"
     stop
    end if
    ! Next case is when the cells containing x1 and x2 are the same or next to each other or separated by one or two cells. In this case we will divide each affected cell into two 
    if ((x1_i <= x2_i) .and. (x1_i > x2_i-4)) then  
     allocate (xmesh(0:N+(x2_i-x1_i+1)), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then 
      print *, "Set1DHOVMesh: Allocation error for variable (xmesh) 40"
      stop
     end if
     xmesh(0) = x_left
     do i=1,x1_i-1
      xmesh(i) = xmesh(i-1)+dx
     end do 
     do i=0,x2_i-x1_i
      xmesh(x1_i+2*i)=xmesh(x1_i-1+2*i) + dx/2.0_DP
      xmesh(x1_i+2*i+1)=xmesh(x1_i+2*i) + dx/2.0_DP
     end do 
     do i=x2_i+2+(x2_i-x1_i),N+(x2_i-x1_i+1)
      xmesh(i) = xmesh(i-1)+dx
     end do  
     N = N + (x2_i-x1_i+1)
    end if 
    ! Next case is when the cells containing x1 and x2 are separated by at least three cells. In this case we only subdivide four cells by half and fourse of the half cells by half 
    if ((x1_i <= x2_i) .and. (x1_i < x2_i-3)) then  
     allocate (xmesh(0:N+8), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then 
      print *, "Set1DHOVMesh: Allocation error for variable (xmesh) 41"
      stop
     end if
     xmesh(0) = x_left
     do i=1,x1_i-1
      xmesh(i) = xmesh(i-1)+dx
     end do 
     xmesh(x1_i)=xmesh(x1_i-1) + dx/2.0_DP
     xmesh(x1_i+1)=xmesh(x1_i) + dx/4.0_DP
     xmesh(x1_i+2)=xmesh(x1_i+1) + dx/4.0_DP
     xmesh(x1_i+3)=xmesh(x1_i+2) + dx/4.0_DP
     xmesh(x1_i+4)=xmesh(x1_i+3) + dx/4.0_DP
     xmesh(x1_i+5)=xmesh(x1_i+4) + dx/2.0_DP
     do i=x1_i+6,x2_i+2
      xmesh(i) = xmesh(i-1)+dx
     end do  
     xmesh(x2_i+3)=xmesh(x2_i+2) + dx/2.0_DP
     xmesh(x2_i+4)=xmesh(x2_i+3) + dx/4.0_DP
     xmesh(x2_i+5)=xmesh(x2_i+4) + dx/4.0_DP
     xmesh(x2_i+6)=xmesh(x2_i+5) + dx/4.0_DP
     xmesh(x2_i+7)=xmesh(x2_i+6) + dx/4.0_DP
     xmesh(x2_i+8)=xmesh(x2_i+7) + dx/2.0_DP
     do i=x2_i+9,N+8
      xmesh(i) = xmesh(i-1)+dx
     end do 
     N = N + 8 
    end if 
  case default 
    allocate (xmesh(0:3), stat=loc_alloc_stat)
     !
    if (loc_alloc_stat >0) then 
    print *, "Set1DHOVMesh: Allocation error for variable (xmesh)"
    stop
    end if
     !
    xmesh = (/ (x_left + (x_right-x_left)/Real(3,DP)*i, i=0,3) /) 
    N=3 ! the new value of (N) in case somebody uses it directly... 
    print *, "Set1DHOVMesh: unsupported choice of non-uniform mesh. x_nonuniform_mesh_type= ", x_nonuniform_mesh_type
end select
!!!!! End Modified Alex 10/18/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
 end if  

end subroutine Set1DHOVspatMesh



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set1DHOVSqp
! This subroutine sets up the Sqp marix. Do not detouch from the main program 
! 
! Sqp,j = \int_{-1}^{1} \varphi_{p,j}(x)\partial_{x}\varphi_{q,j}(x) d x
!
! The values are independent of j
!
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters k is selected!! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Set1DHOVSqp
use common_variables_mod, only: k,Sqp
use spectral_tools_mod
!!!
integer :: loc_alloc_stat ! to keep allocation status
!!!
!  Sqp,j = \int_{-1}^{1} \varphi_{p,j}(x)\partial_{x}\varphi_{q,j}(x) d x
!  first, we prepare the space for Sqp:
allocate (Sqp(0:k,0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVSqp: Allocation error for variable (Sqp)"
     end if 
     !
!  Now we calculate Sqp
call IntLeg_dP_q_P_p (Sqp) ! 

!!! 
end subroutine Set1DHOVSqp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set1DHOVinitftil
!
! This subroutine sets up the initial data for ftil1 and ftil2 arrays. Do not detouch from the main program 
! 
!
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters umesh, xmesh, s,k, max_degree, selected_exact_sol, curr_time are selected!! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Set1DHOVinitftil
use common_variables_mod, only: k,xmesh,N,ftil,curr_time

use DGV_commvar, only: nodes_u,nodes_v,nodes_w 

use spectral_tools_mod

use DGV_sf02, only: f_1D3D
!!!
intrinsic min
!!!
integer :: loc_alloc_stat ! to keep allocation status
!real (DP), dimension (0:s) :: ff1 ! dump variable to compute matrix product
!real (DP), dimension (0:s) :: ff2 ! dump variable to compute matrix product
!!!
!  first, we prepare the space for the solution: ftil(p,j,m):
                                                      ! -- p is the index in the basis functions in x 
                                                      ! -- m is the index in the basis functions in u/nodes in u
                                                      ! -- j is the cell in x
allocate (ftil(0:k,size(nodes_u,1),N), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVinitftil: Allocation error for variable (ftil)"
     end if 
     !
! now we must set up ftil 
! the first step is to calculate the spectal coefficients of the intial data
ftil = DecompLeg1D3D(xmesh,nodes_u,nodes_v,nodes_w,k,f_1D3D,curr_time)	! the main distribution, 
end subroutine Set1DHOVinitftil

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SelectCFL_1DHOV
!
! This subroutine sets the CFL number based on the order of the method. 
! Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters s,k are assigned values
! 
! Most of the varaibles are looked up in the common_varaibles_module
! 
subroutine SelectCFL_1DHOV
use common_variables_mod, only: rk,cfl
!!!
integer (I4B) :: the_degree_of_interest ! some derived number to use in selection of cfl
!!!
the_degree_of_interest = rk  ! the clf will be selected based on the order of approximation
                        ! the actual algorithm is unclear, I will use order of 
                        ! the runge Kutta integrator
!!!
select case (the_degree_of_interest)
  case (0)
  cfl = 1/Real(2,DP)
  case (1)
  cfl = 1/Real(4,DP)
  case (2)
  cfl = 1/Real(8,DP)
  case (3)
  cfl = 1/Real(16,DP)
  case (4)
  cfl = 1/Real(32,DP)
  case (5)
  cfl = 1/Real(64,DP)
  case (6)
  cfl = 1/Real(128,DP)
  case default 
  cfl = 1/Real(128,DP)
end select 
end subroutine SelectCFL_1DHOV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetQ2plus1Altones
!
! This subroutine sets the useful vector of cofficients (1,3,...,2q+1,...,2k+1). 
! Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameter k is assigned values
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetQ2plus1Altones
use common_variables_mod, only: k, q2plus1_vec, altones
integer :: q ! local counter
integer :: loc_alloc_stat ! to keep allocation status
!
allocate (q2plus1_vec(0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Setq2plus1: Allocation error for variable (q2plus1_vec)"
     end if 
     ! 
q2plus1_vec(:) = (/ (2*q+1, q=0,k ) /)
allocate (altones(0:k), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Setq2plus1: Allocation error for variable (altones)"
     end if 
     ! 
altones(:) = (/ ((-1)**q, q=0,k) /)
end subroutine SetQ2plus1Altones

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L2errorEval1D1D
!
! This subroutine evaluates the L_2 error of approximation of the (external) function 
! (ex_sol) and its 1D1D discrete representation by a spectral decomposition (coeff).
!
! EXPECTS 
! xmesh --- mesh points in x
! umesh --- mesh points in u 
! k --- order of the spectral represent in x
! s --- order of the spectral represenataion in u.
! max_deg --- max degree of the basis polynomial (in both x and u)
! t --- coordinate time (used in exact solution)
! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
! ex_sol_fun --- name of the external function that gives the exact solution
!             ! it expects ex_sol_fun(x,w,t), 
!             !  x is scalar, w is a vector, and t is a scalar
! refine_mesh --- coefficient of the mesh refinement : every cell will be divided into refine_mesh^2 parts
! gauss_nodes,gauss_weights -- nodes and weights of the gauss integration formula to use for the error evaluation 
!
! RETURNS
! 
! ans = L2 norm of the f-\hat{f} 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function L2errorEval1D1D(xmesh,umesh,k,s,max_deg,t,coeff,ex_sol_fun,refine_mesh_in,gauss_nodes,gauss_weigths) result (ans)
use spectral_tools_mod  ! this one has the assembling routine
!!
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in variable x
real (DP), dimension (0:), intent (in) :: umesh ! mesh points in variable u
integer (I4B), intent (in) :: k,s,max_deg ! order in x, order in u and max degree of polynomial approximation 
real (DP), intent (in) :: t ! time (used in the exact solution fuction call),
real (DP), dimension (0:,0:,:,:) :: coeff ! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
integer (I4B), intent (in) :: refine_mesh_in ! coefficient of mesh refinement 
real (DP), dimension (:), intent (in) :: gauss_nodes, gauss_weigths ! nodes and weights of guassian formula on [-1,1] interval!
real (DP) :: ans ! the value of the computed L_2 norm of the difference
integer (I4B) :: refine_mesh ! coefficient of mesh refinement
!
integer :: i,j,ci,cj ! local counters
real (DP), dimension (:), allocatable :: xxmesh, uumesh, wweights ! meshes and gauss weights where the spectral representation will be evaluated 
integer :: loc_alloc_stat ! to keep allocation status 
integer :: gaussorder ! local variable to keep the order of gauss integration for error evaluation
real (DP) :: du,dx,ans_tempi,ans_tempi_tempj ! to keep temp calculations
 ! interface block for the dummy function ex_sol_fun: (x-scalar, u-vector)
 ! first variable must be scalar and second must be vector! 
 interface 
   function ex_sol_fun (x, u, t) result (y)
   ! ATTN: NOT CLEAR IF this is legal.... 
    use nrtype ! needs to remind where to get these constatns from? 
   !
    real (DP), intent (in)  :: x     ! the value of variable (x) where the function needs to be evaluated  
    real (DP), dimension (:), intent (in)  :: u     ! vector of values in variable (u) where the function needs to be evaluated  
    real (DP), intent (in) :: t ! the time variable 
    real (DP), dimension (size(u)) :: y ! values of the function  
   end function ex_sol_fun
 end interface
!!!!!!!!
! just in case, check if garbage came in refine_mesh: ---
if (refine_mesh_in < 1) then  
   refine_mesh=1
   print *, "L2errorEval1D1D: supplied refine_mesh parameter is zero or negative. Set refine_mesh=1"
   else 
   refine_mesh=refine_mesh_in
end if 
gaussorder=size(gauss_nodes)
!! allocate space for inter cell meshes ... and nodes
 allocate (xxmesh(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (xxmesh)"
   end if 
allocate (uumesh(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (uumesh)"
   end if 
allocate (wweights(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (xxmesh)"
   end if 
!! 
! we prepare a larger array of weights (dublicates of gauss weights) to be used on subdivided cells
do cj=1,refine_mesh
 wweights(1+(cj-1)*gaussorder:cj*gaussorder) = Real (gauss_weigths, DP)
end do 
! ready to calculate the L_2 error:
ans=0 
do i=1,size(umesh)-1      !loop in cells in u
 ! now each cell [umesh(i-1),umesh(i)] is subdivided into (refine_mesh) parts:
 ! and the new inside-cell meshes are recorded. 
 du = ( umesh(i) - umesh(i-1) ) / Real(refine_mesh,DP)
 do ci=1,refine_mesh
  uumesh(1+(ci-1)*gaussorder:ci*gaussorder) = umesh(i-1) + du*(Real(ci,DP) - .5) + Real(gauss_nodes,DP)*du/2 
 end do
 ans_tempi=0 
 do j=1,size(xmesh)-1      !loop in cells in x
   ! now each cell [xmesh(j-1),xmesh(j)] is subdivided into (refine_mesh)^2 parts:
   ! and the new inside-cell meshes are recorded. 
   dx = ( xmesh(j) - xmesh(j-1) ) / Real(refine_mesh,DP)
   do cj=1,refine_mesh
   xxmesh(1+(cj-1)*gaussorder:cj*gaussorder) = xmesh(j-1) + dx*(Real(cj,DP) - .5) + Real(gauss_nodes,DP)*dx/2 
   end do  
   ! 
   ans_tempi_tempj=0 
   do cj=1,gaussorder*refine_mesh
   ans_tempi_tempj = ans_tempi_tempj + sum( ( (  &
   EvCellLeg1D1D(xxmesh(cj),Real(uumesh,DP),xmesh(j-1),xmesh(j),umesh(i-1),umesh(i),k,s,max_deg,coeff(:,:,i,j)) &
    - ex_sol_fun (xxmesh(cj), uumesh, t) )**2) *wweights)*wweights(cj)
   end do 
   !
   ans_tempi = ans_tempi + ans_tempi_tempj*dx
 end do  ! end loop in (j) --- cells in x
 ans = ans + ans_tempi*du
end do  ! end loop in cells in u
ans=sqrt(ans)/Real(2,DP)                   ! division by 2 comes from the integration
deallocate (uumesh,xxmesh,wweights)
!
end function L2errorEval1D1D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! added 09/05/08 Alex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set1DHOVLegGaussArraysSol
! This subroutine sets up the arrays that will keep gauss nodes, gauss weights and 
! the values of the Legende polynomials evaluated at gaussian nodes multiplied by 
! the Gaussian Weigths. These arrays will be used to evaluate the spectral 
! decomposition of the solution and boundary data
!
! Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameters k,s are selected!! 
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Set1DHOVLegGaussArraysSol
use common_variables_mod, only: k,s,x_gauss_nodes, x_gauss_weights,u_gauss_nodes,u_gauss_weights,& 
                                x_wleg_pol, u_wleg_pol 
use gaussian_mod
use basis_fun_mod
use poly_tool_mod
!!!
integer :: loc_alloc_stat ! to keep allocation status
integer :: i ! local counter
!!! 
!!! Prepare the Gaussian nodes and weights in x variables !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the first variable
  allocate (x_gauss_nodes(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (x_gauss_nodes)"
     end if 
     !
  allocate (x_gauss_weights(1:max(k+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (x_gauss_weights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the first variable    
  if (k == 0) then 
   x_gauss_nodes(:)= (/ 0.0_DP /)
   x_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),x_gauss_nodes,x_gauss_weights)     ! the number of nodes/presition of the quadrature depends on the size of gnodes,gweights
  end if 
!!! Prepare the values of the Legendre polynomials for integration in the first variable 
  allocate (x_wleg_pol(0:k,size(x_gauss_nodes)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (x_wleg_pol)"
     end if 
     !
!!! Evaluate the basis functions on the nodes
  do i=0,k
   x_wleg_pol(i,:) = EvalHorner (LegendrePoly (i), Real(x_gauss_nodes,DP))*Real(x_gauss_weights,DP)
  end do 
!!! End evaluation of basis functions for x variable on gauss nodes  
!!!
!!! Prepare the Gaussian nodes and weights in the second variable (u) !!!!!!!!!!!!!!!!!!!!!!
!!! Allocate arrays for gaussian nodes and weights to inegration in the second variable
  allocate (u_gauss_nodes(1:max(s+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (u_gauss_nodes)"
     end if 
     !
  allocate (u_gauss_weights(1:max(s+1,1)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (u_gauss_weights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the second variable    
  if (s == 0) then 
   u_gauss_nodes(:)= (/ 0.0_DP /)
   u_gauss_weights(:)= (/ 2.0_DP /)  
  else  
   call GauLeg(Real(-1,DP),Real(1,DP),u_gauss_nodes,u_gauss_weights)     ! the number of nodes/prescision of the quadrature depends on the size of gnodes,gweights
  end if
!!! Prepare the values of the Legendre polynomials for integration in the first variable 
  allocate (u_wleg_pol(0:s,size(u_gauss_nodes)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "Set1DHOVLegGaussArraysSol: Allocation error for variable (u_wleg_pol)"
     end if 
     !
!!! Evaluate the basis functions on the nodes and muliply by the weights
  do i=0,s
   u_wleg_pol(i,:) = EvalHorner (LegendrePoly (i), Real(u_gauss_nodes,DP))*Real(u_gauss_weights,DP)
  end do 
!!! End evaluation of basis functions for u variable on gauss nodes  

end subroutine Set1DHOVLegGaussArraysSol

!!!!!!!! end added 09/05/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! added 11/07/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L2NormEval1D1D
!
! This subroutine evaluates the L_2 error of approximation of the (external) function 
! (ex_sol) and its 1D1D discrete representation by a spectral decomposition (coeff).
!
! EXPECTS 
! xmesh --- mesh points in x
! umesh --- mesh points in u 
! k --- order of the spectral represent in x
! s --- order of the spectral represenataion in u.
! max_deg --- max degree of the basis polynomial (in both x and u)
! t --- coordinate time (used in exact solution)
! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
! refine_mesh --- coefficient of the mesh refinement : every cell will be divided into refine_mesh^2 parts
! gauss_nodes,gauss_weights -- nodes and weights of the gauss integration formula to use for the error evaluation 
!
! RETURNS
! 
! ans = L2 norm of the f 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function L2NormEval1D1D(xmesh,umesh,k,s,max_deg,coeff,refine_mesh_in,gauss_nodes,gauss_weigths) result (ans)
use spectral_tools_mod  ! this one has the assembling routine
!!
real (DP), dimension (0:), intent (in) :: xmesh ! mesh points in variable x
real (DP), dimension (0:), intent (in) :: umesh ! mesh points in variable u
integer (I4B), intent (in) :: k,s,max_deg ! order in x, order in u and max degree of polynomial approximation 
real (DP), dimension (0:,0:,:,:) :: coeff ! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
integer (I4B), intent (in) :: refine_mesh_in ! coefficient of mesh refinement 
real (DP), dimension (:), intent (in) :: gauss_nodes, gauss_weigths ! nodes and weights of guassian formula on [-1,1] interval!
real (DP) :: ans ! the value of the computed L_2 norm of the difference
integer (I4B) :: refine_mesh ! coefficient of mesh refinement
!
integer :: i,j,ci,cj ! local counters
real (DP), dimension (:), allocatable :: xxmesh, uumesh, wweights ! meshes and gauss weights where the spectral representation will be evaluated 
integer :: loc_alloc_stat ! to keep allocation status 
integer :: gaussorder ! local variable to keep the order of gauss integration for error evaluation
real (DP) :: du,dx,ans_tempi,ans_tempi_tempj ! to keep temp calculations
!!!!!!!!
! just in case, check if garbage came in refine_mesh: ---
if (refine_mesh_in < 1) then  
   refine_mesh=1
   print *, "L2errorEval1D1D: supplied refine_mesh parameter is zero or negative. Set refine_mesh=1"
   else 
   refine_mesh=refine_mesh_in
end if 
gaussorder=size(gauss_nodes)
!! allocate space for inter cell meshes ... and nodes
 allocate (xxmesh(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (xxmesh)"
   end if 
allocate (uumesh(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (uumesh)"
   end if 
allocate (wweights(1:gaussorder*refine_mesh), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "L2errorEval1D1D: Allocation error for variable (xxmesh)"
   end if 
!! 
! we prepare a larger array of weights (dublicates of gauss weights) to be used on subdivided cells
do cj=1,refine_mesh
 wweights(1+(cj-1)*gaussorder:cj*gaussorder) = Real (gauss_weigths, DP)
end do 
! ready to calculate the L_2 error:
ans=0 
do i=1,size(umesh)-1      !loop in cells in u
 ! now each cell [umesh(i-1),umesh(i)] is subdivided into (refine_mesh) parts:
 ! and the new inside-cell meshes are recorded. 
 du = ( umesh(i) - umesh(i-1) ) / Real(refine_mesh,DP)
 do ci=1,refine_mesh
  uumesh(1+(ci-1)*gaussorder:ci*gaussorder) = umesh(i-1) + du*(Real(ci,DP) - .5) + Real(gauss_nodes,DP)*du/2 
 end do
 ans_tempi=0 
 do j=1,size(xmesh)-1      !loop in cells in x
   ! now each cell [xmesh(j-1),xmesh(j)] is subdivided into (refine_mesh)^2 parts:
   ! and the new inside-cell meshes are recorded. 
   dx = ( xmesh(j) - xmesh(j-1) ) / Real(refine_mesh,DP)
   do cj=1,refine_mesh
   xxmesh(1+(cj-1)*gaussorder:cj*gaussorder) = xmesh(j-1) + dx*(Real(cj,DP) - .5) + Real(gauss_nodes,DP)*dx/2 
   end do  
   ! 
   ans_tempi_tempj=0 
   do cj=1,gaussorder*refine_mesh
   ans_tempi_tempj = ans_tempi_tempj + sum( (&
   EvCellLeg1D1D(xxmesh(cj),Real(uumesh,DP),xmesh(j-1),xmesh(j),umesh(i-1),umesh(i),k,s,max_deg,coeff(:,:,i,j))**2)&
                                            *wweights)*wweights(cj)
   end do 
   !
   ans_tempi = ans_tempi + ans_tempi_tempj*dx
 end do  ! end loop in (j) --- cells in x
 ans = ans + ans_tempi*du
end do  ! end loop in cells in u
ans=sqrt(ans)/Real(2,DP)                   ! division by 2 comes from the integration
deallocate (uumesh,xxmesh,wweights)
!
end function L2NormEval1D1D



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! L2NormEvalHOV1D
!
! This function computes the spectral coefficients of the representation from the spectral-characteristic coefficients
! and calls L_2 norm the solution evaluation.
! 
! Do not detouch from the main program
!
! This subrouting is highly dependent on the variables of the main probram and can not be called before 
! the variables k,s,max_deg,xmesh,umesh,Aml,selected_exact_sol,error_eval_refine_mesh,&
! error_gauss_nodes,error_gauss_weights,curr_time,M,N ... are initialized
!
! Most of the variables are looked up in the common_variables_mod 
!
! USES  
! xmesh --- mesh points in x
! umesh --- mesh points in u 
! k --- order of the spectral represent in x
! s --- order of the spectral represenataion in u.
! max_deg --- max degree of the basis polynomial (in both x and u)
! curr_time --- coordinate time (used in exact solution)
! coeff(p,m,i,j) -- coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
! error_eval_refine_mesh --- coefficient of the mesh refinement : every cell will be divided into refine_mesh^2 parts
! (error_eval_refine_mesh) is set in the common_variables_mod !!!
! error_gauss_nodes -- nodes of the gaussian formula to use for the error evaluation 
! error_gauss_weigths --- weigths of the gaussian formula to use for error evaluation
!
! RETURNS
! 
! ans --- is the value of L2-norm of the spectral solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function L2NormEvalHOV1D(ftil1,ftil2) result (ans)
use common_variables_mod, only: k,s,max_deg,xmesh,umesh,Aml,selected_exact_sol,error_eval_refine_mesh,&
    error_gauss_nodes,error_gauss_weights,curr_time,M,N 
intrinsic min
!!!
real (DP), dimension (0:,0:,:,:), intent (in) :: ftil1,ftil2 ! coeff(p,m,i,j) -- characteristic coefficients of spectral representation
      ! -- p is the index in the basis functions in x 
      ! -- m is the index in the basis functions in u
      ! -- i is the cell in u
      ! -- j is the cell in x
real (DP), dimension (2) :: ans ! the computed values of the norm
real (DP), dimension (0:k,0:s,M,N) :: f  ! auxiliary variable -- coefficients of spectral decomposition 
                                                     ! we will calculate them from the characteristic variables  
                                                     !  f(p,m,i,j)( ftil(p,m,i,j), fti2(p,m,i,j):)
                                                     ! -- p is the index in the basis functions in x 
                                                     ! -- m is the index in the basis functions in u
                                                     ! -- i is the cell in u
                                                     ! -- j is the cell in x
!!!
integer  :: i,j,l,q,m_count    ! local counters
!! 
!!! COMPONENT 1 
   ! convert the characteristic coefficients into non-characteristic, that is
   ! evaluates the coefficients of the spectral representation from the characteristic representation
   f=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f(q,l,i,j)= f(q,l,i,j) + Aml(l,m_count,i)*ftil1(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
   ! Call L_2 Norm Evaluation
   ans(1)=L2NormEval1D1D(xmesh,umesh,k,s,max_deg,f,error_eval_refine_mesh,&
            error_gauss_nodes,error_gauss_weights)
!!! COMPONENT 2
   ! convert the characteristic coefficients into non-characteristic, that is
   ! evaluate the coefficients of the spectral representation from the characteristic representation
   f=0
   do j=1,N
   do i=1,M
   do l=0,s
   do q=0,min(max_deg-l,k)
     do m_count=0,s
     f(q,l,i,j)= f(q,l,i,j) + Aml(l,m_count,i)*ftil2(q,m_count,i,j)
     end do   
   end do
   end do 
   end do 
   end do 
   ! Call L_2 norm evaluation 
   ans(2)=L2NormEval1D1D(xmesh,umesh,k,s,max_deg,f,error_eval_refine_mesh,&
         error_gauss_nodes,error_gauss_weights)
end function L2NormEvalHOV1D

!!!!!!!!!!!!!! end added 11/07/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!  added Alex 05/24/2010

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SetErrTimeEvalArraysHOV1D
!
! This subroutine allocates the array where the error or the norm of the solution and the times 
! when these were evaluatd will be stored . 
! Also, it initiates the vectors of gaussian nodes and weights 
! Do not detouch from the main program 
! 
! This subroutine is highly dependent on the main program. It is created mainly 
! to organize the main program. 
! 
! The subroutine can not be called before the parameter num_eval_error are assigned values 
!
! (error_gauss_order) is set in the common_variables_mod !!
! 
! Most of the varaibles are looked up in the common_varaibles_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetErrTimeEvalArraysHOV1D
use common_variables_mod, only: error_gauss_order,num_eval_error,error_gauss_nodes,error_gauss_weights,err_time_sol
use gaussian_mod
!
integer :: loc_alloc_stat ! to keep allocation status
!
print *, "SetErrTimeEvalArraysHOV1D: Set error_gauss_order=", error_gauss_order
!!! Prepare the Gaussian nodes and weights that will be used in error evaluation
!!! Allocate arrays for gaussian nodes and weights to inegration in the second variable
  allocate (error_gauss_nodes(1:max(error_gauss_order,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetErrTimeEvalArraysHOV1D: Allocation error for variable (error_gauss_nodes)"
     end if 
     !
  allocate (error_gauss_weights(1:max(error_gauss_order,2)), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetErrTimeEvalArraysHOV1D: Allocation error for variable (error_gauss_weights)"
     end if 
     !
!!! Compute the gaussian nodes and weights for integration in the second variable    
  call GauLeg(Real(-1,DP),Real(1,DP),error_gauss_nodes,error_gauss_weights)  
!!!
! allocate space for error recording
allocate (err_time_sol(0:num_eval_error+1,3), stat=loc_alloc_stat)
     !
     if (loc_alloc_stat >0) then 
     print *, "SetErrTimeEvalArraysHOV1D: Allocation error for variable (err_sol)"
     end if 
     !
err_time_sol=0     
end subroutine SetErrTimeEvalArraysHOV1D

!!!!!!!! end added Alex 05/24/2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module misc_setup_mod