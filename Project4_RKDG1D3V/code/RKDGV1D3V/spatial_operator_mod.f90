! This Module implements the spatial operator of the 1D High Order Velocity model (see notes) 
!
! The operations supported: 
! VolTerms  --- evaluates the volume integrals of the 1D HOV model (see notes)
! FluxTermsUpwind  ---  evaluate the flux terms of the 1d HOV models using Upwind flux(see notes)
! FluxTermsSymmetric  ---  evaluate the flux terms of the 1d HOV models using Symmetric Flux (see notes) 
!                 The flux terms are calling outside function responsible for boundary conditions.      
!!!!!!!!!!! 07/28/08 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
module spatial_operator_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
   implicit none
   
contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VolTerms_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
! 
! EXPECTS:
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,l,i,j) -- p is the index in the basis functions in x 
!                           -- l is the index in the basis functions in u
!                           -- i is the cell in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! s -- order of the polynomial interpolation in u
! M  --- umber of cells in u 
! N  --- umber of cells in x 
! 
! USES: (from common_varaibles_mod)
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    Sqp --- matrix of integrated basis functions in x. The exact formula for S is 
!             S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!                          
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function VolTerms_HOV1D3D (ftil, tau, k, M, N) result (f)
use common_variables_mod, only: xmesh, Sqp, q2plus1_vec

!!! Use DGVlib                                 
use DGV_commvar, only: nodes_u         
!!! 

                       
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,i,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- i is the nodes in u
                                   ! -- j is the cell in x
real (DP), intent (in) :: tau                                               ! time step parameter. Is used in time integration. 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- umber of cells in x 

! xmesh --- dimension (0:N) ! mesh points in varaible x 
! alphami --- dimension (0:s,M) ! coefficients of the evolution equations \alpha_{m_count,i}=alphami(m_count,i), m_count-basis function i - cell in u
! Sqp --- dimension (0:k,0:k)   ! matrix of products of basis functions in x: S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
! q2plus1_vec, dimension (0:k)             ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:k, M, N) :: f ! values of the volume integral term
!!! auxiliary variables: 
integer (I4B) :: i, j ! some local counters 
!!! 
intrinsic MatMul
!!!
do j=1,N
  do i=1,M
    f(:,i,j) = nodes_u(i)*q2plus1_vec*MatMul(Sqp,ftil(:,i,j))
  end do
  f(:,:,j)=(tau/(xmesh(j)-xmesh(j-1)))*f(:,:,j) 
end do 
!!!
end function VolTerms_HOV1D3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VolTerms_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
! 
! EXPECTS:
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,l,i,j) -- p is the index in the basis functions in x 
!                           -- l is the index in the basis functions in u
!                           -- i is the cell in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! s -- order of the polynomial interpolation in u
! M  --- umber of cells in u 
! N  --- umber of cells in x 
! 
! USES: (from common_varaibles_mod)
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    Sqp --- matrix of integrated basis functions in x. The exact formula for S is 
!             S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!                          
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function VolTermsOrds_HOV1D3D (ftil, tau, k, M, N, nodes_u) result (f)
use common_variables_mod, only: xmesh, Sqp, q2plus1_vec

                       
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,i,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- i is the nodes in u
                                   ! -- j is the cell in x
real (DP), intent (in) :: tau                                               ! time step parameter. Is used in time integration. 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- umber of cells in x 
real (DP), dimension (:), intent (in) :: nodes_u ! ovalues of ordinates for the components of the transport part
						! note that M should be equal to size(nodes_u)

! xmesh --- dimension (0:N) ! mesh points in varaible x 
! alphami --- dimension (0:s,M) ! coefficients of the evolution equations \alpha_{m_count,i}=alphami(m_count,i), m_count-basis function i - cell in u
! Sqp --- dimension (0:k,0:k)   ! matrix of products of basis functions in x: S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
! q2plus1_vec, dimension (0:k)             ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:k, M, N) :: f ! values of the volume integral term
!!! auxiliary variables: 
integer (I4B) :: i, j ! some local counters 
!!! 
intrinsic MatMul
!!!
do j=1,N
  do i=1,M
    f(:,i,j) = nodes_u(i)*q2plus1_vec*MatMul(Sqp,ftil(:,i,j))
  end do
  f(:,:,j)=(tau/(xmesh(j)-xmesh(j-1)))*f(:,:,j) 
end do 
!!!
end function VolTermsOrds_HOV1D3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VolTermsOMP_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
!
!!! 
! This is a copy of the above subroutines parallellized using OpenMP
!!!
!
! EXPECTS:
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,l,i,j) -- p is the index in the basis functions in x 
!                           -- l is the index in the basis functions in u
!                           -- i is the cell in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! s -- order of the polynomial interpolation in u
! M  --- umber of cells in u 
! N  --- umber of cells in x 
! 
! USES: (from common_varaibles_mod)
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    Sqp --- matrix of integrated basis functions in x. The exact formula for S is 
!             S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!                          
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function VolTermsOMP_HOV1D3D (ftil, tau, k, M, N) result (f)
use common_variables_mod, only: xmesh, Sqp, q2plus1_vec,Num_OMP_threads

!!! Use DGVlib                                 
use DGV_commvar, only: nodes_u         
!!! 

                       
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                   ! characteristic variables  ftil(p,m_count,i,j):
                                   ! -- p is the index in the basis functions in x 
                                   ! -- i is the nodes in u
                                   ! -- j is the cell in x
real (DP), intent (in) :: tau                                               ! time step parameter. Is used in time integration. 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- umber of cells in x 

! xmesh --- dimension (0:N) ! mesh points in varaible x 
! alphami --- dimension (0:s,M) ! coefficients of the evolution equations \alpha_{m_count,i}=alphami(m_count,i), m_count-basis function i - cell in u
! Sqp --- dimension (0:k,0:k)   ! matrix of products of basis functions in x: S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
! q2plus1_vec, dimension (0:k)             ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:size(ftil,1)-1,size(ftil,2),size(ftil,3)) :: f ! values of the volume integral term
!!! auxiliary variables: 
integer (I4B) :: i, j ! some local counters 
!!! 
intrinsic MatMul
!!!
! Open MP Directives 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OpenMP set the number of threads: 
!!!! !$OMP call omp_set_num_threads(Num_OMP_threads)

!$OMP PARALLEL DEFAULT (shared) PRIVATE(j,i) NUM_THREADS(Num_OMP_threads) 
!$OMP DO SCHEDULE(DYNAMIC, 1)
do j=1,N
! start OpenMP loop
  do i=1,M
    f(:,i,j) = nodes_u(i)*q2plus1_vec*MatMul(Sqp,ftil(:,i,j))
  end do
  f(:,:,j)=(tau/(xmesh(j)-xmesh(j-1)))*f(:,:,j) 
! end Open MP loop 
end do
!$OMP END DO 
!$OMP END PARALLEL 
!!!
end function VolTermsOMP_HOV1D3D


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FluxTermsUpwind_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
! 
! EXPECTS: 
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,i,j) -- p is the index in the basis functions in x 
!                           -- i is the index in the basis functions in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- number of cells in x 
!
! USES: (from the common_variables_mod) 
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!    altones  --- a useful vector of (1,-1,1,-1,..., (-1)^k)                      
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FluxTermsUpwind_HOV1D3D (ftil, TilFmiRight, TilFmiLeft, tau, k, M, N) result (f)
use common_variables_mod, only: xmesh, q2plus1_vec, altones

!!! Use DGVlib
use DGV_commvar, only: nodes_u
!!!
use spectral_tools_mod
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                 ! characteristic variables ftil(p,m_count,i,j):
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m_count is the index in the basis functions in u
                                 ! -- j is the cell in x
real (DP), dimension (:), intent (in) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{i}{t,x} 
							     ! at the right and left boundary points, TilFmiRight(m,i), TilFmiLeft(m,i)
                                 ! -- m is the index in the basis functions in u
                                 ! -- i is the cell in u
real (DP), intent (in) :: tau    ! time step parameter. Is used in time integration. 
! xmesh -- dimension (0:size(ftil,4)) ! mesh points in varaible x 
! alphami -- dimension (0:size(ftil,2)-1 , size(ftil,3)) ! coefficients of the evolution equations \alpha_{m,i}=alphami(m,i), m-basis function i - cell in u
! q2plus1_vec, dimension (0:size(ftil,1)-1)  ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
! altones, dimension (0:size(ftil,1)-1)! convenient vector of (1,-1,1,\ldots,(-1)^k) will use in matrix operations 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodeal points in u 
! N  --- umber of cells in x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:k, M, N) :: f ! values of the flux term
!!! auxiliary variables: 
integer (I4B) :: m_count, i, j ! some local counters 
intrinsic MatMul, Sum
!!!
do i=1,M
   ! First we do the interior elements: 
   do j=2,N-1 
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j-1),xmesh(j),ftil(:,i,j)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-2),xmesh(j-1),ftil(:,i,j-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j),xmesh(j+1),ftil(:,i,j+1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-1),xmesh(j),ftil(:,i,j))*altones ) ! left endpoint of the cell
   else 
f(:,i,j) = 0                  
   end if
f(:,i,j) = (tau/(xmesh(j)-xmesh(j-1)))*q2plus1_vec*nodes_u(i)*f(:,i,j)
   end do  !! end loop in j
   ! Now let us take care of the boundary elements 
   ! this will involve the boundary conditions
   ! j=1 --> the left boundary 
   !!!
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,1) = (EvCellLeg1D(xmesh(1),xmesh(0),xmesh(1),ftil(:,i,1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - TilFmiLeft(i)*altones)                           ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,1)= (EvCellLeg1D(xmesh(1),xmesh(1),xmesh(2),ftil(:,i,2)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(0),xmesh(0),xmesh(1),ftil(:,i,1))*altones ) ! left endpoint of the cell
   else 
f(:,i,1) = 0                  
   end if
f(:,i,1) = (tau/(xmesh(1)-xmesh(0)))*q2plus1_vec*nodes_u(i)*f(:,i,1)
   ! end left boundary 
   ! j=N --> the right boundary
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,N) = (EvCellLeg1D(xmesh(N),xmesh(N-1),xmesh(N),ftil(:,i,N)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-2),xmesh(N-1),ftil(:,i,N-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,N) = (TilFmiRight(i) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-1),xmesh(N),ftil(:,i,N))*altones ) ! left endpoint of the cell
   else 
f(:,i,N) = 0                  
   end if
f(:,i,N) = (tau/(xmesh(N)-xmesh(N-1)))*q2plus1_vec*nodes_u(i)*f(:,i,N)
   ! end right boundary 
end do !end loop in i
end function FluxTermsUpwind_HOV1D3D 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FluxTermsUpwindOrds_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
! 
! This is a copy of the above subroutine except the ordinates used in the flux are supplied to the subroutine as a variable.
! 
! EXPECTS: 
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,i,j) -- p is the index in the basis functions in x 
!                           -- i is the index in the basis functions in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- number of cells in x 
!
! USES: (from the common_variables_mod) 
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!    altones  --- a useful vector of (1,-1,1,-1,..., (-1)^k)                      
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FluxTermsUpwindOrds_HOV1D3D (ftil, TilFmiRight, TilFmiLeft, tau, k, M, N, nodes_u) result (f)
use common_variables_mod, only: xmesh, q2plus1_vec, altones

use spectral_tools_mod
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                 ! characteristic variables ftil(p,m_count,i,j):
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m_count is the index in the basis functions in u
                                 ! -- j is the cell in x
real (DP), dimension (:), intent (in) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{i}{t,x} 
							     ! at the right and left boundary points, TilFmiRight(m,i), TilFmiLeft(m,i)
                                 ! -- m is the index in the basis functions in u
                                 ! -- i is the cell in u
real (DP), intent (in) :: tau    ! time step parameter. Is used in time integration. 
! xmesh -- dimension (0:size(ftil,4)) ! mesh points in varaible x 
! alphami -- dimension (0:size(ftil,2)-1 , size(ftil,3)) ! coefficients of the evolution equations \alpha_{m,i}=alphami(m,i), m-basis function i - cell in u
! q2plus1_vec, dimension (0:size(ftil,1)-1)  ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
! altones, dimension (0:size(ftil,1)-1)! convenient vector of (1,-1,1,\ldots,(-1)^k) will use in matrix operations 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodeal points in u 
! N  --- umber of cells in x 
real (DP), dimension (:), intent (in) :: nodes_u ! ordinates are supplied to the subroutine. Check that M=size(nodes_u,1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:k, M, N) :: f ! values of the flux term
!!! auxiliary variables: 
integer (I4B) :: m_count, i, j ! some local counters 
intrinsic MatMul, Sum
!!!
do i=1,M
   ! First we do the interior elements: 
   do j=2,N-1 
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j-1),xmesh(j),ftil(:,i,j)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-2),xmesh(j-1),ftil(:,i,j-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j),xmesh(j+1),ftil(:,i,j+1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-1),xmesh(j),ftil(:,i,j))*altones ) ! left endpoint of the cell
   else 
f(:,i,j) = 0                  
   end if
f(:,i,j) = (tau/(xmesh(j)-xmesh(j-1)))*q2plus1_vec*nodes_u(i)*f(:,i,j)
   end do  !! end loop in j
   ! Now let us take care of the boundary elements 
   ! this will involve the boundary conditions
   ! j=1 --> the left boundary 
   !!!
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,1) = (EvCellLeg1D(xmesh(1),xmesh(0),xmesh(1),ftil(:,i,1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - TilFmiLeft(i)*altones)                           ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,1)= (EvCellLeg1D(xmesh(1),xmesh(1),xmesh(2),ftil(:,i,2)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(0),xmesh(0),xmesh(1),ftil(:,i,1))*altones ) ! left endpoint of the cell
   else 
f(:,i,1) = 0                  
   end if
f(:,i,1) = (tau/(xmesh(1)-xmesh(0)))*q2plus1_vec*nodes_u(i)*f(:,i,1)
   ! end left boundary 
   ! j=N --> the right boundary
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,N) = (EvCellLeg1D(xmesh(N),xmesh(N-1),xmesh(N),ftil(:,i,N)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-2),xmesh(N-1),ftil(:,i,N-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,N) = (TilFmiRight(i) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-1),xmesh(N),ftil(:,i,N))*altones ) ! left endpoint of the cell
   else 
f(:,i,N) = 0                  
   end if
f(:,i,N) = (tau/(xmesh(N)-xmesh(N-1)))*q2plus1_vec*nodes_u(i)*f(:,i,N)
   ! end right boundary 
end do !end loop in i
end function FluxTermsUpwindOrds_HOV1D3D 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FluxTermsUpwindOMP_HOV1D3D() 
! Evaluates the volume terms of the spatial operator in the 1D HOV model (see notes) 
! 
!
!!! 
! This is a copy of the above subroutines parallellized using OpenMP
!!!
! EXPECTS: 
!    ftil --- the main variable --- coefficienst of spactral decomposition in characteristic basis. 
!             ftil(p,i,j) -- p is the index in the basis functions in x 
!                           -- i is the index in the basis functions in u
!                           -- j is the cell in x
!    tau --- time step parameter. This parameter can be used if the spatial operator is used in the time integration. 
!            by setting tau to the actual time step will result in multiplying the volume integral terms by this number. 
!            The tau parameter is introduced to off set the division by $\Delta x_{j}$time integration. 
!            To turn this feature off, set tau=1. This will result in multiplication of the volume terms by 1, which "changes nothing". 
! k -- order of the polynomial interpolation in x 
! M  --- number of nodes in u 
! N  --- number of cells in x 
!
! USES: (from the common_variables_mod) 
!    xmesh -- meshpoints in variable x
!    alphami -- array of coefficents (the characteristic speeds) in the evolution equations
!    q2plus1_vec --- useful (constant) vector, used in the matrix multiplication (1,2,...,2q+1,...,2k+1) 
!    altones  --- a useful vector of (1,-1,1,-1,..., (-1)^k)                      
! RETURNS: 
!    f --- array of the volume terms of the same size as ftil 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function FluxTermsUpwindOMP_HOV1D3D (ftil, TilFmiRight, TilFmiLeft, tau, k, M, N) result (f)
use common_variables_mod, only: xmesh, q2plus1_vec, altones, Num_OMP_threads

!!! Use DGVlib
use DGV_commvar, only: nodes_u
!!!
use spectral_tools_mod
!!! main variables
real (DP), dimension (0:,:,:), intent (in) :: ftil ! main variable -- coefficients of spectral decomposition in 
                                 ! characteristic variables ftil(p,m_count,i,j):
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m_count is the index in the basis functions in u
                                 ! -- j is the cell in x
real (DP), dimension (:), intent (in) :: TilFmiRight, TilFmiLeft ! values of the \tilde{f}_{i}{t,x} 
							     ! at the right and left boundary points, TilFmiRight(m,i), TilFmiLeft(m,i)
                                 ! -- m is the index in the basis functions in u
                                 ! -- i is the cell in u
real (DP), intent (in) :: tau    ! time step parameter. Is used in time integration. 
! xmesh -- dimension (0:size(ftil,4)) ! mesh points in varaible x 
! alphami -- dimension (0:size(ftil,2)-1 , size(ftil,3)) ! coefficients of the evolution equations \alpha_{m,i}=alphami(m,i), m-basis function i - cell in u
! q2plus1_vec, dimension (0:size(ftil,1)-1)  ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
! altones, dimension (0:size(ftil,1)-1)! convenient vector of (1,-1,1,\ldots,(-1)^k) will use in matrix operations 
integer (I4B), intent (in) :: k,M,N
! k -- order of the polynomial interpolation in x 
! M  --- number of nodeal points in u 
! N  --- umber of cells in x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
real (DP), dimension (0:size(ftil,1)-1,size(ftil,2),size(ftil,3)) :: f ! values of the flux term
!!! auxiliary variables: 
integer (I4B) :: m_count, i, j ! some local counters 
intrinsic MatMul, Sum

!!!
! Open MP Directives 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! OpenMP set the number of threads: 
!!!! !$OMP call omp_set_num_threads(Num_OMP_threads)

!$OMP PARALLEL DEFAULT (shared) PRIVATE(i,j) NUM_THREADS(Num_OMP_threads)
!$OMP DO SCHEDULE(DYNAMIC,128) 
!!!
do i=1,M
! begin Open MP loop 
   ! First we do the interior elements: 
   do j=2,N-1 
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j-1),xmesh(j),ftil(:,i,j)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-2),xmesh(j-1),ftil(:,i,j-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,j)= (EvCellLeg1D(xmesh(j),xmesh(j),xmesh(j+1),ftil(:,i,j+1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(j-1),xmesh(j-1),xmesh(j),ftil(:,i,j))*altones ) ! left endpoint of the cell
   else 
f(:,i,j) = 0                  
   end if
f(:,i,j) = (tau/(xmesh(j)-xmesh(j-1)))*q2plus1_vec*nodes_u(i)*f(:,i,j)
   end do  !! end loop in j
   ! Now let us take care of the boundary elements 
   ! this will involve the boundary conditions
   ! j=1 --> the left boundary 
   !!!
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,1) = (EvCellLeg1D(xmesh(1),xmesh(0),xmesh(1),ftil(:,i,1)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - TilFmiLeft(i)*altones)                           ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,1)= (EvCellLeg1D(xmesh(1),xmesh(1),xmesh(2),ftil(:,i,2)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(0),xmesh(0),xmesh(1),ftil(:,i,1))*altones ) ! left endpoint of the cell
   else 
f(:,i,1) = 0                  
   end if
f(:,i,1) = (tau/(xmesh(1)-xmesh(0)))*q2plus1_vec*nodes_u(i)*f(:,i,1)
   ! end left boundary 
   ! j=N --> the right boundary
   if (nodes_u(i) > 0) then               ! makes sense to replace 0 with a small number....
f(:,i,N) = (EvCellLeg1D(xmesh(N),xmesh(N-1),xmesh(N),ftil(:,i,N)) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-2),xmesh(N-1),ftil(:,i,N-1))*altones ) ! left endpoint of the cell
   else if (nodes_u(i) < 0) then          ! makes sense to replace 0 with a small number....
f(:,i,N) = (TilFmiRight(i) &    ! multiplication by ones is omitted. altones are the values of basis functions at the 
             - EvCellLeg1D(xmesh(N-1),xmesh(N-1),xmesh(N),ftil(:,i,N))*altones ) ! left endpoint of the cell
   else 
f(:,i,N) = 0                  
   end if
f(:,i,N) = (tau/(xmesh(N)-xmesh(N-1)))*q2plus1_vec*nodes_u(i)*f(:,i,N)
   ! end right boundary 
! end Open MP loop
end do !end loop in i
!$OMP END DO 
!$OMP END PARALLEL

end function FluxTermsUpwindOMP_HOV1D3D 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SpatialOperatorHOV1D3D_UpwindFlux
! 
! This function  implement the spatial differential operator of the DG 1D HOV model 
! with the upwind flux. 
! 
! Do not detouch form the main program
! 
! Do not call before k,s,M,N,Ainvml,u_gauss_nodes,u_wleg_pol have been initialized
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SpatialOperatorHOV1D3D_UpwindFlux(f1,ftil1,curr_time,tau) 

use common_variables_mod, only: selected_leftbcond,selected_rightbcond,periodic_bc,exact_bc,diffusive_reflection_bc,& 
                                xmesh,k,N,selected_exact_sol
                      
!!!! DGVlib                                
use DGV_sf02
use DGV_commvar, only: nodes_u
!!!!!!!!!!!!
                                
use boundary_cond_mod
use right_side_mod
! main variables 

real (DP), dimension (0:,:,:), intent (out) :: f1 ! values of the discrete transport operator 
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m is the index in the basis functions in u
                                 ! -- j is the cell in x
real (DP), dimension (0:,:,:), intent (in) :: ftil1 ! main variable main variable -- coefficients of spectral decomposition in 
                                                      !  ftil(p,m,j):
                                                      ! -- p (0:k) is the index in the basis functions in x 
                                                      ! -- m (1:Mu*su*Mv*sv*Mw*sw) is the index in the basis functions in u/velocity nodes - comes from DGVlib
                                                      ! -- j is the cell in x
                                                      
real (DP), intent (in) :: curr_time ! the coordinate time
real (DP), intent (in) :: tau ! the time step = dt                                                                  
! auxiliary variables 
real (DP), dimension (:), allocatable  :: TilFmi1Left   ! values of \tilde{f}_{m,i}(t,x) at left boundary
real (DP), dimension (:), allocatable  :: TilFmi1Right  ! values of \tilde{f}_{m,i}(t,x) at right boundary
                                ! at the right and left boundary points, TilFmiRight(m), TilFmiLeft(m)
                                ! -- m is the index in the basis functions/nodes  in u
                                
real (DP), dimension (:,:,:), allocatable :: r1!  temp varaibles to calculate the contribution from the right side
                                 ! (p,m,j)
                                 ! -- p is the index in the basis functions in x 
                                 ! -- m is the index in the basis functions in u
                                 ! -- j is the cell in x
integer :: loc_alloc_stat ! to keep allocation status                                 

!!!!!!!!!!!!!!
! k -- order of the polynomial interpolation in x 
! N  --- umber of cells in x 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! we use constants from common_variables_mod to select between different available types of boundary 
! conditions and different boundary data (if available)
!!!!!!!!!!!
!! Evaluation of the spatial operator:
!! ADD VOLUME TERMS:
f1 = VolTermsOMP_HOV1D3D (ftil1,tau,k,size(ftil1,2),N)
!! ADD FLUX TERMS INCLUDING BCs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Modified Alex 05/15/09  !!!!!!!!!!!!!!!!!
allocate (TilFmi1Left(1:size(ftil1,2)),TilFmi1Right(1:size(ftil1,2)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "SpatialOperatorHOV1D3D_UpwindFlux: Allocation error for variable  (TilFmi1Left,TilFmi1Right)"
 stop
end if
! First, we need to compute the values of the function at the boudaries
if (selected_leftbcond == selected_rightbcond) then 
! if the boundary conditions are the same type on both boundaries 
select case (selected_leftbcond)
 case(periodic_bc) 
     call PeriodicBCsLeg(TilFmi1Left,TilFmi1Right,xmesh,ftil1)
 case(exact_bc)  
   select case (selected_exact_sol)
       case(1)
     call ExactBCs1D3DLeg (TilFmi1Left,TilFmi1Right,xmesh,curr_time,f_1D3D_u_vectors) ! the same as case 2
       case(2)
     call ExactBCs1D3DLeg (TilFmi1Left,TilFmi1Right,xmesh,curr_time,f_1D3D_u_vectors) ! the same as case 1
       case default 
	 print *, "Unsupported type of exact solution"
	 TilFmi1Left = 0
	 TilFmi1Right = 0
   end select
case (diffusive_reflection_bc) 
    print *, "Diffusion Boundary conditions has not been implemented yet. Stop"
    stop
    !call DiffuiveBCsLeftHOV1D(TilFmi1Left,ftil1(:,:,1)) 
    ! call DiffuiveBCsRightHOV1D(TilFmi1Right,ftil1(:,:,N)) 
case default 
		    print *, "BoundaryConditionsTilFmiLeft: Unsupported Type of Boundary Conditions"
		    stop
end select 
else 
print *, "Different BCs on opposite boundaries have not been implemented yet!"
		    stop
end if     
! Now when the solution at the boundary is computed, we can use it in the numerical flux... 
f1 = f1 - FluxTermsUpwindOMP_HOV1D3D (ftil1,TilFmi1Right,TilFmi1Left,tau,k,size(nodes_u,1),N) 
deallocate (TilFmi1Right,TilFmi1Left)
!! ADD THE RIGHT SIDE (COLLISION TERMS)
!! COMMENTED with !* IF FREE MOLECULAR FLOW IS DESIRED (Alex, 05/19/2010) 
! allocate space for the right side values
allocate (r1(0:size(ftil1,1)-1, size(ftil1,2), size(ftil1,3)), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then 
   print *, "SpatialOperatorHOV1D3D_UpwindFlux: Allocation error for variable (r1)"
   end if 
! evaluate the right side   
call RightSide1D3DOMP_HOV1D (r1,ftil1,tau,curr_time)
f1 = f1 + r1
! deallocate temporary variables r1
deallocate (r1)
!!! End Spatial Operator
end subroutine SpatialOperatorHOV1D3D_UpwindFlux

end module spatial_operator_mod 