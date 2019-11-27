! 09/05/08 Alex  
!
!  common_variables_mod.f90
!  
!  This module contains common variables for the 1D HOV model. 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module common_variables_mod
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none
!!!  Global Constants:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Selected Boundary conditions:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), parameter :: periodic_bc = 1    ! constant to denote the periodic boundary conditions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! added 09/07/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), parameter :: exact_bc = 2       ! constant to denote the exact boundary conditions when the 
                                               ! distribution function is known exactly at the boundary
integer (I4B), parameter :: diffusive_reflection_bc = 3  ! constant to denote the diffusive reflective boundary conditions 
!!!!!!!!!!!! end added 09/07/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! Main Variables (the components of the solution, mesh, time and such)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: k_list, s_list, M_list, N_list ! arrays of different values of k,s,M,N   
integer (I4B) :: k,s,M,N, max_deg, rk ! order of approximation and number of cells:
                         ! k -- order of polynomial approximation in x
                         ! s -- order of polynomial approximation in u/number of nodes in nodal DG
                         ! max_deg -- order of combined polynomial approximation 
                         ! M -- number of mesh cells in u
                         ! N -- number of mesh cells in Runge Kutta method
                         ! rk -- order of the Runge Kutta method 
integer (I4B) :: num_save_solution, num_eval_error ! how many records of solution from intial time to final time
                                                   ! how many times to caclulate the error
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! added 09/05/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable ::  x_gauss_nodes, x_gauss_weights ! nodes and weights for calculation of spectral 
                          ! decomposition of the solution/initial data in x
real (DP), dimension (:), allocatable ::  u_gauss_nodes, u_gauss_weights ! nodes and weights for calculation of spectral 
                          ! decomposition of the solution/initial data and bounday data in u
real (DP), dimension (:,:), allocatable :: x_wleg_pol,u_wleg_pol !  arrays to store values of leg_polynomials on 
                          ! gaussian nodes (x_gauss_nodes, u_gauss_nodes) multiplied by the corresponding weights
!!!!!!!!!!!! end added 09/05/08 Alex !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                  
integer (I4B) :: error_eval_refine_mesh=2, error_gauss_order=8 ! coefficient of mesh refinement for error evaluation -- each cell 
                         ! [umesh(i-1),umesh(i)]\times[xmesh(j-1),xmesh(j)] will be refined 
                         ! into (error_eval_refine_mesh)^2 parts
                         ! error_gauss_order is the order of the gauss formula used for error evaluation 
real (DP), dimension (:), allocatable ::  error_gauss_nodes,error_gauss_weights ! nodes and weights for error calculation                          
real (DP), dimension (:,:), allocatable :: err_time_sol ! array to store the error (or the norm) of the solution (whent the exact solution is not available, we will compute norm instead                       


real (DP), dimension (:,:,:), allocatable :: ftil ! \tilde{f}^_{p,j;m}(t)
                                                  ! main variable -- coefficients of spectral decomposition in 
                                                      ! in the nodal dg methd the variables are autopmatically characteristic. characteristic  ftil(p,m,j):
                                                      ! -- p (0:k) is the index in the basis functions in x 
                                                      ! -- m (1:Mu*su*Mv*sv*Mw*sw) is the index in the basis functions in u/velocity nodes - comes from DGVlib
                                                      ! -- j is the cell in x
real (DP), dimension (:,:,:), allocatable :: ftil1,ftil2,ftil3,ftil4,frhs1,frhs2,frhs3,frhs4,frhs5 ! more of the same to be used in multiple time stepping
                                                      
real (DP) :: curr_time, initial_time, final_time ! the coordinate time variables  
real (DP) :: dt ! \Delta t -- the step in time
real (DP) :: cfl ! the cfl number = dt/dx
real (DP) :: x_left, x_right, u_left, u_right   !  boundaries of the intervals in x and u     
logical :: mesh_x_uniform, mesh_u_uniform       ! mesh is uniform in x and u (true or false) 
real (DP), dimension (:), allocatable :: xmesh  ! mesh points in variable x 
real (DP), dimension (:), allocatable :: umesh  ! mesh points in variable u 
! characteristic speeds and related things 
real (DP), dimension (:,:), allocatable :: alphami ! coefficients of the evolution equations \alpha_{m,i}=alpha(m,i), m-basis function i - cell in u
real (DP), dimension (:,:,:), allocatable :: Aml,Ainvml 
             ! alphami(m,i) =\alpha_{m,i} is the array of eigenvalues 
             ! and Aml(:,m,i)=A_{ml,i} is the array of eignevectors of of \hat{T}_{pl,i}=\int_{U_{i}}u\lambda_{l,i}(u)\lambda_{p,i}(u)
             !          i -- is the interval in velocity U_{i}
             !          m -- is the eigenvalue number
             ! Ainvml(:,:,i) is the matrix inverse of Aml(:,:,i)
! Choice of BCs and the sourse function
integer (I4B) :: selected_rightbcond, selected_leftbcond ! variable to denote the choice of BCs
integer (I4B) :: selected_exact_sol ! variable to denote the choice of the exact solution (if available) 

!!! Auxiliary variables (usefull arrays and other things that will not change during the 
!!!  execution of the program) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! integrals of products of basis functions 
real (DP), dimension (:,:), allocatable :: Sqp   ! matrix of products of basis functions in x: S_{qp}=\int_{-1}^{1} P_{p}(y) P'_{q}(y) dy
! values of the basis fuctions

! misc coefficients 
integer (I4B), dimension (:), allocatable  :: q2plus1_vec ! a convenient vector of cofficients (1,3,...,2q+1,...,2k+1) will use in matrix operations 
integer (I4B), dimension (:), allocatable  :: altones     ! a convenient vector of alternating 1 and -1: (1,-1,...,(-1)^q,...,(-1)^k) will use in matrix operations 
!
character (len=50) :: current_solution_dir  ! name of the directory (relative to the local tree?) where the solution is stored
character (len=20) :: current_solution_base_name ! base name for the solution file -- more information about order, mesh, and time is added to it automatically
!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! Added Alex 10/15/08              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! variables and arrays used for the integration of moments

integer (I4B) :: moments_x_gauss_order, moments_u_gauss_order ! order of gausss quadrateures in x and u to calculated (total) moments
integer (I4B) :: moments_refine_x, moments_refine_u ! coefficients of mesh refinement in x and u for the evaluation of moments
real (DP), dimension (:), allocatable :: moments_x_gauss_nodes, moments_u_gauss_nodes     ! arrays to keep gauss nodes for evaluation of moments 
real (DP), dimension (:), allocatable :: moments_x_gauss_weights, moments_u_gauss_weights ! arrays to keep gauss weights for evaluation of moments  
real (DP), dimension (:,:,:), allocatable :: moments_time_sol, moments_time_exact ! errays to keep evaluated total moments in the spectral solution and in the exact solution 
														  ! moments_sol (i,j,k) 
                                                          ! i --- ## of moment
                                                          ! j = 0,1,...  --- # of the moment of saving. 
                                                          ! k = 1,2 --- for each component of the solution
! varaibles related to generation of non-unifom meshes: 

integer (I4B) :: x_nonuniform_mesh_type, u_nonuniform_mesh_type ! parameter to describe non-uniform mesh in x and u 

! supported types of nonuminform mesh: 
! 1 --- mesh is build based on the gauss nodes used for the integration of moments, moments_x_gauss_nodes and 
!       moments_u_gauss_nodes. intervals [x_left,x_right], [u_left,u_right] is divided in subintervals 
!       as to have gauss nodes at centerpoints. Some extra points need to be introduced to make this possible.
!     

real (DP), dimension (:,:,:), allocatable :: rec_dens, rec_avel, rec_tempr ! arrays to keep the record of the macroparameters (in the spectral represetation?)
!!!!!!!!!!!! end added Alex 10/15/08           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!! added Alex 12/25/08 !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP) :: gasR,gasTref,gasalpha,gasmuref ! is the normal gas constant, gas reference temperature, gas reference viscosity and alpha constant
                                           ! --- need to be read from the parameter.dat

!!!! THESE ARE supplementary arrays, so we do not have to recalculate is over and over.... 
real (DP), dimension (:,:), allocatable :: Xnodes,Unodes ! to keep the nodes on the refined cells in x and u
                                                         ! Xnodes(p,j) p -- node, j -- cell in x
                                                         ! Unodes(m,i) m -- node, i -- cell in u 
real (DP), dimension (:,:), allocatable :: XLegPol, ULegPol, XWeightLegPol, UWeightLegPol ! to keep the values of Leg. polynomials on refined cells
									! XWeightLegPol(p,l) p -- nodes, l -- number of polynomial
									! UWeightLegPol(m,i) m -- nodes, i -- number of polynomial 
real (DP), dimension (:), allocatable :: Xnodes_unit, Unodes_unit ! to keep the refined nodes on the "canonical" cell [-1,1]
                                    ! Xnodes_unit(p) p -- nodes
                                    ! Unodes_unit(m) m -- nodes
real (DP), dimension (:), allocatable :: Xweights_unit,Uweights_unit ! to keep the refined weights on the "canonical" cell [-1,1]
                                    ! Xweights_unit(p) p -- nodes
                                    ! Uweights_unit(m) m -- nodes
                                    
!!! End added Alex 12/25/08 !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Added Alex 05/22/09
real (DP) :: T_w_left, T_w_right ! temperature of the left and right wall. Only used with the diffusive BCs
!!! End Added Alex 05/22/09 !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Added Alex 08/27/10

character (len=20) :: restart_time_txt  ! time of last saved solution/restart time in text format
logical :: need_to_restart                ! tells if the simulations need to be restarted 


!!! End Added Alex 08/27/10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!
integer (I2B) :: Num_OMP_threads ! number of threads for OpenMP.


 
end module common_variables_mod
