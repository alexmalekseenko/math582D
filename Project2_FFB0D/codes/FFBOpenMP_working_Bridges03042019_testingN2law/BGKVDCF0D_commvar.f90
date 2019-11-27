! 09/12/2014 Alex Alekseenko, CSUN, alexander.alekseenko@csun.edu
! BGKVDCF0D_commvar.f90 -- module containing definitions of global variables for the 
!						zero dimensional driver (spatially homogeneous solution) 
!						for kinetic equations. The variables will be accessed by the 
!						subroutines and modules of the driver. 
!                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BGKVDCF0D_commvar
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

!!!  Global Constants:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Selected Boundary conditions:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), parameter :: periodic_bc = 1    ! constant to denote the periodic boundary conditions 
integer (I4B), parameter :: exact_bc = 2       ! constant to denote the exact boundary conditions when the 
                                               ! distribution function is known exactly at the boundary
integer (I4B), parameter :: diffusive_reflection_bc = 3  ! constant to denote the diffusive reflective boundary conditions 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Global Variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Variables to describe zero's level of refinement (So far only refinement in the variable u will be 
!! used. All constants in x are the zero level constants. In U we will have levels of refinement. Parameters of 
!! the refinement will be build in, calculated automatically of solicited somehow separately. However, this file will 
!! define constants to create the mesh of level zero refinement. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer (I4B), dimension (:), allocatable :: k_b_list, k_c_list, N_list, d_c_list, d_b_list ! arrays of different values of k_c,k_b,s,M,N   

integer (I4B) :: k_c,k_b,N ! N is the number of physical cells, k_c is order of the highest Legendre polynomials 
						! used in the basis inside the cell  k_d is the number of Gauss-Lobatto nodes used on  
						! the t=const interface of the cell. 
integer (I4B) :: d_c,d_b ! d_c is the order of the Legendre basis in time variable used inside the cell (for the local solutions/test functions).
						! d_b is the number of Gauss-Lobatto nodes on the x=const cell interface
integer (I4B) :: max_deg_xt ! this is the max degree of the two-dimensional basis in x and t. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP), dimension (:), allocatable ::  x_gauss_n, x_gauss_w ! gauss nodes and weights for calculation of spectral coefficients or integrals in u 
real (DP), dimension (:), allocatable ::  t_gauss_n, t_gauss_w ! Gauss nodes and weights for calculation of spectral coefficients or integrals in t

!
real (DP), dimension (:), allocatable ::  x_lobat_n, x_lobat_w ! Gauss-Lobatto nodes and weights for calculation of spectral coefficients or integrals in u 
real (DP), dimension (:), allocatable ::  t_lobat_n, t_lobat_w ! Gauss-Lobatto nodes and weights for calculation of spectral coefficients or integrals in t
!! NOTE: orders of gauss lobatto quadratures in t and in x are controlled by d_b and k_b. 

real (DP) :: x_L, x_R ! left and right endpoints of the interval in X
real (DP), dimension (:), allocatable ::  xmesh_l, xmesh_r ! left and right endpoints of the cells in x (usually, the length of xmesh = N)
real (DP), dimension (:), allocatable ::  xgnodes ! Gauss nodal points for the cells in x. listed for all cells from 
real (DP), dimension (:), allocatable ::  xlnodes ! Gauss-Lobatto nodal points for the interfaces t=const in
 ! cells in x. 
real (DP) :: t_L, t_R ! left and right endpoints of the interval in T
real (DP) :: time ! current time
real (DP) :: dt ! the timestep
real (DP), dimension (:), allocatable ::  tlnodes  ! Gauss-Lobatto nodal points for the interfaces x=const in


logical :: mesh_x_uniform       ! mesh is uniform in x and u (true or false) (applied only to meshes of level 0)

!!!!!!
! Choice of BCs and the source function
integer (I4B) :: selected_rightbcond, selected_leftbcond ! variable to denote the choice of BCs
integer (I4B) :: selected_exact_sol ! variable to denote the choice of the exact solution (if available) 
!!!!!!
! The names of computed solution file and the directory to store it
character (len=50) :: current_solution_dir  ! name of the directory (relative to the local tree?) where the solution is stored
character (len=20) :: current_solution_base_name ! base name for the solution file -- more information about order, mesh, and time is added to it automatically
!!!!
integer (I4B) :: num_save_solution, num_eval_error ! how many records of solution from intial time to final time
                                                   ! how many times to caclulate the error
!!!!
! varaibles related to generation of non-unifom meshes: 

integer (I4B) :: x_nonuniform_mesh_type! parameter to describe non-uniform mesh in x and u 

! supported types of nonuminform mesh: 
! 1 --- the mesh is build based on the gauss nodes used for the integration of moments, moments_x_gauss_nodes and 
!       moments_u_gauss_nodes. intervals [x_left,x_right], [u_left,u_right] is divided in subintervals 
!       as to have gauss nodes at centerpoints. Some extra points need to be introduced to make this possible.
! 2 --- this mesh is to be used with the diffusive boundary conditions. The velocity of the wall (currently only u_{w}=0)
!       will be included in the mesh. Also, the cell near the walls will be 1/8 - 1/4 - 1/2        
! 3 --  for variable u -- this mesh will have small cells surrounding u=0 as prescribed by parameters
!       sml -- small cell levels amd smr -- small cell refinement factor 
! 3 --- This is a non-niform mesh in "x" with cells near wall be 1/4-1/2. Currenlty is not supported for meshes in "u"

! parameters related to the integration of moments in u and integration in x
!
integer (I4B) :: moments_x_gauss_order ! order of gausss quadrateures in x and u to calculated (total) moments
integer (I4B) :: moments_refine_x ! coefficients of mesh refinement in x and u for the evaluation of moments
real (DP), dimension (:), allocatable :: moments_x_gauss_nodes    ! arrays to keep gauss nodes for evaluation of moments 
real (DP), dimension (:), allocatable :: moments_x_gauss_weights! arrays to keep gauss weights for evaluation of moments  

real (DP), dimension (:,:), allocatable :: g_nds_all, g_wts_all ! arrays to store vaules of used gauss nodes and weights. 
                                                                 ! g_nds_all(i,j) i is the number of the node and j is the order of the gauss formula 
integer (I4B) :: g_max_order=10 ! maximum gauss order that will be used ...
                                           
! variables for temperatures on the walls -- will use in Diffusive BCs.
real (DP) :: T_w_left, T_w_right ! temperature of the left and right wall. Only used with the diffusive BCs

! to check whether restart is used.
character (len=20) :: restart_time_txt  ! time of last saved solution/restart time in text format
logical :: need_to_restart                ! tells if the simulations need to be restarted 
 
 
integer (I4B), dimension (:), allocatable :: I1_list ! list of basis functions for which A-array should be computed
                                     
!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable :: f ! the distribution function -- nodal values. 
real (DP), dimension (:), allocatable :: f1,f2,f3,f4,frhs1,frhs2,frhs3,frhs4,frhs5 ! the distribution function -- nodal values. 

!!!!!!!!!!!!!!!!!!!
integer (I2B) :: Num_OMP_threads ! number of threads for OpenMP. 

!!!!!!!!!!!!!!!!!!!
integer (I4B) :: rkmts ! The variable to keep the order of the Runge-Kutta method and the Adams-Bashworth MTS schemes

!!!!!!!!!! These variables are invoked to save the relative L1 error of the deviation of the distribution function from the local Maxwellian
real (DP), dimension (:), allocatable :: L1_a,L1_t ! L1_a keeps the error, L1_t keeps the time moments that correspond to the erro values
integer (I4B) :: L1_count ! This is the counter to keep track of the records
real (DP) :: L1_record_period,L1_next_time_record,L1_err_eval_period,L1_next_time_eval ! some scrap timing variables.

integer (I4B), dimension(:), allocatable :: run_mode_array ! -- this will be used to keep track of run_mode overtime 


!!!!!!!!!! These variables are invoked to save the values of the coefficients in the 
!!!! of te velocity dependent collision frequency
real (DP), dimension (:,:), allocatable :: Cco_rec ! first column keeps the time and the rest will keep the coefficients 
integer (I4B) :: Cco_count ! This is the counter to keep track of the records
real (DP) :: Cco_record_period,Cco_next_time_record,Cco_eval_period,Cco_next_time_eval ! some scrap timing variables.



end module BGKVDCF0D_commvar