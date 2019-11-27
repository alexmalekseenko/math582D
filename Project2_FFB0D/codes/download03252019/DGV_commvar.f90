! 11/12/2014 Alex Alekseenko, CSUN, alexander.alekseenko@csun.edu
!            Credit: Craig Euler, former graduate student at CSUN.
! DGV_commvar.f90 -- module containing definitions of global variables for the library for collision operator. 
!						subroutines and modules of the library will access these variables. 
!						also, subroutines that use the library will access these variables. 
!                       
!!!!!!!!!!!!!!!!!!!!!!!!!!!

module DGV_commvar
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Global Variables:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! Variables to describe zero's level of refinement (The code has subroutines and the capacity to 
!! implement octree refinement in the velocity variable. Parameters of 
!! the refinement will be build in, calculated automatically of solicited somehow separately. 
!! However, the following constants will define the original mesh, that is, mesh on the level zero refinement. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

integer (I4B), dimension (:), allocatable :: su_list, sv_list, sw_list, Mu_list, Mv_list, Mw_list ! arrays 
			! of different values of s,M for each velocity direction. Here s stands for the number of Gauss 
			! interpolation nodes on each interval and M is the number of intervals. 
			! Originally, the program will solisit lists (one, two, three, or more values (up to 20))
			! the lists are then used if several runs with different resolutions are desired. 
			! presently only the first values is used in the driver. 

integer (I4B) :: su,sv,sw,Mu,Mv,Mw ! M is the number of velocity cells of level 0, s is the number of gauss points per level zero cell.
            ! In the current implementation of the library the level zero mesh is uniform. These variables 
            ! keep the parameters of the zero level mesh.

real (DP), dimension (:), allocatable ::  u_gauss_n, u_gauss_w ! Gauss nodes and weights for calculation of spectral coefficients or integrals in x
			! these arrays will keep the values of Gauss nodes and weights that can be used to evaluate gauss integrals on velcotiy cells. 
			! by default, u_gauss_n and u_gauss_w have the dimension s and contains nodes that are used to construct nodes of the velcotiy discretization 
			! on the zero mesh. 

real (DP) :: u_L, u_R, v_L, v_R, w_L, w_R ! left and right endpoints of the interval in U,V,and W 
            ! These variables are used to set up the boundaries of the velocity domain.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The next group of variables enables the velocity cells and the velocity nodal points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable ::   grids_u, grids_v, grids_w 
integer (I4B), dimension (:), allocatable :: grids_cap_u,grids_cap_v,grids_cap_w
     ! grids_cap_u, grids_cap_v and _w contain the total number of 1D meshpoints (including both enpoints, e.g., one cell has two points 
     ! two cells have three and so on. recall that the 3D grid is a cartesian product of the 1d meshes.      
     ! grids_u, grids_v, and grids_w contain 1D meshes for each of the grids. Meshes are following each other, 
     ! to find the appropriate records use the grid capacity array "grids_cap_u/_v/_w
real (DP), dimension (:), allocatable :: cells_lu, cells_lv, cells_lw, cells_ru, cells_rv, cells_rw
integer (I4B), dimension (:), allocatable :: cells_pgrid, cells_cgrid, cells_refu, cells_refv, &
                                      cells_refw, cells_gow, cells_gou, cells_gov
integer (I4B), dimension (:), allocatable :: cells_ugi, cells_vgi, cells_wgi ! relative addresses on the parent grid, currently only works in the case of a single uniform grid .. ... 
     ! cells are the main object of the discretization, they form an ierarchical grid. 
     ! cells_pgrid is the number of the parent grid. 
     ! cells_cgrid is -1 if the cell is not refined or is equal to the number of the child grid. 
     ! cells_lu/_lv/_lw and cells_ru/_rv/_rw are the corners of the cells
     ! cells_refu/_refv/_refw are the refinement coefficients 
     ! cells_gou/_gow /_gov Cell's Gauss order. --- so far will be the same in both u, v, and w 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
real (DP), dimension (:), allocatable ::  nodes_u, nodes_v, nodes_w, nodes_gwts
real (DP), dimension (:), allocatable ::  nodes_ut, nodes_vt, nodes_wt ! u - ubar, v - vbar, w - wbar. t in name comes from tilde. updated when the moments are updated
integer (I4B), dimension (:), allocatable :: nodes_pcell, nodes_ui, nodes_vi, nodes_wi
     ! velocity nodes of the nodal-DG discretization. nodes are listed for each cell in the order cells are listed
     ! nodes_pcell is the number of the parent cell
     ! nodes_u/_v/_w are the coordinates of the node. 
     ! nodes_gwts are the values of products of weights of gaussian quadrature formulas to be used in 3D integration
     ! nodes_ui/_vi/_wi are the values of the local indices that determine addresses within the cell of the 1D nodal 
     ! points in each variable associated with the node

logical :: mesh_u_uniform, mesh_v_uniform, mesh_w_uniform       ! mesh is uniform in x and u (true or false) (applied only to meshes of level 0)
	 ! this is the variable that determines if the zero level mesh is uniform. In the current implementation the value is 
	 ! always true. ATTENTION: used both for primary and for secondary mesh (see blow)		

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These variables are complimentary and are udes in the evaluation of the collision integral by the 
! method of Alekseenko-Josyula. They depend on the velocity nodes and cells and also on the pre-computed operator A
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 integer (I4B), dimension (:), allocatable :: nodes_Ashift,nodes_dui,nodes_dvi,nodes_dwi,nodes_phican    
     ! this array contains the shift in the A index to be used for periodic grids/self similar basis functions
     ! nodes_dui/_dvi/_dwi contains the shift in cell idex for each node. = Shift in the index to go from the cell that has the node to the canonical cell. 
     ! nodes_phican contains the number of the canonical node on the canonical cell that corresponds to the given node.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These are variables of the secondary velocity grid. The secondary grid is currently uniform. 
! The seocondary grid is sufficiently functional -- separate copies of subroutines will be 
! created to work with the secondary grid. 
! The defauls suffix that is dedcated to the secondary grid is II. For examples, nodes_uII designates a
! variable of the secondary grid. 
!  
! Secondary grid will use the same boundaries as the primary grid. Therefore the same u_L, u_R, v_L, v_R, w_L, w_R
!
logical :: SecondaryVelGridInUse	= .false.		! when this flag is true, the secondary 
                                                ! velocity grid is used (to evaluate Boltzmann collision integral on 
                                                ! velocity mesh other then the primary 
                                 !Presently, the parameters of the first (primary) velocity grid and the second (secondary) velocity grid are both taken from the arrays: 
                                 !\texttt{su\_list, sv\_list, sw\_list, Mu\_list, Mv\_list, Mw\_list}. The first element of these arrays gives the parameters fo the 
                                 !primary grid, the second element gives the parameters of the secondary grid. 
!
integer (I4B) :: suII,svII,swII,MuII,MvII,MwII ! M is the number of velocity cells of level 0, s is the number of gauss points per level zero cell.
            ! In the current implementation of the library the level zero mesh is uniform. These variables 
            ! keep the parameters of the zero level mesh.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The next group of variables enables the velocity cells and the velocity nodal points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable ::   grids_uII, grids_vII, grids_wII 
integer (I4B), dimension (:), allocatable :: grids_cap_uII,grids_cap_vII,grids_cap_wII
     ! grids_cap_u, grids_cap_v and _w contain the total number of 1D meshpoints (including both enpoints, e.g., one cell has two points 
     ! two cells have three and so on. recall that the 3D grid is a cartesian product of the 1d meshes.      
     ! grids_u, grids_v, and grids_w contain 1D meshes for each of the grids. Meshes are following each other, 
     ! to find the appropriate records use the grid capacity array "grids_cap_u/_v/_w
real (DP), dimension (:), allocatable :: cells_luII, cells_lvII, cells_lwII, cells_ruII, cells_rvII, cells_rwII
integer (I4B), dimension (:), allocatable :: cells_pgridII, cells_cgridII, cells_refuII, cells_refvII, &
                                      cells_refwII, cells_gowII, cells_gouII, cells_govII
integer (I4B), dimension (:), allocatable :: cells_ugiII, cells_vgiII, cells_wgiII ! relative addresses on the parent grid, currently only works in the case of a single uniform grid .. 
     ! cells are the main object of the discretization, they form an ierarchical grid. 
     ! cells_pgrid is the number of the parent grid. 
     ! cells_cgrid is -1 if the cell is not refined or is equal to the number of the child grid. 
     ! cells_lu/_lv/_lw and cells_ru/_rv/_rw are the corners of the cells
     ! cells_refu/_refv/_refw are the refinement coefficients 
     ! cells_gou/_gow /_gov Cell's Gauss order. --- so far will be the same in both u, v, and w 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
real (DP), dimension (:), allocatable ::  nodes_uII, nodes_vII, nodes_wII, nodes_gwtsII
real (DP), dimension (:), allocatable ::  nodes_utII, nodes_vtII, nodes_wtII ! u - ubar, v - vbar, w - wbar. t in name comes from tilde. updated when the moments are updated
integer (I4B), dimension (:), allocatable :: nodes_pcellII, nodes_uiII, nodes_viII, nodes_wiII
     ! velocity nodes of the nodal-DG discretization. nodes are listed for each cell in the order cells are listed
     ! nodes_pcell is the number of the parent cell
     ! nodes_u/_v/_w are the coordinates of the node. 
     ! nodes_gwts are the values of products of weights of gaussian quadrature formulas to be used in 3D integration
     ! nodes_ui/_vi/_wi are the values of the local indices that determine addresses within the cell of the 1D nodal 
     ! points in each variable associated with the node
integer (I4B), dimension (:), allocatable :: nodes_AshiftII,nodes_duiII,nodes_dviII,nodes_dwiII,nodes_phicanII    
     ! this array contains the shift in the A index to be used for periodic grids/selkf similar basis functions
     ! nodes_dui/_dvi/_dwi contains the shift in cell idex for each node. = Shift in the index to go from the cell that has the node to the canonical cell. 
     ! nodes_phican contains the number of the canonical node on the canonical cell that corresponds to the given node.
integer (I4B), dimension (:), allocatable :: nodes_primcellII ! this arrays is to help with mapping the solution on the primary cell to the solution on the secondary cell
     ! for each node on the secondary mesh it contains the number of the cell on the primary mesh where the node of the secondary mesh falls. 
!!!!!!!!!!!!!!!!!!!!!!!!
! END of Variables of Secondary Velocity Grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The next group of variables is reserved for passing around values of the moments of kinetic solutios. 
! In many cases, the use of these can be avoided, but they are used by some legacy subroutines.
!
!!!!!!!!!!!!! moments of the distribution function f!!!!!!!!!!!!!!!
real (DP) :: nden,ubar,vbar,wbar,temp,temp_u,temp_v,temp_w
real (DP) :: mom3_u,mom3_v,mom3_w,mom4_u,mom4_v,mom4_w,mom5_u,mom5_v,mom5_w,mom6_u,mom6_v,mom6_w
real (DP) :: DfTxMom, Dfu3Mom, Dfu4Mom, Dfu5Mom, Dfu6Mom ! directional temperature moment in the x-direction about (fM - f)

!!!!!!!!!!!!!!!!!!!
real (DP) :: LocDens, LocUbar, LocVbar, LocWbar, LocTempr ! these are the variables to hold the local values of the macroparameters
!!!!!!!!!!!!!!!!!!!

! These are variables related to the linearization of the collision operator
real (DP), dimension (:), allocatable :: fm ! this is a storage to keep the local maxwellian
real (DP), dimension (:,:), allocatable :: fmA ! this is to keep the linearized operator
!!!!!!!!!!!!!!!!!!!                          


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters of the BGK model with velocity-dependent collision frequency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Attention: the number of cases in functions:psi_basis_funcs and kernls_enfrsd_moms must be at least as large as the parameters next: 
integer, parameter :: MaxNumEnforcedMoments=15 ! this is the maximum allowed number of enforced moments.
						! note that MaxNumEnforcedMoments is the total number of enforced moments. Depending on the 
						! choice of user, this may include the conserved moments in this number. However, it is not required by the method 
						! that the concerved quantities are enforced. 
integer, parameter :: MaxNumBFunctionsCollFreq=15 ! this is the maximum allowed number of basis fucntions 
              !in the representation of the collision frequency. 

integer :: Order ! this is the number of the enforced moments. To determine what moments are enforced, 
				 ! the user has to encode/provide the functions located in the module DGV_moment_defs
				 ! The functions are numbered concesutevily. Order determines how many of the tese functions are 
				 ! enforced. 
				 !  
integer :: Order_nu	 ! This is the number of the basis functions in the represntation of the collision frequency.
				 ! The user have to provide/encode the functions. Depending on the implementation of the method, 
				 ! number of the functions in the representation of the collision frequency may be equal or less than the
				 ! the number of enforced moments. If Order_nu<Order, then the coefficients of the collision frequency are found by
				 ! solving the linear least squares problem 
real (DP) :: mft_coeff ! This is the proportionality coefficient of local mean free time to be used to 
                 ! in dermination of the update time period for the relaxation rates of the enforced moments. 
                 ! The relaxation rates for moments are updated in time intervals. The time intervals are computed proportionately 
                 ! to the local mean free time. mft_coeff is the coefficient of the poportionality. Thus mft_coeff=1.0 means that 
                 ! reates will be updated every local mean free time (local mean free tiem is estiamted based on macroparameters and
                 ! gas constants inside the update rates subroutine. Practice shows that mft_coeff=0.5 or less should be used for time-dependent simulations.
                 ! The steady state probelms do not require frequenct updates and mft_coeff can be set at 10s and 100s of the mean free times.   
                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 0D variables that support velocity dependent collision frequency for spatially homogeneous problems(a single spatial cell)
!!!!!!!!			 
real (DP), dimension(:), allocatable :: MoRatesArry0D    ! this is the array to keep the coefficients of the 
				! of the velocity dependent collision frequency in the case of spatially homogeneous problem
logical :: MoRatesArry0DAllocated  = .false.			! determine if the array MoRatesArry0D has been allocated
logical :: MoRatesReliable0D = .false. ! this variable tells if the last computed relaxation rates were obtained from the Bolzmann collision integral
                ! or if the default rates were used for all entries. If at least one entry was computed from the Boltzmann collision integral, 
                ! the value of the fal is true. The values are computed from the Boltzmann integral is the value of signal are higher than then the error indicators. 
                ! otherwise default value of the relaxation rate is used.  


real (DP) :: NuNextUpdateTime0D ! this variable stores the next update time
real (DP) :: NuLastDens0D,NuLastTemp0D  ! these variables store the values of the density and temperature that were used to evaluate the velocity depemdent collision frequency
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mom matrix is used in the new BGK model with velocity dependent collision frequency model \nu(v).
! the coefficients of the velocity dependent collision frequency (VDCF) are determined by solving Mom*Cco = Dphi
! that is the coefficients Cco = Mom^{-1} Dphi
real (DP), dimension (:), allocatable :: Cco		! coefficients for the velocity dependent collision frequency
!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D variables that support velocity dependent collision frequency for multiple spatial cells :	
!!!!!!!!			 
integer (I4B) :: TotNumSpatCells ! This number will store the total number of spatial cells. used for checks
! 
real (DP), dimension(:,:), allocatable :: MoRatesArry1D  ! th`is is the array to keep the coefficients of the 
				! of the velocity dependent collision frequency in the case of 1D, 2D, and 3D problems. 
				! First index is the moments, second index is cells
logical :: MoRatesArry1DAllocated  = .false.
logical, dimension (:), allocatable :: MoRatesReliable1D ! this is the array to keep the flag if the last computed relaxation rates 
                ! for this cell were obtained from the Bolzmann collision integral
                ! or if the default rates were used for all entries. If at least one entry was computed from the Boltzmann collision integral, 
                ! the value stored in the array for this cell is true. The index goes over all cells in space. The values are computed from the Boltzmann integral is the value of signal are higher than then the error indicators. 
                ! otherwise default value of the relaxation rate is used.  

real (DP), dimension (:,:), allocatable :: Cco1D		! coefficients for the velocity dependent collision frequency
                 ! First index is the coefficients, second index is cells    
real (DP), dimension(:), allocatable :: NuNextUpdateTime1D,NuLastDens1D,NuLastTemp1D   ! these are arrays to store the 
               ! values of the next time the coefficients must be updated, the density and the temprature at the last update. Each array
               ! cell corresponds to one spatial point in the solution. 
logical :: Nu1DUpdateArrysAllocated = .false.

               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The names of computed solution file and the directory to store it
character (len=50) :: current_solution_dir  ! name of the directory (relative to the local tree?) where the solution is stored
character (len=20) :: current_solution_base_name ! base name for the solution file -- more information about order, mesh, and time is added to it automatically
!!!!
! varaibles related to generation of non-unifom meshes: 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! These parameters are describing the non-unofirm meshes in the velcity variable. 
! notice that they are currently not implemented -- reserved for future implementation
!  
integer (I4B) :: u_nonuniform_mesh_type, v_nonuniform_mesh_type, w_nonuniform_mesh_type ! parameter to describe non-uniform mesh in x and u 


! supported types of nonuminform mesh : 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! Variables that are used to take galerkin projections and evaluation of moments.
!! Not all are used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: moments_u_gauss_order ! order of gausss quadrateures in x and u to calculated (total) moments
integer (I4B) :: moments_refine_u ! coefficients of mesh refinement in x and u for the evaluation of moments
real (DP), dimension (:), allocatable :: moments_u_gauss_nodes     ! arrays to keep gauss nodes for evaluation of moments 
real (DP), dimension (:), allocatable :: moments_u_gauss_weights ! arrays to keep gauss weights for evaluation of moments  

real (DP), dimension (:,:), allocatable :: g_nds_all, g_wts_all ! arrays to store vaules of used gauss nodes and weights. 
                                                                 ! g_nds_all(i,j) i is the number of the node and j is the order of the gauss formula 
integer (I4B) :: g_max_order=10 ! maximum gauss order that will be used ...
                                           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Parameters specifying the gas
!
!
real (DP), parameter :: kBoltzmann = 1.3806503D-23
real (DP) :: mol_mass = 6.63d-26 ! variable to keep molecular mass 
!!!!!!!!!!!!!!! Parameters of Hard shepre model
real (DP) :: mol_diam !         molecular diameter                            

!!!!!!!!!!!!!!! Parameters of Dimensionless Reduction
! C_inf = termal velocity in m/s 
! L_inf = characteristic length in m
! N_inf = the total number of molecules in the volume
! T_inf =  the normalization for time is selected from the condition T_inf*C_inf = L_inf calculated automatically
! Temp_inf = the characteristic temperature
!!!!!!!!!!!!!!!
real (DP) :: T_inf, C_inf, L_inf, N_inf, Temp_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ES-BGK parameters
! alpha is the corrector term in the ES-BGK distribution for producing a Prandtl number closer to the true value (BGK produces Pr=1 only).
! This value has a range (need to check what this is) specific for keeping the Tensor positive definite (so we dont have exp(+))
real (DP) :: alpha ! factor to adjust the Prandtl number
real (DP) :: gas_viscosity, gas_T_reference, gas_alpha ! from constants from the power viscocity law. gas-viscotiy is the reference visc. 
               ! gas_T_reference is the reference temperature for the power viscocity law and gas_alpha is the exponent of the power viscosity law.
real (DP) :: gasR ! is the normal gas constant
real (DP) :: nu ! collision frequency
logical :: flag_use_VSS_diam ! if this flag is set true, the relaxation rates for moments are computed using the reference diameter of the VSS moment (is used to mimic real gas behavior in HS model)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Variables involved in the Operator A -- the kernel of the collision operator
!
!
character (len=20) :: current_Aoperator_base_name !base name for the files containing the A-operator.
integer (I4B) :: numchnks ! oprator A may be stored in several files -- this variable gived the total number of such files. 
!!!!!!!!!!!!! operator A and the related staff !!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: A_capphi ! A_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p})
                                       ! for each basis function with number i
integer (I4B), dimension (:), allocatable :: A_xi,A_xi1 !for each non-zero entry of A this one keeps the index of velocities that produced that non-zero entries. 
                                       ! We assume that all velocities are indexed with one index
integer (I4B), dimension (:), allocatable :: A_phi ! this is the index of the basis function fro which the non-zero entry was computed.                                        
real (DP), dimension (:), allocatable :: A ! this is the collision information operator A(\xi,\xi_{1},\phi^{j}_{p})
                                !  A(i) i -- is the index of nonzero entires, we need to use other A-arrrays to restore 
                                ! what velocities and wnat basis function this index correcpods  
                                ! for example, A_sind_xi(i) gives the first velocity, A_sind_xi1 gives the second velcity
                                ! and A_phi gives the index of the used basis function. 
real (DP) :: Trad ! cutoff radius for A-array (it is actually diameter of the collisio sphere -- distance between xi and xi1)
real (DP) :: ErrEps, ErrChi, min_sens ! errors in integrals in chi and epsilon and cutoff level for A
! Notice that min_sense will have to uses and they normally require values different in order of magnitude. 
! The historically first use of min_sense is in the subroutine for the evaluation of $A$. These values should be decresing 
! with resolution and essentially be equal to ErrChi.
! The second use is for after the truncation of pre-cmoputed A. 
complex (DP), dimension (:,:,:,:,:), allocatable :: FA_ce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Begin  Variables for Operator A supported on the secondary velocity mesh.
!
character (len=20) :: current_Aoperator_base_nameII !base name for the files containing the A-operator.
integer (I4B) :: numchnksII ! oprator A may be stored in several files -- this variable gived the total number of such files. 
!!!!!!!!!!!!! operator A and the related staff secondary grid!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: A_capphiII ! A_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p})
                                       ! for each basis function with number i
integer (I4B), dimension (:), allocatable :: A_xiII,A_xi1II !for each non-zero entry of A this one keeps the index of velocities that produced that non-zero entries. 
                                       ! We assume that all velocities are indexed with one index
integer (I4B), dimension (:), allocatable :: A_phiII ! this is the index of the basis function fro which the non-zero entry was computed.                                        
real (DP), dimension (:), allocatable :: AII! this is the collision information operator A(\xi,\xi_{1},\phi^{j}_{p})
                                !  A(i) i -- is the index of nonzero entires, we need to use other A-arrrays to restore 
                                ! what velocities and wnat basis function this index correcpods  
                                ! for example, A_sind_xi(i) gives the first velocity, A_sind_xi1 gives the second velcity
                                ! and A_phi gives the index of the used basis function. 
complex(DP), dimension(:,:), allocatable :: FA,FAII !This stores the fourier transform of A. This might be very large depending on the size of
                                ! of A. The last dimension gives the transform for the basis.
complex (DP), dimension(:,:), allocatable :: FAIIRV ! Fourier of A on the secondary grid. Two dimensional because the matrix will be full.
                                            !FAIIRV is the right eigenvectors of FAII. 
complex (DP), dimension(:), allocatable :: FAIIE !A single dimensional array that stores the eigenvalues of FAII

complex (DP), dimension(:,:,:,:,:), allocatable :: FA_poly, FA_polyII ! Stores the fast fourier transform of A but also supports different indices of A.
                                ! It has 5 indices, j',j'',i',i'',i. 
                                ! j',j'' correspond to the cell number
                                ! i',i'' correspond to the basis number within the cell
                                ! i corresponds to the basis that value of A was computed against, the canonical basis.

! End of Variables for Operator A on the secondary grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters and variables used in the Korobov of integration of the the collision bilinear form 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (1:7) :: korob_net_param ! this array will contain 7 numbers: p, a_1, ..., a_6 
                   ! that determine the Korbov net. This array is used for creating the A-array and also for the 
                   ! integration of the bilinear form. 
                   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Variables involved in the Operator Akor -- the kernel of the collision operator prepard for integration
!   using Korobov nets 
! These arrays will only be used to create the Akor arrays. Once they are created they will 
! be saved on the hard drive. To work with the created arrays arrays of pointers will be used. 
!
character (len=20) :: current_AKoroperator_base_name !base name for the files containing the A-operator.
integer (I4B) :: numKorchnks ! oprator Akor may be stored in several files -- this variable gives the total number of files that store a single net. 
integer (I4B) :: numKornets ! oprator Akor may be stored in several files -- this variable gives the total number of korobov nets for all Akor arrays that will be used.
integer (I4B), dimension (:), allocatable :: KorNetsChunks ! this arrays stores the information about the pre-computed values on Akor stored on the hard drive
                  ! and chunnks in each net
                  ! specifically, the length of the array is the number of nets and each records tell how many chunks is in the net
                  ! example:  KorNetsChunks = [3,4,2] will expect to pre-computed values for three nets in total. The Akor for the first net will come in three chunks, 
                  ! for the second net in 4 chunks and for the thirs in 2 chunks. When reading Akor from file, Akor01, Akor02 and Akor03 and so on will be used. 
                  ! when creating values for the operator A, array Akor will be used  
!!!!!!!!!!!!! operator A for korobov integration and the related staff !!!!!!!!!!!!!!!!!!!!!
integer (I4B), dimension (:), allocatable :: Akor_capphi,Akor_capphiII ! Akor_capphi(i) gives the number of non-zero entries in A(\xi,\xi_{1},\varphi^{j}_{p})
                                       ! for each basis function with number i for a given Korobov net
integer (I4B), dimension (:), allocatable :: Akor_k,Akor_kII !for each non-zero entry of Akor this one keeps the index of the Korbov 
                                       ! velocity node. Korobov velocities can be comouted from 
                                       ! p, a_1,\ldots, a_6 by a know formula: on [0,1] \xi_{i}(k)=\{ (a_{i}k)/p \}  that produced that non-zero entries. 
integer (I4B), dimension (:), allocatable :: Akor_phi,Akor_phiII ! this is the index of the basis function for which the non-zero entry was computed.                                        
real (DP), dimension (:), allocatable :: Akor,AkorII ! this is the collision information operator A(\xi(k),\xi_{1}(k),\phi^{j}_{p})
                                !  Akor(i) i -- is the index of nonzero entires, we need to use other A-arrrays to restore 
                                ! what velocities and wnat basis function this index correcpods  
                                ! for example, Akor_k(i) gives the index of the KOrobov node for which Akor(i) is 
                                ! computed. 
                                ! and Akor_phi(i) gives the index of the used basis function. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To evaluate the bilinear collision operator using the Korobov nodes, we will use 
! the floowing arrays of pointers. 
!
type DP1Darray_pointer
 real (DP), dimension (:), pointer :: p
end type DP1Darray_pointer

type I4B1Darray_pointer
 integer (I4B), dimension (:), pointer :: p
end type I4B1Darray_pointer

!
!  Now we define arrays of pointers. Each array will point to a collection of 
!  arrays that will play role of dublicates of Akor, Akor_k, Akor_phi, Akor_capphi 
!  that will correspond to different korobov nets. 
! 
type (DP1Darray_pointer), dimension(15) :: AkorAllNets, AkorAllNetsII
type (I4B1Darray_pointer), dimension(15) :: AkorAllNets_k, AkorAllNets_phi, AkorAllNets_capphi,&
                                        AkorAllNets_kII, AkorAllNets_phiII, AkorAllNets_capphiII,&
                                        korob_net_paramAllNets, korob_net_paramAllNetsII,&
                                        nodes_AkorshiftAllNets, nodes_AkorshiftAllNetsII 
! the pointers will be pointing to regions in the memory to that will store copies of 
!  Akor, Akor_k, Akor_phi, Akor_capphi  for each Korobov net used.
! for the description of the arrays, please read the comments above
!
! End variabels workign with Akor arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
integer (I4B), dimension (:), allocatable :: I1_list ! list of basis functions for which A-array should be computed
integer (I4B) :: I2_from, I2_to ! variables to designate the range of secodn inside loop in velocity node in the subroutine SetA_DGV                                     
                                ! these numbers are only invoked is the I1_list has only one record. If I1_list has mopre than one record the 
                                ! the loop in I2 will go over all nodes...                              
!!!!!!!!!!!!!!!!!!!
integer (I2B) :: Num_OMP_threads ! number of threads for OpenMP. 



!!!!!!!!!!!!!!!!!!!
integer (I4B) :: run_mode ! This is a variable that keeps track in which mode is the simulations: 
                          ! run_mode=0 is the simulation is strongly non-equilibrium. Set decomp_lev < |f_M-f|
                          ! run_mode=1 is when the simulation are close to equilibrium, that is that the local 
                          !            velocity distribution function is close to a maxwellian, however, it is 
                          !            not quite converged to the Maxwellian yet.  linear_lev < |f_M - f| < decomp_lev.
                          !            Max_decomp decides when the solution need to be decomposed into a Maxwell and $$ 
                          !
                          ! run_mode=2 is the linearized regime. In this regime the function is a perturbation of 
                          !            a Maxwellian and the norm of the pertrurbation is small so that quadratic term in the 
                          !            projection of the collision integral can be neglected. |f_{M}-f|<linear_lev, 
                          !            min_linear is the parameter that decides whether the addition is small enough. 
                          ! run_mode=3 is the regime to run BGK model with velocity-dependent collision frequency
                          !
                          !
                          ! run_mode=4 is the regime to run classical ES-BGK or Shakhov or BGK model 
                          !
                         
integer (I4B), dimension(:), allocatable :: run_mode_1D ! this is the array in wich run_mode is stored for each cell. The values are the same 
						  ! as for the variable run_mode

                          
real (DP) :: vel_lev, ES_lev, linear_lev, decomp_lev ! these are parameters detemining at which level linearization or perturbative decomposition starts. 
								            ! make sure that linear_lev is in the order of magnitude smaller than decomp_lev                         


real (DP), parameter :: pi25DT = 3.141592653589793238462643d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! MPI Variables for the evaluation of collision operator !!!!!!!!!
!!!

integer :: MPI_LINEAR_COLL_WORLD
integer :: MPI_LINEAR_COLL_WORLD_GROUP ! Communications and commpnicator group handles for creating a communicator
                                        ! That will contain processes participating in consolidated evaluation of linearied operator.  
integer, dimension(:), allocatable :: MPI_Acopy_UNIVS, MPI_Acopy_UNIVS_GROUP ! This array will hold pointers to universes that will receive a single copy of A_array.
integer, dimension(:), allocatable :: MPI_LinFmA_UNIVS, MPI_LinFmA_UNIVS_GROUP ! This array will hold pointers to universes that will receive a single copy of A_array.

                                          ! the first element will keep the number of receives and the other elements will keep the number of the universes
integer (I4B), dimension (:), allocatable :: procs_nodes_wld ! this array will store nodes that are assigned for a processor 
                                                               ! on which the program is running. the index runs trhought the nodes assigned for the 
                                                               ! evaluation. It is expected that the node knows the necessary peices of arrys A

integer (I4B), dimension (:), allocatable :: lin_proc_nodes_wld ! this array will collect the information about the local workload of the each processor in the 
                                                              ! group of processors responsible for the evaluation of the linearized collision operator
                                                              ! the index runs over assigned nodes. 
integer (I4B) :: My_Acopy_Univ_indx ! these keep the index of the Acopy universe that contain this proceesor   
integer  :: Lin_Self_Send_FmA_Flag ! this flag =0 if the master node of a Acopy Universe also appears as a processor recieving 
									! components of the linearized operator from this universe.   

integer (I4B), dimension (:,:), allocatable :: nodes_proc_groups ! this array will contain the beaking of velocity nodes among the Acopy universes
                     ! format nodes_proc_groups(a,b) a=1 or 2, b=the index of the Acopy universe, a=1,num_Acopies
                     ! nodes_proc_groups(1,b) is the first and nodes_proc_groups(2,b) is the last node assigned to the Acopy Univ #b
                     
integer (I4B), dimension(:), allocatable :: LinFmA_Univs_callorder ! entries of this array determine the order in wich LinFmA universes enter the 
                     ! collective communication. The entires are numbers from 1 to num_Acopies that are permuted in some way. 
                                                               
integer (I4B), dimension (:), allocatable ::  linprocndswld_Acopy_univs, linprocndswld_BfmA_addr, linprocndswld_Acopy_addr
        ! These arrays are located on the group of processors computing linearized solutions
        ! they store local information about what components from what universice the processor is 
        ! expecting as well as the information about the addresss of the components in the 
        ! broadcasted arrays 
        ! receive_proc_lin and send_process_lin
        !! in these arrays the information is in the following format:
        !! (# number of the LinFmA_recv_univ universe where the components of the linearized operator will come from, 
        !!  # total number of nodes for which the componets of the linearizeed operator need to be recieved from this universe
        !!  #node1, #node2, and so on, #nodeN, 
        !! then continue onto the next universe 
        !! the stucture of the array linprocndswld_Acopy_addr is the same as of linprocndswld_Acopy_univs
        !! exception that instead of the nodes, we have their indices in the corresponding Acopy universe's 
        !! Acopy_wrkld array. This array will be useful for the preparation of the linearized operator
        !! linprocndswld_BfmA_addr has the same structure except the records are indices of nodes from array LinFmA_recv_univ
        !! as they show up in the array BfmA
        !!
        !!
integer (I4B), dimension (:), allocatable :: Acopy_wrkld 
        !! This array is very similar to the array procs_nodes_wld, with the exception that is applies to the ]
        !! Acopy_universe rather than to individual processor. 
        !! Specifically, each Acopy universe will have some nodes assigned to it. All processors in this universe 
        !! will have a copy of this array containing the nodes that are assigned to their Acopy universe
        !! no nodes belong to two Acopy universes at the same time and no two processors belong to the same Acopy 
        !! universes, to there will be no ambiguity
integer (I4B), dimension (:), allocatable :: procs_nodes_wld_Acopy_addr
        !! each slave processor will contribute components of the linearized operator 
        !! in each Acopy universe an array will be created (FmA) containing the components of the linearized operator 
        !! the rows of this array will correspond to the nodes listed in the Acopy_wrkld array
        !! now Acopy_wrkld_addr will contain indices of nodes from the procs_nodes_wld on this processor in the 
        !! its A copy Universe's array Acopy_wrkld. The length of the  procs_nodes_wld_Acopy_addr is tha same as      
        !! array procs_nodes_wld. 
       
integer (I4B) :: num_lin_proc ! this variable will keep the numer of processors dedicated to the evaluation of the linearized operator. This number must be less than or equal to the 
                              ! total number of processors in the MPI universe. 
integer (I4B) :: num_Acopies ! this variable will store the desired number of times A is copied in memory between the processor. 
                             ! this is a good approach is A does not take a lot of storage, so say it fits in one node memory, say 16 procs.  
                             ! however there are many nodes and essentially we chunk th enodes between groups of processors. each gorup will have a copy of entire A                              
integer (I4B) :: num_linstep_wait,lin_no_update ! this is the numeber of steps to wait until start computing linear solutio using consolidated linearied collision operator                             
              ! lin_no_update is the counter to see how many times we computed linearied collision operator without updating it.
integer :: Acopy_irank ! MPI variable for universe Acopy. Gives the rank of this processor in the Acopy universe it belongs. Does not make sense on master node i.e. not valid with irank=0!!!                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: pad !Temporary. Need better solution than global variable?
                    
end module DGV_commvar