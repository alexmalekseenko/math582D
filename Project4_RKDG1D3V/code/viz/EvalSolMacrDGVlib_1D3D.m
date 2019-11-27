%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This function is used to evaluate dimensionless macroparamters of density, bulk velocity and
%% temperature as well as some higher order moments. 
%%
%% The subroutine takes a solution that is discretized in the velocity space using Nodal_DG (DGVlib)
%% and by DG in the physical space. The subroutine then performs the integration in the vleocity space 
%% using the meshes of the native velocity discretization. The result is a collection of Space-DG coefficients
%% that corresponds to different Kernels/different moments.
%%
%% The Kernels are coded in an external function. To change what Kernels are computed, change the external function.
%%
%% On the second step, the subroutine reconstructed the graphs of the moment from their DG discretizations. 
%% In 1D spatial variable a DG discretization with Legendre basis is used. -- amke sure this is what is used in the provided solution 
%% 
%% f is the discrete solution
%% xmesh is points where the moments solution need to be evaluated. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xn,xubar,xvbar,xwbar,xmoments] = ... 
     EvalSolMacrDGVlib_1D3D (xpoints,xmesh,nodes_u,nodes_v,nodes_w,nodes_gwts,n_mom,f) % 

% some description from fortran --- may will delete later on 
% real (DP), dimension (nxpoints,n_mom), intent (out) :: xmoments! arrays to store values of computed moments
%                                        ! nxpoints --- is the number of spatial points where the momenrts 
%                                       ! need to be evaluated.                                           
% real (DP), dimension (:), intent (in) :: xpoints
%                                        ! nodes in x where the macroparameters will be evaluated 
%(:) xmesh ! array of mesh points that make cells in variable x for DG discretization
%real (DP), dimension (:), intent (in) :: nodes_u,nodes_v,nodes_w ! u,v, and w components of the velocity nodes 
%real (DP), dimension (:), intent (in) :: nodes_gwts ! values of the gauss weights used in integreation in the velocity variable. 
%integer (I4B), intent (in) :: n_mom ! order in x,
%real (DP), dimension (0:,:,:), intent (in) :: f ! f12(p,m,j) -- coefficients of spectral decomposition (converted from charactristic variables) 
%                                   ! -- p is the index in the basis functions in x 
%                                   ! -- m is the index in the basis functions in u
%                                   ! -- j is the cell in x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First we need to set up some arrays for integration... 
epsilon=0.000000001;
%%% create space for macroparamters
xmoments = zeros(length(xpoints),n_mom);   
xn=zeros(length(xpoints));
xubar=zeros(length(xpoints));
xvbar=zeros(length(xpoints));
xwbar=zeros(length(xpoints));
%%%
k = size(f,1);
N = size(f,3);
M = length(nodes_u);
%%%%%
Nx = length(xmesh); 
%%%%% We will go point to point and we will reconstruct the DG solution in x for every nodal value. The 
%% result will be stored in f_loc. All moments will be calculated from that f_loc.
%%%%%%
for p = 1:length(xpoints) % loop in the points in x where moments need to be evaluated 
    % next we will find the interval on DG mesh for j where particular x point belongs:
    if (xpoints(p) <= xmesh(1) - epsilon) | (xpoints(p) >= xmesh(Nx) + epsilon)
        error('graphing point is outside the discrete mesh') ;
    else 
        % find the intervals on DG mesh for m and n:
        m_i=2;
        % find the interval on DG mesh for m:
        while (xpoints(p) > xmesh(m_i) + epsilon) & (m_i < Nx)
        m_i = m_i+1;
        end
        %% we found m_i such that xxmesh(m)\in [xmesh(m_i-1), xmesh(m_i)]
        % first we evaluate the basis functions in x: 
        xx=2*(xpoints(p)-(xmesh(m_i) + xmesh(m_i-1))/2)/(xmesh(m_i) - xmesh(m_i-1));
        phi=zeros(1,k);
        for pp=0:k-1
            phi(pp+1) = horners_method(legendre_p(pp),xx);
        end
        %%% Next we will reconstruct the value of the solution at every ordinate at point xpoints(p):
        f_loc = nodes_u*0; % set it equal to zero and the right size
        %%%
        for ii=1:M
         for pp=1:k   
          f_loc(ii) = f_loc(ii) + f(pp,ii,m_i-1)*phi(pp);   
         end 
        end     
        %% Next, we compute the value of density and bulk velocity of the moments at the point xpoint(p) in physical space
        n = sum(f_loc.*nodes_gwts);
        ubar=sum(nodes_u.*f_loc.*nodes_gwts);
        vbar=sum(nodes_v.*f_loc.*nodes_gwts);
        wbar=sum(nodes_w.*f_loc.*nodes_gwts);
        %
        ubar=ubar/n; vbar=vbar/n; wbar=wbar/n;
        %% Now we compute the rest of the moments
        xn(p)=n;
        xubar(p)=ubar;
        xvbar(p)=vbar;
        xwbar(p)=wbar;
        for i=1:n_mom 
         xmoments(p,i)=sum(momkernels(nodes_u,nodes_v,nodes_w,i,n,ubar,vbar,wbar).* ...
            f_loc.*nodes_gwts);
        end  
        %% reconstructed the moments at point xpoint(p)        
    end % end if 
end % -- end loop in p --- points x where moments need to be calculated. 

