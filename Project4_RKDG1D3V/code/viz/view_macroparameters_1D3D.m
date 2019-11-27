%%%%%%%
%% This subroutine computes macroaparameters from a 1D3D kinetic solution (produced using DGVlib)
%%%%%%
function view_macroparameters_1D3D; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol_file_name = 'M380_G13_4k5s4mdgr5rk30NU10MU0.00000cfl0.18000000time_1D3Dsol.dat';
fid = fopen(sol_file_name,'rb');
dump = fread(fid,1,'int32');
k=fread(fid,1,'int32');   % degree of the dg polynomials
N=fread(fid,1,'int32');   % number of cells in x
xmesh = fread(fid,N+1,'double'); % mech in x
M=fread(fid,1,'int32');   % number of velocity nodes 
rk=fread(fid,1,'int32');  % order of RK methdo used in time integration
dt=fread(fid,1,'double'); % time step of the temporal discretization

dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
f = fread(fid,(k+1)*N*M,'double');
f = reshape(f,k+1,M,N);
fclose(fid);
xmesh(1)
xmesh(end)
%%%%%%%%%%%%
% Next we read the nodes of the nodal_DG velcity discetization
%%%%%%%%%%%%
exact_sol_file_name = ['M380_G13_5su5sv5sw10MuUU10MvVU10MwWU_nodes.dat'];
fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
nn=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
nodes_u = fread(fid,nn,'double');
nodes_v = fread(fid,nn,'double');
nodes_w = fread(fid,nn,'double');
nodes_gwts = fread(fid,nn,'double');
nodes_pcell = fread(fid,nn,'int32');
nodes_ui = fread(fid,nn,'int32');
nodes_vi = fread(fid,nn,'int32');
nodes_wi = fread(fid,nn,'int32');
fclose(fid);
%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALIZE
x_nodes=[0.0: 0.0005: 0.06];
%x_nodes=[0.0695442299: 0.000000000001:0.0695442300];
%% This one is to produce shifted files:
%x_left = 0;
%x_right = 0.2;
%dx=(x_right-x_left)/1000;
%x_center= 0.06954422998500;
%x_nodes=[x_center:dx:x_center+300*dx];
%x_nodes=[[x_center-300*dx:dx:x_center - dx] x_nodes];


%x_nodes=[ x_center ];
%x = x_center + dx;
%while x < x_right
%x_nodes = [x_nodes, x];
%x = x + dx;
%end 
%x = x_center - dx;
%while x > x_left
%x_nodes = [x x_nodes];
%x = x - dx;
%end 

n_mom = 13

[xn,xubar,xvbar,xwbar,xmoments] = EvalSolMacrDGVlib_1D3D (x_nodes,xmesh,nodes_u,nodes_v,nodes_w,nodes_gwts,n_mom,f);
%%%
%write_macropar_txt02(xmoments, 'M300HSES_4_4k5s4mdgr5rk20NU10MU0.00000cfl1.10030000time_');

figure(2);
hold on;
figure(2);
hold on;
plot(x_nodes, xmoments(:,1),'-r', 'LineWidth', 1, 'MarkerSize', 5);

figure(3);
hold on;
plot(x_nodes, xmoments(:,2),'-r', 'LineWidth', 1, 'MarkerSize', 5);

figure(4);
hold on;
plot(x_nodes, xmoments(:,5),'-r', 'LineWidth', 1, 'MarkerSize', 5);

figure(5);
hold on;
plot(x_nodes, xmoments(:,1).*xmoments(:,2),'-b', 'LineWidth', 1, 'MarkerSize', 5);

%plot(x_nodes, xmoments(:,1),'-b', 'LineWidth', 1, 'MarkerSize', 5);
%plot(x_nodes, (avel-u2)/(u1-u2),'-k', 'LineWidth', 1, 'MarkerSize', 5);
%plot(x_nodes, (tempr-T1)/(T2-T1),'-k', 'LineWidth', 1, 'MarkerSize', 5);

%% for the diffusion conditions ... 
%plot(x_nodes, (dens),'-k', 'LineWidth', 1, 'MarkerSize', 5);
%plot(x_nodes, (avel),'-k', 'LineWidth', 1, 'MarkerSize', 5);
%plot(x_nodes, (tempr),'-k', 'LineWidth', 1, 'MarkerSize', 5);


xlabel('x','FontSize',14);
%ylabel(' (n - n_1)  / (n_2-n_1) ','FontSize', 14, 'Position', [-0.035-.009,  .5]);
%ylabel(' (u - u_2)  / (u_1-u_2) ','FontSize', 14, 'Position', [-0.035-.009,  .5]);
%ylabel(' (T - T_1)  / (T_1-T_2) ','FontSize', 14, 'Position', [-0.035-.009,  .5]);
%ylabel(' n ','FontSize', 14) %, 'Position', [-0.035-.009,  .5]);
%ylabel(' u ','FontSize', 14)%, 'Position', [-0.035-.009,  .5]);
ylabel(' T ','FontSize', 14)%, 'Position', [-0.035-.009,  .5]);


%set(gca, 'XTick', -0.0:0.02:0.13, 'YTick', -0.2:.2:1.2);
%set(gca, 'XTick', -0.0:0.02:0.2);

%set(gca, 'XLim', [-0.03 0.13], 'YLim', [-0.2 1.2]) %, 'ZLim', [0 1]);
%set(gca, 'XLim', [0.0 0.2]) %, 'ZLim', [0 1]);


%write_macropar_txt01(x_nodes,dens,avel,tempr,eflux,mol_mass, k,s,max_deg,N,M,sol_file_name);
%print_macropar_txt01 (x_nodes,dens, avel, tempr,eflux, mol_mass, k,s,max_deg,N,M,sol_file_name);
%[x1,x2]= massa

function write_macropar_txt (x_nodes, dens, avel, tempr, k, s, max_deg, N, M, sol_file_name)

text_file_name = 'sol080909/macroparameters.txt' ;
%%% Prepare the data for writing ... .
macrp = zeros(4,length(x_nodes));
macrp(1,:)= x_nodes;
macrp(2,:)= dens';
macrp(3,:)= avel';
macrp(4,:)= tempr';
%%% end preapring data to be written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fid = fopen(text_file_name,'wt');
fprintf(fid,['Experiment of ' date '\n parameters: \n']);
fprintf(fid,'k= %2i, s=%2i, max_deg=%2i, N=%4i, M=%4i \n', k, s, max_deg, N, M);
fprintf(fid, ['solution file: ' sol_file_name '\n \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
fprintf(fid,'x_nodes \t \t Density \t \t  Average Velocity  \t \t Temperature \n' );
fprintf(fid,'%20.14f \t  %20.14f \t  %20.14f \t %20.14f \n', macrp );
fclose(fid);


function write_macropar_txt01 (x_nodes, dens, avel, tempr, eflux, mol_mass, k, s, max_deg, N, M, sol_file_name)

sol_file_name = sol_file_name(1:length(sol_file_name)-8);
text_file_name = [sol_file_name 'macro.txt'];

%%% Prepare the data for writing ... .
macrp = zeros(6,length(x_nodes));
macrp(1,:)= x_nodes;
macrp(2,:)= dens';
macrp(3,:)= avel';
macrp(4,:)= tempr';
macrp(5,:)= (avel').*(dens')*mol_mass;
macrp(6,:)= eflux';
%%% end preapring data to be written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fid = fopen(text_file_name,'wt');
fprintf(fid,['Experiment of ' date '\n parameters: \n']);
fprintf(fid,'k= %2i, s=%2i, max_deg=%2i, N=%4i, M=%4i \n', k, s, max_deg, N, M);
fprintf(fid, ['solution file: ' sol_file_name '\n \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
fprintf(fid,'x_nodes \t \t Density \t \t  Average Velocity  \t \t Temperature \t \t Mass Flux \t \t Kin.Energy Flux \n' );
fprintf(fid,'%20.14f \t  %20.14f \t  %20.14f \t %20.14f \t %20.14f \t %20.14f \n', macrp );
fclose(fid);

function print_macropar_txt (x_nodes, dens, avel, tempr, k, s, max_deg, N, M, sol_file_name)

%%% Prepare the data for writing ... .
macrp = zeros(4,length(x_nodes));
macrp(1,:)= x_nodes;
macrp(2,:)= dens';
macrp(3,:)= avel';
macrp(4,:)= tempr';
%%% end preapring data to be written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fprintf(1,['Experiment of ' date '\n parameters: \n']);
fprintf(1,'k= %2i, s=%2i, max_deg=%2i, N=%4i, M=%4i \n', k, s, max_deg, N, M);
fprintf(1, ['solution file: ' sol_file_name '\n \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
fprintf(1,'x_nodes \t \t Density \t \t  Average Velocity  \t \t Temperature \n' );
fprintf(1,'%20.14f \t  %20.14f \t  %20.14f \t %20.14f \n', macrp );

function print_macropar_txt01 (x_nodes, dens, avel, tempr, eflux,mol_mass, k, s, max_deg, N, M, sol_file_name)

%%% Prepare the data for writing ... .
macrp = zeros(6,length(x_nodes));
macrp(1,:)= x_nodes;
macrp(2,:)= dens';
macrp(3,:)= avel';
macrp(4,:)= tempr';
macrp(5,:)= (avel').*(dens')*mol_mass;
macrp(6,:)= eflux';
%%% end preapring data to be written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fprintf(1,['Experiment of ' date '\n parameters: \n']);
fprintf(1,'k= %2i, s=%2i, max_deg=%2i, N=%4i, M=%4i \n', k, s, max_deg, N, M);
fprintf(1, ['solution file: ' sol_file_name '\n \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
fprintf(1,'x_nodes \t \t Density \t \t  Ubar  \t \t Temperature \t \t Mass Flux \t \t Kin.Energy Flux \n' );
fprintf(1,'%20.14f \t  %20.14f \t  %20.14f \t %20.14f \t %20.14f \t %20.14f \n', macrp );

function write_macropar_txt02 (xmoments, macro_file_name_start)

text_file_name = [macro_file_name_start 'macro.txt'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fid = fopen(text_file_name,'wt');
fprintf(fid,['Experiment of ' date '\n parameters: \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
mytring = 'x_nodes, n, u_bar, v_bar, w_bar, T, Tu, Tv, (u-u_)^3, (v-v_)^3,'
mytring = [mytring  ' (u-u_)^4, (v-v_)^4, (u-u_)^5, (v-v_)^5, (u-u_)^6, (v-v_)^6 \n'];
fprintf(fid, mytring);

for i=1:size(xmoments,1) 
   mytring = num2str(xmoments(i,1),'%16.14f');
   for j=2:size(xmoments,2)-1 
    mytring = [mytring ', ' num2str(xmoments(i,j),'%16.14f') ];
   end  
   mytring = [mytring ', ' num2str(xmoments(i,end),'%16.14f') ];
   fprintf(fid, mytring );
end  
fclose(fid);



function [d1,d2]=massa
mass=6.63d-26; % for Argon gas

%%% This is for mach 10
rho1=1.0668d-4;   % density in kg/m^3.
d1=rho1/mass;     % is the density n(t,x)
rho2=4.1428d-4;   % density in kg/m^3.
d2=rho2/mass;    

%this is for mach 1.2
%rho1=1.0668d-4;
%d1=rho1/mass;     % is the density n(t,x)
%rho2=1.3837d-4;
%d2=rho2/mass;