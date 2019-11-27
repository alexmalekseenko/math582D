function Astruct
% This function helps to visualize the eigenvalues of operator A ... 

%%%%%%%%%%%% READ THE A-ARRAYS ... 
name='exp07122012_2kc1su1sv1sw3NXU7MuUU7MvVU7MwWU';
exact_sol_file_name = ['sol080909/' name '_grids.dat'];
fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
grids_cap_u = fread(fid,mm,'int32');
grids_cap_v = fread(fid,mm,'int32');
grids_cap_w = fread(fid,mm,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
grids_u = fread(fid,mm,'double'); 
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
grids_v = fread(fid,mm,'double'); 
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');
grids_w = fread(fid,mm,'double'); 
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads CELLS ARRAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exact_sol_file_name = ['sol080909/' name '_cells.dat'];
fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
cells_pgrid = fread(fid,mm,'int32');
cells_cgrid = fread(fid,mm,'int32');
cells_lu = fread(fid,mm,'double'); 
cells_lv = fread(fid,mm,'double');
cells_lw = fread(fid,mm,'double');
cells_ru = fread(fid,mm,'double');
cells_rv = fread(fid,mm,'double');
cells_rw = fread(fid,mm,'double');
cells_refu = fread(fid,mm,'int32');
cells_refv = fread(fid,mm,'int32');
cells_refw = fread(fid,mm,'int32');
cells_gow = fread(fid,mm,'int32');
cells_gou = fread(fid,mm,'int32');
cells_gov = fread(fid,mm,'int32');
fclose(fid);
%%%%%%%%%%%%
% READS NODES ARRAYS
%%%%%%%%%%%%
exact_sol_file_name = ['sol080909/' name '_nodes.dat']; 
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
%%%%%%%%%%%%%
% READS THE SOLUTION A-ARRAYS:
%%%%%%%%%%%%%
exact_sol_file_name = ['sol080909/' name '_Aarrs.dat']; 
fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
mm=fread(fid,1,'int32');
nn=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
A_capphi=fread(fid,mm,'int32');
A=fread(fid,nn,'double');
A_xi=fread(fid,nn,'int32');
A_xi1=fread(fid,nn,'int32');
A_phi=fread(fid,nn,'int32');
fclose(fid);
%%%%%%%%%%%%%%
% END READING THE ARRAYS
%%%%%%%%%%%%%%

%%% REFORMAT A as a matrix, 

Amat=sparse(mm,mm)
for i=1:nn
    Amat(A_xi(i),A_xi1(i))=A(i);
    Amat(A_xi1(i),A_xi(i))=A(i);
end 

%%% take a peak at Amat:

%spy(Amat(1:100,1:100));

%%% Find the eigenvaluse and the eigenvectors of the Amat

[D]=eigs(Amat) %,100,'LM')

[V,D,Flag]=eigs(Amat,10,'LM');

%[D]=eigs(Amat,7^3)

%D - eigenvalues
%V - eigenvectors

%%% Graph the eigenvecotor

%%%%%%% CREATE THE MESH FOR THE PLOT..
u1d=[-3:0.1:3]; % one dimensional arrays
v1d=[-3:0.1:3]; 
w1d=0;
[u2d,v2d]=meshgrid(u1d,v1d);
w2d= 0 + 0*u2d; 

%% Pick the component

for j=1:5
f = (V(:,j));

sol = AsseSolu2D(u2d,v2d,w2d,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);
figure(j)                         
surf(sol)
h=surf(u1d,v1d,sol,'LineStyle', '-', 'LineWidth', 0.5);
colormap(gray);

xlabel('u', 'FontSize', 18, 'FontWeight', 'Demi');
ylabel('v', 'FontSize', 18, 'FontWeight', 'Demi');
% 
%axis([0 1 0 1 0 1]);
% end drawing
set(gca, 'XTick', -3.0:1.5:3.0, 'YTick', -3.0:1.5:3.0);
set(gca, 'XLim', [-3 3], 'YLim', [-3 3]);  %'ZLim', [0 2]
end 