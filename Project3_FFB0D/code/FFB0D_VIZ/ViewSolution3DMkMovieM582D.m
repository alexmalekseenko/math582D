function ViewSolution3DMkMovieM582D

%%%%%%%%%%%%%
bsname = 'M300Trf4dt5D5_1su1sv1sw41MuUU41MvVU41MwWU'; % the common part of the name for grids, cells and nodes files

%%% Next we fill an array with all solutions files names in increasing time order
arrySFNames=[ 'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0002000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0019500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0037000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0054500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0072000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0089500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0107000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0124000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0141500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0159000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0176500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0194000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0211000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0228500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0246000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0263500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0281000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0298500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0315500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0333000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0350500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0368000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0385500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0403000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0420000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0437500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0455000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0472500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0490000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0507500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0594500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0612000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0629000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0646500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0664000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0681500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0699000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0716500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0733500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0751000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0768500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0786000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0803500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0821000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0838500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0855500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0873000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0890500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0908000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0925500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0943000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0960000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0977500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0995000000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.1012500000_sol.dat'; ...
              'M300Trf4dt5D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.1030000000_sol.dat'];

%%% Next we will call the subroutine that will make a graph and capture that graph into a movie
nframes = size(arrySFNames,1)
nframes = 56;

figure(4) %

aviobj = avifile('my3Dmovie.avi','fps',4,'compression','Cinepak'); 

clf;
for i=1:nframes
 exact_sol_file_name= arrySFNames(i,:);  
 h = MakeFrame3DcontourPlotDGV(bsname,exact_sol_file_name);
 view([35, 14]);
 axis([-1 1 -1 1 -.4 1.2]);
 M(i)=getframe;
 aviobj = addframe(aviobj,M(i));
 clf;
end 
aviobj = close(aviobj);
movie(M,2);


function h= MakeFrame3DcontourPlotDGV(bsname,exact_sol_file_name)
%%%%%
%%%%% This function will visualize a 2D section of the solution. 
%%%%% 
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE ONE: READ THE ARRAYS >>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads Grid Arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mesh_file_name = [bsname '_grids.dat'];
fid = fopen(mesh_file_name,'rb');
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
mesh_file_name = [bsname '_cells.dat'];
fid = fopen(mesh_file_name,'rb');
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
mesh_file_name = [bsname '_nodes.dat'];
fid = fopen(mesh_file_name,'rb');
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
% READS THE SOLUTION ARRAY
%%%%%%%%%%%%%

fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
sol_time = fread(fid,1,'double');
sol_dt = fread(fid,1,'double');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
rkmts=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
nn=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
f = fread(fid,nn,'double');
fclose(fid);
%%%%%%%%%%%%%%
% Finished reading DG solution. 
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
% END READING THE ARRAYS
%%%%%%%%%%%%%%

%%%%%%%%%%
%volume pictures? 
%%%%%%%%%%
u=reshape(nodes_u,41,41,41);
v=reshape(nodes_v,41,41,41);
w=reshape(nodes_w,41,41,41);
ffb = reshape(f,41,41,41);

%u=permute(u,[1,2,3]);
%v=permute(v,[3,2,1]);
%w=permute(w,[3,2,1]);
%ffb=permute(ffb,[3,2,1]);

maxf=max(f);
minf=min(f);


%%%%% Set Arrays of slices 
Sx=[-0.02:0.15:1.5];
Sy=[minf:(maxf-minf)/10:maxf]; 

contourslice(v,w,u,ffb,[0],[0],Sx,[minf:(maxf-minf)/10:maxf]);



