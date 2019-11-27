function ViewCells;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exact_sol_file_name = 'sol080909/exp0207_A2d8020_2_2kc7su7sv7sw3NXU5MuUU5MvVU5MwWU_aopr.dat';
fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
m1=fread(fid,1,'int32'); % reads the first dimension 
m2=fread(fid,1,'int32'); % reads the second dimension
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
A2d_patch = fread(fid,m1*m2,'double');
fclose(fid);
A2d_patch=reshape(A2d_patch,m1,m2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now we just need to plot it...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
surf(A2d_patch(1:2:m1,1:2:m2))