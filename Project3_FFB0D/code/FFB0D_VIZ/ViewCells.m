function ViewCells;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads FILE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exact_sol_file_name = 'sol080909/test_2kc2su3sv4sw4NXU3MuUU3MvVU3MwWU_cells.dat';
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
% now we need to visualize this... 
%%%%%%%%%%%%%

[vert, fac]=make_rect(cells_lu(1),cells_ru(1),cells_lv(1),cells_rv(1),cells_lw(1),cells_rw(1));
mm
for i=2:1:mm
[v, f]=make_rect(cells_lu(i),cells_ru(i),cells_lv(i),cells_rv(i),cells_lw(i),cells_rw(i));
f = f + size(vert,1);
vert=cat(1,vert,v);
fac=cat(1,fac,f);
end
patch('Faces',fac,'Vertices',vert,'FaceAlpha',0.1,'FaceColor','m','EdgeColor','b')

function [v,f]=make_rect(xl,xr,yl,yr,zl,zr)
figure(1)
v=[xl,yl,zl;xl,yl,zr;xl,yr,zl;xl,yr,zr; xr,yl,zl;xr,yl,zr;xr,yr,zl;xr,yr,zr];
f=[1 2 6 5; 1 3 4 2;3 4 8 7; 6 8 7 5; 1 3 7 5 ];
%f=[1 2 6 5; 1 3 4 2; 1 3 7 5 ]
%f=[3 4 8 7; 1 3 7 5 ]

