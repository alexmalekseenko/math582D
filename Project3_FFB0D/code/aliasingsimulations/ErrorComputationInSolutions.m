function ViewSolution2D
%%%%%
%%%%% This function will compute errors in the numerical solutions 
%%%%%
%%%%% The function uses mesh and nodes that are saved on the hard drive. 
%%%%% The function computes the solution at those nodes and takes the integral
%%%%%
%%%%% 
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE ZERO: read the integration grids and integration nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integr_dirname = 'RedoM3\';
integr_bsname = 'M3_p36A35L36_2kc1su1sv1sw3NXU21MuUU21MvVU21MwWU';
[inodes_gwts, inodes_u, inodes_v, inodes_w] = read_integration_nodes_weights(integr_dirname,integr_bsname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE ONE: READ THE SOLUTION 1 and evaluate on the integration mesh>>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1dirname = 'RedoM3\';
f1bsname = 'M3_p36A35L36_2kc1su1sv1sw3NXU21MuUU21MvVU21MwWU'; %_time0.0000038500_sol
time1 = '0.0003200000'
f1= eval_solution_integr_nodes(inodes_u,inodes_v,inodes_w,f1dirname,f1bsname,time1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE TWO: READ THE SOLUTION 1 and evaluate on the integration mesh>>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2dirname = 'RedoM3\';
f2bsname = 'M3_p36A35L36_2kc1su1sv1sw3NXU21MuUU21MvVU21MwWU'; %_time0.0000038500_sol
time2 = '0.0003200000'
f2 = eval_solution_integr_nodes(inodes_u,inodes_v,inodes_w,f2dirname,f2bsname,time2);
%%% Use sum of two maxwellians for the second solution
% m1=1
% u1=0
% v1=0
% w1=0
% t1=1
% m2=1
% u2=0
% v2=0
% w2=0
% t2=1
% 
% 
%f2 = eval_double_maxwell_nodes(inodes_u,inodes_v,inodes_w,m1,u1,v1,w1,t1,m2,u2,v2,w2,t2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE THREE: Evaluate the error:
%%%%%%% 
% Check mass: 
mass=sum(f1.*inodes_gwts);
% Check momentum 
ubar=sum(f1.*inodes_gwts.*inodes_u)/mass;
vbar=sum(f1.*inodes_gwts.*inodes_v)/mass;
wbar=sum(f1.*inodes_gwts.*inodes_w)/mass;
% Check temperature
tempr = sum(f1.*inodes_gwts.*((inodes_u-ubar).^2+(inodes_v-vbar).^2+(inodes_w-wbar).^2))/mass/3.0*2.0 % dimensionless temperature
% L2_err
L2_abs_sqr = sum(((f1-f2).^2).*inodes_gwts);
L2_ref_sqr = sum((f2.^2).*inodes_gwts);
L2_rel = sqrt(L2_abs_sqr/L2_ref_sqr);
% 
fprintf(1,['solution 1 ' f1bsname '\n ']);
fprintf(1,['solution 1 mass: %20.16f \t b.vel. %20.16f, %20.16f, %20.16f \t temp. %20.16f \n '], mass, ubar, vbar, wbar, tempr);

fprintf(1,['solution 2 ' f2bsname '\n ']);
fprintf(1,['L2 error  %20.16f \n '], L2_rel);


% check mass Dimensionless Maxwell:


%fMaxw=zeros(length(nodes_u),1);
%for i=1:length(nodes_u) 
%    fMaxw(i) = DimentionlessMaxwell(nodes_u(i),nodes_v(i),nodes_w(i),4.0,0.6123858,0.0,0.0,0.6833286);
%end 
%massMaxw=sum(fMaxw.*nodes_gwts)
%error=sum(abs(f-fMaxw).*nodes_gwts);

%% 2D PLOTS: 
%figure(2)
% section in v=0 through the tip of the maxwellian
%u1d=[-3:0.01:3]; % one dimensional arrays
%v1d=[0.0]; 
%w1d=[0.0];
%[u2d,v2d]=meshgrid(u1d,v1d);
%w2d = 0 + 0*u2d;
%SolSecv0w0 = AsseSolu2D(u2d,v2d,w2d,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
%                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
%                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);
%                         
%SolSecv0w0=(reshape(SolSecv0w0,1,length(u1d)));
%
%Maxwell1Sectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,1.0,1.2247472,0.0,0.0,0.2);
%Maxwell1Sectv0w0=reshape(Maxwell1Sectv0w0,1,length(u1d));
%Maxwell2Sectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,3.0,0.4082448 ,0.0,0.0,0.7333333);
%Maxwell2Sectv0w0=reshape(Maxwell2Sectv0w0,1,length(u1d));

%MaxwellSectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,4.0,0.6123858,0.0,0.0,0.6833286);
%MaxwellSectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,1.0,0.0,0.0,0.0,1.0600000000);
%MaxwellSectv0w0=reshape(MaxwellSectv0w0,1,length(u1d));

%plot(u1d*706.7652406937231,SolSecv0w0,'r',u1d*706.7652406937231,MaxwellSectv0w0,'b');
%plot(u1d*706.7652406937231,SolSecv0w0,'r',u1d*706.7652406937231,Maxwell1Sectv0w0+Maxwell2Sectv0w0,'b');

%write_2arry1Ddata (SolSecv0w0, Maxwell1Sectv0w0+Maxwell2Sectv0w0, [bsname '_secv0w0init.txt'])

%figure(3)
% section in v=0 through the tip of the maxwellian
%u1d=[0.646510]; % one dimensional arrays
%v1d=[-3:0.01:3.0]; 
%w1d=[0.0];
%[u2d,v2d]=meshgrid(u1d,v1d);
%w2d = 0 + 0*u2d;
%SolSecv0w0 = AsseSolu2D(u2d,v2d,w2d,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
%                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
%                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);
                         
%SolSecv0w0=reshape(SolSecv0w0,1,length(v1d));
%
%MaxwellSectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,4.0,0.6123858,0.0,0.0,0.6833286);
%MaxwellSectv0w0=reshape(MaxwellSectv0w0,1,length(v1d));

%plot(v1d*706.7652406937231,SolSecv0w0,'r',v1d*706.7652406937231,MaxwellSectv0w0,'b');



% this function evaluates a 3D maxwellian, using one-dimensional arrays of velocity points u,v,w. It produces an arrays where u 
% correspond to the first index, v - to the second index and w - to the third
function y=DimentionlessMaxwell(u,v,w,n,u_bar,v_bar,w_bar,T)
pi25DT=3.141592653589793238462643d0;
y = zeros(length(u),length(v),length(w));
for iu=1:length(u)
    for iv=1:length(v)
        for iw=1:length(w)
        y(iu,iv,iw)=n*exp( -((u(iu)-u_bar)*(u(iu)-u_bar)+(v(iv)-v_bar)*(v(iv)-v_bar)+(w(iw)-w_bar)*(w(iw)-w_bar))/T)/(sqrt(pi25DT*T)*(pi25DT*T));    
        end 
    end
end

function write_2arry1Ddata (arry1, arry2, text_file_name )

%%% Prepare the data for writing ... .
macrp = zeros(2,length(arry1));
macrp(1,:)= arry1;
macrp(2,:)= arry2;
%%% end preapring data to be written

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW WE WILL WRITE THE TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Record the solution in the file
%% first we include a small signature. I mean state some of the parameters in the header of the file
fid = fopen(text_file_name,'wt');
fprintf(fid,['Experiment of ' date '\n parameters: \n']);
fprintf(fid,['solution file: ' text_file_name '\n \n']);
%% fprintf(fid,'init_time=%6.4f, final_time=%6.4f\n', initial_time, final_time);
%% Now go the moments in three columns
fprintf(fid,'array1 \t array2 \n' );
fprintf(fid,'%20.14f \t  %20.14f \n', macrp );
fclose(fid);

function [nodes_gwts, nodes_u, nodes_v, nodes_w]=read_integration_nodes_weights(dirname,bsname)
%exact_sol_file_name = [dirname bsname '_grids.dat'];
%fid = fopen(exact_sol_file_name,'rb');
%dump = fread(fid,1,'int32');
%mm=fread(fid,1,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');  %
%grids_cap_u = fread(fid,mm,'int32');
%grids_cap_v = fread(fid,mm,'int32');
%grids_cap_w = fread(fid,mm,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%mm=fread(fid,1,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%grids_u = fread(fid,mm,'double'); 
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%mm=fread(fid,1,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%grids_v = fread(fid,mm,'double'); 
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%mm=fread(fid,1,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');
%grids_w = fread(fid,mm,'double'); 
%fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads CELLS ARRAYS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exact_sol_file_name = [dirname bsname '_cells.dat'];
%fid = fopen(exact_sol_file_name,'rb');
%dump = fread(fid,1,'int32');
%mm=fread(fid,1,'int32');
%dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
%dump = fread(fid,1,'int32');  %
%cells_pgrid = fread(fid,mm,'int32');
%cells_cgrid = fread(fid,mm,'int32');
%cells_lu = fread(fid,mm,'double'); 
%cells_lv = fread(fid,mm,'double');
%cells_lw = fread(fid,mm,'double');
%cells_ru = fread(fid,mm,'double');
%cells_rv = fread(fid,mm,'double');
%cells_rw = fread(fid,mm,'double');
%cells_refu = fread(fid,mm,'int32');
%cells_refv = fread(fid,mm,'int32');
%cells_refw = fread(fid,mm,'int32');
%cells_gow = fread(fid,mm,'int32');
%cells_gou = fread(fid,mm,'int32');
%cells_gov = fread(fid,mm,'int32');
%fclose(fid);
%%%%%%%%%%%%
% READS NODES ARRAYS
%%%%%%%%%%%%

exact_sol_file_name = [dirname bsname '_nodes.dat'];
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

function inodes_sol = eval_solution_integr_nodes(in_u, in_v, in_w, dirname, bsname, time)
exact_sol_file_name = [dirname bsname '_grids.dat'];
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
exact_sol_file_name = [dirname bsname '_cells.dat'];
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
exact_sol_file_name = [dirname bsname '_nodes.dat'];
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
% READS THE SOLUTION ARRAY
%%%%%%%%%%%%%
exact_sol_file_name = [dirname bsname '_time' time '_sol.dat'];
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
% END READING THE ARRAYS
%%%%%%%%%%%%%%

%%%%%%% evaluate on the provided mesh 
inodes_sol = AsseSolu2D(in_u,in_v,in_w,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);
                    
function fMaxw = eval_double_maxwell_nodes(in_u, in_v, in_w,m1,u1,v1,w1,t1,m2,u2,v2,w2,t2)
                         
% evaluate Dimensionless Maxwell:


fMaxw=zeros(length(in_u),1);
for i=1:length(in_u) 
    fMaxw(i) = DimentionlessMaxwell(in_u(i),in_v(i),in_w(i),m1,u1,v1,w1,t1) ...
        + DimentionlessMaxwell(in_u(i),in_v(i),in_w(i),m2,u2,v2,w2,t2);
end 
