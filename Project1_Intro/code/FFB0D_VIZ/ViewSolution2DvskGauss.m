function ViewSolution2D
%%%%%
%%%%% This function will visualize a 2D section of the solution. 
%%%%% 
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAGE ONE: READ THE ARRAYS >>>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Reads Grid Arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dirname = ''
bsname = 'M300Trf4dt2D5_1su1sv1sw41MuUU41MvVU41MwWU' %_time0.0960000000_sol
exact_sol_file_name = [bsname '_grids.dat'];
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
exact_sol_file_name = [bsname '_cells.dat'];
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
exact_sol_file_name = [bsname '_nodes.dat'];
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
%exact_sol_file_name = [ 'M300Trf4dt2D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0035800000_sol.dat'];
%exact_sol_file_name = [ 'M300Trf4dt2D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0018400000_sol.dat'];
exact_sol_file_name = [ 'M300Trf4dt2D5_2kc1su1sv1sw3NXU41MuUU41MvVU41MwWU_time0.0070600000_sol.dat'];

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
%% Now will read kGauss
%%%%%%%%%%%%%%
exact_sol_file_name = ['kGssM3S2R258P10dt2D5_1su1sv1sw33MuUU33MvVU33MwWU_time0.0069800000_kGausol.dat'];
%exact_sol_file_name = ['kGssM3S2R258P10dt2D5_1su1sv1sw33MuUU33MvVU33MwWU_time0.0017600000_kGausol.dat'];
%exact_sol_file_name = ['kGssM3S2R258P10dt2D5_1su1sv1sw33MuUU33MvVU33MwWU_time0.0069800000_kGausol.dat'];
%exact_sol_file_name = ['kGssM3S3CR258P10dt1D5_1su1sv1sw21MuUU21MvVU21MwWU_time0.0070100000_kGausol.dat'];

fid = fopen(exact_sol_file_name,'rb');
dump = fread(fid,1,'int32');
sol_time = fread(fid,1,'double');
sol_dt = fread(fid,1,'double');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
sze_gausses=fread(fid,1,'int32');
dump = fread(fid,1,'int32');  % Each new write statement in Fortran adds 64 bits of crap!
dump = fread(fid,1,'int32');  %
gausses=fread(fid,sze_gausses,'double');
fclose(fid);

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

maxf=max(f)
minf=min(f)


%xnodes=[-3:6/41:3];
%x0=(xnodes(1:41)+xnodes(2:42))/2;
%[X,Y,Z] = MESHGRID(x0,x0,x0);

 %u matches to Z,v to X and w to Y
 %(X,Y,Z) matches to (v,w,u)
%break 

Sx=[-0.02:0.15:1.5]

figure(1)
contourslice(v,w,u,ffb,[0],[0],Sx,[minf:(maxf-minf)/10:maxf]);
%%%%%p=patch(isosurface(v,w,u,ffb,(minf+(maxf-minf)*0.7)));

%% assemble kGauss
sol_g = assemble2DKgauss(u,v,w,gausses);
figure(5)
contourslice(v,w,u,sol_g,[0],[0],Sx,[minf:(maxf-minf)/10:maxf]);

% now the difference: 

diff=abs(ffb-sol_g);
min_diff=0; %min(min(min(diff)));
max_diff=max(max(max(diff)))
levs_diff = [min_diff+0.01:(max_diff-min_diff)/12:max_diff];
figure(4)
contourslice(v,w,u,ffb-sol_g,[0],[0],Sx,levs_diff);


break

%%%%%%% CREATE THE MESH FOR THE PLOT..
u1d=[-2:0.08:3]; % one dimensional arrays
v1d=[-2:0.08:3]; 
w1d=0;

[u2d,v2d]=meshgrid(u1d,v1d);

w2d = 0 + 0*u2d; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assemble kGauss Solution %%%%%%%%%%%%%

sol_g = assemble2DKgauss(u2d,v2d,w2d,gausses);

figure(1)
h=surf(u1d,v1d,sol_g,'LineStyle', 'none', 'LineWidth', 0.5, 'EdgeAlpha', 0.2);
colormap(jet);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sol = AsseSolu2D(u2d,v2d,w2d,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);

figure(2)
h=surf(u1d,v1d,sol,'LineStyle', 'none', 'LineWidth', 0.5, 'EdgeAlpha', 0.2);
colormap(jet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
h=surf(u1d,v1d,sol-sol_g,'LineStyle', 'none', 'LineWidth', 0.5, 'EdgeAlpha', 0.2);
colormap(jet);


%set(h,'FaceLighting','flat','FaceColor','interp','AmbientStrength',0.4, ...
%        'EdgeLighting', 'gouraud');
%light1=light('Position',[0.99 0.13 0.02],'Style','infinite');
%light2=light('Position',[0.99 12.13 0.02],'Style','infinite');
%view(130.5,30);

%xlabel('u', 'FontName', 'Courier', 'FontSize', 26, 'FontWeight', 'Bold');
%ylabel('v', 'FontName', 'Courier', 'FontSize', 26, 'FontWeight', 'Bold');
 
%axis([0 1 0 1 0 3]);
%axis off
%set(gca, 'XTick', -3.0:3.0:3.0, 'YTick', -3.0:3.0:3.0, 'ZTick', 0.0:0.1:0.1, 'FontSize', 20, 'FontName', 'Courier');
%set(gca, 'XLim', [-3.0 3.0], 'YLim', [-3.0 3.0], 'ZLim', [0 3],'FontWeight', 'Demi');


% end drawing
%%% Optiopns for Double Constant relaxation 
%%% view(130.5,30);
%%% set(gca, 'XTick', -3.0:3.0:3.0, 'YTick', -3.0:3.0:3.0, 'ZTick', 0.0:0.1:0.1, 'FontSize', 20, 'FontName', 'Courier');
%%% set(gca, 'XLim', [-3.0 3.0], 'YLim', [-3.0 3.0], 'ZLim', [0 0.12],'FontWeight', 'Demi');
%%% light1=light('Position',[0.99 0.13 0.02],'Style','infinite');
%%% light2=light('Position',[0.99 12.13 0.02],'Style','infinite');
%%% Options for Mach 3 relaxation
%%%set(gca, 'XTick', -3.0:3.0:3.0, 'YTick', -3.0:3.0:3.0, 'ZTick', 0.0:1:2.0, 'FontSize', 20, 'FontName', 'Courier');
%%%set(gca, 'XLim', [-3.0 3.0], 'YLim', [-3.0 3.0], 'ZLim', [0 2],'FontWeight', 'Demi');
%%%light1=light('Position',[1.99 5.13 1.5],'Style','infinite');
%%%light2=light('Position',[-0.99 -1.13 1.0],'Style','infinite');
%%%view(-140.5, 30)


% Check mass: 
mass=sum(f.*nodes_gwts);
% Check momentum 
ubar=sum(f.*nodes_gwts.*nodes_u)/mass;
vbar=sum(f.*nodes_gwts.*nodes_v)/mass;
wbar=sum(f.*nodes_gwts.*nodes_w)/mass;
%s

% Check temperature
tempr = sum(f.*nodes_gwts.*((nodes_u-ubar).^2+(nodes_v-vbar).^2+(nodes_w-wbar).^2))/mass/3.0*2.0 % dimensionless temperature
% check mass Dimensionless Maxwell:
fMaxw=zeros(length(nodes_u),1);
for i=1:length(nodes_u) 
    fMaxw(i) = DimentionlessMaxwell(nodes_u(i),nodes_v(i),nodes_w(i),4.0,0.6123858,0.0,0.0,0.6833286);
end 
massMaxw=sum(fMaxw.*nodes_gwts);
error=sum(abs(f-fMaxw).*nodes_gwts);



%{
%% 2D PLOTS: 
%figure(2)
% section in v=0 through the tip of the maxwellian
%u1d=[-3:0.02:3]; % one dimensional arrays
%v1d=[0.0]; 
%w1d=[0.0];
%[u2d,v2d]=meshgrid(u1d,v1d);
%w2d = 0 + 0*u2d;
%SolSecv0w0 = AsseSolu2D(u2d,v2d,w2d,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
%                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
%                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);
                         
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
%                         
%SolSecv0w0=reshape(SolSecv0w0,1,length(v1d));
%
%MaxwellSectv0w0=DimentionlessMaxwell(u1d,v1d,w1d,4.0,0.6123858,0.0,0.0,0.6833286);
%MaxwellSectv0w0=reshape(MaxwellSectv0w0,1,length(v1d));%
%
%plot(v1d*706.7652406937231,SolSecv0w0,'r',v1d*706.7652406937231,MaxwellSectv0w0,'b');%

%}


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