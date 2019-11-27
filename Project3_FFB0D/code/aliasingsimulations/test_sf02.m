function test_sf01

% stationary maxwellian Shock waves initial data: 

% Before the shock wave,
T1=0.0428571429; % Temperature of first stream, 
rho1=1.6094213293*2.0; % density of first stream.
u1=1.2283828258; % bulk velocity of the first stream 

% After the shock wave, 
T2=0.6031580199; %temperature fo the second stream, 
rho2=6.0108772843*2.0; %density of the second stream 
u2=0.3289013189; % Mach2*sqrt(gamma*R*T2)=m/s.

%% MAKE SURE x starts at zero and u1>u2>0 

% exact function: 
uumesh = [-2000:50:3000];
xxmesh = [0:0.1:10];
[X,U]=meshgrid(xxmesh,uumesh);
%
exact_plot = bimodal_maxwellian(R,T1,u1,T2,u2,X,U);
%
figure(5);
%mesh(xxmesh,uumesh,exact_plot, 'EdgeColor','black');
h=surf(xxmesh,uumesh,exact_plot,'LineStyle', 'none', 'LineWidth', 0.5);
%set(h, 'LineStyle', 'none')
colormap(jet)

%mesh(xxmesh,uumesh,exact_plot, 'EdgeColor','black');
%h = get(3,'CurrentAxes')
%set(h, 'View', [0,90])

uumesh=[-2000:50:3000];
xxmesh=[0:0.1:10];
[X,U]=meshgrid(xxmesh,uumesh);

exact_plot= bimodal_maxwellian(R,T1,u1,T2,u2,X,U);
figure(4);
%surf(xxmesh,uumesh,exact_plot);
%h = get(4,'CurrentAxes')
%set(h, 'View', [0,90])
%colorbar

image([0,1500], [0,10], exact_plot, 'CDataMapping','scaled')
h = get(4,'CurrentAxes')
set(h, 'YDir', 'normal')
axis square
colormap(gray)
colorbar


function y = bimodal_maxwellian (R,T1,u1,T2,u2,X,U)
x0=0;
alpha=5*sqrt(2*R.*T1)/u2*(u1-u2)/(u1+u2)/3/sqrt(pi);
N1= 1./(1+exp((X-x0).*alpha/2)); %% assume lambda= !!!
y = maxwellian_dist (T1,R,u1,U).*N1 + maxwellian_dist (T2,R,u2,U).*(1-N1);


%% This function evaluates the 1D reduced maxweillian velocity distribution 
%% from the given average velocity, temperature and the gas constants
%% at velocity point u (possibly mesh matrix) at each point in (x) (possibly mesh matrix)   
%% u_0 -- scalar, 
%% R -scalar


function y = maxwellian_dist (T,R,u_0,u)

beta=2*R.*T;
y = exp(-(u-u_0).*(u-u_0)./max(beta,0.00001))./sqrt(pi*beta);