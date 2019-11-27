function test_sf02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function can be used to make  aplot of initial data for two Maxwellians 
% to see if they fit the domain nicely.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Maxwellian 1:
T1=1.0428571429 % dimensionless temperature  -- ratio to T_{\infyt}, 
d1=1.6094213293 % dimensionless density == ration to initial number density
u1=1.2283828258 % dimensionless bulk velocity  -- ratio to C_{\infty}
v1=0.0 % v1=?
w1=0.0 % w1=?

% Maxwellian 2: 
T2=1.031580199   %  dimensionless temperature  -- ratio to T_{\infyt}, 
d2=6.0108772843   %  dimensionless density == ration to initial number density
u2=0.3289013189   % dimensionless bulk velocity  -- ratio to C_{\infty}
v2=0.0 % v2=?
w2=0.0 % w2=
% 
%


u1d=[-3:0.1:3]; % one dimensional arrays
v1d=[-3:0.1:3]; 
w1d=[-3:0.1:3];

[v3d,w3d,u3d]=meshgrid(v1d,w1d,u1d);

% evaluate the first maxwellian 
f1= maxwellian_dist_3D (d1,T1,u1,v1,w1,u3d,v3d,w3d);

maxf=max(max(max(f1)));
minf=min(min(min(f1)));

%%%%% Set Arrays of slices 
Sx=[-3:0.1:3];  % slices in x
Sf=[minf:(maxf-minf)/10:maxf]; % values to be used in countor plots of the function 
%%%%%%%%%%%%
figure(3)
contourslice(v3d,w3d,u3d,f1,[0],[0],Sx,Sf);

% evaluate the first maxwellian 
f2 = maxwellian_dist_3D (d2,T2,u2,v2,w2,u3d,v3d,w3d);


maxf=max(max(max(f2)));
minf=min(min(min(f2)));

%%%%% Set Arrays of slices 
Sx=[-3:0.2:3];  % slices in x
Sf=[minf:(maxf-minf)/10:maxf]; % values to be used in countor plots of the function 
%%%%%%%%%%%%
figure(4)
contourslice(v3d,w3d,u3d,f2,[0],[0],Sx,Sf);

maxf=max(max(max(f1+f2)));
minf=min(min(min(f1+f2)));

%%%%% Set Arrays of slices 
Sx=[-3:0.2:3];  % slices in x
Sf=[minf:(maxf-minf)/20:maxf]; % values to be used in countor plots of the function 
%%%%%%%%%%%%
figure(5)
contourslice(v3d,w3d,u3d,f2+f1,[0],[0],Sx,Sf);



%%%%%%%%%%%%%%%%%
% maxwellian_dist_3D (n,T,u_0,v_0,w_0,u,v,w)
%
% evaluates dimensionless Maxwellian
%
% n - density
% T - temperature
% u_0,v_0,w_0 -- components of the bulk velocity. 
% 
% u,v,w -- 3D arrays ready to be used to plot 3D functions
%
%%%%%%%%%%%%

function y = maxwellian_dist_3D (n,T,u_0,v_0,w_0,u,v,w)

beta=(n*pi*T)^(-3/2);
y = exp(-((u-u_0).*(u-u_0)+(v-v_0).*(v-v_0)+(w-w_0).*(w-w_0))./max(T,0.00001))./beta;