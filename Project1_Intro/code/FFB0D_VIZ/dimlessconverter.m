function dimlessconverter

kBlzm = 1.380648813e-23; % m^2 kg s^{-2} K^{-1} The Boltzmann constant
amu = 1.66053892173d-27; %kg  standard atomic mass unit = 1/12 of the mass of carbon 12. 
m_arg = 39.948*amu 
m_nit = 14.0067*amu; % per atom
R_arg = kBlzm/m_arg  % 208
R_nit = kBlzm/m_nit/2; % because there are two atoms togeter
d_arg = 188*(1.0d-12)*2 % Van der Waals diameter of argon 

%%%% CASE MACH 3 initial data. Estimate mean time between collisions: 
n_1=6.634184236615457e-6/m_arg;
u_1=967.78; v_1=0; w_1=0;
T_1=300;
%
n_2=1.989989930149256e-5/m_arg;
u_2=322.59; v_2=0; w_2=0;
T_2=1100;
%%% END MACH 3 
%%%% CASE MACH 1.55 initial data. Estimate mean time between collisions: 
%n_1=0.0001067613/m_arg 
%u_1=500.0180209; v_1=0; w_1=0;
%T_1=300;
%
%n_2 = 0.0001899076/m_arg
%u_2 = 281.0975399; v_2=0; w_2=0;
%T_2 = 464.3212637;
%%% END MACH 3 

%%%% CASE MACH 6 initial data. Estimate mean time between collisions: 
%n_1=0.0001067613/m_arg; 
%u_1=2096.849765; v_1=0; w_1=0;
%T_1=300;
%
%n_2 = 0.0003987328/m_arg
%u_2 = 561.4346267; v_2=0; w_2=0;
%T_2 = 4222.106139;
%%% END MACH 3 

% evaluate macroparameters of the mixture
[n,u,v,w,T]=SumTwoStreamsDimensional(n_1,u_1,v_1,w_1,T_1,n_2,u_2,v_2,w_2,T_2,R_arg)
% evaluate the mean collision time:
mci0 = MeanTimeBetweenCollisions(n,1,T,R_arg,d_arg)
mfp = MeanFreePath(n,1,d_arg)

%%%%% CASE 2 constant distributions: Estimate mean time between collisions: 
%n_1=5.0048993e+24; 
%u_1=106.0; v_1=0; w_1=0;
%T_1=300;
%%
%n_2=5.0048993e+24;
%u_2=-106.0; v_2=0; w_2=0;
%T_2=300;
%% evaluate macroparameters of the mixture
%[n,u,v,w,T]=SumTwoStreamsDimensional(n_1,u_1,v_1,w_1,T_1,n_2,u_2,v_2,w_2,T_2,R_arg);
% evaluate the mean collision time:
%mtbc = MeanTimeBetweenCollisions(n,1,T,R_arg,d_arg);

%Reference temperature and density : out of the BLUE pretty much,
T_0=4000 % Kelvins
n_0=1.0d+20 % number density 
L_0=1.0 % domain characteristic length
d_0=d_arg % molecular diameter (reference = argon)

C_inf = sqrt(2*R_arg*T_0)   % characteristic velocity
CapT = L_0/C_inf          % chjaracteristic time

final_time=10*mci0
final_time_dimls=final_time/CapT

[timeA_dls,xA_dls,n1A_dls,u1A_dls,v_dls,w_dls,Tempr1A_dls]=ConvertToDimensionlessAlekseenko(mci0,1.0,u,v,w,n,T,T_0,n_0,L_0,R_arg);
fprintf (1,'Parameters of the sum of streams in Alekseenkos code:\n');
fprintf (1,'time= %14.10f , x= %12.10E , n=%12.10E , u=%12.10f , v=%12.10f , w=%12.10f, Tempr=%12.10f\n', timeA_dls,xA_dls,n1A_dls,u1A_dls,v_dls,w_dls,Tempr1A_dls);
% first steream
[timeA_dls,xA_dls,n1A_dls,u1A_dls,v_dls,w_dls,Tempr1A_dls]=ConvertToDimensionlessAlekseenko(mci0,1.0,u_1,v_1,w_1,n_1,T_1,T_0,n_0,L_0,R_arg);
fprintf (1,'Parameters of the first stream to Put in Alekseenko code:\n');
fprintf (1,'time= %14.10f , x= %12.10E , n=%12.10E , u=%12.10f , v=%12.10f , w=%12.10f, Tempr=%12.10f\n', timeA_dls,xA_dls,n1A_dls,u1A_dls,v_dls,w_dls,Tempr1A_dls);
% second stream
[timeA_dls,xA_dls,n2A_dls,u2A_dls,v_dls,w_dls,Tempr2A_dls]=ConvertToDimensionlessAlekseenko(mci0,1.0,u_2,v_2,w_2,n_2,T_2,T_0,n_0,L_0,R_arg);
fprintf (1,'Parameters of the second stream to Put in Alekseenko code:\n');
fprintf (1,'time= %14.10f , x= %12.10E , n=%12.10E , u=%12.10f , v=%12.10f , w=%12.10f, Tempr=%12.10f\n', timeA_dls,xA_dls,n2A_dls,u2A_dls,v_dls,w_dls,Tempr2A_dls);



function [n,u,v,w,T]=SumTwoStreamsDimensional(n1,u1,v1,w1,T1,n2,u2,v2,w2,T2,R_gas);
n=n1+n2;
%%
u=(n1*u1+n2*u2)/(n1+n2);
v=(n1*v1+n2*v2)/(n1+n2);
w=(n1*w1+n2*w2)/(n1+n2);
%%
T=(T1*n1+T2*n2)/(n1+n2)+ ... 
       ((u1*u1+v1*v1+w1*w1)*n1+(u2*u2+v2*v2+w2*w2)*n2)/(n1+n2)/R_gas/3 - ...
       ((n1*u1+n2*u2)^2+(n1*v1+n2*v2)^2+(n1*w1+n2*w2)^2)/(n1+n2)/(n1+n2)/R_gas/3;

function [time,x,n,u,v,w,Tempr]=ConvertToDimensionalTcheremissine(time_dls,x_dls,u_dls,v_dls,w_dls,n_dls,Tempr_dls,T_0,n_0,R,d)
% Reference quantities  
v_0sqrd=R*T_0;
v_0=sqrt(R*T_0);   % characteristic velocity
lambdasqrd = sqrt(2)*pi*n_0*d^2; % characteristic length 
tau_0 = sqrt(lambdasqrd/v_0sqrd); % chjaracteristic time
% Converting to dimensional:
time=time_dls*tau_0;
x=x_dls*sqrt(lambdasqrd);
u=u_dls*v_0;v=v_dls*v_0;w=w_dls*v_0;
n=n_dls*n_0;
Tempr=Tempr_dls*T_0;

function [time,x,n,u,v,w,Tempr]=ConvertToDimensionalAlekseenko(time_dls,x_dls,u_dls,v_dls,w_dls,n_dls,Tempr_dls,T_0,n_0,L_0,R)
% Reference quantities  
C_inf = sqrt(2*R*T_0);   % characteristic velocity
CapT = L_0/C_inf;          % chjaracteristic time
% Converting to dimensional:
time=time_dls*CapT;
x=x_dls*L_0;
u=u_dls*C_inf; v=v_dls*C_inf; w=w_dls*C_inf;
n=n_dls*n_0/((L_0)^3);
Tempr=Tempr_dls*T_0;

function [time_dls,x_dls,n_dls,u_dls,v_dls,w_dls,Tempr_dls]=ConvertToDimensionlessAlekseenko(time,x,u,v,w,n,Tempr,T_0,n_0,L_0,R)
% Reference quantities  
C_inf = sqrt(2*R*T_0);   % characteristic velocity
CapT = L_0/C_inf;         % chjaracteristic time
% Converting to dimensional:
time_dls=time/CapT;
x_dls = x/L_0;
u_dls = u/C_inf; v_dls=v/C_inf; w_dls=w/C_inf;
n_dls = n/n_0*((L_0)^3);
Tempr_dls=Tempr/T_0;

function mci = MeanTimeBetweenCollisions(n_0,L_0,T_0,R,d) 
%n_0 == total number of molucles in the volume
%L_0 == characteristic size of the volume
%T_0 == reference temperature
%R - specific gas constant 
%d -- molecular diameter
mci=(L_0)^3/(4*n_0*d*d*sqrt(pi*R*T_0));

function mfp = MeanFreePath(n_0,L_0,d) 
%n_0 == total number of molucles in the volume
%L_0 == characteristic size of the volume
%T_0 == reference temperature
%R - specific gas constant 
%d -- molecular diameter
mfp=(L_0)^3/(n_0*d*d*sqrt(2)*pi);