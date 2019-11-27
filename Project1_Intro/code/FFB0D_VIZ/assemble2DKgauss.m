function sol = assemble2DKgauss(u,v,w,gausses);
% this subroutine assembles the k gauss solution from the arrray of Gaussians --- t
% the same format as in the kGauss code dens, bvel, 6 comp to compute inverse stress tensor.

% input --- 3 2d arrays of points (u,v,w) coordinates, 
% input gausses --- array of macroparameters.

nn=size(gausses,1)/10; % number of gaussians

sol = u*0; % create an empty array for the result

for j=1:nn
   %%
   yyy= evalGaussDens(u,v,w,gausses(((j-1)*10+1):j*10));
   sol=sol+yyy; 
end    
  

%%%%% function to evaluate Gaussian Density.
%%%%% L=(a, 0, 0 ; d  b, 0 ; e, m , g)
%%%%% Theta^{-1}=L*L^T
%%%% f(u)_{Gauss} = n*a*b*g*(\pi)^{-3/2}*exp((u-ubar)^T Theta^{-1} (u-ubar))


   function y=evalGaussDens(u,v,w,gauss);
   % \pi: 
   pi25DT = 3.141592653589793238462643;
   n = gauss(1);
   u_0 = gauss(2);  %zeta vector bulk velocity u term
   v_0 = gauss(3);  %zeta vector bulk velocity v term
   w_0 = gauss(4);  %zeta vector bulk velocity w term
   a = gauss(5);  % functions parameterizing L^-1 
   b = gauss(6) ;  
   g = gauss(7);
   d = gauss(8);
   e = gauss(9);
   m = gauss(10);
   
   % Compute the inside of the Gaussian 
   qc1 = a*a*(u-u_0) + a*d*(v-v_0) + a*e*(w-w_0) ;
   qc2 = a*d*(u-u_0) + (d*d + b*b)*(v-v_0) + (d*e+m*b)*(w-w_0) ;
   qc3 = a*e*(u-u_0) + (d*e+m*b)*(v-v_0) + (e*e+m*m+g*g)*(w-w_0) ;
   qq = -qc1.*(u-u_0) - qc2.*(v-v_0) - qc3.*(w-w_0);
   %
   z=sqrt(pi25DT)*(pi25DT);
   y = n*a*b*g*exp(qq)/z;
