% This function implements the kernels of the moments that need to be computed from a kinetic solution. 
% it also includes division by the density and or coefficients, so all is left is just sum it with the kinetic solution and gauss weights.,
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here nodes_u,nodes_v,nodes_w are one or more u,v,w components of the points in velocity space where kernels need to be computed. 
% i is the number of the desired moment
% n,ubar,vbar,wbar -density, and components of the bulk velocity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = momkernels(nodes_u,nodes_v,nodes_w,i,n,ubar,vbar,wbar)

%% We will try to automatically adjust the dimensionality:
%%%%%%%%%%%%%%%%%%%%%%%
switch i
 case 1,
  y = 1 + nodes_u*0.0d0; % kernel density 
 case 2, 
  y = nodes_u/n;			%! u component of the kernel of momentum
 case 3,
  y = nodes_v/n;			%! v component of the kernel of momentum
 case 4,
  y = nodes_w/n;			%! w component of the kernel of momentum
 case 5,  
  y = ((nodes_u-ubar).^2 +(nodes_v-vbar).^2 +(nodes_w-wbar).^2)/3.0*2.0/n;	%! kernel of tempearture, please notice the dimensionless factor of $2/3$
 case 6,
  y = ((nodes_u-ubar).^2)/3.0*2.0/n;                                   % ! directional temperature in direction u
 case 7,
  y = ((nodes_v-vbar).^2)/3.0*2.0/n;                                   % ! directional temperature in direction u
 case 8, 
  y = ((nodes_u-ubar).^3)/n;
 case 9,
  y = ((nodes_v-vbar).^3)/n;
 case 10,
  y = ((nodes_u-ubar).^4)/n;
 case 11,
  y = ((nodes_v-vbar).^4)/n;
 case 12,
  y = ((nodes_u-ubar).^5)/n;
 case 13,
  y = ((nodes_v-vbar).^5)/n;
 case 14,
  y = ((nodes_u-ubar).^6)/n;
 case 15,
  y = ((nodes_v-vbar).^6)/n; 
 otherwise 
  fprintf(1,'can not process i:', i); 
  return
end   