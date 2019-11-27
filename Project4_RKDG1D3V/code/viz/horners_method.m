% Implements Horner's method.
% TAKES 
% a -- array of coefficients of a polynomial in format
%           in format [a_{n}, a_{n-2}, ... , a_{0}]
% x -- point or array of points where needs evaluation 
% RETURNS 
% result of division of polynomial P(x) by (x-x_0): P(x)=Q(x)(x-x_0)+P(x_0) 
% b_{0} --- values of a P(x_0)
% b -- coefficients of the quotient Q(x),  NOTICE that P'(x_0)=Q(x_0)
%
function [b0, b] = horners_method(a, x0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(a,2); % this is to count how many coefficients came in
m = size(x0,2); % this is to count how many values of x_0 came in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if n>0  % in case an non empty array came -- do Horner's method 
  b=[a(1)*ones(m,1)];                     % initialize b as array of just one element b_n=a_n  
  for k=2:1:n-1
    b=[b, a(k) + b(:,k-1).*(x0)'] ; % ***see the explanations below
  end     
  if n>1        % check for the trivial case when  the polynomial is just a constant (one coefficient) 
     b0 = a(n) + b(:,n-1).*(x0)';    % last coefficient
  else 
     b0 = a(n).*ones(m,1); b=[];          % there is only one last coefficient, polinomial Q(x) does not exist 
  end     
    %***
    %a(i)+b(n-1)*x0 calculates the new element according to the formula b_{n-1}=a_{n-1}+b_{n}x0
else   % in case array a is empty 
    b0=0;
    b=[];
end