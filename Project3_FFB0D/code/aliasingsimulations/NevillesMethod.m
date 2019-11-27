function f = NevillesMethod(xnodes,ynodes,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This functio evaluates interpolation Lagrange's polynomial 
% using the method of Lagrange.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=size(xnodes,1);
mtr=zeros(n,n); % get a storage ready... 
%
mtr(:,1) = ynodes;
%
for i=2:n
    for j=1:n-i+1
mtr(j,i)=((xnodes(j+i-1)-x)*mtr(j,i-1)+(x-xnodes(j))*mtr(j+1,i-1))/(xnodes(j+i-1)-xnodes(j)) ;
    end 
end
f=mtr(1,n);

%function test
%
%xnodes=[1:2:10];
%ynodes=sin(xnodes);
%
%n=300;
%x=[1:9/n:10];
%y=0.*x;
%for i=1:n+1
%    y(i) = NevillesMethod(xnodes,ynodes,x(i));
%end
%plot(x,y)