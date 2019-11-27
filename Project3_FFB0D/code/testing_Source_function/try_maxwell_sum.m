function try_maxwell_sum

u=[-3:0.001:3];

figure(1)
f_1=fM(u,1,.4229,0);
plot(u,f_1);

figure(2)
f_2=fM(u,1,0.25,1.1);
plot(u,f_2)

figure(3)
plot(u,f_1,u,f_2,u,f_1+f_2);

function y=fM(u,n,T,ubar)

y = n*exp(-(u-ubar).^2/T);