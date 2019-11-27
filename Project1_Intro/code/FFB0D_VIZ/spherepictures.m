function spherepictures
%% CUBE %%
xL=0.3; xR=1.0
yL=0.3; yR=1.0;
zL=1.1; zR=1.7;
% sphere %%
xc=1.8,yc=1,zc=1 %% coordinates of the center
% direction of g
%gu=1;gv=-1;gw=0.5
gu=1.7;gv=-0.3;gw=0.5
gu=gu/sqrt(gu^2+gv^2+gw^2)/1.1;gv=gv/sqrt(gu^2+gv^2+gw^2)/1.1;gw=gw/sqrt(gu^2+gv^2+gw^2)/1.1;
xiu=xc+gu;xiv=yc+gv;xiw=zc+gw;
xi1u=xc-gu;xi1v=yc-gv;xi1w=zc-gw;
%
%%% DRAWING 
figure(1)
hold on
hidden on
%% The sphere
[x,y,z]=sphere(40);
c = (x-0.5).^2+(y+.5).^2;
%chi = recsupp(x+xc,y+yc,z+zc,xL,xR,yL,yR,zL,zR);
h=surf(x+xc,y+yc,z+zc,c);
axis equal
%colormap(hsv(128))
colormap(lines(128))
set(h, 'AlphaDataMapping','none','FaceAlpha',0.1)
%set(h, 'AlphaDataMapping','none','FaceAlpha','interp','AlphaData',chi)
set(h, 'EdgeAlpha', 0.5, 'Facecolor', 'interp')
set(gca, 'XLim', [-0.1 3.5], 'YLim', [-0.1 2], 'ZLim', [-0.1 2]);  
%%% The cube
%[vert, fac]=make_rect(xL,xR,yL,yR,zL,zR);
%patch('Faces',fac,'Vertices',vert,'FaceAlpha',0.1,'FaceColor','c','EdgeColor','w', 'EdgeAlpha', 1.0, 'LineWidth', 2 )
%%% The vectors
line([0.0 xiu], [0.0 xiv], [0.0 xiw],'LineWidth', 1.5)
line([0.0 xi1u], [0.0 xi1v], [0.0 xi1w],'LineWidth', 1.5)
line([xiu xi1u], [xiv xi1v], [xiw xi1w],'LineWidth', 1.5)
%%%%% ADD SOME MORE BALLSSSSSW
% direction of g
gu=3.0;gv=-.5;gw=-0.0
gu=gu/sqrt(gu^2+gv^2+gw^2)/1.1; gv=gv/sqrt(gu^2+gv^2+gw^2)/1.1; gw=gw/sqrt(gu^2+gv^2+gw^2)/1.1;
xi1u=xiu-2*gu;xi1v=xiv-2*gv;xi1w=xiw-2*gw;
xc=(xi1u+xiu)/2,yc=(xi1v+xiv)/2,zc=(xi1w+xiw)/2;
% the ball:
%% The sphere
[x,y,z]=sphere(40);
c = (x-0.7).^2+(y+.3).^2;
%chi = recsupp(x+xc,y+yc,z+zc,xL,xR,yL,yR,zL,zR);
h=surf(x+xc,y+yc,z+zc,c);
colormap(lines(128))
set(h, 'AlphaDataMapping','none','FaceAlpha',0.1)
set(h, 'EdgeAlpha', 0.5, 'Facecolor', 'interp')
%Vectors
line([0.0 xi1u], [0.0 xi1v], [0.0 xi1w],'LineWidth', 1.5)
line([xiu xi1u], [xiv xi1v], [xiw xi1w],'LineWidth', 1.5)
% direction of g
gu=10;gv=10;gw=0.2
gu=gu/sqrt(gu^2+gv^2+gw^2)/1.1;gv=gv/sqrt(gu^2+gv^2+gw^2)/1.1;gw=gw/sqrt(gu^2+gv^2+gw^2)/1.1;
xi1u=xiu-2*gu;xi1v=xiv-2*gv;xi1w=xiw-2*gw;
xc=(xi1u+xiu)/2,yc=(xi1v+xiv)/2,zc=(xi1w+xiw)/2;
% the ball:
%% The sphere
[x,y,z]=sphere(40);
c = (x-0.7).^2+(y+.3).^2;
%chi = recsupp(x+xc,y+yc,z+zc,xL,xR,yL,yR,zL,zR);
h=surf(x+xc,y+yc,z+zc,c);
colormap(lines(128))
set(h, 'AlphaDataMapping','none','FaceAlpha',0.1)
set(h, 'EdgeAlpha', 0.5, 'Facecolor', 'interp')
%Vectors
line([0.0 xi1u], [0.0 xi1v], [0.0 xi1w],'LineWidth', 1.5)
line([xiu xi1u], [xiv xi1v], [xiw xi1w],'LineWidth', 1.5)



hold off



function chi = recsupp(x,y,z,xL,xR,yL,yR,zL,zR)

chi=zeros(size(x));

[n,m]=size(x);
for j=1:m
    for i=1:n
        if ((x(i,j)>xL) & (x(i,j)<xR) & (y(i,j)>yL) & (y(i,j)<yR) & (z(i,j)>zL) & (z(i,j)<zR)) 
        chi(i,j)=0.8;
        else
        chi(i,j)=0.2;
        end 
    end 
end    

function [v,f]=make_rect(xl,xr,yl,yr,zl,zr)
v=[xl,yl,zl;xl,yl,zr;xl,yr,zl;xl,yr,zr; xr,yl,zl;xr,yl,zr;xr,yr,zl;xr,yr,zr];
%f=[1 2 6 5; 1 3 4 2; 3 4 8 7; 6 8 7 5; 1 3 7 5; 2 4 8 6];
f=[1 2 6 5; 1 3 4 2; 3 4 8 7; 6 8 7 5; 1 3 7 5; 2 4 8 6; 5 7 8 6;1 3 7 5];
%f=[1 2 6 5; 1 3 4 2; 1 3 7 5 ]
%f=[3 4 8 7; 1 3 7 5 ]
