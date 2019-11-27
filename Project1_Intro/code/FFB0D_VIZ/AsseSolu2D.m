function sol = AsseSolu2D(u,v,w,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
                             nodes_ui,nodes_vi,nodes_wi,cells_lu,cells_lv,cells_lw,cells_ru,cells_rv,cells_rw);

%%%%%%%%%%%%%%% Ppeporatory work %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before we do anything, we need to prepare some arrays: 
maxdeg=10;
%
gauss_nodes = zeros(10,10);
gauss_weights = zeros(10,10);
%
gauss_nodes(1,1) = 0;
gauss_weights(1,1) = 2;
%
for i=2:10 
   [xx,ww] = gaussquad(i,-1,1);
   gauss_nodes(i,1:i)=xx';
   gauss_weights(i,1:i)=ww';
end 
%
sol=zeros(size(u)); % this will create the storage for the function values...
% Now we are preparing to evaluate the solution:
n=size(u,1); % take the sizes of the arrays
m=size(u,2); % 
% Main loop: 
for i=1:n
    for j=1:m
% Find the cell where the point in velocity space falls:
celli = locate_cell(u(i,j),v(i,j),w(i,j),grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid);
% cut the portion of the solution that gives collocation points on this cells. Sort it in a 3D array                         
gou=cells_gou(celli);
gov=cells_gov(celli);
gow=cells_gow(celli);
f3d = sortf3d(f,nodes_pcell,nodes_ui,nodes_vi,nodes_wi,gou,gov,gow,celli);
% Now we evaluate the normalized coordinates: 
u_norm = (u(i,j)-(cells_lu(celli)+cells_ru(celli))/2)/(cells_ru(celli)-cells_lu(celli))*2;
v_norm = (v(i,j)-(cells_lv(celli)+cells_rv(celli))/2)/(cells_rv(celli)-cells_lv(celli))*2;
w_norm = (w(i,j)-(cells_lw(celli)+cells_rw(celli))/2)/(cells_rw(celli)-cells_lw(celli))*2;
% Now we will start evaluating the solution at this point. 
% we will assemble the third dimension, first and then second and then first: here we go 
f2d=zeros(gou,gov);
for ii=1:gou
    for jj=1:gov
        % assemble in the third dimension ... W
        ynodes=reshape(f3d(ii,jj,:),gow,1);
        xnodes=reshape(gauss_nodes(gow,1:gow),gow,1);
        f2d(ii,jj) = NevillesMethod(xnodes,ynodes,w_norm);
    end 
end     
f1d=zeros(gou,1);
for ii=1:gou
    % assemble in the second dimension ... V
    ynodes=reshape(f2d(ii,:),gov,1);
    xnodes=reshape(gauss_nodes(gov,1:gov),gov,1);
    f1d(ii)= NevillesMethod(xnodes,ynodes,v_norm);
end     
% and finally, assemble in the first dimension ... U 
xnodes = reshape(gauss_nodes(gou,1:gou),gou,1);
sol(i,j) = NevillesMethod((gauss_nodes(gou,1:gou))',f1d,u_norm);
% DONE WITH THE ASSEMPBILE point (u(i,j),v(i,j),w(i,j))
    end 
end 
% DONE ASSEMBLING...


function celli = locate_cell(u,v,w,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid);
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% first, given a velocity, we need to identify the cell where it came from and its 1D numbers 
% we start looking in course cells first. If a cell is found that has that velocity we check if the cell is refined. 
% if it is, then we look into the corresponding grid and seek the smaller cell that has it and so on. 
%!!!!!!!!

j=1;
%% set up the shift in grids_u/_v/_w corresponding to the first grid (zero shift):
gi_u=0;
gi_v=0;
gi_w=0;
while (j <= size(grids_cap_u,1))
  [ui,vi,wi] = finduiviwi(grids_u(gi_u+1:gi_u+grids_cap_u(j)),grids_v(gi_v+1:gi_v+grids_cap_v(j)), ...
                       grids_w(gi_w+1:gi_w+grids_cap_w(j)),u,v,w);
  if (ui*vi*wi ~= 0)  
     % now we need to find the cell that correspond to this grid and these ui,vi,wi
     celli=1;
     for jj=1:j-1
      celli = celli + (grids_cap_u(jj)-1)*(grids_cap_v(jj)-1)*(grids_cap_w(jj)-1);
     end 
     celli = celli + (ui-2)*(grids_cap_v(j)-1)*(grids_cap_w(j)-1)+(vi-2)*(grids_cap_w(j)-1)+(wi-2);
     %  check if this cell is refined
     if ((cells_refu(celli)>1) | (cells_refv(celli)>1) | (cells_refw(celli)>1)) then  
        j=cells_cgrid(celli)
        % set up the shift in grids_u/_v/_w corresponding to the j'th grid:
        gi_u=0;
        gi_v=0;
        gi_w=0;
        for jj=1:j-1
         gi_u = gi_u + grids_cap_u(jj);
         gi_v = gi_v + grids_cap_v(jj);
         gi_w = gi_w + grids_cap_w(jj);
        end
     else 
        break;
     end
  else 
     celli=0;
     break
  end                      
end
% now "celli" either is the number of the cell or 0 (0 means that the velocity in not on any cell)



function [ui,vi,wi]=finduiviwi(ugrid,vgrid,wgrid,u,v,w)

epsill=10d-11; % a small parameter to help curb the effects of the round off error..  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the grid sizes:
un=size(ugrid,1);
vn=size(vgrid,1);
wn=size(wgrid,1);
% now we start to look for the intervals where coordinates of $(u,v,w)$ belong.
if ((u < ugrid(1)-epsill) | (u > ugrid(un) + epsill)) 
  ui=0;
else
  ui=1;
  while ((u >=ugrid(ui) - epsill) & (ui < un))
  ui=ui+1;
  end 
end
%
if ((v < vgrid(1)-epsill) | (v > vgrid(vn) + epsill)) 
  vi=0;
else
  vi=1;
  while ((v >= vgrid(vi)-epsill) & (vi < vn))
  vi=vi+1;
  end
end
%
if ((w < wgrid(1)-epsill) | (w > wgrid(wn) + epsill))
  wi=0;
else
  wi=1;
  while ((w >= wgrid(wi)-epsill) & (wi<wn))
  wi=wi+1;
  end
end


%%%%%%%%%%%%%%%%%%%%%%
% THis function sets up a convenioent 3D array that stores 
% nodal values of the solution (f) on the cell with number (celli)
%%%%%%%%%%%%%%%%%%%%%% 

function f3d = sortf3d(f,nodes_pcell,nodes_ui,nodes_vi,nodes_wi,gou,gov,gow,celli)

f3d=zeros(gou,gov,gow);
for i=1:size(f,1)
  if (nodes_pcell(i) == celli)
      f3d(nodes_ui(i),nodes_vi(i),nodes_wi(i))=f(i);
  end 
end

