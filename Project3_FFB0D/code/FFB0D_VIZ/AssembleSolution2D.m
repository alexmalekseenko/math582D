function sol = AssembleSolution2D(u,v,w,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
                             cells_refv,cells_refw,cells_cgrid,cells_gou,cells_gov,cells_gow,nodes_pcell,f, ...
                             nodes_ui,nodes_vi,nodes_wi);
sol=0;
%
%celli = locate_cell(u,v,w,grids_cap_u,grids_cap_v,grids_cap_w,grids_u,grids_v,grids_w,cells_refu, ... 
%                             cells_refv,cells_refw,cells_cgrid);
%                         
%f3d = sortf3d(f,nodes_pcell,nodes_ui,nodes_vi,nodes_wi,cells_gou(celli),cells_gov(celli),cells_gow(celli),celli);
%

