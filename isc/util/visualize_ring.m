 function visualize_ring(shape,start_vertex_loc,directions,fg)
if nargin==3,
    fg = 1;
end
u_start     = [1,0,0];
perp = [0;0;0];
wd = 5;
ha = axes;
%% draw net aroudn point
tt = shape.idxs{start_vertex_loc};
for k=1:length(tt),
    vrt_number     = setdiff(unique(shape.TRIV(tt(k),:)),start_vertex_loc);
    Xs=[shape.X(start_vertex_loc),shape.X(vrt_number)'];
    Ys=[shape.Y(start_vertex_loc),shape.Y(vrt_number)'];
    Zs=[shape.Z(start_vertex_loc),shape.Z(vrt_number)'];
    visualize_triangle(fg,[Xs;Ys;Zs],u_start,perp,'r',wd,'k',wd,[]);
end

%% draw angles around it
for k = 1:size(directions,2),
    dr = directions(:,k);
    dr = dr./sqrt(sum(dr.^2));
    visualize_triangle(fg,[Xs;Ys;Zs],u_start,directions(:,k)/2,'r',wd,'k',wd,k); %k); %[]);
end
axis off;
set(ha,'xtick',[],'ytick',[],'ztick',[],'xticklabel',[],'yticklabel',[],'zticklabel',[]);
set(ha,'zcolor',[1,1,1])
set(ha,'ycolor',[1,1,1])
set(ha,'xcolor',[1,1,1])
