function [all_touched,coordinates]= geodesic_triangles(v_source,shoot_dir,shape,scale)

perp = [0;0;0];

[current_triangle,v_loc,o_loc,recon,d1,d2] = project_direction(v_source,shoot_dir,perp,shape);
coordinates = [];

iter =0;
shape.D = -shape.D;

doviz = 0;

all_touched = zeros(400,1);
all_touched_pairs = zeros(400,2);
while 1,
    v       = shape.TRIV(current_triangle,:);
    X(1,:)  = shape.X(v);
    X(2,:)  = shape.Y(v);
    X(3,:)  = shape.Z(v);
    
    if iter >0 && (shape.D(v(1)))<-scale, break; end
    
    iter =iter+1;
    
    if iter==1,
        u_prev =  zeros(1,3);
        u_prev(v_loc) = 0;

        u_prev(o_loc(1)) = d2./(d1+d2);
        u_prev(o_loc(2)) = d1./(d1+d2);
                
        direction_prev    = recon/100;
        if doviz,
            visualize_triangle(vizim,X,u_start,recon*30,'r');            
        end
    else
        triangle_prev = shape.TRIV(previous_triangle,:);
        points        = X_prev; 
        if ~u_prev(1),
            outside =1; 
            within = [2,3];
        elseif ~u_prev(2),
            outside=2; 
            within = [1,3];
        else outside=3; 
            within = [1,2];
        end        
        u_p    = u_prev([within,outside]);
                
        %% vertex only at new triangle
        curr = shape.TRIV(current_triangle,:);
        prev = ((curr==triangle_prev(1))|(curr==triangle_prev(2))|(curr==triangle_prev(3)));
        new_vert  = curr(~prev);
        new_Vert(1,1)  = shape.X(new_vert);
        new_Vert(2,1)  = shape.Y(new_vert);
        new_Vert(3,1)  = shape.Z(new_vert);                               
        
        %% unfolding
        points_us        = [points(:,[within,outside])];
        starting_point   = points_us*(u_p)';    
        
        points_us(:,end) = unfold(points_us,new_Vert);
        u_n              = get_new_barycentric(points_us,starting_point,direction_prev*10);
        
        points_us(:,end) = new_Vert;
        new_direction    = (points_us*u_n') - starting_point;

            idxs_new = [triangle_prev(within),new_vert];
            u_prev = 0*u_prev;
            for id = [1:3],
                %w = find(curr == idxs_new(id));
                %u_prev(w)   = sum(u_n(id);
                u_prev = u_prev + (curr == idxs_new(id)).*u_n(id);
                %u_prev(w)   = sum(u_n.*;
            end
        
        
        direction_prev = new_direction;
    end
    previous_triangle = current_triangle;
    X_prev = X;
    if iter==1,
        coordinates = X(:,v_loc);
    end
    coordinates(:,end+1)   = X_prev*u_prev(:);
    
    [mx,I] = max(u_prev);
    cur_touched = shape.TRIV(previous_triangle,I);
    all_touched(iter,1) = cur_touched;
    all_touched_pairs(iter,:) = [cur_touched previous_triangle];
    
    [current_triangle]=next_triangle_u_start(previous_triangle,u_prev,shape);
    if ...
            isempty(current_triangle) || ...
            ismember([ cur_touched current_triangle ],all_touched_pairs(1:iter,:),'rows') % to prevent infinete loops
        break
    end
 end
all_touched = all_touched(1:max(iter-1,1)); %((iter+1):end,1) = [];

function  [new_vertex_unfolded,vec,bs] = unfold(points_reference,new_Vert)

point_00 = points_reference(:,1);  % (A)
point_0l = points_reference(:,2);  % (B) 
point_xy = points_reference(:,3);  % (C) 

%% Plane Basis

%vec = zeros(3,2); 
%vec(:,1) = point_0l - point_00;
%vec(:,2) = point_xy - point_00;
vec = ([point_0l - point_00,point_xy - point_00]);
bs = vec;

%%           (0,l)
%%%   C''____B____C'
%        \   |\  /
%         \  | \/
%          \ | /\ C (x,y)
%    _______\|/________
%            A (0,0)
%%
%% common vertices: {A,B}
%% previous vertex :C
%% new vertex:      C'

%% d1: d(A,C')
%% d2: d(B,C')
%% l : d(A,B)

%% require:   d(B,C'') = d(B,C') ^ d(A,C'') = d(A,C')
%% (C'' == (x,y)) ->
%% (l - y)^2  + x^2 = d(B,C')
%% y^2 + x^2  = d(A,C')

df = [vec(:,1),new_Vert - point_00,new_Vert - point_0l];
nr = sum(pow_2(df),1);
l = sqrt(nr(1));
d1     = nr(2);
d2     = nr(3);
y      = ((d1 - d2) + l*l)/(2*l);
x      = -sqrt( d1 -  y*y);
%vec = single(vec); l = single(l);
%bs  = vec;
% Gram Schmidt
bs(:,1)  = vec(:,1)./l;
bs(:,2)  = vec(:,2) -  (vec(:,2)'*bs(:,1))*bs(:,1);
bs(:,2)  = bs(:,2)./sqrt(sum(pow_2(bs(:,2))));

new_vertex_unfolded =  point_00 + y*bs(:,1) + x*bs(:,2);

function baryc =  get_new_barycentric(points,start_point,direction)

end_point = start_point + 10*direction;
low_corner_left  = points(:,1);
low_corner_right = points(:,2);
peak             = points(:,3);


%
%            outside  *
%                    / \
%       end point o /   \
%                 |/     \
%       meeting   X       \
%                /|        \
%               / |         \
%   within(1)  *--x----------* within(2)
%                 start_point

%

[ ds_left, prop_diag_left] = ...
    lineline(start_point,end_point,low_corner_left, peak);

[ ds_right, prop_diag_right] = ...
    lineline(start_point,end_point,low_corner_right,peak);

if ~((prop_diag_left>=0)&(prop_diag_left<=1)),
    ds_left  = 10000;
end
if ~((prop_diag_right>=0)&(prop_diag_right<=1)),
    ds_right = 10000;
end
baryc  = zeros(1,3);
if ds_left<ds_right
    %% choose left coordinates
    baryc(2) = 0;
    baryc(1) = 1-prop_diag_left;
    baryc(3) = prop_diag_left;
else
    %% choose right coordinates
    baryc(1) = 0;
    baryc(2) = 1-prop_diag_right;
    baryc(3)   = prop_diag_right;
end


function [coef,recon,resid,vec,d1,d2] = reconstruct_on_triangle(points,shoot_dir)
vec(:,1)  = points(:,1) - points(:,3);
vec(:,2)  = points(:,2) - points(:,3);

design = vec'*vec;
inpr   = vec'*shoot_dir;

coef   = design\inpr;
%if nargout~=6,
%    coef = coef.*(all(coef>0));
%end
%coef
zs = 0;
if any(coef<=0),
    coef = 0*coef;
    zs = 1;
end
%coef = coef.*(all(coef>0));
recon =  vec*coef;
resid  = shoot_dir - recon;
if nargout==6,
    if ~zs
        prj = recon./sqrt(sum(pow_2(recon)));
    else
        prj = 0*recon;
    end
    d1 =  vec(:,1) - sum(vec(:,1).*prj)*prj;
    d2 =  vec(:,2) - sum(vec(:,2).*prj)*prj;
    d1 = sqrt(sum(pow_2(d1)));
    d2 = sqrt(sum(pow_2(d2)));
    if 1==2,
        figure,
        quiver3(0,0,0,shoot_dir(1),shoot_dir(2),shoot_dir(3),'color','r');
        hold on,
        quiver3(0,0,0,vec(1,1),vec(2,1),vec(3,1),'color','b');
        hold on,
        quiver3(0,0,0,vec(1,2),vec(2,2),vec(3,2),'color','g');
        hold on,
        quiver3(0,0,0,resid(1),resid(2),resid(3),'color','k')
        hold on,
        quiver3(0,0,0,recon(1),recon(2),recon(3),'color','y');
        hold on;
        legend({'to reconstruct','basis 1','basis 2','residual','reconstruction'});
        axis([-.01,.01,-.01,.01,-.01,.01]);
        %title(sprintf(' %.2f',coef));
    end
end


function [current_triangle,v_loc,o_loc,recon,d1,d2,dsts] = project_direction(v_source,shoot_dir,perp,shape);

triang_idxs = shape.idxs{v_source};
t_source{1}= triang_idxs(find(shape.TRIV(triang_idxs,1)==v_source));  %% all triangles where the starting node appears as node 1
t_source{2}= triang_idxs(find(shape.TRIV(triang_idxs,2)==v_source));  %% all triangles where the starting node appears as node 2
t_source{3}= triang_idxs(find(shape.TRIV(triang_idxs,3)==v_source));  %% all triangles where the starting node appears as node 3

cnt  = 0;
cnto = 0;
dsts = [];
offsets =[0,-1,1,-2,2,-3,3,-4,4,-5,5]/5;

%% loop over perturbations to direction vector
while 1,
    cnto = cnto  + 1;
    if cnto>11,
        break
    end
    %% loop over all triangles to which starting node is adjacent
    for cn = [1:3],
        switch cn,
            case 1,
                v_loc = [1];
                o_loc = [2,3];
            case 2,
                v_loc = [2];
                o_loc = [1,3];
            case 3,
                v_loc = [3];
                o_loc = [1,2];
        end
        if cnto==12,
            cnto;
        end
        
        for tr_cand = t_source{cn}
            cnt = cnt + 1;
            %% reconstruct vector within triangle and compute
            %% norm of residual from original
            
            triangle = shape.TRIV(tr_cand,:);
            points = [shape.X(triangle),shape.Y(triangle),shape.Z(triangle)]';
            points = points(:,[o_loc,v_loc]);
            shoot_dir_proj  = shoot_dir + offsets(cnto)*perp;
            
            [coef,recon,resid,~]  = reconstruct_on_triangle(points,shoot_dir_proj);
            %% triangle basis
            nrm(cnt)    = sqrt(sum(pow_2(resid)));
            vrt(cnt)    = v_loc;
            ort(cnt,:)  = o_loc;
            cnd(cnt)    = tr_cand;
            rcn(:,cnt)  = recon;
            offs(cnt)   = offsets(cnto);

            if cnto==-1,
                visualize_triangle(cnto,points,[0,0,1],10*shoot_dir);
                visualize_triangle(cnto,points,[0,0,1],shoot_dir_proj);
                title(sprintf(' %.3f %.3f %.3f',coef(1),coef(2),nrm(cnt)));
            end
        end
    end
    cnto;
    if any(nrm<.8),  break; end
end

[~,idx]  = min(nrm);
tr_cand  = cnd(idx);
recon    = rcn(:,idx);
v_loc    = vrt(idx);
o_loc    = ort(idx,:);
offs     = offs(idx);

current_triangle = tr_cand;
triangle = shape.TRIV(tr_cand,:);
points = [shape.X(triangle),shape.Y(triangle),shape.Z(triangle)]';

points = points(:,[o_loc,v_loc]);
[~,recon,~,~,d1,d2] = reconstruct_on_triangle(points,shoot_dir + offs*perp);

function [next_triangle]=next_triangle_u_start(current_triangle,u_on_edge,shape)

%%%% This function returns the number of the next neighbour triangle 
%%% Input: current_triangle-current triangle number
%%%        u_on_edge-barycentrical coordinates of the start point on the edge
%%%        shape- the mesh data structure
%%% Output: neighbour triangle number

zeroed_coordinate=find(abs(u_on_edge)<=eps);
if numel(zeroed_coordinate) ~= 1 % for the rare case that barycentrical coords are on a vertex
    [dummy,zeroed_coordinate] = min(u_on_edge); %#ok<ASGLU>
end

switch zeroed_coordinate
    case 1
        V = shape.TRIV(current_triangle,[2 3]);
    case 2
        V = shape.TRIV(current_triangle,[1 3]);
    case 3
        V = shape.TRIV(current_triangle,[1 2]);
end

triangles_with_vertex1 = shape.idxs{V(1)};
triangles_with_vertex2 = shape.idxs{V(2)};
next_triangle  = [];
for k =1:length(triangles_with_vertex1),
    trv = triangles_with_vertex1(k);
    if trv~=current_triangle,
        if any(triangles_with_vertex2 == trv),
            next_triangle = trv;
            break;
        end
    end
end







if 1==2,
    visualize_triangle(1,points,[0,0,1],recon);
end
if 1==2,
    
    if 1==2,
        visualize_triangle(36,points_reference,[1,0,0],vc(:,1),'r');
        visualize_triangle(36,points_reference,[1,0,0],vc(:,2),'r');
        visualize_triangle(36,trig,[1,0,0],0*vc(:,2),'k');
        visualize_triangle(37,points_reference,[1,0,0],bs(:,1)/10,'r');
        visualize_triangle(37,points_reference,[1,0,0],bs(:,2)/10,'b');
        visualize_triangle(37,trig,[1,0,0],0*vc(:,2),'k');
        
        visualize_triangle(34,[points_reference(:,[1,2]),new_Vert],u_p,direction_prev*10);
        visualize_triangle(33,points_unfolded, u_p,direction_prev*10);
    end
    if 1==2
        %visualize_triangle(vizim,triang,u_p,new_direction*10);
        visualize_triangle(vizim,triang,u_p,direction_prev*10);
        visualize_triangle(vizim,triang,u_p,new_direction*10,'r');
        %visualize_triangle(vizim,triang,u_n,new_direction*10,'r');
        hold on,scatter3(starting_point(1),starting_point(2),starting_point(3),'r','filled');
        %hold on,scatter3(new_ending(1),new_ending(2),new_ending(3),'b','filled');
        idx_within  = triangle_prev(within);idx_outside = triangle_prev(outside); idx_new = new_vert;
        visualize_triangle(21,points_reference,up_eff,recon/20,'b');
        visualize_triangle(21,points_unfolded,up_eff,recon/20,'b');
        visualize_triangle(21,points_unfolded,u_new,recon/20,'k');
        hold on,scatter3(x_start(1),x_start(2),x_start(3),'r','filled')
        hold on,scatter3(new_ending(1),new_ending(2),new_ending(3),'m','filled')
        visualize_triangle(22,points_unfolded,u_new,direction_new*10,'k');
    end
    
    if 1==2,
        visualize_triangle(8,unfolded_triangle,u_prev,shoot_dir);
        visualize_triangle(9,unfolded_triangle,u_unfolded,shoot_dir);
        hold on,scatter3(x_start(1),x_start(2),x_start(3),'r','filled')
        dp = recon/200;
        hold on,scatter3(x_start(1) + dp(1),x_start(2) + dp(2),x_start(3) + dp(3),'r','filled')
    end
    
end
if 1==2,
    v = shape.TRIV(current_triangle,:);
    x = shape.X(v')';
    y = shape.Y(v')';
    z = shape.Z(v')';
    X = [x; y; z];
    
    x_start = X*u_start(:);
    x_end = X*u_end(:);
    start_point = x_start;
    %figure(1+iter),clf;
    %figure(1+iter);
    
    figure(2)
    plot3([points(1,within(1)),points(1,within(2))],[points(2,within(1)),points(2,within(2))],[points(3,within(1)),points(3,within(2))],'r');
    hold on,
    plot3([points(1,within(1)),points(1,outside)],  [points(2,within(1)),points(2,outside)],  [points(3,within(1)),points(3,outside)],'b');
    hold on,
    plot3([points(1,within(2)),points(1,outside)],  [points(2,within(2)),points(2,outside)],  [points(3,within(2)),points(3,outside)],'g');
    hold on,
    plot3([x_start(1),x_end(1)],[x_start(2),x_end(2)],[x_start(3),x_end(3)],'c');
    hold on,
    reconsh = recon0/100;
    quiver3(start_point(1),start_point(2),start_point(3),reconsh(1),reconsh(2),reconsh(3),'k');
    hold on,
    reconsh = recon/1000;
    quiver3(start_point(1),start_point(2),start_point(3),reconsh(1),reconsh(2),reconsh(3),'m');
    hold on,
    scatter3(start_point(1), start_point(2), start_point(3), 20,'c','filled');
    if iter>1,
        hold on,
        scatter3(meeting_vert_left(1), meeting_vert_left(2), meeting_vert_left(3), 20,'y','filled');
        hold on,
        scatter3(meeting_vert_right(1),meeting_vert_right(2),meeting_vert_right(3),20,'m','filled');
        
        hold on,
        scatter3(low_corner_left(1), low_corner_left(2), low_corner_left(3), 10,'r','filled');
        hold on,
        scatter3(low_corner_right(1),low_corner_right(2),low_corner_right(3),10,'g','filled');
        hold on,
        scatter3(peak(1),peak(2),peak(3),10,'b','filled');
    end
    
    legend({'line with point','side1','side2','connection','shooting direction','shooting recon'});
    
    if 1==2,
        scatter3(points(1,within(1)),points(2,within(1)),points(3,within(1)),'r');
        hold on
        scatter3(points(1,within(2)),points(2,within(2)),points(3,within(2)),'b');
    end
end