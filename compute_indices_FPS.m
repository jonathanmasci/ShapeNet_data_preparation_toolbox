function indices = compute_indices_FPS(q,shape,flag_dist)

n = size(shape.X,1);

indices = zeros(q,1);
indices(1) = floor(rand*(n-1))+1;

if strcmp(flag_dist,'Eucl')
    Eucl_dists = pdist2([shape.X,shape.Y,shape.Z],[shape.X,shape.Y,shape.Z]);
end

for i = 2:q
    
    if strcmp(flag_dist,'geod')
        d = compute_geodesic_dist(shape,indices(find(indices~=0)));
    elseif strcmp(flag_dist,'Eucl')
        d = min(Eucl_dists(indices(find(indices~=0)),:),[],1);
    end
    
    [~,idx] = max(d);
    
    indices(i) = idx;
    
end

