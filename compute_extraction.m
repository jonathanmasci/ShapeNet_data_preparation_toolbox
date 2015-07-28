function [M, dists] = compute_extraction(shape, params)

% parameters
rad       = params.rad; % radius used for descriptor construction
flag_dist = params.flag_dist;
nbinsr    = params.nbinsr; 
nbinsth   = params.nbinsth;
fhs       = params.fhs;
fha       = params.fha;

rr        = [1:nbinsr]/nbinsr*rad;
th        = [1:nbinsth]/nbinsth*2*pi;

% compute nearest-neighbor
fprintf('[i] computing nearest-neighbor... ');
start_time   = tic;
shape.idxs   = compute_vertex_face_ring(shape.TRIV');
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);
clear start_time elapsed_time;

% compute signature: just to double check with ISC and our M
fprintf('[i] computing signature... ');
start_time   = tic;
[desc,shape] = signature(shape,'wks');
elapsed_time = toc(start_time);
fprintf('%2.4fs\n',elapsed_time);
clear start_time elapsed_time;

% compute geodesic distances
if isfield(shape, 'dists')
    fprintf('[i] loading geodesic distances directly from shape... \n');
    dists = shape.dists;
else
    fprintf('[i] computing geodesic distances... ');
    if params.geod_th
        fprintf('\n    thresholded')
    end
    start_time = tic;
    idxs_  = compute_indices_FPS(50, shape, 'geod');
    dists_ = compute_geodesic_dist_1vsAll(shape,idxs_,1e+07);
    diam = max(dists_(:));
    fprintf('\n    diameter %f', diam)
    shape = scale_shape(shape,1/diam);
    
    dists = zeros(size(shape.X,1),size(shape.X,1));
    for i = 1:size(shape.X,1)
        if params.geod_th
            dists(i,:) = compute_geodesic_dist_1vsAll(shape,i,params.rad*1.5);
        else
            dists(i,:) = compute_geodesic_dist_1vsAll(shape,i);
        end
    end
    elapsed_time = toc(start_time);
    fprintf('%2.4fs\n',elapsed_time);
    clear start_time elapsed_time;
end

% compute disk for each vertex
M = cell(size(shape.X,1),1);
fprintf('[i] computing disk for each vertex...\n');
start_time = tic;
bad_disks = 0;
for i = 1:size(shape.X,1)
    
    if mod(i,100) == 0
       fprintf('    %d/%d %2.0fs\n',i,size(shape.X,1),toc(start_time));
    end
    
    % compute disks
    shape.D = dists(i,:)';
    
    % make an empty disk in case it fails
    try
        [in_ray,in_ring,shape,geod,directions,ds] = ...
            compute_disk(shape,i,dists,flag_dist,'scales',[0,rr],'N_rays',length(th),'fhs',fhs,'fha',fha);

        % compute descriptor on disk
        areascaling = full(diag(shape.A)); 
        [desc_net_,M_] = get_descriptor_from_net(in_ray,in_ring,desc,areascaling); 
    
        % sanity check with ISC
        % shuffle dimensions so that fastest index corresponds 
        % to bins and then to rings
        desc_net_shuffled = permute(desc_net_, [1,3,2]);
        desc_net = M_ * desc;
        if norm(reshape(desc_net_shuffled,100,80)' - desc_net) > 1e-04
            error('[e] descriptors not matching');
        end
    catch
        warning('something went wrong, making an empty disk')
        M_ = sparse(zeros(80,size(shape.X,1)));
        size(M_)
        bad_disks = bad_disks + 1;
    end
    M{i} = M_;
    clear M_;
    
end
fprintf('\n%i bad disks for this shape\n\n', bad_disks)
