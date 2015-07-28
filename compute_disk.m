function [in_ray,in_ring,shape,rays_from_point,directions,ds] ...
    = compute_disk(shape,start_vertex,dists,flag_dist,varargin)

% default parameters -----------------------------------------------------
% determine rings
scale_max = 80;
scale_min = 2;
nscales   = 25;
scales    = logspace(log(scale_min)/log(10),log(scale_max)/log(10),nscales+1);
% determine rays
N_rays = 8;
% ------------------------------------------------------------------------

%
expand_varargin;

%
tt         = shape.idxs{start_vertex};
directions = get1ring(shape,start_vertex,tt,N_rays);

% geodesic building phase
for ind_geod = [1:N_rays]
    direction = directions(:,ind_geod);
    [touched{ind_geod},rays_from_point{ind_geod}] = geodesic_triangles(start_vertex,direction,shape,max(scales));
end

% establish coordinate system
ds = zeros(N_rays,length(shape.D));
fm = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));
for ray = 1 : N_rays
    closest_to_ray = unique(touched{ray});
    if strcmp(flag_dist,'fmm')
        src = Inf(size(shape.X,1),1);
        src(closest_to_ray) = 0;
        tmp = fastmarchmex('march', fm, double(src), double(1e7));

%         tmp       = compute_geodesic_dist(shape,closest_to_ray);

        ds(ray,:) = tmp';
    elseif strcmp(flag_dist,'min')
        ds(ray,:) = min(dists(closest_to_ray,:),[],1);
    else
        error('[e] distance computation not supported');
    end
end
fastmarchmex('deinit', fm);
ds    = single(ds);
D     = single(max(shape.D',.00001));
scmax = scales(end) + 3;

% missing: 
% normalize in_{ray,scale}_soft to sum to 1 at each pixel

in_ray  = assign_to_rays(ds,D,single(fha),scmax);
in_ring = assign_to_scales(D,scales,fhs,scmax);

end

% ---------------------------------------------------- auxiliary functions
function in_ray = assign_to_rays(ds,D,soft_an,scmax)
Nrays = size(ds,1);
do_soft_or = soft_an>0;
wt  = find(D<scmax);

if do_soft_or,
    in_ray_pt = exp(-pow_2(ds(:,wt)./repmat(D(wt),[Nrays,1]))/soft_an);
    nrm = (max(sum(in_ray_pt,1),eps));
    in_ray_pt       = in_ray_pt./(repmat(nrm,[Nrays,1]));
    in_ray          = zeros(size(ds),'single');
    in_ray(:,wt)    = in_ray_pt;
else
    [~,idx] = min(ds,[],1);
    in_ray = ones(size(ds),'single');
    for ray = 1:Nrays
        in_ray(ray,(idx==ray)) = single(1);
    end
end
end

function in_ring = assign_to_scales(D,scales,fh,scmax)
nscales = length(scales)-1;
scs = [];
soft_sc = fh>0;
in_ring =  zeros(nscales,length(D),'single');
wt = find(D<scmax);
for sc = 1:nscales,
    if soft_sc
        scale_mid = (scales(sc) + scales(sc+1))/2;
        scs = [scs,scale_mid];
        if sc==1,
            in_ring(sc,wt) =  1./(1+exp((D(wt) - scale_mid)./(fh*scale_mid)));
        else
            dm = diff([scales(sc+1),scales(sc)]);
            in_ring(sc,wt) =  exp(-single(pow_2((D(wt) - scale_mid))./pow_2(fh*dm)));
        end
    else
        in_ring(sc,wt) =  double((D(wt)<=scales(sc+1))&(D(wt)>scales(sc)));
    end
end

if soft_sc,
    in_ring(:,wt) = in_ring(:,wt)./repmat(max(sum(in_ring(:,wt),1),eps),[nscales,1]);
end
end
