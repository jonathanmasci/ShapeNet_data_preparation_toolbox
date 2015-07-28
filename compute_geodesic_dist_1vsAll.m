function D = compute_geodesic_dist_1vsAll(shape,idx,dist_th)

n = size(shape.X,1);

if nargin < 3
    dist_th = 1e7;
end

D = zeros(length(idx),n);

f = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));
for i = 1 : length(idx)
    src = Inf(n,1);
    src(idx(i)) = 0;
    d = fastmarchmex('march', f, double(src), double(dist_th));
    D(i,:) = d(:)';
end

fastmarchmex('deinit', f);

