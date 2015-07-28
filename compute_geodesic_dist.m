function d = compute_geodesic_dist(shape,idx,th)
if nargin < 3
    th = 1e7;
end
n = size(shape.X,1);

f = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));
src = Inf(n,1);
src(idx) = 0;
d = fastmarchmex('march', f, double(src), double(th));
fastmarchmex('deinit', f);
