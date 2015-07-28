function show_shape(shape,fld,idx);
if (nargin==1)|(isempty(fld)),
    fld = ones(size(shape.X));
    if nargin==3,
        fld(idx) = 0;
        fld = 1-fld;
    end
end


ha = axes;
hf = trisurf(shape.TRIV,shape.X,shape.Y,shape.Z,fld);
axis image;
axis off; shading interp; lighting phong; camlight head;