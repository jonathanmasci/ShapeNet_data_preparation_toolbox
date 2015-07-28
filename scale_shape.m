function shape_scaled = scale_shape(shape, s)

shapeV = [shape.X shape.Y shape.Z];

meanV = mean([
    max(shape.X) , min(shape.X)
    max(shape.Y) , min(shape.Y)
    max(shape.Z) , min(shape.Z)],2)';

shapeV = bsxfun(@minus, shapeV, meanV) * s;
shapeV = bsxfun(@plus, shapeV, meanV);

shape_scaled = shape;
shape_scaled.TRIV = shape.TRIV;
shape_scaled.X = shapeV(:,1);
shape_scaled.Y = shapeV(:,2);
shape_scaled.Z = shapeV(:,3);