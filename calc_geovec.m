function [geoVec, basis] = calc_geovec(evecs, evals, geovec_params)

sampleLocs = bsxfun(@minus, evals(:), geovec_params.evalSamples) / geovec_params.dEval;

basis = cubicB_spline(sampleLocs);

geoVec = evecs.^2 * basis;

% L2 normalization
if geovec_params.doNormalize    
    geoVec = bsxfun(@rdivide, geoVec, sqrt(sum(geoVec.^2, 2)));
end


function b = cubicB_spline(t)

b = zeros(size(t));

idx4 = find(t >= 0 & t < 1);
idx3 = find(t >= 1 & t < 2);
idx2 = find(t >= 2 & t < 3);
idx1 = find(t >= 3 & t <= 4);

b(idx4) = t(idx4).^3 / 6;
b(idx3) = ( -3*(t(idx3)-1).^3  +3*(t(idx3)-1).^2  +3*(t(idx3)-1)  + 1 ) / 6;
b(idx2) = ( +3*(t(idx2)-2).^3  -6*(t(idx2)-2).^2                  + 4 ) / 6;
b(idx1) = (   -(t(idx1)-3).^3  +3*(t(idx1)-3).^2  -3*(t(idx1)-3)  + 1 ) / 6;
