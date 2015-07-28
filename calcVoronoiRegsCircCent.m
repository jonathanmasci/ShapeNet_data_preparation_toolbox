function [VorRegsVertices] = calcVoronoiRegsCircCent(Tri, Vertices)

%% Preps.:
A1 = Vertices(Tri(:,1), :);
A2 = Vertices(Tri(:,2), :);
A3 = Vertices(Tri(:,3), :);
a = A1 - A2;    % Nx3
b = A3 - A2;    % Nx3
c = A1 - A3;    % Nx3
M1 = 1/2*(A2 + A3);    % Nx3
M2 = 1/2*(A1 + A3);    % Nx3
M3 = 1/2*(A2 + A1);    % Nx3
N = size(A1, 1);

%% Circumcenters calculation
O = zeros(size(A1));

obtuseAngMat = [(dot(a, b, 2) < 0), (dot(-b, c, 2) < 0), (dot(-c, -a, 2) < 0)];
obtuseAngInds = any(obtuseAngMat, 2);
O(obtuseAngInds, :) = ...
    M1(obtuseAngInds, :).*(obtuseAngMat(obtuseAngInds, 1)*[1 1 1]) + ...
    M2(obtuseAngInds, :).*(obtuseAngMat(obtuseAngInds, 2)*[1 1 1]) + ...
    M3(obtuseAngInds, :).*(obtuseAngMat(obtuseAngInds, 3)*[1 1 1]);

OM3 = -repmat(dot(c, a, 2), 1, 3).*b + repmat(dot(b, a, 2), 1, 3).*c;
OM1 = -repmat(dot(c, b, 2), 1, 3).*a + repmat(dot(a, b, 2), 1, 3).*c;
M1M3 = M1 - M3;
tmp = M3 + OM3.*repmat(dot(cross(M1M3, OM1, 2), cross(OM3, OM1, 2), 2), 1, 3)./...
    repmat(dot(cross(OM3, OM1, 2), cross(OM3, OM1, 2), 2), 1, 3);
O(not(obtuseAngInds), :) = tmp(not(obtuseAngInds), :);

%% Voronoi Regions calculation (for each vertex in each triangle.
VorRegs = zeros(N, 3);
% For all the triangles do (though the calculation is correct for
% non-obtuse triangles only:
VorRegs(:,1) = calcArea(A1, M3, O) + calcArea(A1, M2, O);
VorRegs(:,2) = calcArea(A2, M1, O) + calcArea(A2, M3, O);
VorRegs(:,3) = calcArea(A3, M2, O) + calcArea(A3, M1, O);
% % For obtuse triangles:
% TriA = calcArea(A1, A2, A3);
% VorRegs(obtuseAngInds, :) = (1/4*ones(sum(obtuseAngInds), 3) + ...
%     1/4*obtuseAngMat(obtuseAngInds, :)).*repmat(TriA(obtuseAngInds), [1 3]);

%% Voronoi Regions per Vertex
M = size(Vertices, 1);
% VorRegsVertices = zeros(M, 1);
VorRegsVertices = sparse(M, M);
for k = 1:M
%     VorRegsVertices(k) = sum(VorRegs(Tri == k)); 
% as I understand - at diagonal areas of Voronois cells
    VorRegsVertices(k, k) = sum(VorRegs(Tri == k)); 
    
    %% UPD 12.08.2012 by Artiom
     VorRegsVertices(k, k) = max(  VorRegsVertices(k, k), 1e-7 );
end

end

%% --------------------------------------------------------------------- %%
function [area_tri] = calcArea(A, B, C)
% Calculate areas of triangles

% Calculate area of each triangle
% area_tri = cross(B - A, C - A, 2);
% area_tri = 1/2*sqrt(sum(area_tri.^2, 2));

area_tri = 1/2*sqrt(sum((B - A).^2, 2).*sum((C - A).^2, 2) - dot(B - A, C - A, 2).^2);

end