function [M, DiagS] = calcLB(shape)
% The L-B operator matrix is computed by DiagS^-1*M.

% Calculate the weights matrix M
M = calcCotMatrixM1([shape.X, shape.Y, shape.Z], shape.TRIV);

% Calculate the diagonal of matrix S 
DiagS = calcVoronoiRegsCircCent(shape.TRIV, [shape.X, shape.Y, shape.Z]);
%%
DiagS = abs( DiagS );
%%

end


% ----------------------------------------------------------------------- %
function [M] = calcCotMatrixM1(Vertices, Tri)

N = size(Vertices, 1);
M = sparse(N, N);

v1 = Vertices(Tri(:, 2), :) - Vertices(Tri(:, 1), :); %v1 = v1./repmat(normVec(v1), 1, 3);
v2 = Vertices(Tri(:, 3), :) - Vertices(Tri(:, 1), :); %v2 = v2./repmat(normVec(v2), 1, 3);
v3 = Vertices(Tri(:, 3), :) - Vertices(Tri(:, 2), :); %v3 = v3./repmat(normVec(v3), 1, 3);

% cot1 = dot( v1, v2, 2)./normVec(cross( v1, v2, 2)); %cot1(cot1 < 0) = 0;
% cot2 = dot(-v1, v3, 2)./normVec(cross(-v1, v3, 2)); %cot2(cot2 < 0) = 0;
% cot3 = dot(-v2, -v3, 2)./normVec(cross(-v2, -v3, 2)); %cot3(cot3 < 0) = 0;
tmp1 = dot( v1,  v2, 2); cot1 = tmp1./sqrt(normVec(v1).^2.*normVec(v2).^2 - (tmp1).^2); clear tmp1;
tmp2 = dot(-v1,  v3, 2); cot2 = tmp2./sqrt(normVec(v1).^2.*normVec(v3).^2 - (tmp2).^2); clear tmp2;
tmp3 = dot(-v2, -v3, 2); cot3 = tmp3./sqrt(normVec(v2).^2.*normVec(v3).^2 - (tmp3).^2); clear tmp3;

for k = 1:size(Tri, 1)
    M(Tri(k, 1), Tri(k, 2)) = M(Tri(k, 1), Tri(k, 2)) + cot3(k);
    M(Tri(k, 1), Tri(k, 3)) = M(Tri(k, 1), Tri(k, 3)) + cot2(k);
    M(Tri(k, 2), Tri(k, 3)) = M(Tri(k, 2), Tri(k, 3)) + cot1(k);
end
M = 0.5*(M + M'); % here she does the normalization (comment - Artiom)
    
% inds = sub2ind([N, N], [Tri(:, 2); Tri(:, 1); Tri(:, 1)], [Tri(:, 3); Tri(:, 3); Tri(:, 2)]);
% M(inds) = M(inds) + [cot1; cot2; cot3];
% inds = sub2ind([N, N], [Tri(:, 3); Tri(:, 3); Tri(:, 2)], [Tri(:, 2); Tri(:, 1); Tri(:, 1)]);
% M(inds) = M(inds) + [cot1; cot2; cot3];
% M = 0.5*(M + M');
% % M(M < 0) = 0;

M = M - diag(sum(M, 2)); % making it Laplacian

    function normV = normVec(vec)
        normV = sqrt(sum(vec.^2, 2));
    end
%     function normalV = normalizeVec(vec)
%         normalV = vec./repmat(normVec(vec), 1, 3);
%     end

end

% ----------------------------------------------------------------------- %
function [M] = calcCotMatrixM(Vertices, Tri) %#ok<DEFNU>

N = size(Vertices, 1);
[transmat] = calcTransmat(N, Tri);

% Calculate the matrix M, when {M}_ij = (cot(alpha_ij) + cot(beta_ij))/2
% [transrow, transcol] = find(triu(transmat,1) > 0);
[transrow, transcol] = find((triu(transmat,1) > 0) | (triu(transmat',1) > 0));
M = sparse(N, N);

for k = 1:length(transrow)
    
    P = transrow(k);
    Q = transcol(k);
    S = transmat(P,Q);
    R = transmat(Q,P);
    
    %%
%     u1 = Vertices(Q, :) - Vertices(R, :); u1 = u1./norm(u1);
%     v1 = Vertices(P, :) - Vertices(R, :); v1 = v1./norm(v1);
%     u2 = Vertices(P, :) - Vertices(S, :); u2 = u2./norm(u2);
%     v2 = Vertices(Q, :) - Vertices(S, :); v2 = v2./norm(v2);
%     M(P,Q) = -1/2*(dot(u1, v1)/norm(cross(u1, v1)) + dot(u2, v2)/norm(cross(u2, v2)));

    tmp1 = 0;
    tmp2 = 0;
    
    if (R ~= 0)
        u1 = Vertices(Q, :) - Vertices(R, :); u1 = u1./norm(u1);
        v1 = Vertices(P, :) - Vertices(R, :); v1 = v1./norm(v1);
        tmp1 = dot(u1, v1)/norm(cross(u1, v1));
    end

    if (S ~= 0)
        u2 = Vertices(P, :) - Vertices(S, :); u2 = u2./norm(u2);
        v2 = Vertices(Q, :) - Vertices(S, :); v2 = v2./norm(v2);
        tmp2 = dot(u2, v2)/norm(cross(u2, v2));
    end
    
    M(P,Q) = -1/2*(tmp1 + tmp2);
    %%
    
end

M = 0.5*(M + M');
M = M - diag(sum(M, 2));

end


% ----------------------------------------------------------------------- %
function [transmat] = calcTransmat(N, Tri)

% Calculation of the map of all the connected vertices: for each i,j,
% transmat(i,j) equals to the third vertex of the triangle which connectes
% them; if the vertices aren't connected - transmat(i,j) = 0.
transmat = sparse(N, N);
transmat(sub2ind(size(transmat), Tri(:,1), Tri(:,2))) = Tri(:,3);
transmat(sub2ind(size(transmat), Tri(:,2), Tri(:,3))) = Tri(:,1);
transmat(sub2ind(size(transmat), Tri(:,3), Tri(:,1))) = Tri(:,2);

end