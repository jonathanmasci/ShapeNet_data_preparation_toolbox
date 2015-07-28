function [features,shape] = signature(shape,type)

K = 100;            % number of eigenfunctions
alpha = 2;          % log scalespace basis

T1 = [5:0.5:16];    % time scales for HKS
T2 = [1:0.2:20];    % time scales for SI-HKS
Omega = 2:20;       % frequencies for SI-HKS


% compute cotan Laplacian
if strcmp(computer,'MACI64') || strcmp(computer,'PCWIN64')
	[shape.W shape.A] = mshlp_matrix(shape);
	shape.A = spdiags(shape.A,0,size(shape.A,1),size(shape.A,1));
elseif strcmp(computer,'GLNXA64')
    [shape.W shape.A] = calcLB(shape);
else
    error('[e] load the appropriate mex file for this architecture');
end

% compute eigenvectors/values
try
    [shape.evecs,shape.evals] = eigs(shape.W,shape.A,K,'SM');
catch
    [shape.evecs,shape.evals] = eigs(shape.W,shape.A,K,-1e-8);
end
shape.evals = -diag(shape.evals);

switch type
    case 'sihks'
        features = sihks(shape.evecs,shape.evals,alpha,T2,Omega);
    case 'hks'
        features  = hks(shape.evecs,shape.evals,alpha.^T1);
    case 'wks'
        params_.neig     = 100;
        params_.N        = 100;
        params_.variance = 6;     
        features = wks(shape,params_);
end
