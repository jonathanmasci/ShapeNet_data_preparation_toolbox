% Computes heat kernel signature H_t(x,x), where H_t(x,y) is the heat kernel
%
% Usage:  desc = hks(evecs,evals,T)
%
% Input:  evecs   - (n x k) Laplace-Beltrami eigenvectors arranged as columns 
%         evals   - (k x 1) corresponding Laplace-Beltrami eigenvalues 
%         T       - (1 x t) time values
%
% Output: desc    - (n x t) matrix with the values of HKS for different t as columns
%
% (c) Michael Bronstein 2012    http://www.inf.usi.ch/bronstein
%
% M. M. Bronstein, I. Kokkinos, "Scale-invariant heat kernel signatures for non-rigid 
% shape recognition", CVPR 2010. 

function desc = hks(evecs,evals,T)

desc = zeros(size(evecs,1),length(T));

for k = 1:length(T)
    desc(:,k) = sum(evecs.^2 .* repmat(exp(-T(k)*evals(:))',[size(evecs,1) 1]),2);
end


