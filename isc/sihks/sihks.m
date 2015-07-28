% Computes scale-covariant and scale-invariant heat kernel signature (SI-HKS)
%
% Usage:  [sc,si] = sihks(evecs,evals,T,Omega)
%
% Input:  evecs   - (n x k) Laplace-Beltrami eigenvectors arranged as columns 
%         evals   - (k x 1) corresponding Laplace-Beltrami eigenvalues 
%         alpha   - log scalespace basis
%         T       - (1 x t) time values (used as alpha.^T)
%         Omega   - (1 x w) frequency samples
%
% Output: sc      - (n x t) matrix with the values of the scale-covariant HKS 
%                   for different t as columns
%         si      - (n x w) matrix with the values of the scale-invariant HKS 
%                   for different omega as columns
%
% (c) Michael Bronstein 2012    http://www.inf.usi.ch/bronstein
%
% M. M. Bronstein, I. Kokkinos, "Scale-invariant heat kernel signatures for non-rigid 
% shape recognition", CVPR 2010. 

function [si sc] = sihks(evecs,evals,alpha,T,Omega)

sc = zeros(size(evecs,1),length(T));

for k = 1:length(T)
    sc(:,k) =  - log(alpha)*sum(evecs.^2 .* repmat(alpha^T(k)*evals(:)'.*exp(-alpha^T(k)*evals(:))',[size(evecs,1) 1]),2) ./ ...
        sum(evecs.^2 .* repmat(exp(-alpha^T(k)*evals(:))',[size(evecs,1) 1]),2);
end

si = abs(fft(sc,[],2));
si = si(:,Omega);
