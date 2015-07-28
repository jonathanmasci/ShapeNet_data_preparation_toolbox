function [evecs,evals,area] = calc_lbo(shape, n)

[W, A] = calcLB(shape);
area = diag(A);

[evecs,evals] = eigs(W, A, n, -1e-5, struct('disp', 0));
evals = abs(diag(evals));

[evals, perm] = sort(evals);
evecs = evecs(:, perm);
