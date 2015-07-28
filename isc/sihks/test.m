load shapes

K = 100;            % number of eigenfunctions
alpha = 2;          % log scalespace basis

T1 = [5:0.5:16];    % time scales for HKS
T2 = [1:0.2:20];    % time scales for SI-HKS
Omega = 2:20;       % frequencies for SI-HKS

for k = 1:2
    % compute cotan Laplacian
    [shape{k}.W shape{k}.A] = mshlp_matrix(shape{k});
    shape{k}.A = spdiags(shape{k}.A,0,size(shape{k}.A,1),size(shape{k}.A,1));

    % compute eigenvectors/values
    [shape{k}.evecs,shape{k}.evals] = eigs(shape{k}.W,shape{k}.A,K,'SM');
    shape{k}.evals = -diag(shape{k}.evals);

    % compute descriptors
    shape{k}.hks   = hks(shape{k}.evecs,shape{k}.evals,alpha.^T1);
    [shape{k}.sihks, shape{k}.schks] = sihks(shape{k}.evecs,shape{k}.evals,alpha,T2,Omega); 
end

i = [959 7106 43365];

% show examples of HKS and SI-HKS at a few points
figure(1)

subplot(1,4,1),
hold on
trisurf(shape{1}.TRIV,shape{1}.X,shape{1}.Y,shape{1}.Z), axis image, shading interp, 
trisurf(shape{2}.TRIV,shape{2}.X+150,shape{2}.Y,shape{2}.Z), axis image, shading interp, 
view([-10 20]), colormap([1 1 1]*0.9), lighting phong, camlight
SYM = {'b.','g.','r.'};
for k = 1:length(i)
    plot3(shape{1}.X(i(k)),shape{1}.Y(i(k)),shape{1}.Z(i(k)),SYM{k})
    plot3(shape{2}.X(i(k))+150,shape{2}.Y(i(k)),shape{2}.Z(i(k)),SYM{k})
end
axis off

subplot(1,4,2), hold on
plot(T1,log(shape{1}.hks(i,:)))
plot(T1,log(shape{2}.hks(i,:)),'--')
title('HKS (log-log)')
xlabel('\tau')
axis square

subplot(1,4,3), hold on
plot(T2,shape{1}.schks(i,:))
plot(T2,shape{2}.schks(i,:),'--')
title('Scale-covariant HKS')
xlabel('\tau')
axis square

subplot(1,4,4), hold on
plot(Omega(:),shape{1}.sihks(i,:))
plot(Omega(:),shape{2}.sihks(i,:),'--')
title('Scale-invariant HKS')
xlabel('\omega')
axis square



%show one component of HKS and SI-HKS at all points
figure(2)

subplot(1,2,1)
hold on
trisurf(shape{1}.TRIV,shape{1}.X,shape{1}.Y,shape{1}.Z, shape{1}.hks(:,20)), axis image, shading interp, 
trisurf(shape{2}.TRIV,shape{2}.X+150,shape{2}.Y,shape{2}.Z, shape{2}.hks(:,20)), axis image, shading interp, 
view([-10 20]),
title('HKS')
axis off

subplot(1,2,2)
hold on
trisurf(shape{1}.TRIV,shape{1}.X,shape{1}.Y,shape{1}.Z, shape{1}.sihks(:,1)), axis image, shading interp, 
trisurf(shape{2}.TRIV,shape{2}.X+150,shape{2}.Y,shape{2}.Z, shape{2}.sihks(:,1)), axis image, shading interp, 
view([-10 20]),
title('SI-HKS')
axis off

