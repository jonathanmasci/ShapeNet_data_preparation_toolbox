function [Vray,TRIV_] = get1ring(shape,i,tt,K,visualize)
if nargin==4,
    visualize = 0;
end

XX = [shape.X shape.Y shape.Z];

TRIV_ = find_adj_triang(shape.TRIV,tt,i);
% rng(42,'twister')
TRIV_ = circshift(TRIV_,randi(numel(tt)));% do random re-order of 1-ring, so directions won't depent on it

V1 = XX(TRIV_(:,2),:) - XX(TRIV_(:,1),:);
V2 = XX(TRIV_(:,3),:) - XX(TRIV_(:,1),:);

V1n = V1 ./ repmat(sqrt(sum(V1.^2,2)),[1 3]);
V2n = V2 ./ repmat(sqrt(sum(V2.^2,2)),[1 3]);


tang = acos(dot(V1n,V2n,2));
ttang = sum(tang);
pang = cumsum(tang/ttang*2*pi);
subang = [0; pang(1:end-1)];

%% <Iasonas slack>
sang = [0:K-1]'/K*2*pi + .5/K*2*pi;
%% </ Iasonas slack>
tray = sum(~(repmat(sang,[1 length(pang)]) <= repmat(pang',[length(sang) 1])),2)+1;

offang = (sang - subang(tray));
alpha = offang ./ (tang(tray)/ttang*2*pi);

% rotation angles
rotang  = alpha.*tang(tray);
rotaxis = cross(V1n,V2n);
rotaxis = rotaxis(tray,:);

V1   = V1./repmat(sqrt(sum(V1.^2,2)),[1,3]);
Vray = zeros(K,3);
for t = 1:K
    Vray(t,:) = (rotmtx(rotaxis(t,:),rotang(t))*V1(tray(t),:)')';
end

directions_0 = Vray; 
dist_tips = dist2(directions_0,directions_0);
dist_tips = dist_tips + 1e10*eye(size(dist_tips)); 

ordered     = [];
ordered(1)  = 1;
for k=1:K-1,
    [m,ordered(k+1)] = min(dist_tips(ordered(k),:));
    dist_tips(ordered,ordered(k+1)) = 1e10;
    dist_tips(ordered(k+1),ordered) = 1e10;
end
Vray = directions_0(ordered,:)';
if visualize 
    visualize_ring(shape,i,tt,Vray)
end

function R = rotmtx(u,th)

ux = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];
ut = repmat(u(:),[1 3]).*repmat(u(:),[1 3])';

R = eye(3)*cos(th) + ux*sin(th) + ut*(1-cos(th));

function [TRIV_,ADJ] = find_adj_triang(TRIV,tt,i)

% reorder vertices in triangle
%tt = find(TRIV(:,1)==i | TRIV(:,2)==i | TRIV(:,3)==i);
TRIV_ = TRIV(tt,:);
[~,ii] = sort(TRIV_~=i,2);
TRIV_ = TRIV_(sub2ind(size(TRIV_),repmat([1:size(TRIV_,1)]',[1 3]),ii));

% reorder triangles
ADJ = spalloc(size(TRIV_,1),size(TRIV_,1),2*size(TRIV_,1));
for t = 1:size(TRIV_,1)
    ADJ(t,find(TRIV_(:,2)==TRIV_(t,2) | TRIV_(:,3)==TRIV_(t,2) | TRIV_(:,2)==TRIV_(t,3) | TRIV_(:,3)==TRIV_(t,3))) = 1;
end
ADJ = ADJ - diag(diag(ADJ));

%%%%%%%%%%%%%%%%%%%%%%
% treat boundaries !!!
%%%%%%%%%%%%%%%%%%%%%%
try
    
t = [1];
idx = setdiff(1:size(TRIV_,1),t);
while ~isempty(idx)
    ic = setdiff(find(ADJ(t(end),:)),t);
    t = [t ic(1)];
    idx = setdiff(1:size(TRIV_,1),t);
end

TRIV_ = TRIV_(t,:);

% reorder orientations
if size(TRIV_,1)>1
    
    for t = 2:size(TRIV_,1)
        if find(TRIV_(t,2:3)==TRIV_(t-1,end))==2
            TRIV_(t,[2 3]) = TRIV_(t,[3 2]);
        end
    end
    
end

catch
    
    
end
