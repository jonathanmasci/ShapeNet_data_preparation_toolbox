close all;

rootdir  = fileparts(mfilename('fullpath'));
addpath(fullfile(rootdir,'util'));
addpath(fullfile(rootdir,'sihks'));
addpath(genpath(fullfile('..','..','..','graph_CNN')));

load(fullfile(rootdir,'shapes','0001.scale.1.mat'),'shape')
%name = 'mesh000_fixed_subsampled_3000';
%load(fullfile('~','Dropbox','shape_cnn','dataset','SCAPE',[name,'.mat']));

%%--------------------------------------------------------------------
%% SIHKS signature
%%--------------------------------------------------------------------
fprintf('preprocessing');
shape.idxs    = compute_vertex_face_ring(shape.TRIV');
ndesc         = 1;
[desc,shape] = signature(shape,'sihks');
fprintf('.');

%%--------------------------------------------------------------------
%% ISC settings
%%--------------------------------------------------------------------
rad     = 8;   % radius used for descriptor construction
nbinsr  = 5;    % number of rings
nbinsth = 16;   % number of rays

rr      = [1:nbinsr]/nbinsr*rad;
th      = [1:nbinsth]/nbinsth*2*pi;

fhs     = 2;         %% factor determining hardness of scale quantization
fha     = .01;       %% factors determining hardness of angle quantization

shape.f_dns      = fastmarchmex('init', int32(shape.TRIV-1), double(shape.X(:)), double(shape.Y(:)), double(shape.Z(:)));

[~,vertex] = max(shape.Z);
shape                                 = fast_marching(vertex,shape,'vertex',0,1,shape.f_dns);
[in_ray,in_ring,shp,geod,directions,ds]  = get_net(shape,vertex,'scales',[0,rr],'N_rays',length(th),'fhs',fhs,'fha',fha);

% shape.Av = full(diag(shape.A));

[desc_net,M]                          = get_descriptor_from_net(in_ray,in_ring,desc,shape.Av); 

clear desc_net
tmp_ = reshape( M * desc(:,1), nbinsr, nbinsth );
desc_net(1,:,:) = tmp_;

fastmarchmex('deinit', shape.f_dns);

%% end of code, visualization of results 
a = -82.5000;
v = 18;

dind = 1;
figure; clf;  show_shape(shape,desc(:,dind));
hold on,scatter3(geod{1}(1,1),geod{1}(2,1),geod{1}(3,1),'filled','SizeData',150,'Cdata',[1,0,0])
for k=[1:length(geod)]
    hold on;
    h = plot3(geod{k}(1,:),geod{k}(2,:),geod{k}(3,:));
    set(h,'Color',[0 0 0],'LineWidth',2);
end
for r = [1:length(rr)]
    plot_ring(shape,rr(r));
end
title('Net around vertex')
view(-82.50,18)

figure(2); clf; plot_polarhist(squeeze(desc_net(dind,:,:)),rr,th,0);



