rng(42,'twister')
addpath(genpath('isc'))

%% Compute LBO
nLBO = 300;
extract_lbo('data/train/shapes/', 'data/train/lbo', nLBO);
extract_lbo('data/test/shapes/', 'data/test/lbo', nLBO);

%% Compute GEOVEC
nGEOVEC = 150;
geovec_params = estimate_geovec_params('data/train/lbo', nGEOVEC);
extract_geovec('data/train/lbo', 'data/train/geovec', geovec_params);
extract_geovec('data/test/lbo', 'data/test/geovec', geovec_params);

%% Compute patch operator (disk)
patch_params.rad          = 0.01;    % disk radius
patch_params.flag_dist    = 'fmm';   % possible choices: 'fmm' or 'min'
patch_params.nbinsr       = 5;       % number of rings
patch_params.nbinsth      = 16;      % number of rays
patch_params.fhs          = 0.2;     % factor determining hardness of scale quantization
patch_params.fha          = 0.01;    % factor determining hardness of angle quantization
patch_params.geod_th      = true;
extract_patch_operator('data/train/shapes', 'data/test/disk', patch_params);
extract_patch_operator('data/test/shapes', 'data/test/disk', patch_params);
