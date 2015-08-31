function extract_geovec(srcpath, dstpath, geovec_params)

if ~exist(dstpath, 'dir')
    mkdir(dstpath);
end

fnames = dir(fullfile(srcpath, '*.mat'));
parfor i = 1 : length(fnames)
    if exist(fullfile(dstpath, fnames(i).name), 'file')
        fprintf('%s already processed, skipping\n', fnames(i).name)
        continue
    end
    fprintf('Processing %s\n', fnames(i).name)
    tmp = load(fullfile(srcpath, fnames(i).name));
    [desc, ~] = calc_geovec(tmp.Phi, tmp.Lambda, geovec_params);
    parsave(fullfile(dstpath, fnames(i).name), desc);
end
end

function parsave(fn, desc)
save(fn, 'desc', '-v7.3');
end