function extract_geovec(srcpath, dstpath, geovec_params)

fnames = dir(fullfile(srcpath, '*.mat'));
parfor i = 1 : length(fnames)
    fprintf('Processing %s\n', fnames(i).name)
    tmp = load(fullfile(srcpath, fnames(i).name));
    [desc, ~] = calc_geovec(tmp.Phi, tmp.Lambda, geovec_params);
    parsave(fullfile(dstpath, fnames(i).name), desc);
end
end

function parsave(fn, desc)
save(fn, 'desc', '-v7.3');
end