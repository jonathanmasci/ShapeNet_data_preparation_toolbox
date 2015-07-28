function extract_lbo(srcpath, dstpath, nLBO)

fnames = dir(fullfile(srcpath, '*.mat'));
parfor i = 1 : length(fnames)
    fprintf('Processing %s\n', fnames(i).name)
    tmp = load(fullfile(srcpath, fnames(i).name));
    [Phi, Lambda, A] = calc_lbo(tmp.shape, nLBO);
    parsave(fullfile(dstpath, fnames(i).name), Phi, Lambda, A);
end
end

function parsave(fn, Phi, Lambda, A)
save(fn, 'Phi', 'Lambda', 'A', '-v7.3')
end