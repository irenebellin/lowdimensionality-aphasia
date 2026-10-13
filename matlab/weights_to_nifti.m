function vol = weights_to_nifti(w, voxelIdx, templateFile, outFile, weightThreshold)
% WEIGHTS_TO_NIFTI  Map voxel weights back to the brain, normalize, threshold, smooth, save.
%
%   weights_to_nifti(w, voxelIdx, templateFile, outFile, weightThreshold)
%
%   Steps (as in Supplementary S1): scale to [-1, 1], set |value| < threshold
%   to zero (default 0.05), then Gaussian smoothing (smooth3, 5-voxel window).

if nargin < 5 || isempty(weightThreshold), weightThreshold = 0.05; end
ref  = niftiread(templateFile);
info = niftiinfo(templateFile);

vol = zeros(size(ref));
vol(voxelIdx) = w;
vol = vol ./ max(abs(vol(:)));
vol(abs(vol) < weightThreshold) = 0;
vol = smooth3(vol, 'gaussian', 5);

vol = single(vol);
info.Datatype = 'single';
info.BitsPerPixel = 32;
info.MultiplicativeScaling = 1;
info.AdditiveOffset = 0;

[folder, base] = fileparts(outFile);
base = regexprep(base, '\.nii$', '');
niftiwrite(vol, fullfile(folder, base), info, 'Compressed', true);
end
