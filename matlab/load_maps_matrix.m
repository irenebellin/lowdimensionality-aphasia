function [X, voxelIdx, names] = load_maps_matrix(mapDir, pattern, binarizeAt)
% LOAD_MAPS_MATRIX  Read one NIfTI per subject into a subjects x voxels matrix.
%
%   [X, voxelIdx, names] = load_maps_matrix(mapDir, pattern, binarizeAt)
%
%   mapDir     folder with the maps (all in the same space, e.g. MNI 2 mm)
%   pattern    file pattern, e.g. '*MNI.nii.gz'
%   binarizeAt []  keep the values as they are (e.g. disconnection probability)
%              0.5 binarize at this value (value >= binarizeAt becomes 1)
%              For binary lesion masks use 0.5 as well (mask == 1).
%
%   X          n x nUsedVoxels (only voxels altered in at least one subject)
%   voxelIdx   linear indices (into the volume) of the columns of X
%   names      file names, in the order of the rows of X

files = dir(fullfile(mapDir, pattern));
assert(~isempty(files), 'No files matching %s in %s', pattern, mapDir);
n = numel(files);
names = {files.name}';
vox = cell(n, 1);
val = cell(n, 1);
for i = 1:n
    v = niftiread(fullfile(files(i).folder, files(i).name));
    v = double(v);
    if ~isempty(binarizeAt), v = double(v >= binarizeAt); end
    vox{i} = find(v ~= 0);
    val{i} = v(vox{i});
end
voxelIdx = unique(vertcat(vox{:}));
X = zeros(n, numel(voxelIdx));
for i = 1:n
    [~, loc] = ismember(vox{i}, voxelIdx);
    X(i, loc) = val{i};
end
end
