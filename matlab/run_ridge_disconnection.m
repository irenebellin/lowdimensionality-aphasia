%% Ridge regression of a behavioral score on disconnection maps
% Same pipeline as run_ridge_lesion.m. Disconnection maps (e.g. BCBtoolkit
% probability maps) must be in the same space, one file per patient.
%
% Set binarizeAt = [] to use the probability values as they are, or a value
% (e.g. 0.5) to binarize. Use the same choice as in the original analysis.

disconDir    = '/path/to/disconnection_maps_MNI';
disconFiles  = '*disconnection.nii.gz';
binarizeAt   = [];
scoresFile   = '/path/to/scores.csv';
scoreColumn  = 'PC1';
templateFile = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask_dil.nii';
outFile      = 'ridge_PC1_disconnection.nii.gz';

[X, voxelIdx, names] = load_maps_matrix(disconDir, disconFiles, binarizeAt);

T  = readtable(scoresFile);
ID = string(T.ID);
y  = nan(numel(names), 1);
for i = 1:numel(names)
    j = find(startsWith(string(names{i}), ID));
    assert(numel(j) == 1, 'File %s does not match exactly one ID.', names{i});
    y(i) = T.(scoreColumn)(j);
end
assert(~any(isnan(y)), 'Missing scores.');

opts = struct('correction', 'maxT');   % 'none' reproduces the uncorrected analysis
res  = ridge_pca_permutation(X, y, opts);

weights_to_nifti(res.voxelWeights, voxelIdx, templateFile, outFile);
save('ridge_disconnection_results.mat', 'res', 'names');
