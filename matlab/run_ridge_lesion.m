%% Ridge regression of a behavioral score on lesion maps
% Edit the paths below. Lesion masks must be binary, in the same space
% (MNI152 2 mm), one file per patient. scores.csv needs a column "ID" and a
% column with the score; each lesion file name must START with the patient ID.

lesionDir    = '/path/to/lesions_MNI';
lesionFiles  = '*MNI.nii.gz';
scoresFile   = '/path/to/scores.csv';
scoreColumn  = 'PC1';
templateFile = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask_dil.nii';
outFile      = 'ridge_PC1_lesion.nii.gz';

[X, voxelIdx, names] = load_maps_matrix(lesionDir, lesionFiles, 0.5);

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
save('ridge_lesion_results.mat', 'res', 'names');
