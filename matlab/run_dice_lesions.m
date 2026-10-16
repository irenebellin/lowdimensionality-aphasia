%% Dice analysis, lesions (Figure 5A)
% Edit the paths. groups.csv (kept private, it links patients to syndromes) needs two
% columns: ID (e.g. DIS_006) and syndrome (one of cfg.syndromeNames).

cfg.mode          = 'lesion';
cfg.null          = 'shuffle_voxels';      % or 'shuffle_labels'
cfg.nPerm         = 1000;
cfg.pc1File       = '/path/to/ridge_regression_lesion_PC1.nii';
cfg.brainMaskFile = '/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz';
cfg.syndromeNames = {'Motor Transcortical', 'Wernicke', 'Broca', 'Conduction', 'Dyslexia', ...
                     'Anomic', 'Global', 'Numerical System', 'Mixed Transcortical'};

groupMapFolder    = '/path/to/group_lesion_maps/';     % real group maps, one per syndrome
cfg.syndromeMapFiles = fullfile(groupMapFolder, {'r_transcorticale_motoria.nii', 'r_wernicke.nii', ...
    'r_broca.nii', 'r_conduzione.nii', 'r_dislessia.nii', 'r_anomica.nii', 'r_globale.nii', ...
    'r_numero.nii', 'r_transcorticale_mista.nii'});

lesionFolder = '/path/to/individual_lesions_MNI_resampled/';
G = readtable('/path/to/groups.csv');
cfg.subjectMapFiles = cell(1, numel(cfg.syndromeNames));
for i = 1:numel(cfg.syndromeNames)
    ids = cellstr(string(G.ID(strcmp(string(G.syndrome), cfg.syndromeNames{i}))));
    cfg.subjectMapFiles{i} = fullfile(lesionFolder, strcat('r_', ids, '_lesion_MNI.nii'));
end

res = dice_permutation_analysis(cfg);
save('dice_lesion_results.mat', 'res');
writetable(res.table, 'dice_lesion_results.csv');
plot_dice_results(res, [-0.02 0.4]);
