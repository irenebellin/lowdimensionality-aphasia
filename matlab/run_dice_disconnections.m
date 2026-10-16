%% Dice analysis, disconnections (Figure 5B)
% Edit the paths. groups.csv (kept private) needs two columns: ID (e.g. DIS_006) and
% syndrome (one of cfg.syndromeNames).

cfg.mode          = 'disconnection';
cfg.null          = 'random_tracts';       % or 'shuffle_labels'
cfg.nPerm         = 1000;
cfg.pc1File       = '/path/to/ridge_regression_disconnection_PC1.nii';
cfg.brainMaskFile = '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz';
cfg.syndromeNames = {'Motor Transcortical', 'Wernicke', 'Broca', 'Conduction', 'Dyslexia', ...
                     'Anomic', 'Global', 'Numerical System', 'Mixed Transcortical'};

groupMapFolder    = '/path/to/group_disconnection_maps/';   % real group maps, one per syndrome
cfg.syndromeMapFiles = fullfile(groupMapFolder, {'transcorticale_motoria.nii', 'wernicke.nii', ...
    'broca.nii', 'conduzione.nii', 'dislessia.nii', 'anomica.nii', 'globale.nii', ...
    'numero.nii', 'transcorticale_mista.nii'});

cfg.probabilityCsv = '/path/to/results_DIS_bcb/probability.csv';   % BCBtoolkit output
cfg.tractFolder    = '/path/to/BCBToolKit/Tracts';
cfg.idSuffix       = '_lesion_MNI';

G = readtable('/path/to/groups.csv');
cfg.subjectIds      = cell(1, numel(cfg.syndromeNames));
cfg.subjectMapFiles = cell(1, numel(cfg.syndromeNames));      % only needed for 'shuffle_labels'
disconFolder = '/path/to/results_DIS_bcb/';
for i = 1:numel(cfg.syndromeNames)
    ids = cellstr(string(G.ID(strcmp(string(G.syndrome), cfg.syndromeNames{i}))));
    cfg.subjectIds{i}      = ids;
    cfg.subjectMapFiles{i} = fullfile(disconFolder, strcat(ids, '_lesion_MNI.nii'));
end

res = dice_permutation_analysis(cfg);
save('dice_disconnection_results.mat', 'res');
writetable(res.table, 'dice_disconnection_results.csv');
plot_dice_results(res, [-0.02 0.6]);
