function res = dice_permutation_analysis(cfg)
% DICE_PERMUTATION_ANALYSIS  Overlap (Dice) between the PC1-related map and the
% group lesion / disconnection map of each aphasia syndrome, tested against a
% null distribution with one-sided permutation tests.
%
%   res = dice_permutation_analysis(cfg)
%
% For each syndrome:
%   1. The PC1-related map is thresholded to the same number of voxels as the real
%      group map (top-ranking voxels), and the observed Dice coefficient between
%      the two masks is computed.
%   2. A null distribution of Dice values is generated (cfg.nPerm iterations):
%        'shuffle_voxels' (lesion mode): every patient's lesion is moved to random
%             voxels of the brain (lesion size kept), the simulated group map is
%             thresholded to the same number of voxels, and compared with the
%             PC1 mask.
%        'random_tracts' (disconnection mode): a random set of tracts, as many as
%             the unique tracts disconnected (probability > 0.5) in the group, is
%             compared with the real group disconnection map. The size of the
%             random mask is NOT matched to the real map.
%        'shuffle_labels' (both modes): syndrome labels are shuffled across all
%             patients, the maps of the patients assigned to the syndrome are
%             summed, thresholded to the same number of voxels, and compared with
%             the PC1 mask.
%   3. One-sided permutation p-value: p = (1 + #{null >= observed}) / (1 + nPerm),
%      with Holm correction across syndromes.
%
% cfg fields (required)
%   .mode             'lesion' | 'disconnection'
%   .pc1File          NIfTI of the PC1-related map (ridge regression weights)
%   .brainMaskFile    NIfTI brain mask in the same space as all maps
%   .syndromeNames    1 x S cellstr
%   .syndromeMapFiles 1 x S cellstr, real group maps (lesion or disconnection)
%   .subjectMapFiles  1 x S cell; each a cellstr with the maps of the patients of
%                     the syndrome (needed for 'shuffle_voxels' and 'shuffle_labels')
%   for 'random_tracts' also:
%   .probabilityCsv   BCBtoolkit table: first column = patient file name, other
%                     columns = probability of disconnection of each tract
%   .tractFolder      folder with the tract maps (<tract name>.nii.gz)
%   .subjectIds       1 x S cell; each a cellstr with patient identifiers
% cfg fields (optional; default)
%   .null             'shuffle_voxels' (lesion) | 'random_tracts' (disconnection)
%   .nPerm            1000
%   .seed             1
%   .rankBy           'signed' | 'absolute'                  ('signed')
%                     'signed' ranks voxels by their value (original analysis);
%                     'absolute' ranks by absolute value.
%   .tieBreak         'index' | 'random'                     ('index')
%                     voxels with equal values are ranked by position ('index',
%                     original analysis) or at random ('random').
%   .tractThreshold   0.5
%   .idSuffix         '_lesion_MNI'   appended to the identifier to find the patient
%                                     in the probability table
%
% Output (struct): dice_real, dice_nulls (nPerm x S), pUpper, pLower, pHolm,
%   syndromeNames, cfg.

cfg = set_defaults(cfg);
rng(cfg.seed);
S = numel(cfg.syndromeNames);

%% Brain space: all maps are handled as vectors over the voxels of the brain mask
mask = niftiread(cfg.brainMaskFile) > 0;
idx  = find(mask);
nB   = numel(idx);

pc1 = double(niftiread(cfg.pc1File));
assert(isequal(size(pc1), size(mask)), 'PC1 map and brain mask must have the same size.');
pc1 = pc1(idx);
fprintf('PC1 map: %d voxels with negative values, %d with positive values.\n', nnz(pc1 < 0), nnz(pc1 > 0));
if nnz(pc1 < 0) > 0 && strcmp(cfg.rankBy, 'signed')
    warning(['The PC1 map has negative values and voxels are ranked by their signed value: ' ...
             'only the highest positive values are selected. Use rankBy = ''absolute'' to rank by absolute value.']);
end
rankVal = pc1;
if strcmp(cfg.rankBy, 'absolute'), rankVal = abs(pc1); end

%% Patient maps (loaded once)
needSubjects = any(strcmp(cfg.null, {'shuffle_voxels', 'shuffle_labels'}));
if needSubjects
    subjectMaps = cell(1, S);
    for g = 1:S
        files = cfg.subjectMapFiles{g};
        subjectMaps{g} = false(numel(files), nB);
        for s = 1:numel(files)
            v = niftiread(files{s}) > 0;
            subjectMaps{g}(s, :) = v(idx)';
        end
    end
    if strcmp(cfg.null, 'shuffle_labels')
        allMaps   = vertcat(subjectMaps{:});
        allLabels = repelem(1:S, cellfun(@(m) size(m, 1), subjectMaps))';
    end
end

%% Tract masks (loaded once) for the disconnection null
if strcmp(cfg.null, 'random_tracts')
    prob       = readtable(cfg.probabilityCsv);
    tractNames = prob.Properties.VariableNames(2:end);
    nT = numel(tractNames);
    tractMasks = false(nT, nB);
    for t = 1:nT
        f = fullfile(cfg.tractFolder, [tractNames{t} '.nii.gz']);
        if ~isfile(f), warning('Tract file not found: %s', f); continue; end
        v = niftiread(f) > cfg.tractThreshold;
        tractMasks(t, :) = v(idx)';
    end
    probIds = string(prob{:, 1});
end

%% Main loop over syndromes
dice = @(a, b) 2 * nnz(a & b) / (nnz(a) + nnz(b));
res.dice_real  = nan(S, 1);
res.dice_nulls = nan(cfg.nPerm, S);

for i = 1:S
    a = niftiread(cfg.syndromeMapFiles{i}) > 0;
    a = a(idx);                                   % real group map
    N = nnz(a);
    if N == 0
        warning('Empty group map for %s', cfg.syndromeNames{i});
        continue;
    end
    m = top_n(rankVal, N, cfg.tieBreak);          % PC1 mask, same size as the group map
    res.dice_real(i) = dice(m, a);

    nullVals = zeros(cfg.nPerm, 1);
    switch cfg.null
        case 'shuffle_voxels'
            L = subjectMaps{i};
            sizes = sum(L, 2);
            for p = 1:cfg.nPerm
                freq = zeros(nB, 1);
                for s = 1:numel(sizes)
                    ix = randperm(nB, sizes(s));
                    freq(ix) = freq(ix) + 1;
                end
                nullVals(p) = dice(top_n(freq, N, cfg.tieBreak), m);
            end
        case 'shuffle_labels'
            for p = 1:cfg.nPerm
                sel  = allLabels(randperm(numel(allLabels))) == i;
                freq = sum(allMaps(sel, :), 1)';
                nullVals(p) = dice(top_n(freq, N, cfg.tieBreak), m);
            end
        case 'random_tracts'
            ids = string(cfg.subjectIds{i}(:)) + cfg.idSuffix;
            used = [];
            for s = 1:numel(ids)
                r = find(probIds == ids(s));
                if isempty(r), warning('Patient %s not found in the table', ids(s)); continue; end
                used = union(used, find(prob{r, 2:end} > cfg.tractThreshold));
            end
            k = numel(used);
            fprintf('%s: %d unique tracts above %.1f\n', cfg.syndromeNames{i}, k, cfg.tractThreshold);
            for p = 1:cfg.nPerm
                rm = any(tractMasks(randperm(nT, k), :), 1)';
                nullVals(p) = dice(rm, a);
            end
    end
    res.dice_nulls(:, i) = nullVals;
    fprintf('%s: Dice = %.4f | permutation p = %.4f\n', cfg.syndromeNames{i}, ...
        res.dice_real(i), (1 + sum(nullVals >= res.dice_real(i))) / (1 + cfg.nPerm));
end

%% Permutation p-values (one-sided: overlap greater than chance) and Holm correction
res.pUpper = ((1 + sum(res.dice_nulls >= res.dice_real', 1)) / (1 + cfg.nPerm))';
res.pLower = ((1 + sum(res.dice_nulls <= res.dice_real', 1)) / (1 + cfg.nPerm))';
res.pHolm  = holm(res.pUpper);
res.syndromeNames = cfg.syndromeNames;
res.cfg = cfg;

T = table(cfg.syndromeNames(:), res.dice_real, mean(res.dice_nulls, 1)', std(res.dice_nulls, 0, 1)', ...
    res.pUpper, res.pHolm, 'VariableNames', {'syndrome', 'dice', 'null_mean', 'null_sd', 'p_upper', 'p_holm'});
disp(T);
res.table = T;
end

%% ----------------------------------------------------------------- helpers
function cfg = set_defaults(cfg)
if ~isfield(cfg, 'null') || isempty(cfg.null)
    if strcmp(cfg.mode, 'lesion'), cfg.null = 'shuffle_voxels'; else, cfg.null = 'random_tracts'; end
end
def = struct('nPerm', 1000, 'seed', 1, 'rankBy', 'signed', 'tieBreak', 'index', ...
             'tractThreshold', 0.5, 'idSuffix', '_lesion_MNI');
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(cfg, f{i}) || isempty(cfg.(f{i})), cfg.(f{i}) = def.(f{i}); end
end
end

function m = top_n(values, N, tieBreak)
% Logical mask of the N highest values (ties broken by position or at random).
values = values(:);
if strcmp(tieBreak, 'random')
    span = max(values) - min(values);
    if span == 0, span = 1; end
    values = values + 1e-6 * span * rand(size(values));
end
[~, ix] = sort(values, 'descend');
m = false(size(values));
m(ix(1:N)) = true;
end

function adj = holm(p)
[ps, ix] = sort(p(:));
mm = numel(ps);
a = min(1, cummax(ps .* (mm:-1:1)'));
adj = zeros(mm, 1);
adj(ix) = a;
end
