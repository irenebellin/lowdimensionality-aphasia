function res = ridge_pca_permutation(X, y, opts)
% RIDGE_PCA_PERMUTATION  PCA + ridge regression with leave-one-out
% cross-validation and permutation testing.
%
%   res = ridge_pca_permutation(X, y, opts)
%
% Inputs
%   X     n x p matrix (subjects x voxels): binary lesion masks, or
%         disconnection maps, restricted to voxels altered in >= 1 subject.
%   y     n x 1 behavioral score (e.g. PC1 of the OCS/NIHSS scores).
%   opts  (optional) struct with fields
%         .varianceKept    cumulative variance of the PCA to keep    (0.95)
%         .lambdas         ridge penalties to try          (logspace(-2,4,200))
%         .nPerm           number of permutations                    (10000)
%         .alpha           two-sided alpha                           (0.05)
%         .correction      'none' | 'maxT' | 'fdr'                   ('maxT')
%                          'none' : each component is compared with its own
%                                   2.5th/97.5th permutation percentiles
%                                   (as in the original analysis, uncorrected)
%                          'maxT' : family-wise error control across components
%                                   (|beta| compared with the (1-alpha) percentile
%                                   of the maximum |beta| over components)
%                          'fdr'  : Benjamini-Hochberg on per-component
%                                   permutation p-values
%         .backprojection  'original' | 'unstandardized'            ('original')
%                          'original'       : voxel weights = loadings * beta
%                          'unstandardized' : betas divided by the SD of each
%                                   component score before back-projection
%         .seed            random seed                                    (1)
%
% Output (struct)
%   lambda, r2Loo (leave-one-out cross-validated R2 at the chosen lambda),
%   r2Grid, pR2 (permutation p-value of the cross-validated R2), beta,
%   pComp, pAdj, sig (logical, significant components), nSig, nComp,
%   explained, voxelWeights (p x 1), yPredLoo.
%
% Notes
%   - Leave-one-out predictions are computed in closed form from the hat
%     matrix of the ridge model; they are identical to refitting the model n
%     times, without an intercept, on the z-scored component scores.
%   - Regression coefficients are those of the model fitted on all subjects.
%   - lambda is selected with the same cross-validation that is then reported,
%     so r2Loo is slightly optimistic (nested cross-validation would remove this).

if nargin < 3, opts = struct(); end
opts = set_defaults(opts);
rng(opts.seed);

y = y(:);
n = size(X, 1);
assert(numel(y) == n, 'X and y must contain the same number of subjects.');

%% 1) PCA on the imaging matrix
[coeff, score, ~, ~, explained] = pca(X, 'Centered', true);
k = find(cumsum(explained) >= 100 * opts.varianceKept, 1, 'first');
k = min(k, size(score, 2));
compSD = std(score(:, 1:k));
assert(all(compSD > 0), 'A retained component has zero variance.');
Z  = zscore(score(:, 1:k));
yz = zscore(y);

%% 2) Ridge regression: leave-one-out cross-validation over lambda
[U, S, V] = svd(Z, 'econ');
s   = diag(S);
Uty = U' * yz;
nL  = numel(opts.lambdas);
r2Grid = zeros(nL, 1);
for l = 1:nL
    lam  = opts.lambdas(l);
    d    = s.^2 ./ (s.^2 + lam);
    yhat = U * (d .* Uty);
    h    = sum((U.^2) .* d', 2);                 % diagonal of the hat matrix
    loo  = yz - (yz - yhat) ./ (1 - h);          % leave-one-out prediction
    r2Grid(l) = 1 - sum((yz - loo).^2) / sum((yz - mean(yz)).^2);
end
[r2Loo, iBest] = max(r2Grid);
lambda = opts.lambdas(iBest);
d    = s.^2 ./ (s.^2 + lambda);
h    = sum((U.^2) .* d', 2);
yPredLoo = yz - (yz - U * (d .* Uty)) ./ (1 - h);
beta = V * ((s ./ (s.^2 + lambda)) .* Uty);

%% 3) Permutation test (behavioral scores shuffled across subjects)
betaNull = zeros(k, opts.nPerm);
r2Null   = zeros(opts.nPerm, 1);
for p = 1:opts.nPerm
    yp  = yz(randperm(n));
    Uyp = U' * yp;
    betaNull(:, p) = V * ((s ./ (s.^2 + lambda)) .* Uyp);
    loo = yp - (yp - U * (d .* Uyp)) ./ (1 - h);
    r2Null(p) = 1 - sum((yp - loo).^2) / sum((yp - mean(yp)).^2);
end
pR2   = (1 + sum(r2Null >= r2Loo)) / (1 + opts.nPerm);
pComp = (1 + sum(abs(betaNull) >= abs(beta), 2)) / (1 + opts.nPerm);

%% 4) Which components are significant (with or without correction)
switch lower(opts.correction)
    case 'none'
        lo  = prctile(betaNull, 100 * opts.alpha / 2, 2);
        hi  = prctile(betaNull, 100 * (1 - opts.alpha / 2), 2);
        sig  = (beta < lo) | (beta > hi);
        pAdj = pComp;
    case 'maxt'
        maxNull = max(abs(betaNull), [], 1);
        sig  = abs(beta) > prctile(maxNull, 100 * (1 - opts.alpha));
        pAdj = (1 + sum(maxNull >= abs(beta), 2)) / (1 + opts.nPerm);
    case 'fdr'
        [ps, ix] = sort(pComp);
        m = numel(ps);
        below = find(ps <= (1:m)' / m * opts.alpha, 1, 'last');
        sig = false(k, 1);
        if ~isempty(below), sig(ix(1:below)) = true; end
        pAdj = nan(k, 1);
        adj  = min(1, cummin(ps(end:-1:1) .* m ./ (m:-1:1)'));
        pAdj(ix) = adj(end:-1:1);
    otherwise
        error('opts.correction must be ''none'', ''maxT'' or ''fdr''.');
end

%% 5) Back-projection of the significant components to voxel space
betaSig = beta;
betaSig(~sig) = 0;
if strcmpi(opts.backprojection, 'unstandardized')
    betaSig = betaSig ./ compSD(:);
end
voxelWeights = coeff(:, 1:k) * betaSig;

res = struct('lambda', lambda, 'r2Loo', r2Loo, 'r2Grid', r2Grid, 'pR2', pR2, ...
    'beta', beta, 'pComp', pComp, 'pAdj', pAdj, 'sig', sig, 'nSig', sum(sig), ...
    'nComp', k, 'explained', explained(1:k), 'voxelWeights', voxelWeights, ...
    'yPredLoo', yPredLoo, 'opts', opts);

fprintf('Components kept: %d | lambda = %.4g | LOOCV R2 = %.3f | permutation p = %.4f\n', ...
    k, lambda, r2Loo, pR2);
fprintf('Significant components (%s): %d of %d\n', opts.correction, sum(sig), k);
end

function opts = set_defaults(opts)
def = struct('varianceKept', 0.95, 'lambdas', logspace(-2, 4, 200), ...
    'nPerm', 10000, 'alpha', 0.05, 'correction', 'maxT', ...
    'backprojection', 'original', 'seed', 1);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts, f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
end
