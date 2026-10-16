function plot_dice_results(res, yLimits, starsFrom)
% PLOT_DICE_RESULTS  Boxplots of the null Dice distributions with the observed Dice
% of each syndrome (red dots) and significance stars (as in Figure 5).
%
%   plot_dice_results(res, [-0.02 0.4])              stars from uncorrected p
%   plot_dice_results(res, [-0.02 0.4], 'holm')      stars from Holm-corrected p
%
% Stars: * p < 0.05, ** p < 0.01, *** p < 0.001 (one-sided permutation tests).

if nargin < 2 || isempty(yLimits), yLimits = [-0.02, 1.05 * max(res.dice_real)]; end
if nargin < 3, starsFrom = 'uncorrected'; end
S = numel(res.syndromeNames);
if strcmp(starsFrom, 'holm'), p = res.pHolm; else, p = res.pUpper; end

figure('Color', 'w'); hold on;
boxplot(res.dice_nulls, 'Labels', res.syndromeNames, 'Colors', 'k', 'Symbol', '');
set(gca, 'XTickLabelRotation', 30, 'FontSize', 15);
ylabel('Dice Coefficient');
ylim(yLimits);
title(['Dice Coefficients: Aphasic Syndrome vs PC1 - ' res.cfg.mode]);

colors = lines(S);
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:numel(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), colors(numel(boxes) - j + 1, :), ...
          'FaceAlpha', 0.4, 'EdgeColor', 'none');
end
for i = 1:S
    plot(i, res.dice_real(i), 'r.', 'MarkerSize', 30);
end
for i = 1:S
    if p(i) < 0.001, sig = '***'; elseif p(i) < 0.01, sig = '**'; elseif p(i) < 0.05, sig = '*'; else, sig = ''; end
    if ~isempty(sig)
        y = max(res.dice_nulls(:, i)) + 0.04 * diff(yLimits);
        text(i, y, sig, 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'k');
    end
end
end
