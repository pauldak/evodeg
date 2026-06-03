% STATS_EXTRACT  Computes paired t-tests and effect sizes for evolutionary slope comparisons.
%
% Reads the deep evolutionary analysis output table and computes:
%   - Paired t-tests (global vs degron slopes) for structural, sequence, and functional
%     data types (Figures 5C, 5D, 5E), with Cohen's d effect size
%   - One-way ANOVA with Tukey post-hoc comparisons for area-between-curve metrics (Figure 5F)
%
% Requires on disk: evolutionary_rates_summary.xlsx
%   Columns: Type (Structural/Sequence/Functional), Slopes_global_total,
%            Slopes_degron_total, Areas_Between_1, Areas_Between_2, Areas_Between_3
% Produces: statistics printed to console

T = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'evolutionary_rates_summary.xlsx'));

% Figure 5C - Structural
idx = strcmp(T.Type, 'Structural');
[~, p, ~, stats] = ttest(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
d = mean(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx)) / std(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
fprintf('Structural: t(%d) = %.2f, p = %g, d = %.2f\n', stats.df, stats.tstat, p, d);

% Figure 5D - Sequence
idx = strcmp(T.Type, 'Sequence');
[~, p, ~, stats] = ttest(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
d = mean(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx)) / std(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
fprintf('Sequence: t(%d) = %.2f, p = %g, d = %.2f\n', stats.df, stats.tstat, p, d);

% Figure 5E - Functional
idx = strcmp(T.Type, 'Functional');
[~, p, ~, stats] = ttest(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
d = mean(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx)) / std(T.Slopes_global_total(idx) - T.Slopes_degron_total(idx));
fprintf('Functional: t(%d) = %.2f, p = %g, d = %.2f\n', stats.df, stats.tstat, p, d);

%%
idx = strcmp(T.Type, 'Structural');  % pick any one type
areas = [T.Areas_Between_1(idx), T.Areas_Between_2(idx), T.Areas_Between_3(idx)];
group = [ones(720,1); 2*ones(720,1); 3*ones(720,1)];
[p, tbl, stats] = anova1(areas(:), group, 'off');
eta_sq = tbl{2,2} / tbl{4,2};
fprintf('F(%d, %d) = %.2f, p = %g, eta2 = %.4f\n', tbl{2,3}, tbl{3,3}, tbl{2,5}, tbl{2,6}, eta_sq);
[c,m] = multcompare(stats, 'Display', 'off');
for i = 1:size(c,1)
    fprintf('Group %d vs %d: p = %g\n', c(i,1), c(i,2), c(i,6));
end