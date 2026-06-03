% DEMO  Reproduces the main statistical figures from the evodeg pipeline
% using a representative subset of the pre-computed coupling data.
%
% No external dependencies beyond standard MATLAB toolboxes.
% No PyMOL, no QCDPred, no AlphaFold downloads required.
%
% Produces 6 figures:
%   Figure 1  - RMSD distributions per organism (box plots, degron vs global)
%   Figure 2  - Mean RMSD per organism (paired scatter, degron vs global)
%   Figure 3  - Sequence alignment score distributions (box plots)
%   Figure 4  - Mean alignment score per organism (paired scatter)
%   Figure 5  - Successful vs failed structural coupling (stacked bar + chi-square)
%   Figure 6  - Evolutionary rate of change: structural, sequence, functional
%
% Required toolboxes: Statistics and Machine Learning Toolbox
%
% Usage:
%   Run this script from the repository root directory.
%   Demo data files must be present in the demo_data/ subfolder.

clearvars; close all; clc;

% =========================================================================
% 1. CONFIGURATION
% =========================================================================
demo_dir = fullfile(fileparts(mfilename('fullpath')), 'demo_data');

organism_list = {
    'Saccharomyces_cerevisiae', ...
    'Amphimedon_queenslandica', ...
    'Caenorhabditis_elegans', ...
    'Nematostella_vectensis', ...
    'Drosophila_melanogaster', ...
    'Danio_rerio', ...
    'Homo_sapiens'
};
organism_labels = strrep(organism_list, '_', ' ');

% Exclude S. cerevisiae (index 1) - it is the reference organism and its
% self-comparison RMSD is trivially zero.
comparison_idx    = 2:length(organism_list);
comparison_list   = organism_list(comparison_idx);
comparison_labels = organism_labels(comparison_idx);
n_orgs            = length(comparison_list);

fprintf('Loading demo data from: %s\n', demo_dir);

% =========================================================================
% 2. LOAD DATA
% =========================================================================
global_tbl = readtable(fullfile(demo_dir, 'global_coupling_table_degron_alignment.xlsx'));
deg_tbl    = readtable(fullfile(demo_dir, 'global_deg_coupling_table_degron_alignment.xlsx'));
reg_global = readtable(fullfile(demo_dir, 'global_coupling_table_regular_alignment.xlsx'));
reg_deg    = readtable(fullfile(demo_dir, 'global_deg_coupling_table_regular_alignment.xlsx'));
evo_tbl    = readtable(fullfile(demo_dir, 'evolutionary_rates_summary.xlsx'));

fprintf('Loaded %d global rows, %d degron rows.\n', height(global_tbl), height(deg_tbl));

% =========================================================================
% 3. EXTRACT PER-ORGANISM VECTORS (Figures 1-4)
%    Separate failed couplings (RMSD == 99) from valid ones.
%    Perform Wilcoxon rank-sum test (unpaired, per organism) comparing
%    degron vs global interval RMSD and alignment scores.
% =========================================================================
rmsd_global_vec   = {};
rmsd_deg_vec      = {};
align_global_vec  = {};
align_deg_vec     = {};
rmsd_global_mean  = zeros(1, n_orgs);
rmsd_deg_mean     = zeros(1, n_orgs);
align_global_mean = zeros(1, n_orgs);
align_deg_mean    = zeros(1, n_orgs);
rmsd_sig          = zeros(1, n_orgs);
align_sig         = zeros(1, n_orgs);

fprintf('\n%-30s  %6s  %6s  %8s  |  %6s  %6s  %8s\n', ...
    'Organism', 'N_gl', 'N_dg', 'p_RMSD', 'N_gl', 'N_dg', 'p_align');

for i = 1:n_orgs
    org = comparison_list{i};

    % Global rows for this organism (valid only)
    g_rows = global_tbl(strcmp(global_tbl.organism_name, org) & ...
                        global_tbl.coupling_rmsd ~= 99, :);
    % Degron rows for this organism (valid only)
    d_rows = deg_tbl(strcmp(deg_tbl.organism_name, org) & ...
                     deg_tbl.coupling_rmsd ~= 99, :);

    rmsd_global_vec{i}  = g_rows.coupling_rmsd;
    rmsd_deg_vec{i}     = d_rows.coupling_rmsd;
    align_global_vec{i} = g_rows.algmnt_score;
    align_deg_vec{i}    = d_rows.algmnt_score;

    rmsd_global_mean(i)  = mean(g_rows.coupling_rmsd);
    rmsd_deg_mean(i)     = mean(d_rows.coupling_rmsd);
    align_global_mean(i) = mean(g_rows.algmnt_score);
    align_deg_mean(i)    = mean(d_rows.algmnt_score);

    % Wilcoxon rank-sum tests (unpaired)
    [p_rmsd,  ~] = ranksum(g_rows.coupling_rmsd, d_rows.coupling_rmsd);
    [p_align, ~] = ranksum(g_rows.algmnt_score,  d_rows.algmnt_score);

    rmsd_sig(i)  = sig_level(p_rmsd);
    align_sig(i) = sig_level(p_align);

    fprintf('%-30s  %6d  %6d  %8.4g  |  %6d  %6d  %8.4g\n', ...
        strrep(org,'_',' '), height(g_rows), height(d_rows), ...
        p_rmsd, height(g_rows), height(d_rows), p_align);
end

% Wilcoxon signed-rank on paired means (across organisms)
[p_rmsd_paired,  ~] = signrank(rmsd_deg_mean,  rmsd_global_mean);
[p_align_paired, ~] = signrank(align_deg_mean, align_global_mean);
fprintf('\nPaired signed-rank (means across organisms):\n');
fprintf('  RMSD:      p = %.4g\n', p_rmsd_paired);
fprintf('  Alignment: p = %.4g\n', p_align_paired);

% =========================================================================
% 4. FIGURE 1 - RMSD box plots per organism
% =========================================================================
figure('Name', 'Figure 1: RMSD distributions', 'Color', 'w');
plot_grouped_boxes(rmsd_deg_vec, rmsd_global_vec, ...
    'Structural RMSD — Degron vs Global intervals', ...
    'RMSD (Å)', comparison_labels, rmsd_sig);

% =========================================================================
% 5. FIGURE 2 - Mean RMSD paired scatter
% =========================================================================
figure('Name', 'Figure 2: Mean RMSD per organism', 'Color', 'w');
plot_paired_scatter(rmsd_deg_mean, rmsd_global_mean, ...
    'Mean RMSD (Å)', p_rmsd_paired);

% =========================================================================
% 6. FIGURE 3 - Alignment score box plots
% =========================================================================
figure('Name', 'Figure 3: Alignment score distributions', 'Color', 'w');
plot_grouped_boxes(align_deg_vec, align_global_vec, ...
    'Sequence Alignment Score — Degron vs Global intervals', ...
    'Alignment Score (AU)', comparison_labels, align_sig);

% =========================================================================
% 7. FIGURE 4 - Mean alignment score paired scatter
% =========================================================================
figure('Name', 'Figure 4: Mean alignment score per organism', 'Color', 'w');
plot_paired_scatter(align_deg_mean, align_global_mean, ...
    'Sequence Alignment Score (AU)', p_align_paired);

% =========================================================================
% 8. FIGURE 5 - Successful vs failed coupling (stacked bar + chi-square)
% =========================================================================
figure('Name', 'Figure 5: Coupling success rate', 'Color', 'w');

% Classify rows: good match = valid RMSD and global_rmsd / global_rmsd_projections < 2
classify = @(t) (t.coupling_rmsd ~= 99) & ...
                (t.global_rmsd ./ t.global_rmsd_projections < 2);

deg_good  = sum(classify(reg_deg));
deg_bad   = height(reg_deg) - deg_good;
glob_good = sum(classify(reg_global));
glob_bad  = height(reg_global) - glob_good;

deg_pct  = 100 * [deg_good,  deg_bad]  / height(reg_deg);
glob_pct = 100 * [glob_good, glob_bad] / height(reg_global);

data = [deg_pct; glob_pct];
b = bar(data, 'stacked');
b(1).FaceColor = [0.25 0.55 0.25];
b(2).FaceColor = [0.75 0.25 0.25];
set(gca, 'XTickLabel', {'Degron Regions', 'Global Regions'}, 'FontSize', 11);
ylabel('Percentage of Intervals (%)');
title('Successful vs Failed Structural Coupling');
ylim([0 120]);
legend({'Successful Match', 'Failed Match'}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');
grid on;

% Chi-square test of independence
contingency = [deg_good, deg_bad; glob_good, glob_bad];
[~, chi2, p_chi2] = crosstab(contingency(:,1), contingency(:,2));
stars = sig_stars(p_chi2);
ymax = max(sum(data, 2)) + 5;
line([1 2], [ymax ymax], 'Color', 'k', 'LineWidth', 1.5);
text(1.5, ymax + 2, stars, 'HorizontalAlignment', 'center', ...
    'FontSize', 16, 'FontWeight', 'bold');
fprintf('\nFigure 5 — Chi-square p = %.4g  (%s)\n', p_chi2, stars);

% =========================================================================
% 9. FIGURE 6 - Evolutionary rate of change (structural / sequence / functional)
% =========================================================================
figure('Name', 'Figure 6: Evolutionary rate of change', 'Color', 'w');

types  = {'Structural', 'Sequence', 'Functional'};
colors = [0.0000 0.4470 0.7410;   % blue  - global
          0.8500 0.3250 0.0980];  % orange - degron

fprintf('\nFigure 6 — Paired t-tests (global vs degron slopes):\n');
for k = 1:3
    subplot(1, 3, k);
    idx   = strcmp(evo_tbl.Type, types{k});
    glob  = evo_tbl.Slopes_global_total(idx);
    deg   = evo_tbl.Slopes_degron_total(idx);
    allD  = [glob; deg];
    grp   = [ones(size(glob)); 2*ones(size(deg))];
    boxplot(allD, grp, 'Labels', {'Global', 'Degron'}, 'Symbol', '');
    hold on;
    boxes = findobj(gca, 'Tag', 'Box');
    for b = 1:length(boxes)
        patch(get(boxes(b),'XData'), get(boxes(b),'YData'), ...
              colors(b,:), 'FaceAlpha', 0.5, 'EdgeColor', colors(b,:));
    end
    plot(1:2, [mean(glob), mean(deg)], 'kd', ...
        'MarkerFaceColor', 'k', 'MarkerSize', 6);
    [~, p, ~, st] = ttest(glob - deg);
    d = mean(glob - deg) / std(glob - deg);
    ymax = max(allD);
    line([1 2], [ymax*1.05 ymax*1.05], 'Color', 'k', 'LineWidth', 1.2);
    text(1.5, ymax*1.08, sig_stars(p), 'HorizontalAlignment', 'center', 'FontSize', 11);
    ylabel('Rate of Change (AU)');
    title(types{k});
    set(gca, 'FontSize', 10);
    grid on; hold off;
    fprintf('  %s: t(%d) = %.2f, p = %.4g, d = %.2f\n', ...
        types{k}, st.df, st.tstat, p, d);
end
sgtitle('Evolutionary Rate of Change: Global vs Degron Regions', 'FontSize', 13);

fprintf('\nDemo complete. 6 figures produced.\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function plot_grouped_boxes(x_data, y_data, ttl, ylbl, org_names, sig_arr)
% Grouped box plots: x_data = degron (blue), y_data = global (orange).
    n = length(org_names);
    colors = [0.0000 0.4470 0.7410;   % blue  - degron
              0.8500 0.3250 0.0980];  % orange - global

    pos_deg  = (1:n) - 0.18;
    pos_glob = (1:n) + 0.18;

    all_deg  = vertcat(x_data{:});
    all_glob = vertcat(y_data{:});
    lbl_deg  = repelem(org_names, cellfun(@numel, x_data))';
    lbl_glob = repelem(org_names, cellfun(@numel, y_data))';

    hold on;
    boxplot(all_deg,  lbl_deg,  'Positions', pos_deg,  'Widths', 0.3, ...
        'Symbol', '', 'Colors', colors(1,:));
    boxplot(all_glob, lbl_glob, 'Positions', pos_glob, 'Widths', 0.3, ...
        'Symbol', '', 'Colors', colors(2,:));

    for i = 1:n
        plot(pos_deg(i),  mean(x_data{i}), 'k+', 'MarkerSize', 6, 'LineWidth', 1);
        plot(pos_glob(i), mean(y_data{i}), 'k+', 'MarkerSize', 6, 'LineWidth', 1);
    end

    all_vals = [all_deg; all_glob];
    Q3  = prctile(all_vals, 90);
    ymax = Q3 + 1.5 * iqr(all_vals);
    ymax = max(all_vals(all_vals <= ymax));

    for i = 1:n
        if sig_arr(i) == 0; continue; end
        line([pos_deg(i) pos_glob(i)], [ymax*1.08 ymax*1.08], ...
            'Color', 'k', 'LineWidth', 1);
        text(i, ymax*1.11, repmat('*',1,sig_arr(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
    end

    set(gca, 'XTick', 1:n, 'XTickLabel', org_names, ...
        'XTickLabelRotation', 20, 'FontSize', 9);
    ylabel(ylbl);
    title(ttl);
    ylim([-0.5 ymax*1.2]);
    grid on; hold off;

    % Manual legend
    hd = plot(nan, nan, 's', 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));
    hg = plot(nan, nan, 's', 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', colors(2,:));
    legend([hd hg], {'Degron intervals', 'Global intervals'}, 'Location', 'northeast');
end


function plot_paired_scatter(x_data, y_data, ylbl, p_val)
% Paired scatter: each dot = one organism mean; left = degron, right = global.
    colors = [0.0000 0.4470 0.7410;
              0.8500 0.3250 0.0980];
    n = length(x_data);
    hold on;
    jitter = 0.04;
    scatter(ones(n,1) + (rand(n,1)-0.5)*jitter, x_data, 60, ...
        'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', 'k');
    scatter(2*ones(n,1) + (rand(n,1)-0.5)*jitter, y_data, 60, ...
        'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'k');
    for i = 1:n
        plot([1 2], [x_data(i) y_data(i)], 'k-', 'LineWidth', 0.5, 'Color', [0.6 0.6 0.6]);
    end
    errorbar(1, mean(x_data), std(x_data)/sqrt(n), 'k', 'LineWidth', 2);
    errorbar(2, mean(y_data), std(y_data)/sqrt(n), 'k', 'LineWidth', 2);
    plot(1, mean(x_data), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot(2, mean(y_data), 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 8);
    plot([1 2], [mean(x_data) mean(y_data)], 'k--', 'LineWidth', 2);
    yl = ylim;
    ypos = yl(2) * 1.05;
    line([1 2], [ypos ypos], 'Color', 'k', 'LineWidth', 1.2);
    text(1.5, ypos*1.03, sig_stars(p_val), 'HorizontalAlignment', 'center', 'FontSize', 11);
    ylim([yl(1) ypos*1.12]);
    xlim([0.5 2.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Degron regions', 'Global regions'}, 'FontSize', 10);
    ylabel(ylbl);
    grid on; hold off;
end


function level = sig_level(p)
% Convert p-value to integer significance level (0-3).
    if     p < 0.001; level = 3;
    elseif p < 0.01;  level = 2;
    elseif p < 0.05;  level = 1;
    else;             level = 0;
    end
end


function s = sig_stars(p)
% Convert p-value to significance star string.
    if     p < 0.001; s = '***';
    elseif p < 0.01;  s = '**';
    elseif p < 0.05;  s = '*';
    else;             s = 'ns';
    end
end
