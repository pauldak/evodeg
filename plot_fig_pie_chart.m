% PLOT_FIG_PIE_CHART  Visualises evolutionary slope and area-between-curve comparisons.
%
% Reads evolutionary_rates_summary.xlsx and produces:
%   - Box plots of per-protein slopes (global vs degron) for structural, sequence,
%     and functional conservation, with paired t-tests and Cohen's d effect sizes
%   - Box plot of area-between-curves for all three conservation types, with
%     one-way ANOVA and Tukey post-hoc pairwise comparisons
%
% Requires on disk: evolutionary_rates_summary.xlsx
%   Columns: Type (Structural/Sequence/Functional), Slopes_global_total,
%            Slopes_degron_total, Areas_Between_1, Areas_Between_2, Areas_Between_3
% Produces:         4 figure windows; Cohen's d and ANOVA results printed to console

close all; clear all;

% ================== SETTINGS ==================
is_title = 1;  % Toggle subplot titles ON/OFF
T = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'evolutionary_rates_summary.xlsx'));

%% Filter rows by type
idx_structural = strcmp(T.Type, 'Structural');
idx_sequence   = strcmp(T.Type, 'Sequence');
idx_func       = strcmp(T.Type, 'Functional');

%% Extract the variables
% Slopes (Global vs Degron)
data_struct_global = T.Slopes_global_total(idx_structural);
data_struct_degs   = T.Slopes_degron_total(idx_structural);

data_seq_global    = T.Slopes_global_total(idx_sequence);
data_seq_degs      = T.Slopes_degron_total(idx_sequence);

data_func_global   = T.Slopes_global_total(idx_func);
data_func_degs     = T.Slopes_degron_total(idx_func);

% Areas Between 1–3
area_structural = T.Areas_Between_1;
area_sequence   = T.Areas_Between_2;
area_func       = T.Areas_Between_3;

%% Labels
groupLabels = {'Global', 'Degron'};
changeLabels = {'Structural', 'Sequence', 'Functional'};

% ================== COLORS ==================
colors = [
    0.0000 0.4470 0.7410;  % blue
    0.8500 0.3250 0.0980;  % orange
];

% ================== PLOT 1: Structural ==================
figure();
plot_box_with_stats(data_struct_global, data_struct_degs, groupLabels, 'Structural change Between ortholog pairs', colors)
% ================== PLOT 2: Sequence ==================
figure();
plot_box_with_stats(data_seq_global, data_seq_degs, groupLabels, 'Sequence change between ortholog pairs', colors)
% ================== PLOT 3: Functional ==================
figure();
plot_box_with_stats(data_func_global, data_func_degs, groupLabels, 'Functional change between ortholog pairs', colors)

% ================== EFFECT SIZES (Cohen's d) ==================
d_struct = mean(data_struct_global - data_struct_degs) / std(data_struct_global - data_struct_degs);
d_seq    = mean(data_seq_global - data_seq_degs) / std(data_seq_global - data_seq_degs);
d_func   = mean(data_func_global - data_func_degs) / std(data_func_global - data_func_degs);
fprintf('Cohen''s d — Structural: %.3f\n', d_struct);
fprintf('Cohen''s d — Sequence:   %.3f\n', d_seq);
fprintf('Cohen''s d — Functional: %.3f\n', d_func);
% ================== PLOT 4: Area Between Graphs (ANOVA + Tukey) ==================
figure();
data = [area_structural, area_sequence, area_func];
boxplot(data, 'Labels', changeLabels);
hold on;
% Means
means = mean(data);
plot(1:3, means, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
ylabel('Area Between Graphs');
if is_title
    title('Area between graphs');
end
% ===== ANOVA + Tukey =====
[p, tbl, stats] = anova1(data, [], 'off');
c = multcompare(stats, 'Display', 'off');   % Tukey table
eta_sq = tbl{2,2} / tbl{4,2};
fprintf('Eta-squared (ANOVA): %.4f\n', eta_sq);
yMax = max(data(:));
% ---- Apply custom colors ----
areaColors = [
    0.6353 0.0784 0.1843;   % red
    0.4627 0.3294 0.8392;   % purple
    0.2039 0.5961 0.2706;   % dark green
];
boxes = findobj(gca, 'Tag', 'Box');
for i = 1:length(boxes)
    patch(get(boxes(i),'XData'), get(boxes(i),'YData'), ...
          areaColors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', areaColors(i,:));
end
% --- Draw Tukey pairwise bars ---
pairY = yMax * 1.1;      % height for bars
stepY = yMax * 0.08;     % spacing between bars
pairs = [1 2; 1 3; 2 3];
for i = 1:3
    p_val = c(i,6);          % column 6 = p-value
    s = sigstar(p_val);      % convert to stars
    x1 = pairs(i,1);
    x2 = pairs(i,2);
    % plot horizontal bar
    line([x1 x2], [pairY pairY], 'Color','k','LineWidth',1.2);
    % text label
    text(mean([x1 x2]), pairY + stepY/4, s, 'HorizontalAlignment','center', 'FontSize', 11);
    pairY = pairY + stepY;   % move next bar up
end
set(gca, 'FontSize', 11);
ylim([min(data(:))  yMax*1.55]);
grid on;
hold off;

function plot_box_with_stats(data1, data2, labels, titleStr, colors)
    % Combine data
    allData   = [data1; data2];
    groupsNum = [ones(size(data1)); 2*ones(size(data2))];
    % Boxplot
    boxplot(allData, groupsNum, 'Labels', labels);
    hold on;
    % ---- Apply custom colors ----
    boxes = findobj(gca, 'Tag', 'Box');
    for i = 1:length(boxes)
        patch(get(boxes(i),'XData'), get(boxes(i),'YData'), ...
              colors(i,:), 'FaceAlpha', 0.5, 'EdgeColor', colors(i,:));
    end
    % Means
    means = [mean(data1), mean(data2)];
    plot(1:2, means, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    ylabel('Rate of Change (AU)');
    title(titleStr);
    % Significance test
    [h, p] = ttest(data1 - data2);  % Paired t-test  % Two-sample t-test
    yMax = max(allData);
    line([1 2], [yMax*1.05 yMax*1.05], 'Color','k','LineWidth',1.2);

    if p < 0.001
        sig = '***';
        disp(p);
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end
    set(gca, 'FontSize', 11);
    text(1.5, yMax*1.08, sprintf('%s', sig), ...
         'HorizontalAlignment', 'center', 'FontSize', 11);

    grid on;
    hold off;
end

function stars = sigstar(p)
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    elseif p < 0.05
        stars = '*';
    else
        stars = 'ns';
    end
end