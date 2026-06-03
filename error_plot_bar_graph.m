% ERROR_PLOT_BAR_GRAPH  Stacked bar plot of successful vs failed interval couplings.
%
% Classifies each entry in the global coupling tables as a successful match
% (coupling_rmsd ~= 99 AND global_rmsd / global_rmsd_projections < 2) or a failed match,
% and plots the proportions as stacked bars for degron vs global regions.
% Performs a Chi-square test of independence on the 2x2 contingency table.
%
% Requires on disk: global_deg_coupling_table_regular_alignment.xlsx
%                   global_coupling_table_regular_alignment.xlsx
% Produces:         figure window with significance annotation
close all;
clear all;
degron_table = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'global_deg_coupling_table_regular_alignment.xlsx'));
global_table = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'global_coupling_table_regular_alignment.xlsx'));

% Define a function to classify matches
function match_labels = classify_matches(tbl)
    match_labels = strings(height(tbl),1);
    for i = 1:height(tbl)
        if tbl.coupling_rmsd(i) == 99
            match_labels(i) = "bad_match";
        elseif (tbl.global_rmsd(i) / tbl.global_rmsd_projections(i)) >= 2
            match_labels(i) = "bad_match";
        else
            match_labels(i) = "good_match";
        end
    end
end

% Apply classification to each population
degron_labels = classify_matches(degron_table);
global_labels = classify_matches(global_table);

% Count matches as percentages
total_degrons = length(degron_labels);
total_global = length(global_labels);
% Reverse column order so good_match is first (bottom of stack)
degron_counts = (100 / total_degrons) * [sum(degron_labels == "good_match"), sum(degron_labels == "bad_match")];
global_counts = (100 / total_global) * [sum(global_labels == "good_match"), sum(global_labels == "bad_match")];

% Combine
data = [degron_counts; global_counts];

% Plot as stacked bars
figure;
bar(data, 'stacked');
set(gca, 'XTickLabel', {'Degron Regions', 'Global Regions'});
ylabel('Percentage of Matches (%)');
title('Degron Vs Global Regions, Successful and Unsuccessful Matches');
ylim([0 120]);  % Cap at 100% with a bit of headroom
grid on;

% Chi-square test of independence using crosstab
group_labels = [repmat("Degron", height(degron_table), 1); repmat("Global", height(global_table), 1)];
match_labels = [degron_labels; global_labels];

[~, chi2, p_value, labels] = crosstab(group_labels, match_labels);

% Display p-value
disp(['Chi-square p-value: ', num2str(p_value, '%.4g')]);
% Optional: Annotate p-value on plot
if p_value < 0.001
    stars = '***';
elseif p_value < 0.01
    stars = '**';
elseif p_value < 0.05
    stars = '*';
else
    stars = 'n.s.';
end

y_max = max(sum(data, 2)) + 5;
bar_positions = [1, 2];

hold on;
line([bar_positions(1), bar_positions(2)], [y_max, y_max], 'Color', 'k', 'LineWidth', 1.5);
line([bar_positions(1), bar_positions(1)], [y_max-1, y_max], 'Color', 'k', 'LineWidth', 1.2);
line([bar_positions(2), bar_positions(2)], [y_max-1, y_max], 'Color', 'k', 'LineWidth', 1.2);

text(mean(bar_positions), y_max + 1.5, stars, ...
    'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
legend({'Successful Match', 'Unsuccessful Match'}, 'Location', 'northoutside', 'Orientation', 'horizontal');

hold off;
