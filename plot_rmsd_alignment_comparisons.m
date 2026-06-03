% BAR_PLOTS_NO_SC  Statistical analysis and visualization of RMSD and alignment scores,
% excluding S. cerevisiae from all comparisons.
%
% Identical to bar_plots.m but drops index 1 (S. cerevisiae) from all plots and tests.
% S. cerevisiae is the structural reference and always has RMSD = 0 against itself;
% excluding it avoids distortion of cross-species comparisons.
%
% Requires loaded: global_coupling_table, global_deg_coupling_table (set in addresses.m),
%                  organism_list, organism_list_txt, number_of_organisms (from folder_paths.mat)
% Produces:        statistic_test.txt, figure windows
clearvars -except hogs_array; clc; close all;
load('folder_paths.mat');
cd(matlab_folder);

rmsd_mean_vector = [];
align_mean_vector = [];
prcnt_bad_projection = [];
prcnt_99 = [];

rmsd_deg_mean_vector = [];
align_deg_mean_vector = [];
prcnt_deg_bad_projection = [];
prcnt_deg_99 = [];
%all merged
all_rmsd_deg_vector = [];
all_align_deg_vector = [];
all_rmsd_global_vector = [];
all_align_global_vector = [];
%non means
rmsd_deg_vector = {};
rmsd_global_vector = {};
align_deg_vector = {};
align_global_vector = {};



rmsd_significants = [];
align_significant = [];
results_filename = 'statistic_test.txt';
fid = fopen(results_filename, 'w');
% saving data
p_rmsd_per_org = zeros(number_of_organisms, 1);
p_align_per_org = zeros(number_of_organisms, 1);
%
for i = 1: number_of_organisms
    crnt_organism = organism_list{i};
    org_indx = contains(global_coupling_table.organism_name, crnt_organism);
    org_sub_tble = global_coupling_table(org_indx, :);
    
    idx_99 = find(org_sub_tble.coupling_rmsd == 99);
    idx_not_99 = find(org_sub_tble.coupling_rmsd ~= 99);
    org_sub_tble_only_99 = org_sub_tble(idx_99, :);
    org_sub_tble_not_99 = org_sub_tble(idx_not_99, :);
    %rmsd scroe
    crnt_average_rmsd = mean(org_sub_tble_not_99.coupling_rmsd);
    rmsd_mean_vector = [rmsd_mean_vector, crnt_average_rmsd];
    rmsd_global_vector = [rmsd_global_vector, org_sub_tble_not_99.coupling_rmsd];
    all_rmsd_global_vector = [all_rmsd_global_vector; org_sub_tble_not_99.coupling_rmsd];
    %align data
    crnt_average_align = mean(org_sub_tble_not_99.algmnt_score);
    align_mean_vector = [align_mean_vector, crnt_average_align];
    align_global_vector = [align_global_vector, org_sub_tble_not_99.algmnt_score];
    all_align_global_vector = [all_align_global_vector; org_sub_tble_not_99.algmnt_score];
    % false projection
    global_rmsd_projections = org_sub_tble_not_99.global_rmsd_projections;
    global_rmsd = org_sub_tble_not_99.global_rmsd;
    diff_rmsd_proj_global = global_rmsd_projections ./ global_rmsd;
    number_of_false_align = sum(diff_rmsd_proj_global > 2);
    prcnt_crnt_bad_projection = (100 / size(org_sub_tble_not_99,1)) * number_of_false_align;
    prcnt_bad_projection = [prcnt_bad_projection, prcnt_crnt_bad_projection];
    % 99 cases
    crnt_prcnt_99 = (100 / size(org_sub_tble, 1)) * size(org_sub_tble_only_99, 1);
    prcnt_99 = [prcnt_99 ,crnt_prcnt_99];
    
    
    %% degron version
    org_deg_indx = contains(global_deg_coupling_table.organism_name, crnt_organism);
    org_deg_sub_tble = global_deg_coupling_table(org_deg_indx, :);
    
    idx_deg_99 = find(org_deg_sub_tble.coupling_rmsd == 99);
    idx_deg_not_99 = find(org_deg_sub_tble.coupling_rmsd ~= 99);
    org_deg_sub_tble_only_99 = org_deg_sub_tble(idx_deg_99, :);
    org_deg_sub_tble_not_99 = org_deg_sub_tble(idx_deg_not_99, :);
    %rmsd deg scroe
    crnt_average_deg_rmsd = mean(org_deg_sub_tble_not_99.coupling_rmsd);
    rmsd_deg_mean_vector = [rmsd_deg_mean_vector, crnt_average_deg_rmsd];
    rmsd_deg_vector = [rmsd_deg_vector, org_deg_sub_tble_not_99.coupling_rmsd];
    all_rmsd_deg_vector = [all_rmsd_deg_vector; org_deg_sub_tble_not_99.coupling_rmsd];
    %align data
    crnt_deg_average_align = mean(org_deg_sub_tble_not_99.algmnt_score);
    align_deg_mean_vector = [align_deg_mean_vector, crnt_deg_average_align];
    align_deg_vector = [align_deg_vector, org_deg_sub_tble_not_99.algmnt_score];
    all_align_deg_vector = [all_align_deg_vector; org_deg_sub_tble_not_99.algmnt_score];
    % false projection
    global_deg_rmsd_projections = org_deg_sub_tble_not_99.global_rmsd_projections;
    global_deg_rmsd = org_deg_sub_tble_not_99.global_rmsd;
    diff_deg_rmsd_proj_global = global_deg_rmsd_projections ./ global_deg_rmsd;
    number_deg_of_false_align = sum(diff_deg_rmsd_proj_global > 2);
    prcnt_deg_crnt_bad_projection = (100 / size(org_deg_sub_tble_not_99,1)) * number_deg_of_false_align;
    prcnt_deg_bad_projection = [prcnt_deg_bad_projection, prcnt_deg_crnt_bad_projection];
    % 99 cases
    crnt_deg_prcnt_99 = (100 / size(org_deg_sub_tble, 1)) * size(org_deg_sub_tble_only_99, 1);
    prcnt_deg_99 = [prcnt_deg_99 ,crnt_deg_prcnt_99];

    %% plotting population comparisons mann whitney test
    fprintf(fid, '%s\n', crnt_organism);
    a = org_sub_tble_not_99.coupling_rmsd;
    b = org_deg_sub_tble_not_99.coupling_rmsd;
    [p1, h1] = ranksum(a, b);
    p_rmsd_per_org(i) = p1;
    if p1 < 0.001
        h1 = 3;
    elseif p1 < 0.01
        h1 = 2;
    elseif p1 < 0.05
        h1 = 1;
    else
        h1 = 0;
    end
    rmsd_significants = [rmsd_significants, h1];
    fprintf(fid, 'coupling RMSD p value - %s, hypothesys rejection - %s\n', num2str(p1), num2str(h1));
    %alignment significance test
    c = org_sub_tble_not_99.algmnt_score;
    d = org_deg_sub_tble_not_99.algmnt_score;
    [p2, h2] = ranksum(c, d);
    p_align_per_org(i) = p2;
    fprintf(fid, 'Alignment p value - %s, hypothesys rejection - %s\n', num2str(p2), num2str(h2));
    align_significant = [align_significant, h2];
    
    %%
end
%all average RMSD significance test
all_non_99_idx = find(global_coupling_table.coupling_rmsd ~= 99);
all_non_99_table = global_coupling_table(all_non_99_idx, :);
all_deg_non_99_idx = find(global_deg_coupling_table.coupling_rmsd ~= 99);
all_deg_non_99_table = global_deg_coupling_table(all_deg_non_99_idx, :);
all_rmsd = all_non_99_table.coupling_rmsd;
all_deg_rmsd = all_deg_non_99_table.coupling_rmsd;
[p3, h3] = ranksum(all_rmsd, all_deg_rmsd);
fprintf(fid, 'RMSD all p value - %s, hypothesys rejection - %s\n', num2str(p3), num2str(h3));
%all average alignment significance test
all_alignment = all_non_99_table.algmnt_score;
all_deg_alignment = all_deg_non_99_table.algmnt_score;
[p4, h4] = ranksum(all_alignment, all_deg_alignment);
fprintf(fid, 'Alignment all p value - %s, hypothesys rejection - %s\n', num2str(p4), num2str(h4));

fclose(fid);
% ================== STATS FOR FIGURE LEGENDS ==================

% --- 3A: RMSD per organism (unpaired) ---
fprintf('\n========== Figure 3A (RMSD — per organism) ==========\n');
fprintf('Test: Wilcoxon rank-sum (unpaired)\n');
fprintf('Overall: N_global = %d, N_degron = %d, p = %e\n', ...
    length(all_rmsd_global_vector), length(all_rmsd_deg_vector), p3);
for i = 1:number_of_organisms
    fprintf('%s: N_global = %d, N_degron = %d, p = %e\n', ...
        organism_list{i}, length(rmsd_global_vector{i}), length(rmsd_deg_vector{i}), p_rmsd_per_org(i));
end

% --- 3B: RMSD means across species (paired) ---
fprintf('\n========== Figure 3B (RMSD — paired means) ==========\n');
fprintf('Test: Wilcoxon signed-rank (paired)\n');
[p_3B, ~] = signrank(rmsd_deg_mean_vector, rmsd_mean_vector);
fprintf('N = %d organism pairs, p = %e\n', length(rmsd_mean_vector), p_3B);
fprintf('Mean degron = %.4f, Mean global = %.4f\n', ...
    mean(rmsd_deg_mean_vector(2:end)), mean(rmsd_mean_vector(2:end)));

% --- 3C: Sequence alignment per organism (unpaired) ---
fprintf('\n========== Figure 3C (Sequence Alignment — per organism) ==========\n');
fprintf('Test: Wilcoxon rank-sum (unpaired)\n');
fprintf('Overall: N_global = %d, N_degron = %d, p = %e\n', ...
    length(all_align_global_vector), length(all_align_deg_vector), p4);
for i = 1:number_of_organisms
    fprintf('%s: N_global = %d, N_degron = %d, p = %e\n', ...
        organism_list{i}, length(align_global_vector{i}), length(align_deg_vector{i}), p_align_per_org(i));
end

% --- 3D: Sequence alignment means across species (paired) ---
fprintf('\n========== Figure 3D (Sequence Alignment — paired means) ==========\n');
fprintf('Test: Wilcoxon signed-rank (paired)\n');
[p_3D, ~] = signrank(align_deg_mean_vector, align_mean_vector);
fprintf('N = %d organism pairs, p = %e\n', length(align_mean_vector), p_3D);
fprintf('Mean degron = %.4f, Mean global = %.4f\n', ...
    mean(align_deg_mean_vector(2:end)), mean(align_mean_vector(2:end)));
%%
% exporting table
my_table = table(organism_list', rmsd_mean_vector', rmsd_deg_mean_vector', align_mean_vector', align_deg_mean_vector', ...
    prcnt_bad_projection', prcnt_deg_bad_projection', prcnt_99', prcnt_deg_99');
my_table.Properties.VariableNames = ...
    {'organism_list', 'rmsd_mean_vector', 'rmsd_deg_mean_vector', 'align_mean_vector', 'align_deg_mean_vector', ...
                        'prcnt_bad_projection', 'prcnt_deg_bad_projection', 'prcnt_99', 'prcnt_deg_99'};
%writetable(my_table, 'results.xlsx');
    

%multiple bars - RMSD (excluding SC at index 1)
figure();
crnt_ylabel = 'Å';
plot_title = 'RMSD of all pairwise comparisons';
multi_box_plots(rmsd_deg_vector(2:end), rmsd_global_vector(2:end), plot_title, crnt_ylabel, rmsd_significants(2:end), organism_list_txt(2:end));

%single bars - RMSD
figure();
crnt_ylabel = 'RMSD (Å)';
plot_title = 'Mean RMSD of all pairwise comparisons';
plot_scatter_org(rmsd_deg_mean_vector(2:end), rmsd_mean_vector(2:end), plot_title, crnt_ylabel, h3);

%multiple bars - Alignment (excluding SC at index 1)
figure();
crnt_ylabel = 'Alignment Score (AU)';
plot_title = 'Degron Alignment score of all pairwise comparisons';
align_significant(1) = 0;  % S. cerevisiae self-alignment is excluded from significance testing.
multi_box_plots(align_deg_vector(2:end), align_global_vector(2:end), plot_title, crnt_ylabel, align_significant(2:end), organism_list_txt(2:end));

%single bars - Alignment
figure();
crnt_ylabel = 'Sequence Alignment (AU)';
plot_title = 'Mean degron functionality based alignment';
plot_scatter_org(align_deg_mean_vector(2:end), align_mean_vector(2:end), plot_title, crnt_ylabel, h4);

%single bars - Wrong Projection
figure();
crnt_ylabel = 'Percent (%)';
plot_title = 'Prcnt of wrong projection - RMSD based';
plot_bar_1(prcnt_deg_bad_projection, prcnt_bad_projection, plot_title, crnt_ylabel, 0);

%single bars - failed coupling
figure();
crnt_ylabel = 'Percent (%)';
plot_title = 'Prcnt of failed coupling';
plot_bar_1(prcnt_deg_99, prcnt_99, plot_title, crnt_ylabel, 0);
a=5;




function plot_bar_1(x_data, y_data, plot_title, crnt_ylabel, is_significant)
    % Create a bar plot
    all_organism_deg_rmsd_average = mean(x_data(2:end));
    all_organism_rmsd_average = mean(y_data(2:end));
    bar_data = [all_organism_deg_rmsd_average, all_organism_rmsd_average];
    names = categorical({'Degron regions', 'Entire Protein'});
    names = reordercats(names, {'Degron regions', 'Entire Protein'});
    % Create the bar plot with specified bar width (0.4)
    bar(names, bar_data, 0.4);
    title(plot_title);
    ylabel(crnt_ylabel);
    if is_significant
        % Calculate the x position for the center of the line
        x_position = [1, 2];
        x_center = mean(x_position);
        % Calculate the y position for the line and asterisk
        y_position = max(bar_data) * 1.1; % Adjust the vertical position as needed
        % Add a line connecting the two bars
        line(x_position, [y_position, y_position], 'Color', 'k', 'LineWidth', 1);
        % Add an asterisk above the center of the line
        text(x_center, y_position * 1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 16);
        ylim([0 (y_position * 1.2)]);
    else
        ylim([0 (max(bar_data) * 1.2)]);
    end
end


function multi_box_plots(x_data, y_data, plot_title, crnt_ylabel, significance_arr, org_names)
    % x_data, y_data: cell arrays of per-organism vectors
    % org_names: cell array of organism label strings (same length as x_data)
    n_orgs = length(org_names);

    % Build grouped data for boxplot
    data_set_degs = repelem(org_names, cellfun(@numel, x_data))';
    data_set_global = repelem(org_names, cellfun(@numel, y_data))';

    all_deg_data = vertcat(x_data{:});
    all_global_data = vertcat(y_data{:});
    allData = [all_deg_data; all_global_data];

    positions_deg = (1:n_orgs) - 0.18;
    positions_glob = (1:n_orgs) + 0.18;

    colors = [
        0.0000 0.4470 0.7410;  % blue - Degron
        0.8500 0.3250 0.0980;  % orange - Global
    ];

    hold on;
    % Plot degron boxes
    boxplot(all_deg_data, data_set_degs, 'Positions', positions_deg, ...
        'Widths', 0.3, 'Symbol', '', 'Colors', colors(1,:));


    % Plot global boxes
    boxplot(all_global_data, data_set_global, 'Positions', positions_glob, ...
        'Widths', 0.3, 'Symbol', '', 'Colors', colors(2,:));


    % Plot mean as + markers for each organism
    for i = 1:n_orgs
        mean_deg = mean(x_data{i});
        mean_glob = mean(y_data{i});
        plot(positions_deg(i), mean_deg, 'k+', 'MarkerSize', 6, 'LineWidth', 1);
        plot(positions_glob(i), mean_glob, 'k+', 'MarkerSize', 6, 'LineWidth', 1);
    end

    set(gca, 'XTick', 1:n_orgs, 'XTickLabel', org_names, 'FontSize', 10);
    ylabel(crnt_ylabel);
    title(plot_title);

    % IQR-based y-axis scaling (same as plot_box_1)
    IQR_val = iqr(allData);
    Q3 = prctile(allData, 90);
    upper_limit = Q3 + 1.5 * IQR_val;
    yMax = max(allData(allData <= upper_limit));

    % Find max whisker height for significance marker placement
    whisker_lines = findobj(gca, 'Tag', 'Upper Whisker');
    max_whisker = max(cellfun(@max, get(whisker_lines, 'YData')));

    %% drawing significance markers
    for i = 1:n_orgs
        crnt_significant = significance_arr(i);
        if crnt_significant == 0
            continue
        end
        % Calculate the x position for the center of the line
        x_position = [i - 1/4, i + 1/4];
        x_center = mean(x_position);
        % Place line and asterisk 10% above the max whisker
        y_position = max_whisker * 1.10;
        % Add a line connecting the two bars
        line(x_position, [y_position, y_position], 'Color', 'k', 'LineWidth', 1);
        % Add asterisks above the center of the line based on significance level
        sig_label = repmat('*', 1, crnt_significant);  % 1='*', 2='**', 3='***'
        text(x_center, y_position * 1.03, sig_label, 'HorizontalAlignment', 'center', 'FontSize', 16);
    end
    ylim([-1 max_whisker * 1.25]);
    grid on;
    hold off;
end


function plot_scatter_org(x_data, y_data, plot_title, crnt_ylabel, is_significant)
    % x_data = degron (per organism means)
    % y_data = global (per organism means)
    % assumes same length (paired by organism)
    colors = [
        0.0000 0.4470 0.7410;  % blue - Degron
        0.8500 0.3250 0.0980;  % orange - Global
    ];
    x_data = x_data(:);
    y_data = y_data(:);
    % Remove NaNs consistently (paired removal)
    valid_idx = ~isnan(x_data) & ~isnan(y_data);
    x_data = x_data(valid_idx);
    y_data = y_data(valid_idx);
    n = length(x_data);
    hold on;
    % Scatter points with small jitter
    jitter = 0.05;
    scatter(ones(n,1) + (rand(n,1)-0.5)*jitter, x_data, ...
        60, 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor','k');

    scatter(2*ones(n,1) + (rand(n,1)-0.5)*jitter, y_data, ...
        60, 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor','k');

    % Connect paired points
    for i = 1:n
        %plot([1 2], [x_data(i) y_data(i)], 'k-', 'LineWidth', 0.8);
    end
    % Mean ± SEM
    mean_x = mean(x_data);
    mean_y = mean(y_data);
    sem_x = std(x_data)/sqrt(n);
    sem_y = std(y_data)/sqrt(n);

    errorbar(1, mean_x, sem_x, 'k', 'LineWidth',2);
    errorbar(2, mean_y, sem_y, 'k', 'LineWidth',2);

    plot(1, mean_x, 'kd', 'MarkerFaceColor','k','MarkerSize',8);
    plot(2, mean_y, 'kd', 'MarkerFaceColor','k','MarkerSize',8);

    % Axis formatting
    xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],...
            'XTickLabel',{'Degron regions','Entire Protein'},...
            'FontSize',9);

    ylabel(crnt_ylabel);
    %title(plot_title);
    grid on;

    % ----- Statistical test (paired Wilcoxon) -----
    [p, h] = signrank(x_data, y_data);

    yl = ylim;
    yPos = yl(2) * 1.05;
    line([1 2], [yPos yPos], 'Color','k','LineWidth',1.2);

    if p < 0.001
        sig = '***';
    elseif p < 0.01
        sig = '**';
    elseif p < 0.05
        sig = '*';
    else
        sig = 'ns';
    end

    text(1.5, yPos*1.02, sig, ...
        'HorizontalAlignment','center','FontSize',11);

    ylim([yl(1) yPos*1.1]);
    plot([1 2], [mean_x mean_y], 'k--', 'LineWidth', 2);
    hold off;

end