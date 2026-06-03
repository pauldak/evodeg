% EVO_STRUCT_ALIGN  Builds all-vs-all pairwise evolutionary distance matrices from
% per-HOG coupling tables and plots structural, sequence, and functional divergence trends.
%
% Reads all per-HOG coupling Excel files from a directory, extracts average RMSD and
% alignment scores per organism pair, and assembles 7x7 distance matrices for:
%   structural (RMSD, regular alignment), sequence (alignment score, regular alignment),
%   and functional (alignment score, QCDPred degron alignment) comparisons.
% Both global and degron-specific populations are tracked. Calls plot_evolutionary_trends
% to visualise divergence from H. sapiens across organisms.
%
% Requires on disk: directory of per-HOG coupling .xlsx files with filenames encoding
%   zero_organism ($name$) and alignment_type (@type@) via regex tokens
% Requires loaded:  folder_paths.mat (organism_list, organism_list_evo_ordered,
%                   number_of_organisms, organism_list_txt)
% Produces:         figure with 3 subplots (structural, sequence, functional divergence)
function evo_struct_align()
    clc; clear all; close all;
    load('folder_paths.mat');
    % Initialise result matrices (7x7, one entry per organism pair)
    db.structural_mat_global = zeros(number_of_organisms, number_of_organisms);
    db.sequence_mat_global = zeros(number_of_organisms, number_of_organisms);
    db.functional_mat_global = zeros(number_of_organisms, number_of_organisms);
    db.structural_mat_deg = zeros(number_of_organisms, number_of_organisms);
    db.sequence_mat_deg = zeros(number_of_organisms, number_of_organisms);
    db.functional_mat_deg = zeros(number_of_organisms, number_of_organisms);
    %
    db.structural_mat_global_cell = cell(number_of_organisms, number_of_organisms);
    db.sequence_mat_global_cell = cell(number_of_organisms, number_of_organisms);
    db.functional_mat_global_cell = cell(number_of_organisms, number_of_organisms);
    db.structural_mat_deg_cell = cell(number_of_organisms, number_of_organisms);
    db.sequence_mat_deg_cell = cell(number_of_organisms, number_of_organisms);
    db.functional_mat_deg_cell = cell(number_of_organisms, number_of_organisms);

    % Path to the multi-evolution coupling tables directory.
    % Update this to the absolute path of multi_evolution_tables/ in the repository
    % before running, e.g.: directory = fullfile(pwd, 'multi_evolution_tables\');
    directory = 'D:\multi_evolution\all_tables\';
    files = dir([directory '*.xlsx']);
    path_list = fullfile({files.folder}, {files.name});
    %path_list = {[asa_library 'profile_O13329.csv']}; %debugging
    files_count = size(path_list, 2);
    for i = 1: files_count
        %readtable
        [filepath, name, ext] = fileparts(path_list{i}); % bonus
        crnt_table = readtable(path_list{i});
        crnt_table_filtered = crnt_table(crnt_table.global_rmsd_projections_by_global_rmsd < 2, :);
        %extract info from filename
        zero_organism   = regexp(name, '\$(.*?)\$', 'tokens', 'once');
        alignment_type  = regexp(name, '@(.*?)@', 'tokens', 'once');
        zero_organism  = zero_organism{1};
        alignment_type = alignment_type{1};
        if is_multi_evo
            zero_organism_index = find(strcmp(organism_list_evo_ordered, zero_organism));
        else
            zero_organism_index = find(strcmp(organism_list, zero_organism));
        end

        if contains(name, 'global_coupling')
            population = 'global';
        elseif contains(name, 'global_deg_coupling')
            population = 'degron';
        end
        %averaging
        summary_table = average_table(crnt_table_filtered);
        cell_array_grouped_raw = get_raw_grouped_data(crnt_table_filtered);
        %assigining into matrix
        if strcmp(alignment_type, 'regular_align') & strcmp(population, 'global')
            db.structural_mat_global(:,zero_organism_index) = summary_table.coupling_rmsd;
            db.sequence_mat_global(:,zero_organism_index) = summary_table.algmnt_score;
            % cell
            % Column layout of cell_array_grouped_raw: 1=organism_name, 2=coupling_rmsd,
            % 3=global_rmsd_projections, 4=global_rmsd, 5=algmnt_score, 6=ratio
            db.structural_mat_global_cell(:,zero_organism_index) = cell_array_grouped_raw(:,2); %coupling_rmsd
            db.sequence_mat_global_cell(:,zero_organism_index) = cell_array_grouped_raw(:,5); %algmnt_score

        elseif strcmp(alignment_type, 'regular_align') & strcmp(population, 'degron')
            db.structural_mat_deg(:,zero_organism_index) = summary_table.coupling_rmsd;
            db.sequence_mat_deg(:,zero_organism_index) = summary_table.algmnt_score;
            %cell
            db.structural_mat_deg_cell(:,zero_organism_index) = cell_array_grouped_raw(:,2); %coupling_rmsd
            db.sequence_mat_deg_cell(:,zero_organism_index) = cell_array_grouped_raw(:,5); %algmnt_score

        elseif strcmp(alignment_type, 'degron_align') & strcmp(population, 'global')
            db.functional_mat_global(:,zero_organism_index) = summary_table.algmnt_score;
            % cell
            db.functional_mat_global_cell(:,zero_organism_index) = cell_array_grouped_raw(:,5); %algmnt_score

        elseif strcmp(alignment_type, 'degron_align') & strcmp(population, 'degron')
            db.functional_mat_deg(:,zero_organism_index) = summary_table.algmnt_score;
            %cell
            db.functional_mat_deg_cell(:,zero_organism_index) = cell_array_grouped_raw(:,5); %algmnt_score
        end
        disp(['finished table ' num2str(i) ' out of ' num2str(files_count)]);
    end
    plot_evolutionary_trends(db)
end

function grouped_data = get_raw_grouped_data(crnt_table_filtered)
        load('folder_paths.mat');
    % List of numeric columns to average
    cols_to_hold = {
        'coupling_rmsd', ...
        'global_rmsd_projections', ...
        'global_rmsd', ...
        'algmnt_score', ...
        'global_rmsd_projections_by_global_rmsd'
    };
    [organisms, ~, group_idx] = unique(crnt_table_filtered.organism_name, 'stable');
    num_cols = numel(cols_to_hold);
    grouped_data = cell(number_of_organisms, num_cols + 1);
    grouped_data(:,1) = organism_list';
    for i = 1:number_of_organisms
        rows = (group_idx == i);
        for j = 1:num_cols
            col_name = cols_to_hold{j};
            grouped_data{i, j+1} = crnt_table_filtered{rows, col_name};  % Extract as array
        end
    end
end

function summary_table = average_table(crnt_table_filtered)
    load('folder_paths.mat');
    % List of numeric columns to average
    cols_to_average = {
        'coupling_rmsd', ...
        'global_rmsd_projections', ...
        'global_rmsd', ...
        'algmnt_score', ...
        'global_rmsd_projections_by_global_rmsd'
    };
    summary_table = groupsummary( ...
        crnt_table_filtered, ...               % Table to summarize
        'organism_name', ...                   % Grouping variable
        'mean', ...                            % Aggregation method
        cols_to_average ...                    % Columns to aggregate
    );
    summary_table.Properties.VariableNames = strrep(summary_table.Properties.VariableNames, 'mean_', '');
    %sorting
    [~, sort_idx] = ismember(summary_table.organism_name, organism_list_evo_ordered);
    [~, desired_order] = sort(sort_idx);
    summary_table = summary_table(desired_order, :);

end

function plot_evolutionary_trends(db)
    load('folder_paths.mat');
    % Index of Homo sapiens (assumed last in list)
    human_idx = find(strcmp(organism_list, 'Homo_sapiens'));
    % Organism names to compare (exclude Homo sapiens)
    comparison_orgs = strrep(organism_list, '_', ' ');
    %comparison_orgs(human_idx) = [];
    % Indices for plotting (all except Homo sapiens)
    indices = 1:length(organism_list);
    %indices(human_idx) = [];
    % Extract distances to Homo sapiens from each matrix
    struct_global_vals = db.structural_mat_global(human_idx, indices)';
    struct_deg_vals    = db.structural_mat_deg(human_idx, indices)';

    seq_global_vals = db.sequence_mat_global(human_idx, indices)';
    seq_deg_vals    = db.sequence_mat_deg(human_idx, indices)';

    func_global_vals = db.functional_mat_global(human_idx, indices)';
    func_deg_vals    = db.functional_mat_deg(human_idx, indices)';
    % Categories
    types = {'Structural', 'Sequence', 'Functional'};
    data_global = {struct_global_vals, seq_global_vals, func_global_vals};
    data_degron = {struct_deg_vals, seq_deg_vals, func_deg_vals};
    evolution_time_scale = {'1275', '...', '...', '570', '...', '...', '0'};
    for k = 1:3
        subplot(3, 1, k);
        plot(1:length(comparison_orgs), data_global{k}, '-o', 'LineWidth', 2);
        hold on;
        plot(1:length(comparison_orgs), data_degron{k}, '-o', 'LineWidth', 2);
        xticks(1:length(comparison_orgs));
        xticklabels(evolution_time_scale);
        xtickangle(45);
        ylabel('Divergence to Human');
        title([types{k} ' Divergence']);
        legend({'Global', 'Degron'}, 'Location', 'northwest');
        grid on;
        box off;
    end
 end