% DEGRON_SCND_STRUCT_COMPARISON  Compares structurally coupled degron intervals against
% QCDPred-predicted degron regions, with secondary structure parsed directly from CIF files.
%
% Input:
%   overlap_value - Minimum fractional overlap to count an interval as matching a degron
%                   (default 0.25; the value passed as argument is overridden internally)
%
% Output:
%   optimezed_diff - Absolute difference between mean real and mean random overlap percentages
%
% Functionally equivalent to degron_alignment_comparison, but derives secondary structure
% assignments at runtime by parsing AlphaFold CIF files directly via
% parseSecondaryStructureFromCIF, rather than reading pre-computed Excel tables.
% Uses a corrected significance annotation function that takes raw counts rather than
% percentages. Exports scnd_strct_real.xlsx and scnd_strct_rand.xlsx.
%
% Requires on disk: global_deg_coupling_table_regular_alignment.xlsx
%                   data_all_degron_profiles/profile_<id>.csv
%                   all_cif_files/<id>.cif
%                   all_pdb_files/<id>.pdb
function optimezed_diff = degron_scnd_struct_comparison(overlap_value)
    overlap_value = 0.25;
    close all; clearvars -except overlap_value;
    load('folder_paths.mat'); % should define main_folder
    addpath(main_folder);   
    table_path = fullfile(source_data_folder, 'global_deg_coupling_table_regular_alignment.xlsx');
    coupling_table = readtable(table_path);
    table_length = size(coupling_table, 1);
    org_index_map = containers.Map(organism_list, 1:number_of_organisms);
    % Initialize arrays for parallel accumulation
    total_counts_tmp = zeros(table_length, 1);
    intersect_real_tmp = zeros(table_length, 1);
    scnd_strct_real = {};
    scnd_strct_rand = {};
    intersect_rand_tmp = zeros(table_length, 1);
    org_idx_tmp = zeros(table_length, 1);

    % To ensure main_folder is available inside parfor, you can:
    main_folder_const = parallel.pool.Constant(main_folder);

    for i = 1:table_length
        % Use main_folder from the constant
        mf = main_folder_const.Value;
        crnt_protein_id = coupling_table.pritein_id{i};
        crnt_organism = coupling_table.organism_name{i};
        crnt_coupling_interval = coupling_table.coupling_interval{i};
        crnt_ratio = coupling_table.coupling_rmsd(i) / coupling_table.global_rmsd_projections(i);
        crnt_is_valid_alignment = ~((crnt_ratio > 2) || (coupling_table.coupling_rmsd(i) == 99));
        if ~crnt_is_valid_alignment
            continue
        end
        org_idx = org_index_map(crnt_organism);
        org_idx_tmp(i) = org_idx;
        total_counts_tmp(i) = 1;
        crnt_csv_profile_path = [mf 'data_all_degron_profiles\profile_' crnt_protein_id];
        try
            degron_intervals = get_degron_intervals(crnt_csv_profile_path);
        catch
            continue
        end
        [isIntersecting_real, matched_degron] = ...
                            check_interval_overlap(degron_intervals, crnt_coupling_interval, overlap_value);
        if isIntersecting_real
            intersect_real_tmp(i) = 1;
            try
                tmp_scnd_strct_real = get_scnd_struct(crnt_protein_id, matched_degron);
                scnd_strct_real = [scnd_strct_real; tmp_scnd_strct_real];
            catch Error
                continue
            end
        end
        random_interval = randomize_single_interval(crnt_coupling_interval, crnt_csv_profile_path);
        isIntersecting_random = check_interval_overlap(crnt_coupling_interval, random_interval, overlap_value);
        if isIntersecting_random
            intersect_rand_tmp(i) = 1;
        else
            try 
                nums = sscanf(strrep(random_interval, '_', ' '), '%d')';
                tmp_scnd_strct_rand = get_scnd_struct(crnt_protein_id, nums);
                scnd_strct_rand = [scnd_strct_rand; tmp_scnd_strct_rand];
            catch
                continue
            end
        end
    end

    % Aggregate results per organism
    total_counts_arr = zeros(1, number_of_organisms);
    intersect_counts_real_arr = zeros(1, number_of_organisms);
    intersect_counts_random_arr = zeros(1, number_of_organisms);
    for idx = 1:table_length
        if org_idx_tmp(idx) == 0
            continue
        end
        org_i = org_idx_tmp(idx);
        total_counts_arr(org_i) = total_counts_arr(org_i) + total_counts_tmp(idx);
        intersect_counts_real_arr(org_i) = intersect_counts_real_arr(org_i) + intersect_real_tmp(idx);
        intersect_counts_random_arr(org_i) = intersect_counts_random_arr(org_i) + intersect_rand_tmp(idx);
    end
    % Prepare data for plotting and stats
    real_percents = (intersect_counts_real_arr ./ total_counts_arr) * 100;
    random_percents = (intersect_counts_random_arr ./ total_counts_arr) * 100;

    plot_real_vs_random(organism_list, real_percents, random_percents, intersect_counts_real_arr, total_counts_arr,...
        intersect_counts_random_arr);

    fprintf('Organism-wise overlap summary (Real vs Random):\n');
    for k = 1:number_of_organisms
        fprintf('%s: Real %.1f%%, Random %.1f%% (%d total)\n', ...
            organism_list{k}, real_percents(k), random_percents(k), total_counts_arr(k));
    end
    optimezed_diff = abs(mean(real_percents) - mean(random_percents));
    disp(['diff ' num2str(optimezed_diff)]);
    %export tables
    writetable(scnd_strct_real, 'secondary_structure_degron_matched.xlsx');
    writetable(scnd_strct_rand, 'secondary_structure_random_control.xlsx');
end

function plot_real_vs_random(organisms, real_vals, random_vals, intersect_counts_real_arr, total_counts_arr,...
                                intersect_counts_random_arr)
    % Calculate confidence intervals for real and random values
    n = length(organisms);
    real_CI = zeros(n, 2);
    random_CI = zeros(n, 2);
    for i = 1:n
        k_real = intersect_counts_real_arr(i);
        n_real = total_counts_arr(i);
        [low, high] = binomial_confidence_interval(k_real, n_real);
        real_CI(i,:) = [low, high] * 100;

        k_rand = intersect_counts_random_arr(i);
        n_rand = total_counts_arr(i);
        [low, high] = binomial_confidence_interval(k_rand, n_rand);
        random_CI(i,:) = [low, high] * 100;
    end

    figure('Color','w'); hold on; box on;

    x_real = (1:n) - 0.15;
    x_random = (1:n) + 0.15;

    errorbar(x_real, real_vals, real_vals' - real_CI(:,1), real_CI(:,2) - real_vals', ...
        'o', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5, 'CapSize', 10, 'MarkerFaceColor', [0 0.4470 0.7410]);

    errorbar(x_random, random_vals, random_vals' - random_CI(:,1), random_CI(:,2) - random_vals', ...
        'o', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5, 'CapSize', 10, 'MarkerFaceColor', [0.8500 0.3250 0.0980]);

    for i = 1:n
        plot([x_real(i), x_random(i)], [real_vals(i), random_vals(i)], 'k-', 'LineWidth', 0.3, 'LineStyle', '--');
    end

    % Significance using actual counts (fixed)
    add_significance_annotations(x_real, x_random, ...
        intersect_counts_real_arr, intersect_counts_random_arr, total_counts_arr);

    yl1 = yline(mean(real_vals),'-',['~' num2str(round(mean(real_vals),0)) '%'],'LineWidth',0.2, 'Color', [0.5 0.5 0.5]);
    yl2 = yline(mean(random_vals),'-',['~' num2str(round(mean(random_vals),0)) '%'],'LineWidth',0.2, 'Color', [0.5 0.5 0.5]);

    organisms_labels = {'Saccharomyces Cerevisiae', 'Amphimedon Queenslandica', 'Caenorhabditis Elegans',...
                'Nematostella Vectensis', 'Drosophila Melanogaster', 'Danio Rerio', 'Homo Sapiens'};
    set(gca, 'XTick', 1:n, 'XTickLabel', organisms_labels, 'XTickLabelRotation', 45);
    xlim([0.5, n + 0.5]);
    ylim([0 115]);
    ylabel('% Intervals Matching');
    title('QCDPred and Structural Aligned Degrons Coupling Overlap');
    legend({'Degron Intervals', 'Randomized Intervals'}, 'Location', 'northeast');

    grid on;
    set(gca, 'FontSize', 12);
end

function degron_intervals = get_degron_intervals(csv_profile)
    data = readtable(csv_profile);
    matches = data.resi(data.logit_smooth > 0.85);
    if isempty(matches)
        intervals = [];
        return;
    end
    raw_intervals = [matches - 8, matches + 8];
    raw_intervals = sortrows(raw_intervals, 1);
    merged = [];
    current = raw_intervals(1, :);
    for i = 2:size(raw_intervals, 1)
        if raw_intervals(i,1) <= current(2)
            current(2) = max(current(2), raw_intervals(i,2));
        else
            merged = [merged; current];
            current = raw_intervals(i, :);
        end
    end
    merged = [merged; current];
    degron_intervals = merged;
end

function [isIntersecting, matched_degron] = check_interval_overlap(degron_intervals, interval_str, min_overlap_fraction)
% CHECK_INTERVAL_OVERLAP checks if an interval overlaps degron intervals by at least min_overlap_fraction
%
% Inputs:
%   degron_intervals       - Nx2 matrix of intervals OR a string like '14_25'
%   interval_str           - string, e.g., '14_25'
%   min_overlap_fraction   - fraction between 0 and 1 (e.g., 0.3 for 30% overlap)
%
% Output:
%   isIntersecting         - true if at least one overlap meets threshold, false otherwise

    if nargin < 3
        disp('no min_overlap_fraction');  % Default: any overlap counts
        return % Default: 30% overlap required
    end

    % --- Handle case where first input is a string like '14_25'
    if ischar(degron_intervals) || isstring(degron_intervals)
        tokens = regexp(degron_intervals, '\d+_\d+', 'match');
        if isempty(tokens)
            error('Invalid degron interval format. Expected format like "14_25".');
        end
        nk = sscanf(tokens{1}, '%d_%d');
        degron_intervals = [nk(1), nk(2)];  % Convert to 1x2 interval
    end

    % --- Parse the query interval from interval_str
    tokens = regexp(interval_str, '\d+_\d+', 'match');
    if isempty(tokens)
        error('Invalid interval format. Expected format like "14_25".');
    end
    nk = sscanf(tokens{1}, '%d_%d');
    query_start = nk(1);
    query_end = nk(2);
    query_len = query_end - query_start;
    if query_len <= 0
        error('Invalid interval: end must be greater than start.');
    end

    % --- Check each degron interval for overlap
    isIntersecting = false;
    matched_degron = false;
    for i = 1:size(degron_intervals, 1)
        degron_start = degron_intervals(i,1);
        degron_end = degron_intervals(i,2);
        % Calculate overlap range
        overlap_start = max(query_start, degron_start);
        overlap_end   = min(query_end, degron_end);
        overlap_len   = overlap_end - overlap_start;
        if overlap_len > 0  % real overlap
            overlap_fraction = overlap_len / query_len;
            if overlap_fraction >= min_overlap_fraction
                isIntersecting = true;
                matched_degron =  degron_intervals(i,:);
                return;
            end
        end
    end
end

function [low, high] = binomial_confidence_interval(k, n)
    if n == 0
        low = NaN; high = NaN;
        return;
    end
    alpha = 0.05;
    z = norminv(1 - alpha/2);
    phat = k / n;
    denom = 1 + z^2 / n;
    center = (phat + z^2/(2*n)) / denom;
    halfwidth = z * sqrt((phat*(1 - phat) + z^2/(4*n)) / n) / denom;
    low = max(0, center - halfwidth);
    high = min(1, center + halfwidth);
end

function random_interval_str = randomize_single_interval(coupling_interval_str, csv_path)
% RANDOMIZE_SINGLE_INTERVAL randomizes a single n_k interval while preserving its length
%
% Inputs:
%   coupling_interval_str - String in the format 'start_end' (e.g., '10_40')
%   csv_path              - Path to profile_XXXXXX.csv (without .csv extension)
%
% Output:
%   random_interval_str   - String in format 'start_end' for randomized interval

    % Load the profile CSV to get the length
    T = readtable([csv_path, '.csv']);
    profile_length = height(T);

    % Parse original interval
    tokens = regexp(coupling_interval_str, '(\d+)_(\d+)', 'tokens');
    if isempty(tokens)
        error('Invalid coupling interval format: %s', coupling_interval_str);
    end
    nums = str2double(tokens{1});
    original_start = nums(1);
    original_end = nums(2);

    % Compute length
    interval_length = original_end - original_start;

    % Generate valid random start index
    max_start = profile_length - interval_length;
    if max_start < 1
        error('Protein is too short to randomize interval of length %d', interval_length);
    end
    rand_start = randi([1, max_start]);
    rand_end = rand_start + interval_length;

    % Return as string
    random_interval_str = sprintf('%d_%d', rand_start, rand_end);
end

function p_values = add_significance_annotations(x_real, x_random, real_vals, random_vals)

    n = length(real_vals);
    p_values = zeros(n,1);
    max_val = 0;
    for i = 1:n
        % Chi-squared test for 2x2 table based on counts
        observed = [real_vals(i), random_vals(i)];
        total = sum(observed);
        if total == 0
            p_values(i) = NaN;
            continue;
        end

        % Use binomial test or chi-squared approximation
        % Using Chi-squared test of proportions
        [~, p] = chi2gof([1, 2], 'Freq', observed, 'Expected', mean(observed) * [1, 1], 'Emin', 1);
        p_values(i) = p;

        % Add significance stars
        if p < 0.001
            star_str = '***';
        elseif p < 0.01
            star_str = '**';
        elseif p < 0.05
            star_str = '*';
        else
            continue;
        end

        % Find y-position for placing the star
        tmp_max = max(real_vals(i), random_vals(i));
        if tmp_max > max_val
            max_val = tmp_max;
        end
        y_pos = max_val + 5; % add offset

        % x position centered between the two points
        x_pos = mean([x_real(i), x_random(i)]);

        text(x_pos, y_pos, star_str, 'FontSize', 16, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Color', 'k');
        % text(x_pos, y_pos, '$\left\{$', ...
        %         'Interpreter', 'latex', 'FontSize', 20, ...
        %         'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
        
    end
end

function match_table = get_scnd_struct(crnt_protein_id, degron_interval)
    % Load folder path
    load('folder_paths.mat');
    local_path_cif = [main_folder 'all_cif_files\'];
    local_path_pdb = [main_folder 'all_pdb_files\'];
    filename_cif = [local_path_cif crnt_protein_id '.cif'];
    filename_pdf = [local_path_pdb crnt_protein_id '.pdb'];
    if ~isfile(filename_cif)
        error('CIF file "%s" not found.', filename_cif);
    end
    try
        crnt_pdb = pdbread(filename_pdf);
        protein_length = crnt_pdb.Sequence.NumOfResidues;
    catch
        protein_length = 0;
    end
    secStructTable = parseSecondaryStructureFromCIF(filename_cif);
    scnd_struct_table = generate_full_table(secStructTable, protein_length);
    %narrowing data
    start_res = degron_interval(1);
    end_res   = degron_interval(2);
    degron_residues = start_res:end_res;
    relevant_segment = scnd_struct_table(degron_residues, :);
    % putting all into a histogram table
    [uniqueTypes, ~, idx] = unique(relevant_segment.SecondaryStructure);
    counts = accumarray(idx, 1);
    match_table = table(uniqueTypes, counts, 'VariableNames', {'Key', 'Value'});
    % expanding the table
    match_table.protein_name = repmat({crnt_protein_id}, height(match_table), 1);
    match_table.protein_length = repmat(protein_length, height(match_table), 1);
    match_table.start_res = repmat(start_res, height(match_table), 1);
    match_table.end_res = repmat(end_res, height(match_table), 1);
end

function secStructTable = parseSecondaryStructureFromCIF(cifFile)
    % Reads a CIF file and extracts secondary structure info into a table.

    % Read file lines
    fid = fopen(cifFile, 'r');
    rawLines = {};
    tline = fgetl(fid);
    while ischar(tline)
        rawLines{end+1} = strtrim(tline); %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    % Initialize structure map
    residueStructMap = containers.Map('KeyType', 'int32', 'ValueType', 'char');

    % Regular expression to match: chain resName resNum ...
    expr = '(?<chain>\w+)\s+(?<resName>\w+)\s+(?<resNum>\d+)\s+\w+\s+\w+\s+\d+\s+(?<ss>\w+)';

    for i = 1:length(rawLines)
        line = rawLines{i};
        tokens = regexp(line, expr, 'names');
        if ~isempty(tokens)
            resNumStart = str2double(tokens(1).resNum);
            ssType = tokens(1).ss;

            % Some lines are ranges, e.g., A LEU 37 ... A ILE 53 ...
            tokensEnd = regexp(line, expr, 'once');
            if length(tokens) > 1
                resNumEnd = str2double(tokens(2).resNum);
            else
                resNumEnd = resNumStart;
            end

            for r = resNumStart:resNumEnd
                residueStructMap(r) = ssType;
            end
        end
    end

    % Build the output table
    allResidues = sort(cell2mat(residueStructMap.keys));
    residueNames = cell(length(allResidues),1);
    ssLabels = cell(length(allResidues),1);

    for i = 1:length(allResidues)
        resNum = allResidues(i);
        % Try to extract resName from raw lines (inefficient but ok for now)
        for j = 1:length(rawLines)
            if contains(rawLines{j}, sprintf(' %d ', resNum))
                % crude match, improve if needed
                parts = split(strtrim(rawLines{j}));
                idx = find(strcmp(parts, num2str(resNum)), 1);
                if idx > 1
                    residueNames{i} = parts{idx - 1};
                else
                    residueNames{i} = 'UNK';
                end
                break;
            end
        end
        ssLabels{i} = residueStructMap(resNum);
    end

    % Output table
    secStructTable = table(allResidues(:), residueNames, ssLabels, ...
        'VariableNames', {'ResidueNumber', 'ResidueName', 'SecondaryStructure'});
end

function scnd_struct_table = generate_full_table(partial_table, protein_length)    
    % Create columns
    ResidueNumber = (1:protein_length)';  % Column vector of residue indices from 1 to protein_length.
    ResidueName = repmat("unknown", protein_length, 1);          % String array of "unknown"
    SecondaryStructure = repmat("unknown", protein_length, 1);   % String array of "unknown"
    
    % Create the table
    scnd_struct_table = table(ResidueNumber, ResidueName, SecondaryStructure);
    for i = 2:height(partial_table)
            idx = partial_table.ResidueNumber(i);  % Get the ResidueNumber from partial_table
            
            % Check if idx is within the valid range
            if idx >= 1 && idx <= protein_length
                % Replace the corresponding row in scnd_struct_table
                scnd_struct_table.ResidueName(idx) = partial_table.ResidueName(i);
                scnd_struct_table.SecondaryStructure(idx) = partial_table.SecondaryStructure(i);
            end
        end
end
