% GENERATE_EVO_TIME_VS_PARAMETERS  Fits linear slopes of structural and sequence conservation
% against evolutionary divergence time across a hierarchical subset of organisms.
%
% Extracts RMSD and alignment scores from hogs_array for four hierarchically ordered species
% (S. cerevisiae, A. queenslandica, N. vectensis, H. sapiens). Filters out high-variance
% proteins (std > 1 across species), optionally normalises values, and fits a linear slope
% per protein against evolutionary time (MYA). Calls ribbon_plot to visualise results.
%
% Requires in workspace: hogs_array
% Requires loaded:       folder_paths.mat
% Produces:              ribbon plot figures; mean, std, and CV of slopes printed to console

    clc;
    clearvars -except hogs_array
    load('folder_paths.mat');
    % load([matlab_folder 'hogs_array_regular_align.mat']); % Uncomment to load a pre-computed hogs_array (~70 GB).
    number_of_hogs = length(hogs_array);
    variableTypes = {'double', 'double', 'double', 'double', 'double', 'double', 'double'};
    table_names = {'intervals_rmsd_scores', 'alignment_score', 'deg_intervals_rmsd_scores', 'degron_alignment_score'};
    organism_list = organism_list([1, 2, 4, 3, 5, 6, 7]); % reordering by evolution time
    %creating tables
    for i = 1:length(table_names)
        eval(['table_' table_names{i} ' = table(''Size'', [0, number_of_organisms], ''VariableNames'', organism_list, ''VariableTypes'', variableTypes);']);
    end
    for i =1: number_of_hogs
        crnt_hog = hogs_array{i};
        number_of_segments = length(crnt_hog.(['crnt_' zero_organism '_intervals']));
        %creating empty cell arrays for tmp data collection
        for p = 1:length(table_names)
            eval(['tmp_mtrx_'  table_names{p} ' = [];']);
        end
        for k = 1: length(organism_list)
            crnt_organism = organism_list{k};

            for p = 1:length(table_names)
                try
                    eval(['crnt_'  table_names{p} ' = crnt_hog.crnt_'  crnt_organism '_' table_names{p} ''';']);
                catch Error
                    tmp_size = length(eval(['crnt_hog.crnt_' zero_organism '_' table_names{p}]));
                    eval(['crnt_'  table_names{p} ' = nan(1,' num2str(tmp_size) ')'';']);
                end
            end
            % updating tmp matrices
            for p = 1:length(table_names)
                try
                    eval(['tmp_mtrx_'  table_names{p} ' = [tmp_mtrx_'  table_names{p} ', crnt_' table_names{p} '];']);
                catch Error_b
                    tmp_size = length(eval(['crnt_hog.crnt_' zero_organism '_' table_names{p}]));
                    eval(['crnt_' table_names{p} ' = nan(1, ' num2str(tmp_size) ')'';']);
                    eval(['tmp_mtrx_'  table_names{p} ' = [tmp_mtrx_'  table_names{p} ', crnt_' table_names{p} '];']);
                    continue
                end
            end
        end
        %updating to tables
        for p = 1:length(table_names)
            try
                eval(['table_'  table_names{p} ' = [table_'  table_names{p} '; array2table(tmp_mtrx_' table_names{p}...
                                            ', ''VariableNames'', table_'  table_names{p} '.Properties.VariableNames)];']);
            catch Error_b
                disp(Error_b)
                continue
            end
                
        end
    end

    %initiation_variables
    relevant_hirarchial_organisms =  {'Saccharomyces_cerevisiae', 'Amphimedon_queenslandica',...
        'Nematostella_vectensis', 'Homo_sapiens'};
    time_array = [0,517,560,846];
    scaled_time_array = scale_time_arrat(time_array);
    %calling analyzis
    array_intervals_rmsd_scores = set_table_data(table_intervals_rmsd_scores, 0, relevant_hirarchial_organisms, scaled_time_array);
    array_deg_intervals_rmsd_scores = set_table_data(table_deg_intervals_rmsd_scores, 0, relevant_hirarchial_organisms, scaled_time_array);
    array_alignment_score = set_table_data(table_alignment_score, 1, relevant_hirarchial_organisms, scaled_time_array);
    array_degron_alignment_score = set_table_data(table_degron_alignment_score, 1, relevant_hirarchial_organisms, scaled_time_array);
    %calling    
    ribbon_plot(array_intervals_rmsd_scores(:, 1:end-1), array_deg_intervals_rmsd_scores(:, 1:end-1), time_array);
    ribbon_plot(array_alignment_score, array_degron_alignment_score);
    

function data_array = set_table_data(table_data, is_normalize_max_min, relevant_hirarchial_organisms, time_array)
    % filter non hirarchial organism
    relevant_hirarchial_organisms_b = relevant_hirarchial_organisms(2:end);
    table_data_hirarchial = table_data(:, relevant_hirarchial_organisms);
    %filter those with high stdv
    if is_normalize_max_min == 0
        for q = 1: size(table_data_hirarchial, 1)
            crnt_row = table_data_hirarchial(q, relevant_hirarchial_organisms_b);
            crnt_stdv = std(crnt_row{:, :});
            if crnt_stdv > 1
                nanCellArray = cell(1, length(relevant_hirarchial_organisms_b));
                nanCellArray(:) = {NaN};
                table_data_hirarchial(q, relevant_hirarchial_organisms_b) = nanCellArray;
            end
        end
    end
    %preparing for ploting
    data_array = table2array(table_data_hirarchial);
    %convert 99 values
    data_array(data_array == 99) = NaN;
    % cleaning nans
    data_array = data_array(~any(isnan(data_array), 2), :);
    % min max normalization
    if is_normalize_max_min == 1
        for q = 1: size(data_array, 1)
            max_value = max(data_array(q, :));
            min_value = min(data_array(q, :));
            number_of_columns = size(data_array, 2);
            for r = 1: number_of_columns
                crnt_value = data_array(q, r);
                tmp_new_value = (crnt_value - min_value) / (max_value - min_value);
                data_array(q, r) = tmp_new_value;
            end
        end
    end
    data_array = get_slopes_column(data_array, time_array);
    disp_results(data_array);
end


function new_data_table = get_slopes_column(data_array, time_array)
    [numRows, ~] = size(data_array); 
    slopes = zeros(numRows, 1); 
    %time_array = log10(time_array); time_array(1) = 0;
    for i = 1:numRows
        yValues = data_array(i, :);
        p = polyfit(time_array, yValues, 1);  
        slopes(i) = p(1);
    end
    new_data_table = [data_array, slopes];
end

function disp_results(data_array)
    tmp_mean = mean(data_array(:, end));
    tmp_std = std(data_array(:, end));
    tmp_cv = tmp_std/tmp_mean;
    disp(['test results']);
    disp(['average - ' num2str(tmp_mean)]);
    disp(['stdv - ' num2str(tmp_std)]);
    disp(['cv - ' num2str(tmp_cv)]);
end

function scaled_time_array = scale_time_arrat(time_array)
    max_value = max(time_array);
    min_value = min(time_array);
    number_of_columns = length(time_array);
    for r = 1: number_of_columns
        crnt_value = time_array(r);
        tmp_new_value = (crnt_value - min_value) / (max_value - min_value);
        time_array(r) = tmp_new_value;
        scaled_time_array = time_array;
    end
end