% GET_COUPLING_DATA  Maps reference organism intervals to structurally equivalent regions
% in each comparison species using pairwise Ca distance matrices.
%
% Input:
%   ds - Data structure with sliding-window intervals, degron intervals, Ca distance
%        matrices, and organism metadata
%
% Output:
%   ds - Updated data structure; for each comparison organism:
%        ds.crnt_<organism>_intervals: coupled sliding-window interval strings
%        ds.crnt_<organism>_degrons_intervals: coupled degron interval strings
%
% For each reference interval boundary, searches within a ±5 residue window of the
% Ca distance matrix to find the median closest residue in the comparison protein,
% then expands by gap_coupling (8 aa) on each side. Intervals with no valid match
% receive the sentinel '-1_-1'.
function ds = get_coupling_data(ds)
   % Maps reference organism intervals to structurally equivalent regions in each comparison species
   % using the pairwise Ca distance matrix. Sliding-window and degron intervals are processed
   % independently rather than parameterised by type, as this function defines the initial coupling.
    zero_org_sliding_window_intervals =  ds.(['crnt_' ds.zero_organism '_intervals']);
    zero_org_degrons_intervals = ds.(['crnt_' ds.zero_organism '_degrons_intervals']);
    for i = 1: ds.number_of_organisms
        crnt_organism = ds.organism_list{i};
        term_zero_rganism = strcmp(crnt_organism, ds.zero_organism);
        term_no_data = strcmp(ds.(['crnt_' crnt_organism '_protein_id']), 'NA');
        if term_zero_rganism || term_no_data
            continue
        end
        %case all
        number_of_sliding_window_zero_intervals = length(zero_org_sliding_window_intervals);
        for j= 1: number_of_sliding_window_zero_intervals
            crnt_interval = zero_org_sliding_window_intervals{j};
            try
                coupler_interval = get_coupler_interval(crnt_organism, crnt_interval, ds);
            catch Error
                continue
            end
            % Intervals are stored as 'start_end' strings to maintain a consistent format across ds fields.
            interval_converted = [num2str(coupler_interval(1)) '_' num2str(coupler_interval(2))];
            ds.(['crnt_' crnt_organism '_intervals']) = [ds.(['crnt_' crnt_organism '_intervals']) interval_converted];
        end
        %case degrons
        number_of_degron_zero_intervals = length(zero_org_degrons_intervals);
        for k = 1: number_of_degron_zero_intervals
            crnt_interval = zero_org_degrons_intervals{k};
            try
                coupler_interval = get_coupler_interval(crnt_organism, crnt_interval, ds);
            catch Error
                continue
            end 
            interval_converted = [num2str(coupler_interval(1)) '_' num2str(coupler_interval(2))];
            ds.(['crnt_' crnt_organism '_degrons_intervals']) = ...
               [ds.(['crnt_' crnt_organism '_degrons_intervals']) interval_converted];
        end
     end
end



function coupler_interval = get_coupler_interval(crnt_organism, crnt_interval, ds)
    dist_matrix = ds.(['crnt_' crnt_organism '_zero_dist_matrix']);
    first_aa_position_zero = str2double(extractBefore(crnt_interval, '_'));
    last_aa_position_zero = str2double(extractAfter(crnt_interval, '_'));
    searching_intervals_gap = [-5 -3 0 3 5];
    first_aa_position_zero_search = searching_intervals_gap + first_aa_position_zero;
    last_aa_position_zero_search = searching_intervals_gap + last_aa_position_zero;
    %first position
    crnt_frsr_matching_aa_position = find_median_fit(first_aa_position_zero_search, dist_matrix);
    if isnan(crnt_frsr_matching_aa_position)
        coupler_interval = [ds.wrong_index_indicator ds.wrong_index_indicator];
        return
    end
    %last position
    crnt_lst_matching_aa_position = find_median_fit(last_aa_position_zero_search, dist_matrix);
    
    function median_fit = find_median_fit(searching_interval, dist_matrix)
        number_of_iterations = length(searching_interval);
        matching_interval = [];
        for k = 1: number_of_iterations
           crnt_matrix_row = searching_interval(k);
            if crnt_matrix_row < 1 || crnt_matrix_row > size(dist_matrix, 1)
                continue
            end
            [tmp_val, crnt_min_value_postion] = min(dist_matrix(crnt_matrix_row,:));
            if tmp_val > ds.minimal_angstrom_dist
                % Filtering by single-residue distance is intentionally omitted: endpoint
                % distances do not reliably represent the interval-average RMSD.
            end
            matching_interval = [matching_interval, crnt_min_value_postion];    
        end
            median_fit = round(median(matching_interval));
    end
    %concatenating two intervals
    crnt_frst_coupler = crnt_frsr_matching_aa_position - ds.gap_coupling;
    crnt_last_coupler = crnt_lst_matching_aa_position + ds.gap_coupling;
    if crnt_frst_coupler < 1
        crnt_frst_coupler = 1;
    end
    if crnt_last_coupler > size(dist_matrix, 2)
        crnt_last_coupler = size(dist_matrix, 2);
    end
    term_nan = isnan(crnt_frst_coupler) || isnan(crnt_last_coupler);
    term_good_order = crnt_last_coupler < crnt_frst_coupler;
    if  term_nan || term_good_order
        coupler_interval = [ds.wrong_index_indicator ds.wrong_index_indicator];
    else
        coupler_interval = [crnt_frst_coupler crnt_last_coupler];
    end
end
     