% GET_DEGRONS_INTERVALS  Identifies and merges degron regions from QCDPred profiles.
%
% Input:
%   ds - Data structure with prediction_threshold and the reference organism profile table
%
% Output:
%   ds - Updated data structure; ds.crnt_<zero_organism>_degrons_intervals populated
%        with interval strings in 'start_end' format
%
% Residues with logit_smooth >= prediction_threshold (0.85) are flagged; each flagged
% position is expanded by 8 residues in both directions; contiguous expanded regions
% are merged into single intervals. Degron intervals are defined on the reference
% organism (S. cerevisiae) only.
function ds = get_degrons_intervals(ds)
    amount_of_iterations = length(ds.organism_list);
    for j = 1: amount_of_iterations
        crnt_rmsd_av_vector = [];
        try
            crnt_organism = ds.organism_list{j};
            %%
            crnt_organism = ds.zero_organism; % Degron intervals are identified on the reference organism only.
            %%
            cnrt_csv_profile = ds.(['crnt_' crnt_organism '_profile_table']);
            crnt_prediction_arr = cnrt_csv_profile.logit_smooth;
            crnt_binaries = crnt_prediction_arr;
            for i =1: length(crnt_prediction_arr)
                if crnt_prediction_arr(i) >= ds.prediction_threshold
                    crnt_binaries(i-8:i+8) = 1;
                elseif crnt_binaries(i) == 1
                    continue;
                else
                    crnt_binaries(i) = 0;
                end
            end
        catch
            continue
        end
        get_ranges_from_bin_array()
        crnt_pdb_path = ds.(['crnt_' crnt_organism '_pdb_path']);
    end 
    
    
    function get_ranges_from_bin_array()
        % Initialize variables to store the ranges
        ranges = {};
        start_index = 0;
        % Iterate through the binary array
        for k = 1:length(crnt_binaries)
            if crnt_binaries(k) == 1
                if start_index == 0
                    start_index = k;  % Start of a new range
                end
            else
                if start_index > 0
                    ranges{end+1} = [start_index, k-1];  % End of the current range
                    start_index = 0;
                end
            end
        end
        % Check for the last range if it extends to the end
        if start_index > 0
            ranges{end+1} = [start_index, length(crnt_binaries)];
        end
        % Convert ranges to the desired format
        cell_array = cell(1, numel(ranges));
        for k = 1:numel(ranges)
            cell_array{k} = sprintf('%d_%d', ranges{k}(1), ranges{k}(2));
        end
        ds.(['crnt_' crnt_organism '_degrons_intervals']) = cell_array;
    end
end