% GET_INTERVAL_RMSD_SCORES  Computes RMSD for each coupled interval across all organisms.
%
% Inputs:
%   ds           - Data structure with interval lists, Ca distance matrices, and alignment paths
%   analyze_type - 'all' for sliding-window intervals; 'deg' for degron intervals
%
% Output:
%   ds - Updated data structure; for each organism:
%        ds.crnt_<organism>_<type>intervals_rmsd_scores: interval-level RMSD array
%        ds.crnt_<organism>_<type>global_projections_zero_rmsd: RMSD re-computed using
%          the global alignment regions as a cross-reference metric
%
% Intervals flagged as '-1_-1' or with missing alignment files receive the
% sentinel value ds.nan_alignment_score (99).
function ds = get_interval_rmsd_scores(ds, analyze_type)
    load('folder_paths.mat');
    if strcmp(analyze_type, 'all')
        zero_intervals = ds.(['crnt_' ds.zero_organism '_intervals']);
        number_of_iterations = length(zero_intervals);
    elseif strcmp(analyze_type, 'deg')
        zero_degs = ds.(['crnt_' ds.zero_organism '_degrons_intervals']);
        number_of_iterations = length(zero_degs);
    end
    %% non degron regions    
    for j = 1: length(ds.organism_list)
        crnt_organism = ds.organism_list{j};
        try
            if strcmp(analyze_type, 'all')
                crnt_organism_intervals = ds.(['crnt_' crnt_organism '_intervals']);
                interval_char = '_';
            elseif strcmp(analyze_type, 'deg')
                crnt_organism_intervals = ds.(['crnt_' crnt_organism '_degrons_intervals']);
                interval_char = '_deg_';
            end
        catch Error
            continue
        end
        crnt_prot_id = ds.(['crnt_' crnt_organism '_protein_id']);    
        for t = 1: number_of_iterations
            crnt_interval = crnt_organism_intervals{t};
            crnt_folder_path = ds.(['crnt_' crnt_organism '_folder_path']);
            try
                if strcmp(crnt_interval, ds.nan_interval) %-1_-1 case
                    error(ds.nan_interval);
                    continue
                end
                    %getting aligned regions
                    crnt_interval_aln_path = [crnt_folder_path crnt_prot_id interval_char crnt_interval '_int_tx.aln'];
                    aligned_regions = extract_aligned_regions(crnt_interval_aln_path);
                    ds.(['crnt_' crnt_organism interval_char crnt_interval '_int_zero_aligned_regions']) = aligned_regions;               
                    %getting rmsd               
                    crnt_dist_mat = ds.(['crnt_' crnt_organism interval_char crnt_interval '_int_zero_dist_mat']);
                    crnt_rmsd_interval = get_rmsd_from_aa_array(crnt_dist_mat, aligned_regions);
                    %measuring the new alignment as projection of the global
                    global_tx_aligned_regions = ds.(['crnt_' crnt_organism '_aligned_regions']);
                    crnt_rmsd_global_interval = get_rmsd_from_aa_array(crnt_dist_mat, global_tx_aligned_regions);
                    %adding RMSD data to array
                    ds.(['crnt_' crnt_organism interval_char 'intervals_rmsd_scores']) = ...
                            [ds.(['crnt_' crnt_organism interval_char 'intervals_rmsd_scores']), crnt_rmsd_interval];
                    %adding projection data to array
                    ds.(['crnt_' crnt_organism interval_char 'global_projections_zero_rmsd']) =...
                            [ds.(['crnt_' crnt_organism interval_char 'global_projections_zero_rmsd']), crnt_rmsd_global_interval];
            catch Error
                ds.(['crnt_' crnt_organism interval_char 'intervals_rmsd_scores']) = ...
                        [ds.(['crnt_' crnt_organism interval_char 'intervals_rmsd_scores']), ds.nan_alignment_score];
                ds.(['crnt_' crnt_organism interval_char 'global_projections_zero_rmsd']) =...
                        [ds.(['crnt_' crnt_organism interval_char 'global_projections_zero_rmsd']), ds.nan_alignment_score];
                continue
            end
        end    
    end
end