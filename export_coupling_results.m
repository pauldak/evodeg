% EXPORT_COUPLING_RESULTS  Writes per-organism interval coupling data to CSV files.
%
% Inputs:
%   ds           - Data structure with computed RMSD scores, alignment scores, and interval lists
%   analyze_type - 'all' for sliding-window intervals; 'deg' for degron intervals
%
% Output:
%   CSV file written to: data_folder/<gene_group>/<organism>/coupling_<protein_id>_<type>.csv
%   Columns: zero_org_interval, coupling_interval, coupling_rmsd,
%            global_rmsd_projections, global_rmsd, algmnt_score
%
% Only executes when is_exporting_local_files = 1.
function export_coupling_results(ds, analyze_type)
    load('folder_paths.mat');
    if ~is_exporting_local_files
        return
    end
    for i = 1: length(ds.organism_list)
        try
            crnt_organism = ds.organism_list{i};
            crnt_coupling_name = ['coupling_' ds.(['crnt_' crnt_organism '_protein_id'])];
            if strcmp(analyze_type, 'all')
                zero_org_interval = ds.(['crnt_' ds.zero_organism '_intervals'])';
                coupling_interval = ds.(['crnt_' crnt_organism '_intervals'])';
                coupling_rmsd = ds.(['crnt_' crnt_organism '_intervals_rmsd_scores'])';
                global_rmsd_projections = ds.(['crnt_' crnt_organism '_global_projections_zero_rmsd'])';
                global_rmsd = repelem(ds.(['crnt_' crnt_organism '_zero_global_rmsd_score']) ,length(zero_org_interval))';
                algmnt_score = ds.(['crnt_' crnt_organism '_alignment_score'])';
            elseif strcmp(analyze_type, 'deg')
                zero_org_interval = ds.(['crnt_' ds.zero_organism '_degrons_intervals'])';
                coupling_interval = ds.(['crnt_' crnt_organism '_degrons_intervals'])';
                coupling_rmsd = ds.(['crnt_' crnt_organism '_deg_intervals_rmsd_scores'])';
                global_rmsd_projections = ds.(['crnt_' crnt_organism '_deg_global_projections_zero_rmsd'])';
                global_rmsd = repelem(ds.(['crnt_' crnt_organism '_zero_global_rmsd_score']) ,length(zero_org_interval))';
                algmnt_score = ds.(['crnt_' crnt_organism '_degron_alignment_score'])';
            end
            full_file_path = [data_folder ds.crnt_protein_group '\' crnt_organism '\' crnt_coupling_name '_' analyze_type  '.csv'];
            crnt_table = table(zero_org_interval, coupling_interval, coupling_rmsd, global_rmsd_projections,...
                                global_rmsd, algmnt_score);
            writetable(crnt_table, full_file_path);
        catch Error
            continue
        end
    end
end
