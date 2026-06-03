% GET_ALL_CRNT_VARIABLES  Populates the data structure with per-organism file paths and sequences.
%
% Inputs:
%   ds               - Data structure initialised by initiate_variables
%   data_table_index - Row index into data_table for the current protein group
%
% Output:
%   ds - Updated data structure; for each organism the following fields are added:
%        crnt_<organism>_protein_id, _fasta_path, _pdb_path, _protein_sequence,
%        _protein_length, _profile_path, _profile_table, _folder_path,
%        _global_pdtx_path, _global_aln_path, and empty result arrays for
%        intervals, RMSD scores, and alignment scores
%
% Organisms with missing files or an 'NA' protein ID entry are skipped silently.
function ds = get_all_crnt_variables(ds, data_table_index)
    load('folder_paths.mat');
    ds.crnt_protein_group = data_table.Gene_Group{data_table_index};
    for j = 1: ds.number_of_organisms
        try
            crnt_organism = ds.organism_list{j};
            ds.(['crnt_' crnt_organism '_protein_id']) = data_table.(crnt_organism){data_table_index};
            ds.(['crnt_' crnt_organism '_fasta_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\' ds.(['crnt_' crnt_organism '_protein_id']) '.fasta'];
            ds.(['crnt_' crnt_organism '_pdb_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\' ds.(['crnt_' crnt_organism '_protein_id']) '.pdb'];
            ds.(['crnt_' crnt_organism '_protein_sequence']) = ...
                    fastaread(ds.(['crnt_' crnt_organism '_fasta_path'])).Sequence;
            ds.(['crnt_' crnt_organism '_protein_length']) = ...
                    length(ds.(['crnt_' crnt_organism '_protein_sequence']));
            ds.(['crnt_' crnt_organism '_profile_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\' ds.(['crnt_' crnt_organism '_protein_id']) '.csv'];
            ds.(['crnt_' crnt_organism '_profile_table']) = ...
                    readtable(ds.(['crnt_' crnt_organism '_profile_path']));
            ds.(['crnt_' crnt_organism '_intervals']) = {};
            ds.(['crnt_' crnt_organism '_degrons_intervals']) = {};
            ds.(['crnt_' crnt_organism '_alignment_score']) = [];
            ds.(['crnt_' crnt_organism '_degron_alignment_score']) = [];
            ds.(['crnt_' crnt_organism '_intervals_rmsd_scores']) = [];
            ds.(['crnt_' crnt_organism '_deg_intervals_rmsd_scores']) = [];
            ds.(['crnt_' crnt_organism '_deg_global_projections_zero_rmsd']) = [];
            ds.(['crnt_' crnt_organism '_global_projections_zero_rmsd']) = [];
            ds.(['crnt_' crnt_organism '_zero_dist_matrix']) = [];
            ds.(['crnt_' crnt_organism '_folder_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\'];
            ds.(['crnt_' crnt_organism '_global_pdtx_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\' ds.(['crnt_' crnt_organism '_protein_id']) '_tx.pdb'];
            ds.(['crnt_' crnt_organism '_global_aln_path']) = ...
                [data_folder ds.crnt_protein_group '\' crnt_organism '\' ds.(['crnt_' crnt_organism '_protein_id']) '_tx.aln'];
            ds.(['crnt_' crnt_organism '_zero_global_rmsd_score']) = [];
            ds.(['crnt_' crnt_organism '_aligned_regions']) = [];
        catch
            continue
        end
    end
end