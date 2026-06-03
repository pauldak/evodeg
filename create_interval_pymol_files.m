% CREATE_INTERVAL_PYMOL_FILES  Writes per-organism PyMOL scripts for interval-level structural alignment.
%
% Input:
%   ds - Data structure with protein group metadata, organism folder paths,
%        sliding-window intervals, and degron intervals for all species
%
% Output:
%   ds - Returned unchanged; output written to disk as .pml scripts
%
% For each organism, generates a .pml script that loads the globally-aligned protein pair,
% selects each coupled interval, performs interval-level alignment, and saves the
% transformed structure (_int_tx.pdb), alignment file (_int_tx.aln), and image (_int_tx.png).
% Only executes when is_exporting_local_files = 1.
function ds = create_interval_pymol_files(ds)
    load('folder_paths.mat');
    if ~is_exporting_local_files
        return
    end
    number_of_hogs = size(data_table, 1);
    %creating pymol file
    crnt_zero_pdb_path = ds.(['crnt_' ds.zero_organism '_pdb_path']);
    zero_intervals = ds.(['crnt_' ds.zero_organism '_intervals']);
    zero_degs = ds.(['crnt_' ds.zero_organism '_degrons_intervals']);
    number_of_zero_intervals = length(zero_intervals);
    number_of_zero_degs = length(zero_degs);
        %% non degron regions
        for k = 1: length(organism_list)
            crnt_organism = organism_list{k};
            crnt_organism_folder = [data_folder ds.crnt_protein_group '\' crnt_organism '\'];
            %initiating pymol file
            pymol_script_filename = [crnt_organism_folder crnt_organism '_pdbtx_interval.pml'];
            interval_fid = fopen(pymol_script_filename, 'w');
            %collecting intervals
            try
                crnt_organism_intervals = ds.(['crnt_' crnt_organism '_intervals']);
            catch
                continue
            end
            crnt_prot_id = ds.(['crnt_' crnt_organism '_protein_id']);
            for t = 1: number_of_zero_intervals
                zero_interval = zero_intervals{t};
                crnt_interval = crnt_organism_intervals{t};
                try
                    crnt_pdb_path = ds.(['crnt_' crnt_organism '_global_pdtx_path']);
                    add_interval(zero_interval, crnt_interval, crnt_zero_pdb_path, crnt_pdb_path,...
                                                                            crnt_organism_folder, crnt_prot_id, interval_fid);
                catch Error
                    continue
                end
            end
            %% degron regions
            try
                crnt_organism_degs_intervals = ds.(['crnt_' crnt_organism '_degrons_intervals']);
            catch Error
                continue
            end
            for s = 1: number_of_zero_degs
                zero_degs_interval = zero_degs{s};
                crnt_degs_interval = crnt_organism_degs_intervals{s};
                try
                    crnt_prot_id_edit = [crnt_prot_id '_deg'];
                    crnt_pdb_path = ds.(['crnt_' crnt_organism '_global_pdtx_path']);
                    add_interval(zero_degs_interval, crnt_degs_interval, crnt_zero_pdb_path, crnt_pdb_path,...
                                                                            crnt_organism_folder, crnt_prot_id_edit, interval_fid);
                catch Error
                    continue
                end
            end
            fclose(interval_fid);
        end
end

function add_interval(zero_interval, crnt_interval, crnt_zero_pdb_path, crnt_pdb_path,...
                            crnt_organism_folder, crnt_prot_id, interval_fid)
    zero_interval_edit = strrep(zero_interval, '_', '-');
    crnt_interval_edit = strrep(crnt_interval, '_', '-');
    %writing comands to pml file
    fprintf(interval_fid, 'reinitialize\n');
    fprintf(interval_fid, 'load %s, zero_prot\n', crnt_zero_pdb_path);
    fprintf(interval_fid, 'load %s, crnt_prot\n', crnt_pdb_path);
    fprintf(interval_fid, 'select zero_prot_int, resi %s and chain A and zero_prot\n', zero_interval_edit);
    fprintf(interval_fid, 'select crnt_prot_int, resi %s and chain A and crnt_prot\n', crnt_interval_edit);
    fprintf(interval_fid, 'color yellow, zero_prot_int\n');
    fprintf(interval_fid, 'color red, crnt_prot_int\n');
    fprintf(interval_fid, 'align crnt_prot_int, zero_prot_int, object="crnt_align"\n');
    fprintf(interval_fid, 'zoom complete=1\n');
    fprintf(interval_fid, 'save %s%s_%s_int_tx.pdb, crnt_prot\n', crnt_organism_folder, crnt_prot_id, crnt_interval);
    fprintf(interval_fid, 'save %s%s_%s_int_tx.aln, crnt_align\n', crnt_organism_folder, crnt_prot_id, crnt_interval);
    fprintf(interval_fid, 'png %s%s_%s_int_tx, dpi=900\n', crnt_organism_folder, crnt_prot_id, crnt_interval);
    %next line is removed due to 650KB/file
    %fprintf(interval_fid, 'save %s%s_%s_int_tx.pse\n', crnt_organism_folder, crnt_prot_id, crnt_interval);
end