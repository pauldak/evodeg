% PDB_SUPERPOSITION  Main pipeline entry point. Executes first_run and second_run.
%
% Supports two modes controlled by is_multi_evo in addresses.m:
%
%   Standard mode (is_multi_evo = 0):
%     S. cerevisiae is the reference organism. Runs first_run then second_run.
%     Total runtime approximately 35 minutes (excluding file I/O).
%
%   Multi-evolution mode (is_multi_evo = 1):
%     Any species can be the reference organism (set zero_organism in addresses.m).
%     If full_independant_cycle_run = 1, additionally runs pdbtx_all_protein_initation
%     and delete_empty_hogs before the main loop. Adds a pause after PyMOL to allow
%     all structure files to be written before second_run begins.
%     Total runtime approximately 3 hours per reference organism (with file export).
%
% File output is controlled by is_exporting_local_files in addresses.m.
% Set is_exporting_local_files = 0 for a dry run (no files written to disk).
function pdb_superposition()
    % Run addresses() then load parameters before executing the pipeline.
    clear all; clc; close all; addresses(); load('folder_paths.mat');
    warning(warning_status, 'all');
    if is_multi_evo && full_independant_cycle_run
        fprintf('Started TX initiation (~20 min), %s\n', datestr(now,'HH:MM'));
        pdbtx_all_protein_initation();
        pause(1200);
        fprintf('Deleting empty HOGs, %s\n', datestr(now,'HH:MM'));
        delete_empty_hogs();
    end
    hogs_array = first_run(); %appx 3 min
    hogs_array = second_run(hogs_array); %appx 12 min
    if is_exporting_local_files
        combine_coupling_table();
    end
    try
        fprintf('Saving hogs_array.mat (~20 min for 70 GB), %s\n', datestr(now,'HH:MM'));
        save('hogs_array.mat', 'hogs_array', '-v7.3'); %large file, takes appx 20 min to save, appx 70GB
        fprintf('Save complete, %s\n', datestr(now,'HH:MM'));
    catch Error
        disp(Error);
    end
end

% FIRST_RUN  Phase 1: interval definition, structural coupling, and PyMOL alignment.
%
% Output:
%   hogs_array - Cell array (n x 2); column 1 holds the ds struct per protein group,
%                column 2 holds the gene group ID string
%
% For each protein group (parallelised): initialises ds, loads file paths and sequences,
% skips groups where the reference organism has no PDB file, computes global Ca distance
% matrices and whole-protein RMSD, generates 32 aa sliding-window and degron intervals,
% maps coupled intervals to each comparison species via the distance matrix, and writes
% PyMOL scripts for interval-level structural alignment.
% After the parallel loop (when is_exporting_local_files = 1): combines all .pml files
% and executes PyMOL in batch mode. In multi-evolution mode, pauses to allow PyMOL
% to finish writing all structure files before second_run begins.
function hogs_array = first_run()
    load('folder_paths.mat');
    warning(warning_status, 'all');
    %initiating variables
    amount_of_files = size(data_table, 1);
    hogs_array = cell(amount_of_files, 2);
    hogs_array(:,2) = data_table.Gene_Group;
    % Uncomment the line below to regenerate global PyMOL transformation scripts (5-10 min).
    %pdbtx_all_protein_initation()
    parfor i = 1: amount_of_files %parallel
        warning('off', 'all');
        %% initiation
        ds = initiate_variables();
        ds = get_all_crnt_variables(ds, i);
        if ~isfile(ds.(['crnt_' ds.zero_organism '_pdb_path']))
            continue
        end
        ds = get_dist_matrix(ds);
        ds = get_global_rmsd_scores(ds);
        %% getting intervals - both degron and non-degron regions
        ds = amino_acid_sliding_window(ds, ds.segment_size , ds.segment_size - 1); % the -1 is for sliding window of a single aa
        ds = get_degrons_intervals(ds);
        % creating distances matrices
        ds = get_coupling_data(ds);
        create_interval_pymol_files(ds);
        hogs_array{i} = ds;
        disp(['i ' num2str(i)]);
    end
    if is_exporting_local_files
        pymol_script_filename = combine_all_pml(data_folder);
        % Execute global structural alignment via PyMOL (batch mode). Required once; approximately 2.5 hours for the full dataset.
        run_pymol_comnds(pymol_script_filename);
        if is_multi_evo
            fprintf('Pausing for PyMOL file writing (~90 min), %s\n', datestr(now,'HH:MM'));
            pause(5000);
            fprintf('Resuming, %s\n', datestr(now,'HH:MM'));
        end
    end
end

% SECOND_RUN  Phase 2: interval-level RMSD and sequence alignment scoring.
%
% Input/Output:
%   hogs_array - Cell array from first_run; each ds struct is updated in place
%
% For each protein group (parallelised): skips empty ds entries (groups skipped in
% first_run), computes Ca distance matrices for each coupled interval pair using the
% interval-aligned PDB files, calculates interval RMSD and global-projection RMSD for
% both sliding-window and degron intervals, computes local sequence alignment scores
% using the scoring matrix set by align_to_use, and exports per-organism CSV files.
function hogs_array = second_run(hogs_array)
    load('folder_paths.mat');
    %after aquring basic data - analyzing it
    hogs_size = length(hogs_array);
    parfor k = 1: hogs_size
        ds = hogs_array{k};
        if isempty(ds)
            continue
        end
        ds = dist_matrix_interval_iterator(ds);
        ds = get_interval_rmsd_scores(ds, 'all');
        ds = get_interval_rmsd_scores(ds, 'deg');
        ds = set_alignment_score(ds, 'all');
        ds = set_alignment_score(ds, 'deg');
        export_coupling_results(ds, 'all');
        export_coupling_results(ds, 'deg');
        hogs_array{k} = ds;
        disp(['k ' num2str(k)]);
    end
end