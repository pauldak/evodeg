% ADDRESSES  Defines all project paths and parameters; saves them to folder_paths.mat.
%
% No inputs or outputs. Run once before each session, or whenever paths or mode change.
%
% Variables saved to folder_paths.mat:
%
%   Paths:
%     source_data_folder       - repo data/ directory (auto-resolved, no user action needed)
%     main_folder, proteomes_folder, orthodb_data_folder, data_folder,
%     data_table_path, matlab_folder
%
%   Data tables:
%     data_table               - HOG-to-protein-ID mapping (hog_protein_id_table.xlsx)
%     global_coupling_table    - Pre-computed sliding-window coupling results
%     global_deg_coupling_table - Pre-computed degron coupling results
%
%   Organism lists:
%     zero_organism            - Reference species for the current run
%     organism_list            - All 7 species in standard order
%     organism_list_txt        - Human-readable version of organism_list
%     organism_list_evo_ordered     - Species in evolutionary-time order
%                                     (Nematostella and C. elegans swapped vs standard)
%     organism_list_evo_ordered_txt - Human-readable version of the above
%     number_of_organisms      - 7
%
%   Pipeline control flags:
%     is_multi_evo             - 0 = standard SC pipeline (default)
%                                1 = multi-evolution mode; set zero_organism to the
%                                    desired reference species before running
%     is_exporting_local_files - 0 = dry run (read only); 1 = write files to disk
%     align_to_use             - 1 = standard BLOSUM alignment
%                                2 = QCDPred-derived custom matrix (default)
%     full_independant_cycle_run - 1 = run pdbtx_all_protein_initation + delete_empty_hogs
%                                      at pipeline startup (multi-evolution only)
%     warning_status           - Warning suppression level passed to warning()
%                                ('off' recommended for clean parallel execution)
function addresses()
    % Defines all project-wide folder paths, organism lists, and coupling tables.
    % Saves all variables to folder_paths.mat for access across the pipeline.

    % source_data_folder: pre-computed tables and reference files bundled with the repo.
    % Resolves automatically to the data/ subdirectory next to this file.
    source_data_folder = fullfile(fileparts(mfilename('fullpath')), 'data');

    % main_folder: root of the user's local project (raw data, PDB files, proteomes, etc.).
    % Update this path to match the local installation before running.
    main_folder = ['G:\My_Drive\LabData\Shir\evolution_project\'];
    proteomes_folder = [main_folder 'proteomes\'];
    orthodb_data_folder = [main_folder 'orthodb_data\'];
    data_folder = [main_folder 'data_folder\'];
    data_table_path = fullfile(source_data_folder, 'hog_protein_id_table.xlsx');
    matlab_folder = [main_folder 'matlab\'];
    %variables
    data_table = readtable(data_table_path);
    % Alternative: coupling tables computed with the standard BLOSUM alignment matrix.
        %global_coupling_table = readtable(fullfile(source_data_folder, 'global_coupling_table_regular_alignment.xlsx'));
        %global_deg_coupling_table = readtable(fullfile(source_data_folder, 'global_deg_coupling_table_regular_alignment.xlsx'));
    % Active: coupling tables computed with the QCDPred-derived alignment matrix.
        global_coupling_table = readtable(fullfile(source_data_folder, 'global_coupling_table_degron_alignment.xlsx'));
        global_deg_coupling_table = readtable(fullfile(source_data_folder, 'global_deg_coupling_table_degron_alignment.xlsx'));
    %
    % Pipeline mode flag.
    % 0 = standard pipeline (S. cerevisiae as reference organism).
    % 1 = multi-evolution mode (change zero_organism below and re-run addresses()
    %     before each pipeline execution with a different reference species).
    is_multi_evo = 0;

    % Reference organism. In standard mode keep as Saccharomyces_cerevisiae.
    % In multi-evolution mode set to the desired reference species before running.
    zero_organism = 'Saccharomyces_cerevisiae';

    organism_list = {'Saccharomyces_cerevisiae', 'Amphimedon_queenslandica', 'Caenorhabditis_elegans',...
                            'Nematostella_vectensis', 'Drosophila_melanogaster', 'Danio_rerio', 'Homo_sapiens'};
    [organism_list_txt] = strrep(organism_list,'_',' ');
    number_of_organisms = length(organism_list);

    % Evolutionary-time-ordered organism list - Nematostella and C. elegans swapped
    % relative to organism_list to reflect structural divergence order from S. cerevisiae.
    % Used by evo_struct_align.m and plot_structural_evotime_axis.m when is_multi_evo = 1.
    organism_list_evo_ordered = {'Saccharomyces_cerevisiae', 'Amphimedon_queenslandica', ...
        'Nematostella_vectensis', 'Caenorhabditis_elegans', ...
        'Drosophila_melanogaster', 'Danio_rerio', 'Homo_sapiens'};
    organism_list_evo_ordered_txt = strrep(organism_list_evo_ordered, '_', ' ');

    % Sequence alignment scoring matrix.
    % 1 = standard BLOSUM matrix; 2 = QCDPred-derived custom matrix (default).
    align_to_use = 2;

    % Controls file output: set to 1 to write output files to disk, 0 for a read-only dry run.
    % A full run with is_exporting_local_files = 1 takes approximately 3 hours.
    is_exporting_local_files = 0;

    % Warning suppression flag (passed to warning() at pipeline start).
    warning_status = 'off';

    % When 1, pdb_superposition runs pdbtx_all_protein_initation() and delete_empty_hogs()
    % before the main pipeline loop. Only relevant when is_multi_evo = 1.
    full_independant_cycle_run = 0;

    save folder_paths...
        source_data_folder...
        main_folder...
        proteomes_folder...
        orthodb_data_folder...
        data_folder...
        data_table_path...
        matlab_folder...
        data_table...
        global_coupling_table...
        global_deg_coupling_table...
        zero_organism...
        organism_list...
        organism_list_txt...
        organism_list_evo_ordered...
        organism_list_evo_ordered_txt...
        number_of_organisms...
        is_exporting_local_files...
        is_multi_evo...
        align_to_use...
        warning_status...
        full_independant_cycle_run
end

    