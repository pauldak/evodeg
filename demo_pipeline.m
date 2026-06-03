% DEMO_PIPELINE  Runs the full MATLAB analysis pipeline on 2 pre-computed HOG datasets
% and reproduces the main figures from the manuscript.
%
% This demo covers Phases 1-4 of the pipeline. The only step NOT reproduced is the
% PyMOL structural alignment (Phase 0 / end of Phase 1), because it requires a local
% PyMOL installation. All PyMOL outputs (_tx.pdb, _tx.aln, _int_tx.pdb, _int_tx.aln)
% are provided pre-computed in demo_data/hog_data/.
%
% Required toolboxes: Bioinformatics, Statistics and Machine Learning
%
% Expected directory structure under demo_data/hog_data/:
%   <HOG_ID>/
%     <organism_name>/
%       <protein_id>.fasta              UniProt canonical sequence
%       <protein_id>.csv               QCDPred per-residue degron score profile
%       <protein_id>.pdb               AlphaFold structure (original)
%       <protein_id>_tx.pdb            Globally-aligned structure (PyMOL output)
%       <protein_id>_tx.aln            CLUSTAL global alignment (PyMOL output)
%       <protein_id>_<intv>_int_tx.pdb Interval-aligned structure (PyMOL output)
%       <protein_id>_<intv>_int_tx.aln Interval CLUSTAL alignment (PyMOL output)
%       <protein_id>_deg_<intv>_int_tx.pdb  (same for degron intervals)
%       <protein_id>_deg_<intv>_int_tx.aln
%
% Also required in demo_data/:
%   penalty_matrix_for_matlab.xlsx     QCDPred substitution matrix
%
% Usage:
%   Run from the repository root directory. Produces folder_paths.mat in the
%   working directory (demo-specific, safe to delete afterwards).

clearvars; close all; clc;

% =========================================================================
% 1. LOCATE DEMO DATA
% =========================================================================
repo_root    = fileparts(mfilename('fullpath'));
demo_root    = fullfile(repo_root, 'demo_data');
hog_data_dir = fullfile(demo_root, 'hog_data');
penalty_path = fullfile(demo_root, 'penalty_matrix_for_matlab.xlsx');

if ~exist(hog_data_dir, 'dir')
    error('demo_data/hog_data/ not found. Upload the HOG demo files first.');
end

fprintf('=== evodeg Pipeline Demo ===\n');
fprintf('HOG data: %s\n\n', hog_data_dir);

% =========================================================================
% 2. AUTO-DETECT HOG IDs AND PROTEIN IDs FROM DIRECTORY STRUCTURE
% =========================================================================
organism_list = {
    'Saccharomyces_cerevisiae', ...
    'Amphimedon_queenslandica', ...
    'Caenorhabditis_elegans', ...
    'Nematostella_vectensis', ...
    'Drosophila_melanogaster', ...
    'Danio_rerio', ...
    'Homo_sapiens'
};
organism_list_txt = strrep(organism_list, '_', ' ');
number_of_organisms = length(organism_list);
zero_organism = 'Saccharomyces_cerevisiae';

hog_entries = dir(hog_data_dir);
hog_entries = hog_entries([hog_entries.isdir] & ...
              ~ismember({hog_entries.name}, {'.', '..'}));
hog_ids = {hog_entries.name};
n_hogs  = length(hog_ids);

if n_hogs == 0
    error('No HOG subdirectories found in demo_data/hog_data/.');
end
fprintf('Found %d HOG(s): %s\n\n', n_hogs, strjoin(hog_ids, ', '));

% Build data_table: rows = HOGs, columns = Gene_Group + organisms
col_names = [{'Gene_Group'}, organism_list];
table_data = cell(n_hogs, length(col_names));

for i = 1:n_hogs
    hog_id = hog_ids{i};
    table_data{i, 1} = hog_id;
    for j = 1:number_of_organisms
        org      = organism_list{j};
        org_dir  = fullfile(hog_data_dir, hog_id, org);
        fastas   = dir(fullfile(org_dir, '*.fasta'));
        if ~isempty(fastas)
            % Protein ID is the filename without extension
            table_data{i, j+1} = fastas(1).name(1:end-6);
        else
            table_data{i, j+1} = 'NA';
        end
    end
end

data_table = cell2table(table_data, 'VariableNames', col_names);
fprintf('Data table:\n');
disp(data_table);

% =========================================================================
% 3. SET UP DEMO WORKSPACE (replaces addresses.m for demo purposes)
%    Saves a demo-specific folder_paths.mat used by pipeline functions.
% =========================================================================
main_folder      = demo_root;
data_folder      = [hog_data_dir filesep];
matlab_folder    = demo_root;
proteomes_folder = demo_root;
orthodb_data_folder = demo_root;
data_table_path  = '';

% Dummy coupling tables (not used in pipeline, only in analysis scripts)
global_coupling_table     = table();
global_deg_coupling_table = table();

is_exporting_local_files = 0;  % dry run — no file writes during demo

save('folder_paths.mat', ...
    'main_folder', 'data_folder', 'matlab_folder', 'proteomes_folder', ...
    'orthodb_data_folder', 'data_table_path', 'data_table', ...
    'global_coupling_table', 'global_deg_coupling_table', ...
    'organism_list', 'organism_list_txt', 'number_of_organisms', ...
    'zero_organism', 'is_exporting_local_files');

fprintf('Demo folder_paths.mat saved.\n\n');

% =========================================================================
% 4. LOAD QCDPred PENALTY MATRIX
% =========================================================================
if isfile(penalty_path)
    penalty_matrix = table2array(readtable(penalty_path));
    fprintf('Penalty matrix loaded (%dx%d).\n', size(penalty_matrix));
else
    penalty_matrix = [];
    fprintf('WARNING: penalty_matrix_for_matlab.xlsx not found in demo_data/.\n');
    fprintf('         Sequence alignment will use default BLOSUM50.\n');
end

% =========================================================================
% 5. PHASE 1 - Per HOG: intervals, global RMSD, coupling
% =========================================================================
fprintf('\n--- Phase 1: Interval definition and coupling ---\n');
hogs_array = cell(n_hogs, 1);

for i = 1:n_hogs
    fprintf('\nHOG %d/%d: %s\n', i, n_hogs, hog_ids{i});

    %% 5a. Initialise data structure
    ds = initiate_variables();
    if ~isempty(penalty_matrix)
        ds.qcdpred_penalty_matrix = penalty_matrix;
    end
    ds.crnt_protein_group = hog_ids{i};

    %% 5b. Populate organism-specific file paths and sequences
    ds = get_all_crnt_variables(ds, i);

    %% 5c. Compute Ca pairwise distance matrices (uses _tx.pdb files)
    ds = get_dist_matrix(ds);
    n_valid = sum(cellfun(@(f) ~isempty(ds.(['crnt_' f '_zero_dist_matrix'])), ...
                          organism_list));
    fprintf('  Distance matrices computed for %d organisms.\n', n_valid);

    %% 5d. Compute global RMSD from CLUSTAL alignment files (_tx.aln)
    ds = get_global_rmsd_scores(ds);

    %% 5e. Generate 32 aa sliding-window intervals on S. cerevisiae
    ds = amino_acid_sliding_window(ds, ds.segment_size, ds.segment_size - 1);
    n_intervals = length(ds.(['crnt_' zero_organism '_intervals']));
    fprintf('  Sliding windows: %d intervals.\n', n_intervals);

    %% 5f. Identify degron intervals from QCDPred profiles
    ds = get_degrons_intervals(ds);
    n_degrons = length(ds.(['crnt_' zero_organism '_degrons_intervals']));
    fprintf('  Degron intervals: %d.\n', n_degrons);

    %% 5g. Map intervals to each comparison species (structural coupling)
    ds = get_coupling_data(ds);

    hogs_array{i} = ds;
end

% =========================================================================
% 6. PHASE 2 - Per HOG: interval RMSD and sequence alignment
% =========================================================================
fprintf('\n--- Phase 2: Interval RMSD and alignment scoring ---\n');

for k = 1:n_hogs
    fprintf('\nHOG %d/%d: %s\n', k, n_hogs, hog_ids{k});
    ds = hogs_array{k};

    %% 6a. Compute Ca distance matrices for each coupled interval pair
    ds = dist_matrix_interval_iterator(ds);

    %% 6b. Interval-level RMSD (sliding windows)
    ds = get_interval_rmsd_scores(ds, 'all');

    %% 6c. Interval-level RMSD (degron windows)
    ds = get_interval_rmsd_scores(ds, 'deg');

    %% 6d. Sequence alignment scores (sliding windows)
    ds = set_alignment_score(ds, 'all');

    %% 6e. Sequence alignment scores (degron windows)
    ds = set_alignment_score(ds, 'deg');

    hogs_array{k} = ds;
    fprintf('  Phase 2 complete for %s.\n', hog_ids{k});
end

% =========================================================================
% 7. PHASE 3 - Assemble coupling tables from results
% =========================================================================
fprintf('\n--- Phase 3: Building coupling tables ---\n');

global_tbl = build_coupling_table(hogs_array, organism_list, zero_organism, 'all');
deg_tbl    = build_coupling_table(hogs_array, organism_list, zero_organism, 'deg');

fprintf('Global coupling table: %d rows.\n', height(global_tbl));
fprintf('Degron coupling table: %d rows.\n', height(deg_tbl));

% =========================================================================
% 8. PHASE 4 - Statistical analysis and figures
% =========================================================================
fprintf('\n--- Phase 4: Statistical analysis and figures ---\n');

comparison_idx    = 2:number_of_organisms;
comparison_list   = organism_list(comparison_idx);
comparison_labels = organism_list_txt(comparison_idx);
n_orgs            = length(comparison_list);

rmsd_global_vec  = {};   rmsd_deg_vec  = {};
align_global_vec = {};   align_deg_vec = {};
rmsd_global_mean  = zeros(1, n_orgs);
rmsd_deg_mean     = zeros(1, n_orgs);
align_global_mean = zeros(1, n_orgs);
align_deg_mean    = zeros(1, n_orgs);
rmsd_sig   = zeros(1, n_orgs);
align_sig  = zeros(1, n_orgs);

fprintf('\n%-30s  %6s  %6s  %8s\n', 'Organism', 'N_gl', 'N_dg', 'p_RMSD');

for i = 1:n_orgs
    org = comparison_list{i};

    g = global_tbl(strcmp(global_tbl.organism_name, org) & ...
                   global_tbl.coupling_rmsd ~= 99, :);
    d = deg_tbl(strcmp(deg_tbl.organism_name, org) & ...
                deg_tbl.coupling_rmsd ~= 99, :);

    if isempty(g) || isempty(d); continue; end

    rmsd_global_vec{i}  = g.coupling_rmsd;
    rmsd_deg_vec{i}     = d.coupling_rmsd;
    align_global_vec{i} = g.algmnt_score;
    align_deg_vec{i}    = d.algmnt_score;
    rmsd_global_mean(i)  = mean(g.coupling_rmsd);
    rmsd_deg_mean(i)     = mean(d.coupling_rmsd);
    align_global_mean(i) = mean(g.algmnt_score);
    align_deg_mean(i)    = mean(d.algmnt_score);

    [p_r, ~] = ranksum(g.coupling_rmsd, d.coupling_rmsd);
    [p_a, ~] = ranksum(g.algmnt_score,  d.algmnt_score);
    rmsd_sig(i)  = sig_level(p_r);
    align_sig(i) = sig_level(p_a);

    fprintf('%-30s  %6d  %6d  %8.4g\n', strrep(org,'_',' '), height(g), height(d), p_r);
end

[p_rmsd_paired,  ~] = signrank(rmsd_deg_mean,  rmsd_global_mean);
[p_align_paired, ~] = signrank(align_deg_mean, align_global_mean);
fprintf('\nPaired signed-rank — RMSD: p = %.4g   Alignment: p = %.4g\n', ...
    p_rmsd_paired, p_align_paired);

% Figure 1 - RMSD box plots
figure('Name', 'Figure 1: RMSD distributions (pipeline output)', 'Color', 'w');
plot_grouped_boxes(rmsd_deg_vec, rmsd_global_vec, ...
    'Structural RMSD — Degron vs Global (pipeline output)', ...
    'RMSD (Å)', comparison_labels, rmsd_sig);

% Figure 2 - Mean RMSD paired scatter
figure('Name', 'Figure 2: Mean RMSD per organism', 'Color', 'w');
plot_paired_scatter(rmsd_deg_mean, rmsd_global_mean, 'Mean RMSD (Å)', p_rmsd_paired);

% Figure 3 - Alignment score box plots
figure('Name', 'Figure 3: Alignment scores (pipeline output)', 'Color', 'w');
plot_grouped_boxes(align_deg_vec, align_global_vec, ...
    'Sequence Alignment Score — Degron vs Global (pipeline output)', ...
    'Alignment Score (AU)', comparison_labels, align_sig);

% Figure 4 - Mean alignment paired scatter
figure('Name', 'Figure 4: Mean alignment score per organism', 'Color', 'w');
plot_paired_scatter(align_deg_mean, align_global_mean, ...
    'Sequence Alignment Score (AU)', p_align_paired);

fprintf('\nPipeline demo complete. 4 figures produced from %d HOG(s).\n', n_hogs);
fprintf('To see figures from the full dataset, run demo.m.\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function tbl = build_coupling_table(hogs_array, organism_list, zero_organism, analyze_type)
% Assembles a flat coupling table from all HOGs and organisms in hogs_array.
    gene_group   = {};
    protein_id   = {};
    organism_name = {};
    zero_org_interval   = {};
    coupling_interval   = {};
    coupling_rmsd       = [];
    global_rmsd_projections = [];
    global_rmsd         = [];
    algmnt_score        = [];

    if strcmp(analyze_type, 'all')
        intv_field  = '_intervals';
        rmsd_field  = '_intervals_rmsd_scores';
        proj_field  = '_global_projections_zero_rmsd';
        algn_field  = '_alignment_score';
    else
        intv_field  = '_degrons_intervals';
        rmsd_field  = '_deg_intervals_rmsd_scores';
        proj_field  = '_deg_global_projections_zero_rmsd';
        algn_field  = '_degron_alignment_score';
    end

    for i = 1:length(hogs_array)
        ds  = hogs_array{i};
        hog = ds.crnt_protein_group;
        try
            zero_intvs = ds.(['crnt_' zero_organism intv_field]);
        catch
            continue
        end
        for j = 1:length(organism_list)
            org = organism_list{j};
            if strcmp(org, zero_organism); continue; end
            try
                prot_id  = ds.(['crnt_' org '_protein_id']);
                org_intvs = ds.(['crnt_' org intv_field]);
                rmsd_arr  = ds.(['crnt_' org rmsd_field]);
                proj_arr  = ds.(['crnt_' org proj_field]);
                grmsd     = ds.(['crnt_' org '_zero_global_rmsd_score']);
                algn_arr  = ds.(['crnt_' org algn_field]);
            catch
                continue
            end
            n = min([length(zero_intvs), length(org_intvs), ...
                     length(rmsd_arr), length(proj_arr), length(algn_arr)]);
            for t = 1:n
                gene_group{end+1}   = hog;
                protein_id{end+1}   = prot_id;
                organism_name{end+1} = org;
                zero_org_interval{end+1} = zero_intvs{t};
                coupling_interval{end+1} = org_intvs{t};
                coupling_rmsd(end+1)     = rmsd_arr(t);
                global_rmsd_projections(end+1) = proj_arr(t);
                global_rmsd(end+1)       = grmsd;
                algmnt_score(end+1)      = algn_arr(t);
            end
        end
    end

    tbl = table(gene_group', protein_id', organism_name', ...
                zero_org_interval', coupling_interval', ...
                coupling_rmsd', global_rmsd_projections', ...
                global_rmsd', algmnt_score', ...
        'VariableNames', {'gene_group','pritein_id','organism_name', ...
                          'zero_org_interval','coupling_interval', ...
                          'coupling_rmsd','global_rmsd_projections', ...
                          'global_rmsd','algmnt_score'});
end


function plot_grouped_boxes(x_data, y_data, ttl, ylbl, org_names, sig_arr)
% Grouped box plots: x_data = degron (blue), y_data = global (orange).
    non_empty = ~cellfun(@isempty, x_data) & ~cellfun(@isempty, y_data);
    x_data   = x_data(non_empty);
    y_data   = y_data(non_empty);
    org_names = org_names(non_empty);
    sig_arr   = sig_arr(non_empty);
    n = length(org_names);
    if n == 0; title('No data available'); return; end

    colors = [0.0000 0.4470 0.7410;
              0.8500 0.3250 0.0980];
    pos_deg  = (1:n) - 0.18;
    pos_glob = (1:n) + 0.18;

    all_deg  = vertcat(x_data{:});
    all_glob = vertcat(y_data{:});
    lbl_deg  = repelem(org_names, cellfun(@numel, x_data))';
    lbl_glob = repelem(org_names, cellfun(@numel, y_data))';

    hold on;
    boxplot(all_deg,  lbl_deg,  'Positions', pos_deg,  'Widths', 0.3, ...
        'Symbol', '', 'Colors', colors(1,:));
    boxplot(all_glob, lbl_glob, 'Positions', pos_glob, 'Widths', 0.3, ...
        'Symbol', '', 'Colors', colors(2,:));
    for i = 1:n
        plot(pos_deg(i),  mean(x_data{i}), 'k+', 'MarkerSize', 6, 'LineWidth', 1);
        plot(pos_glob(i), mean(y_data{i}), 'k+', 'MarkerSize', 6, 'LineWidth', 1);
    end

    all_vals = [all_deg; all_glob];
    ymax = prctile(all_vals, 90) + 1.5 * iqr(all_vals);
    ymax = max(all_vals(all_vals <= ymax));

    for i = 1:n
        if sig_arr(i) == 0; continue; end
        line([pos_deg(i) pos_glob(i)], [ymax*1.08 ymax*1.08], 'Color','k','LineWidth',1);
        text(i, ymax*1.11, repmat('*',1,sig_arr(i)), ...
            'HorizontalAlignment','center','FontSize',14);
    end
    set(gca,'XTick',1:n,'XTickLabel',org_names,'XTickLabelRotation',20,'FontSize',9);
    ylabel(ylbl); title(ttl);
    ylim([-0.5 ymax*1.2]); grid on; hold off;
    hd = plot(nan,nan,'s','MarkerFaceColor',colors(1,:),'MarkerEdgeColor',colors(1,:));
    hg = plot(nan,nan,'s','MarkerFaceColor',colors(2,:),'MarkerEdgeColor',colors(2,:));
    legend([hd hg],{'Degron intervals','Global intervals'},'Location','northeast');
end


function plot_paired_scatter(x_data, y_data, ylbl, p_val)
% Paired scatter: each dot = one organism mean.
    valid = x_data ~= 0 & y_data ~= 0;
    x_data = x_data(valid);
    y_data = y_data(valid);
    n = length(x_data);
    if n == 0; title('No data'); return; end

    colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980];
    hold on;
    jitter = 0.04;
    scatter(ones(n,1)+(rand(n,1)-.5)*jitter, x_data, 60, ...
        'MarkerFaceColor',colors(1,:),'MarkerEdgeColor','k');
    scatter(2*ones(n,1)+(rand(n,1)-.5)*jitter, y_data, 60, ...
        'MarkerFaceColor',colors(2,:),'MarkerEdgeColor','k');
    for i = 1:n
        plot([1 2],[x_data(i) y_data(i)],'Color',[0.6 0.6 0.6],'LineWidth',0.5);
    end
    errorbar(1,mean(x_data),std(x_data)/sqrt(n),'k','LineWidth',2);
    errorbar(2,mean(y_data),std(y_data)/sqrt(n),'k','LineWidth',2);
    plot(1,mean(x_data),'kd','MarkerFaceColor','k','MarkerSize',8);
    plot(2,mean(y_data),'kd','MarkerFaceColor','k','MarkerSize',8);
    plot([1 2],[mean(x_data) mean(y_data)],'k--','LineWidth',2);
    yl = ylim; ypos = yl(2)*1.05;
    line([1 2],[ypos ypos],'Color','k','LineWidth',1.2);
    text(1.5,ypos*1.03,sig_stars(p_val),'HorizontalAlignment','center','FontSize',11);
    ylim([yl(1) ypos*1.12]); xlim([0.5 2.5]);
    set(gca,'XTick',[1 2],'XTickLabel',{'Degron regions','Global regions'},'FontSize',10);
    ylabel(ylbl); grid on; hold off;
end


function level = sig_level(p)
    if p < 0.001; level = 3; elseif p < 0.01; level = 2;
    elseif p < 0.05; level = 1; else; level = 0; end
end

function s = sig_stars(p)
    if p < 0.001; s = '***'; elseif p < 0.01; s = '**';
    elseif p < 0.05; s = '*'; else; s = 'ns'; end
end
