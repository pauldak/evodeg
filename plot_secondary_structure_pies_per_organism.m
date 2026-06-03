% PLOT_STRUCTURAL_REGION_PIES_WITH_COUNTS_B  Per-organism pie charts with 5-group
% secondary structure classification.
%
% Reads secondary_structure_random_per_organism.xlsx and secondary_structure_degron_per_organism.xlsx (per-organism secondary structure counts),
% maps detailed structure labels to 5 groups (Alpha Helix, Beta Sheet, Turn, Bend,
% Unstructured), and plots two multi-panel figures with per-organism pies and sample sizes.
%
% Requires on disk: secondary_structure_random_per_organism.xlsx, secondary_structure_degron_per_organism.xlsx
%   Structure: columns = organisms; rows = secondary structure types + 'total_count'
% Produces:         2 figure windows (global and degron structural composition, per organism)
%
% This is the current production version of the per-organism pie chart.
% Supersedes plot_structural_region_pies.m.
function plot_secondary_structure_pies_per_organism()
    close all; clc; clear all;

    % ------------------- LOAD TABLES -------------------
    rand_tbl = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'secondary_structure_random_per_organism.xlsx'), 'PreserveVariableNames', true);
    degs_tbl = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'secondary_structure_degron_per_organism.xlsx'), 'PreserveVariableNames', true);

    % Get organism columns
    orgs = rand_tbl.Properties.VariableNames(3:end);

    % Separate totals
    rand_total_row = strcmpi(rand_tbl.Key, 'total_count');
    degs_total_row = strcmpi(degs_tbl.Key, 'total_count');
    rand_totals = table2array(rand_tbl(rand_total_row, 3:end));
    degs_totals = table2array(degs_tbl(degs_total_row, 3:end));

    % Remove total rows
    rand_tbl(rand_total_row, :) = [];
    degs_tbl(degs_total_row, :) = [];

    % ------------------- GROUPING -------------------
    groups = ["Alpha Helix", "Beta Sheet", "Turn", "Bend", "Unstructured"];

    rand_grp = map_desc_to_group(rand_tbl.Description);
    degs_grp = map_desc_to_group(degs_tbl.Description);

    rand_data = table2array(rand_tbl(:, 3:end));
    degs_data = table2array(degs_tbl(:, 3:end));

    % Accumulate per organism, per group
    n_orgs = numel(orgs);
    rand_grouped = zeros(numel(groups), n_orgs);
    degs_grouped = zeros(numel(groups), n_orgs);
    for i = 1:n_orgs
        rand_grouped(:, i) = accum_by_group(rand_grp, rand_data(:, i), groups);
        degs_grouped(:, i) = accum_by_group(degs_grp, degs_data(:, i), groups);
    end

    % ------------------- COLORS -------------------
    cmap_full = [
        0.20 0.45 0.75;   % Alpha Helix  - blue
        0.85 0.33 0.10;   % Beta Sheet   - orange
        0.49 0.18 0.56;   % Turn         - purple
        0.47 0.67 0.19;   % Bend         - green
        0.30 0.30 0.30];  % Unstructured - gray

    % Remove groups empty across ALL organisms in both conditions
    keep = any(rand_grouped > 0, 2) | any(degs_grouped > 0, 2);
    groups       = groups(keep);
    rand_grouped = rand_grouped(keep, :);
    degs_grouped = degs_grouped(keep, :);
    cmap         = cmap_full(keep, :);

    % ------------------- PLOT -------------------
    plot_multi_pie(rand_grouped, groups, orgs, rand_totals, cmap, 'Global Structural Regions');
    plot_multi_pie(degs_grouped, groups, orgs, degs_totals, cmap, 'Degron Structural Regions');
end


% ==========================================================
% Main plotting function
% ==========================================================
function plot_multi_pie(data, groups, orgs, totals, cmap, fig_title)
    n_orgs = numel(orgs);
    n_cols = ceil(n_orgs / 2);
    n_rows = 2;

    figure('Color', 'w', 'Position', [100, 100, min(300*n_cols, 1600), 650], 'Name', fig_title);
    t = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, fig_title, 'FontSize', 14, 'FontWeight', 'bold');

    legend_handles = [];  % collect patch handles for legend

    for i = 1:n_orgs
        nexttile;
        counts = data(:, i);

        % Skip organisms with no data
        if sum(counts) == 0
            axis off;
            continue;
        end

        h = pie(counts);
        patch_handles = format_pie(h, counts, cmap);  % totals NOT passed here

        if isempty(legend_handles)
            legend_handles = patch_handles;  % grab from first valid pie only
        end

        org_name = strrep(orgs{i}, '_', ' ');
        title({org_name; ' '}, 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
        text(0, -1.25, sprintf('n = %d', totals(i)), ...
            'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'normal');
    end

    % Hide unused tiles
    for i = n_orgs+1 : n_rows*n_cols
        nexttile;
        axis off;
    end

    % ---- Legend using real patch handles ----
    lgd = legend(legend_handles, groups, 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south';
    lgd.FontSize    = 9;
    lgd.Box         = 'off';
end


% ==========================================================
% Helper: group mapping  (identical to original)
% ==========================================================
function grp = map_desc_to_group(desc)
    d   = lower(string(desc));
    grp = strings(size(d));
    grp(contains(d,"unstructured") | contains(d,"disorder"))                                 = "Unstructured";
    grp(grp=="" & (contains(d,"beta")  | contains(d,"strand") | contains(d,"sheet")))        = "Beta Sheet";
    grp(grp=="" & (contains(d,"alpha") | contains(d,"helix")  | contains(d,"polyproline")))  = "Alpha Helix";
    grp(grp=="" & contains(d,"turn"))                                                         = "Turn";
    grp(grp=="")                                                                              = "Bend";
end


% ==========================================================
% Helper: accumulate counts  (identical to original)
% ==========================================================
function counts = accum_by_group(grp, vals, groups)
    counts = zeros(size(groups));
    for i = 1:numel(groups)
        counts(i) = sum(vals(grp == groups(i)));
    end
end


% ==========================================================
% Helper: beautify pie  - returns patch handles for legend
% ==========================================================
function patch_handles = format_pie(h, counts, cmap)
    patches = findobj(h, 'Type', 'Patch');
    texts   = findobj(h, 'Type', 'Text');

    % Apply fixed colors + white borders
    for i = 1:numel(patches)
        patches(i).FaceColor = cmap(i, :);
        patches(i).EdgeColor = 'w';
        patches(i).LineWidth = 0.8;
    end

    % FIX: divide by sum(counts), NOT by external total
    total = sum(counts);
    perc  = counts ./ total * 100;

    for i = 1:numel(texts)
        texts(i).String     = sprintf('%.1f%%', perc(i));
        texts(i).FontSize   = 9;
        texts(i).FontWeight = 'bold';
    end

    patch_handles = patches;  % return for legend
end