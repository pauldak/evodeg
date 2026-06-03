% PLOT_STRUCTURAL_REGION_PIES_AVG_VRS_B  Averaged pie charts with 5-group classification
% and chi-squared test of independence.
%
% Reads secondary_structure_random_summary.xlsx and secondary_structure_degron_summary.xlsx, maps detailed secondary structure labels into
% 5 canonical groups (Alpha Helix, Beta Sheet, Turn, Bend, Unstructured), plots two
% averaged pie charts (global vs degron), and reports a chi-squared test statistic.
%
% Requires on disk: secondary_structure_random_summary.xlsx, secondary_structure_degron_summary.xlsx
%   Columns: Key, Description, count
% Produces:         1 figure window with 2 subplots + chi-squared p-value annotation
%
% This is the current production version of the averaged pie chart.
% Supersedes plot_structural_region_pies_avg.m.
function plot_secondary_structure_pies_averaged()
    close all; clc;
% ------------------- LOAD TABLES -------------------
    rand_tbl = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'secondary_structure_random_summary.xlsx'), 'PreserveVariableNames', true);
    degs_tbl = readtable(fullfile(fileparts(mfilename('fullpath')), 'data', 'secondary_structure_degron_summary.xlsx'), 'PreserveVariableNames', true);
    rand_tbl(strcmpi(rand_tbl.Key,'total'),:) = [];
    degs_tbl(strcmpi(degs_tbl.Key,'total'),:) = [];
% ------------------- GROUPING -------------------
    groups = ["Alpha Helix","Beta Sheet","Turn","Bend","Unstructured"];
    rand_grp = map_desc_to_group(rand_tbl.Description);
    degs_grp = map_desc_to_group(degs_tbl.Description);
    rand_counts = accum_by_group(rand_grp, rand_tbl.count, groups);
    degs_counts = accum_by_group(degs_grp, degs_tbl.count, groups);
    rand_total = sum(rand_counts);
    degs_total = sum(degs_counts);
% ------------------- CHI-SQUARED TEST -------------------
    observed = [rand_counts; degs_counts];
    [~, chi2_p, chi2_stat] = chi2test_2xN(observed);
% ------------------- REMOVE EMPTY GROUPS -------------------
    keep = (rand_counts + degs_counts) > 0;
    groups = groups(keep);
    rand_counts = rand_counts(keep);
    degs_counts = degs_counts(keep);
% ------------------- COLORS -------------------
    cmap = [
        0.20 0.45 0.75;   % Alpha - blue
        0.85 0.33 0.10;   % Beta - orange
        0.49 0.18 0.56;   % Turn - purple
        0.47 0.67 0.19;   % Bend - green
        0.30 0.30 0.30];  % Unstructured - gray
    cmap = cmap(keep,:);
% ------------------- FIGURE -------------------
    figure('Color','w','Position',[200 200 1100 550]);
    t = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
% ------------------- RANDOM -------------------
    nexttile
    h1 = pie(rand_counts);
    format_pie(h1, rand_counts, rand_total, cmap);
    title({'Random Structural Regions'; ' '},'FontSize',11,'FontWeight','bold');
    text(0, -1.25, sprintf('n = %d', rand_total), 'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','normal');
% ------------------- DEGRON -------------------
    nexttile
    h2 = pie(degs_counts);
    format_pie(h2, degs_counts, degs_total, cmap);
    title({'Degron Structural Regions', ' '},'FontSize',11,'FontWeight','bold');
    text(0, -1.25, sprintf('n = %d', degs_total), 'HorizontalAlignment','center', 'FontSize',10, 'FontWeight','normal');
% ------------------- P-VALUE ANNOTATION -------------------
    if chi2_p < 0.001
        p_str = sprintf('\\chi^2 = %.1f, p < 0.001', chi2_stat);
    else
        p_str = sprintf('\\chi^2 = %.1f, p = %.3f', chi2_stat, chi2_p);
    end
    annotation('textbox', [0.3 0.02 0.4 0.06], ...
        'String', p_str, ...
        'HorizontalAlignment','center', 'FontSize', 10, ...
        'EdgeColor','none', 'FitBoxToText','on', ...
        'Interpreter','tex');
% ------------------- LEGEND -------------------
    lgd = legend(groups,'Orientation','horizontal');
    lgd.Layout.Tile = 'south';
    lgd.FontSize = 9;
    lgd.Box = 'off';
end
% =========================
% Helper: group mapping
% =========================
function grp = map_desc_to_group(desc)
    d = lower(string(desc));
    grp = strings(size(d));
    grp(contains(d,"unstructured") | contains(d,"disorder")) = "Unstructured";
    grp(grp=="" & (contains(d,"beta") | contains(d,"strand") | contains(d,"sheet"))) = "Beta Sheet";
    grp(grp=="" & (contains(d,"alpha") | contains(d,"helix") | contains(d,"polyproline"))) = "Alpha Helix";
    grp(grp=="" & contains(d,"turn")) = "Turn";
    grp(grp=="") = "Bend";
end
% =========================
% Helper: accumulate counts
% =========================
function counts = accum_by_group(grp, vals, groups)
    counts = zeros(size(groups));
    for i = 1:numel(groups)
        counts(i) = sum(vals(grp==groups(i)));
    end
end
% =========================
% Helper: beautify pie
% =========================
function format_pie(h, counts, total, cmap)
    patches = findobj(h,'Type','Patch');
    texts   = findobj(h,'Type','Text');
    for i = 1:numel(patches)
        patches(i).FaceColor = cmap(i,:);
        patches(i).EdgeColor = 'w';
        patches(i).LineWidth = 0.8;
    end
    perc = counts ./ total * 100;
    for i = 1:numel(texts)
        texts(i).String = sprintf('%.1f%%', perc(i));
        texts(i).FontSize = 10;
        texts(i).FontWeight = 'bold';
    end
end
% =========================
% Helper: chi-squared test for 2×N table
% =========================
function [h, p, chi2stat] = chi2test_2xN(observed)
    row_sum = sum(observed, 2);
    col_sum = sum(observed, 1);
    grand   = sum(observed(:));
    expected = row_sum * col_sum / grand;
    chi2stat = sum((observed - expected).^2 ./ expected, 'all');
    df = (size(observed,1)-1) * (size(observed,2)-1);
    p = 1 - chi2cdf(chi2stat, df);
    h = p < 0.05;
end