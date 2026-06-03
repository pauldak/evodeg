function plot_structural_evotime_axis()
% PLOT_STRUCTURAL_EVOTIME_AXIS
% Creates a reference-independent evo-time axis from a global RMSD matrix
%
% INPUTS:
%   rmsd_mat        - NxN symmetric matrix of global RMSD values
%   organism_names  - cell array of organism names (length N)
%
% OUTPUT:
%   A figure showing organisms ordered along a structural evo-time axis
%
% NOTE:
%   This is NOT a phylogenetic tree.
%   Evo-time is defined as mean structural divergence to all other organisms.
    load('folder_paths.mat');
    rmsd_mat = readmatrix(fullfile(fileparts(mfilename('fullpath')), 'data', 'pairwise_global_rmsd_matrix.xlsx'));
    % -------- compute consensus distance --------
    % Mean distance of each organism to all others (excluding self)
    mean_rmsd = zeros(number_of_organisms,1);
    for i = 1:number_of_organisms
        others = rmsd_mat(i, :);
        others(i) = NaN; % exclude self
        mean_rmsd(i) = nanmean(others);
    end

    % -------- derive consensus order --------
    % Smaller mean RMSD = structurally closer to others
    [mean_rmsd_sorted, idx] = sort(mean_rmsd, 'ascend');
    ordered_names = organism_list_evo_ordered_txt(idx);

    % -------- plotting --------
    figure('Color','w','Position',[300 200 350 500]);
    hold on;

    y = 1:number_of_organisms;
    x = mean_rmsd_sorted;

    % Plot points
    plot(x, y, 'ks', 'MarkerFaceColor','k', 'MarkerSize',6);

    % Connect with line (visual aid only)
    plot(x, y, '-', 'Color',[0.5 0.5 0.5]);

    % Labels
    set(gca, 'YTick', y, ...
             'YTickLabel', ordered_names, ...
             'YDir','reverse', ...
             'FontSize',10);

    xlabel('Mean global structural RMSD');
    title({'Structure derived consensus evo time axis', ...
           '(global regions RMSD)'}, ...
           'FontWeight','normal');

    box off;
    grid on;
    hold off;

end
