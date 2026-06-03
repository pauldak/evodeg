% TEST_OF_NEW_PLOT  Compares three methods for ordering organisms by structural divergence.
%
% Reads rmsd_mat.xlsx (7x7 global RMSD matrix) and applies three ordering algorithms
% to the 6 metazoan species (excluding S. cerevisiae):
%   1. Mean distance: sort by average RMSD to all other organisms
%   2. Seriation: greedy nearest-neighbour traversal from the most divergent organism
%   3. Fiedler ordering: spectral ordering via the second eigenvector of the graph Laplacian
% Prints a comparison table of the three orderings to the console.
%
% Requires on disk: rmsd_mat.xlsx (update path before running)
% Produces:         comparison table printed to console
%
% Development/exploration script. The chosen ordering is used in evo_struct_align.m
% and generate_evo_time_vs_parameters.m.

%% Consensus evo-time ordering from 7x7 RMSD matrix

% --- Input matrix (rows = from organism, cols = reference organism) ---
organism_names = {'Sc','Aq','Nv','Ce','Dm','Dr','Hs'};
full_names = {'Saccharomyces cerevisiae', 'Amphimedon queenslandica', ...
              'Nematostella vectensis', 'Caenorhabditis elegans', ...
              'Drosophila melanogaster', 'Danio rerio', 'Homo sapiens'};

M = readmatrix(fullfile(fileparts(mfilename('fullpath')), 'data', 'pairwise_global_rmsd_matrix.xlsx'));

%% Step 1: Symmetrize the matrix
D = (M + M') / 2;

%% Exclude Sc (index 1), work with 6 metazoans only
D_full = (M + M') / 2;
D6 = D_full(2:end, 2:end);  % 6x6 metazoan matrix
n = 6;
full_names = {'Amphimedon queenslandica', ...
              'Nematostella vectensis', 'Caenorhabditis elegans', ...
              'Drosophila melanogaster', 'Danio rerio', 'Homo sapiens'};

%% Mean distance
mean_dist = sum(D6, 2) / (n - 1);
[sorted1, idx1] = sort(mean_dist, 'ascend');

%% Seriation (greedy, start from most divergent)
visited = false(1, n);
[~, start] = max(mean_dist);
order2 = start;
visited(start) = true;
for step = 1:(n-1)
    dists = D6(order2(end), :);
    dists(visited) = Inf;
    [~, next] = min(dists);
    order2 = [order2, next];
    visited(next) = true;
end
order2 = fliplr(order2);

%% Fiedler
sigma = median(D6(D6 > 0));
W = exp(-D6.^2 / (2 * sigma^2));
W(logical(eye(n))) = 0;
Deg = diag(sum(W, 2));
L = Deg - W;
[V, eigvals] = eig(L);
[~, eig_order] = sort(diag(eigvals));
fiedler = V(:, eig_order(2));
[~, ce_idx] = max(mean_dist);  % Ce should be at high end
if fiedler(ce_idx) < mean(fiedler)
    fiedler = -fiedler;
end
[sorted3, idx3] = sort(fiedler, 'ascend');

%% Compare
fprintf('\n=== 6-Metazoan Ordering (Sc excluded) ===\n');
fprintf('%-8s %-35s %-35s %-35s\n', 'Rank', 'Mean Dist', 'Seriation', 'Fiedler');
for i = 1:n
    fprintf('%-8d %-35s %-35s %-35s\n', i, ...
        full_names{idx1(i)}, full_names{order2(i)}, full_names{idx3(i)});
end