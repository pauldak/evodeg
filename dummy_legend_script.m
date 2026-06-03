% DUMMY_LEGEND_SCRIPT  Code snippet that adds a manual legend to an existing figure.
%
% Creates invisible NaN plot handles with the correct marker and line styles,
% then attaches a legend for Degrons, Global regions, and Mean line.
% Intended to be pasted into a figure-generating script where the automatic
% legend is missing or incorrectly generated.
%
% Not a standalone script - paste into the relevant plotting context.

    % Same color scheme
    colors = [
        0.0000 0.4470 0.7410;  % blue - Degron
        0.8500 0.3250 0.0980;  % orange - Global
    ];

    hold on;

    % Dummy handles (NaN so nothing is drawn)
    hDeg = plot(nan, nan, 'o', ...
        'MarkerFaceColor', colors(1,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8);

    hGlob = plot(nan, nan, 'o', ...
        'MarkerFaceColor', colors(2,:), ...
        'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8);

    hMean = plot(nan, nan, 'k--', ...
        'LineWidth', 2);

    legend([hDeg, hGlob, hMean], ...
        {'Degrons', 'Global regions', 'Mean line'}, ...
        'Location', 'best');

    hold off;