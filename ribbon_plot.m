% RIBBON_PLOT  Plots two data populations as ribbon plots with percentile confidence bands.
%
% Inputs:
%   rmsd_pop1  - Matrix for population 1 (rows = proteins/samples, columns = time points)
%   rmsd_pop2  - Matrix for population 2 (rows = proteins/samples, columns = time points)
%   time_array - Row vector of time points (e.g. evolutionary divergence time in MYA)
%
% Produces a figure with layered filled regions representing the 5th-95th, 25th-75th,
% and median percentile bands, overlaid with Gaussian-smoothed mean lines.
% Population 1 is shown in blue; population 2 in red.
function ribbon_plot(rmsd_pop1, rmsd_pop2, time_array)
    mean_pop1 = nanmean(rmsd_pop1, 1);
    mean_pop2 = nanmean(rmsd_pop2, 1);

    smoothed_pop1 = smoothdata(mean_pop1, 'gaussian', 5);  % Gaussian smoothing with a window of 5
    smoothed_pop2 = smoothdata(mean_pop2, 'gaussian', 5);
    % Calculate the percentiles for the confidence intervals
    conf_50_pop1 = prctile(rmsd_pop1, [5, 50, 95], 1);  % 0.5, 0.8, 0.95
    conf_50_pop2 = prctile(rmsd_pop2, [5, 50, 95], 1);
    
    % Create the plot
    figure;
    hold on;

    % Plot the ribbons for Population 1
    fill([time_array, fliplr(time_array)], [conf_50_pop1(1, :), fliplr(conf_50_pop1(3, :))], ...
        'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none'); % Light blue for 95% CI
    fill([time_array, fliplr(time_array)], [conf_50_pop1(2, :), fliplr(conf_50_pop1(3, :))], ...
        'b', 'FaceAlpha', 0.35, 'EdgeColor', 'none'); % Blue for 80% CI
    fill([time_array, fliplr(time_array)], [conf_50_pop1(2, :), fliplr(conf_50_pop1(1, :))], ...
        'b', 'FaceAlpha', 0.55, 'EdgeColor', 'none'); % Dark blue for 50% CI
    
    % Plot the ribbons for Population 2
    fill([time_array, fliplr(time_array)], [conf_50_pop2(1, :), fliplr(conf_50_pop2(3, :))], ...
        [0.8, 0.2, 0.2], 'FaceAlpha', 0.15, 'EdgeColor', 'none'); % Light red for 95% CI
    fill([time_array, fliplr(time_array)], [conf_50_pop2(2, :), fliplr(conf_50_pop2(3, :))], ...
        [0.8, 0.2, 0.2], 'FaceAlpha', 0.35, 'EdgeColor', 'none'); % Red for 80% CI
    fill([time_array, fliplr(time_array)], [conf_50_pop2(2, :), fliplr(conf_50_pop2(1, :))], ...
        [0.8, 0.2, 0.2], 'FaceAlpha', 0.55, 'EdgeColor', 'none'); % Dark red for 50% CI
    
    % Plot the smoothed curves (fit lines)
    plot(time_array, smoothed_pop1, 'b-', 'LineWidth', 2);  % Population 1
    plot(time_array, smoothed_pop2, 'r-', 'LineWidth', 2);  % Population 2
    
    % Refined labels and legend
    xlabel('Time (Million Years Ago, MYA)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('RMSD Change', 'FontSize', 12, 'FontWeight', 'bold');
    legend({'Population 1 (95% CI)', 'Population 1 (80% CI)', 'Population 1 (50% CI)', ...
            'Population 2 (95% CI)', 'Population 2 (80% CI)', 'Population 2 (50% CI)', ...
            'Smoothed Population 1', 'Smoothed Population 2'}, ...
            'Location', 'NorthEast', 'FontSize', 10, 'Box', 'off', 'FontWeight', 'bold');
    
    % Set a clean background and grid
    set(gca, 'Box', 'off', 'LineWidth', 1.5, 'FontSize', 12);
    grid on;
    
    % Refine the axis limits and ticks for better presentation
    axis([min(time_array), max(time_array), min(min(conf_50_pop1(1, :)), min(conf_50_pop2(1, :))), ...
        max(max(conf_50_pop1(3, :)), max(conf_50_pop2(3, :)))]);
    hold off;
end