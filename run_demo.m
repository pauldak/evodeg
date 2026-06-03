% RUN_DEMO  Code Ocean entry point - runs demo.m and saves all figures to disk.
%
% When executed inside a Code Ocean capsule, figures are saved to /results/.
% When run locally, figures are saved to a results/ subfolder in the repo root.
%
% Called by the capsule run script. Can also be run directly in MATLAB:
%   run_demo

addpath(genpath(fileparts(mfilename('fullpath'))));

% Run the demo (creates 6 figures)
demo;

% Determine output directory
if exist('/results', 'dir')
    out_dir = '/results';           % Code Ocean capsule
else
    out_dir = fullfile(fileparts(mfilename('fullpath')), 'results');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
end

% Save every open figure as a high-resolution PNG
figs = findall(0, 'Type', 'figure');
fprintf('\nSaving %d figures to %s\n', length(figs), out_dir);

for k = 1:length(figs)
    fig  = figs(k);
    name = fig.Name;
    % Sanitise filename
    name = strrep(name, ' ', '_');
    name = regexprep(name, '[^a-zA-Z0-9_\-]', '');
    name = regexprep(name, '_+', '_');
    if isempty(name)
        name = sprintf('figure_%d', k);
    end
    out_path = fullfile(out_dir, [name '.png']);
    saveas(fig, out_path);
    fprintf('  Saved: %s\n', out_path);
end

fprintf('Done.\n');
