% RUN_PYMOL_COMNDS  Executes a PyMOL script in batch (headless) mode.
%
% Input:
%   pymol_file - Path to the .pml script file to execute
%
% Update pymol_adrs to match the local PyMOL installation path before use.
% Raises an error if PyMOL exits with a non-zero status code.
function run_pymol_comnds(pymol_file)
    pymol_adrs = '"C:\Users\SEA\AppData\Local\Schrodinger\PyMOL2\PyMOLWin.exe"';
    % Execute the PyMOL script
    pymol_command = [pymol_adrs ' -cq ' pymol_file];
    status = system(pymol_command);
    if status == 0
        disp('PyMOL opened with the PML script.');
    else
        error('Failed to open PyMOL with the PML script.');
    end
end