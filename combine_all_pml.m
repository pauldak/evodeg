% COMBINE_ALL_PML  Concatenates all .pml PyMOL scripts under a directory tree into one file.
%
% Input:
%   data_folder - Root directory searched recursively for .pml files
%
% Output:
%   pymol_script_filename - Path to the combined output file
%                           ('create_pdbtx_interval.pml' in the working directory)
function pymol_script_filename = combine_all_pml(data_folder)
    pymol_script_filename = 'create_pdbtx_interval.pml';
    % Combine all .pml files found in folderPath and its subfolders into a single file named outputFileName.
    % Ensure the folder path ends with a slash
    if ~endsWith(data_folder, '/')
        data_folder = [data_folder '/'];
    end
    % Search for all .pml files in the folder and its subfolders
    pmlFiles = dir(fullfile(data_folder, '**', '*.pml'));
    % Open the output file for writing
    outputFileID = fopen(pymol_script_filename, 'w');
    if outputFileID == -1
        error('Cannot open output file: %s', pymol_script_filename);
    end
    % Loop through each .pml file and append its content to the output file
    for k = 1:length(pmlFiles)
        % Construct full path to the .pml file
        pmlFilePath = fullfile(pmlFiles(k).folder, pmlFiles(k).name);
        % Open the .pml file for reading
        fileID = fopen(pmlFilePath, 'r');
        if fileID == -1
            warning('Cannot open file: %s', pmlFilePath);
            continue; % Skip this file and move to the next
        end
        % Read the content of the file
        fileContent = fread(fileID, '*char')';
        % Write the content to the output file
        fwrite(outputFileID, fileContent);        
        % Close the .pml file
        fclose(fileID);
    end
    % Close the output file
    fclose(outputFileID);    
    fprintf('All .pml files have been combined into %s\n', pymol_script_filename);
end