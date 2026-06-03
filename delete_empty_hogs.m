% DELETE_EMPTY_HOGS  Removes HOG folders that contain no PDB file for the reference organism.
%
% Used during multi-evolution pipeline runs to clean up HOG directories where
% the zero organism had no valid AlphaFold structure. Called by pdb_superposition.m
% before the main parallel loop when running in multi-evolution mode.
%
% Requires loaded: folder_paths.mat (data_folder, zero_organism)
% WARNING: permanently deletes folders and their entire contents.
function delete_empty_hogs()
    % Load the folder paths from 'folder_paths.mat'
    load('folder_paths.mat'); % This will load 'data_folder' and 'zero_organism'
    x_folders = dir(data_folder);
    x_folders = x_folders([x_folders.isdir]);
    for i = 1:length(x_folders)
        x_folder_path = fullfile(data_folder, x_folders(i).name);
        if ~ismember(x_folders(i).name, {'.', '..'})
            zero_organism_folder = fullfile(x_folder_path, zero_organism);
            if isfolder(zero_organism_folder)
                pdb_files = dir(fullfile(zero_organism_folder, '*.pdb'));
                if isempty(pdb_files)
                    disp(['Deleting folder: ', x_folder_path]); % Display the folder being deleted (for tracking)
                    rmdir(x_folder_path, 's'); % Delete the folder and its contents
                end
            end
        end
    end
end