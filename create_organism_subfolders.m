% ADD_SUBFOLDERS_TMP_CODE  Creates an evo_comparison subdirectory inside each gene group folder.
%
% One-time setup script. Iterates over all subdirectories of data_folder and creates
% an 'evo_comparison' subfolder in each, if it does not already exist.
%
% Update parentFolder to match the local data_folder path before running.

parentFolder = 'G:\My_Drive\LabData\Shir\evolution_project\data_folder';
items = dir(parentFolder);
for k = 1:length(items)
    itemName = items(k).name;
    if items(k).isdir && ~strcmp(itemName, '.') && ~strcmp(itemName, '..')
        subfolderPath = fullfile(parentFolder, itemName);
        newFolderPath = fullfile(subfolderPath, 'evo_comparison');
        if ~exist(newFolderPath, 'dir')
            mkdir(newFolderPath);
        else
        end
    end
end