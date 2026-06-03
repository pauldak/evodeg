% ORTHODB_WINDOWS_FOLDER_GENERATOR  Creates the output directory structure for all protein groups.
%
% One-time setup script. For each OrthoDB FASTA file in orthodb_data_folder, creates a
% corresponding directory under data_folder with one subdirectory per organism.
%
% Requires on disk: orthodb_data_folder with one .fasta file per gene group (HOG ID)
% Requires loaded:  folder_paths.mat
% Produces:         nested directories: data_folder/<gene_group>/<organism>/
load('folder_paths.mat');

files = dir([orthodb_data_folder '*.fasta']);
%get a list of all fasta files of the pdb results
path_list = fullfile({files.folder}, {files.name});
%path_list = {[asa_library 'profile_O13329.csv']}; %debugging
files_count = size(path_list, 2);

for i = 1: files_count
    [filepath, name, ext] = fileparts(path_list{i}); % bonus
    create_windows_folders(name)
end

function create_windows_folders(folder_name)
    load('folder_paths.mat');
    full_path_name = [data_folder folder_name];
    %subfolders
    subfolder1f = [full_path_name '\' 'Amphimedon_queenslandica'];
    subfolder2f = [full_path_name '\' 'Caenorhabditis_elegans'];
    subfolder3f = [full_path_name '\' 'Danio_rerio'];
    subfolder4f = [full_path_name '\' 'Drosophila_melanogaster'];
    subfolder5f = [full_path_name '\' 'Homo_Sapiens'];
    subfolder6f = [full_path_name '\' 'Nematostella_vectensis'];
    subfolder7f = [full_path_name '\' 'Saccharomyces_cerevisiae'];
    mkdir(full_path_name);
    mkdir(subfolder1f);
    mkdir(subfolder2f);
    mkdir(subfolder3f);
    mkdir(subfolder4f);
    mkdir(subfolder5f);
    mkdir(subfolder6f);
    mkdir(subfolder7f);
end