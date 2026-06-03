% CREATE_TABLE_FOR_EVO_ALIGN  Builds a flat interval table from hogs_array for evolutionary comparison.
%
% Iterates over all protein groups in hogs_array and assembles a table (evo_table) where
% each row corresponds to one reference interval, with columns for the corresponding
% coupled interval in each organism. Calls create_pml_file_evo_test to generate PyMOL
% visualisation scripts for selected interval pairs.
%
% Requires in workspace: hogs_array
% Requires loaded:       folder_paths.mat
% Requires on disk:      evolution_table_column_names.txt (for initiate_table)
% Produces:              evo_table in workspace; PyMOL scripts via create_pml_file_evo_test

%load("G:\My_Drive\LabData\Shir\evolution_project\matlab\hogs_array_regular_align.mat");
clearvars -except hogs_array
load("G:\My_Drive\LabData\Shir\evolution_project\matlab\folder_paths.mat");
number_of_hogs = length(hogs_array);
evo_table = initiate_table();
test_group = '_intervals';
%test_group = '_degrons_intervals';
seperator = '_';
first_index = 0;
end_index = 0;
for i = 1: number_of_hogs
    crnt_hog = hogs_array{i};
    crnt_hog_id = crnt_hog.crnt_protein_group;
    number_of_zero_intervals = length(crnt_hog.(['crnt_' zero_organism test_group])');
    first_index = 1 + end_index;
    end_index = first_index + number_of_zero_intervals - 1;
    for j= 1:number_of_organisms
        crnt_organisem = organism_list{j};
        try
            crnt_global.(crnt_organisem).all =  crnt_hog.(['crnt_' crnt_organisem test_group])';       
        catch
            crnt_global.(crnt_organisem).all = repmat({'-2_-2'}, 1, number_of_zero_intervals);
        end
        evo_table.(['interval_' crnt_organisem])(first_index:end_index) = crnt_global.(crnt_organisem).all;
        evo_table.('hog_id')(first_index:end_index) = crnt_hog_id;
        %add stuff
    end
end
create_pml_file_evo_test(hogs_array, evo_table)








