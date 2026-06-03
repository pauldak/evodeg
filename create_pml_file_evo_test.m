% CREATE_PML_FILE_EVO_TEST  Generates PyMOL scripts for pairwise evolutionary interval comparison.
%
% Inputs:
%   hogs_array - Cell array of ds structs from the main pipeline
%   evo_table  - Table with interval assignments per protein group and organism
%                (produced by create_table_for_evo_align)
%
% Note: This function is incomplete. Several variables referenced in the function body
% (ds, k, zero_interval, crnt_interval, interval_fid) are undefined in this scope.
% The function requires further development before it can be used.
function create_pml_file_evo_test(hogs_array, evo_table)
    load("G:\My_Drive\LabData\Shir\evolution_project\matlab\folder_paths.mat");
    number_iterations = size(evo_table, 1);
    numbers_of_comparisons = 21;
    for i = 1: number_iterations
            crnt_hog_id = evo_table.hog_id{i};
            crnt_hog = hogs_array(find(strcmp(hogs_array(:,2), crnt_hog_id)), 1);
            %
            var_names = evo_table.Properties.VariableNames;
            strct_idx = find(contains(var_names, 'strct', 'IgnoreCase', true));
            total_counts = length(strct_idx);
            %
            for j = 1:total_counts
                crnt_header_idx = strct_idx(j);
                crnt_header = var_names{crnt_header_idx};
                crnt_first_roganism = extractBefore(crnt_header, '__');
                crnt_second_roganism = extractAfter(crnt_header, '__');
                crnt_second_roganism = crnt_second_roganism(1:end-6); %removing 'strct'
                interval_first_organism = evo_table.(['interval_' crnt_first_roganism]){j};
                interval_second_organism = evo_table.(['interval_' crnt_second_roganism]){j};
                %check if measurments are relevant
                minimal_alignment_dist = 20;
                maximal_alignent_dist = 40;
                first_interval_dist = abs(str2double(extractAfter(interval_first_organism, '_'))) - ...
                                                abs(str2double(extractBefore(interval_first_organism, '_')));
                second_interval_dist = abs(str2double(extractAfter(interval_second_organism, '_'))) - ...
                                                abs(str2double(extractBefore(interval_second_organism, '_')));
                term_first = (first_interval_dist < minimal_alignment_dist) || (first_interval_dist > maximal_alignent_dist);
                term_second = (second_interval_dist < minimal_alignment_dist) || (second_interval_dist > maximal_alignent_dist);

                if term_first || term_second
                    crnt_strct_align = 99;
                else
                    crnt_organism = organism_list{k};
                    crnt_organism_folder = [data_folder ds.crnt_protein_group '\' crnt_organism '\'];
                    %initiating pymol file
                    pymol_script_filename = [crnt_organism_folder crnt_organism '_pdbtx_interval.pml'];
                    interval_fid = fopen(pymol_script_filename, 'w');
                    add_interval(zero_interval, crnt_interval, crnt_zero_pdb_path, crnt_pdb_path,...
                        crnt_organism_folder, crnt_prot_id, interval_fid);
                end




            end

    
            add_interval(zero_interval, crnt_interval, crnt_zero_pdb_path, crnt_pdb_path,...
                                    crnt_organism_folder, crnt_prot_id, interval_fid);
    
    end

end
