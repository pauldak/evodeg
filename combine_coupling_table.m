% Aggregates all per-organism coupling CSV files into global Excel summary tables.
clear all; clc; close all;
load('folder_paths.mat');
amount_of_files = size(data_table, 1);
%table intiation
global_coupling_table =  initation_variables_data_gathering();
global_deg_coupling_table =  initation_variables_data_gathering();

for i = 1: amount_of_files
    ds = initiate_variables();
    ds.crnt_protein_group = data_table.Gene_Group{i};
    for j = 1: ds.number_of_organisms
        % Get the protein file name for the current organism.
        ds.cnrt_organism = ds.organism_list{j};
        crnt_folder_path = [data_folder ds.crnt_protein_group '\' ds.cnrt_organism '\'];
        files = dir([crnt_folder_path '*.fasta']);
        try
            ds.crnt_protein = files.name([1: end - 6]);
            crnt_coupling_path = [crnt_folder_path 'coupling_' ds.crnt_protein '_all.csv'];
            crnt_degron_coupling_path = [crnt_folder_path 'coupling_' ds.crnt_protein '_deg.csv'];
            crnt_coupling_table = readtable(crnt_coupling_path);
            crnt_degron_coupling_table = readtable(crnt_degron_coupling_path);
            %adding data to global tables
            global_coupling_table = add_data_to_main_tables(crnt_coupling_table, global_coupling_table, ds);
            global_deg_coupling_table = add_data_to_main_tables(crnt_degron_coupling_table, global_deg_coupling_table, ds);
        catch Error
            continue
        end
    end
    disp(i);
end
full_file_path = [main_folder 'matlab\global_coupling_table.xlsx'];
writetable(global_coupling_table, full_file_path);
full_file_path_deg = [main_folder 'matlab\global_deg_coupling_table.xlsx'];
writetable(global_deg_coupling_table, full_file_path_deg);

function crnt_table =  initation_variables_data_gathering()  
    gene_group = {};
    pritein_id = {};
    organism_name = {};
    zero_org_interval = {};
    coupling_interval = {};
    coupling_rmsd = [];
    global_rmsd_projections = [];
    global_rmsd = [];
    algmnt_score = [];
    crnt_table = table(gene_group, pritein_id, organism_name, ...
                                zero_org_interval, coupling_interval, coupling_rmsd, global_rmsd_projections,...
                                            global_rmsd, algmnt_score);
end

function global_table = add_data_to_main_tables(crnt_table, global_table, ds)
    % Gene group identifier column.
    crnt_gene_group_arr = cell(1, length(crnt_table.zero_org_interval)); 
    crnt_gene_group_arr(:) = {ds.crnt_protein_group};
    gene_group = crnt_gene_group_arr';
    % protein ID
    crnt_protein_arr = cell(1, length(crnt_table.zero_org_interval)); 
    crnt_protein_arr(:) = {ds.crnt_protein};
    pritein_id = crnt_protein_arr';
    % organism name
    crnt_organism_name_arr = cell(1, length(crnt_table.zero_org_interval)); 
    crnt_organism_name_arr(:) = {ds.cnrt_organism};
    organism_name = crnt_organism_name_arr';
    % other
    zero_org_interval = crnt_table.zero_org_interval;
    coupling_interval = crnt_table.coupling_interval;
    coupling_rmsd = crnt_table.coupling_rmsd;
    algmnt_score = crnt_table.algmnt_score;
    global_rmsd_projections = crnt_table.global_rmsd_projections;
    global_rmsd = crnt_table.global_rmsd;
    
    tmp_table = table(gene_group, pritein_id, organism_name, zero_org_interval, coupling_interval, coupling_rmsd, global_rmsd_projections,...
                                            global_rmsd, algmnt_score);
                
    global_table = vertcat(global_table ,tmp_table);
end
