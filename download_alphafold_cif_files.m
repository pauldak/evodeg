% DIC_GENERATOR_SINGLE_USE  Downloads AlphaFold CIF files for a list of UniProt IDs.
%
% One-time script. Reads a list of UniProt accession IDs from uniprot_id.xlsx and
% downloads the corresponding AlphaFold v4 CIF files to all_cif_files/.
%
% Requires on disk: all_cif_files/uniprot_id.xlsx (column: protein_id)
% Requires loaded:  folder_paths.mat
% Produces:         all_cif_files/<id>.cif for each UniProt ID
load('folder_paths.mat');
uniprot_table = readtable("G:\My_Drive\LabData\Shir\evolution_project\all_cif_files\uniprot_id.xlsx");
uniprot_ids = uniprot_table.protein_id;
number_of_ids = length(uniprot_ids);
local_path = [main_folder 'all_cif_files\'];
for i = 1:number_of_ids
    crnt_id = uniprot_ids{i};
    tmp_url = ['https://alphafold.ebi.ac.uk/files/AF-' crnt_id '-F1-model_v4.cif'];
    tmp_name = [local_path crnt_id '.cif'];
    websave(tmp_name, tmp_url);
end

