% Downloads PDB, CIF, FASTA, and degron profile files for each ortholog identified in the OrthoDB dataset.
load('folder_paths.mat');
global  gene_group crnt_index my_table

files = dir([orthodb_data_folder '*.fasta']);
%get a list of all fasta files of the pdb results
path_list = fullfile({files.folder}, {files.name});
%path_list = {[asa_library 'profile_O13329.csv']}; %debugging
files_count = size(path_list, 2);

my_table = cell2table(cell(files_count, 8));
organism_names = {'Gene_Group', 'Saccharomyces_cerevisiae', 'Amphimedon_queenslandica', 'Nematostella_vectensis',...
   'Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Danio_rerio', 'Homo_sapiens'};
my_table.Properties.VariableNames = organism_names;


for i = 1: files_count
    [filepath, gene_group, ext] = fileparts(path_list{i}); % bonus
    blast_cif_downloader2([filepath '\' gene_group ext])
    crnt_index  = i ;
end
%writetable(my_table, [main_folder_path 'hog_protein_id_table.xlsx']);

function blast_cif_downloader2(orthodb_file_path)
    global organism_name gene_group crnt_index my_table
    all_lines_ortho = readlines(orthodb_file_path);
    number_of_lines = length(all_lines_ortho);
    try
        for j = 0: (number_of_lines / 3)
            try
                crnt_count = 2 + j * 3;
                crnt_sequence = all_lines_ortho{crnt_count};
                organism_name = get_organism_name(all_lines_ortho{crnt_count - 1});
                proteome_path = get_proteome_path();
                crnt_fasta_line = findSequenceInFasta(proteome_path, crnt_sequence(1:55));
                uniprot_id = get_uniprot_id(crnt_fasta_line);
                save_pdb(uniprot_id);
                save_cif(uniprot_id)
                fasta_path = save_fasta(uniprot_id);
                save_degron_profile(fasta_path, uniprot_id);
                if ~isempty(my_table.(organism_name){crnt_index})
                    my_table.(organism_name){crnt_index} = 'XXX';
                else
                    my_table.(organism_name){crnt_index} = uniprot_id;
                end
                my_table.Gene_Group{crnt_index} = gene_group;
            catch
                continue
            end
        end
    catch
        return
    end
end

function oraginsm_name = get_organism_name(crnt_line)
    pattern = '"organism_name":"([^"]+)"';
    match = regexp(crnt_line, pattern, 'tokens', 'once');
    % Extract the organism name
    if ~isempty(match)
        oraginsm_name = match{1};
    else
        disp('Organism name not found.');
        oraginsm_name = '';
    end
end

function proteome_path = get_proteome_path()
	load('folder_paths.mat');
    global organism_name
    if contains(organism_name, 'Saccharomyces cerevisiae')
        proteome_path = [main_folder 'proteomes\uniprotkb_Saccharomyces_cerevisiae.fasta'];
        organism_name = 'Saccharomyces_cerevisiae';
    elseif contains(organism_name, 'Danio rerio')
        proteome_path = [main_folder 'proteomes\uniprotkb_Danio_rerio.fasta'];
        organism_name = 'Danio_rerio';
    elseif contains(organism_name, 'Caenorhabditis elegans')
        proteome_path = [main_folder 'proteomes\uniprotkb_Caenorhabditis_elegans.fasta'];
        organism_name = 'Caenorhabditis_elegans';
    elseif contains(organism_name, 'Homo sapiens')
        proteome_path = [main_folder 'proteomes\uniprotkb_Homo_sapiens.fasta'];
        organism_name = 'Homo_sapiens';
    elseif contains(organism_name, 'Drosophila melanogaster')
        proteome_path = [main_folder 'proteomes\uniprotkb_Drosophila_melanogaster.fasta'];
        organism_name = 'Drosophila_melanogaster';
    elseif contains(organism_name, 'Nematostella vectensis')
        proteome_path = [main_folder 'proteomes\uniprotkb_Nematostella_vectensis.fasta'];
        organism_name = 'Nematostella_vectensis';
    elseif contains(organism_name, 'Amphimedon queenslandica')
        proteome_path = [main_folder 'proteomes\uniprotkb_Amphimedon_queenslandica.fasta'];
        organism_name = 'Amphimedon_queenslandica';
    end
end

function fastaLine = findSequenceInFasta(fastaFile, targetSequence)
    % Initialize the result
    fastaLine = [];
    % Open the FASTA file for reading
    fid = fopen(fastaFile, 'r');
    if fid == -1
        error('Unable to open the FASTA file.');
    end
    % Initialize variables to hold header and sequence
    header = '';
    sequence = '';
    inSequence = false;
    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        % Check for header line
        if line(1) == '>'
            % If we were in the previous sequence, check if it matches the target
            if inSequence && contains(sequence, targetSequence)
                fastaLine = sprintf('%s\n%s', header, sequence);
                fclose(fid);  % Close the file
                return;
            end

            % Update the header for the new sequence
            header = line;
            sequence = '';
            inSequence = false;
        else
            % Append the line to the current sequence
            sequence = [sequence, line];
            inSequence = true;
        end
    end
    % Check the last sequence in the file
    if inSequence && contains(sequence, targetSequence)
        fastaLine = sprintf('%s\n%s', header, sequence);
    end
    % Close the file
    fclose(fid);
    % If the target sequence was not found, fastaLine remains empty
end

function uniprot_id = get_uniprot_id(crnt_fasta_line)
    % Initialize the result
    uniprot_id = [];
    % Define the pattern for extracting the accession
    pattern = '([A-Z0-9]+)\|';
    % Use regexp to search for the pattern in the line
    match = regexp(crnt_fasta_line, pattern, 'tokens', 'once');
    % Extract the accession if a match is found
    if ~isempty(match)
        uniprot_id = match{1};
    end
end

function save_pdb(uniprot_id)
    % Downloads and saves the AlphaFold PDB structure file for the given UniProt accession.
    load('folder_paths.mat');
    global gene_group organism_name
    save_path = [data_folder gene_group '\' organism_name '\' uniprot_id '.pdb'];
    website_path = ['https://alphafold.ebi.ac.uk/files/AF-' uniprot_id '-F1-model_v4.pdb'];
    curl_command = ['curl -o ' save_path ' ' website_path];
    system(curl_command); 
end

function save_cif(uniprot_id)
    % Downloads and saves the AlphaFold CIF structure file for the given UniProt accession.
    load('folder_paths.mat');
    global gene_group organism_name
    save_path = [data_folder gene_group '\' organism_name '\' uniprot_id '.cif'];
    website_path = ['https://alphafold.ebi.ac.uk/files/AF-' uniprot_id '-F1-model_v4.cif'];
    curl_command = ['curl -o ' save_path ' ' website_path];
    system(curl_command); 
end

function save_path = save_fasta(uniprot_id)
    load('folder_paths.mat');
    global organism_name gene_group
    new_organism_name = strrep(organism_name, ' ', '_');
    save_path = [data_folder gene_group '\'  new_organism_name '\' uniprot_id '.fasta'];
    website_path = ['https://rest.uniprot.org/uniprotkb/' uniprot_id '.fasta'];
    curl_command = ['curl -o ' save_path ' ' website_path];
    system(curl_command); 
end

function save_degron_profile(fasta_path, uniprot_id)
    load('folder_paths.mat');
    global organism_name gene_group
    addpath 'G:\My_Drive\LabData\Shir\degron_lib\MATLAB\yeast_proteom\';
    save_path = [data_folder gene_group '\' organism_name '\' uniprot_id '.csv'];
    prob_calculator_with_fasta_input(fasta_path, save_path);
end
