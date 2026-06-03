% ORTHO_DB_SEARCH_RESULT_EXTRACTER  Downloads OrthoDB ortholog FASTA files for all gene groups.
%
% One-time data acquisition script. Reads OrthoDB search result group IDs from
% search_results.txt, then downloads a multi-species FASTA file for each group
% from the OrthoDB v11 API (species: H. sapiens, D. rerio, D. melanogaster,
% N. vectensis, C. elegans, A. queenslandica, S. cerevisiae).
%
% Requires on disk: orthodb_data/search_results.txt (OrthoDB search result export)
% Produces:         orthodb_data/<group_id>.fasta for each gene group

global main_folder_path
main_folder_path = 'G:\Dropbox\my_documents\School\Academic\PhD\Year_A\Biological_Systematics_90721\final_project\';
filename = [main_folder_path 'orthodb_data\search_results.txt'];
groupValues = extractGroupValues(filename);
number_of_results = length(groupValues);

for i = 1: number_of_results
    crnt_group_id = groupValues{i};
    save_search_results(crnt_group_id);
end

function groupValues = extractGroupValues(filename)
    % Read the content of the file
    fileID = fopen(filename, 'r');
    fileContent = textscan(fileID, '%s', 'Delimiter', '\n');
    fclose(fileID);
    % Initialize a cell array to store the extracted values
    groupValues = {}; 
    % Iterate through the file content and extract values followed by "Group"
    for i = 1:length(fileContent{1})
        line = fileContent{1}{i};
        % Use regular expression to match values followed by "Group"
        match = regexp(line, '(\w+) at Eukaryota level', 'tokens');
        if ~isempty(match)
            groupValues = [groupValues, match{1}{1}];
        end
    end
end

function save_search_results(ortho_group_id)
    global main_folder_path
    save_path = [main_folder_path 'orthodb_data\'  ortho_group_id '.fasta'];
    website_path = ['https://data.orthodb.org/v11/fasta?id=' ortho_group_id '&species=9606_0,7955_0,7227_0,45351_0,6239_0,400682_0,559292_0'];
    %curl_command = ['curl -o ' save_path ' ' website_path];
    website_data = webread(website_path);
    %saving to file   
    filename_out = save_path;      
    fid = fopen(filename_out, 'w'); %// open file to writing
    fprintf(fid, '%s\n', string(website_data)); %// print string to file
    fclose(fid); %// don't forget to close the file
    %system(curl_command); 
end
