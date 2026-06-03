% GET_SCND_STRUCTURE_FULL_TABLE_SINGLE_USE  Pre-computes per-residue secondary structure tables.
%
% One-time preprocessing script. For each UniProt ID listed in uniprot_id.xlsx, reads the
% corresponding AlphaFold CIF file, extracts secondary structure assignments per residue
% (helix, sheet, loop), and writes a full-length table to data_all_scnd_struct/<id>.xlsx.
% These tables are consumed by degron_alignment_comparison.m.
%
% Requires on disk: all_cif_files/uniprot_id.xlsx  (column: protein_id)
%                   all_cif_files/<id>.cif          (AlphaFold CIF files)
%                   all_pdb_files/<id>.pdb          (for protein sequence length)
% Requires loaded:  folder_paths.mat
% Produces:         data_all_scnd_struct/<id>.xlsx for each UniProt ID

load('folder_paths.mat');
uniprot_table = readtable("G:\My_Drive\LabData\Shir\evolution_project\all_cif_files\uniprot_id.xlsx");
uniprot_ids = uniprot_table.protein_id;
number_of_ids = length(uniprot_ids);
local_path = [main_folder 'all_cif_files\'];
parfor i = 1:number_of_ids
    crnt_id = uniprot_ids{i};
    get_scnd_structure_full_table(crnt_id)
end

function get_scnd_structure_full_table(crnt_protein_id)
    % Load folder path
    load('folder_paths.mat');
    local_path_cif = [main_folder 'all_cif_files\'];
    local_path_pdb = [main_folder 'all_pdb_files\'];
    filename_cif = [local_path_cif crnt_protein_id '.cif'];
    filename_pdf = [local_path_pdb crnt_protein_id '.pdb'];
    if ~isfile(filename_cif)
        error('CIF file "%s" not found.', filename_cif);
    end
    try
        crnt_pdb = pdbread(filename_pdf);
        protein_length = crnt_pdb.Sequence.NumOfResidues;
    catch
        protein_length = 0;
    end
    secStructTable = parseSecondaryStructureFromCIF(filename_cif);
    scnd_struct_table = generate_full_table(secStructTable, protein_length);
    writetable(scnd_struct_table, [main_folder 'data_all_scnd_struct\' crnt_protein_id '.xlsx']);
end

function secStructTable = parseSecondaryStructureFromCIF(cifFile)
    % Reads a CIF file and extracts secondary structure info into a table.

    % Read file lines
    fid = fopen(cifFile, 'r');
    rawLines = {};
    tline = fgetl(fid);
    while ischar(tline)
        rawLines{end+1} = strtrim(tline); %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);

    % Initialize structure map
    residueStructMap = containers.Map('KeyType', 'int32', 'ValueType', 'char');

    % Regular expression to match: chain resName resNum ...
    expr = '(?<chain>\w+)\s+(?<resName>\w+)\s+(?<resNum>\d+)\s+\w+\s+\w+\s+\d+\s+(?<ss>\w+)';

    for i = 1:length(rawLines)
        line = rawLines{i};
        tokens = regexp(line, expr, 'names');
        if ~isempty(tokens)
            resNumStart = str2double(tokens(1).resNum);
            ssType = tokens(1).ss;

            % Some lines are ranges, e.g., A LEU 37 ... A ILE 53 ...
            tokensEnd = regexp(line, expr, 'once');
            if length(tokens) > 1
                resNumEnd = str2double(tokens(2).resNum);
            else
                resNumEnd = resNumStart;
            end

            for r = resNumStart:resNumEnd
                residueStructMap(r) = ssType;
            end
        end
    end

    % Build the output table
    allResidues = sort(cell2mat(residueStructMap.keys));
    residueNames = cell(length(allResidues),1);
    ssLabels = cell(length(allResidues),1);

    for i = 1:length(allResidues)
        resNum = allResidues(i);
        % Try to extract resName from raw lines (inefficient but ok for now)
        for j = 1:length(rawLines)
            if contains(rawLines{j}, sprintf(' %d ', resNum))
                % crude match, improve if needed
                parts = split(strtrim(rawLines{j}));
                idx = find(strcmp(parts, num2str(resNum)), 1);
                if idx > 1
                    residueNames{i} = parts{idx - 1};
                else
                    residueNames{i} = 'UNK';
                end
                break;
            end
        end
        ssLabels{i} = residueStructMap(resNum);
    end

    % Output table
    secStructTable = table(allResidues(:), residueNames, ssLabels, ...
        'VariableNames', {'ResidueNumber', 'ResidueName', 'SecondaryStructure'});
end

function scnd_struct_table = generate_full_table(partial_table, protein_length)    
    % Create columns
    ResidueNumber = (1:protein_length)';  % Column vector of residue indices from 1 to protein_length.
    ResidueName = repmat("unknown", protein_length, 1);          % String array of "unknown"
    SecondaryStructure = repmat("unknown", protein_length, 1);   % String array of "unknown"
    
    % Create the table
    scnd_struct_table = table(ResidueNumber, ResidueName, SecondaryStructure);
    for i = 2:height(partial_table)
            idx = partial_table.ResidueNumber(i);  % Get the ResidueNumber from partial_table
            
            % Check if idx is within the valid range
            if idx >= 1 && idx <= protein_length
                % Replace the corresponding row in scnd_struct_table
                scnd_struct_table.ResidueName(idx) = partial_table.ResidueName(i);
                scnd_struct_table.SecondaryStructure(idx) = partial_table.SecondaryStructure(i);
            end
        end
end

