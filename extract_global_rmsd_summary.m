% TMP_SCRIPT_FOR_GLOBAL_RMSD_TO_ALL_SPECIES  Extracts global RMSD and sequence alignment
% scores across all species from hogs_array into a flat summary table.
%
% Development script. Iterates over hogs_array and collects whole-protein RMSD
% and global sequence alignment scores for every protein group and species pair.
%
% Note: The nwalign call on line 35 is incomplete (4 output variables assigned but
% only one used). This script requires review before use in production.
%
% Requires in workspace: hogs_array, organism_list, zero_organism
% Produces:              table T in the workspace

num_hogs = length(hogs_array);  % Calculate the number of hogs based on the size of hogs_array
num_organisms = length(organism_list);  % Number of organisms (species)

protein_groups = {};
global_rmsds = [];
species_names = {};
seq_align = [];

for i = 1:num_hogs
    protein_group = hogs_array{i}.crnt_protein_group;
    for j = 1:num_organisms
        species_name = organism_list{j};
        rmsd_field_name = ['crnt_' species_name '_zero_global_rmsd_score'];
        crnt_zero_seq_name =  ['crnt_' zero_organism '_protein_sequence'];
        crnt_prot_seq_name =  ['crnt_' species_name '_protein_sequence'];
        crnt_zero_seq = hogs_array{i}.(crnt_zero_seq_name);
        try
            crnt_rmsd = hogs_array{i}.(rmsd_field_name);
        catch
            crnt_rmsd = NaN;
        end
        if isempty(crnt_rmsd)
            crnt_rmsd = NaN;
        end
        %sequence alignment
        try 
            crnt_prot_seq = hogs_array{i}.(crnt_prot_seq_name);
        catch
            crnt_prot_seq = [];
            crnt_align_score = NaN;
        end
        if ~isempty(crnt_prot_seq)
            crnt_align_score = nwalign(crnt_zero_seq, crnt_prot_seq);
            [alignedSeq1, alignedSeq2, crnt_align_score2, alignment] = nwalign(crnt_prot_seq,crnt_zero_seq)
        end

        protein_groups{end + 1} = protein_group; % Add protein group
        global_rmsds(end + 1) = crnt_rmsd; % Add the RMSD value
        species_names{end + 1} = species_name; % Add the species name
        seq_align(end + 1) = crnt_align_score;
    end
end
T = table(protein_groups', global_rmsds',seq_align', species_names', ...
    'VariableNames', {'crnt_protein_group', 'crnt_global_rmsd','crnt_seq_align', 'crnt_species'});