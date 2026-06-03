% SET_ALIGNMENT_SCORE  Computes local sequence alignment scores for all coupled interval pairs.
%
% Inputs:
%   ds           - Data structure with protein sequences, interval lists, and scoring matrix
%   analyze_type - 'all' for sliding-window intervals; 'deg' for degron intervals
%
% Output:
%   ds - Updated data structure:
%        analyze_type 'all': ds.crnt_<organism>_alignment_score appended
%        analyze_type 'deg': ds.crnt_<organism>_degron_alignment_score appended
%
% Uses MATLAB's localalign with the QCDPred-derived substitution matrix stored in
% ds.qcdpred_penalty_matrix. Intervals flagged as '-1_-1' receive the sentinel
% score ds.nan_alignment_score (99).
% SET_ALIGNMENT_SCORE  Computes local sequence alignment scores for all coupled interval pairs.
%
% Inputs:
%   ds           - Data structure with protein sequences, interval lists, and scoring matrix
%   analyze_type - 'all' for sliding-window intervals; 'deg' for degron intervals
%
% Output:
%   ds - Updated data structure:
%        analyze_type 'all': ds.crnt_<organism>_alignment_score appended
%        analyze_type 'deg': ds.crnt_<organism>_degron_alignment_score appended
%
% Alignment matrix is controlled by align_to_use in folder_paths.mat:
%   1 = standard BLOSUM matrix; 2 = QCDPred-derived custom matrix (default).
% Intervals flagged as '-1_-1' receive the sentinel score ds.nan_alignment_score (99).
function ds = set_alignment_score(ds, analyze_type)
    load('folder_paths.mat');
    if strcmp(analyze_type, 'all')
        crnt_intervals = ds.(['crnt_' ds.zero_organism '_intervals']);
    elseif strcmp(analyze_type, 'deg')
        crnt_intervals = ds.(['crnt_' ds.zero_organism '_degrons_intervals']);
    end
    number_of_intervals = length(crnt_intervals);
    for j = 1: length(ds.organism_list)
        for i = 1:number_of_intervals
            try
                crnt_organism = ds.organism_list{j};
                if strcmp(analyze_type, 'all')
                    crnt_zero_organism_interval = ds.(['crnt_' ds.zero_organism '_intervals']){i};
                    crnt_organism_interval = ds.(['crnt_' crnt_organism '_intervals']){i};
                elseif strcmp(analyze_type, 'deg')
                    crnt_zero_organism_interval = ds.(['crnt_' ds.zero_organism '_degrons_intervals']){i};
                    crnt_organism_interval = ds.(['crnt_' crnt_organism '_degrons_intervals']){i};
                end
                if strcmp(crnt_organism_interval, ds.nan_interval) %case '-1_-1'
                   ds = add_alignment_score_to_ds(crnt_organism, analyze_type, ds.nan_alignment_score, ds);
                   continue;
                end
                crnt_fasta_seq = ds.(['crnt_' crnt_organism '_protein_sequence']);
                %interpreter of interval string
                parts = strsplit(crnt_organism_interval, '_');
                start_val = str2double(parts{1});
                end_val = str2double(parts{2});
                crnt_fasta_seq_partial = crnt_fasta_seq(start_val : end_val);
                % all again for the zero organism
                crnt_zero_fasta_seq = ds.(['crnt_' ds.zero_organism '_protein_sequence']);
                %interpreter of interval string
                parts = strsplit(crnt_zero_organism_interval, '_');
                zero_start_val = str2double(parts{1});
                zero_end_val = str2double(parts{2});
                
                crnt_zero_fasta_seq_partial = crnt_zero_fasta_seq(zero_start_val : zero_end_val);
                %alignment
                if align_to_use == 1
                    crnt_align_struct = localalign(crnt_fasta_seq_partial, crnt_zero_fasta_seq_partial);
                else
                    crnt_align_struct = localalign(crnt_fasta_seq_partial, crnt_zero_fasta_seq_partial, 'ScoringMatrix', ds.qcdpred_penalty_matrix);
                end


                ds = add_alignment_score_to_ds(crnt_organism, analyze_type, crnt_align_struct.Score, ds);
            catch
                continue
            end
        end
    end
end

function ds = add_alignment_score_to_ds(crnt_organism, analyze_type, score, ds)
    if strcmp(analyze_type, 'all')
        ds.(['crnt_' crnt_organism '_alignment_score']) = ...
        [ds.(['crnt_' crnt_organism '_alignment_score']), score];
    elseif strcmp(analyze_type, 'deg')
        ds.(['crnt_' crnt_organism '_degron_alignment_score']) = ...
        [ds.(['crnt_' crnt_organism '_degron_alignment_score']), score];
    end
end