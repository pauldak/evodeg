% AMINO_ACID_SLIDING_WINDOW  Generates sliding-window intervals on the reference organism.
%
% Inputs:
%   ds            - Data structure with protein metadata and pipeline parameters
%   interval_size - Window length in amino acids (typically 32)
%   overlap       - Residues shared between consecutive windows (typically 31,
%                   producing a single amino acid step)
%
% Output:
%   ds - Updated data structure; ds.crnt_<zero_organism>_intervals populated
%        with interval strings in 'start_end' format
function ds = amino_acid_sliding_window(ds, interval_size, overlap)
    number_of_iterations = length(ds.organism_list);
    number_of_iterations = 1; % Intervals are generated for the reference organism only; get_coupling_data maps them to all other species.
    for j = 1: number_of_iterations
        try
            crnt_organism = ds.organism_list{j};
            %%
            crnt_organism = ds.zero_organism; % Intervals are defined on the reference organism only.
            % Coupled intervals for each comparison species are derived by get_coupling_data
            % using the pairwise Ca distance matrix, rather than independent sliding windows.
            %%
            cnrt_protein_length = ds.(['crnt_' crnt_organism '_protein_length']);
            ds.(['crnt_' crnt_organism '_intervals']) = generate_Intervals();
        catch
            continue
        end
    end
    function intervals = generate_Intervals()
        intervals = cell(1, 0);  % Initialize an empty cell array to store intervals
        start = 1;
        while start <= cnrt_protein_length
            % Calculate the end of the current interval
            finish = start + interval_size - 1;
            % Check if the finish of the current interval exceeds the number
            if finish > cnrt_protein_length
                finish = cnrt_protein_length;
                break
            end
            % Add the interval to the cell array
            intervals{end+1} = [num2str(start) '_' num2str(finish)];
            % Calculate the start of the next interval with overlap
            start = start + overlap - 1;
        end
    end
end